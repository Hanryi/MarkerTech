"""
Created on August 31, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 MarkerCount.py [bc.txt] [mut.tab] [bam.dir] [fq.dir] [o_csv.dir]
# [bc.txt]: 8bp barcode sequences separated with '\n'
# [mut.tab]: table containing mutation info (CHROM, POS, REF, ALT)
#                1. CHROM: contig of reference
#                2. POS: position of the first base in mutation identifier
#                3. REF: bases in reference sequence
#                4. ALT: mutant bases in reads
# [bam.dir]: directory of bam
# [fq.dir]: directory of fq
# [o_csv.dir]: a directory used to save tmp csv (Num, Frequency, UMI)
#              1. Num: number of sample which has a record at the site
#              2. Frequency: fraction of the mutation in the position
#              3. UMI: number of reads with unique UMI marker

@Function
1. Deduplicate reads with same UMI marker for each loci
2. Count mutant SNP and InDel frequency among reads for each loci
"""

import os
import re
import sys

import openpyxl as opx
import pandas as pd
import pysam


def xlsx_reader(mut_file):
    """ Read mutation args by line from xlsx file.

    :param mut_file: str
        xlsx path

    :return: 2d list
        mutation matrix composed of strings
    """
    mut_mat = []
    mut_book = opx.load_workbook(mut_file, read_only=True)
    mut_sht = mut_book[mut_book.sheetnames[0]]
    # read_only=True only allow rows attr but not columns attr
    mut_row = mut_sht.rows
    next(mut_row)
    for row in mut_row:
        mut_arr = []
        for x in row:
            mut_arr.append(str(x.value))
        mut_mat.append(mut_arr)
    return mut_mat


def text_reader(mut_file):
    """ Read mutation args by line from text file.

    :param mut_file: str
        path of csv or txt (separator=',')

    :return: 2d list
        mutation matrix composed of strings
    """
    mut_mat = []
    mut_obj = open(mut_file)
    next(mut_obj)
    for x in mut_obj:
        mut_mat.append(x.strip().split(','))
    mut_obj.close()
    return mut_mat


def seq_roster(fq_1, fq_2):
    """ Associate serial name and its sequence from PE fq files.

    :param fq_1: (str)
        sample_1.fq filepath
    :param fq_2: (str)
        sample_2.fq filepath

    :return: (dict)
        use serial name as key and its sequence as key value
    """
    f_fq = open(fq_1)
    r_fq = open(fq_2)

    # Init
    roster = {}
    line = 1
    while True:
        try:
            f_row = next(f_fq).strip()
            r_row = next(r_fq).strip()

            # Focus on identity line
            if (line + 3) % 4 == 0:
                pair_idx = 0
                seq_pair = []

                # Search the identifier index in reverse
                for letter in range(len(f_row) - 1, -1, -1):
                    if f_row[letter] == '1' and r_row[letter] == '2':
                        pair_idx = letter
                        break

                # Add pair of sequences in next line
                line += 1
                seq_pair.append(next(f_fq).strip())
                seq_pair.append(next(r_fq).strip())

                # The query_name decided by identifier index
                roster[f_row[1: pair_idx - 1]] = seq_pair

            # Lead to next cycle
            line += 1

        except StopIteration:
            break

    f_fq.close()
    r_fq.close()
    return roster


def mut_locator(seg_seq, leader_pos, mut_len, pos_arr):
    """ Locate the sequence in mutation area.

    The return sequence can be Ref, Alt or any other seq. But only the Ref and
    Alt will be used for calculate frequency (Freq = Alt / (Ref + Alt)) after
    UMI deduplication.

    :param seg_seq: str
        segment sequence in bam
    :param leader_pos: int
        position index of the first base in the mutation identifier
    :param mut_len: int
        maximum length between Ref and Alt
    :param pos_arr: tuple list
        each binary tuple shows the binary position of current base that first
        int element means reads position and second int element means reference
        position

    :return: str
        mutation sequence in the given mutation length
    """
    # Locate the position of leader base
    mut_leader = None
    for pos in pos_arr:
        if pos[0] is not None and pos[1] == leader_pos:
            mut_leader = pos_arr.index(pos)
            break
    # Avoid invalid leader base (just in case)
    if mut_leader is None:
        return None

    # Mutation sequence
    mut_seq = ''
    mut_arr = pos_arr[mut_leader: mut_leader + mut_len]
    for mut in mut_arr:
        if mut[0] is None:
            continue
        mut_seq += seg_seq[mut[0]]
    return mut_seq


# Args
barcode = sys.argv[1]
mutFile = sys.argv[2]
bamDir = sys.argv[3]
fqDir = sys.argv[4]
csvDir = sys.argv[5]

# Definition of regex and variable
UmiRegex = re.compile("ATAGCGACGCGTTTCAAC([ACGT]{6})")
statePanel = {}

""" Mutation matrix (n * 4)
"""
filetype = mutFile.split('.')[-1]
mutMat = []
# TODO(154702913@qq.com): Add VCF format reader.
if filetype.lower() not in ["xlsx", "csv", "txt"]:
    print("Invalid filetype")
if filetype.lower() == "xlsx":
    mutMat = xlsx_reader(mutFile)
if filetype.lower() == "csv" or filetype.lower() == "txt":
    mutMat = text_reader(mutFile)

# Measure sample number by barcodes
sample = 0
with open(barcode) as barcodeFile:
    for bc in barcodeFile:
        if bc.strip() == '':
            continue
        sample += 1

for spl in range(sample):
    # Load bam file
    bamObj = pysam.AlignmentFile(
        os.path.join(bamDir, str(spl) + ".sorted.bam"), 'rb'
    )
    # Check index file
    if not bamObj.check_index():
        continue

    # Load fastq file
    seqRoster = seq_roster(
        os.path.join(fqDir, str(spl) + "_1.fq"),
        os.path.join(fqDir, str(spl) + "_2.fq")
    )

    for mutArr in mutMat:
        # Sequentially check each key of loci
        location = "_".join(mutArr)
        if location not in statePanel.keys():
            statePanel[location] = [[], [], []]

        # Mutation args
        vcfPos = int(mutArr[1])
        vcfRef = mutArr[2]
        vcfAlt = mutArr[3]
        mutLen = max(len(vcfRef), len(vcfAlt))

        # UMI container (Ref and Alt)
        umiPanel = [[], []]

        # Fetch interval endpoints are not equal
        segPile = bamObj.fetch(mutArr[0], vcfPos - 1, vcfPos + len(vcfRef))
        for seg in segPile:
            """ Mutation sequence.
            """
            mutSeq = mut_locator(
                seg.query_sequence, vcfPos - 1, mutLen, seg.get_aligned_pairs()
            )
            # Filter Ref or Alt sequence
            if mutSeq != vcfRef and mutSeq != vcfAlt or mutSeq is None:
                continue
            mutState = 0
            if mutSeq == vcfAlt:
                mutState = 1

            """ UMI calling.
            """
            # Name --> Fq
            seqPair = seqRoster[seg.query_name]
            fSeq = seqPair[0]
            rSeq = seqPair[1]
            # Fq --> UMI
            umiCase = UmiRegex.findall(fSeq) + UmiRegex.findall(rSeq)
            if len(umiCase) != 1:
                continue
            umi = umiCase[0]

            # SNP: UMI
            if umi not in sum(umiPanel, []):
                umiPanel[mutState].append(umi)

        # Frequency
        umiPool = sum(umiPanel, [])
        if len(umiPool) == 0:
            continue
        freq = len(umiPanel[1]) / len(umiPool)

        # Append to dataframe
        statePanel[location][0].append(spl)
        statePanel[location][1].append(freq)
        statePanel[location][2].append(len(umiPool))

for loci in statePanel.keys():
    frame = pd.DataFrame({
        "Num": statePanel[loci][0],
        "Frequency": statePanel[loci][1],
        "UMI": statePanel[loci][2]
    })
    frame.to_csv(os.path.join(csvDir, loci) + ".csv", index=False, sep=',')
