"""
Created on August 31, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 AltFraction.py <bc.txt> <mut.tab> <bam.dir> <umi.dir> <o_csv.dir>
# <bc.txt>: 8bp barcode sequences separated with '\n'
# <mut.tab>: table containing mutation info (CHROM, POS, REF, ALT)
#                1. CHROM: contig of reference
#                2. POS: position of the first base in mutation identifier
#                3. REF: bases in reference sequence
#                4. ALT: mutant bases in reads
# <bam.dir>: directory of bam files
# <umi.dir>: directory of fq files containing umi sequence in fastq format
# <o_csv.dir>: a directory used to save tmp csv (Num, Frequency, UMI)
#              1. Num: number of sample which has a record at the site
#              2. Frequency: fraction of the mutation in the position
#              3. UMI: number of reads with unique UMI marker

@Function
1. Deduplicate reads with same UMI marker for each loci
2. Count mutant SNP and InDel frequency among reads for each loci
"""

import os
import sys

import openpyxl as opx
import pandas as pd
import pysam


class MutMatrix:
    """ Mutation Matrix class.

    This class mainly provide a method for external calls to update mutation
    matrix which deals with multi-mutation sites.
    """
    def __init__(self, mut_file: str):
        self.path = mut_file
        self.type = mut_file.split('.')[-1].lower()

    def update_mutation(self):
        """ Update each matrix in tensor and merge into a matrix.

        :return: (2d list)
            updated mutation matrix
        """
        update_mat = []
        chr_tensor = self.sort_mutation()
        for mat in chr_tensor:
            update_mat += self.update_sub_matrix(mat)
        return update_mat

    def sort_mutation(self):
        """ Divide sorted matrix into tensor by chromosome names.

        :return: (3d list)
            sorted tensor
        """
        mut_matrix = []
        if self.type not in ["csv", "txt", "xlsx"]:
            print("Invalid filetype")
        if self.type == "xlsx":
            mut_matrix = self.xlsx_reader()
        if self.type == "csv" or self.type == "txt":
            mut_matrix = self.text_reader()
        # Sort raw matrix
        mut_matrix.sort(key=lambda info: (info[0], eval(info[1])))

        # Divide matrix into tensor by chromosome names
        chr_tensor = []
        chr_name = ""
        for vec in mut_matrix:
            # Divide list regions
            if vec[0] != chr_name:
                chr_name = vec[0]
                chr_tensor.append([])
                setattr(self, "overlap_" + chr_name, set())
            chr_tensor[-1].append(vec)
        return chr_tensor

    def update_sub_matrix(self, sorted_mat):
        """ Merge possible overlapping sites in a sorted sub_matrix.

        This class method scan the overlapping sites, update the REF and ALT
        sequence and set the following attributes:
        ---------
        self.overlap_<chr>: save row index avoiding use them again. (update)
        self."_".join(<vec>): length of union sequence. (dynamically created)
        ---------

        :return: (2d list)
            updated matrix, each row output a frequency table
        """
        # No double mutation
        if len(sorted_mat) == 1:
            return sorted_mat

        chr_attr = "overlap_" + sorted_mat[0][0]
        update_mat = []
        for i in range(len(sorted_mat) - 1):
            mut_size = len(sorted_mat[i][2])
            # Scan overlapped mutation sites
            if eval(sorted_mat[i][1]) + mut_size >= eval(sorted_mat[i + 1][1]):
                # Add overlapped index to set
                setattr(
                    self, chr_attr,
                    getattr(self, chr_attr) | {i, i + 1}
                )
                union = self.merge_seq(
                    sorted_mat[i][2], sorted_mat[i + 1][2],
                    eval(sorted_mat[i][1]), eval(sorted_mat[i + 1][1])
                )
                head_len = eval(sorted_mat[i + 1][1]) - eval(sorted_mat[i][1])
                merge_vec = [
                    sorted_mat[i][0], sorted_mat[i][1],
                    sorted_mat[i][3] + union[len(sorted_mat[i][2]):],
                    union[: head_len] + sorted_mat[i + 1][3] +
                    union[head_len + len(sorted_mat[i + 1][2]):]
                ]
                # Update mutation vector and add attribute
                update_mat.append(merge_vec)
                setattr(self, "_".join(merge_vec), len(union))

            if i not in getattr(self, chr_attr):
                update_mat.append(sorted_mat[i])
        # Update the last mutation vector
        if len(sorted_mat) - 1 not in getattr(self, chr_attr):
            update_mat.append(sorted_mat[-1])
        return update_mat

    @staticmethod
    def merge_seq(seq_1, seq_2, pos_1, pos_2) -> str:
        """ Merge two overlapping sequences.

        :param seq_1: (str)
            overlapping region is at the end of the seq_1
        :param seq_2: (str)
            overlapping region is at the head of the seq_2
        :param pos_1: (int)
            first base position of seq_1
        :param pos_2: (int)
            first base position of seq_2, and
            pos_2 > pos_1

        :return: (str)
            merged sequence
        """
        if pos_1 + len(seq_1) >= pos_2 + len(seq_2):
            return seq_1
        if pos_1 + len(seq_1) < pos_2 + len(seq_2):
            return seq_1[: pos_2 - pos_1] + seq_2

    def xlsx_reader(self):
        """ Read mutation args by line from a xlsx file.

        :return: (2d list)
            mutation matrix composed of strings
        """
        mut_mat = []
        mut_book = opx.load_workbook(self.path, read_only=True)
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

    def text_reader(self):
        """ Read mutation args by line from a text file (separator=',').

        :return: (2d list)
            mutation matrix composed of strings
        """
        mut_mat = []
        mut_obj = open(self.path)
        next(mut_obj)
        for x in mut_obj:
            mut_mat.append(x.strip().split(','))
        mut_obj.close()
        return mut_mat


def umi_roster(umi_fq):
    """ Associate serial name and its UMI sequence from a fq file.

    The key of dict is a substring which starts from the character behind '@'
    to the character in two bases front of the pair-end sequencing file
    identifier. The read IDs all are valid, read with invalid ID has been
    filtered out before.

    :param umi_fq: (str)
        sample.fq filepath

    :return: (dict)
        a dict links serial name and its UMI sequence like {name: sequence,}
    """
    roster = {}
    line = 1

    fq_obj = open(umi_fq)
    while True:
        try:
            fq_row = next(fq_obj).strip()
            if (line + 3) % 4 == 0:
                # Link read name and its sequence
                roster[fq_row[: -2]] = next(fq_obj).strip()
                line += 1
            line += 1

        except StopIteration:
            break
    fq_obj.close()
    return roster


def mut_locator(seg_seq, leader_pos, mut_len, pos_arr):
    """ Locate the sequence in mutation area.

    The return sequence can be Ref, Alt or any other seq. But only the Ref and
    Alt will be used for calculate frequency (Freq = Alt / (Ref + Alt)) after
    UMI deduplication.

    :param seg_seq: (str)
        segment sequence in bam
    :param leader_pos: (int)
        position index of the first base in the mutation identifier
    :param mut_len: (int)
        maximum length between Ref and Alt
    :param pos_arr: (tuple list)
        each binary tuple shows the binary position of current base that first
        int element means reads position and second int element means reference
        position

    :return: (str)
        mutation sequence in the given mutation length
    """
    # Locate the position of leader base
    mut_leader = None
    mut_follow = None
    for pos in pos_arr:
        if pos[0] is not None and pos[1] == leader_pos:
            mut_leader = pos_arr.index(pos)
        if pos[0] is not None and pos[1] == leader_pos + mut_len:
            mut_follow = pos_arr.index(pos)
            break

    # Mutation sequence
    mut_seq = ''
    mut_arr = pos_arr[mut_leader: mut_follow]
    for mut in mut_arr:
        if mut[0] is None:
            continue
        mut_seq += seg_seq[mut[0]]
    return mut_seq


barcode = sys.argv[1]
mutFile = sys.argv[2]
bamDir = sys.argv[3]
umiDir = sys.argv[4]
# Output
csvDir = sys.argv[5]

mutObj = MutMatrix(mutFile)
# Read mutation table as updated matrix
mutMat = mutObj.update_mutation()

# Measure sample number by barcodes
sample = 0
with open(barcode) as barcodeFile:
    for bc in barcodeFile:
        if bc.strip() == '':
            continue
        sample += 1

statePanel = {}
for spl in range(sample):
    # Load bam file
    bamObj = pysam.AlignmentFile(
        os.path.join(bamDir, str(spl) + ".sorted.bam"), 'rb'
    )
    # Check index file
    if not bamObj.check_index():
        continue

    # Load UMI fastq file
    seqRoster = umi_roster(os.path.join(umiDir, str(spl) + ".fq"))
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
        # Numerical rewriting
        if hasattr(mutObj, location):
            mutLen = getattr(mutObj, location)
        # UMI container (Ref and Alt)
        umiPanel = [[], []]

        # Fetch interval endpoints are not equal
        segPile = bamObj.fetch(mutArr[0], vcfPos - 1, vcfPos + len(vcfRef))
        for seg in segPile:
            """ Mutation state.
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

            """ UMI deduplication.
            """
            umiSeq = seqRoster[seg.query_name]
            if umiSeq not in sum(umiPanel, []):
                umiPanel[mutState].append(umiSeq)

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
