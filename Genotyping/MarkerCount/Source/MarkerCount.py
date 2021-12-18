"""
Created on August 31, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 MarkerCount.py <barcode.txt> <i_mut.table> <i_bam.dir> <i_fq.dir> <o_csv.dir>
# <barcode.txt>: barcode sequence separated with '\n'
# <i_mut.table>: table containing mutation info (CHROM, POS, REF)
#                1. CHROM: contig of reference, usually equals to the number of chromosome
#                2. POS: SNP position in chromosome
#                3. REF: base character in reference genome
# <i_bam.dir>: directory of bam
# <i_fq.dir>: directory of fq
# <o_csv.dir>: directory used to save result csv file for each site (Num, Frequency, UMI)
#              1. Num: number of sample which has a record on the site
#              2. Frequency: appearance probability of the mut base provided by user
#              3. UMI: number of reads with unique UMI marker

@Function
1. Deduplicate UMI marker for each loci in bam file
2. Allow user inputs mutation info with xlsx, txt or csv format file
"""

import os
import re
import sys

import openpyxl as opx
import pandas as pd
import pysam

barcode = sys.argv[1]
mutFile = sys.argv[2]
bamDir = sys.argv[3]
fqDir = sys.argv[4]
csvDir = sys.argv[5]

# Definition of variable and function
baseMap = {"A": 0, "C": 1, "G": 2, "T": 3}
umiObj = re.compile("ATAGCGACGCGTTTCAAC([ACGT]{6})")
mergeFrame = {}

# Measure for sample number by barcodes
sample = 0
with open(barcode, 'r') as barcodeFile:
    for bc in barcodeFile:
        if bc.strip() == '':
            continue
        sample += 1


def xlsx_reader(mut_file):
    """
    Iterate row array in mutation file
    (read_only=True only allow rows attr but not columns attr)

    :param mut_file: str
        file path of xlsx

    :return: list
        iterable and each element in array is str type
    """
    mut_mat = []
    mut_book = opx.load_workbook(mut_file, read_only=True)
    mut_sht = mut_book[mut_book.sheetnames[0]]
    mut_row = mut_sht.rows
    # Iterate in generator
    next(mut_row)
    for row in mut_row:
        mut_arr = []
        for x in row:
            mut_arr.append(str(x.value))
        mut_mat.append(mut_arr)
    return mut_mat


def text_reader(mut_file):
    """
    Iterate row array in mutation file

    :param mut_file: str
        file path of csv or txt(separate should be ',')

    :return: list
        iterable and each element in array is str type
    """
    mut_mat = []
    mut_obj = open(mut_file, 'r')
    next(mut_obj)
    for x in mut_obj:
        mut_mat.append(x.strip().split(','))
    mut_obj.close()
    return mut_mat


def seq_name(fq_1, fq_2):
    """
    Build mapping between serial name and its sequence for a PE fq

    :param fq_1: (str)
        sample_1.fq filepath
    :param fq_2: (str)
        sample_2.fq filepath

    :return: (dict)
        mapping for each reads serial name and its sequence
    """
    seq_alias = {}
    f_fq = open(fq_1)
    r_fq = open(fq_2)
    f_row = '1'
    r_row = '2'
    line = 1

    while f_row and r_row:

        try:

            f_row = next(f_fq).strip()
            r_row = next(r_fq).strip()

            # Focus on serial name
            if (line + 3) % 4 == 0:

                # PE index and seq pair initialize avoiding generate memory
                pe_idx = 0
                seq_pair = []

                # Update the pe_idx
                for letter in range(len(f_row) - 1, -1, -1):
                    if f_row[letter] == '1' and r_row[letter] == '2':
                        pe_idx = letter
                        break

                # Generate by next(generator)
                line += 1
                # Update the seq_pair
                seq_pair.append(next(f_fq).strip())
                seq_pair.append(next(r_fq).strip())

                # The key value should equal to query_name in bam
                seq_alias[f_row[1: pe_idx - 1]] = seq_pair

            # Generate for next() in next cycle
            line += 1

        except StopIteration:
            break

    f_fq.close()
    r_fq.close()

    return seq_alias


# Build iterable object for 3 filetypes
filetype = mutFile.split('.')[-1]
mutation = []
if filetype.lower() not in ["xlsx", "csv", "txt"]:
    print("Invalid filetype")
if filetype.lower() == "xlsx":
    mutation = xlsx_reader(mutFile)
if filetype.lower() == "csv" or filetype.lower() == "txt":
    mutation = text_reader(mutFile)

for spl in range(sample):

    bamObj = pysam.AlignmentFile(
        os.path.join(bamDir, str(spl) + ".sorted.bam"), 'rb'
    )
    # Check index file
    if not bamObj.check_index():
        continue

    # Sorted file for search
    seqAlias = seq_name(
        os.path.join(fqDir, str(spl) + "_1.fq"),
        os.path.join(fqDir, str(spl) + "_2.fq")
    )

    for mutArr in mutation:

        # Build dataframe
        location = "_".join(mutArr)
        if location not in mergeFrame.keys():
            mergeFrame[location] = [[], [], []]

        # Fetch mutate positions in segment iterator
        snp = int(mutArr[1])
        segPile = bamObj.fetch(mutArr[0], snp - 1, snp)
        umiPile = [[] for _ in 'ACGT']

        for seg in segPile:

            # Search SNP base by position number array corresponding to sequence
            refPos = seg.get_reference_positions()
            if snp - 1 not in refPos:
                continue
            # SNP
            pos = [x[0] for x in seg.get_aligned_pairs() if snp - 1 == x[1]][0]
            snpBase = seg.query_sequence[pos]
            if snpBase == 'N':
                continue

            # Name --> Seq(fq)
            seqPair = seqAlias[seg.query_name]
            fSeq = seqPair[0]
            rSeq = seqPair[1]

            # Seq --> UMI
            mergeUmi = umiObj.findall(fSeq) + umiObj.findall(rSeq)
            if len(mergeUmi) != 1:
                continue
            umi = mergeUmi[0]

            # SNP: UMI
            if umi not in sum(umiPile, []):
                umiPile[baseMap[snpBase]].append(umi)

        # Frequency
        umiSet = sum(umiPile, [])
        if len(umiSet) == 0:
            continue
        freq = len(umiPile[baseMap[mutArr[2]]]) / len(umiSet)

        # Append to dataframe
        mergeFrame[location][0].append(spl)
        mergeFrame[location][1].append(freq)
        mergeFrame[location][2].append(len(umiSet))

for loci in mergeFrame.keys():
    frame = pd.DataFrame({
        "Num": mergeFrame[loci][0],
        "Frequency": mergeFrame[loci][1],
        "UMI": mergeFrame[loci][2]
    })
    frame.to_csv(os.path.join(csvDir, loci) + ".csv", index=False, sep=',')
