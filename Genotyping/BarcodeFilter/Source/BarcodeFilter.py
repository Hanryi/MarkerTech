"""
Created on August 28, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 BarcodeFilter.py <i_barcode.txt> <i_1.fq> <i_2.fq> <o_fq.dir>
# <i_barcode.txt>: file containing 8bp barcodes
# <i_1.fq>: forward sequencing file
# <i_2.fq>: reverse sequencing file
# <o_fq.dir>: directory for saving each pair of fq files

@Function
Divide samples from sequencing files by barcode sequences
"""

import os
import sys

barcodeFile = open(sys.argv[1], 'r')
barcode = [bc.strip() for bc in barcodeFile]
barcodeFile.close()

fFq = open(sys.argv[2], 'r')
rFq = open(sys.argv[3], 'r')
fqDir = sys.argv[4]

# Definition of variables and functions
fReads = []
rReads = []
fRow = '1'
rRow = '2'
sample = 0
line = 0
cache = 0
emptySpl = list(range(len(barcode)))

# For case count
bothMatched = 0
onceMatched = 0
noneMatched = 0


def num_case(f_row, r_row, line_ID):
    """
    Select the case of once "1"-"2" occur in F/R ID row respectively
    lack of memory property for only one pair IDs

    :param f_row: (str)
        forward reads ID
    :param r_row: (str)
        reverse reads ID
    :param line_ID: (int)
        line number of reads ID

    :return: boolean
    """
    # Case of unequal length
    if len(f_row) != len(r_row):
        print(f"Mismatched: unequal length in {line_ID} line ({f_row}, {r_row})")
        return False

    # Case of equal length
    case = 0
    for letter in range(len(f_row)):
        # Sub-condition of next "if"
        if f_row[letter] == '1' and r_row[letter] == '2':
            case += 1
            continue
        if f_row[letter] != r_row[letter]:
            print(f"Mismatched: unequal letters in {line_ID} line ({f_row}, {r_row})")
            return False

    # Sub-condition of first "if" in "for"
    if case == 1:
        # Regular output
        return True
    if case != 1:
        print(f"Mismatched: 1-2 occurs {case} times in {line_ID} line ({f_row}, {r_row})")
        return False


while fRow and rRow:

    try:

        # Iterate by row and touch exception
        fRow = next(fFq)
        rRow = next(rFq)
        line += 1

        # May be removed in future
        # Check if reads ID are pair-matched in F/R sequencing files (optional)
        if (line + 3) % 4 == 0 and not num_case(fRow, rRow, line):
            break

        # Limited cache variables and recycle them when they stored 4 elements
        if (line + 3) % 4 == 0:

            # Receive active cache state and turn off it after file writen
            if cache == 1:
                # Write reads into correspond file
                fq_1 = open(os.path.join(fqDir, str(sample) + "_1.fq"), 'a')
                fq_2 = open(os.path.join(fqDir, str(sample) + "_2.fq"), 'a')
                fq_1.write(''.join(fReads))
                fq_2.write(''.join(rReads))
                fq_1.close()
                fq_2.close()

                cache = 0

            fReads = []
            rReads = []

        fReads.append(fRow)
        rReads.append(rRow)

        # Focus on sequence
        if (line + 2) % 4 == 0:

            if fRow[: 8] in barcode and rRow[: 8] in barcode:
                bothMatched += 1
                continue

            # Output reads
            if fRow[: 8] in barcode or rRow[: 8] in barcode:
                onceMatched += 1

                if fRow[: 8] in barcode:
                    sample = barcode.index(fRow[: 8])
                if rRow[: 8] in barcode:
                    sample = barcode.index(rRow[: 8])
                emptySpl[sample] = -1

                # Active file writing action
                cache = 1

            if fRow[: 8] not in barcode and rRow[: 8] not in barcode:
                noneMatched += 1

    except StopIteration:
        break

fFq.close()
rFq.close()

# Create empty file for not found samples
for i in emptySpl:
    if i == -1:
        continue
    empty_1 = open(os.path.join(fqDir, str(i) + "_1.fq"), 'a')
    empty_2 = open(os.path.join(fqDir, str(i) + "_2.fq"), 'a')
    empty_1.close()
    empty_2.close()

# Terminal output
reads = bothMatched + onceMatched + noneMatched
print(
    f"bothMatched: {bothMatched} ({100 * bothMatched / reads:.2f}%), "
    f"onceMatched: {onceMatched} ({100 * onceMatched / reads:.2f}%), "
    f"noneMatched: {noneMatched} ({100 * noneMatched / reads:.2f}%)"
)
