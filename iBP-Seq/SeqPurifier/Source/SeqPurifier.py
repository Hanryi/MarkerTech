"""
Created on August 28, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 SeqPurifier.py <bc.txt> <cln_1.fq> <cln_2.fq> <o_fq.dir> <o_umi.dir>
# <bc.txt>: file containing 8bp barcodes
# <cln_1.fq>: clean data of forward sequencing
# <cln_2.fq>: clean data of reverse sequencing
# <o_fq.dir>: directory for storing pair of fq which removed additional seq
# <o_umi.dir>: directory for storing fq which contains only UMI fragments

@Function
1. Divide mixed sequencing data to several pairs of fq by barcode and bridge
2. Separate UMI sequence in reads from genome sequence
"""

import os
import sys


def valid_tag(f_row, r_row, line_ID):
    """ Judge the validity of current ID pair.

    Valid ID should be judged by identifier which only allowed occur once as
    "1"-"2" in F/R ID respectively.

    :param f_row: (str)
        forward reads ID
    :param r_row: (str)
        reverse reads ID
    :param line_ID: (int)
        line number of reads ID

    :return: (bool)
    """
    # Case of unequal length
    if len(f_row) != len(r_row):
        print(f"Mismatched: unequal length in {line_ID} line "
              f"({f_row}, {r_row})")
        return False

    # Case of equal length
    case = 0
    identifier = 0
    for letter in range(len(f_row)):
        if f_row[letter] == '1' and r_row[letter] == '2':
            identifier = letter
            case += 1
            continue
        # Other mismatches
        if f_row[letter] != r_row[letter]:
            print(f"Mismatched: unequal letters in {line_ID} line "
                  f"({f_row}, {r_row})")
            return False

    # Valid ID
    if case == 1:
        return True
    if case != 1:
        print(f"Mismatched: 1-2 occurs {case} times in {line_ID} line "
              f"({f_row}, {r_row})")
        return False


def sample_num(read_seq, tolerance=1):
    """ Identify the number of sample by barcode sequence and bridge sequence.

    Exact match for barcode and fuzzy match for bridge which allows at most one
    different base in the same position.

    :param read_seq: (str)
        sequence in the second line of read
    :param tolerance: (int, default=1)
        Maximum fault tolerance of base in the specific sequence

    :return:
        number of sample or None
    """
    # Count for different base
    bridge_seq = read_seq[8: 8 + len(BRIDGE)]
    difference = 0
    for base in range(len(BRIDGE)):
        if bridge_seq[base] != BRIDGE[base]:
            difference += 1
    # Identify special barcode within the allowed tolerance
    bc_seq = read_seq[: 8]
    if bc_seq in barcode and difference <= tolerance:
        return barcode.index(bc_seq)
    return None


def fq_assign(num, tag, f_read, r_read) -> None:
    """ Split genome sequence and output 3 fastq files.

    Assign umi and genome sequence from tagged read and another clean genome
    sequence to three different fastq files in two fixed directories.

    :param num: (int)
        number of current sample
    :param tag: (int)
        the value of tag is only 1 or 2, which represents the one of Pair-End
        sequencing file with the equal file tag containing a special barcode,
        bridge and UMI sequence.
    :param f_read: (List)
        contain 4 entries corresponding to the 4 lines of current forward read
    :param r_read: (List)
        contain 4 entries corresponding to the 4 lines of current reverse read

    :return: None
    """
    # Check the equality of sequence and quality length
    if len(f_read[1]) != len(f_read[3]) or len(r_read[1]) != len(r_read[3]):
        print(f"Length Error: sequence in {line - 3} line and "
              f"quality in {line - 1} line are not equal")
        return None

    # Recognize tagged read by tag state
    tag_read = [f_read, r_read][tag - 1]
    cln_read = [f_read, r_read][2 - tag]
    # Write the UMI sequence and clean sequence from tagged read to fastq
    with open(os.path.join(fqDir, str(num) + f"_{tag}.fq"), 'a') as fq_tag:
        fq_tag.write('\n'.join([
            tag_read[0], tag_read[1][8 + len(BRIDGE) + 6:],
            tag_read[2], tag_read[3][8 + len(BRIDGE) + 6:]
        ]) + '\n')
    with open(os.path.join(umiDir, str(num) + ".fq"), 'a') as fq_umi:
        fq_umi.write('\n'.join([
            tag_read[0], tag_read[1][8 + len(BRIDGE): 8 + len(BRIDGE) + 6],
            tag_read[2], tag_read[3][8 + len(BRIDGE): 8 + len(BRIDGE) + 6]
        ]) + '\n')
    # Write another clean read to fastq
    with open(os.path.join(fqDir, str(num) + f"_{3 - tag}.fq"), 'a') as fq_cln:
        fq_cln.write('\n'.join(cln_read) + '\n')
    return None


barcodeFile = open(sys.argv[1])
barcode = [bc.strip() for bc in barcodeFile]
barcodeFile.close()

fFq = open(sys.argv[2])
rFq = open(sys.argv[3])
# Output
fqDir = sys.argv[4]
umiDir = sys.argv[5]

# Init
fRead = []
rRead = []
sample = 0
line = 0
cache = 0
emptySpl = list(range(len(barcode)))
BRIDGE = "ATAGCGACGCGTTTCAAC"

# For case count
bothMatched = 0
onceMatched = 0
noneMatched = 0

while True:
    try:
        # Iterate by row and touch exception
        fRow = next(fFq).strip()
        rRow = next(rFq).strip()
        line += 1

        if (line + 3) % 4 == 0:
            # Check if reads ID are valid
            if not valid_tag(fRow, rRow, line):
                continue

            # Receive active cache state and turn off it after file writen
            if cache > 0:
                # Fastq output
                fq_assign(sample, cache, fRead, rRead)

                # Limited cache variables
                cache = 0
            fRead = []
            rRead = []

        # Sample special sequence identify
        if (line + 2) % 4 == 0:
            fCase = sample_num(fRow)
            rCase = sample_num(rRow)

            if fCase is not None and rCase is not None:
                bothMatched += 1
                continue
            if fCase is None and rCase is None:
                noneMatched += 1
                continue

            # Single-ended matching
            onceMatched += 1
            if fCase is not None and rCase is None:
                sample = fCase
                cache = 1
            if fCase is None and rCase is not None:
                sample = rCase
                cache = 2
            emptySpl[sample] = -1

        fRead.append(fRow)
        rRead.append(rRow)

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
    emptyUMI = open(os.path.join(umiDir, str(i) + ".fq"), 'a')
    emptyUMI.close()

# Terminal output
reads = bothMatched + onceMatched + noneMatched
print(
    f"Double-ended match: {bothMatched} ({100 * bothMatched / reads:.2f}%), "
    f"Single-ended match: {onceMatched} ({100 * onceMatched / reads:.2f}%), "
    f"None-ended match: {noneMatched} ({100 * noneMatched / reads:.2f}%)"
)
