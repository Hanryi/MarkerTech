"""
Created on August 28, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 SeqPurifier.py [*BP] <bc.txt> <_1.fq> <_2.fq> <o_fq.d> <o_umi.d>
# [*BP]: a string option (BP or mBP or iBP)
# <bc.txt>: file containing 8bp barcodes
# <_1.fq>: clean data of forward sequencing
# <_2.fq>: clean data of reverse sequencing
# <o_fq.d>: directory for storing pair of fq which removed additional seq
# <o_umi.d>: (optional)
#            directory for storing fq which contains only UMI fragments

@Function
1. Divide mixed sequencing data to several pairs of fq by barcode and bridge
2. Separate UMI sequence in reads from genome sequence
"""

import os
import sys
from typing import Optional


class Read:
    """ Sequencing Read class.

    This class has a variable attribute named BRIDGE and has a class method
    named sample_num(). After instantiate, each object has an init object
    attribute named idLocation which can bind a value by object method named
    and be called by the object method named fq_reorg().
    """
    plan = sys.argv[1]

    BRIDGE = None
    if plan == "BP" or "iBP":
        BRIDGE = "ATAGCGACGCGTTTCAAC"
    if plan == "mBP":
        BRIDGE = "ATAGCGACG"

    YFP = "CAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAG"
    EndR = "AAACGCCAGCAATAATTACG"

    BC_LEN = 8
    UMI_LEN = 6

    def __init__(self):
        self.idLocation = None
        self.labelEnd = 0

    def valid_tag(self, f_id, r_id):
        """ Judge the validity of current ID pair.

        Valid ID should be judged by identifier which only allowed occur once as
        "1"-"2" in F/R ID respectively.

        :param f_id: (str)
            forward reads ID
        :param r_id: (str)
            reverse reads ID

        :return: (bool)
        """
        # Case of unequal length
        if len(f_id) != len(r_id):
            print(f"Mismatched: unequal length in {line} line")
            return False

        # Case of equal length
        case = 0
        identifier = 0
        for letter in range(len(f_id)):
            if f_id[letter] == '1' and r_id[letter] == '2':
                identifier = letter
                case += 1
                continue
            # Other mismatches
            if f_id[letter] != r_id[letter]:
                print(f"Mismatched: unequal letters in {line} line")
                return False

        # Valid ID
        if case == 1:
            # Bind the identifier to current class attribute
            self.idLocation = identifier
            return True
        if case != 1:
            print(f"Mismatched: 1-2 occurs {case} times in {line} line")
            return False

    def fq_reorg(self, num: int, tag: int, f_read: list, r_read: list) -> None:
        """ Reorganize genome sequence and output 3 fastq files.

        Assign umi and genome sequence from tagged read and another clean
        genome sequence to three different fastq files in two directories.

        :param num: (int)
            number of current sample
        :param tag: (int)
            the value of tag is only 1 or 2, which represents the one of
            Pair-End sequencing file with the equal file tag containing a
            special barcode, bridge and UMI sequence.
        :param f_read: (list)
            contain 4 entries which are the 4 lines of forward read
        :param r_read: (list)
            contain 4 entries which are the 4 lines of reverse read

        :return: None
        """
        # Check the equality of sequence and quality length
        if not self.equal_seq(f_read) or not self.equal_seq(r_read):
            return None

        # Recognize tagged read by tag state
        tag_read = [f_read, r_read][tag - 1]
        cln_read = [f_read, r_read][2 - tag]
        # Write the clean read to fastq
        with open(fastq_file(fqDir, num, f"_{3 - tag}.fq"), 'a') as fq_cln:
            fq_cln.write('\n'.join(cln_read) + '\n')
        if self.plan == "BP":
            self.labelEnd -= self.UMI_LEN
        # Write the clean sequence from tagged read to fastq
        with open(fastq_file(fqDir, num, f"_{tag}.fq"), 'a') as fq_tag:
            fq_tag.write('\n'.join([
                tag_read[0], tag_read[1][self.labelEnd:],
                tag_read[2], tag_read[3][self.labelEnd:]
            ]) + '\n')
        if self.plan == "BP":
            return None

        # Write the UMI sequence to fastq
        with open(fastq_file(umiDir, num, ".fq"), 'a') as fq_umi:
            # Replace read ID with name
            fq_umi.write('\n'.join([
                tag_read[0][1: self.idLocation + 1],
                tag_read[1][self.labelEnd - self.UMI_LEN: self.labelEnd],
                tag_read[2],
                tag_read[3][self.labelEnd - self.UMI_LEN: self.labelEnd]
            ]) + '\n')
        return None

    def sample_num(self, read_seq: str) -> Optional[int]:
        """ Identify the number of sample by barcode, bridge and UMI sequence.

        Exact match for barcode and fuzzy match for bridge. Moreover, "N"
        (unknown base in sequencing) is not allowed in UMI sequence.

        :param read_seq: (str)
            sequence in the second line of read

        :return:
            number of sample or None
        """
        if read_seq.count(self.BRIDGE) != 1:
            return None
        barcode_len = 0
        if self.plan == "BP":
            barcode_len = self.BC_LEN + len(self.YFP) + self.BC_LEN

        flank_seq = read_seq.split(self.BRIDGE)
        self.labelEnd += len(flank_seq[0]) + len(self.BRIDGE) + self.UMI_LEN
        # Barcode check
        if len(flank_seq[0]) < barcode_len:
            return None
        bc_seq = flank_seq[0][-barcode_len:]

        # BP byroad
        if self.plan == "BP":
            if bc_seq in barcode:
                return barcode.index(bc_seq)
            if bc_seq not in barcode:
                return None

        # UMI check
        if len(flank_seq[1]) < self.UMI_LEN:
            return None
        umi_seq = flank_seq[1][: self.UMI_LEN]
        if bc_seq in barcode and umi_seq.count("N") <= 0:
            return barcode.index(bc_seq)
        return None

    @staticmethod
    def equal_seq(single_read):
        """ Check the equality of sequence and quality length.

        :param single_read: (list)
            contain 4 entries corresponding to the 4 lines of input read

        :return: Boolean
        """
        if len(single_read[1]) == len(single_read[3]):
            return True
        print(f"Error: quality does not correspond sequence {single_read[0]}")
        return False


def fastq_file(fq_dir, num_prefix, suffix):
    """ Build fastq filepath by input parameters.

    :param fq_dir: (str)
        a directory of fastq path
    :param num_prefix: (int)
        a sample number
    :param suffix: (str)
        pair-end file tag such as _1 or _2 and file format suffix, such as .fq

    :return: (str)
        a relative path of fastq file with full filename
    """
    return os.path.join(fq_dir, str(num_prefix) + suffix)


def jump_iter(f_iter, r_iter, count, step=2):
    """ Iterate over objects recursively.

    :param f_iter: (object)
        a forward iterable object
    :param r_iter: (object)
        a reverse iterable object
    :param count: (int)
        current line number
    :param step: (int, default=2)
        number of iterations

    :return: (object)
        recursive final two iterators and one changed count number
    """
    next(f_iter)
    next(r_iter)

    count += 1
    # Step by step
    step -= 1
    if step == 0:
        return f_iter, r_iter, count
    return jump_iter(f_iter, r_iter, count, step)


barcodeFile = open(sys.argv[2])
barcode = [bc.strip() for bc in barcodeFile]
barcodeFile.close()
if sys.argv[1] == "mBP":
    barcode = [bc.split(',')[0] for bc in barcode]

fFq = open(sys.argv[3])
rFq = open(sys.argv[4])
# Output
fqDir = sys.argv[5]
umiDir = None
if len(sys.argv) == 7:
    umiDir = sys.argv[6]

# Init
fRead = []
rRead = []
sample = 0
line = 0
cache = 0
emptySpl = list(range(len(barcode)))
thisRead = None

# For case count
bothMatched = 0
onceMatched = 0
noneMatched = 0

while True:
    try:
        # Iterate by row until touch exception
        fRow = next(fFq).strip()
        rRow = next(rFq).strip()
        line += 1

        if (line + 3) % 4 == 0:
            # Limited cache variables
            cache = 0
            fRead = []
            rRead = []

            thisRead = Read()
            # Check ID validity and bind identifier to the global object
            if not thisRead.valid_tag(fRow, rRow):
                continue

        # Sample special sequence identify
        if (line + 2) % 4 == 0 and thisRead.idLocation is not None:
            fCase = thisRead.sample_num(fRow)
            rCase = thisRead.sample_num(rRow)

            if fCase is not None and rCase is not None:
                bothMatched += 1
                fFq, rFq, line = jump_iter(fFq, rFq, line)
                continue
            if fCase is None and rCase is None:
                noneMatched += 1
                fFq, rFq, line = jump_iter(fFq, rFq, line)
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

        # Fastq output
        if line % 4 == 0 and cache > 0:
            thisRead.fq_reorg(sample, cache, fRead, rRead)

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

    # Output UMI sequence for iBP and mBP
    if umiDir is None:
        break
    emptyUMI = open(os.path.join(umiDir, str(i) + ".fq"), 'a')
    emptyUMI.close()

# Terminal output
reads = bothMatched + onceMatched + noneMatched
print(
    f"Double-ended match: {bothMatched} ({100 * bothMatched / reads:.2f}%), "
    f"Single-ended match: {onceMatched} ({100 * onceMatched / reads:.2f}%), "
    f"None-ended match: {noneMatched} ({100 * noneMatched / reads:.2f}%)"
)
