"""
Created on February 7, 2022

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 MethFraction.py <umi.fq> <m.bam> <m.tab> <o_m.csv>
# <umi.fq>: fastq file of fastq
# <m.bam>: a binary file after <m_1/2.fq> aligned to reference genome
# <m.tab>: a table file containing methylation info (e.g. chr, pos, context)
#          1. chr: chromosome name
#          2. pos: methylation position
#          3. context: methylation type
# <o_m.csv>: a table file (Chr, Pos, Context, Active)
#            1. Chr: chromosome name
#            2. Pos: methylation position
#            3. Context: methylation type
#            4. Active: active fraction

@Function
Calculate active fraction for each methylation site after UMI deduplication

@Note
Active = C / (C + T)
"""

import sys

import pysam


def generate_triple(tab_path):
    """ Generate main methylation info as triple.

    The main methylation info containing chr, pos and context.

    :param tab_path: (Str)
        file path of a readable table

    :return: (Generator)
        each time the generated triple format all are str type
    """
    meth_tab = open(tab_path)
    next(meth_tab)
    for line in meth_tab:
        meth_record = line.split("\t", 4)
        yield [meth_record[0], meth_record[1], meth_record[3]]
    meth_tab.close()


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


# Input args
seqRoster = umi_roster(sys.argv[1])
bamObj = pysam.AlignmentFile(sys.argv[2], 'rb')
if not bamObj.check_index():
    # Abnormal exit
    raise SystemExit(1)
triple = generate_triple(sys.argv[3])

# Output arg
separator = ","
header = ["Chr", "Pos", "Context", "Active"]
with open(sys.argv[4], 'w') as methActiveTab:
    methActiveTab.write(separator.join(header) + "\n")

while True:
    try:
        # Iterate methylation by generator
        methTriple = next(triple)
        position = int(methTriple[1])
        segPile = bamObj.fetch(methTriple[0], position - 1, position)

        umiPanel = [[], []]
        for seg in segPile:
            # mBP-Seq type
            mutSeq = mut_locator(
                seg.query_sequence, position - 1, 1, seg.get_aligned_pairs()
            )
            if mutSeq != "C" and mutSeq != "T" or mutSeq is None:
                continue
            mutState = 0
            if mutSeq == "T":
                mutState = 1

            # UMI deduplication
            umiSeq = seqRoster[seg.query_name]
            if umiSeq not in sum(umiPanel, []):
                umiPanel[mutState].append(umiSeq)

        # Calculate active
        umiPool = sum(umiPanel, [])
        if len(umiPool) == 0:
            continue
        active = len(umiPanel[0]) / len(umiPool)

        # Real-time output
        methTriple.append(str(active))
        with open(sys.argv[4], 'a') as methActiveTab:
            methActiveTab.write(separator.join(methTriple) + "\n")
    except StopIteration:
        break
