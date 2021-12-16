"""
Updated on December 09, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3.x PrimerDesigner.py <i_int> <i_text> <i_str_1> <i_str_2> <<i_num_1>...<i_num_8>>
# <i_int>: strategy num
# <i_text>: readable text file. Seq and its SNP index interspersed by '\n'
# <i_str_1>: user account path
# <i_str_2>: filename of tag txt and xlsx is composed by username and time
# <<i_num_1>...<i_num_8>>: 8 num get from default settings

@Function
Containing 2 primer design strategies.
    Each strategy will return the best choice from design result according
    a series of params given by command line
"""

import os
import sys

import primer3
import pandas as pd


def primer_designer(template_name, template_seq, keep_area, len_arr, melt_arr):
    """
    Designer for primers

    :param template_name: Str
        name of the template sequence
    :param template_seq: Str
        template sequence for designing primers
    :param keep_area: Tuple
        pair of index of snp in target_seq
    :param len_arr: List
        containing lens of primers and PCR product
    :param melt_arr: List
        containing Tm value of primers

    :return: Dict
        containing all details of a best primer design strategy
    """
    # Constant length value (118 = 150 - 32)
    left_primer_area_length = 150
    right_primer_area_length = 118

    # Select start index
    left_primer_start = max(
        keep_area[0] - left_primer_area_length,
        0,
    )
    # right_primer_end = min()
    right_primer_start = min(
        keep_area[1] - 2 + right_primer_area_length,
        len(tempSeq) - 1,
    )
    # Parameter configuration
    seq_args = {
        # Identify the source of the chosen primers (str)
        'SEQUENCE_ID': template_name,
        # Target sequence (5'->3') and it cannot span several lines (str)
        'SEQUENCE_TEMPLATE': template_seq,
        # Subset of template used to design primers ([int]; default: all)
        'SEQUENCE_INCLUDED_REGION': [0, len(template_seq)],
        # Internals that primer pairs must flank ([int]; default: None)
        'SEQUENCE_TARGET': [keep_area[0] - 1, keep_area[1] - keep_area[0]],
        # Primers interval which -1 means no limit ([[int],...]; default: None)
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [
            [left_primer_start, left_primer_area_length, -1, -1],
            [-1, -1, right_primer_start, right_primer_area_length]
        ],
    }
    global_args = {
        # ---Design Strategy
        # Num of paired primers to return sorted by "quality" (int; default: 5)
        'PRIMER_NUM_RETURN': 1,
        # Pick an internal primer (bool; default: 0)
        'PRIMER_PICK_INTERNAL_OLIGO': 0,

        # ---Length
        # Optimum length of primers (int; default: 20)
        'PRIMER_OPT_SIZE': len_arr[0],
        # Minimum length of primers (int; default: 18)
        'PRIMER_MIN_SIZE': len_arr[1],
        # Maximum length of primers (int; default: 27)
        'PRIMER_MAX_SIZE': len_arr[2],
        # Product size range list ([int]; default: 0)
        'PRIMER_PRODUCT_SIZE_RANGE': [len_arr[3], len_arr[4]],

        # ---Tm (Celsius)
        # Optimum temperature of primers (float; default: 60.0)
        'PRIMER_OPT_TM': melt_arr[0],
        # Minimum temperature of primers (float; default: 57.0)
        'PRIMER_MIN_TM': melt_arr[1],
        # Maximum temperature of primers (float; default: 63.0)
        'PRIMER_MAX_TM': melt_arr[2],
        # Maximum difference between pair of Tm (float; default 100.0)
        'PRIMER_PAIR_MAX_DIFF_TM': 10,

        # Print secondary structures (boolean; default: 0)
        'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT': 1,

        # ---GC Percent
        # Minimum GC percent of primers (float; default: 20.0)
        'PRIMER_MIN_GC': 20.0,
        # Optimum GC percent of primers (float; default: 50.0)
        'PRIMER_OPT_GC_PERCENT': 50.0,
        # Maximum GC percent of primers (float; default: 80.0)
        'PRIMER_MAX_GC': 80.0,

        # ---Base Feature
        # Maximum repeat length (AAAAAA) of a mononucleotide  (int; default: 5)
        'PRIMER_MAX_POLY_X': 5,
        # Maximum number of N in any primer (int; default: 0)
        'PRIMER_MAX_NS_ACCEPTED': 0,
    }

    return primer3.bindings.designPrimers(seq_args, global_args)


# Constant string
BRIDGE = "ATAGCGACGCGTTTCAAC"

# Args
# strategy = eval(sys.argv[1])
# filename = os.path.join(sys.argv[3], sys.argv[4])
# lenArr = list(map(int, sys.argv[5: 10]))
# meltArr = list(map(float, sys.argv[10: 13]))
strategy = 3
filename = 'testSeq(1)'
lenArr = [20, 18, 25, 280, 320]
meltArr = [60, 50, 65]

# Template file
# tempTxt = sys.argv[2]
tempTxt = "testSeq(1).txt"
primerAttrMat = []
attrName = [
    'Name',
    'Primer (L)', 'Primer (R)',
    '5`pos\n(L)', '5`pos\n(R)',
    'Length\n(L)', 'Length\n(R)',
    'Tm (L)', 'Tm (R)',
    'GC (L)', 'GC (R)',
    'Product\nSize',
]
primerAttrMat.append(attrName)

with open(tempTxt) as tempFile:
    # Init
    tempName = ''
    tempSeq = ''
    keepArea = []

    for line in tempFile:

        if line[0] == '>':
            tempName = line[1:].strip()
            tempSeq = ''
            continue

        if line[0] == '#':
            # Check for name and sequence length
            if tempName == '' or len(tempSeq) < 150:
                continue

            # Check for area annotation
            note = line[1:].strip()
            if note == '':
                continue
            if note.count('-') > 1 or note[0] == '-' or note[-1] == '-':
                continue
            keepArea = [int(i) for i in line[1:].strip().split('-')]

            # Expand single base
            if keepArea[0] > len(tempSeq):
                continue
            if len(keepArea) == 1 and keepArea[0] < len(tempSeq):
                keepArea.append(keepArea[0])

            # All len == 2
            if keepArea[0] > keepArea[1] or keepArea[1] > len(tempSeq):
                continue

            # Primer design
            # print(tempName, len(tempSeq), keepArea)
            primerAttr = primer_designer(
                tempName, tempSeq, keepArea, lenArr, meltArr
            )
            mainAttr = [
                tempName,
                primerAttr['PRIMER_LEFT_0_SEQUENCE'],
                primerAttr['PRIMER_RIGHT_0_SEQUENCE'],
                primerAttr['PRIMER_LEFT_0'][0],
                primerAttr['PRIMER_RIGHT_0'][0],
                primerAttr['PRIMER_LEFT_0'][1],
                primerAttr['PRIMER_RIGHT_0'][1],
                round(primerAttr['PRIMER_LEFT_0_TM'], 2),
                round(primerAttr['PRIMER_RIGHT_0_TM'], 2),
                round(primerAttr['PRIMER_LEFT_0_GC_PERCENT'], 2),
                round(primerAttr['PRIMER_RIGHT_0_GC_PERCENT'], 2),
                primerAttr['PRIMER_PAIR_0_PRODUCT_SIZE'],
            ]

            # The only difference
            if strategy == 3:
                mainAttr[2] = BRIDGE + "NNNNNN" + mainAttr[2]
            if strategy == 2:
                mainAttr[2] = BRIDGE + mainAttr[2]

            primerAttrMat.append(mainAttr)
            continue
            # break

        # Attach sequence string
        tempSeq += line.strip()

# Build dataframe
primerAttrFrame = pd.DataFrame(primerAttrMat)

# Export as excel
frameWriter = pd.ExcelWriter(filename + '.xlsx', engine='xlsxwriter')
primerAttrFrame.to_excel(
    frameWriter, sheet_name='Sheet1', header=False, index=False
)

# # Print all column
# pd.set_option('display.max_columns', None)
# print(primerAttrFrame)

# Customize style
primerBook = frameWriter.book
primerSheet = frameWriter.sheets['Sheet1']
headerStyle = primerBook.add_format({
    'valign': 'vcenter',
    'align': 'center',
    'color': '#2f5b66',
    # 'fg_color': '#F4B084',
    'font_name': 'Times New Roman',
    'bold':  True,
    'text_wrap': True,
})
defaultCellStyle = primerBook.add_format({
    'valign': 'vcenter',
    'align': 'center',
})
# Set row style include header row
primerSheet.set_row(0, 28, cell_format=headerStyle)
for row in range(1, primerAttrFrame.shape[0]):
    primerSheet.set_row(row, 16)
# Set column style
primerSheet.set_column('A:L', cell_format=defaultCellStyle)
# Set column width (pixel)
primerSheet.set_column('A:A', 7.64)
primerSheet.set_column('B:B', 21.55)
primerSheet.set_column('C:C', 47.45)
primerSheet.set_column('D:E', 6.09)
primerSheet.set_column('F:K', 7.18)
primerSheet.set_column('L:L', 8.27)
# Save as file
frameWriter.save()
