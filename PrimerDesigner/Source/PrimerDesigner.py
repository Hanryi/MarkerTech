"""
Updated on December 18, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 PrimerDesigner.py <i_int> <i_fa> <o_str> <i_num>(9)
# <i_int>: strategy num (BP-Seq: 2; iBP-Seq: 3)
# <i_fa>: readable text file
# <o_str>: path of xlsx with suffix
# <i_num>(9): 9 nums get from default settings

@Function
Design specific primers in batches for BP-Seq and iBP-Seq
"""

import sys

import primer3
import pandas as pd


def primer_attr_matrix(temp_name, temp_seq, keep_area, len_arr, melt_arr):
    """ Design primers and return main attr list from json format.

    :param temp_name: Str
        name of the template sequence
    :param temp_seq: Str
        template sequence for designing primers
    :param keep_area: List
        interval endpoints in temp_seq (len == 2)
    :param len_arr: List
        length args of primers and PCR product
    :param melt_arr: List
        Tm args of primers

    :return: List
        main attrs of one pair of primers
    """
    # Constant length value (iBP: 118 = 150 - 32, BP: 36 = 150 - 114)
    left_primer_area_length = 150
    right_primer_area_length = 118
    if strategy == 2:
        right_primer_area_length = 36

    # Select start index
    left_primer_start = max(
        keep_area[1] - left_primer_area_length,
        0,
    )
    right_primer_start = min(
        keep_area[0] - 2 + right_primer_area_length,
        len(tempSeq) - 1,
    )

    # Parameter configuration
    seq_args = {
        # Identify the source of the chosen primers (str)
        'SEQUENCE_ID': temp_name,
        # Target sequence (5'->3') and it cannot span several lines (str)
        'SEQUENCE_TEMPLATE': temp_seq,
        # Subset of template used to design primers ([int]; default: all)
        'SEQUENCE_INCLUDED_REGION': [0, len(temp_seq)],
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

        # ---Primer Length
        # Minimum length of primers (int; default: 18)
        'PRIMER_MIN_SIZE': len_arr[0],
        # Optimum length of primers (int; default: 20)
        'PRIMER_OPT_SIZE': len_arr[1],
        # Maximum length of primers (int; default: 27)
        'PRIMER_MAX_SIZE': len_arr[2],

        # ---Production Length (int, default: 0)
        # Optimum length of production (int; default: 0)
        'PRIMER_PRODUCT_OPT_SIZE': len_arr[4],
        # Production size range list ([int]; default: 0)
        'PRIMER_PRODUCT_SIZE_RANGE': [
            [len_arr[3], len_arr[5]],
        ],

        # ---Tm (Celsius)
        # Minimum temperature of primers (float; default: 57.0)
        'PRIMER_MIN_TM': melt_arr[0],
        # Optimum temperature of primers (float; default: 60.0)
        'PRIMER_OPT_TM': melt_arr[1],
        # Maximum temperature of primers (float; default: 63.0)
        'PRIMER_MAX_TM': melt_arr[2],
        # Maximum difference between pair of Tm (float; default 100.0)
        'PRIMER_PAIR_MAX_DIFF_TM': 10,

        # Print secondary structures (boolean; default: 0)
        'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT': 0,

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
    # Design primer
    primer_attr = primer3.bindings.designPrimers(seq_args, global_args)

    # Organize main attr
    main_attr = [
        temp_name,
        primer_attr['PRIMER_LEFT_0_SEQUENCE'],
        primer_attr['PRIMER_RIGHT_0_SEQUENCE'],
        primer_attr['PRIMER_LEFT_0'][0],
        primer_attr['PRIMER_RIGHT_0'][0],
        primer_attr['PRIMER_LEFT_0'][1],
        primer_attr['PRIMER_RIGHT_0'][1],
        round(primer_attr['PRIMER_LEFT_0_TM'], 2),
        round(primer_attr['PRIMER_RIGHT_0_TM'], 2),
        round(primer_attr['PRIMER_LEFT_0_GC_PERCENT'], 2),
        round(primer_attr['PRIMER_RIGHT_0_GC_PERCENT'], 2),
        primer_attr['PRIMER_PAIR_0_PRODUCT_SIZE'],
    ]
    # BRIDGE inherits from global variable
    main_attr[2] = BRIDGE + main_attr[2]

    return main_attr


# Template file
tempTxt = sys.argv[2]
attrTab = sys.argv[3]
# Args
strategy = eval(sys.argv[1])
lenArr = list(map(int, sys.argv[4: 10]))
meltArr = list(map(float, sys.argv[10:]))

# Add header of primer result
primerAttrMat = []
attrName = [
    'Name',
    'Primer (L)', 'Primer (R)',
    "5'pos\n(L)", "5'pos\n(R)",
    'Length\n(L)', 'Length\n(R)',
    'Tm (L)', 'Tm (R)',
    'GC (L)', 'GC (R)',
    'Product\nSize',
]
primerAttrMat.append(attrName)

# Different variable
BRIDGE = "ATAGCGACGCGTTTCAAC"
sheetName = 'iBP-Primers'
if strategy == 3:
    for product in range(3, 6):
        lenArr[product] -= 32
    BRIDGE += "NNNNNN"
if strategy == 2:
    for product in range(3, 6):
        lenArr[product] -= 114
    sheetName = 'BP-Primers'

# Init
tempName = ''
keepArea = []
tempSeq = ''

tempFile = open(tempTxt)
for line in tempFile:
    # Attach sequence based on valid header
    if line[0] != '>':
        if tempName != '' and keepArea != [] and keepArea[0] <= keepArea[1]:
            tempSeq += line.strip()
        continue

    if len(tempSeq) >= 150 and len(tempSeq) >= keepArea[1]:
        # Design primer for previous valid record
        mainAttr = primer_attr_matrix(
            tempName, tempSeq, keepArea, lenArr, meltArr
        )
        primerAttrMat.append(mainAttr)

    # Entry point
    tempName = ''
    keepArea = []
    tempSeq = ''

    # Description line
    headerVec = line[1:].strip().split('|')
    if len(headerVec) >= 2:
        # Sequence name
        tempName = headerVec[0].strip()

        # Area annotation
        note = headerVec[1].replace(' ', '')
        if note == '' or note.count('-') > 1 or '-' in [note[0], note[-1]]:
            continue
        keepArea = [int(i) for i in note.split('-')]

        # Expand the single base position to a interval
        if len(keepArea) == 1:
            keepArea.append(keepArea[0])
tempFile.close()

# Pick primer for the last sequence
if len(tempSeq) >= 150 and len(tempSeq) >= keepArea[1]:
    mainAttr = primer_attr_matrix(
        tempName, tempSeq, keepArea, lenArr, meltArr
    )
    primerAttrMat.append(mainAttr)

# Export dataframe as excel
primerAttrFrame = pd.DataFrame(primerAttrMat)
frameWriter = pd.ExcelWriter(attrTab, engine='xlsxwriter')
primerAttrFrame.to_excel(
    frameWriter, sheet_name=sheetName, header=False, index=False
)

# # Print all column
# pd.set_option('display.max_columns', None)
# print(primerAttrFrame)

""" Customize style.
"""
primerBook = frameWriter.book
primerSheet = frameWriter.sheets[sheetName]
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
# Most chars in each column
charLen = [0 for _ in range(6)]
aimCol = [0, 1, 2, 3, 4, 11]
for item in primerAttrMat[1:]:
    for col in aimCol:
        if len(str(item[col])) + 2 > charLen[aimCol.index(col)]:
            charLen[aimCol.index(col)] = len(str(item[col])) + 2
# Set column width (pixel)
primerSheet.set_column('A:A', max(8.64, charLen[0]))
primerSheet.set_column('B:B', max(21.64, charLen[1]))
primerSheet.set_column('C:C', max(21.64, charLen[2]))
primerSheet.set_column('D:E', max(7, charLen[3], charLen[4]))
primerSheet.set_column('F:K', 7.64)
primerSheet.set_column('L:L', max(9, charLen[5]))
# Save as file
frameWriter.save()
