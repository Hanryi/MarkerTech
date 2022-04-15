"""
Created on April 14, 2022

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 MethReference.py <bc.txt> <lines.fa>
# <bc.txt>: a barcode file, contains two columns without header (barcode,
#           refType) and the separator of columns should be ','
#           1. barcode: 8bp barcode sequence
#           2. refType: combination of inbred line name (sep='/')
# <lines.fa>: a reference sequence file

@Function
1. Create a directory for each line to save the divided reference file
2. Divide the input reference file into several files according lines
   combination in the barcode file
"""

import os
import sys


class Reference:
    """ Make sub-directory for each inbred line combination, each sub-directory
    saves the reference sequence of all the inbred line in combination.
    Moreover, the class attribute saves the {inbred line: segment number} and
    {inbred line: inbred line sequences} as detail information.
    """
    line_set = {}
    line_seq = {}

    def __init__(self, bc_file: str, ref_file: str):
        self.bcFile = bc_file
        self.refFile = ref_file
        self.refDir = os.path.abspath(os.path.dirname(ref_file))

    def inb_line_type(self):
        """ Filter for unique inbred line combination types.

        :return: (list of list)
            each unique combination type save as a sub-list
        """
        line_type = [t.strip().split(',')[1] for t in open(self.bcFile)]
        type_unit = []
        for i in line_type:
            if sorted(i.split('/')) not in type_unit:
                type_unit.append(sorted(i.split('/')))
        return type_unit

    def sub_seq_dir(self, lineage: str) -> str:
        """ Make directory for each inbred line in the working directory.

        :param lineage: (str)
            name of inbred line

        :return: (str)
            path of sub-directory
        """
        sub_dir = os.path.join(self.refDir, lineage)
        if os.path.exists(sub_dir):
            return sub_dir

        os.mkdir(sub_dir)
        return sub_dir

    def inb_line_library(self):
        """ Build a reference sequence library as combination unit.

        :return: None
        """
        lineage = None
        for i in open(self.refFile):
            if i[0] == ">":
                line_tuple = i[1:].strip().split('-')
                lineage = line_tuple[1]
                # Name uniqueness check (non-functional test)
                if lineage not in self.line_set.keys():
                    self.line_set[lineage] = []
                if line_tuple[0] in self.line_set[lineage]:
                    print(f"DescriptionError: sequence name conflict")
                    break
                self.line_set[lineage].append(line_tuple[0])

                # Build reference sequence library
                if lineage not in self.line_seq.keys():
                    self.line_seq[lineage] = []
            self.line_seq[lineage].append(i)
        return None

    def inb_line_assign(self):
        """ Generate corresponding types of reference sequence files.

        :return: None
        """
        type_unit = self.inb_line_type()
        self.inb_line_library()

        # Assign reference sequence on demand
        for type_comb in type_unit:
            sub_dir = self.sub_seq_dir('_'.join(type_comb))
            for type_ele in type_comb:
                seq_append(
                    os.path.join(sub_dir, '_'.join(type_comb) + '.fa'),
                    "".join(self.line_seq[type_ele])
                )
        return None


def seq_append(file_path: str, text):
    """ Write to file as append.

    :param file_path: (str)
        path of file
    :param text: (str)
        file content

    :return: None
    """
    sub_ref = open(file_path, 'a')
    sub_ref.write(text)
    sub_ref.close()


refObj = Reference(sys.argv[1], sys.argv[2])
refObj.inb_line_assign()
