"""
Created on October 15, 2023

@Author
Han Rui (hanry.rui@gmail.com)

@Usage
python3 categorize_barcode.py barcode.txt -i1 R1.fq -i2 R2.fq -out outputDir

@Function
Categorize NGS data by barcode sequence.
"""
import os
import argparse


def command():
    parser = argparse.ArgumentParser(
        description="Categorize the NGS data for each sample "
                    "according to barcode sequences.",
        # Clean default groups
        add_help=False
    )
    # Add to positional arguments group by default
    parser.add_argument("barcode", help="barcode list")
    # Create a required group behind the positional arguments group
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('-i1', help="fastq_1", required=True)
    required_group.add_argument('-i2', help="fastq_2", required=True)
    required_group.add_argument(
        '-out', help="output directory", required=True
    )
    # Proactively declare the position of optional group
    optional_group = parser.add_argument_group('optional arguments')
    optional_group.add_argument(
        "-h", "--help", action='help', help="show this help message and exit"
    )
    return parser.parse_args()


class Fastq:
    def __init__(self, barcode_file: str):
        self.read_len = 4
        self.barcode_len = 8
        with open(barcode_file) as barcode_obj:
            self.barcode = [b.strip() for b in barcode_obj]

        # Exception capture
        if len(set(self.barcode)) != len(self.barcode):
            print("Error: Redundant sequence in barcodes.")

    def read(self, fq):
        read_unit = []
        fq_obj = open(fq)
        for line in fq_obj:
            read_unit.append(line)
            if len(read_unit) == self.read_len:
                yield read_unit
                read_unit = []
        fq_obj.close()

    def index(self, seq_1: str, seq_2: str) -> int:
        marker_1 = seq_1[: self.barcode_len]
        marker_2 = seq_2[: self.barcode_len]
        if marker_1 in self.barcode and marker_2 not in self.barcode:
            return self.barcode.index(marker_1)
        if marker_1 not in self.barcode and marker_2 in self.barcode:
            return self.barcode.index(marker_2)


if __name__ == '__main__':
    cmd = command()
    fastq = Fastq(cmd.barcode)
    valid = []
    # Iterator of read
    read_1 = fastq.read(cmd.i1)
    read_2 = fastq.read(cmd.i2)
    while True:
        try:
            read_arr_1 = next(read_1)
            read_arr_2 = next(read_2)
            barcode_idx = fastq.index(read_arr_1[1], read_arr_2[1])
            if barcode_idx is None:
                continue
            valid.append(barcode_idx)
            # Write reads
            output_1 = os.path.join(cmd.out, f"{barcode_idx}_1.fq")
            output_2 = os.path.join(cmd.out, f"{barcode_idx}_2.fq")
            with open(output_1, 'a') as out_1:
                out_1.write(''.join(read_arr_1))
            with open(output_2, 'a') as out_2:
                out_2.write(''.join(read_arr_2))
        except StopIteration:
            break
    # Touch empty files
    for item in range(len(fastq.barcode)):
        if item not in valid:
            open(os.path.join(cmd.out, f"{item}_1.fq"), 'a').close()
            open(os.path.join(cmd.out, f"{item}_2.fq"), 'a').close()
