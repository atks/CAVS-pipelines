#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2022 Adrian Tan <adrian_tan@nparks.gov.sg>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permsit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import click
import re
from queue import PriorityQueue

@click.command()
@click.option('-s', '--sample_file', required=True, help='sample file')
@click.option('-o', '--output_fasta_file', required=True, help='output contig fasta file')
def main(sample_file, output_fasta_file):
    """
    Rename contigs to specify sample ID and kraken2 annotation

    e.g. rename_contig.py -s vibrio.sa -i contigs.fa -o renamed.contigs.fasta
   """
    print(f"sample file: {sample_file}")
    print(f"output FASTA file: {output_fasta_file}")

    with open(output_fasta_file, 'w') as out_file:
        with open(sample_file, 'r') as sample_file:
            for line in sample_file:
                #print(line)
                if not line.startswith('#'):
                    sample_id, annotation, read1_fastq, read2_fastq, assembled_contigs_file = line.rstrip().split('\t')
                    print(assembled_contigs_file)
                    no = 0
                    with open(assembled_contigs_file, 'r') as input_fasta_file:
                        block = []
                        for fasta_line in input_fasta_file:
                            if fasta_line.startswith('>'):
                                if block:
                                    out_file.write(''.join(block) + '\n')
                                    block = []
                                no += 1
                                new_header =  f'>{sample_id}_contig_{no} {annotation} {fasta_line.rstrip().lstrip(">")}\n'
                                #print(new_header)
                                out_file.write(new_header)
                            else:
                                block.append(fasta_line.strip())
                        #last record
                        if block:
                            out_file.write(''.join(block) + '\n')
                            block = []

if __name__ == '__main__':
    main()
