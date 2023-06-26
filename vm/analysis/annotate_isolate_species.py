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
@click.option('-o', '--annotated_sample_file', required=True, help='new annotated sample file')
def main(sample_file, annotated_sample_file):
    """
    Annotate isolate species from kraken annotation

    e.g. annotate_isolate_species.py -s vibrio.sa -o vibrio.annotated.sa
   """
    print(f"sample  file: {sample_file}")
    print(f"annotated sample  file: {annotated_sample_file}")

    outfile = {sample_file}

    with open(annotated_sample_file, 'w') as out_file:
        out_file.write(f'#sample-id\tannotation\tread1_fastq\tread2_fastq\tassembled_contigs\n')
        with open(sample_file, 'r') as file:
            for line in file:
                if not line.startswith('#'):
                    sample_id, read1_fasta, read2_fasta, contig_fasta = line.rstrip().split('\t')
                    run_no, sample_no, species, id = sample_id.split('_')
                    analysis_dir = f'/net/singapura/vm/hts/ilm{run_no}/analysis/{sample_no}_{species}_{id}/kraken2_result'
                    kraken_file = f'{analysis_dir}/report.txt'
                    n_s = 0
                    n = 0
                    annotation = ''
                    with open(kraken_file, 'r') as file:
                        q = PriorityQueue()
                        for line in file:
                            percent, clade_no, taxon_no, rank, tax_id, tax_name = [field.strip() for field in line.rstrip().split('\t')]
                            n += 1
                            percent = float(percent)
                            if rank=='S':
                                n_s += 1
                                # print(percent)
                                # print(clade_no)
                                # print(taxon_no)
                                # print(rank)
                                # print(tax_id)
                                # print(tax_name)
                                q.put((-percent, tax_name))
                        # print(f'n_s = {n_s}')
                        while not q.empty():
                            t = q.get()
                            species = t[1]
                            percent = -t[0]
                            if species.startswith("Vibrio") and percent>1:
                                if annotation=='':
                                    annotation += f"{species}:{percent}"
                                else:
                                    annotation += f"|{species}:{percent}"
                    out_file.write(f'{sample_id}\t{annotation}\t{read1_fasta}\t{read2_fasta}\t{contig_fasta}\n')

#  Sample Report Output Format
#  Like Kraken 1, Kraken 2 offers two formats of sample-wide results. Kraken 2's standard sample report format is tab-delimited with one line per taxon. The fields of the output, from left-to-right, are as follows:
#  1. Percentage of fragments covered by the clade rooted at this taxon
#  2. Number of fragments covered by the clade rooted at this taxon
#  3. Number of fragments assigned directly to this taxon
#  4. A rank code, indicating
#       (U)nclassified,
#       (R)oot,
#       (D)omain,
#       (K)ingdom,
#       (P)hylum,
#       (C)lass,
#       (O)rder,
#       (F)amily,
#       (G)enus, or
#       (S)pecies.
#     Taxa that are not at any of these 10 ranks have a rank code that is
#     formed by using the rank code of the closest ancestor rank with a number
#     indicating the distance from that rank. E.g., "G2" is a rank code indicating
#     a taxon is between genus and species and the grandparent taxon is at the genus rank.
#  5. NCBI taxonomic ID number
#  6. Indented scientific name
#   1.65  31100   31100   U       0       unclassified
#  98.35  1859194 9       R       1       root
#  98.35  1859106 27      R1      131567    cellular organisms
#  98.33  1858749 159     D       2           Bacteria
#  98.29  1857899 62      P       1224          Proteobacteria
#  98.28  1857810 6002    C       1236            Gammaproteobacteria
#  97.73  1847477 0       O       135623            Vibrionales
#  97.73  1847477 1976    F       641                 Vibrionaceae
#  97.53  1843636 62018   G       662                   Vibrio
#  92.05  1739936 98729   G1      717610                  Vibrio harveyi group
#  80.50  1521696 1521696 S       669                       Vibrio harveyi
#   1.92  36269   36269   S       696485                    Vibrio owensii
#   1.62  30678   25083   S       680                       Vibrio campbellii
#   0.21  3914    3914    S1      338187                      Vibrio campbellii ATCC BAA-1116
#   0.09  1681    1681    S1      1224742                     Vibrio campbellii CAIM 519 = NBRC 15631

if __name__ == '__main__':
    main()
