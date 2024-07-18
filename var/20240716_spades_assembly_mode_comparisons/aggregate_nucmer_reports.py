#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2024 Adrian Tan <adrian_tan@nparks.gov.sg>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
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


@click.command()
@click.argument("report_files", nargs=-1)
@click.option("-t", "--tag", required=True, help="tag")
@click.option("-o", "--output_file", required=True, help="output data text file")
def main(report_files, tag, output_file):
    """
    aggregate_nucmer_reports.py

    e.g. aggregate_nucmer_reports.py 1.report 2.report -t a_b -o a_b.txt
    """
    print("\t{0:<20} :   {1:<10}".format("report_files", report_files))
    print("\t{0:<20} :   {1:<10}".format("tag", tag))
    print("\t{0:<20} :   {1:<10}".format("output_file", output_file))

    samples = []
    idx = 0
    #aggregate statistics from report
    # 0:meta_viral/59_2_29_A66_23-11_ASFV_blood.contigs.fasta /meta/59_2_29_A66_23-11_ASFV_blood.contigs.fasta
    # 1: NUCMER
    # 2:
    # 3:                             [REF]                [QRY]
    # 4: [Sequences]
    # 5: TotalSeqs                         51               264808
    # 6: AlignedSeqs             48(94.1176%)         490(0.1850%)
    # 7: UnalignedSeqs             3(5.8824%)     264318(99.8150%)
    # 8:
    # 9: [Bases]
    #10: TotalBases                    148294             99993229
    #11: AlignedBases        138783(93.5864%)      208972(0.2090%)
    #12: UnalignedBases         9511(6.4136%)   99784257(99.7910%)
    #13:
    #14: [Alignments]
    #15: 1-to-1                           138                  138
    #16: TotalLength                   134376               134364
    #17: AvgLength                   973.7391             973.6522
    #18: AvgIdentity                  99.6538              99.6538
    #19:
    #20: M-to-M                           556                  556
    #21: TotalLength                   227699               227687
    #22: AvgLength                   409.5306             409.5090
    #23: AvgIdentity                  95.9291              95.9291

    #name
    #seq_type: isolate_meta, metaviral_isolate
    #total_seqs
    # aligned_seqs,
    # unaligned_seqs,
    # total_bases,
    #  aligned_bases,
    #  unaligned_bases,
    # oto_total_length, oto_avg_length, avg_identity, m_to_m, total_length, avg_length, avg_identity

    for f in report_files:
        line_no = 0
        with open(f, "r") as file:
            for line in file:
                if line_no == 1:
                elif line_no == 1:

                if not line.startswith("#"):
                    idx += 1
                    sample_id, contigs_file = line.rstrip().split("\t")
                    samples.append(Sample(idx, sample_id, contigs_file))



if __name__ == "__main__":
    main()  # type: ignore
