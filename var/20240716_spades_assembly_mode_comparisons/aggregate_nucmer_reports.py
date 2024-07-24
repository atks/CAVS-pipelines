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
import re

@click.command()
@click.argument("report_files", nargs=-1)
@click.option("-o", "--output_file", required=True, help="output data text file")
def main(report_files, output_file):
    """
    aggregate_nucmer_reports.py

    e.g. aggregate_nucmer_reports.py 1.report 2.report -t a_b -o a_b.txt
    """
    #print("\t{0:<20} :   {1:<10}".format("report_files", report_files))
    print("\t{0:<20} :   {1:<10}".format("output_file", output_file))

    samples = []
    idx = 0
    #aggregate statistics from report
    # 0:meta_viral/59_2_29_A66_23-11_ASFV_blood.contigs.fasta /meta/59_2_29_A66_23-11_ASFV_blood.contigs.fasta
    # 1: NUCMER
    # 2:
    # 3:                             [REF]                [QRY]
    # 4: [Sequences]                                               - Sequence-centric stats.
    # 5: TotalSeqs                         51               264808 - Total number of input sequences.
    # 6: AlignedSeqs             48(94.1176%)         490(0.1850%) - Number of input sequences with at least one alignment.
    # 7: UnalignedSeqs             3(5.8824%)     264318(99.8150%) - Number of input sequences with no alignment.
    # 8:
    # 9: [Bases]                                                   - Base-pair-centric stats.
    #10: TotalBases                    148294             99993229 - Total number of bases in the input sequences.
    #11: AlignedBases        138783(93.5864%)      208972(0.2090%) - Total number of bases contained within an alignment.
    #12: UnalignedBases         9511(6.4136%)   99784257(99.7910%) - Total number of unaligned bases. This is a rough
    #                                                                measure for the amount of "unique" sequence in the
    #                                                                reference and query.
    #13:
    #14: [Alignments]                                              - Alignment-centric stats.
    #15: 1-to-1                           138                  138 - Number of alignment blocks comprising the 1-
    #                                                                mapping of reference to query. This is a sub
    #                                                                the M-to-M mapping, with repeats removed.
    #16: TotalLength                   134376               134364 - Total length of 1-to-1 alignment blocks.
    #17: AvgLength                   973.7391             973.6522 - Average length of 1-to-1 alignment blocks.
    #18: AvgIdentity                  99.6538              99.6538 - Average identity of 1-to-1 alignment blocks.
    #19:
    #20: M-to-M                           556                  556 - Number of alignment blocks comprising the
    #                                                                many-to-many mapping of reference to query. The
    #                                                                M-to-M mapping represents the smallest set of
    #                                                                alignments that maximize the coverage of both
    #                                                                reference and query. This is a superset of the 1
    #                                                                mapping.
    #21: TotalLength                   227699               227687 - Total length of M-to-M alignment blocks.
    #22: AvgLength                   409.5306             409.5090 - Average length of M-to-M alignment blocks.
    #23: AvgIdentity                  95.9291              95.9291 - Average identity of M-to-M alignment blocks.

    #name
    #seq_type: isolate_meta, metaviral_isolate
    #total_seqs, aligned_seqs,
    # unaligned_seqs,
    # total_bases,
    #  aligned_bases,
    #  unaligned_bases,
    # oto_total_length, oto_avg_length, avg_identity, m_to_m, total_length, avg_length, avg_identity
    file_no = 0

    with open(output_file, "w") as ofile:

        #print header
        output_line = f"file_no\ttype\ttotal_seqs\taligned_seqs\tunaligned_seqs"
        output_line += f"\ttotal_bases\taligned_bases\tunaligned_bases"
        output_line += f"\t1to1_total_length\t1to1_aligned_length\t1to1_aligned_identity"
        output_line += f"\tmtom_total_length\tmtom_aligned_length\tmtom_aligned_identity\n"

        # output_line = f"file_no\tref_total_seqs\tqry_total_seqs\tref_aligned_seqs\tqry_aligned_seqs\tref_unaligned_seqs\tqry_unaligned_seqs"
        # output_line += f"\tref_total_bases\tqry_total_bases\tref_aligned_bases\tqry_aligned_bases\tref_unaligned_bases\tqry_unaligned_bases"
        # output_line += f"\tref_1to1_total_length\tqry_1to1_total_length\tref_1to1_aligned_length\tqry_1to1_aligned_length\tref_1to1_aligned_identity\tqry_1to1_aligned_identity"
        # output_line += f"\tref_mtom_total_length\tqry_mtom_total_length\tref_mtom_aligned_length\tqry_mtom_aligned_length\tref_mtom_aligned_identity\tqry_mtom_aligned_identity\n"
        ofile.write(output_line)

        for f in report_files:
            with open(f, "r") as file:
                ref_total_seqs = ""
                qry_total_seqs = ""
                ref_aligned_seqs = ""
                qry_aligned_seqs = ""
                ref_unaligned_seqs = ""
                qry_unaligned_seqs = ""

                ref_total_bases = ""
                qry_total_bases = ""
                ref_aligned_bases = ""
                qry_aligned_bases = ""
                ref_unaligned_bases = ""
                qry_unaligned_bases = ""

                ref_total_seqs = ""
                qry_total_seqs = ""
                ref_aligned_seqs = ""
                qry_aligned_seqs = ""
                ref_unaligned_seqs = ""
                qry_unaligned_seqs = ""

                ref_1to1_total_length = ""
                qry_1to1_total_length = ""
                ref_1to1_aligned_length = ""
                qry_1to1_aligned_length = ""
                ref_1to1_aligned_identity = ""
                qry_1to1_aligned_identity = ""

                ref_mtom_total_length = ""
                qry_mtom_total_length = ""
                ref_mtom_aligned_length = ""
                qry_mtom_aligned_length = ""
                ref_mtom_aligned_identity = ""
                qry_mtom_aligned_identity = ""

                line_no = 0
                for line in file:
                    if line_no in [5, 6, 7, 10, 11, 12, 16, 17, 18, 21, 22, 23]:
                        pass
                        #print(f"{line_no} {line}")
                    if line_no == 5:
                        m = re.search(r"TotalSeqs\s+(\d+)\s+(\d+)", line)
                        if m is not None:
                            ref_total_seqs = m.group(1)
                            qry_total_seqs = m.group(2)
                    elif line_no == 6:
                        m = re.search(r"AlignedSeqs\s+(\d+)\(.+\)\s+(\d+)\(.+\)", line)
                        if m is not None:
                            ref_aligned_seqs = m.group(1)
                            qry_aligned_seqs = m.group(2)
                    elif line_no == 7:
                        m = re.search(r"UnalignedSeqs\s+(\d+)\(.+\)\s+(\d+)\(.+\)", line)
                        if m is not None:
                            ref_unaligned_seqs = m.group(1)
                            qry_unaligned_seqs = m.group(2)
                    elif line_no == 10:
                        m = re.search(r"TotalBases\s+(\d+)\s+(\d+)", line)
                        if m is not None:
                            ref_total_bases = m.group(1)
                            qry_total_bases = m.group(2)
                    elif line_no == 11:
                        m = re.search(r"AlignedBases\s+(\d+)\(.+\)\s+(\d+)\(.+\)", line)
                        if m is not None:
                            ref_aligned_bases = m.group(1)
                            qry_aligned_bases = m.group(2)
                    elif line_no == 12:
                        m = re.search(r"UnalignedBases\s+(\d+)\(.+\)\s+(\d+)\(.+\)", line)
                        if m is not None:
                            ref_unaligned_bases = m.group(1)
                            qry_unaligned_bases = m.group(2)
                    elif line_no == 16:
                        m = re.search(r"TotalLength\s+(\d+)\s+(\d+)", line)
                        if m is not None:
                            ref_1to1_total_length = m.group(1)
                            qry_1to1_total_length = m.group(2)
                    elif line_no == 17:
                        m = re.search(r"AvgLength\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
                        if m is not None:
                            ref_1to1_aligned_length = m.group(1)
                            qry_1to1_aligned_length = m.group(2)
                    elif line_no == 18:
                        m = re.search(r"AvgIdentity\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
                        if m is not None:
                            ref_1to1_aligned_identity = m.group(1)
                            qry_1to1_aligned_identity = m.group(2)
                    elif line_no == 21:
                        m = re.search(r"TotalLength\s+(\d+)\s+(\d+)", line)
                        if m is not None:
                            ref_mtom_total_length = m.group(1)
                            qry_mtom_total_length = m.group(2)
                    elif line_no == 22:
                        m = re.search(r"AvgLength\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
                        if m is not None:
                            ref_mtom_aligned_length = m.group(1)
                            qry_mtom_aligned_length = m.group(2)
                    elif line_no == 23:
                        m = re.search(r"AvgIdentity\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
                        if m is not None:
                            ref_mtom_aligned_identity = m.group(1)
                            qry_mtom_aligned_identity = m.group(2)
                    line_no += 1

            output_line = f"{file_no}\tref\t{ref_total_seqs}\t{ref_aligned_seqs}\t{ref_unaligned_seqs}"
            output_line += f"\t{ref_total_bases}\t{ref_aligned_bases}\t{ref_unaligned_bases}"
            output_line += f"\t{ref_1to1_total_length}\t{ref_1to1_aligned_length}\t{ref_1to1_aligned_identity}"
            output_line += f"\t{ref_mtom_total_length}\t{ref_mtom_aligned_length}\t{ref_mtom_aligned_identity}\n"
            ofile.write(output_line)

            output_line = f"{file_no}\tquery\t{qry_total_seqs}\t{qry_aligned_seqs}\t{qry_unaligned_seqs}"
            output_line += f"\t{qry_total_bases}\t{qry_aligned_bases}\t{qry_unaligned_bases}"
            output_line += f"\t{qry_1to1_total_length}\t{qry_1to1_aligned_length}\t{qry_1to1_aligned_identity}"
            output_line += f"\t{qry_mtom_total_length}\t{qry_mtom_aligned_length}\t{qry_mtom_aligned_identity}\n"


            # output_line = f"{file_no}\t{ref_total_seqs}\t{qry_total_seqs}\t{ref_aligned_seqs}\t{qry_aligned_seqs}\t{ref_unaligned_seqs}\t{qry_unaligned_seqs}"
            # output_line += f"\t{ref_total_bases}\t{qry_total_bases}\t{ref_aligned_bases}\t{qry_aligned_bases}\t{ref_unaligned_bases}\t{qry_unaligned_bases}"
            # output_line += f"\t{ref_1to1_total_length}\t{qry_1to1_total_length}\t{ref_1to1_aligned_length}\t{qry_1to1_aligned_length}\t{ref_1to1_aligned_identity}\t{qry_1to1_aligned_identity}"
            # output_line += f"\t{ref_mtom_total_length}\t{qry_mtom_total_length}\t{ref_mtom_aligned_length}\t{qry_mtom_aligned_length}\t{ref_mtom_aligned_identity}\t{qry_mtom_aligned_identity}\n"
            ofile.write(output_line)

            file_no += 1

if __name__ == "__main__":
    main()  # type: ignore
