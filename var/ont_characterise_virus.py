#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2023 Adrian Tan <adrian_tan@nparks.gov.sg>
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

import sys
import os
import click
import subprocess
import re


@click.command()
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-p",
    "--primer_fasta_file",
    required=True,
    show_default=True,
    help="primers fasta file",
)
@click.option(
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
def main(working_dir, primer_fasta_file, reference_fasta_file):
    """
    Extracts amplicon from a reference sequence file and a pair of primers
    from reads

    e.g. characterise_virus -p primers.fasta -r ref.fasta
    """

    # read reference sequences
    seq = ""
    with open(reference_fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                pass
            else:
                seq += line.rstrip()

    # read primer sequences
    seq1 = ""
    seq2 = ""
    no_seq = 0
    with open(primer_fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                no_seq += 1
            else:
                if no_seq == 1:
                    seq1 += line.rstrip()
                elif no_seq == 2:
                    seq2 += line.rstrip()
                else:
                    exit("Only 2 sequences are expected")

    # reverse_complement
    seq1_forward = seq1.upper()
    seq1_reverse = (
        seq1_forward[::-1]
        .replace("A", "t")
        .replace("T", "a")
        .replace("G", "c")
        .replace("C", "g")
        .upper()
    )
    seq2_forward = seq2.upper()
    seq2_reverse = (
        seq2_forward[::-1]
        .replace("A", "t")
        .replace("T", "a")
        .replace("G", "c")
        .replace("C", "g")
        .upper()
    )

    # write out to separate files
    output_dir = f"{working_dir}/extract_amplicon_output"
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    output_file = f"{output_dir}/p1_forward.fasta"
    with open(output_file, "w") as file:
        file.write(">primer1_forward\n")
        file.write(seq1_forward)
        file.write("\n")

    output_file = f"{output_dir}/p1_reverse.fasta"
    with open(output_file, "w") as file:
        file.write(">primer1_reverse\n")
        file.write(seq1_reverse)
        file.write("\n")

    output_file = f"{output_dir}/p2_forward.fasta"
    with open(output_file, "w") as file:
        file.write(">primer2_forward\n")
        file.write(seq2_forward)
        file.write("\n")

    output_file = f"{output_dir}/p2_reverse.fasta"
    with open(output_file, "w") as file:
        file.write(">primer2_reverse\n")
        file.write(seq2_reverse)
        file.write("\n")

    # invoke smith waterman alignment for each direction
    water = "/usr/local/emboss-6.6.0/bin/water"
    p1_forward_water_alignment_file = f"{output_dir}/p1_forward.water"
    p1_reverse_water_alignment_file = f"{output_dir}/p1_reverse.water"
    p2_forward_water_alignment_file = f"{output_dir}/p2_forward.water"
    p2_reverse_water_alignment_file = f"{output_dir}/p2_reverse.water"

    # result = subprocess.run(
    # [f'{water} {output_dir}/p1_forward.fasta {reference_fasta_file} output.water -gapopen 10 -gapextend 0.5'], capture_output=True)
    result = subprocess.run(
        [
            water,
            f"{output_dir}/p1_forward.fasta",
            reference_fasta_file,
            p1_forward_water_alignment_file,
            "-gapopen",
            "10",
            "-gapextend",
            "0.5",
        ],
        capture_output=True,
    )
    result = subprocess.run(
        [
            water,
            f"{output_dir}/p1_reverse.fasta",
            reference_fasta_file,
            p1_reverse_water_alignment_file,
            "-gapopen",
            "10",
            "-gapextend",
            "0.5",
        ],
        capture_output=True,
    )
    result = subprocess.run(
        [
            water,
            f"{output_dir}/p2_forward.fasta",
            reference_fasta_file,
            p2_forward_water_alignment_file,
            "-gapopen",
            "10",
            "-gapextend",
            "0.5",
        ],
        capture_output=True,
    )
    result = subprocess.run(
        [
            water,
            f"{output_dir}/p2_reverse.fasta",
            reference_fasta_file,
            p2_reverse_water_alignment_file,
            "-gapopen",
            "10",
            "-gapextend",
            "0.5",
        ],
        capture_output=True,
    )

    # examine alignments
    p1_forward_alignment = parse_water_alignment(p1_forward_water_alignment_file)
    p1_reverse_alignment = parse_water_alignment(p1_reverse_water_alignment_file)
    p2_forward_alignment = parse_water_alignment(p2_forward_water_alignment_file)
    p2_reverse_alignment = parse_water_alignment(p2_reverse_water_alignment_file)

    # for good match - extract amplicon, report length
    # yes I know the issue here.  need to ensure it is the right pair and they are consistent with one another to amplify
    # later fix.
    best_p1_alignment = (
        p1_forward_alignment
        if p1_forward_alignment.identity > p1_reverse_alignment.identity
        else p1_reverse_alignment
    )
    best_p2_alignment = (
        p2_forward_alignment
        if p2_forward_alignment.identity > p2_reverse_alignment.identity
        else p2_reverse_alignment
    )

    # check overlap
    amplicon_size = 0
    overlap = (
        best_p1_alignment.end >= best_p2_alignment.beg
        and best_p1_alignment.beg <= best_p2_alignment.end
    )

    amplicon_beg = min(best_p1_alignment.beg, best_p2_alignment.beg)
    amplicon_end = max(best_p1_alignment.end, best_p2_alignment.end)
    amplicon_size = amplicon_end - amplicon_beg + 1

    print(f"Amplicon stats")
    print(f'primer 1  : {seq1[0:20]}{"..." if len(seq1)>20 else ""} ({len(seq1)}bp)')
    print(f'primer 2  : {seq2[0:20]}{"..." if len(seq2)>20 else ""} ({len(seq2)}bp)')
    print(f'reference : {seq[0:30]}{"..." if len(seq)>30 else ""} ({len(seq)}bp)')
    print(f"overlap: {overlap}")
    print(f"size: {amplicon_size}")
    print(f"locus: {best_p1_alignment.rseq}:{amplicon_beg}-{amplicon_end}")
    print("==========")
    best_p1_alignment.print()
    print("==========")
    best_p2_alignment.print()
    print("==========")

    # p1_forward_alignment.print()
    # p1_reverse_alignment.print()
    # p2_forward_alignment.print()
    # p2_reverse_alignment.print()


# Aligned_sequences: 2
# 1: primer2_reverse
# 2: AY261360.1
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 19
# Identity:      15/19 (78.9%)
# Similarity:    15/19 (78.9%)
# Gaps:           2/19 (10.5%)
# Score: 56.5
def parse_water_alignment(file):
    with open(file, "r") as file:
        qseq = ""
        rseq = ""
        length = 0
        identity = 0
        gaps = 0
        score = 0
        beg = 10000000000
        end = -1
        for line in file:
            if line.startswith("# 1:"):
                m = re.search("\# 1: (.+)", line)
                qseq = m.group(1)
            elif line.startswith("# 2:"):
                m = re.search("\# 2: (.+)", line)
                rseq = m.group(1)
            elif line.startswith("# Length"):
                m = re.search("(\d+)", line)
                length = m.group(1)
            elif line.startswith("# Identity"):
                m = re.search("(\d+)\/", line)
                identity = m.group(1)
            elif line.startswith("# Gaps"):
                m = re.search("(\d+)\/", line)
                gaps = m.group(1)
            elif line.startswith("# Score"):
                m = re.search("([\d\.]+)", line)
                score = m.group(1)
            elif len(rseq) != 0 and line.startswith(rseq):
                m = re.search(f"{rseq}\s+([\d\.]+) [^\d]+ ([\d\.]+)", line)
                beg = int(m.group(1)) if int(m.group(1)) < beg else beg
                end = int(m.group(2)) if int(m.group(2)) > end else end
            else:
                pass

        alignment = Alignment(qseq, rseq, length, identity, gaps, score, beg, end)
        return alignment


class Alignment(object):
    def __init__(self):
        self.qseq = ""
        self.rseq = ""
        self.length = 0
        self.identity = 0
        self.gaps = 0
        self.score = 0
        self.beg = 0
        self.end = 0

    def __init__(self, qseq, rseq, length, identity, gaps, score, beg, end):
        self.qseq = qseq
        self.rseq = rseq
        self.length = length
        self.identity = identity
        self.gaps = gaps
        self.score = score
        self.beg = beg
        self.end = end

    def print(self):
        print(f"qseq      : {self.qseq}")
        print(f"rseq      : {self.rseq}")
        print(f"length    : {self.length}")
        print(f"identity  : {self.identity}")
        print(f"gaps      : {self.gaps}")
        print(f"score     : {self.score}")
        print(f"beg       : {self.beg}")
        print(f"end       : {self.end}")


if __name__ == "__main__":
    main()
