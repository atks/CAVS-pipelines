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

import sys
import os
import click
import subprocess
import re
from shutil import copy2


@click.command()
@click.option(
    "-o",
    "--output_dir",
    default=f"{os.getcwd()}/extract_amplicon_output",
    show_default=True,
    help="output directory",
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
def main(output_dir, primer_fasta_file, reference_fasta_file):
    """
    Extracts amplicon from a reference sequence file and a pair of primers

    e.g. extract_amplicon -p primers.fasta -r ref.fasta
    """

    # version
    version = "1.0.0"

    # programs
    water = "/usr/local/emboss-6.6.0/bin/water"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/extract_amplicon.log")

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
    seq1_reverse = reverse_complement(seq1_forward)
    seq2_forward = seq2.upper()
    seq2_reverse = reverse_complement(seq2_forward)

    # write out to separate files
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    output_file = f"{output_dir}/p1_forward.fasta"
    with open(output_file, "w") as file:
        file.write(">primer1_forward\n")
        file.write(seq1_forward)
        file.write("\n")

    output_file = f"{output_dir}/p1_reverse.fasta"
    with open(output_file, "w") as file:
        file.write(">primer1_reverse\n")
        if seq1_reverse is not None:
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
        if seq2_reverse is not None:
            file.write(seq2_reverse)
        file.write("\n")

    # invoke smith waterman alignment for each direction
    primer_fasta_file = f"{output_dir}/p1_forward.fasta"
    water_alignment_file = f"{output_dir}/p1_forward.water"
    cmd = f"{water} {primer_fasta_file} {reference_fasta_file}  {water_alignment_file} -gapopen 10 -gapextend 0.5"
    tgt = f"{primer_fasta_file}.OK"
    desc = f"P1F sequence alignment"
    mpm.run(cmd, tgt, desc)
    p1_forward_alignment = parse_water_alignment(water_alignment_file)

    primer_fasta_file = f"{output_dir}/p1_reverse.fasta"
    water_alignment_file = f"{output_dir}/p1_reverse.water"
    cmd = f"{water} {primer_fasta_file} {reference_fasta_file}  {water_alignment_file} -gapopen 10 -gapextend 0.5"
    tgt = f"{primer_fasta_file}.OK"
    desc = f"P1R sequence alignment"
    mpm.run(cmd, tgt, desc)
    p1_reverse_alignment = parse_water_alignment(water_alignment_file)

    primer_fasta_file = f"{output_dir}/p2_forward.fasta"
    water_alignment_file = f"{output_dir}/p2_forward.water"
    cmd = f"{water} {primer_fasta_file} {reference_fasta_file}  {water_alignment_file} -gapopen 10 -gapextend 0.5"
    tgt = f"{primer_fasta_file}.OK"
    desc = f"P2F sequence alignment"
    mpm.run(cmd, tgt, desc)
    p2_forward_alignment = parse_water_alignment(water_alignment_file)

    primer_fasta_file = f"{output_dir}/p2_reverse.fasta"
    water_alignment_file = f"{output_dir}/p2_reverse.water"
    cmd = f"{water} {primer_fasta_file} {reference_fasta_file} {water_alignment_file} -gapopen 10 -gapextend 0.5"
    tgt = f"{primer_fasta_file}.OK"
    desc = f"P2R sequence alignment"
    mpm.run(cmd, tgt, desc)
    p2_reverse_alignment = parse_water_alignment(water_alignment_file)

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

    # copy files to trace
    copy2(__file__, trace_dir)

    # write log file
    mpm.print_log()

def reverse_complement(seq):
    return (
        seq[::-1]
        .upper()
        .replace("A", "t")
        .replace("T", "a")
        .replace("G", "c")
        .replace("C", "g")
        .replace("M", "k")
        .replace("R", "y")
        .replace("W", "w")
        .replace("S", "s")
        .replace("Y", "r")
        .replace("K", "m")
        .replace("V", "b")
        .replace("H", "d")
        .replace("D", "h")
        .replace("B", "v")
        .replace("N", "n")
        .upper()
    )

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
        align = ""
        for line in file:
            if line.startswith("# 1:"):
                m = re.search(r"\# 1: (.+)", line)
                if m is not None:
                    qseq = m.group(1)
            elif line.startswith("# 2:"):
                m = re.search(r"\# 2: (.+)", line)
                if m is not None:
                    rseq = m.group(1)
            elif line.startswith("# Length"):
                m = re.search(r"(\d+)", line)
                if m is not None:
                   length = m.group(1)
            elif line.startswith("# Identity"):
                m = re.search(r"(\d+)\/", line)
                if m is not None:
                    identity = int(m.group(1))
            elif line.startswith("# Gaps"):
                m = re.search(r"(\d+)\/", line)
                if m is not None:
                    gaps = m.group(1)
            elif line.startswith("# Score"):
                m = re.search(r"([\d\.]+)", line)
                if m is not None:
                    score = m.group(1)
            elif len(rseq) != 0 and line.startswith(rseq[:13]):
                m = re.search(fr"{rseq[:13]}\s+([\d\.]+) [^\d]+ ([\d\.]+)", line)
                if m is not None:
                    beg = int(m.group(1)) if int(m.group(1)) < beg else beg
                    end = int(m.group(2)) if int(m.group(2)) > end else end

            if not line.startswith("#") and len(line) > 1:
                align += line

        return Alignment(qseq, rseq, length, identity, gaps, score, beg, end, align)


class Alignment(object):
    def __init__(self, qseq, rseq, length, identity, gaps, score, beg, end, align):
        self.qseq = qseq
        self.rseq = rseq
        self.length = length
        self.identity = identity
        self.gaps = gaps
        self.score = score
        self.beg = beg
        self.end = end
        self.align = align

    def print(self):
        print(f"qseq      : {self.qseq}")
        print(f"rseq      : {self.rseq}")
        print(f"length    : {self.length}")
        print(f"identity  : {self.identity}")
        print(f"gaps      : {self.gaps}")
        print(f"score     : {self.score}")
        print(f"beg       : {self.beg}")
        print(f"end       : {self.end}")
        print(f"\n{self.align}")

class MiniPipeManager(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log_msg = []

    def run(self, cmd, tgt, desc):
        try:
            if os.path.exists(tgt):
                self.log(f"{desc} -  already executed")
                self.log(cmd)
                return
            else:
                self.log(f"{desc}")
                subprocess.run(cmd, shell=True, check=True)
                subprocess.run(f"touch {tgt}", shell=True, check=True)
                self.log(cmd)
        except subprocess.CalledProcessError as e:
            self.log(f" - failed")
            exit(1)

    def log(self, msg):
        print(msg)
        self.log_msg.append(msg)

    def print_log(self):
        self.log(f"\nlogs written to {self.log_file}")
        with open(self.log_file, "w") as f:
            f.write("\n".join(self.log_msg))



if __name__ == "__main__":
    main() # type: ignore[arg-type]
