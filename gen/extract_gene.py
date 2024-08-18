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
from shutil import copy2

@click.command()
@click.option(
    "-o",
    "--output_dir",
    default=os.getcwd(),
    show_default=True,
    help="output directory",
)
@click.option(
    "-g",
    "--gene_fasta_file",
    required=True,
    show_default=True,
    help="gene fasta file",
)
@click.option(
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
@click.option(
    "-h",
    "--extracted_gene_fasta_header",
    required=False,
    default="gene",
    show_default=True,
    help="extracted gene FASTA header",
)
def main(
    output_dir, gene_fasta_file, reference_fasta_file, extracted_gene_fasta_header
):
    """
    Extracts gene from a reference sequence file based on a gene sequence

    e.g. extract_gene -g gene.fasta -r ref.fasta
    """

    # create working directory
    # write out to separate files
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")
        exit(1)

    # version
    version = "1.0.0"

    # programs
    water = "/usr/local/emboss-6.6.0/bin/water"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/extract_amplicon.log")
    mpm.set_ignore_targets(True)
    
    # read reference sequences
    seq = ""
    with open(reference_fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                pass
            else:
                seq += line.rstrip()

    # read gene sequence
    gene = ""
    with open(gene_fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                pass
            else:
                gene += line.rstrip()

    # reverse_complement
    gene_fwd = gene.upper()
    gene_rev = (
        gene[::-1]
        .replace("A", "t")
        .replace("T", "a")
        .replace("G", "c")
        .replace("C", "g")
        .upper()
    )

    # write out to separate files
    output_fasta_file = f"{output_dir}/gene_fwd.fasta"
    cmd = f"echo '>gene_fwd\n{gene_fwd}' > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Create forward gene FASTA file"
    mpm.run(cmd, tgt, desc)

    output_fasta_file = f"{output_dir}/gene_rev.fasta"
    cmd = f"echo '>gene_rev\n{gene_rev}' > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Create reverse gene FASTA file"
    mpm.run(cmd, tgt, desc)

    # invoke smith waterman alignment for each direction
    input_fasta_file = f"{output_dir}/gene_fwd.fasta"
    alignment_file = f"{output_dir}/gene_fwd_ref.water"
    cmd = f"{water} {input_fasta_file} {reference_fasta_file} {alignment_file} -gapopen 10 -gapextend 0.5"
    tgt = f"{alignment_file}.OK"
    desc = f"Align forward gene to reference"
    mpm.run(cmd, tgt, desc)

    input_fasta_file = f"{output_dir}/gene_rev.fasta"
    alignment_file = f"{output_dir}/gene_rev_ref.water"
    cmd = f"{water} {input_fasta_file} {reference_fasta_file} {alignment_file} -gapopen 10 -gapextend 0.5"
    tgt = f"{alignment_file}.OK"
    desc = f"Align reverse gene to reference"
    mpm.run(cmd, tgt, desc)

    # examine alignments
    gene_fwd_alignment = parse_water_alignment(f"{output_dir}/gene_fwd_ref.water")
    gene_rev_alignment = parse_water_alignment(f"{output_dir}/gene_rev_ref.water")


    print(f"fwd score: {gene_fwd_alignment.score}")
    print(f"rev score: {gene_rev_alignment.score}")
    # find best alignment.
    best_alignment = (
        gene_fwd_alignment
        if gene_fwd_alignment.score > gene_rev_alignment.score
        else gene_rev_alignment
    )

    print("best alignment")
    best_alignment.print()

    # write out extracted gene sequence
    output_fasta_file = f"{output_dir}/extracted_gene.fasta"
    gene_seq = seq[best_alignment.beg - 1 : best_alignment.end]
    print(f"orig: {gene_seq}")
    if best_alignment.qseq == "gene_rev":
        gene_seq = (
            gene_seq[::-1]
            .replace("A", "t")
            .replace("T", "a")
            .replace("G", "c")
            .replace("C", "g")
            .upper()
        )
    cmd = f"echo '>{extracted_gene_fasta_header}\n{gene_seq}' > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Extract gene and save in FASTA file"
    mpm.run(cmd, tgt, desc)

    print(f"Extracted gene stats")
    print(f'gene  : {gene[0:20]}{"..." if len(gene)>20 else ""} ({len(gene)}bp)')
    print(f'reference : {seq[0:30]}{"..." if len(seq)>30 else ""} ({len(seq)}bp)')
    print(f"best alignment  : {best_alignment.qseq}")
    print(f"size: {best_alignment.length}")
    print(f"locus: {best_alignment.qseq}:{best_alignment.beg}-{best_alignment.end}")
    print("==========")
    gene_fwd_alignment.print()
    print("==========")
    gene_rev_alignment.print()
    print("==========")

    # copy files to trace
    copy2(__file__, trace_dir)

    # write log file
    mpm.print_log()


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
                m = re.search(r"\# 1: (.+)", line)
                qseq = m.group(1)
            elif line.startswith("# 2:"):
                m = re.search(r"\# 2: (.+)", line)
                rseq = m.group(1)
            elif line.startswith("# Length"):
                m = re.search(r"(\d+)", line)
                length = m.group(1)
            elif line.startswith("# Identity"):
                m = re.search(r"(\d+)\/", line)
                identity = m.group(1)
            elif line.startswith("# Gaps"):
                m = re.search(r"(\d+)\/", line)
                gaps = m.group(1)
            elif line.startswith("# Score"):
                m = re.search(r"([\d\.]+)", line)
                score = m.group(1)
            elif len(rseq) != 0 and line.startswith(rseq[:13]):
                m = re.search(fr"{rseq[:13]}\s+([\d\.]+) [^\d]+ ([\d\.]+)", line)
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


class MiniPipeManager(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log_msg = []
        self.ignore_targets = False

    def run(self, cmd, tgt, desc):
        try:
            if not self.ignore_targets and os.path.exists(tgt):
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

    def set_ignore_targets(self, ignore_targets):
        self.ignore_targets = ignore_targets

    def log(self, msg):
        print(msg)
        self.log_msg.append(msg)

    def print_log(self):
        self.log(f"\nlogs written to {self.log_file}")
        with open(self.log_file, "w") as f:
            f.write("\n".join(self.log_msg))



if __name__ == "__main__":
    main()
