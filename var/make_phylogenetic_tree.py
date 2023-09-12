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
    "-f",
    "--fasta_file",
    required=True,
    show_default=True,
    help="sequence panel fasta file",
)
@click.option(
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
def main(working_dir, sequence_panel_fasta_file):
    """
    Generates a phylogenetic tree from a panel of sequences and reference sequences

    e.g. make_phylogenetic_tree.py -s ndv.fasta -r ndv_ref.fasta
    """

    # create working directory
    output_dir = f"{working_dir}/phylo"
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    # mafft ndv_classII.fasta > ndv_classII.msa.fasta
    # raxml-ng  --msa ndv_classII.msa.fasta --model GTR+G --prefix  ndv_classII --bootstrap > phylo.log
    # raxml-ng --consense MRE --tree alltrees.nw --prefix consMRE
    # combine sequence files
    output_fasta_file = f"{output_dir}/gene_fwd.fasta"
    cmd = f"echo '>gene_fwd\n{gene_fwd}' > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Create forward gene FASTA file"
    run(cmd, tgt, desc)

    output_fasta_file = f"{output_dir}/gene_rev.fasta"
    cmd = f"echo '>gene_rev\n{gene_rev}' > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Create reverse gene FASTA file"
    run(cmd, tgt, desc)

    # invoke smith waterman alignment for each direction
    water = "/usr/local/emboss-6.6.0/bin/water"

    input_fasta_file = f"{output_dir}/gene_fwd.fasta"
    alignment_file = f"{output_dir}/gene_fwd_ref.water"
    cmd = f"{water} {input_fasta_file} {reference_fasta_file} {alignment_file} -gapopen 10 -gapextend 0.5"
    tgt = f"{alignment_file}.OK"
    desc = f"Align forward gene to reference"
    run(cmd, tgt, desc)

    # examine alignments
    gene_fwd_alignment = parse_water_alignment(f"{output_dir}/gene_fwd_ref.water")
    gene_rev_alignment = parse_water_alignment(f"{output_dir}/gene_rev_ref.water")

    # find best alignment.
    best_alignment = (
        gene_fwd_alignment
        if gene_fwd_alignment.score > gene_rev_alignment.score
        else gene_rev_alignment
    )

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
    run(cmd, tgt, desc)

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
            elif len(rseq) != 0 and line.startswith(rseq[:13]):
                m = re.search(f"{rseq[:13]}\s+([\d\.]+) [^\d]+ ([\d\.]+)", line)
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


def run(cmd, tgt, desc):
    try:
        if os.path.exists(tgt):
            print(f"{desc} -  already executed")
            return
        else:
            print(f"{cmd}")
            subprocess.run(cmd, shell=True, check=True)
            subprocess.run(f"touch {tgt}", shell=True, check=True)
            print(f"{desc} -  successfully executed")
    except subprocess.CalledProcessError as e:
        print(f" - failed")
        exit(1)


if __name__ == "__main__":
    main()
