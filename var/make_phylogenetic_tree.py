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
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-f",
    "--fasta_file",
    required=False,
    show_default=True,
    help="sequence panel fasta file",
)
@click.option(
    "-r",
    "--ref_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
def main(working_dir, fasta_file, ref_fasta_file):
    """
    Generates a phylogenetic tree from a panel of sequences and reference sequences

    e.g. make_phylogenetic_tree.py -s lsdv.fasta -r lsdv_ref.fasta
    """

    output_dir = f"{working_dir}/phylo"
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # if fasta_file not empty, compare sequences to reference
    combined_fasta_file = f"{output_dir}/combined.fasta"
    if fasta_file is None:
        fasta_file = ""
    combined_fasta_file = f"{output_dir}/combined.fasta"
    cmd = f"cat {fasta_file} {ref_fasta_file} | seqkit replace -p \"[\s;:,\(\)\']\" -r \"_\"  > {combined_fasta_file}"
    tgt = f"{combined_fasta_file}.OK"
    desc = f"Create combined FASTA file for multiple sequence alignement"
    run(cmd, tgt, desc)

    #programs
    mafft = "/usr/local/mafft-7.490/bin/mafft"
    raxml = "/usr/local/raxml-ng-1.1.0/raxml-ng"

    # perform multiple sequence alignment
    input_fasta_file = f"{output_dir}/{ref_fasta_file}"
    if fasta_file is not None:
        input_fasta_file = f"{output_dir}/combined.fasta"
    output_fasta_file = f"{output_dir}/msa.fasta"
    log = f"{output_dir}/msa.log"
    cmd = f"{mafft} {input_fasta_file} > {output_fasta_file} 2>{log}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Multiple sequence alignment"
    run(cmd, tgt, desc)

    # construct phylogenetic tree with bootstrap
    input_msa_fasta_file = f"{output_dir}/msa.fasta"
    log = f"{output_dir}/tree.log"
    cmd = f"cd {output_dir}; {raxml} --threads 10 --msa {input_msa_fasta_file} --model GTR+G --prefix tree --bootstrap > {log}"
    tgt = f"{output_dir}/construct_trees.OK"
    desc = f"Constructing phylogenetic tree"
    run(cmd, tgt, desc)

    # construct consensus tree
    input_msa_fasta_file = f"{output_dir}/msa.fasta"
    cmd = f"cd {output_dir}; {raxml} --consense MRE --tree tree.raxml.bootstraps --prefix consMRE "
    tgt = f"{output_dir}/consensus_tree.OK"
    desc = f"Constructing consensus tree"
    run(cmd, tgt, desc)

    #copy files to trace
    copy2(__file__, trace_dir)

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
