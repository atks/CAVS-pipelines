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
    default="",
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
@click.option(
    "-p",
    "--prefix",
    required=True,
    default="tree",
    show_default=True,
    help="for RAXML file naming",
)
def main(working_dir, fasta_file, ref_fasta_file, prefix):
    """
    Generates a phylogenetic tree from a panel of sequences and reference sequences

    e.g. make_phylogenetic_tree.py -s lsdv.fasta -r lsdv_ref.fasta
    """
    output_dir = f"{os.getcwd()}/phylo"
    if working_dir != "":
        output_dir = os.path.abspath(working_dir)
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # version
    version = "1.0.0"

    # programs
    mafft = "/usr/local/mafft-7.490/bin/mafft"
    raxml = "/usr/local/raxml-ng-1.1.0/raxml-ng"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/make_phylogenetic_tree.log")

    # log text
    log_text = ""

    # if fasta_file not empty, compare sequences to reference
    desc = (
        f"Generating combined FASTA file with clean IDs for multiple sequence alignment"
    )
    if fasta_file is None:
        fasta_file = ""
        desc = f"Generating FASTA file with clean IDs from reference FASTA only for multiple sequence alignment"
    combined_fasta_file = f"{output_dir}/combined.fasta"
    cmd = fr'cat {fasta_file} {ref_fasta_file} | {seqkit} replace -p "[\s;:,\(\)\']" -r "_"  > {combined_fasta_file}'
    tgt = f"{combined_fasta_file}.OK"
    mpm.run(cmd, tgt, desc)

    # perform multiple sequence alignment
    input_fasta_file = f"{output_dir}/{ref_fasta_file}"
    if fasta_file is not None:
        input_fasta_file = f"{output_dir}/combined.fasta"
    output_fasta_file = f"{output_dir}/msa.fasta"
    log = f"{output_dir}/msa.log"
    cmd = f"{mafft} {input_fasta_file} > {output_fasta_file} 2>{log}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Multiple sequence alignment"
    mpm.run(cmd, tgt, desc)

    # construct phylogenetic tree with bootstrap
    input_msa_fasta_file = f"{output_dir}/msa.fasta"
    log = f"{output_dir}/construct_trees.log"
    cmd = f"cd {output_dir}; {raxml} --threads 10 --msa {input_msa_fasta_file} --model GTR+G --prefix {prefix} --bootstrap > {log}"
    tgt = f"{output_dir}/construct_trees.OK"
    desc = f"Constructing phylogenetic tree"
    mpm.run(cmd, tgt, desc)

    # construct consensus tree
    input_msa_fasta_file = f"{output_dir}/msa.fasta"
    log = f"{output_dir}/consensus_tree.log"
    cmd = f"cd {output_dir}; {raxml} --consense MRE --tree {prefix}.raxml.bootstraps --redo --prefix {prefix} > {log}"
    tgt = f"{output_dir}/consensus_tree.OK"
    desc = f"Constructing consensus tree"
    mpm.run(cmd, tgt, desc)

    # copy files to trace
    copy2(__file__, trace_dir)

    # write log file
    mpm.print_log()


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
