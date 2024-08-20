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
import sys
import re
import sys
from shutil import copy2
from datetime import datetime


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="sequence_orientation_evaluation_pipeline.mk",
    help="make file name",
)
@click.option("-f", "--fasta_file", required=True, help="FASTA file containing the sequences to be checked")
@click.option("-r", "--ref_fasta_file", required=True, help="reference sequence file to compare against")
@click.option(
    "-w",
    "--working_dir",
    default= f"{os.getcwd()}/sequence_orientation_evaluation",
    show_default=True,
    help="working directory",
)
def main(make_file, fasta_file, ref_fasta_file, working_dir):
    """
    Compares the orientation of a set of sequences against a reference sequence.
    This is motivated by phylogenetic tree building where the multiple sequence alignment has to be in the same orientation

    e.g. generate_sequence_orientation_evaluation_pipeline -f panel.fasta -i ref.fasta -w verify_orientation
    """
    working_dir = os.path.abspath(working_dir)
    fasta_file = os.path.abspath(fasta_file)
    prefix_fasta_file = os.path.basename(fasta_file).split(".")[0]
    ref_fasta_file = os.path.abspath(ref_fasta_file)

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("fasta_file", fasta_file))
    print("\t{0:<20} :   {1:<10}".format("ref_fasta_file", ref_fasta_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    #version
    version = "1.0.0"

    # initialize
    pg = PipelineGenerator(make_file)

    # programs
    seqkit = "/usr/local/seqkit-2.20/seqkit"
    compare_sequence_orientation = "/home/atks/programs/CAVS-pipelines/gen/compare_sequence_orientation.py"

    ids = []

    #read through the fasta file and record the IDs
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                m = re.search(r"(?<=\>)(.*?)(?=\s)", line)
                if m is not None:
                    ids.append(m.group(0))

    # create directories in destination folder directory
    trace_dir = f"{working_dir}/trace"
    split_dir = f"{working_dir}/fasta_split"
    orientation_dir = f"{working_dir}/orientation"
    try:
        os.makedirs(trace_dir, exist_ok=True)
        os.makedirs(split_dir, exist_ok=True)
        os.makedirs(orientation_dir, exist_ok=True)
        for id in ids:
            os.makedirs(f"{orientation_dir}/{id}", exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # split fasta file
    # seqkit split --by-id ../36seq_asfv_ref_p72_genotype.fasta
    output_dir = f"{working_dir}/fasta_split"
    log = f"{output_dir}/fasta_split.log"
    dep = f""
    tgt = f"{working_dir}/fasta_split.OK"
    cmd = f"{seqkit} split -f --by-id {fasta_file} -O {output_dir} > {log}"
    pg.add(tgt, dep, cmd)

    #perform pairwise orientation checking
    compare_sequence_orientation_dep = ""
    compare_sequence_orientation_reports = ""
    for id in ids:
        query_fasta_file = f"{split_dir}/{prefix_fasta_file}.id_{id}.fasta"
        output_dir = f"{orientation_dir}/{id}"
        log = f"{output_dir}/compare_sequence_orientation.log"
        dep = f"{working_dir}/fasta_split.OK"
        tgt = f"{output_dir}/{id}.orientation_check.OK"
        compare_sequence_orientation_dep += f" {tgt}"
        compare_sequence_orientation_reports += f" {output_dir}/result.txt"
        cmd = f"{compare_sequence_orientation} -q {query_fasta_file} -r {ref_fasta_file} -o {output_dir} > {log}"
        pg.add(tgt, dep, cmd)

    #compile orientation detection results
    # seqkit split --by-id ../36seq_asfv_ref_p72_genotype.fasta
    output_text_file = f"{working_dir}/orientation_report.txt"
    dep = f"{compare_sequence_orientation_dep}"
    tgt = f"{output_text_file}.OK"
    cmd = f"cat {compare_sequence_orientation_reports} | sort | uniq > {output_text_file}"
    pg.add(tgt, dep, cmd)

    #fix orientation of sequences


    # write make file
    print("Writing pipeline")
    pg.write()

    #copy files to trace
    copy2(__file__, trace_dir)
    copy2(make_file, trace_dir)

class PipelineGenerator(object):
    def __init__(self, make_file):
        self.make_file = make_file
        self.tgts = []
        self.deps = []
        self.cmds = []
        self.clean_cmd = ""

    def add_srun(self, tgt, dep, cmd, cpu):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(f"srun --mincpus {cpu} {cmd}")

    def add(self, tgt, dep, cmd):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(cmd)

    def add_clean(self, cmd):
        self.clean_cmd = cmd

    def write(self):
        with open(self.make_file, "w") as f:
            f.write("SHELL:=/bin/bash\n")
            f.write(".DELETE_ON_ERROR:\n\n")
            f.write("all : ")
            for i in range(len(self.tgts)):
                f.write(f"{self.tgts[i]} ")
            f.write("\n\n")

            for i in range(len(self.tgts)):
                f.write(f"{self.tgts[i]} : {self.deps[i]}\n")
                f.write(f"\t{self.cmds[i]}\n")
                f.write(f"\ttouch {self.tgts[i]}\n\n")

            if self.clean_cmd != "":
                f.write(f"clean : \n")
                f.write(f"\t{self.clean_cmd}\n")

if __name__ == "__main__":
    main() # type: ignore
