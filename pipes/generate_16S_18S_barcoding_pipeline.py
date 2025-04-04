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
    default="16S_18S_barcoding_pipeline.mk",
    help="make file name",
)
@click.option("-n", "--fastq_file", required=True, help="Nanopore fastq reads to be classified")
@click.option(
    "-w",
    "--working_dir",
    default= f"{os.getcwd()}/16S_18S_barcoding",
    show_default=True,
    help="working directory",
)
def main(make_file, fastq_file, working_dir):
    """
    Classifies reads by 16S and 18S barcoding

    e.g. generate_16S_18S_barcoding_pipeline.py -n nanopore.fastq.gz
    """
    working_dir = os.path.abspath(working_dir)
    ref_fasta_file = "/net/singapura/vfp/ref/16s_18s/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
    fastq_file = os.path.abspath(fastq_file)

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("fastq_file", fastq_file))
    print("\t{0:<20} :   {1:<10}".format("ref_fasta_file", ref_fasta_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    #version
    version = "1.0.0"

    # initialize
    pg = PipelineGenerator(make_file)

    # programs
    seqkit = "/usr/local/seqkit-2.20/seqkit"

    # create directories in destination folder directory
    trace_dir = f"{working_dir}/trace"
    try:
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    ##########
    # minimap2
    ##########

    # build index
    output_dir = f"{working_dir}/fasta_split"
    log = f"{output_dir}/fasta_split.log"
    dep = f""
    tgt = f"{working_dir}/fasta_split.OK"
    cmd = f"{seqkit} split -f --by-id {fastq_file} -O {output_dir} > {log}"
    pg.add(tgt, dep, cmd)

    # align

    # process alignment

    # prepare text file for krona tools

    # generate krona plot

    ##########
    # blast
    ##########

    ##########
    # qiime 2
    ##########

    ##########
    # kraken 2
    ##########



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
