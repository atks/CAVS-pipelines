#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2022 Adrian Tan <adrian_tan@nparks.gov.sg>
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
    default="analyse_insert_length.mk",
    required=True,
    help="make file",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
@click.option(
    "-w", "--working_dir", default=os.getcwd(), required=True, help="working directory"
)
def main(make_file, working_dir, sample_file):
    """
    Aligns reads to reference and estimate insert length

    e.g. generate_analyse_insert_length_pipeline -s ilm38.fastq.txt
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                sample_id, fastq1, fastq2 = line.rstrip().split("\t")
                samples.append(Sample(sample_id, fastq1, fastq2))

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    log_dir = f"{working_dir}/log"
    bam_dir = f"{working_dir}/bam"
    stats_dir = f"{working_dir}/stats"
    try:
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # create index
    bwa = "/usr/local/bwa-0.7.17/bwa"
    tgt = f"{log_dir}/bam_index.OK"
    ref_fasta_file = f"{ref_dir}/NC_039223.1.fasta"
    dep = ""
    cmd = f"{bwa} index -a bwtsw {ref_fasta_file}"
    pg.add(tgt, dep, cmd)

    for idx, sample in enumerate(samples):
        # align reads
        bwa = "/usr/local/bwa-0.7.17/bwa"
        samtools = "/usr/local/samtools-1.14/samtools"
        bam_file = f"{bam_dir}/{sample.id}.bam"
        log = f"{log_dir}/{sample.id}.bam.log"
        err = f"{log_dir}/{sample.id}.bam.err"
        tgt = f"{log_dir}/{sample.id}.bam.OK"
        dep = f"{log_dir}/bam_index.OK"
        cmd = f"{bwa} mem -M {ref_fasta_file} {sample.fastq1} {sample.fastq2} | {samtools} sort - | {samtools} view - -o {bam_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # get stats
        stats_file = f"{stats_dir}/{sample.id}.stats.txt"
        log = f"{log_dir}/{sample.id}.stats.txt.log"
        err = f"{log_dir}/{sample.id}.stats.txt.err"
        tgt = f"{log_dir}/{sample.id}.stats.txt.OK"
        dep = f"{log_dir}/{sample.id}.bam.OK"
        cmd = f"{samtools} stats {bam_file} > {stats_file}"
        pg.add(tgt, dep, cmd)

    # write make file
    print("Writing pipeline")
    pg.write()


class PipelineGenerator(object):
    def __init__(self, make_file):
        self.make_file = make_file
        self.tgts = []
        self.deps = []
        self.cmds = []
        self.clean_cmd = ""

    def add(self, tgt, dep, cmd):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(cmd)

    def add_clean(self, cmd):
        self.clean_cmd = cmd

    def print(self):
        print(".DELETE_ON_ERROR:")
        for i in range(len(self.tgts)):
            print(f"{self.tgts[i]} : {self.deps[i]}")
            print(f"\t{self.cmds[i]}")
            print(f"\ttouch {self.tgts[i]}")

    def write(self):
        with open(self.make_file, "w") as f:
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


class Sample(object):
    def __init__(self):
        self.id = ""
        self.fastq1 = ""
        self.fastq2 = ""

    def __init__(self, id, fastq1, fastq2):
        self.id = id
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"id     : {self.id}")
        print(f"fastq1 : {self.fastq1}")
        print(f"fastq2 : {self.fastq2}")


if __name__ == "__main__":
    main()
