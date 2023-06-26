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

import os
import click
import sys


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_ndv_consensus_analysis.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(make_file, working_dir, sample_file, reference_fasta_file):
    """
    NDV consensus and coverage calculations

    e.g. generate_ndv_coverage_consensus_analysis_pipeline
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # create directories in destination folder directory
    log_dir = f"{working_dir}/log"
    try:
        os.makedirs(log_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    align_and_consensus = "/home/atks/programs/cavs/var/align_and_consensus.py"

    # initialize
    pg = PipelineGenerator(make_file)

    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                id, ilm_fastq1, ilm_fastq2, ont_fastq = line.rstrip().split("\t")
                output_dir = f"{working_dir}/{id}"
                dep = ""
                tgt = f"{log_dir}/{id}.OK"
                log = f"{log_dir}/{id}.log"
                err = f"{log_dir}/{id}.err"
                cmd = ""
                if ont_fastq == "n/a":
                    cmd = f"{align_and_consensus} -s {id} -w {output_dir} -r {reference_fasta_file} -1 {ilm_fastq1} -2 {ilm_fastq2} > {log} 2> {err}"
                else:
                    cmd = f"{align_and_consensus} -s {id} -w {output_dir} -r {reference_fasta_file} -1 {ilm_fastq1} -2 {ilm_fastq2} -n {ont_fastq} > {log} 2> {err}"

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
        self.seq_tech = ""
        self.fastq1 = ""
        self.fastq2 = ""

    def __init__(self, id, seq_tech, fastq1, fastq2):
        self.id = id
        self.seq_tech = seq_tech
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"id         : {self.id}")
        print(f"seq tech   : {self.seq_tech}")
        print(f"fastq1     : {self.fastq1}")
        print(f"fastq2     : {self.fastq2}")


if __name__ == "__main__":
    main()
