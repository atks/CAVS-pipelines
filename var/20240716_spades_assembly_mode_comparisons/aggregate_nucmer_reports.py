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
from shutil import copy2


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_assembly_comparisons_analysis.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(make_file, working_dir, sample_file):
    """
    Generates assembly comparisons statistics.

    e.g. generate_assembly_comparisons_pipeline.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # create directories in destination folder directory
    pairwise_alignments_dir = f"{working_dir}/pairwise_alignments"
    trace_dir = f"{working_dir}/trace"
    try:
        os.makedirs(pairwise_alignments_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    nucmer = "/usr/local/mummer-4.0.0/bin/nucmer"
    dnadiff = "/usr/local/mummer-4.0.0/bin/dnadiff"

    # initialize
    pg = PipelineGenerator(make_file)

    samples = []
    idx = 0
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                idx += 1
                sample_id, contigs_file = line.rstrip().split("\t")
                samples.append(Sample(idx, sample_id, contigs_file))

    samtools_multiqc_dep = ""

    #aggregate statistics from report
    # NUCMER
    #
    #                             [REF]                [QRY]
    # [Sequences]
    # TotalSeqs                         51               264808
    # AlignedSeqs             48(94.1176%)         490(0.1850%)
    # UnalignedSeqs             3(5.8824%)     264318(99.8150%)
    #
    # [Bases]
    # TotalBases                    148294             99993229
    # AlignedBases        138783(93.5864%)      208972(0.2090%)
    # UnalignedBases         9511(6.4136%)   99784257(99.7910%)
    #
    # [Alignments]
    # 1-to-1                           138                  138
    # TotalLength                   134376               134364
    # AvgLength                   973.7391             973.6522
    # AvgIdentity                  99.6538              99.6538
    #
    # M-to-M                           556                  556
    # TotalLength                   227699               227687
    # AvgLength                   409.5306             409.5090
    # AvgIdentity                  95.9291              95.9291



    pg.add_clean(f"rm -rf {pairwise_alignments_dir}")

    # write make file
    print("Writing pipeline")
    pg.write()

    # copy files to trace
    copy2(__file__, trace_dir)
    copy2(make_file, trace_dir)
    copy2(sample_file, trace_dir)

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
    def __init__(self, idx, id, contigs_file):
        self.idx = idx
        self.id = id
        self.contigs_file = contigs_file

    def print(self):
        print(f"idx           : {self.idx}")
        print(f"id            : {self.id}")
        print(f"contigs_file  : {self.contigs_file}")

if __name__ == "__main__":
    main()  # type: ignore
