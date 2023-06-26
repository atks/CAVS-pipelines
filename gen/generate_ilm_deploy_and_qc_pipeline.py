#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2021 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import argparse
import subprocess
import textwrap
import gzip
import re
import sys
from shutil import copy2
from datetime import datetime


def main():
    # programs
    fastqc = "/usr/local/FastQC-0.11.9/fastqc"
    multiqc = "/usr/local/bin/multiqc"

    # options
    parser = argparse.ArgumentParser(
        description="Moves illumina fastq files to a destination and performs QC",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
           usage: generate_ilm_deploy_and_qc_pipeline -m <make_file> -d <dest_dir>

           """
        ),
    )
    parser.add_argument(
        "-m", "--make_file", help="make file name", type=str, default="deploy_and_qc.mk"
    )
    parser.add_argument(
        "-d", "--dest_dir", help="destination directory", type=str, required=True
    )
    parser.add_argument("-s", "--sample", help="sample file", type=str, required=True)
    args = parser.parse_args()

    if not os.path.isabs(args.dest_dir):
        sys.exit("Please input an absolute path for destination directory")

    print("Options")
    print(f"destination: {args.dest_dir}")
    print(f"sample file: {args.sample}")

    # read sample file
    SAMPLES = []
    file = open(args.sample, "r")
    for line in file:
        if not line.startswith("#"):
            id, fastq1, fastq2 = line.rstrip().split("\t")
            SAMPLES.append(Sample(id, fastq1, fastq2))
    file.close()

    # create directories and copy files
    dest_dir = args.dest_dir
    ymd = datetime.today().strftime("%Y%m%d")
    print(f"YMD={ymd}")
    os.makedirs(dest_dir, exist_ok=True)
    print(dest_dir)

    for sample in SAMPLES:
        os.makedirs(f"{dest_dir}/{sample.id}", exist_ok=True)
        os.makedirs(f"{dest_dir}/{sample.id}/{ymd}_fastqc_results", exist_ok=True)
        copy2(sample.fastq1, f"{dest_dir}/{sample.id}")
        print(f"Copied {sample.fastq1}")
        copy2(sample.fastq2, f"{dest_dir}/{sample.id}")
        print(f"Copied {sample.fastq2}")
        # update paths
        sample.fastq1 = f"{dest_dir}/{sample.id}/{os.path.basename(sample.fastq1)}"
        sample.fastq2 = f"{dest_dir}/{sample.id}/{os.path.basename(sample.fastq2)}"
        # sample.print()

    # initialize
    pg = PipelineGenerator(args.make_file)
    multiqc_dep = ""

    for sample in SAMPLES:
        # fastqc
        fastqc_dir = f"{dest_dir}/{sample.id}/{ymd}_fastqc_results"

        log = f"{fastqc_dir}/fastqc1.log"
        err = f"{fastqc_dir}/fastqc1.err"
        tgt = f"{fastqc_dir}/fastqc1.OK"
        dep = ""
        cmd = f"{fastqc} {sample.fastq1} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        multiqc_dep += f" {tgt}"

        log = f"{fastqc_dir}/fastqc1.log"
        err = f"{fastqc_dir}/fastqc1.err"
        tgt = f"{fastqc_dir}/fastqc2.OK"
        dep = ""
        cmd = f"{fastqc} {sample.fastq2} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        multiqc_dep += f" {tgt}"

    # multiqc
    multiqc_dir = f"{dest_dir}/{ymd}_multiqc_results"
    os.makedirs(f"{multiqc_dir}", exist_ok=True)

    log = f"{multiqc_dir}/multiqc.log"
    err = f"{multiqc_dir}/multiqc.err"
    tgt = f"{multiqc_dir}/multiqc.OK"
    dep = multiqc_dep
    cmd = f"{multiqc} {dest_dir} -o {multiqc_dir} > {log} 2> {err}"
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

    def add(self, tgt, dep, cmd):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(cmd)

    def print(self):
        print(".DELETE_ON_ERROR:")
        for i in range(len(self.tgts)):
            print(f"{self.tgts[i]} : {self.deps[i]}")
            print(f"\t{self.cmds[i]}")
            print(f"\ttouch {self.tgts[i]}")

    def write(self):
        f = open(self.make_file, "w")
        f.write(".DELETE_ON_ERROR:\n\n")
        f.write("all : ")
        for i in range(len(self.tgts)):
            f.write(f"{self.tgts[i]} ")
        f.write("\n\n")

        for i in range(len(self.tgts)):
            f.write(f"{self.tgts[i]} : {self.deps[i]}\n")
            f.write(f"\t{self.cmds[i]}\n")
            f.write(f"\ttouch {self.tgts[i]}\n\n")
        f.close()


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


main()
