#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2026 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import gzip
import sys
import click
import re
import subprocess


@click.command()
@click.option(
    "-m",
    "--make_file",
    default="download_sra_fastq.mk",
    show_default=True,
    help="make file name",
)
@click.option(
    "-s", "--sequence_run_id_file", required=True, help="list of sequence run IDs"
)
@click.option(
    "-o", "--output_dir", default=os.getcwd(), show_default=True, help="output fastq directory",
)
def main(make_file, sequence_run_id_file, output_dir):
    """
    generate_download_sra_fastq_pipeline -d viral

    e.g.  generate_download_sra_fastq_pipeline -s  -m download_gb.mk -o /home/atks/downloads

    """

    ids = []
    with open(sequence_run_id_file, "r") as file:
        for line in file:
            id = line.rstrip()
            ids.append(id)
   
    output_dir = os.path.abspath(output_dir)

    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("sequence run id file", sequence_run_id_file))
    print("\t{0:<20} :   {1:<10}".format("output directory", output_dir))
    print("\n")
    print(
        "Please invoked pipeline with not more than 8 jobs running concurrently due to NCBI restrictions"
    )
    print("i.e. make -f download_sra_fastq.mk -j 8 -k")

    #programs
    fastq_dump = "/usr/local/sratoolkit.3.3.0/bin/fastq-dump"


    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {output_dir} cannot be created")

    # generate make file
    print("Generating pipeline")
    pg = PipelineGenerator(make_file)

    # download each file
    for id in ids:
        err = f"{output_dir}/{id}.err"
        tgt = f"{output_dir}/{id}.OK"
        dep = ""
        cmd = f"{fastq_dump} {id} --split-files --gzip -O {output_dir} 2> {err}"
        pg.add(tgt, dep, cmd)

    # clean files
    cmd = f"rm {output_dir}/*.gz  {output_dir}/*.OK {output_dir}/*.err"
    pg.add_clean(cmd)

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


if __name__ == "__main__":
    main() # type: ignore
