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
import re

@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_salmonella_serotyping.mk",
    help="make file name",
)
@click.option(
    "-o",
    "--output_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(make_file, output_dir, sample_file):
    """
    Generates assemblies and run raw reads on
        a. seroseq
        b. seroseq2
        c. SISTR
        d. MOST

    e.g. vm_generate_salmonella_serotyping_pipeline.py -s illu13.sa
    """
    log_dir = f"{output_dir}/log"

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("output_dir", output_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, fastq1, fastq2 = line.rstrip().split("\t")
                samples.append(
                    Sample(index, sample_id, fastq1, fastq2)
                )

    # create directories in destination folder directory
    new_dir = ""
    try:
        os.makedirs(log_dir, exist_ok=True)

        for sample in samples:
            new_dir = f"{output_dir}/{sample.id}"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{output_dir}/{sample.id}/spades_assembly"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{output_dir}/{sample.id}/seqsero2"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{output_dir}/{sample.id}/mlst"
            os.makedirs(new_dir, exist_ok=True)

    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # programs
    spades = "/usr/local/SPAdes-3.15.5/bin/spades.py"
    seqsero2 = "/usr/local/SeqSero2/bin/SeqSero2_package.py"
    mlst = "/usr/local/mlst-2.23.0/bin/mlst"

    # analyze
    for idx, sample in enumerate(samples):

        # spades
        out_dir = f"{output_dir}/{sample.id}/spades_assembly"
        log = f"{out_dir}/run.log"
        err = f"{out_dir}/run.err"
        tgt = f"{log_dir}/{sample.id}_spades_assembly_contigs.OK"
        dep = ""
        cmd = f'{spades} -o {out_dir} --isolate -1 {sample.fastq1} -2 {sample.fastq2} > {log} 2> {err}'
        pg.add(tgt, dep, cmd)

        # serovar and antigen
        # /usr/local/SeqSero2/bin/SeqSero2_package.py -d efg -n 23_1704 -t 2 -i /net/singapura/vm/hts/illu13/13_1_salmonella_23_1704_R1.fastq.gz  /net/singapura/vm/hts/illu13/13_1_salmonella_23_1704_R2.fastq.gz
        out_dir = f"{output_dir}/{sample.id}/seqsero2"
        log = f"{out_dir}/run.log"
        err = f"{out_dir}/run.err"
        tgt = f"{log_dir}/{sample.id}.seqsero2.OK"
        dep = ""
        cmd = f'{seqsero2} -d {out_dir} -n {sample.id} -t 2 -i {sample.fastq1} {sample.fastq2} > {log} 2> {err}'
        pg.add(tgt, dep, cmd)

        # sequence typing
        # mlst illu13/*/spades_assembly/contigs.fasta --json out.json --
        out_dir = f"{output_dir}/{sample.id}/mlst"
        input_contig_fasta_file = f"{output_dir}/{sample.id}/spades_assembly/contigs.fasta"
        log = f"{out_dir}/run.log"
        err = f"{out_dir}/run.err"
        tgt = f"{log_dir}/{sample.id}.mlst.OK"
        dep = f"{log_dir}/{sample.id}_spades_assembly_contigs.OK"
        cmd = f'{mlst} {input_contig_fasta_file} --json typing.json --scheme senterica_achtman_2  --nopath > {log} 2> {err}'
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
    def __init__(self): # type: ignore
        self.idx = ""
        self.id = ""
        self.fastq1 = ""
        self.fastq2 = ""

    def __init__(self, idx, id, fastq1, fastq2):
        self.idx = idx
        self.id = id
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"index   : {self.idx}")
        print(f"id      : {self.id}")
        print(f"fastq1  : {self.fastq1}")
        print(f"fastq2 : {self.fastq2}")


if __name__ == "__main__":
    main() # type: ignore