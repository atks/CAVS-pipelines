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


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_crocodylus_porosus_mito_analysis.mk",
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
    ND4, ATP6 and 16S gene mitochondria analysis

    e.g. generate_crocodylus_porosus_mito_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                sample_id, gene, fastq1, fastq2 = line.rstrip().split("\t")
                samples.append(Sample(sample_id, gene, fastq1, fastq2))

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    bam_dir = f"{working_dir}/bam"
    stats_dir = f"{working_dir}/stats"
    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    #################
    # reference files
    #################
    id = "AJ810453.1"
    output_genbank_file = f"{ref_dir}/{id}.genbank"
    dep = ""
    tgt = f"{output_genbank_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format genbank > {output_genbank_file}"
    pg.add(tgt, dep, cmd)

    output_fasta_file = f"{ref_dir}/{id}.fasta"
    dep = ""
    tgt = f"{output_fasta_file}.OK"
    cmd = f"efetch -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    #########
    # mapping
    #########
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.16/bin/samtools"
    seqtk = "/usr/local/seqtk-1.3/seqtk"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

    # create bwa indices
    input_fasta_file = f"{ref_dir}/AJ810453.1.fasta"
    output_bwt_file = f"{input_fasta_file}.bwt"
    log = f"{output_bwt_file}.log"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_bwt_file}.OK"
    cmd = f"{bwa} index -a bwtsw {input_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    for sample in samples:
        input_fastq_file1 = sample.fastq1
        input_fastq_file2 = sample.fastq2
        ref_fasta_file = f"{ref_dir}/AJ810453.1.fasta"
        output_bam_file = f"{bam_dir}/{sample.id}.bam"
        log = f"{output_bam_file}.log"
        dep = f"{output_bwt_file}.OK"
        tgt = f"{output_bam_file}.OK"
        cmd = f"set -o pipefail; {bwa} mem -t 2 -M {ref_fasta_file} {input_fastq_file1} {input_fastq_file2} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
        pg.add(tgt, dep, cmd)

        input_bam_file = output_bam_file
        dep = f"{input_bam_file}.OK"
        tgt = f"{input_bam_file}.bai.OK"
        cmd = f"{samtools} index {input_bam_file}"
        pg.add(tgt, dep, cmd)

    # clean
    pg.add_clean(f"rm -fr {ref_dir} {bam_dir} {stats_dir}")

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


class Sample(object):
    def __init__(self):
        self.id = ""
        self.gene = ""
        self.fastq1 = ""
        self.fastq2 = ""

    def __init__(self, id, gene, fastq1, fastq2):
        self.id = id
        self.gene = gene
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"id      : {self.id}")
        print(f"gene    : {self.gene}")
        print(f"fastq1  : {self.fastq1}")
        print(f"fastq2  : {self.fastq2}")


if __name__ == "__main__":
    main()
