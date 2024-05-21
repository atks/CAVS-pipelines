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
    default="run_ndv_coverage_analysis.mk",
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
    NDV coverage calculations

    e.g. generate_ndv_coverage_analysis_pipeline
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    samples = []
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                seq_tech = "ilm"
                id, fastq1, fastq2 = line.rstrip().split("\t")
                run_id, sample_id, rest = fastq1.split("_", 2)
                if fastq2 == "n/a":
                    seq_tech = "ont"
                    fastq2 = ""
                else:
                    fastq2 = f"/net/singapura/var/hts/{seq_tech}{run_id}/{fastq2}"
                fastq1 = f"/net/singapura/var/hts/{seq_tech}{run_id}/{fastq1}"

                samples.append(Sample(id, seq_tech, fastq1, fastq2))

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    log_dir = f"{working_dir}/log"
    try:
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(f"{working_dir}/data", exist_ok=True)
        os.makedirs(f"{working_dir}/bams", exist_ok=True)
        os.makedirs(f"{working_dir}/stats", exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    # data directory
    # NDV reference sequence NC_039223.1.fasta
    ref_ndv_fasta = f"{working_dir}/NC_039223.1.fasta"

    # mapping with minimap2
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.15/samtools"

    # bwa index -a bwtsw NC_039223.1.fa
    # minimap2 -d HQ878327.mmi HQ878327.fasta
    bwa_ref_fasta_file = f"{working_dir}/bwa_ref/NC_039223.1.fasta"
    minimap2_ref_mmi_file = f"{working_dir}/minimap2_ref/NC_039223.1.mmi"

    for sample in samples:
        # sample.print()
        if sample.seq_tech == "ont":
            output_bam_file = f"{working_dir}/bams/{sample.id}.minimap2.bam"
            log = f"{log_dir}/{sample.id}.minimap2.log"
            err = f"{log_dir}/{sample.id}.minimap2.err"
            dep = f""
            tgt = f"{log_dir}/{sample.id}.minimap2.bam.OK"
            cmd = f"{minimap2} -ax map-ont {minimap2_ref_mmi_file} {sample.fastq1} | {samtools} sort - | {samtools} view -o {output_bam_file} > {log} 2> {err}"
            pg.add(tgt, dep, cmd)

            input_bam_file = f"{working_dir}/bams/{sample.id}.minimap2.bam"
            output_txt_file = f"{working_dir}/stats/{sample.id}.minimap2.coverage.txt"
            dep = f"{log_dir}/{sample.id}.minimap2.bam.OK"
            tgt = f"{log_dir}/{sample.id}.minimap2.coverage.txt.OK"
            cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
            pg.add(tgt, dep, cmd)
        else:
            output_bam_file = f"{working_dir}/bams/{sample.id}.bwa.bam"
            log = f"{log_dir}/{sample.id}.bwa.log"
            err = f"{log_dir}/{sample.id}.bwa.err"
            dep = f""
            tgt = f"{log_dir}/{sample.id}.bwa.bam.OK"
            cmd = f"{bwa} mem -t 2 -M {bwa_ref_fasta_file} {sample.fastq1} {sample.fastq2} | {samtools} sort - | {samtools} view -o {output_bam_file} > {log} 2> {err}"
            pg.add(tgt, dep, cmd)

            input_bam_file = f"{working_dir}/bams/{sample.id}.bwa.bam"
            output_txt_file = f"{working_dir}/stats/{sample.id}.bwa.coverage.txt"
            dep = f"{log_dir}/{sample.id}.bwa.bam.OK"
            tgt = f"{log_dir}/{sample.id}.bwa.coverage.txt.OK"
            cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
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
