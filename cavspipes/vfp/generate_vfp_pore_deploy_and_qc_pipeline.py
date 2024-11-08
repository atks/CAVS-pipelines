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
import re
import sys
from shutil import copy2
from datetime import datetime

@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="pore_deploy_and_qc.mk",
    help="make file name",
)
@click.option("-r", "--run_id", required=True, help="Run ID")
@click.option("-i", "--nanopore_dir", required=True, help="nanopore directory")
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-q", "--qscore", help="Minimum passing read Q-score", required=False, default=7
)
@click.option(
    "-l", "--len", help="Minimum passing read length", required=False, default=20
)
@click.option("-x", "--memory", help="Memory for fastqc", required=False, default=2048)
@click.option("-y", "--model", help="Model for Dorado calling", required=False, default="dna_r10.4.1_e8.2_400bps_sup@v5.0.0")
@click.option("-k", "--kit", default="SQK-LSK114", show_default=True, help="Kit ID")
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(
    make_file,
    run_id,
    nanopore_dir,
    working_dir,
    qscore,
    len,
    memory,
    model,
    kit,
    sample_file,
):
    """
    Moves Oxford Nanopore Technology fastq files to a destination and performs QC

    e.g. generate_vfp_pore_deploy_and_qc_pipeline -r pore6 -i raw -s pore6.sa
    """
    dest_dir = working_dir + "/" + run_id
    nanopore_dir = os.path.abspath(nanopore_dir)

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<20} :   {1:<10}".format("nanopore_dir", nanopore_dir))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("minumum qscore", qscore))
    print("\t{0:<20} :   {1:<10}".format("minumum length", len))
    print("\t{0:<20} :   {1:<10}".format("memory", memory))
    print("\t{0:<20} :   {1:<10}".format("model", model))
    print("\t{0:<20} :   {1:<10}".format("kit", kit))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("dest_dir", dest_dir))

    # version
    version = "1.0.0"

    # programs
    dorado = "/usr/local/dorado-0.8.1/bin/dorado"
    #dna_r10.4.1_e8.2_400bps_hac@v4.4.0
    #dna_r10.4.1_e8.2_400bps_sup@v4.3.0
    #dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
    dorado_basecall_model = f"/usr/local/dorado-0.8.1/models/{model}"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    fastqc = f"/usr/local/FastQC-0.12.1/fastqc --adapters /usr/local/FastQC-0.12.1/Configuration/adapter_list.nanopore.txt --memory {memory}"
    multiqc = "/usr/local/bin/multiqc"
    nanoplot = "/usr/local/NanoPlot-1.43/bin/NanoPlot"
    ft = "/usr/local/cavstools-0.0.1/ft"

    run = Run(run_id)
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, barcode = line.rstrip().split("\t")
                run.add_sample(index, sample_id, barcode)

    # create directories in destination folder directory
    analysis_dir = f"{dest_dir}/analysis"
    log_dir = f"{working_dir}/log"
    aux_dir = f"{working_dir}/aux"
    bam_dir = f"{working_dir}/bam"
    fastq_dir = f"{working_dir}/fastq"
    trace_dir = f"{dest_dir}/trace"
    try:
        os.makedirs(analysis_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(aux_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
        for sample in run.samples:
            sample_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(sample_dir, exist_ok=True)
            os.makedirs(f"{sample_dir}/fastqc_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/nanoplot_result", exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # base call
    # dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.2.0  pod5s/ > calls.bam
    output_bam_file = f"{working_dir}/bam/basecalls.bam"
    pod5_files_dir = nanopore_dir
    log = f"{log_dir}/dorado_base_caller.log"
    err = f"{log_dir}/dorado_base_caller.err"
    tgt = f"{log_dir}/dorado_base_caller.OK"
    dep = ""
    cmd = f"{dorado} duplex -r {dorado_basecall_model} {pod5_files_dir} > {output_bam_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    #trim and filter and extract only simple single reads and duplex reads
    input_bam_file = f"{bam_dir}/basecalls.bam"
    output_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
    tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz.OK"
    dep = f"{log_dir}/dorado_base_caller.OK"
    cmd = f"{dorado} trim {input_bam_file} | {samtools} bam2fq -T \"*\" | grep -iP \"dx:i:0|dx:i:1\" -A 3 | grep -vP \"^--$$\"  | {ft} filter -q 7 -l 20 -o {output_fastq_file} "
    pg.add(tgt, dep, cmd)

    # symbolic link for fastqc
    fastqc_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"

    src_fastq = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
    dst_fastq = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}.fastq.gz"
    tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastq.gz.OK"
    dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz.OK"
    cmd = f"ln -sf {src_fastq} {dst_fastq}"
    pg.add(tgt, dep, cmd)

    # fastqc
    input_fastq_file = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}.fastq.gz"
    log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc.log"
    err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc.err"
    tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}_fastqc.OK"
    dep = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastq.gz.OK"
    cmd = f"{fastqc} {input_fastq_file} -o {fastqc_dir} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    #nanoplot
    input_fastq_file = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}.fastq.gz"
    nanoplot_dir =f"{analysis_dir}/{sample.idx}_{sample.id}/nanoplot_result"
    log = f"{log_dir}/{sample.idx}_{sample.id}_nanoplot.log"
    err = f"{log_dir}/{sample.idx}_{sample.id}_nanoplot.err"
    tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}_nanoplot.OK"
    dep = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastq.gz.OK"
    cmd = f"{nanoplot} --fastq {input_fastq_file} -o {nanoplot_dir} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # write make file
    print("Writing pipeline")
    pg.write()

    #copy files to trace
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


class Sample(object):
    def __init__(self, idx, id, barcode):
        self.idx = idx
        self.padded_idx = f"{idx:02}"
        self.id = id
        self.barcode = barcode

    def print(self):
        print(f"index   : {self.idx}")
        print(f"id      : {self.id}")
        print(f"barcode : {self.barcode}")


class Run(object):
    def __init__(self, id):
        m = re.match(r"\D+(\d+)", id)
        if m is not None:
            self.idx = m.group(1)
        else:
            self.idx = 999
        self.id = id
        self.samples = []

    def add_sample(self, idx, sample_id, barcode):
        self.samples.append(Sample(idx, sample_id, barcode))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id}")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")

if __name__ == "__main__":
    main()  # type: ignore