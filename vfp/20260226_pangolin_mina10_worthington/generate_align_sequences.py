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
    default = "generate_align_worthington2024_sequences.mk",
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
@click.option("-g", "--genome_reference_fasta_file", required=True, help="genome reference FASTA file")
def main(make_file, working_dir, sample_file, genome_reference_fasta_file):
    """
    Moves fastq files to a destination and performs QC

    e.g. generate_align_sequences.py -s worthington2024.sa -g /usr/local/ref/vfp/ManJav1.0_HiC.fasta
    """
    
    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, run = line.rstrip().split("\t")
                fastq1 = f"/net/singapura/vfp/datasets/worthington2024/fastq/{run}_1.fastq.gz"
                fastq2 = f"/net/singapura/vfp/datasets/worthington2024/fastq/{run_id}_2.fastq.gz"
                samples.append(Sample(sample_id, fastq1, fastq2))


    if make_file == "":
        make_file = f"{working_dir}/{run_id}_novogene_wgs_deploy_and_qc.mk"

    print("\t{0:<28} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<28} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<28} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<28} :   {1:<10}".format("genome reference FASTA file", genome_reference_fasta_file))
    
    # create directories in destination folder directory
    log_dir = f"{working_dir}/log"
    ref_dir = f"{working_dir}/ref"
    bam_dir = f"{working_dir}/bam"
    stats_dir = f"{working_dir}/analysis"
    try:
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(f"{stats_dir}/coverage", exist_ok=True)
        os.makedirs(f"{stats_dir}/general_stats", exist_ok=True)
        os.makedirs(f"{stats_dir}/flag_stats", exist_ok=True)
        os.makedirs(f"{stats_dir}/idx_stats", exist_ok=True)
        os.makedirs(f"{stats_dir}/all", exist_ok=True)
            
    except OSError as error:
        print(f"{error.filename} cannot be created")

    #version
    version = "1.0.0"

    #programs
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    plot_bamstats = "/usr/local/samtools-1.17/bin/plot-bamstats"
    multiqc = "docker run  -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc "
    seqkit = "/usr/local/seqkit-2.10.1/seqkit"
        
    # initialize
    pg = PipelineGenerator(make_file)

    samtools_multiqc_dep = ""
    
    # copy reference
    dst_genome_reference_fasta_file = f"{ref_dir}/genome_reference.fasta"
    dep = ""
    tgt = f"{log_dir}/genome_reference.fasta.OK"
    cmd = f"ln -s {genome_reference_fasta_file} {dst_genome_reference_fasta_file}"
    pg.add(tgt, dep, cmd)

    # construct reference
    log = f"{log_dir}/genome_reference.bwa_index.log"
    dep = f""
    tgt = f"{log_dir}/genome_reference.bwa_index.OK"
    cmd = f"{bwa} index -a bwtsw {dst_genome_reference_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    for sample in samples:
        
        #  align
        output_bam_file = f"{bam_dir}/{sample.id}.bam"
        log = f"{log_dir}/{sample.id}.align.log"
        sort_log = f"{log_dir}/{sample.id}.align.sort.log"
        dep = f"{log_dir}/genome_reference.bwa_index.OK"
        tgt = f"{log_dir}/{sample.id}.bam.OK"
        cmd = f"{bwa} mem -t 2 -M {dst_genome_reference_fasta_file} {sample.fastq1} {sample.fastq2} 2> {log} | {samtools} view -h | {samtools} sort -o {output_bam_file} 2> {sort_log}"
        pg.add(tgt, dep, cmd)

        #  index
        input_bam_file = f"{bam_dir}/{sample.id}.bam"
        dep = f"{log_dir}/{sample.id}.bam.OK"
        tgt = f"{log_dir}/{sample.id}.bam.bai.OK"
        cmd = f"{samtools} index {input_bam_file}"
        pg.add(tgt, dep, cmd)

        # coverage
        input_bam_file = f"{bam_dir}/{sample.id}.bam"
        output_stats_file = f"{stats_dir}//coverage/{sample.id}.coverage.txt"
        dep = f"{log_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{log_dir}/{sample.id}.coverage.stats.OK"
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # stats
        output_stats_file = f"{stats_dir}/general_stats/{sample.id}.txt"
        dep = f"{log_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{log_dir}/{sample.id}.stats.OK"
        cmd = f"{samtools} stats {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # flag stats
        output_stats_file = f"{stats_dir}/flag_stats/{sample.id}.txt"
        dep = f"{log_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{log_dir}/{sample.id}.flag.stats.OK"
        cmd = f"{samtools} flagstat {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        #  idx stats
        output_stats_file = f"{stats_dir}/idx_stats/{sample.id}.txt"
        dep = f"{log_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{log_dir}/{sample.id}.idx.stats.OK"
        cmd = f"{samtools} idxstats {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # plot samtools stats
        input_stats_file = f"{stats_dir}/general_stats/{sample.id}.txt"
        dep = f"{log_dir}/{sample.id}.stats.OK"
        tgt = f"{log_dir}/{sample.id}.plot_bamstats.OK"
        cmd = f"{plot_bamstats} -p  {stats_dir}/plot_bamstats/plot {input_stats_file}"
        pg.add(tgt, dep, cmd)
        
    # plot samtools multiqc results
    analysis = "samtools"
    output_dir = f"{stats_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = samtools_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {stats_dir}; {multiqc} . -m samtools -f -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
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

    def add_srun(self, tgt, dep, cmd, cpu):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(f"srun --mincpus {cpu} {cmd}")

    def add_srun_blastdb(self, tgt, dep, cmd, cpu):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(f"srun --mincpus {cpu} --export=ALL,BLASTDB=/db/blast/nt {cmd}")

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

    def __init__(self, id, fastq1, fastq2):
        self.id = id
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"id             : {self.id}")
        print(f"fastq1         : {self.fastq1}")
        print(f"fastq2         : {self.fastq2}")
        



if __name__ == "__main__":
    main() # type: ignore
