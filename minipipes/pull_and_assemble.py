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

import sys
import os
import click
import subprocess


@click.command()
@click.option(
    "-o",
    "--output_dir",
    default=os.path.join(os.getcwd(), "pull_assemble"),
    show_default=True,
    help="output directory",
)
@click.option(
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
@click.option(
    "-s",
    "--sample_id",
    required=True,
    default="assembly",
    show_default=True,
    help="assembled sequence sample ID",
)
@click.option(
    "-1",
    "--input_ilm_read1_fastq_files",
    required=False,
    default="",
    help="Illumina read 1 fastq files",
)
@click.option(
    "-2",
    "--input_ilm_read2_fastq_files",
    required=False,
    default="",
    help="Illumina read 2 fastq files",
)
def main(
    input_ilm_read1_fastq_files,
    input_ilm_read2_fastq_files,
    sample_id,
    output_dir,
    reference_fasta_file,
):
    """
    Aligns all fastq files to a reference sequence file, perform assembly on reads.  Perform this iteratively.

    e.g. pull_and_assemble -1 r1.fastq.gz -2 r2.fastq.gz -r ref.fasta
    """
    illumina = len(input_ilm_read1_fastq_files) != 0

    # version
    version = "1.0.0"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/align_and_consense.log")

    # programs
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    spades = "/usr/local/Spades-4.0.0/bin/spades.py"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

    # make directories
    output_dir = os.path.abspath(output_dir)
    fastq_dir = os.path.join(output_dir, "fastq")
    bam_dir = os.path.join(output_dir, "bam")
    fasta_dir = os.path.join(output_dir, "fasta")
    stats_dir = os.path.join(output_dir, "stats")
    consensus_dir = os.path.join(output_dir, "consensus")
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(consensus_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    #  combine fastq files
    read1_files = input_ilm_read1_fastq_files.split(",")
    output_read1_fastq_file = os.path.join(fastq_dir, "ilm.r1.fastq.gz")
    cmd = f"zcat {' '.join(read1_files)} | gzip > {output_read1_fastq_file}"
    tgt = f"{output_read1_fastq_file}.OK"
    desc = f"Combining read 1 fastq files"
    mpm.run(cmd, tgt, desc)

    #  combine fastq files
    read2_files = " ".join(input_ilm_read2_fastq_files.split(","))
    output_read2_fastq_file = os.path.join(fastq_dir, "ilm.r2.fastq.gz")
    cmd = f"zcat {read2_files} | gzip > {output_read2_fastq_file}"
    tgt = f"{output_read2_fastq_file}.OK"
    desc = f"Combining read 2 fastq files"
    mpm.run(cmd, tgt, desc)

    #  construct reference
    cmd = f"{bwa} index -a bwtsw {reference_fasta_file}"
    tgt = f"{reference_fasta_file}.bwa_index.OK"
    desc = f"Construct bwa reference"
    mpm.run(cmd, tgt, desc)

    #  align
    output_bam_file = os.path.join(bam_dir, "ilm.bam")
    log = os.path.join(bam_dir, "bwa_mem.log")
    cmd = f"{bwa} mem -t 2 -M {reference_fasta_file} {output_read1_fastq_file} {output_read2_fastq_file} 2> {log}| {samtools} view -hF4 | {samtools} sort -o {output_bam_file}"
    tgt = f"{output_bam_file}.OK"
    desc = f"Align to reference with bwa mem"
    mpm.run(cmd, tgt, desc)

    #  index
    input_bam_file = os.path.join(bam_dir, "ilm.bam")
    cmd = f"{samtools} index {input_bam_file}"
    tgt = f"{input_bam_file}.bai.OK"
    desc = f"Index bam file"
    mpm.run(cmd, tgt, desc)

    #  coverage
    input_bam_file = os.path.join(bam_dir, "ilm.bam")
    output_stats_file = os.path.join(stats_dir, "ilm.stats.txt")
    cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
    tgt = f"{output_stats_file}.OK"
    desc = f"Illumina coverage statistics"
    mpm.run(cmd, tgt, desc)

    #  consensus
    input_bam_file = os.path.join(bam_dir, "ilm.bam")
    output_fasta_file = os.path.join(fasta_dir, "ilm.consensus.fasta")
    cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Illumina consensus"
    mpm.run(cmd, tgt, desc)




class MiniPipeManager(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log_msg = []

    def run(self, cmd, tgt, desc):
        try:
            if os.path.exists(tgt):
                self.log(f"{desc} -  already executed")
                self.log(cmd)
                return
            else:
                self.log(f"{desc}")
                self.log(cmd)
                subprocess.run(cmd, shell=True, check=True)
                subprocess.run(f"touch {tgt}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            self.log(f" - failed")
            exit(1)

    def log(self, msg):
        print(msg)
        self.log_msg.append(msg)

    def print_log(self):
        self.log(f"\nlogs written to {self.log_file}")
        with open(self.log_file, "w") as f:
            f.write("\n".join(self.log_msg))


if __name__ == "__main__":
    main() # type: ignore
