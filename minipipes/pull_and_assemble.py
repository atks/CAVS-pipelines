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
    mpm = MiniPipeManager(f"{output_dir}/pull_and_assemble.log")

    # programs
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    spades = "/usr/local/SPAdes-4.0.0/bin/spades.py"
    seqtk = "/usr/local/seqtk-1.4/seqtk"
    iteration = 0

    reference_fasta_file = os.path.abspath(reference_fasta_file)

    # make directories
    output_dir = os.path.abspath(output_dir)
    ref_dir = os.path.join(output_dir, "ref")
    fastq_dir = os.path.join(output_dir, "fastq")
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    ref_fasta_file_basename = os.path.basename(reference_fasta_file)
    last_ref_fasta_file = reference_fasta_file

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

    input_fastq_file1 = os.path.join(fastq_dir, "ilm.r1.fastq.gz")
    input_fastq_file2 = os.path.join(fastq_dir, "ilm.r2.fastq.gz")

    improved_assembly = True
    while improved_assembly:

        # make directories
        iteration_dir = os.path.join(output_dir, f"iter_{iteration}")
        ref_dir = os.path.join(iteration_dir, "ref")
        assembly_dir = f"{iteration_dir}/assembly"
        try:
            os.makedirs(iteration_dir, exist_ok=True)
            os.makedirs(ref_dir, exist_ok=True)
            os.makedirs(assembly_dir, exist_ok=True)
        except OSError as error:
            print(f"{error.filename} cannot be created")

        #  copy reference file to ref_dir
        output_fasta_file = os.path.join(ref_dir, ref_fasta_file_basename)
        cmd = f"cp {last_ref_fasta_file} {output_fasta_file}"
        tgt = f"{output_fasta_file}.OK"
        desc = f"Step {iteration}: Copy reference sequence to ref directory"
        mpm.run(cmd, tgt, desc)

        #  construct reference
        reference_fasta_file = os.path.join(ref_dir, ref_fasta_file_basename)
        log = os.path.join(ref_dir, "bwa_index.log")
        cmd = f"{bwa} index -a bwtsw {reference_fasta_file} 2> {log}"
        tgt = f"{reference_fasta_file}.bwa_index.OK"
        desc = f"Step {iteration}: Construct bwa reference"
        mpm.run(cmd, tgt, desc)

        #  align
        output_bam_file = os.path.join(iteration_dir, "ilm.bam")
        log = os.path.join(iteration_dir, "bwa_mem.log")
        cmd = f"{bwa} mem -t 2 -M {reference_fasta_file} {input_fastq_file1} {input_fastq_file2} 2> {log}| {samtools} view -hF4 | {samtools} sort -o {output_bam_file}"
        tgt = f"{output_bam_file}.OK"
        desc = f"Step {iteration}: Align to reference with bwa mem"
        mpm.run(cmd, tgt, desc)

        #  index
        input_bam_file = os.path.join(iteration_dir, "ilm.bam")
        cmd = f"{samtools} index {input_bam_file}"
        tgt = f"{input_bam_file}.bai.OK"
        desc = f"Step {iteration}: Index bam file"
        mpm.run(cmd, tgt, desc)

        #  coverage
        input_bam_file = os.path.join(iteration_dir, "ilm.bam")
        output_stats_file = os.path.join(iteration_dir, "ilm.stats.txt")
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        tgt = f"{output_stats_file}.OK"
        desc = f"Step {iteration}: Illumina coverage statistics"
        mpm.run(cmd, tgt, desc)

        #  consensus
        input_bam_file = os.path.join(iteration_dir, "ilm.bam")
        output_fasta_file = os.path.join(iteration_dir, "ilm.consensus.fasta")
        cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
        tgt = f"{output_fasta_file}.OK"
        desc = f"Step {iteration}: Illumina consensus"
        mpm.run(cmd, tgt, desc)

        # extract reads IDs that were aligned
        input_bam_file = os.path.join(iteration_dir, "ilm.bam")
        id_text_file = f"{iteration_dir}/aligned_read_IDs.txt"
        cmd = f"{samtools} view {input_bam_file} | cut -f1 | sort | uniq > {id_text_file}"
        tgt = f"{id_text_file}.OK"
        desc = f"Step {iteration}: Extract Aligned Illumina reads IDs"
        mpm.run(cmd, tgt, desc)

        # extract fastq paired reads
        aligned_fastq_file1 = os.path.join(iteration_dir, "aligned.r1.fastq.gz")
        cmd = f"{seqtk} subseq {input_fastq_file1} {id_text_file} | gzip > {aligned_fastq_file1 } "
        tgt = f"{aligned_fastq_file1}.OK"
        desc = f"Step {iteration}: Extract Read 1"
        mpm.run(cmd, tgt, desc)

        aligned_fastq_file2 = os.path.join(iteration_dir, "aligned.r2.fastq.gz")
        tgt = f"{ aligned_fastq_file2}.OK"
        cmd = f"{seqtk} subseq {input_fastq_file2} {id_text_file} | gzip > {aligned_fastq_file2 } "
        desc = f"Step {iteration}: Extract Read 2"
        mpm.run(cmd, tgt, desc)

        # assemble
        log = f"{assembly_dir}/assembly.log"
        err = f"{assembly_dir}/assembly.err"
        tgt = f"{assembly_dir}/assembly.OK"
        cmd = f"{spades} -1 {aligned_fastq_file1} -2 {aligned_fastq_file2} -o {assembly_dir} --threads 10 --isolate  > {log} 2> {err}"
        desc = f"Step {iteration}: Illumina assembly"
        mpm.run(cmd, tgt, desc)

        # update last reference file used
        last_ref_fasta_file = f"{assembly_dir}/contigs.fasta"

        iteration += 1

        if iteration == 3:
            improved_assembly = False

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
