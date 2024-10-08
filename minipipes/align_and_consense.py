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

import sys
import os
import click
import subprocess
from shutil import copy2

@click.command()
@click.option(
    "-o",
    "--output_dir",
    default=os.path.join(os.getcwd(), "align_consensus"),
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
    default="consensus",
    show_default=True,
    help="consensus sequence sample ID",
)
@click.option(
    "-m",
    "--ont_model",
    required=False,
    #r941_min_hac_g507
    #r1041_e82_400bps_sup_v5.0.0
    default="r1041_e82_400bps_sup_v4.3.0",
    help="ONT Machine Model",
)
@click.option(
    "-n",
    "--input_ont_fastq_files",
    required=False,
    default="",
    help="ONT fastq files",
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
    input_ont_fastq_files,
    input_ilm_read1_fastq_files,
    input_ilm_read2_fastq_files,
    ont_model,
    sample_id,
    output_dir,
    reference_fasta_file,
):
    """
    Aligns all fastq files to a reference sequence file and generates a consensus sequence.

    e.g. align_and_consensus -1 r1.fastq.gz -2 r2.fastq.gz -n ont.fastq.qz  -r ref.fasta
    """
    illumina = len(input_ilm_read1_fastq_files) != 0
    nanopore = len(input_ont_fastq_files) != 0

    print("\t{0:<20} :   {1:<10}".format("sample ID", sample_id))
    print("\t{0:<20} :   {1:<10}".format("output directory", output_dir))
    print("\t{0:<20} :   {1:<10}".format("reference fasta file", reference_fasta_file))
    if illumina:
        print("\tIllumina reads")
        print("\t{0:<20} :   {1:<10}".format("fastq1", input_ilm_read1_fastq_files))
        print("\t{0:<20} :   {1:<10}".format("fastq2", input_ilm_read1_fastq_files))
    if nanopore:
        print("\tNanopore reads")
        print("\t{0:<20} :   {1:<10}".format("ont fastq", input_ont_fastq_files))

    # version
    version = "1.1.0"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/align_and_consense.log")

    # programs
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    #"docker run -t -v  `pwd`:`pwd` -w `pwd` fischuu/quast quast.py "
    activate_medaka_virtualenv = "source /usr/local/medaka-1.12.1/bin/activate"
    medaka = "/usr/local/medaka-1.12.1/bin/medaka"
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

    # process illumina reads
    if illumina:
        print("Processing illumina reads")

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

    # process nanopore reads
    if nanopore:
        print("Processing nanopore reads")

        #  combine fastq files
        fastq_files = " ".join(input_ont_fastq_files.split(","))
        output_fastq_file = os.path.join(fastq_dir, "ont.fastq.gz")
        cmd = f"zcat {fastq_files} | gzip > {output_fastq_file}"
        tgt = f"{output_fastq_file}.OK"
        desc = f"Combining nanopore fastq files"
        mpm.run(cmd, tgt, desc)

        #  construct reference
        cmd = f"{minimap2} -d {reference_fasta_file}.mmi {reference_fasta_file}"
        tgt = f"{reference_fasta_file}.minimap2_index.OK"
        desc = f"Construct minimap2 reference"
        mpm.run(cmd, tgt, desc)

        #  align
        minimap2_ref_mmi_file = f"{reference_fasta_file}.mmi"
        input_fastq_file = os.path.join(fastq_dir, "ont.fastq.gz")
        output_bam_file = os.path.join(bam_dir, "ont.bam")
        cmd = f"{minimap2} -ax map-ont {minimap2_ref_mmi_file} {input_fastq_file}  | {samtools} view -hF4 | {samtools} sort -o {output_bam_file}"
        tgt = f"{output_bam_file}.OK"
        desc = f"Align to reference with minimap2"
        mpm.run(cmd, tgt, desc)

        #  index
        input_bam_file = os.path.join(bam_dir, "ont.bam")
        cmd = f"{samtools} index {input_bam_file}"
        tgt = f"{input_bam_file}.bai.OK"
        desc = f"Index bam file"
        mpm.run(cmd, tgt, desc)

        #  coverage
        input_bam_file = os.path.join(bam_dir, "ont.bam")
        output_stats_file = os.path.join(stats_dir, "ont.stats.txt")
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        tgt = f"{output_stats_file}.OK"
        desc = f"Nanopore coverage statistics"
        mpm.run(cmd, tgt, desc)

        #  consensus
        input_bam_file = os.path.join(bam_dir, "ont.bam")
        output_hdf_file = os.path.join(bam_dir, "ont.hdf")
        cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model {ont_model}"
        tgt = f"{output_hdf_file}.OK"
        desc = f"Oxford Nanopore consensus contigs"
        mpm.run(cmd, tgt, desc)

        #  stitch
        input_hdf_file = os.path.join(bam_dir, "ont.hdf")
        output_fasta_file = os.path.join(fasta_dir, "ont.consensus.fasta")
        cmd = f"{medaka} stitch --fill_char N {input_hdf_file} {reference_fasta_file} {output_fasta_file} "
        tgt = f"{output_fasta_file}.OK"
        desc = f"Oxford Nanopore consensus assembly/scaffold"
        mpm.run(cmd, tgt, desc)

    # combine illumina and nanopore reads
    if illumina and nanopore:
        print("Combine illumina and nanopore reads")

        # merge bams
        input_bam_file1 = os.path.join(bam_dir, "ilm.bam")
        input_bam_file2 = os.path.join(bam_dir, "ont.bam")
        output_bam_file = os.path.join(bam_dir, "ilm_ont.bam")
        cmd = (
            f"{samtools} merge -f {output_bam_file} {input_bam_file1} {input_bam_file2}"
        )
        tgt = f"{output_bam_file}.OK"
        desc = f"Merge Illumina and Nanopore alignments"
        mpm.run(cmd, tgt, desc)

        # index bam
        input_bam_file = os.path.join(bam_dir, "ilm_ont.bam")
        cmd = f"{samtools} index {input_bam_file}"
        tgt = f"{output_bam_file}.bai.OK"
        desc = f"Index merged alignments"
        mpm.run(cmd, tgt, desc)

        # coverage
        input_bam_file = os.path.join(bam_dir, "ilm_ont.bam")
        output_stats_file = os.path.join(stats_dir, "ilm_ont.stats.txt")
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        tgt = f"{output_stats_file}.OK"
        desc = f"Compute coverage of merged alignments"
        mpm.run(cmd, tgt, desc)

        # consensus
        input_bam_file = os.path.join(bam_dir, "ilm_ont.bam")
        output_fasta_file = os.path.join(fasta_dir, "ilm_ont.consensus.fasta")
        cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
        tgt = f"{output_fasta_file}.OK"
        desc = f"Illumina and Nanopore consensus"
        mpm.run(cmd, tgt, desc)

    # copy out consensus
    input_fasta_file = ""
    if illumina and nanopore:
        input_fasta_file = os.path.join(fasta_dir, "ilm_ont.consensus.fasta")
    elif illumina:
        input_fasta_file = os.path.join(fasta_dir, "ilm.consensus.fasta")
    else:
        input_fasta_file = os.path.join(fasta_dir, "ont.consensus.fasta")

    output_fasta_file = os.path.join(consensus_dir, f"{sample_id}.fasta")
    cmd = f"{seqkit} replace  -p '^.+$' -r {sample_id} {input_fasta_file} -o {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Final consensus"
    mpm.run(cmd, tgt, desc)

    # copy files to trace
    copy2(__file__, trace_dir)

    # write log file
    mpm.print_log()

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
