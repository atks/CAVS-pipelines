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
    help="reference fasta file",
)
@click.option(
    "-s",
    "--sample_id",
    required=False,
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
    default="r1041_e82_400bps_hac_v4.3.0",
    show_default=True,
    help="ONT Machine Model e.g. r941_min_hac_g507, r1041_e82_400bps_sup_v5.0.0",
)
@click.option(
    "--in_situ_ref",
    is_flag=True,
    default=False,
    show_default=True,
    help="Build reference databases in the location of the supplied reference FASTA file",
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
    in_situ_ref,
    ont_model,
    sample_id,
    output_dir,
    reference_fasta_file,
):
    """
    Aligns all fastq files to a reference sequence file and generates a consensus sequence.

    e.g. align_and_consense.py -1 r1.fastq.gz -2 r2.fastq.gz -n ont.fastq.qz  -r ref.fasta
    """
    ex_situ_ref = not in_situ_ref

    illumina = len(input_ilm_read1_fastq_files) != 0
    nanopore = len(input_ont_fastq_files) != 0

    reference_fasta_file = os.path.abspath(reference_fasta_file)

    print("\t{0:<20} :   {1:<10}".format("sample ID", sample_id))
    print("\t{0:<20} :   {1:<10}".format("output directory", output_dir))
    print("\t{0:<20} :   {1:<10}".format("reference fasta file", reference_fasta_file))
    print("\t{0:<20} :   {1}".format("in situ ref db", in_situ_ref))
    if illumina:
        print("\tIllumina reads")
        print("\t{0:<20} :   {1:<10}".format("fastq1", input_ilm_read1_fastq_files))
        print("\t{0:<20} :   {1:<10}".format("fastq2", input_ilm_read2_fastq_files))
    if nanopore:
        print("\tNanopore reads")
        print("\t{0:<20} :   {1:<10}".format("ont fastq", input_ont_fastq_files))
        print("\t{0:<20} :   {1:<10}".format("ONT model", ont_model))

    # version
    version = "1.1.1"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/align_and_consense.log")

    # programs
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    #ontresearch/medaka:sha3486abaab0d3b90351617eb8622acf2028edb154 medaka v1.12.0
    #ontresearch/medaka:sha447c70a639b8bcf17dc49b51e74dfcde6474837b medaka v2.0.0
    medaka = "docker run  -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` ontresearch/medaka:sha3486abaab0d3b90351617eb8622acf2028edb154 medaka"
    seqkit = "/usr/local/seqkit-2.10.1/bin/seqkit"

    # make directories
    output_dir = os.path.abspath(output_dir)
    fastq_dir = os.path.join(output_dir, "fastq")
    bam_dir = os.path.join(output_dir, "bam")
    fasta_dir = os.path.join(output_dir, "fasta")
    stats_dir = os.path.join(output_dir, "stats")
    consensus_dir = os.path.join(output_dir, "consensus")
    trace_dir = f"{output_dir}/trace"
    ref_dir = os.path.abspath(os.path.dirname(reference_fasta_file))
    if ex_situ_ref:
        ref_dir = f"{output_dir}/ref"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(consensus_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
        if ex_situ_ref:
            os.makedirs(ref_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    if ex_situ_ref:
        copy2(reference_fasta_file, ref_dir)
        reference_fasta_file = os.path.join(ref_dir, os.path.basename(reference_fasta_file))

    # process illumina reads
    if illumina:
        print("Processing illumina reads")

        read1_files = input_ilm_read1_fastq_files.split(",")
        output_read1_fastq_file = os.path.join(fastq_dir, "ilm.r1.fastq.gz")
        #  combine fastq files
        if len(read1_files)>1:
            cmd = f"zcat {' '.join(read1_files)} | gzip > {output_read1_fastq_file}"
            tgt = f"{output_read1_fastq_file}.OK"
            desc = f"Combining read 1 fastq files"
            mpm.run(cmd, tgt, desc)
        else:
            cmd = f"cp {read1_files[0]} {fastq_dir}/ilm.r1.fastq.gz"
            tgt = f"{output_read1_fastq_file}.OK"
            desc = f"Copying read 1 fastq file"
            mpm.run(cmd, tgt, desc)

        read2_files = input_ilm_read2_fastq_files.split(",")
        output_read2_fastq_file = os.path.join(fastq_dir, "ilm.r2.fastq.gz")
        #  combine fastq files
        if len(read2_files)>1:
            cmd = f"zcat {' '.join(read2_files)} | gzip > {output_read2_fastq_file}"
            tgt = f"{output_read2_fastq_file}.OK"
            desc = f"Combining read 2 fastq files"
            mpm.run(cmd, tgt, desc)
        else:
            cmd = f"cp {read2_files[0]} {fastq_dir}/ilm.r2.fastq.gz"
            tgt = f"{output_read2_fastq_file}.OK"
            desc = f"Copying read 2 fastq file"
            mpm.run(cmd, tgt, desc)

        #  construct reference
        log = os.path.join(ref_dir, f"{reference_fasta_file}.bwa_index.log")
        cmd = f"{bwa} index -a bwtsw {reference_fasta_file} 2> {log}"
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
        fastq_files = input_ont_fastq_files.split(",")
        output_fastq_file = os.path.join(fastq_dir, "ont.fastq.gz")
        if len(fastq_files)>1:
            cmd = f"zcat {' '.join(fastq_files)} | gzip > {output_fastq_file}"
            tgt = f"{output_fastq_file}.OK"
            desc = f"Combining nanopore fastq files"
            mpm.run(cmd, tgt, desc)
        else:
            cmd = f"cp {fastq_files[0]} {fastq_dir}/ont.fastq.gz"
            tgt = f"{output_fastq_file}.OK"
            desc = f"Copying nanopore fastq file"
            mpm.run(cmd, tgt, desc)

        #  construct reference
        log = os.path.join(ref_dir, f"{reference_fasta_file}.minimap2_index.log")
        cmd = f"{minimap2} -d {reference_fasta_file}.mmi {reference_fasta_file} 2> {log}"
        tgt = f"{reference_fasta_file}.minimap2_index.OK"
        desc = f"Construct minimap2 reference"
        mpm.run(cmd, tgt, desc)

        #  align
        minimap2_ref_mmi_file = f"{reference_fasta_file}.mmi"
        minimap2_log_file = os.path.join(bam_dir, "minimap2.log")
        input_fastq_file = os.path.join(fastq_dir, "ont.fastq.gz")
        output_bam_file = os.path.join(bam_dir, "ont.bam")
        cmd = f"{minimap2} -ax map-ont {minimap2_ref_mmi_file} {input_fastq_file} 2> {minimap2_log_file}  | {samtools} view -hF4 | {samtools} sort -o {output_bam_file}"
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
        log = os.path.join(bam_dir, "medaka_consensus.log")
        cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model {ont_model} 2> {log}"
        tgt = f"{output_hdf_file}.OK"
        desc = f"Oxford Nanopore consensus contigs"
        mpm.run(cmd, tgt, desc)

        #  stitch
        input_hdf_file = os.path.join(bam_dir, "ont.hdf")
        output_fasta_file = os.path.join(fasta_dir, "ont.consensus.fasta")
        log = os.path.join(bam_dir, "medaka_stitch.log")
        cmd = f"{medaka} stitch --fill_char N {input_hdf_file} {reference_fasta_file} {output_fasta_file} 2> {log}"
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

    # write log file
    mpm.print_log()
    #print(" ".join(sys.argv))

    # copy files to trace
    copy2(__file__, trace_dir)
    subprocess.run(f'echo {" ".join(sys.argv)} > {trace_dir}/cmd.txt', shell=True, check=True)

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
