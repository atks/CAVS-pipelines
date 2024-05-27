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


@click.command()
@click.option(
    "-w",
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
    default="sample1",
    show_default=True,
    help="sample ID",
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
    sample_id,
    output_dir,
    reference_fasta_file,
):
    """
    Aligns all fastq files to a reference sequence file and generates a consensus sequence.

    e.g. align_and_consensus -1   -2   -n  -r ref.fasta
    """
    illumina = len(input_ilm_read1_fastq_files) != 0
    nanopore = len(input_ont_fastq_files) != 0

    # programs
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    medaka = "/home/atks/.miniconda3/bin/medaka"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

    # make directories
    fastq_dir = os.path.join(output_dir, "fastq")
    bam_dir = os.path.join(output_dir, "bam")
    fasta_dir = os.path.join(output_dir, "fasta")
    stats_dir = os.path.join(output_dir, "stats")
    consensus_dir = os.path.join(output_dir, "consensus")
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(consensus_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    # process illumina reads
    if illumina:
        print("Processing illumina reads")

        #  combine fastq files
        read1_files = input_ilm_read1_fastq_files.split(",")
        output_read1_fastq_file = os.path.join(fastq_dir, "ilm.r1.fastq.gz")
        cmd = f"zcat {' '.join(read1_files)} | gzip > {output_read1_fastq_file}"
        tgt = f"{output_read1_fastq_file}.OK"
        desc = f"Combining read 1 fastq files"
        run(cmd, tgt, desc)

        #  combine fastq files
        read2_files = " ".join(input_ilm_read2_fastq_files.split(","))
        output_read2_fastq_file = os.path.join(fastq_dir, "ilm.r2.fastq.gz")
        cmd = f"zcat {read2_files} | gzip > {output_read2_fastq_file}"
        tgt = f"{output_read2_fastq_file}.OK"
        desc = f"Combining read 2 fastq files"
        run(cmd, tgt, desc)

        #  construct reference
        cmd = f"{bwa} index -a bwtsw {reference_fasta_file}"
        tgt = f"{reference_fasta_file}.bwa_index.OK"
        desc = f"Construct bwa reference"
        run(cmd, tgt, desc)

        #  align
        output_bam_file = os.path.join(bam_dir, "ilm.bam")
        cmd = f"{bwa} mem -t 2 -M {reference_fasta_file} {output_read1_fastq_file} {output_read2_fastq_file} | {samtools} view -hF4 | {samtools} sort -o {output_bam_file}"
        tgt = f"{output_bam_file}.OK"
        desc = f"Align to reference with bwa"
        run(cmd, tgt, desc)

        #  index
        input_bam_file = os.path.join(bam_dir, "ilm.bam")
        cmd = f"{samtools} index {input_bam_file}"
        tgt = f"{input_bam_file}.bai.OK"
        desc = f"Index bam file"
        run(cmd, tgt, desc)

        #  coverage
        input_bam_file = os.path.join(bam_dir, "ilm.bam")
        output_stats_file = os.path.join(stats_dir, "ilm.stats.txt")
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        tgt = f"{output_stats_file}.OK"
        desc = f"Illumina coverage statistics"
        run(cmd, tgt, desc)

        #  consensus
        input_bam_file = os.path.join(bam_dir, "ilm.bam")
        output_fasta_file = os.path.join(fasta_dir, "ilm.consensus.fasta")
        cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
        tgt = f"{output_fasta_file}.OK"
        desc = f"Illumina consensus"
        run(cmd, tgt, desc)

    # process nanopore reads
    if nanopore:
        print("Processing nanopore reads")

        #  combine fastq files
        fastq_files = " ".join(input_ont_fastq_files.split(","))
        output_fastq_file = os.path.join(fastq_dir, "ont.fastq.gz")
        cmd = f"zcat {fastq_files} | gzip > {output_fastq_file}"
        tgt = f"{output_fastq_file}.OK"
        desc = f"Combining nanopore fastq files"
        run(cmd, tgt, desc)

        #  construct reference
        cmd = f"{minimap2} -d {reference_fasta_file}.mmi {reference_fasta_file}"
        tgt = f"{reference_fasta_file}.minimap2_index.OK"
        desc = f"Construct minimap2 reference"
        run(cmd, tgt, desc)

        #  align
        minimap2_ref_mmi_file = f"{reference_fasta_file}.mmi"
        input_fastq_file = os.path.join(fastq_dir, "ont.fastq.gz")
        output_bam_file = os.path.join(bam_dir, "ont.bam")
        cmd = f"{minimap2} -ax map-ont {minimap2_ref_mmi_file} {input_fastq_file}  | {samtools} view -hF4 | {samtools} sort -o {output_bam_file}"
        tgt = f"{output_bam_file}.OK"
        desc = f"Align to reference with minimap2"
        run(cmd, tgt, desc)

        #  index
        input_bam_file = os.path.join(bam_dir, "ont.bam")
        cmd = f"{samtools} index {input_bam_file}"
        tgt = f"{input_bam_file}.bai.OK"
        desc = f"Index bam file"
        run(cmd, tgt, desc)

        #  coverage
        input_bam_file = os.path.join(bam_dir, "ont.bam")
        output_stats_file = os.path.join(stats_dir, "ont.stats.txt")
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        tgt = f"{output_stats_file}.OK"
        desc = f"Nanopore coverage statistics"
        run(cmd, tgt, desc)

        #  consensus
        input_bam_file = os.path.join(bam_dir, "ont.bam")
        output_hdf_file = os.path.join(bam_dir, "ont.hdf")
        cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model r941_min_hac_g507"
        tgt = f"{output_hdf_file}.OK"
        desc = f"Oxford Nanopore consensus contigs"
        run(cmd, tgt, desc)

        #  consensus
        input_hdf_file = os.path.join(bam_dir, "ont.hdf")
        output_fasta_file = os.path.join(fasta_dir, "ont.consensus.fasta")
        cmd = f"{medaka} stitch --fill_char N {input_hdf_file} {reference_fasta_file} {output_fasta_file} "
        tgt = f"{output_fasta_file}.OK"
        desc = f"Oxford Nanopore consensus assembly/scaffold"
        run(cmd, tgt, desc)

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
        run(cmd, tgt, desc)

        # index bam
        input_bam_file = os.path.join(bam_dir, "ilm_ont.bam")
        cmd = f"{samtools} index {input_bam_file}"
        tgt = f"{output_bam_file}.bai.OK"
        desc = f"Index merged alignments"
        run(cmd, tgt, desc)

        # coverage
        input_bam_file = os.path.join(bam_dir, "ilm_ont.bam")
        output_stats_file = os.path.join(stats_dir, "ilm_ont.stats.txt")
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        tgt = f"{output_stats_file}.OK"
        desc = f"Compute coverage of merged alignments"
        run(cmd, tgt, desc)

        # consensus
        input_bam_file = os.path.join(bam_dir, "ilm_ont.bam")
        output_fasta_file = os.path.join(fasta_dir, "ilm_ont.consensus.fasta")
        cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
        tgt = f"{output_fasta_file}.OK"
        desc = f"Illumina and Nanopore consensus"
        run(cmd, tgt, desc)

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
    run(cmd, tgt, desc)


def run(cmd, tgt, desc):
    try:
        if os.path.exists(tgt):
            print(f"{desc} -  already executed")
            return
        else:
            print(f"{cmd}")
            subprocess.run(cmd, shell=True, check=True)
            subprocess.run(f"touch {tgt}", shell=True, check=True)
            print(f"{desc} -  successfully executed")
    except subprocess.CalledProcessError as e:
        print(f" - failed")
        exit(1)


if __name__ == "__main__":
    main() # type: ignore
