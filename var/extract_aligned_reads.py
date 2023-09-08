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
    "-d",
    "--output_fastq_dir",
    required=True,
    show_default=True,
    help="output fastq directory",
)
@click.option(
    "-o",
    "--output_fastq_root_name",
    required=True,
    show_default=True,
    help="output fastq root name",
)
@click.option(
    "-b",
    "--input_bam_file",
    required=True,
    default="",
    help="Aligned reads bam file",
)
@click.option(
    "-1",
    "--input_ilm_read1_fastq_file",
    required=False,
    default="",
    help="Illumina read 1 fastq file",
)
@click.option(
    "-2",
    "--input_ilm_read2_fastq_file",
    required=False,
    default="",
    help="Illumina read 2 fastq file",
)
def main(
    input_bam_file,
    input_ilm_read1_fastq_file,
    input_ilm_read2_fastq_file,
    output_fastq_dir,
):
    """
    Extracted aligned reads from a bam file and extracts from fastq file.

    e.g. extract_aligned_reads -b aligned.bam -1 r1.fastq.gz -2 r2.fastq.gz -o extracted_aligned_reads
    """

    # programs
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    seqtk = "/usr/local/seqtk-1.4/seqtk"

    # extract aligned read IDs from bam
    output_id_txt_file = f"{output_fastq_root}.id.txt"
    cmd = f"{samtools} view -F 4  {input_bam_file} | cut -f1 | sort | uniq  > {output_id_txt_file}"
    tgt = f"{output_id_txt_file}.OK"
    desc = f"Extract aligned read IDs"
    run(cmd, tgt, desc)

    # extract files from fastq files
    output_fastq_file1 = f"{output_fastq_dir}/{output_fastq_root_name}_r1.fastq.gz"
    cmd = f"{seqtk} subseq {input_ilm_read1_fastq_file} {output_id_txt_file} | gzip > {output_fastq_file1}"
    tgt = f"{output_fastq_file1}.OK"
    desc = f"Extract R1 fastq reads"
    run(cmd, tgt, desc)

    output_fastq_file1 = f"{output_fastq_dir}/{output_fastq_root_name}_r2.fastq.gz"
    cmd = f"{seqtk} subseq {input_ilm_read2_fastq_file} {output_id_txt_file} | gzip > {output_fastq_file2}"
    tgt = f"{output_fastq_file2}.OK"
    desc = f"Extract R2 fastq reads"
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
    main()
