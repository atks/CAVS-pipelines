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
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
def main(
    reference_fasta_file,
):
    """
    Aligns all fastq files to a reference sequence file and generates a consensus sequence.

    e.g. index reference file for bwa and minimap2
    """
    # programs
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    bwa = "/usr/local/bwa-0.7.17/bwa"

    #  construct bwa reference
    cmd = f"{bwa} index -a bwtsw {reference_fasta_file}"
    tgt = f"{reference_fasta_file}.bwa_index.OK"
    desc = f"Construct bwa reference"
    run(cmd, tgt, desc)

    #  construct minimap2 reference
    cmd = f"{minimap2} -d {reference_fasta_file}.mmi {reference_fasta_file}"
    tgt = f"{reference_fasta_file}.minimap2_index.OK"
    desc = f"Construct minimap2 reference"
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
