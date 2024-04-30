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
import sys

@click.command()
@click.argument("analysis_dir", required=True)
def main(analysis_dir):
    """
    Takes in an analysis directory from a CAVS data freeze and inspects the fastqc_result directory to ensure run completion.
    This is required as fastqc does not always proper output a failure status.

    e.g. check_fastqc_runs.py.py analysis
    """
    fastq_path = ""

    # 01_wild_e_coli_756-21_R1.fastq.gz
    # 01_wild_e_coli_756-21_R1_fastqc.html
    # 01_wild_e_coli_756-21_R1_fastqc.zip
    # 01_wild_e_coli_756-21_R2.fastq.gz
    # 01_wild_e_coli_756-21_R2_fastqc.html
    # 01_wild_e_coli_756-21_R2_fastqc.zip

    dirs = os.listdir(analysis_dir)
    new_dirs = []
    for dirname in dirs:
        if dirname != "all":
            index, sample_id = dirname.split("_", 1)
            padded_index = index.zfill(2)
            new_dirs.append((int(index), padded_index, dirname, sample_id))

    new_dirs.sort(key=lambda x: x[0])

    all_good = True
    for index, padded_index, dirname, sample_id in new_dirs:
        # print(f"Checking {dirname}")
        file = f"{analysis_dir}/{dirname}/fastqc_result/{padded_index}_{sample_id}_R1_fastqc.html"
        if not os.path.isfile(file):
            print(f"{file} missing", file=sys.stderr)
            all_good = False
        file = f"{analysis_dir}/{dirname}/fastqc_result/{padded_index}_{sample_id}_R2_fastqc.html"
        if not os.path.isfile(file):
            print(f"{file} missing", file=sys.stderr)
            all_good = False
        file = f"{analysis_dir}/{dirname}/fastqc_result/{padded_index}_{sample_id}_R1_fastqc.zip"
        if not os.path.isfile(file):
            print(f"{file} missing", file=sys.stderr)
            all_good = False
        file = f"{analysis_dir}/{dirname}/fastqc_result/{padded_index}_{sample_id}_R2_fastqc.zip"
        if not os.path.isfile(file):
            print(f"{file} missing", file=sys.stderr)
            all_good = False

    if all_good:
        exit(0)
    else:
        exit(1)

if __name__ == "__main__":
    main()  # type: ignore
