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
import re
from shutil import copy2

@click.command()
@click.option(
    "-o",
    "--output_fasta_file",
    required=True,
    default="",
    show_default=True,
    help="output directory",
)
@click.option(
    "-s",
    "--sequence_orientation_report_file",
    required=True,
    show_default=True,
    help="sequence orientation report file",
)
def main(output_fasta_file, sequence_orientation_report_file):
    """
    Normalise orientation of sequences from an orientation report file.

    e.g. normalise_sequence_orientation -s orientation_report.txt -o normalised.fasta
    """

    # version
    version = "1.0.0"

    # programs
    seqtk = "/usr/local/seqtk-1.4/seqtk"

    # initialize
    mpm = MiniPipeManager("")
    mpm.set_ignore_targets(True)

    # sequence orientation report file
    index = 0
    with open(sequence_orientation_report_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                file_name, qseq, rseq, qry_len,fwd_align_len, rev_align_len, fwd_score, rev_score, fwd_similarity, rev_similarity, same_orientation = line.rstrip().split("\t")
                if index == 0:
                    if same_orientation == "True":
                        cmd = f"{seqtk} seq -A {file_name} > {output_fasta_file}"
                        mpm.run_lite(cmd)
                    else:
                        cmd = f"{seqtk} seq -Ar {file_name} > {output_fasta_file}"
                        mpm.run_lite(cmd)
                else:
                    if same_orientation == "True":
                        cmd = f"{seqtk} seq -A {file_name} >> {output_fasta_file}"
                        mpm.run_lite(cmd)
                    else:
                        cmd = f"{seqtk} seq -Ar {file_name} >> {output_fasta_file}"
                        mpm.run_lite(cmd)
                index += 1

class MiniPipeManager(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log_msg = []
        self.ignore_targets = False

    def run(self, cmd, tgt, desc):
        try:
            if not self.ignore_targets and os.path.exists(tgt):
                self.log(f"{desc} -  already executed")
                self.log(cmd)
                return
            else:
                self.log(f"{desc}")
                subprocess.run(cmd, shell=True, check=True)
                subprocess.run(f"touch {tgt}", shell=True, check=True)
                self.log(cmd)
        except subprocess.CalledProcessError as e:
            self.log(f" - failed")
            exit(1)

    def set_ignore_targets(self, ignore_targets):
        self.ignore_targets = ignore_targets

    def log(self, msg):
        print(msg)
        self.log_msg.append(msg)

    def print_log(self):
        self.log(f"\nlogs written to {self.log_file}")
        with open(self.log_file, "w") as f:
            f.write("\n".join(self.log_msg))

    def run_lite(self, cmd):
        try:
            print(cmd)
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(" - failed")
            exit(1)

if __name__ == "__main__":
    main() # type: ignore