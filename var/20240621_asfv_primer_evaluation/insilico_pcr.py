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
    "--output_dir",
    default=f"{os.getcwd()}/insilico_pcr",
    show_default=True,
    help="output directory",
)
@click.option(
    "-p",
    "--primer_fasta_file",
    required=True,
    show_default=True,
    help="primers fasta file",
)
@click.option(
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
def main(output_dir, primer_fasta_file, reference_fasta_file):
    """
    Performs in silico PCR on a sample fasta file from a set of degenerate primers

    e.g. insilico_pcr.py -p primers.fasta -r sample.fasta
    """

    # version
    version = "1.0.0"

    # programs
    bowtie2_build = "/usr/local/bowtie2-2.4.5/bowtie2-build"
    bowtie2 = "/usr/local/bowtie2-2.4.5/bowtie2"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/extract_amplicon.log")

    # build bowtie2 data base from sample fasta file

    # read primer sequence pairs
    # primers must be named in pairs and be labelled forward and reverse
    # >primer1_fwd
    # >primer1_rev
    #
    # suppose you have
    #
    # a pair of flavivirus forward primers
    # TACAACATgATggggAARAgAgARAA / TACAACATgATgggMAAACgYgARAA
    #
    # and a flavivirus reverse primer
    # gTgTCCCAGCCNgCKgTRTCRTC
    #
    # the you should generate a fasta file with the following content
    #
    # >flavi1_fwd
    # TACAACATgATggggAARAgAgARAA
    # >flavi1_rev
    # gTgTCCCAGCCNgCKgTRTCRTC
    # >flavi2_fwd
    # TACAACATgATgggMAAACgYgARAA
    # >flavi2_rev
    # gTgTCCCAGCCNgCKgTRTCRTC

    primers = dict()
    with open(primer_fasta_file, "r") as file:
        for line in file:
            fwd_primer_name, fwd_orientation = line.strip().lstrip(">").split("_")
            fwd_seq = next(file).strip()
            rev_primer_name, rev_orientation = next(file).strip().lstrip(">").split("_")
            rev_seq = next(file).strip()
            if fwd_orientation != "fwd" or rev_orientation != "rev":
                print("primer names must be labeled fwd and rev")
                exit(1)
            if fwd_primer_name != rev_primer_name:
                print("Primer names must be the same")
                exit(1)
            primers[fwd_primer_name] = PrimerSet(fwd_primer_name, fwd_seq, rev_seq)
            primers[fwd_primer_name].print()

    # write out to separate files
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # copy files to trace
    copy2(__file__, trace_dir)

    # write log file
    mpm.print_log()

class PrimerSet(object):
    def __init__(self, name, fwd_primer, rev_primer):
        self.name = name
        self.fwd_primer = fwd_primer
        self.rev_primer = rev_primer

    def print(self):
        print(f"name       : {self.name}")
        print(f"fwd_primer : {self.fwd_primer}")
        print(f"rev_primer : {self.rev_primer}")

def reverse_complement(seq):
    return (
        seq[::-1]
        .upper()
        .replace("A", "t")
        .replace("T", "a")
        .replace("G", "c")
        .replace("C", "g")
        .replace("M", "k")
        .replace("R", "y")
        .replace("W", "w")
        .replace("S", "s")
        .replace("Y", "r")
        .replace("K", "m")
        .replace("V", "b")
        .replace("H", "d")
        .replace("D", "h")
        .replace("B", "v")
        .replace("N", "n")
        .upper()
    )


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
                subprocess.run(cmd, shell=True, check=True)
                subprocess.run(f"touch {tgt}", shell=True, check=True)
                self.log(cmd)
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
    main() # type: ignore[arg-type]
