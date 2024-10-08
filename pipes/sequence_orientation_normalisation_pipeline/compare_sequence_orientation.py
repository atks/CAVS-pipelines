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
    default=f"{os.getcwd()}/compare_sequence_orientation",
    show_default=True,
    help="output directory",
)
@click.option(
    "-q",
    "--query_fasta_file",
    required=True,
    show_default=True,
    help="query fasta file",
)
@click.option(
    "-r",
    "--reference_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
def main(
    output_dir, query_fasta_file, reference_fasta_file
):
    """
    Compare orientation of query sequence against reference sequence

    e.g. compare_sequence_orientation -q query.fasta -r ref.fasta
    """

    # create working directory
    # write out to separate files
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")
        exit(1)

    # version
    version = "1.0.0"

    # programs
    stretcher = "/usr/local/emboss-6.6.0/bin/stretcher"
    seqtk = "/usr/local/seqtk-1.4/seqtk"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/compare_sequence_orientation.log")
    mpm.set_ignore_targets(True)

    # read reference sequence
    seq = ""
    query_len = 0
    with open(query_fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                pass
            else:
                query_len += len(line.rstrip())

    # write out to separate files
    output_fasta_file = f"{output_dir}/query_fwd.fasta"
    cmd = f"{seqtk} seq -A {query_fasta_file} > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Create forward sequence FASTA file"
    mpm.run(cmd, tgt, desc)

    output_fasta_file = f"{output_dir}/query_rev.fasta"
    cmd = f"{seqtk} seq -Ar {query_fasta_file} > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"Create reverse completed sequence FASTA file"
    mpm.run(cmd, tgt, desc)

    # invoke needleman-wunsch alignment for each direction
    input_fasta_file = f"{output_dir}/query_fwd.fasta"
    alignment_file = f"{output_dir}/query_fwd_ref.stretcher"
    log_file = f"{output_dir}/query_fwd_ref.stretcher.log"
    cmd = f"{stretcher} {input_fasta_file} {reference_fasta_file} {alignment_file} 2> {log_file}"
    tgt = f"{alignment_file}.OK"
    desc = f"Align forward sequence to reference"
    mpm.run(cmd, tgt, desc)

    input_fasta_file = f"{output_dir}/query_rev.fasta"
    alignment_file = f"{output_dir}/query_rev_ref.stretcher"
    log_file = f"{output_dir}/query_rev_ref.stretcher.log"
    cmd = f"{stretcher} {input_fasta_file} {reference_fasta_file} {alignment_file} 2> {log_file}"
    tgt = f"{alignment_file}.OK"
    desc = f"Align reverse sequence to reference"
    mpm.run(cmd, tgt, desc)

    # examine alignments
    query_fwd_alignment = parse_stretcher_alignment(f"{output_dir}/query_fwd_ref.stretcher")
    query_rev_alignment = parse_stretcher_alignment(f"{output_dir}/query_rev_ref.stretcher")

    # find best alignment.
    in_same_orientation = (
        True
        if query_fwd_alignment.score > query_rev_alignment.score
        else False
    )

    result_file  = f"{output_dir}/result.txt"
    with open(result_file, "w") as file:
        file.write(f"#file\tqseq\trseq\tqry_len\tfwd_align_len\trev_align_len\tfwd_score\trev_score\tfwd_similarity\trev_similarity\tsame_orientation\n")
        fwd_identity_percent = float(query_fwd_alignment.identity)/float(query_fwd_alignment.length)*100
        rev_identity_percent = float(query_rev_alignment.identity)/float(query_rev_alignment.length)*100
        file.write(f"{query_fasta_file}\t{query_fwd_alignment.qseq}\t{query_fwd_alignment.rseq}\t{query_len}\t{query_fwd_alignment.length}\t{query_rev_alignment.length}\t{query_fwd_alignment.score}\t{query_rev_alignment.score}\t{fwd_identity_percent:.2f}\t{rev_identity_percent:.2f}\t{in_same_orientation}\n")


    print(f"Compare sequence orientation stats")
    print(f'In same orientation: {in_same_orientation}')
    print("==========")
    query_fwd_alignment.print()
    print("==========")
    query_rev_alignment.print()
    print("==========")

    # copy files to trace
    copy2(__file__, trace_dir)

    # write log file
    mpm.print_log()

#=======================================
#
# Aligned_sequences: 2
# 1: FR682468.2
# 2: AF270707.1
# Matrix: EDNAFULL
# Gap_penalty: 16
# Extend_penalty: 4
#
# Length: 417
# Identity:     402/417 (96.4%)
# Similarity:   402/417 (96.4%)
# Gaps:           0/417 ( 0.0%)
# Score: 1950
#
#
#=======================================
def parse_stretcher_alignment(file):
    with open(file, "r") as file:
        qseq = ""
        rseq = ""
        length = 0
        identity = 0
        gaps = 0
        score = 0
        beg = 10000000000
        end = -1
        for line in file:
            if line.startswith("# 1:"):
                m = re.search(r"\# 1: (.+)", line)
                if m is not None:
                    qseq = m.group(1)
                else:
                    exit("Cannot parse qseq")
            elif line.startswith("# 2:"):
                m = re.search(r"\# 2: (.+)", line)
                if m is not None:
                    rseq = m.group(1)
                else:
                    exit("Cannot parse rseq")
            elif line.startswith("# Length"):
                m = re.search(r"(\d+)", line)
                if m is not None:
                    length = int(m.group(1))
                else:
                    exit("Cannot parse length")
            elif line.startswith("# Identity"):
                m = re.search(r"(\d+)\/", line)
                if m is not None:
                    identity = int(m.group(1))
                else:
                    exit("Cannot parse identity")
            elif line.startswith("# Gaps"):
                m = re.search(r"(\d+)\/", line)
                if m is not None:
                    gaps = int(m.group(1))
                else:
                    exit("Cannot parse gaps")
            elif line.startswith("# Score"):
                m = re.search(r"([\d\.]+)", line)
                if m is not None:
                    score = float(m.group(1))
                else:
                    exit("Cannot parse score")
            else:
                pass

        alignment = Alignment(qseq, rseq, length, identity, gaps, score)
        return alignment

class Alignment(object):
    def __init__(self, qseq, rseq, length, identity, gaps, score):
        self.qseq = qseq
        self.rseq = rseq
        self.length = length
        self.identity = identity
        self.gaps = gaps
        self.score = score

    def print(self):
        print(f"qseq      : {self.qseq}")
        print(f"rseq      : {self.rseq}")
        print(f"length    : {self.length}")
        print(f"identity  : {self.identity}")
        print(f"gaps      : {self.gaps}")
        print(f"score     : {self.score}")

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


if __name__ == "__main__":
    main() # type: ignore