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
import re


@click.command()
@click.argument("input_fasta_file")
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-r",
    "--ref_fasta_file",
    required=True,
    show_default=True,
    help="reference fasta file",
)
def main(working_dir, input_fasta_file, ref_fasta_file):
    """
    Extracts similar subsequences from a FASTA file

    e.g. extract_asfv_p72_segment.py 44genomes.complete.asfv.fasta -r MT851941.1.fasta
    """

    # write out to separate files
    output_dir = f"{working_dir}/extract_segment_output"
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    # read reference sequences
    sequences = parse_fasta(input_fasta_file)
    ref_seq = parse_fasta(ref_fasta_file)[0]
    ref_seq.reverse_complement()
    write_fasta(ref_seq, f"{output_dir}/rc_seq.fasta")
    rc_ref_fasta_file = f"{output_dir}/rc_seq.fasta"

    for seq in sequences:
        target_fasta_file = f"{output_dir}/{seq.chrom}.fasta"
        output_alignment_file = f"{output_dir}/{seq.chrom}.water"
        output_rc_alignment_file = f"{output_dir}/rc_{seq.chrom}.water"
        write_fasta(seq, f"{output_dir}/{seq.chrom}.fasta")

        # print(f'processing {seq.chrom}')
        # print(f'\t{ref_fasta_file}')
        # print(f'\t{target_fasta_file}')
        # print(f'\t{output_alignment_file}')

        # print("==============")
        alignment = get_water_alignment(
            ref_fasta_file, target_fasta_file, output_alignment_file
        )
        # alignment.print()
        # print("++++++++++++++")
        rc_alignment = get_water_alignment(
            rc_ref_fasta_file, target_fasta_file, output_rc_alignment_file
        )
        # rc_alignment.print()
        # print(
        #     f'{alignment.rseq}\t{alignment.get_similarity_score()}\t{rc_alignment.get_similarity_score()}')

        # print("==============")
        print(f"{seq.desc}\n{seq.seq[rc_alignment.beg-1:rc_alignment.end]}")


def get_water_alignment(ref_fasta_file, target_fasta_file, output_alignment_file):
    water = "/usr/local/emboss-6.6.0/bin/water"
    result = subprocess.run(
        [
            water,
            ref_fasta_file,
            target_fasta_file,
            output_alignment_file,
            "-gapopen",
            "10",
            "-gapextend",
            "0.5",
        ],
        capture_output=True,
    )
    if result:
        return parse_water_alignment(output_alignment_file)


# Aligned_sequences: 2
# 1: primer2_reverse
# 2: AY261360.1
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 19
# Identity:      15/19 (78.9%)
# Similarity:    15/19 (78.9%)
# Gaps:           2/19 (10.5%)
# Score: 56.5


def parse_water_alignment(file):
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
                m = re.search("\# 1: (.+)", line)
                qseq = m.group(1)
            elif line.startswith("# 2:"):
                m = re.search("\# 2: (.+)", line)
                rseq = m.group(1)
            elif line.startswith("# Length"):
                m = re.search("(\d+)", line)
                length = int(m.group(1))
            elif line.startswith("# Identity"):
                m = re.search("(\d+)\/", line)
                identity = m.group(1)
            elif line.startswith("# Gaps"):
                m = re.search("(\d+)\/", line)
                gaps = m.group(1)
            elif line.startswith("# Score"):
                m = re.search("([\d\.]+)", line)
                score = m.group(1)
            elif len(rseq) != 0 and line.startswith(rseq):
                m = re.search(f"{rseq}\s+([\d\.]+) [^\d]+ ([\d\.]+)", line)
                beg = int(m.group(1)) if int(m.group(1)) < beg else beg
                end = int(m.group(2)) if int(m.group(2)) > end else end
            else:
                pass

        alignment = Alignment(qseq, rseq, length, identity, gaps, score, beg, end)
        return alignment


class Alignment(object):
    def __init__(self):
        self.qseq = ""
        self.rseq = ""
        self.length = 0
        self.identity = 0
        self.gaps = 0
        self.score = 0
        self.beg = 0
        self.end = 0

    def __init__(self, qseq, rseq, length, identity, gaps, score, beg, end):
        self.qseq = qseq
        self.rseq = rseq
        self.length = length
        self.identity = identity
        self.score = score
        self.gaps = gaps
        self.beg = beg
        self.end = end

    def get_similarity_score(self):
        return int(self.identity) / int(self.length) * 100

    def get_similarity(self):
        return f"{int(self.identity)/int(self.length)*100:.2f}% ({self.identity}/{self.length})"

    def print(self):
        print(f"qseq      : {self.qseq}")
        print(f"rseq      : {self.rseq}")
        print(f"length    : {self.length}")
        print(f"identity  : {self.identity}")
        print(
            f"identity  : {int(self.identity)/int(self.length)*100:.2f}% ({self.identity}/{self.length})"
        )
        print(f"gaps      : {self.gaps}")
        print(f"score     : {self.score}")
        print(f"beg       : {self.beg}")
        print(f"end       : {self.end}")


def parse_fasta(file):
    sequences = []
    desc = ""
    seq = ""
    with open(file, "r") as file:
        for line in file:
            if line.startswith(">"):
                if desc != "":
                    sequences.append(FastaEntry(desc, seq))
                desc = line.rstrip()
                seq = ""
            else:
                seq += line.rstrip()
        if desc != "":
            sequences.append(FastaEntry(desc, seq))
    return sequences


def write_fasta(sequence, file):
    with open(file, "w") as file:
        file.write(sequence.desc)
        file.write("\n")
        file.write(sequence.seq)


class FastaEntry(object):
    def __init__(self):
        self.desc = ""
        self.chrom = ""
        self.seq = ""

    def __init__(self, desc, seq):
        self.desc = desc
        self.chrom = desc.lstrip(">").split(" ")[0]
        self.seq = seq

    def reverse_complement(self):
        rc_seq = (
            self.seq[::-1]
            .replace("A", "t")
            .replace("T", "a")
            .replace("G", "c")
            .replace("C", "g")
            .upper()
        )
        self.seq = rc_seq

    def print(self):
        print(f"desc  : {self.desc}")
        print(f"chrom : {self.chrom}")
        print(f"seq   : {self.seq}")


class Record(object):
    def __init__(self):
        self.id = ""
        self.origin = ""
        self.serotype = ""

    def __init__(self, id, origin, serotype):
        self.id = id
        self.origin = origin
        self.serotype = serotype

    def print(self):
        print(f"id        : {self.id}")
        print(f"origin    : {self.origin}")
        print(f"serotype  : {self.serotype}")


if __name__ == "__main__":
    main()
