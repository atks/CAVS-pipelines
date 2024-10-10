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
import re
from shutil import copy2

@click.command()
@click.argument("input_tbl_file")
@click.option(
    "-o",
    "--output_tbl_file",
    default="filtered.tbl",
    required=False,
    show_default=True,
    help="output genbank submission table file",
)
@click.option(
    "-r",
    "--ref_fasta_file",
    required=True,
    help="reference fasta file",
)
def main(
    input_tbl_file,
    output_tbl_file,
    ref_fasta_file
):
    """
    Filters prokka genbank submission file.

    e.g.filter_prokka_tbl.py in.tbl
    """
    print("\t{0:<20} :   {1:<10}".format("input tbl file", input_tbl_file))
    print("\t{0:<20} :   {1:<10}".format("output tbl file", output_tbl_file))
    print("\t{0:<20} :   {1:<10}".format("reference fasta file", ref_fasta_file))

    # version
    version = "0.0.1"

    # get reference sequence
    ref_id = ""
    ref_seq = ""

    with open(output_tbl_file, "w") as e:
        with open(ref_fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    ref_id = line.rstrip().split(" ")[0][1:]
                else:
                    ref_seq += line.rstrip()

        ref_len = len(ref_seq)
        print("")
        print(f"reference sequence")
        print(f"header: {ref_id}")
        print(f"length: {ref_len}")

        # parse input tbl file
        feature_table = FeatureTable()
        feature_table.parse(input_tbl_file)

        print("")
        print(f"feature table parsed")
        print(f"no. records: {len(feature_table.features)}")

        # filter features
        no_genes = 0
        no_stop_codons = 0
        no_fwd_strand = 0
        no_rev_strand = 0
        for feature in feature_table.features:
            #delete ab initio prediction annotation
            is_gene = False

            for annotation in feature.annotations:
                if annotation.tag=="inference" and annotation.description.startswith("ab initio prediction"):
                    feature.annotations.remove(annotation)
                if annotation.tag=="gene":
                    is_gene = True

            if is_gene:
                no_genes += 1

                #check stop codon
                codon = ""
                is_stop_codon = False
                if feature.end>feature.beg:
                    no_fwd_strand += 1
                    seq = ref_seq[feature.beg-1:feature.end]
                    codon = ref_seq[feature.end-3:feature.end]
                    if codon=="TAA" or codon=="TAG" or codon=="TGA":
                        no_stop_codons += 1
                        is_stop_codon = True
                else:
                    no_rev_strand += 1
                    seq = ref_seq[feature.end-1:feature.beg]
                    codon = ref_seq[feature.end-1:feature.end+2]
                    if codon=="TTA" or codon=="CTA" or codon=="TCA":
                        no_stop_codons += 1
                        is_stop_codon = True

                #check to ensure stop codon is present
                if not is_stop_codon:
                    seq = ref_seq[feature.end-1:feature.beg]
                    print(f"seq: {seq}")
                    print(f"codon: {codon}")

                #edit gene annotation locus_tag
                for annotation in feature.annotations:
                    if annotation.tag=="locus_tag":
                        annotation.description = f"EHIMFBMF_{no_genes}"

                feature.write(e)


    print(f"no. gene records: {no_genes}")
    print(f"\tno. forward strand: {no_fwd_strand}")
    print(f"\tno. reverse strand: {no_rev_strand}")
    print(f"\tno. stop codons: {no_stop_codons}")
class FeatureTable(object):
    def __init__(self):
        self.features = []

    def parse(self, file):
        feature_text = ""
        with open(file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                if re.search(r"^\d", line):
                    if len(feature_text)!=0:
                        feature = Feature()
                        feature.parse(feature_text)
                        self.features.append(feature)
                        feature_text = ""
                    feature_text = line
                else:
                    feature_text += line

            if len(feature_text)!=0:
                feature = Feature()
                feature.parse(feature_text)
                self.features.append(feature)
                feature_text = ""

class Feature(object):
    def __init__(self):
        self.beg = 0
        self.end = 0
        self.feature = ""
        self.annotations = []

    def parse(self, text):
        #print(text)
        lines = text.split("\n")
        beg, end, self.feature = lines[0].split("\t")
        self.beg = int(beg)
        self.end = int(end)
        for line in lines[1:]:
            if len(line)==0:
                continue
            tag, description = line.strip().split("\t")
            annotation = Annotation(tag, description)
            self.annotations.append(annotation)

    def print(self):
        print(f"{self.beg}\t{self.end}\t{self.feature}")
        for annotation in self.annotations:
            annotation.print()

    def write(self, file):
        file.write(f"{self.beg}\t{self.end}\t{self.feature}\n")
        for annotation in self.annotations:
            annotation.write(file)

class Annotation(object):
    def __init__(self, tag, description):
        self.tag = tag
        self.description = description

    def print(self):
        print(f"\t\t\t{self.tag}\t{self.description}")

    def write(self, file):
        file.write(f"\t\t\t{self.tag}\t{self.description}\n")

if __name__ == "__main__":
    main() # type: ignore
