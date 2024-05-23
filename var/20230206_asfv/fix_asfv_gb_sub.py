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
@click.argument("input_prokka_tbl_file")
@click.option(
    "-r",
    "--ref_proteins_faa_file",
    required=True,
    show_default=True,
    help="reference proteins fasta file",
)
@click.option(
    "-o",
    "--output_annotated_tbl_file",
    required=True,
    show_default=True,
    help="output annotated tbl file",
)
@click.option(
    "-d",
    "--debug",
    is_flag=True,
    default=False,
    help="output debug statements",
)
def main(
    input_prokka_tbl_file, ref_proteins_faa_file, output_annotated_tbl_file, debug
):
    """
    Fix prokka ASFV annotations

    e.g. fix_asfv_gb_sub.py prokka.tbl -r ref.faa -o annotated.tbl
    """

    # read reference protein file
    ref_proteins = {}
    id = ""
    gene = ""
    sequence = ""
    with open(ref_proteins_faa_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                m = re.match(r">(\S+) ~~~(.+)~~~", line)
                id = m.group(1)
                gene = m.group(2)
            else:
                sequence = line.strip()
                ref_proteins[id] = Protein(id, gene, sequence)

    if debug:
        print(f"No of protein records read: {len(ref_proteins)}")

    #    CDS             complement(9490..10107)
    #                      /gene="MGF 110-5L-6L CDS"
    #                      /codon_start=1
    #                      /protein_id="CAD2068349.1"
    #                      /translation="MLVIFLGILGLLASQVSSQLVGQLRPTEEPPEEELEYWCAYMES
    #                      CQFCWDCQDGTCINKIDGSVIYKNEYVKSCLVSRWLDKCMYDLDKGIYHTMNCNQVLG
    #                      LPNQPAGQLHPTDNPPQEELEYWCTYTENCKFCWNCQNGLCEGKLENTTILENEYVQS
    #                      CIVSRWLNKCMYDLGQGIHHVMACSEPKPWNPYKILKREWKENNS"

    # 186926	187996	CDS
    # 			inference	ab initio prediction:Prodigal:002006
    # 			inference	similar to AA sequence:FR682468.2.gb:CAD2068525.1
    # 			locus_tag	KEICJGFC_00173
    # 			note
    # 			product	hypothetical protein
    # 188830	188633	CDS
    # 			inference	ab initio prediction:Prodigal:002006
    # 			locus_tag	KEICJGFC_00174
    # 			product	hypothetical protein

    first_record = True
    beg = 0
    end = 0
    name = ""
    gene = ""
    protein_id = ""
    expected_plen = 0
    n = 1
    no_len_discrepancy = 0
    with open(output_annotated_tbl_file, "w") as of:
        with open(input_prokka_tbl_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    if debug:
                        print(line.strip())
                else:
                    if not line.startswith("\t"):
                        # print out record
                        if not first_record:
                            if debug:
                                print("#===============================")
                                print("#output  last record")
                            feature = Feature(
                                n, beg, end, name, gene, protein_id, expected_plen
                            )
                            if feature.plen != feature.expected_plen:
                                no_len_discrepancy += 1
                                # feature.print()
                            of.write(f"{feature.beg}\t{feature.end}\t{feature.name}\n")
                            of.write(f"\t\t\tid\tprot_{feature.n}\n")
                            of.write(f"\t\t\tgene\t{feature.gene}\n")
                            of.write(f"\t\t\tprotein_id\t{feature.protein_id}\n")
                            if debug:
                                feature.print()

                            beg = 0
                            end = 0
                            name = ""
                            gene = ""
                            protein_id = ""
                            expected_plen = 0
                            n += 1
                            if debug:
                                print("#===============================")

                        first_record = False

                        if debug:
                            print("processing record")
                        beg, end, name = line.strip().split("\t")
                        beg = int(beg)
                        end = int(end)
                        if debug:
                            print(f"\tread in {beg} {end} {features}")
                        # if protein_id in ref_proteins:
                        #     ref_protein = ref_proteins[protein_id]
                        #     if ref_protein.sequence == name:
                        #         print(
                        #             f"Matched {protein_id} {ref_protein.gene} {ref_protein.sequence}"
                        #         )
                        #     else:
                        #         print(
                        #             f"Unmatched {protein_id} {ref_protein.gene} {ref_protein.sequence}"
                        #         )
                        #         print(f"  {name}")
                        # else:
                        #     print(f"Missing {protein_id}")

                    else:
                        if debug:
                            print(f"\t{line.strip()}")
                        if re.search("AA sequence", line.strip()):
                            m = re.match(r".+AA sequence:[^:]+:(.+)", line.strip())
                            protein_id = m.group(1)
                            # print(f"\t\tprotein_id {protein_id}")
                            if protein_id in ref_proteins:
                                gene = ref_proteins[protein_id].gene
                                expected_plen = len(ref_proteins[protein_id].sequence)

                        # ref_proteins[id] = Protein(id, gene, sequence)
            if debug:
                print("#===============================")
            if debug:
                print("#output last record")
            feature = Feature(n, beg, end, name, gene, protein_id, expected_plen)
            of.write(f"{feature.beg}\t{feature.end}\t{feature.name}\n")
            of.write(f"\t\t\tid\tprot_{feature.n}\n")
            of.write(f"\t\t\tgene\t{feature.gene}\n")
            of.write(f"\t\t\tprotein_id\t{feature.protein_id}\n")
            if feature.plen != feature.expected_plen:
                no_len_discrepancy += 1
                # print("===========")
                # feature.print()
            if debug:
                feature.print()
            if debug:
                print("#===============================")
            print(f"no of length discrepancies:  {no_len_discrepancy}")


class Protein(object):
    def __init__(self):
        self.id = ""
        self.gene = ""
        self.sequence = ""

    def __init__(self, id, gene, sequence):
        self.id = id
        self.gene = gene
        self.sequence = sequence

    def print(self):
        print(f"id        : {self.id}")
        print(f"gene      : {self.gene}")
        print(f"sequence  : {self.sequence}")


class Feature(object):
    def __init__(self):
        self.n = 0
        self.beg = 0
        self.end = 0
        self.name = ""
        self.gene = ""
        self.protein_id = ""
        self.plen = 0
        self.expected_plen = 0

    def __init__(self, n, beg, end, name, gene, protein_id, expected_plen):
        self.n = n
        self.beg = beg
        self.end = end
        self.name = name
        self.gene = gene
        self.protein_id = protein_id
        self.plen = (abs(beg - end) + 1) / 3 - 1
        self.expected_plen = expected_plen

    def print(self):
        print(f"n            : {self.n}")
        print(f"beg          : {self.beg}")
        print(f"end          : {self.end}")
        print(f"name         : {self.name}")
        print(f"gene         : {self.gene}")
        print(f"protein_id   : {self.protein_id}")
        print(f"plen          : {self.plen}")
        print(f"expected plen : {self.expected_plen}")


if __name__ == "__main__":
    main()
