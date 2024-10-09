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

    # version
    version = "0.0.1"

    # get reference sequence
    ref_id = ""
    ref_seq = ""
    with open(ref_fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                ref_id = line.rstrip().split(" ")[0][1:]
            else:
                ref_seq += line.rstrip()

    ref_len = len(ref_seq)
    print(f"header: {ref_id}")
    print(f"length: {ref_len}")

    with open(input_tbl_file, "r") as f:
        lines = f.readlines()
        for line in f:
            if line.startswith("LOCUS"):
                break

if __name__ == "__main__":
    main() # type: ignore
