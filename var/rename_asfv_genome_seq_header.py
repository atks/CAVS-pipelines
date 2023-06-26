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

import os.path
import click
import re


@click.command()
@click.argument("input_fasta_file", nargs=1)
@click.option("-n", "--rename_file", required=True, help="reference fasta file")
def main(input_fasta_file, rename_file):
    """
    Rename FASTA headers

    e.g. rename_asfv_genome_seq_header.py target.fa -n name.txt
    """

    headers = dict()
    with open(rename_file) as file:
        for line in file:
            acc, country, year, header = line.strip().split("\t")
            headers[header] = f">{acc}_{country}_{year}".replace(" ", "_")

    with open(input_fasta_file) as file:
        for line in file:
            if line.startswith(">"):
                print(headers[line.strip()])
            else:
                print(line, end="")


if __name__ == "__main__":
    main()
