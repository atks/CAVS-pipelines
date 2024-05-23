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
import sys
import click
import re


@click.command()
@click.argument("genbank_file", nargs=1)
@click.option("-o", "--bed_file", required=False, default="", help="output bed file")
def main(genbank_file, bed_file):
    """
    Extracts genes from a genbank file

    e.g. genbank2genebed.py ref.genbank -o out.bed
    """

    with open(genbank_file) as ifile:
        if bed_file == "":
            ofile = sys.stdout
        else:
            ofile = open(bed_file, "w")

        gene_name_line = False
        beg = -1
        end = -1
        id = ""
        for line in ifile:
            line = line.strip()
            if line.startswith("VERSION"):
                m = re.match("VERSION\s+(\S+)", line)
                id = m.group(1)
            elif line.startswith("ORIGIN"):
                break
            elif line.startswith("gene"):
                # ofile.write(f'\t{line}\n')
                m = re.match("[^\d]+(\d+)\.\.>?(\d+).*", line)
                beg = int(m.group(1))
                end = int(m.group(2))
                gene_name_line = True
            elif gene_name_line:
                # ofile.write(f'\t\t{line}\n')
                m = re.match('/gene="(.+)"', line)
                gene = m.group(1)
                if beg <= end:
                    ofile.write(f"{id}\t{beg-1}\t{end}\t{gene}\n")
                else:
                    ofile.write(f"\t\tOOPS::{gene}\t{beg-1}\t{end}\n")
                gene_name_line = False
                beg = -1
                end = -1
            else:
                pass


if __name__ == "__main__":
    main()
