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

import os.path
import click
import re


@click.command()
@click.argument("input_genbank_files", nargs=-1)
@click.option(
    "-w",
    "--write",
    required=False,
    default=False,
    show_default=True,
    help="reference fasta file",
)
def main(input_genbank_files, write):
    """
    Extract country IDs from genbank files.

    e.g. genbank_extract_country_year.py
    """
    n = 0
    n_country_annotated = 0
    n_date_annotated = 0

    for genbank_file in input_genbank_files:
        n += 1
        acc = ""
        header = ""
        with open(genbank_file) as file:
            for line in file:
                if line.startswith(">"):
                    header = line.strip()
                    acc = header.lstrip(">").split(" ")[0]
                    break
        genbank_file = input_fasta_file.replace("fasta", "seq")
        country = "no country"
        date = "no date"
        year = "no year"
        if os.path.isfile(genbank_file):
            with open(genbank_file) as file:
                for line in file:
                    if line.startswith("ORIGIN"):
                        break
                    else:
                        line = line.strip()
                        m = re.match(r'/country="(.+)"', line)
                        if m:
                            n_country_annotated += 1
                            country = m.group(1)
                            continue
                        m = re.match(r"JOURNAL   Submitted \((.+)\).+", line)
                        if m:
                            n_date_annotated += 1
                            date = m.group(1)
                            year = date.split("-")[2]
                            continue

        print(f"{acc}\t{country}\t{year}\t{header}")

    print(f"no fasta files       : {n}")
    print(f"no country annotated : {n_country_annotated}")
    print(f"no date annotated    : {n_date_annotated}")

    #  JOURNAL   Submitted (15-MAY-2007) Chapman D.A., Microbiology, Institute Of
    #             Animal Health, Pirbright Laboratory, Ash Road, Pirbright, Woking,
    #             Surrey, GU24 0NF, UNITED KINGDOM

    # FEATURES             Location/Qualifiers
    #      source          1..182284
    #                      /organism="African swine fever virus Benin 97/1"
    #                      /mol_type="genomic DNA"
    #                      /isolate="Benin 97/1"
    #                      /db_xref="taxon:443876"
    #                      /country="Benin"


if __name__ == "__main__":
    main()  # type: ignore
