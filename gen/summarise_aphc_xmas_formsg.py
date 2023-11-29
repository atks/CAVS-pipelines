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

import os
import click
import random
import math
import sys


@click.command()
@click.argument("file_name")
def main(file_name):
    """
    summarise aphc xmas 2024 formSG results

    e.g. summarise_aphc_xmas_formsg.py
    """
    # print(f'processing "{file_name}"')

    print(f"#giving\t#name\t#email\t#section\t#wish")

    # read phrase file
    phrases = []
    line_no = 0
    no_participating = 0
    no_skipping = 0
    section_no = {
        "CAR": 0,
        "CWR": 0,
        "VAR": 0,
        "VM": 0,
        "VFP": 0,
        "PQA": 0,
        "PBPS": 0,
        "PHLA": 0,
    }

    with open(file_name, "r") as file:
        for line in file:
            line_no += 1
            if line_no > 6:
                line = line.rstrip()
                # print(line)
                tokens = line.split(",")
                l = len(tokens)
                attendance = tokens[4]
                name = tokens[3]
                section = tokens[-1]
                if attendance == "Yes":
                    if not (
                        name == "Adrian" or (name == "Adrian Tan" and section == "PHLA")
                    ):
                        wish = ",".join(tokens[5:-2])
                        wish = wish.rstrip('"').lstrip('"')
                        section_no[section] += 1
                        email = tokens[-2]
                        print(f"YES\t{name}\t{email}\t{section}\t{wish}")
                        no_participating += 1
                else:
                    no_skipping += 1
                    if name == "Ng Oi Wing":
                        section = "VAR"
                    if name == "Amy Chan":
                        section = "VFP"
                    print(f"NO\t{name}\tn/a\t{section}\tn/a")

    # randomise assignments

    print(f"\n\n", file=sys.stderr)
    print(f"Summary", file=sys.stderr)
    print(f"No respondants:   {no_participating+no_skipping}", file=sys.stderr)
    print(f"No participating: {no_participating}", file=sys.stderr)
    print(f"No skipping:      {no_skipping}", file=sys.stderr)

    print(f"car:      {section_no['CAR']}", file=sys.stderr)
    print(f"cwr:      {section_no['CWR']}", file=sys.stderr)
    print(f"var:      {section_no['VAR']}", file=sys.stderr)
    print(f"vm:       {section_no['VM']}", file=sys.stderr)
    print(f"vfp:      {section_no['VFP']}", file=sys.stderr)
    print(f"pqa:      {section_no['PQA']}", file=sys.stderr)
    print(f"pbps:     {section_no['PBPS']}", file=sys.stderr)
    print(f"phla:     {section_no['PHLA']}", file=sys.stderr)


if __name__ == "__main__":
    main()
