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


@click.command()
@click.argument("file_name")
def main(file_name):
    """
    summarise aphc xmas 2024 formSG results

    e.g. summarise_aphc_xmas_formsg.py
    """
    print(f'processing "{file_name}"')

    # read phrase file
    phrases = []
    line_no = 0
    no_participating = 0
    no_skipping = 0
    section_no = dict(car=0, cwr=0, var=0, vm=0, vfp=0, pqa=0, pbps=0, phla=0)
    with open(file_name, "r") as file:
        for line in file:
            line_no += 1
            if line_no > 6:
                line = line.rstrip()
                # print(line)
                tokens = line.split(",")
                l = len(tokens)
                attendance = tokens[4]
                if attendance == "Yes":
                    wish = ",".join(tokens[5:-2])
                    wish = wish.rstrip('"').lstrip('"')
                    name = tokens[3]
                    section = tokens[-1]
                    section_no[section] += 1
                    email = tokens[-2]
                    if not (
                        name == "Adrian" or (name == "Adrian Tan" and section == "PHLA")
                    ):
                        print(f"{name}\t{email}\t{section}\t{wish}")
                        no_participating += 1
                else:
                    no_skipping += 1

    print(f"Summary\n")
    print(f"No respondants:   {no_participating+no_skipping}")
    print(f"No participating: {no_participating}")
    print(f"No skipping:      {no_skipping}")

    print(f"car:      {section_no[car]}")


if __name__ == "__main__":
    main()
