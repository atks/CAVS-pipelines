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
import re
import fnmatch


@click.command()
@click.argument("ilm_sample_sheet", required=True)
@click.option("-o", "--output_csv_file", required=True, help="Output Sample Sheet name")
def main(ilm_sample_sheet, output_csv_file):
    """
    Create new sample sheet file without adaptor trimming

    e.g. make_ilm_sample_sheet_without_adaptor_trimming.py SampleSheet.csv -o notrim_SampleSheet.csv
    """

    # Sample sheet setting for adapter trimming
    # [Settings]
    # adapter,CTGTCTCTTATACACATCT

    in_data = False

    os.makedirs(os.path.dirname(output_csv_file), exist_ok=True)
    with open(output_csv_file, "w") as ofile:
        with open(ilm_sample_sheet, "r") as ifile:
            for line in ifile:
                if not line.startswith("adapter"):
                    if line.startswith("[Data]"):
                        in_data = True
                    if in_data and not line.startswith("Sample_ID,"):
                        ofile.write(line.replace("_", "-"))
                    else:
                        ofile.write(line)


if __name__ == "__main__":
    main()  # type: ignore
