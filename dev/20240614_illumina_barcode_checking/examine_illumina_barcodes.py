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

import os
import click
import sys

@click.command()
@click.argument("fastqc_dir", required=True)
@click.argument("sample_id", required=True)
def main(fastqc_dir, sample_id):
    """
    Takes in an FASTQC analysis directory and sample ID and inspects the fastqc_result directory to ensure run completion.
    This is required as fastqc does not always properly output a failure status.

    e.g. check_fastqc_run.py /home/atks/hts/illu14/illu14/analysis/10_staphylococcus_pseudintermedius_416-18/fastqc 10_staphylococcus_pseudintermedius_416-18_R1
    """

    all_good = True

    file = f"{fastqc_dir}/{sample_id}_fastqc.html"
    if not os.path.isfile(file):
        print(f"{file} missing", file=sys.stderr)
        all_good = False
    file = f"{fastqc_dir}/{sample_id}_fastqc.zip"
    if not os.path.isfile(file):
        print(f"{file} missing", file=sys.stderr)
        all_good = False

    if all_good:
        exit(0)
    else:
        exit(1)

if __name__ == "__main__":
    main()  # type: ignore
