#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2021 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import argparse
import textwrap
import re


def main():
    cwd = os.getcwd()

    parser = argparse.ArgumentParser(
        description="Extract STR plots",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
           usage: extract_str_plots <sample_id> -marker <marker> -v <input_vcf_file> -s <sample_phenotype_file> -f <sample_file_list>

           <sample_id> maybe be multiple IDs delimited by a comma or all

           output will be  <sample-id>_<marker>.pdf
           e.g. BW002_ZFXY.pdf
                BW002_blue.pdf
                BW002_all.pdf
           """
        ),
    )
    parser.add_argument(
        "sample_id",
        type=str,
        help="input directory containing fsa.csv files from GeneticAnalyzer",
    )
    parser.add_argument(
        "-v",
        "--vcf_file",
        type=str,
        required=True,
        help="input directory containing fsa.csv files from GeneticAnalyzer",
    )
    parser.add_argument(
        "-m", "--marker", type=str, default="all", help="marker or dye channel or all"
    )
    parser.add_argument(
        "-s",
        "--sample_phenotype_file",
        type=str,
        required=True,
        help="sample file containing phenotype",
    )
    parser.add_argument(
        "-f",
        "--sample_file_list",
        type=str,
        required=True,
        help="sample file containing file locations",
    )
    args = parser.parse_args()

    for arg in vars(args):
        print("\t{0:<24} :   {1:<10}".format(arg, getattr(args, arg) or ""))


main()
