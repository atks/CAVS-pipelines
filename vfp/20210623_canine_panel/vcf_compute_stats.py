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
        description="Compute stats from dog VCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
           usage: vcf_compute_stats <vcf_file> -o <stats_file>
           """
        ),
    )
    parser.add_argument("vcf_file", type=str, help="input VCF file")
    parser.add_argument(
        "-o", "--output_stats_file", type=str, default="-", help="output stats file"
    )
    args = parser.parse_args()

    for arg in vars(args):
        print("\t{0:<20} :   {1:<10}".format(arg, getattr(args, arg) or ""))

    vcf = open(args.vcf_file, "r")

    SAMPLES = []

    strs_stats = open("strs_call.txt", "w")
    strs_stats.write(f"locus\tcall\tmulti\tmissingness\tweak\n")
    for line in vcf:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                (
                    chrom,
                    pos1,
                    id,
                    ref,
                    alt,
                    qual,
                    filter,
                    info,
                    format,
                    samples,
                ) = line.rstrip().split("\t", 9)

                col = 0
                for sample in samples.split("\t"):
                    SAMPLES.append(Stats(sample))
                    col += 1

        else:
            chrom, pos1, id, ref, alt, qual, filter, info, format, gts = line.split(
                "\t", 9
            )

            gts = gts.split("\t")
            n = 0
            c = 0
            m = 0
            x = 0
            y = 0

            col = 0
            for gt in gts:
                if gt != "./." and gt != "?/?" and gt != "N/D":
                    c += 1
                    SAMPLES[col].c += 1
                if gt == "./.":
                    y += 1
                    SAMPLES[col].y += 1
                if gt == "N/D":
                    x += 1
                    SAMPLES[col].x += 1
                if gt == "?/?":
                    m += 1
                    SAMPLES[col].m += 1
                n += 1
                SAMPLES[col].n += 1
                col += 1
            call_rate = round(c / n * 100, 2)
            amb_rate = round(m / n * 100, 2)
            nd_rate = round(x / n * 100, 2)
            y_rate = round(y / n * 100, 2)

            strs_stats.write(f"{id}\t{call_rate}\t{amb_rate}\t{nd_rate}\t{y_rate}\n")

    samples_stats = open("samples_call.txt", "w")
    samples_stats.write(f"locus\tcall\tmulti\tmissingness\tweak\n")
    for sample in SAMPLES:
        call_rate = round(sample.c / sample.n * 100, 2)
        amb_rate = round(sample.m / sample.n * 100, 2)
        nd_rate = round(sample.x / sample.n * 100, 2)
        y_rate = round(sample.y / sample.n * 100, 2)
        samples_stats.write(
            f"{sample.id}\t{call_rate}\t{amb_rate}\t{nd_rate}\t{y_rate}\n"
        )


class Stats(object):
    def __init__(self, sample_id):
        self.id = sample_id
        self.n = 0
        self.c = 0
        self.m = 0
        self.x = 0
        self.y = 0


main()
