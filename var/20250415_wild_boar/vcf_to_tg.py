#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2025 Adrian Tan <adrian_tan@nparks.gov.sg>
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

@click.command()
@click.argument("vcf_file")
@click.option("-o", "--output_dir", required=True, help="output directory")
def main(vcf_file, output_dir):
    """
    Convert VCF file to tg format.

    e.g. vcf_to_tg.py
    """
    print("\t{0:<20} :   {1:<10}".format("vcf file", vcf_file))
    print("\t{0:<20} :   {1:<10}".format("output directory", output_dir))


    out_tg_file = os.path.join(output_dir, os.path.basename(vcf_file).replace(".vcf", ".tg"))

    # read VCF file, obtain master matrix of data
    data = []
    samples = []
    no_variants = 0
    no_samples = 0
    with open(out_tg_file, "w") as tg_file:
        with open(vcf_file, "r") as file:
            for line in file:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        samples = "\t".join(line.rstrip().split("\t")[9:])
                        tg_file.write(f"snp-id\t{samples}\n")
                else:
                    chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")
                    tg_file.write(f"{id}")
                    for geno in genotypes:
                        if geno.startswith("./."):
                            tg_file.write("\t-1")
                        else:
                            geno = geno.split(":")[0].split("/")
                            tg_file.write(f"\t{int(geno[0])+int(geno[1])}")
                    tg_file.write("\n")


if __name__ == "__main__":
    main() # type: ignore
