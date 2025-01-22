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
@click.option("-s", "--sample_call_rate_cutoff", default=0.9, help="sample call rate cutoff")
@click.option("-v", "--variant_call_rate_cutoff", default=0.9, help="variant call rate cutoff")
@click.argument("vcf_file")
def main(vcf_file, sample_call_rate_cutoff, variant_call_rate_cutoff):
    """
    Filter VCF file obtained from ddRADSeq.

    e.g. filter_ddradseq_vfc.py
    """
    print("\t{0:<20} :   {1:<10}".format("vcf file", vcf_file))
    print("\t{0:<20} :   {1:<10}".format("minimum sample call rate", sample_call_rate_cutoff))
    print("\t{0:<20} :   {1:<10}".format("minimum variant call rate", variant_call_rate_cutoff))

    # read VCF file, obtain master matrix of data
    data = []
    samples = []
    variants = []
    no_variants = 0
    no_samples = 0
    with open(vcf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    print(line)
                    samples = line.rstrip().split("\t")[9:]
                    no_samples = len(samples)
            else:
                chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")

                no_variants +=1

    print(f"no samples : {no_samples}")
    print(f"no variants : {no_variants}")


    filtered_samples = []
    filtered_variants = []

    change = False
    iter_no = 0

    while change:

        #compute call rates for samples
        for i in filtered_samples:
            if filtered_samples[i]:
                pass

        #filter samples


        #compute call rates and MAF for variants


        #filter variants


        #check for change
        change = False
        iter_no += 1

        #quick summary
        print(f"no samples : {no_samples}")
        print(f"no variants : {no_variants}")

    #compute final call rates and maf, report

    #quick summary
    print(f"no samples : {no_samples}")
    print(f"no variants : {no_variants}")

class Variant(object):
    def __init__(self, id, chrom, pos, ref, alt):
        self.id = id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.genotypes = []

    def add_genotype(self, genotype):
        self.genotypes.append(genotype)

if __name__ == "__main__":
    main() # type: ignore
