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
from random import sample

@click.command()
@click.argument("vcf_file")
@click.option("-s", "--no_subsets", default=100, help="number of subsets to generate")
@click.option("-n", "--subset_size", default=2000, help="size of subsets")
def main(vcf_file, no_subsets, subset_size):
    """
    Convert VCF file to relpair format.

    e.g. vcf_to_relpair.py
    """
    print("\t{0:<20} :   {1:<10}".format("vcf file", vcf_file))
    print("\t{0:<20} :   {1:<10}".format("no subsets", no_subsets))
    print("\t{0:<20} :   {1:<10}".format("subset size", subset_size))

    # read VCF file, obtain master matrix of data
    data = []
    samples = []
    vcf_hdr = ""
    no_variants = 0
    no_samples = 0
    with open(vcf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.rstrip().split("\t")[9:]
                    no_samples = len(samples)
                else:
                    vcf_hdr += line
            else:
                chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")
                ns, dp, ad, af = info.split(";")
                h,f = af.split("=")
                data.append(Variant(id, chrom, pos, ref, alt, float(f), genotypes))
                no_variants +=1

    for i in range(no_subsets):
        relpair_loc_file = vcf_file.replace(".vcf", ".loc")
        relpair_ped_file = vcf_file.replace(".vcf", ".ped")
        relpair_ctl_file = vcf_file.replace(".vcf", ".ctl")
        relpair_out_file = vcf_file.replace(".vcf", ".out")

        #subsample
        subset_indices = sorted(sample(range(no_variants), subset_size))

        with open(relpair_ctl_file, "w") as ctl_file:
            ctl_file.write(f"{relpair_loc_file}\n")
            ctl_file.write(f"{relpair_ped_file}\n")
            ctl_file.write(f"{relpair_out_file}\n")
            ctl_file.write("all\n")
            ctl_file.write("n\n")
            ctl_file.write("n\n")
            ctl_file.write("F\n")
            ctl_file.write("M\n")
            ctl_file.write("2\n")
            ctl_file.write("0.01\n")
            ctl_file.write("1\n")
            ctl_file.write("10.0\n")

        with open(relpair_loc_file, "w") as loc_file:
            #write locus header
            for i in subset_indices:
                loc_file.write(f"MARKER{i} AUTOSOME 2 0 0.0\n")
                loc_file.write(f"R {1-data[i].af:.2f}\n")
                loc_file.write(f"A {data[i].af}\n")

        with open(relpair_ped_file, "w") as ped_file:
            ped_file.write(f"(I2,1X,A8)\n")
            ped_file.write(f"(3A8,2A1,A3,{subset_size-1}(1X,A3))\n")
            ped_file.write(f"{no_samples} FAMILY1\n")
            for j in range(no_samples):
                ped_file.write(f"{samples[j]}                ")
                sample_line = "F"
                for i in subset_indices:
                    gt = data[i].genotypes[j].gt
                    if gt == -1:
                        sample_line += "    "
                    elif gt == 0:
                        sample_line += " R/R"
                    elif gt == 1:
                        sample_line += " R/A"
                    elif gt == 2:
                        sample_line += " A/A"
                ped_file.write(f"{sample_line}\n")



class Variant(object):
    def __init__(self, id, chrom, pos, ref, alt, af, genotypes):
        self.id = id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.af = af
        self.genotypes = []
        for g in genotypes:
            self.add_genotype(g)
        self.ts = self.is_ts(ref,alt)
        #print(f"tv: {self.tv} | {ref} | {alt}")

    def is_ts(self, ref, alt):
        return (ref == "A" and alt == "G") or \
        (ref == "G" and alt == "A") or \
        (ref == "C" and alt == "T") or \
        (ref == "T" and alt == "C")

    def add_genotype(self, genotype):
        if genotype == "./.":
            self.genotypes.append(Genotype(-1, -1, -1, -1, "-1,-1,-1"))
        else:
            gt, dp, ad, gq, gl = genotype.split(":")
            if gt == "./.":
                print(genotype)
                self.genotypes.append(Genotype(-1, -1, -1, -1, "-1,-1,-1"))
            else:
                if gt == "0/0":
                    gt = 0
                elif gt == "0/1":
                    gt = 1
                elif gt == "1/1":
                    gt = 2
                self.genotypes.append(Genotype(gt, int(dp), ad, int(gq), gl))

class Genotype(object):
    def __init__(self, gt, dp, ad, gq, gl):
        self.gt = gt
        self.dp = dp
        if ad == -1:
            self.ad = -1
        else:
            self.ad = int(ad.split(",")[1])
        self.gq = gq
        self.gl = gl.split(",")

if __name__ == "__main__":
    main() # type: ignore