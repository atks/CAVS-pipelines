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
import random

@click.command()
@click.option("-i", "--input_vcf_file", required=True, help="input VCF file", type=str)
@click.option("-o", "--output_vcf_file", required=True, help="output VCF file", type=str)
def main(input_vcf_file, output_vcf_file):
    """
    Generate a panmictic population based on a VCF file.
    
    e.g. simulate_panmictic_pop.py
    """
    print("\t{0:<20} :   {1:<10}".format("input vcf file", input_vcf_file))
    print("\t{0:<20} :   {1:<10}".format("output vcf file", output_vcf_file))

    # read VCF file, obtain master matrix of data
    data = []
    samples = []
    vcf_hdr = ""
    no_variants = 0
    no_samples = 0
    with open(input_vcf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.rstrip().split("\t")[9:]
                    no_samples = len(samples)
                
                vcf_hdr += line
            else:
                chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")
                data.append(Variant(id, chrom, pos, ref, alt, info))
                no_variants +=1

    #write to file
    with open(output_vcf_file, "w") as file:
        file.write(vcf_hdr)
        for i in range(no_variants):
            ac_info = data[i].info.split(";")[4].split("=")[1]
            an_info = data[i].info.split(";")[5].split("=")[1]
            af = float(ac_info) / float(an_info)

            GENOTYPES = ""
            INFO_NS = 0
            INFO_DP = 0
            INFO_AD0 = 0
            INFO_AD1 = 0
            INFO_AF = 0
            INFO_AC = 0
            INFO_AN = 0
            for j in range(no_samples):
                #draw alleles
                gt = 0 
                if random.random() < af:
                    gt += 1
                if random.random() < af:
                    gt += 1

                #GT:DP:AD:GQ:GL
                GT = "./."
                DP = 10
                GQ = 30
                AD0 = 0
                AD1 = 0
                GL = "0,0,0"
                if gt == 0:
                    GT = "0/0"
                    AD0 = 10
                    AD1 = 0
                    GL = "0,30,60"
                elif gt == 1:
                    GT = "0/1"
                    AD0 = 5
                    AD1 = 5 
                    GL = "30,0,30"
                elif gt == 2:
                    GT = "1/1"    
                    AD0 = 0
                    AD1 = 10
                    GL = "60,30,0"

                GENOTYPES +=  f"\t{GT}:{DP}:{AD0},{AD1}:{GQ}:{GL}"
            
                INFO_NS += 1
                INFO_DP += DP
                INFO_AD0 += AD0
                INFO_AD1 += AD1
                INFO_AC += gt
                INFO_AN += 2
                
            INFO_AF = float(INFO_AC) / float(INFO_AN)

            INFO = f"NS={INFO_NS};DP={INFO_DP};AD={INFO_AD0},{INFO_AD1};AF={INFO_AF:.2f};AC={INFO_AC};AN={INFO_AN}"
            file.write(f"{data[i].chrom}\t{data[i].pos}\t{data[i].id}\t{data[i].ref}\t{data[i].alt}\t.\tPASS\t{INFO}\tGT:DP:AD:GQ:GL{GENOTYPES}\n")


      



      
class Variant(object):
    def __init__(self, id, chrom, pos, ref, alt, info):
        self.id = id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

if __name__ == "__main__":
    main() # type: ignore
