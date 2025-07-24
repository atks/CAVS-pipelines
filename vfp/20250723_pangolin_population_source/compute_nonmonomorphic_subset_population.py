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
@click.option("-n", "--subpop_size", default=5, help="subpopulation size", type=int)
@click.option("-s", "--sample_file", help="sample file", type=str)
@click.argument("vcf_file")
def main(vcf_file, subpop_size, sample_file):
    """
    Compute number of monomorphic sites in subpopulations from a VCF file

    e.g. compute_nonmonomorphic_subset_population.py
    """
    print("\t{0:<20} :   {1:<10}".format("vcf file", vcf_file))
    print("\t{0:<20} :   {1:<10}".format("subpopulation size", subpop_size))
    print("\t{0:<20} :   {1:<10}".format("sample file", sample_file))
    
    samples = {}

    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                sample_id, rest = line.strip().split("\t", 1)
                samples[sample_id] = 1  

    print(f"Number of samples: {len(samples)}")

    # read VCF file, obtain master matrix of data
    data = []
    subsamples_idx = []
    no_variants = 0
    with open(vcf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    sample_ids = line.rstrip().split("\t")[9:]
                    for idx,id in enumerate(sample_ids):
                        #print(f"sample {idx}: {id}")
                        if id in samples:
                            subsamples_idx.append(idx)
                    # print(f"INITsubsamples_idx length: {len(subsamples_idx)}")        
            else:
                chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")
                subsample_genotypes = []
                for idx in subsamples_idx:
                    subsample_genotypes.append(genotypes[idx])
                data.append(Variant(id, chrom, pos, ref, alt, subsample_genotypes))
                no_variants +=1

    print(f"Number of variants: {no_variants}")

    c = [[0],[1]]
    d = [l + [2] for l in c]
    print(f"c: {c}")
    print(f"d: {d}")

    # print(f"data length: {len(data)}")
    # print(f"subsamples_idx length: {len(subsamples_idx)}")
    # permutate all subsamples
    # for variant in data:
    #   variant.print()
    def subsets(n, k):
        print(f"recursive step n={n} k={k}")
        if k > n:
            return []
        elif k==n:
            return [list(range(n))]
        elif k == 0:
            return []
        else:            
            return subsets(n-1, k) + [l + [n-1] for l in subsets(n-1, k-1)]

    subsets_list = subsets(2, 1)
    print(f"subsets_list length: {len(subsets_list)}") 
    print(subsets_list)
    # for l in subsets_list:
    #     print(l)
    
class Variant(object):
    def __init__(self, id, chrom, pos, ref, alt, genotypes):
        self.id = id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
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
            #print(genotype)
            gt, dp, ad, gq, gl = genotype.split(":")
            if gt == "./.":
                self.genotypes.append(Genotype(-1, -1, -1, -1, "-1,-1,-1"))
            else:
                if gt == "0/0":
                    gt = 0
                elif gt == "0/1":
                    gt = 1
                elif gt == "1/1":
                    gt = 2
                self.genotypes.append(Genotype(gt, int(dp), ad, int(gq), gl))
    
    def print(self):
        print(f"{self.id}", end="")
        print(f"\t{self.chrom}", end="")
        print(f"\t{self.pos}", end="")
        print(f"\t{self.ref}", end="")
        print(f"\t{self.alt}", end="")
        for g in self.genotypes:
            g.print()
        print("")    
            
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

    def print(self):
        print(f"\t{self.gt}", end="")

if __name__ == "__main__":
    main() # type: ignore
