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
@click.option("-s", "--sample_file", default=0.9, help="sample file")
@click.argument("vcf_file")
def main(vcf_file, subpop_size, sample_file):
    """
    Compute number of monomorphic sites in subpopulations from a VCF file

    e.g. compute_nonmonomorphic_subset_population.py
    """
    print("\t{0:<20} :   {1:<10}".format("vcf file", vcf_file))
    print("\t{0:<20} :   {1:<10}".format("subpopulation size", subpop_size))
    print("\t{0:<20} :   {1:<10}".format("sample file", sample_file))
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

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
                if no_variants %10000 ==0:
                    print(f"adding variant {no_variants}")
                chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")
                data.append(Variant(id, chrom, pos, ref, alt, genotypes))
                no_variants +=1

    #no change in size
    sample_c = [0]*no_samples
    variant_c = [0]*no_variants
    variant_ac = [0]*no_variants

    #reducing size
    #pointer to subset of samples and variants
    filtered_samples = range(no_samples)
    filtered_variants = range(no_variants)

    new_filtered_samples = []
    new_filtered_variants = []
    finalized_samples = []
    finalized_variants = []

    ts = 0
    tv = 0
    for i in filtered_variants:
        if data[i].ts:
            ts += 1
        else:
            tv += 1
    
    
        print(f"iteration {iter_no}")
        change = False

        if iter_no == 2:
            print("updating sample call rate cut off to 0.9 from second iteration onwards")
            sample_call_rate_cutoff = 0.9
            print("updating variant call rate cut off to 0.9 from second iteration onwards")
            variant_call_rate_cutoff = 0.9

        #compute call rates for samples
        no_filtered_variants = len(filtered_variants)

        for j in filtered_samples:
            sample_c[j] = 0

        for i in filtered_variants:
            for j in filtered_samples:
                if data[i].genotypes[j].gt != -1:
                    sample_c[j] += 1

        #filter samples
        with open(f"{output_dir}/sample_call_rate_iter_{iter_no}.txt", "w") as file:
            new_filtered_samples.clear()
            file.write("#sample\tsample_call_rate\n")
            for j in filtered_samples:
                sample_call_rate = float(sample_c[j])/no_filtered_variants
                file.write(f"{samples[j]}\t{sample_call_rate}\n")
                #print(f"sample call rate: {sample_call_rate}")
                if sample_call_rate >= sample_call_rate_cutoff:
                    new_filtered_samples.append(j)

        change = len(new_filtered_samples) != len(filtered_samples)
        filtered_samples = new_filtered_samples.copy()
        no_filtered_samples = len(filtered_samples)

        print(f"no filtered samples : {no_filtered_samples}")

        with open(f"{output_dir}/snp_call_rate_iter_{iter_no}.txt", "w") as file:
            with open(f"{output_dir}/maf_iter_{iter_no}.txt", "w") as maf_file:
                file.write("#variant_call_rate\n")
                maf_file.write("#variant_maf\n")
                for i in filtered_variants:
                    #compute variant call rates
                    variant_c[i] = 0
                    variant_ac[i] = 0
                    for j in filtered_samples:
                        if data[i].genotypes[j].gt != -1:
                            variant_c[i] += 1
                            variant_ac[i] += data[i].genotypes[j].gt

                    #filter variants
                    variant_call_rate = float(variant_c[i])/no_filtered_samples
                    if variant_c[i] == 0:
                        variant_maf = 0
                    else:
                        variant_af = float(variant_ac[i])/(variant_c[i]*2)
                    variant_maf = min(variant_af, 1-variant_af)
                    #print(f"variant call rate: {variant_call_rate: .2f} | {variant_maf: .2f}")
                    file.write(f"{variant_call_rate}\n")
                    maf_file.write(f"{variant_maf}\n")
                    if variant_call_rate >= variant_call_rate_cutoff and variant_maf >= variant_maf_cutoff:
                        new_filtered_variants.append(i)
                        if data[i].ts:
                            ts += 1
                        else:
                            tv += 1

    
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
