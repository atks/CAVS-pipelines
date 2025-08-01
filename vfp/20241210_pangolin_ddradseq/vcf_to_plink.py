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
@click.option("-s", "--sample_file", required=True, help="sample file")
@click.option("-o", "--output_plink_file_base_name", default="plink", help="plink file base name - two files .ped and .map")
def main(vcf_file, sample_file, output_plink_file_base_name):
    """
    Convert VCF file to plink format.

    e.g. vcf_to_plink.py
    """
    print("\t{0:<20} :   {1:<10}".format("vcf file                   ", vcf_file))
    print("\t{0:<20} :   {1:<10}".format("sample file                ", sample_file))
    print("\t{0:<20} :   {1:<10}".format("output plink file base name", output_plink_file_base_name))

    # Family ID
    # Individual ID
    # Paternal ID
    # Maternal ID
    # Sex (1=male; 2=female; other=unknown)
    # Phenotype
    fam_file = f"{output_plink_file_base_name}.fam"
    
    # chromosome (1-22, X, Y or 0 if unplaced)
    # rs# or snp identifier
    # Genetic distance (morgans)
    # Base-pair position (bp units)
    map_file = f"{output_plink_file_base_name}.map"

    # family ID
    # individual ID
    # snp ID
    # allele 1 of this genotype
    # allele 2 of this genotype             
    lgen_file = f"{output_plink_file_base_name}.lgen"
   
    no_variants = 0
    no_samples = 0
    
    SAMPLES = {}

    with open(sample_file, "r") as sa:
        for line in sa:
            if not line.startswith("#"):
                sample_id, latitude, longitude, sex = line.rstrip().split("\t")
                SAMPLES[sample_id] = Sample(sample_id, latitude, longitude, sex)

    print("{0:<20} :   {1:<10}".format("no. of samples in sample file", len(SAMPLES)))

    samples = []
    families = []
    with open(lgen_file, "w") as lgen:
        with open(map_file, "w") as map:    
            with open(vcf_file, "r") as vcf:
                for line in vcf:
                    if line.startswith("#"):
                        if line.startswith("#CHROM"):
                            samples = line.rstrip().split("\t")[9:]
                            no_samples = len(samples)
                            print("{0:<20} :   {1:<10}".format("no. of samples in vcf file", no_samples))
                            with open(fam_file, "w") as fam:
                                for idx, sample in enumerate(samples):
                                    # Family ID
                                    # Individual ID
                                    # Paternal ID
                                    # Maternal ID
                                    # Sex (1=male; 2=female; other=unknown)
                                    # Phenotype
                                    family_id = f"FAM{idx+1:2d}"
                                    sample_id = sample
                                    families.append(family_id)
                                    sex_code = -1
                                    if sample_id in SAMPLES:
                                        if SAMPLES[sample_id].sex == "M":
                                            sex_code = 1
                                        elif SAMPLES[sample_id].sex == "F":
                                            sex_code = 2
                                    fam.write(f"{family_id}\t{sample_id}\t-1\t-1\t{sex_code}\t0\n")
                    else:
                        print(f"Processing variant record line: {line.rstrip()}")
                        no_variants += 1
                        if no_variants % 1000 == 0:
                            print(f"\r{no_variants} variants processed", end="")
                    
                        chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")
                        
                        # chromosome (1-22, X, Y or 0 if unplaced)
                        # rs# or snp identifier
                        # Genetic distance (morgans) 
                        # Base-pair position (bp units)
                        map.write(f"{chrom} {id} 0 {pos}\n")
                        
                        # family ID
                        # individual ID
                        # snp ID
                        # allele 1 of this genotype
                        # allele 2 of this genotype 
                        for idx, genotype in enumerate(genotypes):
                            sample_id = samples[idx]
                            if genotype == "./.":
                                lgen.write(f"{families[idx]}\t{samples[idx]}\t{id}\t0 0\n")
                            else:
                                gt, dp, ad, gq, gl = genotype.split(":")
                                if gt == "./.":
                                   lgen.write(f"{families[idx]}\t{samples[idx]}\t{id}\t0 0\n")
                                else:
                                    allele1 = 0
                                    allele2 = 0
                                    if gt == "0/0":
                                        allele1 = ref
                                        allele2 = ref
                                    elif gt == "0/1":
                                        allele1 = ref
                                        allele2 = alt
                                    elif gt == "1/1":
                                        allele1 = alt
                                        allele2 = alt
                                    lgen.write(f"{families[idx]}\t{samples[idx]}\t{id}\t{allele1}\t{allele2}\n")

class Sample:
    def __init__(self, sample_id, latitude, longitude, sex):
        self.sample_id = sample_id
        self.latitude = latitude
        self.longitude = longitude
        self.sex = sex

if __name__ == "__main__":
    main() # type: ignore