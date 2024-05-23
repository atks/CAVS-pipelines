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
        description="Genotype canine VNTRs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
           usage: abi_vntr_genotype_caller -s <sample_file> -o <output_vcf_file>
           """
        ),
    )
    parser.add_argument(
        "-s",
        "--sample_file",
        type=str,
        required=True,
        help="input directory containing fsa.csv files from GeneticAnalyzer",
    )
    parser.add_argument(
        "-o", "--output_vcf_file", type=str, required=True, help="output VCF file"
    )
    args = parser.parse_args()

    for arg in vars(args):
        print("\t{0:<20} :   {1:<10}".format(arg, getattr(args, arg) or ""))

    gmp = GeneMapperParser(args.output_vcf_file)

    #######################
    # aggregate sample files
    #######################
    sa_file = open(args.sample_file, "r")
    for line in sa_file:
        if not line.startswith("#"):
            id, file = line.rstrip().split("\t")
            file = file.strip('"')
            print(f"{id} => {file}")
            gmp.parse(file)
            gmp.genotype(id)
            gmp.write_to_vcf()


class GeneMapperParser(object):
    def __init__(self, vcf_file_name):
        self.vcf_file_name = vcf_file_name
        self.markers = []
        self.markers.append(
            VNTR("FH2309", "1", 85772974, 85773377, "GAAA", 339, 428, "Blue")
        )
        self.markers.append(
            VNTR("PEZ21", "2", 36438658, 36438751, "AAAT", 83, 103, "Yellow")
        )
        self.markers.append(
            VNTR("FH3377", "3", 78748898, 78749090, "GAAAA", 183, 305, "Yellow")
        )
        self.markers.append(
            VNTR("FH2107", "3", 83830247, 83830574, "GAAA", 291, 426, "Yellow")
        )
        self.markers.append(
            VNTR("PEZ17", "4", 71904833, 71905038, "GAAA", 190, 225, "Blue")
        )
        self.markers.append(
            VNTR("FH2004", "11", 32161381, 32161621, "AAAG", 232, 326, "Green")
        )
        self.markers.append(
            VNTR("FH2054", "12", 37914504, 37914739, "GATA", 139, 177, "Yellow")
        )
        self.markers.append(
            VNTR("PEZ05", "12", 60326434, 60326541, "TTTA", 92, 117, "Green")
        )
        self.markers.append(
            VNTR("FH2017", "15", 37914470, 37914741, "AGGT/AGAT/GATA", 256, 276, "Blue")
        )
        self.markers.append(
            VNTR("FH2088", "15", 53905651, 53905779, "TTTA/TTCA", 94, 138, "Red")
        )
        self.markers.append(
            VNTR("PEZ02", "17", 13276076, 13276209, "GGAA", 104, 145, "Blue")
        )
        self.markers.append(
            VNTR("FH3313", "19", 24606038, 24606459, "GAAA", 340, 446, "Red")
        )
        self.markers.append(
            VNTR("FH2001", "23", 50961325, 50961475, "GATA", 118, 160, "Green")
        )
        self.markers.append(
            VNTR("FH2010", "24", 5196383, 5196605, "ATGA", 221, 243, "Red")
        )
        self.markers.append(
            VNTR("PEZ16", "27", 10305692, 10305995, "GAAA", 280, 332, "Red")
        )
        self.markers.append(
            VNTR("vWF.X", "27", 41977918, 41978074, "AGGAAT", 151, 187, "Red")
        )
        self.markers.append(
            VNTR("FH2361", "29", 19723594, 19723782, "GAAA", 322, 439, "Green")
        )
        self.markers.append(
            VNTR("FH2328", "33", 19158127, 19158477, "GAAA", 171, 213, "Green")
        )
        self.markers.append(
            VNTR(
                "ZFX/Y", "XY", 19690330, 19748404, "GTTTTTGCCAGTCTGA", 159, 164, "Blue"
            )
        )
        self.signals = []
        self.sample = []

    def parse(self, file_name):
        file = open(file_name, "r")
        self.signals.clear()
        for line in file:
            if not self.is_header(line):
                self.parse_line(line)

    def is_header(self, line):
        return True if line.startswith('"Dye"') or line.startswith("Dye") else False

    def is_informative(self, tokens):
        for val in tokens:
            if val == " ":
                return False
        return True

    def parse_line(self, line):
        tokens = line.strip().split(",")
        tokens = list(map(lambda x: x.strip('"'), tokens))
        if self.is_informative(tokens):
            self.signals.append(
                Record(
                    tokens[0],
                    tokens[1],
                    tokens[2],
                    tokens[3],
                    tokens[4],
                    tokens[5],
                    tokens[6],
                    tokens[7],
                    tokens[8],
                    tokens[9],
                    tokens[10],
                    tokens[11],
                    tokens[12],
                )
            )
        # self.signals[-1].print_lite()

    def genotype(self, id):
        for vntr in self.markers:
            vntr.genotype(id, self.signals)

    def write_to_vcf(self):
        vcf_file = open(self.vcf_file_name, "w")
        vcf_file.write("##VCF4.3\n")
        vcf_file.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for indiv in self.markers[0].samples:
            vcf_file.write(f"\t{indiv.id}")
        vcf_file.write("\n")
        for vntr in self.markers:
            vcf_file.write(
                f"{vntr.chrom}\t{vntr.pos1}\t{vntr.locus}\t{vntr.motif}\t<STR>\t.\t.\tDYE={vntr.dye};MIN_SIZE={vntr.min_size};MAX_SIZE={vntr.max_size}\tGT"
            )
            for indiv in vntr.samples:
                vcf_file.write(f"\t{indiv.gt}")
            vcf_file.write("\n")
        vcf_file.close()


# "Dye","Sample File Name","Size","Height","Area in Point","Area in BP","Data Point","Begin Point","Begin BP","End Point","End BP","Width in Point","Width in BP"
# CHROM POS ID REF ALT QUAL FILTER INFO
# 20 2 . TCG T . PASS DP=100
# GT
##genotype
##raw data from genemapper?
class Record(object):
    def __init__(
        self,
        dye,
        sample_file_name,
        size,
        height,
        area_in_point,
        area_in_bp,
        data_point,
        begin_point,
        begin_bp,
        end_point,
        end_bp,
        width_in_point,
        width_in_bp,
    ):
        self.dye = dye
        self.sample_file_name = sample_file_name
        self.size = float(size)
        self.height = float(height)
        self.area_in_point = float(area_in_point)
        self.area_in_bp = float(area_in_bp)
        self.data_point = float(data_point)
        self.begin_point = float(begin_point)
        self.begin_bp = float(begin_bp)
        self.end_point = float(end_point)
        self.end_bp = float(end_bp)
        self.width_in_point = float(width_in_point)
        self.width_in_bp = float(width_in_bp)

    def print(self):
        print(
            f"{self.dye}|{self.sample_file_name}|{self.size}|{self.height}|"
            + f"{self.area_in_point}|{self.area_in_bp}|{self.data_point}|"
            + f"{self.begin_point}|{self.begin_bp}|{self.end_point}|{self.end_bp}|"
            + f"{self.width_in_point}|{self.width_in_bp}"
        )

    def print_lite(self):
        print(f"{self.dye}|{self.size}|{self.height}")


class Individual(object):
    def __init__(self, id):
        self.id = id
        self.extended_range = False
        self.gt_info = ""
        self.alleles = []
        self.gt = ""


class VNTR(object):
    def __init__(self, locus, chrom, pos1, end1, motif, min_size, max_size, dye):
        self.locus = locus
        self.chrom = chrom
        self.pos1 = pos1
        self.end1 = end1
        self.motif = motif
        self.min_size = min_size
        self.max_size = max_size
        self.dye = dye

        # genotype related
        self.samples = []

    #        self.extended_range = False
    #        self.alleles = []
    #        self.gt = ''

    def genotype(self, id, signals):
        indiv = Individual(id)

        self.print_lite()
        #        for rec in signals:
        #            if rec.dye == self.dye:
        #                rec.print_lite()
        print("=====")

        records_exists = False
        for rec in signals:
            if (
                rec.dye == self.dye
                and rec.size >= self.min_size
                and rec.size <= self.max_size
            ):
                records_exists = True
                #                if rec.height >=1000 and rec.height < 4000:
                if rec.height >= 1000:
                    indiv.alleles.append(round(rec.size, 1))
                print("\tgeno: ", end="")
                rec.print_lite()
        if not records_exists:
            records_exists = False
            for rec in signals:
                if (
                    rec.dye == self.dye
                    and rec.size >= self.min_size - 2
                    and rec.size <= self.max_size + 2
                ):
                    records_exists = True
                    if rec.height >= 1000:
                        indiv.alleles.append(round(rec.size, 1))
                    print("\textended geno: ", end="")
                    rec.print_lite()
            if not records_exists:
                indiv.gt_info = "ND"
                print("Cannot find data")
            else:
                extended_range = True
        if len(indiv.alleles) == 1:
            if self.locus == "ZFX/Y":
                indiv.gt = f"0/0"
            else:
                indiv.gt = f"{indiv.alleles[0]}/{indiv.alleles[0]}"
        elif len(indiv.alleles) == 2:
            if self.locus == "ZFX/Y":
                indiv.gt = f"0/1"
            else:
                indiv.gt = f"{indiv.alleles[0]}/{indiv.alleles[1]}"
        elif len(indiv.alleles) > 2:
            print("More than 2 alleles")
            indiv.gt = f"?/?"
        else:
            indiv.gt = "./."

        if indiv.gt_info == "ND":
            indiv.gt = "N/D"

        self.samples.append(indiv)

        print(f"\tCALLED GENOTYPE : {indiv.gt}  ")
        print(f"\t       ALLELES  : {indiv.alleles}  ")

    def print(self):
        print(f"locus: {self.locus}")
        print(f"chrom: {self.chrom}")
        print(f"motif: {self.motif}")
        print(f"pos1: {self.pos1}")
        print(f"motif: {self.motif}")
        print(f"min_size: {self.min_size}")
        print(f"max_size: {self.max_size}")
        print(f"dye: {self.dye}")

    def print_lite(self):
        print(f"locus: {self.locus}")
        print(f"{self.dye} : {self.min_size}-{self.max_size}")


main()
