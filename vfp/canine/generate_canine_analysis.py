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
        description="Generate canine analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
           usage: generate_canine_analysis -m <make_file> -d <database>
           """
        ),
    )
    parser.add_argument(
        "-m",
        "--make_file",
        help="make file name",
        type=str,
        default="download_refseq_db.mk",
    )
    parser.add_argument(
        "-i",
        "--input_directory",
        help="input directory containing fsa.csv files from GeneticAnalyzer",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o", "--output_directory", help="output directory", type=str, required=True
    )
    args = parser.parse_args()

    for arg in vars(args):
        print("\t{0:<20} :   {1:<10}".format(arg, getattr(args, arg) or ""))

    try:
        os.makedirs(args.output_directory, exist_ok=True)
    except OSError as error:
        print(f"Directory {args.output_directory} cannot be created")

    pg = PipelineGenerator(args.make_file)
    gmp = GeneMapperParser(args.output_directory + os.sep + "dog.vcf")

    #######################
    # aggregate sample files
    #######################
    for subdir, dirs, files in os.walk(args.input_directory):
        print(f"SUBDIR\t{subdir}")
        for file in files:
            if re.search("\.fsa\.csv$", file) != None:
                print(f"\t\t" + subdir + os.sep + file)
                # process one whole STR panel for each dog
                gmp.parse(subdir + os.sep + file)

                gmp.genotype()

                exit()

    ################
    # write make file
    ################
    print("Writing pipeline")
    pg.write()


class GeneMapperParser(object):
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.markers = []
        self.markers.append(VNTR("PEZ02", "17", "1", "GGAA", 104, 145, "Blue"))
        self.markers.append(VNTR("ZFX/Y", "XY", "2", "-", 159, 164, "Blue"))
        self.markers.append(VNTR("PEZ17", "4", "3", "GAAA", 190, 225, "Blue"))
        self.markers.append(
            VNTR("FH2017", "15", "4", "AGGT/AGAT/GATA", 256, 276, "Blue")
        )
        self.markers.append(VNTR("FH2309", "1", "5", "GAAA", 339, 428, "Blue"))
        self.markers.append(VNTR("PEZ05", "12", "6", "TTTA", 92, 117, "Green"))
        self.markers.append(VNTR("FH2001", "23", "7", "GATA", 118, 160, "Green"))
        self.markers.append(VNTR("FH2328", "33", "8", "GAAA", 171, 213, "Green"))
        self.markers.append(VNTR("FH2004", "11", "9", "AAAG", 232, 326, "Green"))
        self.markers.append(VNTR("FH2361", "29", "10", "GAAA", 322, 439, "Green"))
        self.markers.append(VNTR("PEZ21", "2", "11", "AAAT", 83, 103, "Yellow"))
        self.markers.append(VNTR("FH2054", "12", "12", "GATA", 139, 177, "Yellow"))
        self.markers.append(VNTR("FH3377", "3", "13", "GAAAA", 183, 305, "Yellow"))
        self.markers.append(VNTR("FH2107", "3", "14", "GAAA", 291, 426, "Yellow"))
        self.markers.append(VNTR("FH2088", "15", "15", "TTTA/TTCA", 94, 138, "Red"))
        self.markers.append(VNTR("vWF.X", "27", "16", "AGGAAT", 151, 187, "Red"))
        self.markers.append(VNTR("FH2010", "24", "17", "ATGA", 221, 243, "Red"))
        self.markers.append(VNTR("PEZ16", "27", "18", "GAAA", 280, 332, "Red"))
        self.markers.append(VNTR("FH3313", "19", "19", "GAAA", 340, 446, "Red"))
        self.signals = []
        self.sample = []

    def parse(self, file_name):
        file = open(file_name, "r")
        for line in file:
            if not self.is_header(line):
                self.parse_line(line)

    def is_header(self, line):
        return True if line.startswith('"Dye"') else False

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
        self.signals[-1].print_lite()

    def genotype(self):
        for vntr in self.markers:
            vntr.genotype(self.signals)


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
    def __init__(self, locus, chrom, pos1, motif, min_size, max_size, dye):
        self.extended_range = False
        self.alleles = []
        self.gt = ""


class VNTR(object):
    def __init__(self, locus, chrom, pos1, motif, min_size, max_size, dye):
        self.locus = locus
        self.chrom = chrom
        self.pos = pos1
        self.motif = motif
        self.min_size = min_size
        self.max_size = max_size
        self.dye = dye

        # genotype related
        self.samples = []

        self.extended_range = False
        self.alleles = []
        self.gt = ""

    def genotype(self, signals):
        self.print_lite()
        records_exists = False
        for rec in signals:
            if (
                rec.dye == self.dye
                and rec.size >= self.min_size
                and rec.size <= self.max_size
            ):
                records_exists = True
                if rec.height >= 1000 and rec.height < 4000:
                    self.alleles.append(rec.size)
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
                    print("\textended geno: ", end="")
                    rec.print_lite()
            if not records_exists:
                print("Cannot find data")
            else:
                extended_range = True
        if len(self.alleles) == 1:
            self.gt = f"{self.alleles[0]}/{self.alleles[0]}"
        elif len(self.alleles) == 2:
            self.gt = f"{self.alleles[0]}/{self.alleles[1]}"
        else:
            self.gt = "./."

        print(f"\tCALLED GENOTYPE : {self.gt}  ")
        print(f"\t       ALLELES  : {self.alleles}  ")

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


class PipelineGenerator(object):
    def __init__(self, make_file):
        self.make_file = make_file
        self.tgts = []
        self.deps = []
        self.cmds = []

    def add(self, tgt, dep, cmd):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(cmd)

    def print(self):
        print(".DELETE_ON_ERROR:")
        for i in range(len(self.tgts)):
            print(f"{self.tgts[i]} : {self.deps[i]}")
            print(f"\t{self.cmds[i]}")
            print(f"\ttouch {self.tgts[i]}")

    def write(self):
        f = open(self.make_file, "w")
        f.write(".DELETE_ON_ERROR:\n\n")
        f.write("all : ")
        for i in range(len(self.tgts)):
            f.write(f"{self.tgts[i]} ")
        f.write("\n\n")

        for i in range(len(self.tgts)):
            f.write(f"{self.tgts[i]} : {self.deps[i]}\n")
            f.write(f"\t{self.cmds[i]}\n")
            f.write(f"\ttouch {self.tgts[i]}\n\n")
        f.close()


main()
