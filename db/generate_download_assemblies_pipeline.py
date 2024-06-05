#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2023 Adrian Tan <adrian_tan@nparks.gov.sg>
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

import click
import os
import re
import subprocess
import sys


@click.command()
@click.option(
    "-m",
    "--make_file",
    default="download_assemblies.mk",
    show_default=True,
    help="make file name",
)
@click.option("-s", "--species_file", required=True, help="list of species information")
@click.option(
    "-o",
    "--output_dir",
    default=os.getcwd(),
    show_default=True,
    help="output directory",
)
def main(make_file, species_file, output_dir):
    """
    Download assemblies

    e.g.  generate_download_assemblies_pipeline -s species.txt -m download_assemblies.mk
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("species file", species_file))
    print("\t{0:<20} :   {1:<10}".format("output dir", output_dir))

    animals = []
    with open(species_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                id, binomial_name, acc_id = line.rstrip().split("\t")
                if acc_id != "n/a":
                    animals.append(Animal(id, binomial_name, acc_id))

    print(f"Will download {len(animals)} genomes from NCBI datasets")

    # generate make file
    print("Generating pipeline")
    pg = PipelineGenerator(make_file)

    # datasets download genome accession GCF_014570535.1 --filename pangolin_genome.zip

    datasets = "/usr/local/ncbi-cli-14.13.0/datasets"

    for animal in animals:
        output_file = f'{animal.binomial_name.replace(" " , "_").lower()}.zip'
        tgt = f"{output_file}.OK"
        dep = ""
        cmd = f"{datasets} download genome accession {animal.acc_id} --filename {output_file}"
        pg.add(tgt, dep, cmd)

    # clean files
    cmd = f"rm -fr {output_dir}/*.OK  {output_dir}/*.err"
    pg.add_clean(cmd)

    # write make file
    print("Writing pipeline")
    pg.write()


class PipelineGenerator(object):
    def __init__(self, make_file):
        self.make_file = make_file
        self.tgts = []
        self.deps = []
        self.cmds = []
        self.clean_cmd = ""

    def add(self, tgt, dep, cmd):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(cmd)

    def add_clean(self, cmd):
        self.clean_cmd = cmd

    def write(self):
        with open(self.make_file, "w") as f:
            f.write("SHELL:=/bin/bash\n")
            f.write(".DELETE_ON_ERROR:\n\n")
            f.write("all : ")
            for i in range(len(self.tgts)):
                f.write(f"{self.tgts[i]} ")
            f.write("\n\n")

            for i in range(len(self.tgts)):
                f.write(f"{self.tgts[i]} : {self.deps[i]}\n")
                f.write(f"\t{self.cmds[i]}\n")
                f.write(f"\ttouch {self.tgts[i]}\n\n")

            if self.clean_cmd != "":
                f.write(f"clean : \n")
                f.write(f"\t{self.clean_cmd}\n")


class Animal(object):
    def __init__(self):
        self.id = ""
        self.binomial_name = ""
        self.acc_id = ""

    def __init__(self, id, binomial_name, acc_id):
        self.id = id
        self.binomial_name = binomial_name
        self.acc_id = acc_id

    def print(self):
        print(f"id             : {self.id}")
        print(f"binomial_name  : {self.binomial_name}")
        print(f"acc_id         : {self.acc_id}")


if __name__ == "__main__":
    main()
