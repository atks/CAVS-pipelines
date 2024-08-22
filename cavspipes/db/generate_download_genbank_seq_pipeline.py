#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2022 Adrian Tan <adrian_tan@nparks.gov.sg>
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
    default="download_genbank_sequences.mk",
    show_default=True,
    help="make file name",
)
@click.option(
    "-s", "--sequence_id_file", required=True, help="list of sequence accension IDs"
)
@click.option(
    "-t",
    "--download_type",
    default="fasta",
    type=click.Choice(["fasta", "genbank"]),
    help="download type - fasta or genbank",
)
@click.option(
    "-o",
    "--output_dir",
    default=os.getcwd(),
    show_default=True,
    help="output directory",
)
def main(make_file, sequence_id_file, download_type, output_dir):
    """
    Download genbank sequences

    e.g.  generate_download_genbank_seq_pipeline -s id.txt -m download_gb_seq.mk -o /home/atks/downloads
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("sequence ID file", sequence_id_file))
    print("\t{0:<20} :   {1:<10}".format("download type", download_type))
    print("\t{0:<20} :   {1:<10}".format("output dir", output_dir))

    ext = "fasta"
    if download_type == "genbank":
        ext = "genbank"

    release_number = subprocess.run(
        ["curl", f"https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number"],
        text=True,
        capture_output=True,
    ).stdout.rstrip()

    ids = []
    with open(sequence_id_file, "r") as file:
        for line in file:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            if "-" in line:
                print(line, end="\t")
                result = re.search(r"(\D+)(\d+)-(\D+)(\d+)", line)
                if result is None:
                    print("invalid range")
                    exit()
                prefix1 = result.group(1)
                num1 = result.group(2)
                prefix2 = result.group(3)
                num2 = result.group(4)
                if prefix1 != prefix2:
                    print("prefix not the same")
                    exit()
                if len(num1) != len(num2):
                    print("numeric length not the same")
                    exit()
                prefix = prefix1
                nlen = len(num1)
                num1 = int(num1)
                num2 = int(num2)
                if num1 >= num2:
                    print("num1 > num2")
                    exit()
                # generate sequence
                no_ids = 0
                for i in range(num1, num2 + 1):
                    str_i = f"{i}"
                    id = f"{prefix}{str_i.zfill(nlen)}"
                    ids.append(id)
                    no_ids += 1
                print(f"generated {no_ids} IDs")
            else:
                ids.append(line)
    print(
        f"Will download {len(ids)} {download_type} records from GenBank release {release_number}"
    )

    # generate make file
    print("Generating pipeline")
    pg = PipelineGenerator(make_file)

    efetch = "/usr/local/edirect-17.0/efetch"
    is_file_empty =  "/usr/local/cavspipes-1.1.0/is_file_empty"

    for id in ids:
        output_sequence_file = f"{output_dir}/{id}.{ext}"
        log = ""
        err = f"{output_sequence_file}.err"
        tgt = f"{output_sequence_file}.OK"
        dep = ""
        cmd = f"{efetch} -db nuccore -id {id} -format {download_type} > {output_sequence_file} 2> {err}"
        cmd += f"; {is_file_empty} {output_sequence_file}"
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

    def print(self):
        print(".DELETE_ON_ERROR:")
        for i in range(len(self.tgts)):
            print(f"{self.tgts[i]} : {self.deps[i]}")
            print(f"\t{self.cmds[i]}")
            print(f"\ttouch {self.tgts[i]}")

    def write(self):
        with open(self.make_file, "w") as f:
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


if __name__ == "__main__":
    main() # type: ignore
