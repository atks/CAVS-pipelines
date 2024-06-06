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
    default="download_uniprot_sequences.mk",
    show_default=True,
    help="make file name",
)
@click.option(
    "-o",
    "--output_dir",
    default=os.getcwd(),
    show_default=True,
    help="output directory",
)
def main(make_file, output_dir):
    """
    Download genbank sequences

    e.g.  generate_download_genbank_seq_pipeline -s id.txt -m download_gb_seq.mk -o /home/atks/downloads
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("output dir", output_dir))

    release_notes = subprocess.run(
        [
            "curl",
            f"https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.release_note",
        ],
        text=True,
        capture_output=True,
    ).stdout.rstrip()

    m = re.search(r"Release: ([\d_]+),", release_notes)

    print(release_notes)

    if m is not None:
        release_number = m.group(1)
    print(f"UniRef {release_number}")

    # https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
    # https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
    # https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
    # https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

    # generate make file
    print("Generating pipeline")
    pg = PipelineGenerator(make_file)

    # download
    input_fasta_file = (
        "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta"
    )
    output_uniprot_file = f"uniref90.{release_number}.fasta.gz"
    log = f"{output_uniprot_file}.log"
    err = f"{output_uniprot_file}.err"
    dep = ""
    tgt = f"{output_uniprot_file}.OK"
    cmd = f"wget -c {input_fasta_file} -P {output_dir} -O {output_uniprot_file} > {log} 2> {err}"
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
    main() # type: ignore[arg-type]
