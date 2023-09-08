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
    default="download_genbank_database.mk",
    show_default=True,
    help="make file name",
)
@click.option(
    "-d", "--database", default="vrl", show_default=True, help="database to download"
)
@click.option(
    "-o",
    "--output_directory",
    default=os.getcwd(),
    show_default=True,
    help="output directory, database files will be downloaded to <out_dir>/<db_release>/<db>",
)
def main(make_file, database, output_directory):
    """
    Download genbank database

    e.g.  generate_download_genbank_db_pipeline -d vrl -m download_gb.mk -o /home/atks/downloads

    \b
    database bct|con|env|est|gss|htc|htg
             inv|mam|pat|phg|pln|pri|rod
             sts|syn|tsa|una|vrl|vrt
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("database", database))
    print("\t{0:<20} :   {1:<10}".format("output_directory", output_directory))

    for database in database.split(","):
        if database not in [
            "bct",
            "con",
            "env",
            "est",
            "gss",
            "htc",
            "htg",
            "inv",
            "mam",
            "pat",
            "phg",
            "pln",
            "pri",
            "rod",
            "sts",
            "syn",
            "tsa",
            "una",
            "vrl",
            "vrt",
        ]:
            print("error : database not valid\n", file=sys.stderr)
            with click.Context(main) as ctx:
                click.echo(main.get_help(ctx))
            exit()

    # check on database validity
    # PRI - primate sequences
    # ROD - rodent sequences
    # MAM - other mammalian sequences
    # VRT - other vertebrate sequences
    # INV - invertebrate sequences
    # PLN - plant, fungal, and algal sequences
    # BCT - bacterial sequences
    # VRL - viral sequences
    # PHG - bacteriophage sequences
    # SYN - synthetic sequences
    # UNA - unannotated sequences
    # ENV - environmental sampling sequences
    # TSA - transcriptome shotgun data
    # PAT - patent sequences
    # high throughput
    # EST - EST sequences (expressed sequence tags)
    # STS - STS sequences (sequence tagged sites)
    # GSS - GSS sequences (genome survey sequences)
    # HTG - HTG sequences (high-throughput genomic sequences)
    # HTC - unfinished high-throughput cDNA sequencing

    release_number = subprocess.run(
        ["curl", f"https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number"],
        text=True,
        capture_output=True,
    ).stdout.rstrip()
    output_dir = f"{output_directory}/{release_number}/{database}"
    print(f"\nDatabase will be downloaded to {output_dir}")
    output_dir = f"{output_directory}/{release_number}/{database}"

    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {output_dir} cannot be created")

    # get listing of files to download
    print("Getting directory listings")
    out = subprocess.run(
        ["curl", f"https://ftp.ncbi.nlm.nih.gov/genbank/"],
        text=True,
        capture_output=True,
    )

    files = []
    for line in out.stdout.splitlines():
        m = re.search(f".+(gb{database}\d+\.seq\.gz)<", line)
        if m != None:
            files.append(m.group(1))
        else:
            pass

    # sort files?

    print(f"Downloading {len(files)} files.")

    # generate make file
    print("Generating pipeline")
    pg = PipelineGenerator(make_file)

    concat_fasta_file_list = ""
    concat_fasta_file_OK_list = ""

    for file_name in files:
        file_core = file_name.split(".")[0]

        # download
        output_genbank_file = f"{output_dir}/{file_name}"
        log = f"{output_genbank_file}.log"
        err = f"{output_genbank_file}.err"
        dep = ""
        tgt = f"{output_genbank_file}.OK"
        cmd = f"wget -c https://ftp.ncbi.nlm.nih.gov/genbank/{file_name} -P {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # convert to fasta
        gb2fasta = "/home/atks/programs/cavspipes/gen/gb2fasta"
        input_genbank_file = f"{output_dir}/{file_name}"
        output_fasta_file = f"{output_dir}/{file_core}.fasta.gz"
        log = f"{output_fasta_file}.err"
        err = f"{output_fasta_file}.err"
        dep = f"{input_genbank_file}.OK"
        tgt = f"{output_fasta_file}.OK"
        cmd = f"{gb2fasta} {input_genbank_file} -o {output_fasta_file} > {log} 2> {err}"
        concat_fasta_file_list += f" {output_fasta_file}"
        concat_fasta_file_OK_list += f" {output_fasta_file}.OK"
        pg.add(tgt, dep, cmd)

    # combine into one file
    output_fasta_file = f"{output_dir}/genbank.{release_number}.{database}.fasta.gz"
    err = f"{output_fasta_file}.err"
    tgt = f"{output_fasta_file}.OK"
    dep = concat_fasta_file_OK_list
    cmd = f"zcat {concat_fasta_file_list} | gzip > {output_fasta_file} 2> {err}"
    pg.add(tgt, dep, cmd)

    # get headers
    input_fasta_file = f"{output_dir}/genbank.{release_number}.{database}.fasta.gz"
    output_text_file = f"{output_dir}/genbank.{release_number}.{database}.id.txt"
    log = f"{output_text_file}.log"
    err = f"{output_text_file}.err"
    tgt = f"{output_text_file}.OK"
    dep = f"{input_fasta_file}.OK"
    cmd = f'zcat {input_fasta_file} | grep -P "^>" > {output_text_file} 2> {err}'
    pg.add(tgt, dep, cmd)

    # clean files
    cmd = f"rm -fr {output_dir}/gbvrl*.gz  {output_dir}/*.OK {output_dir}/*.err {output_dir}/*.log"
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
    main()
