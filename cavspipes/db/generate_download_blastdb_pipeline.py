#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2024 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import gzip
import sys
import click
import re
import subprocess


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    help="make file name",
)
@click.option(
    "-d", "--database", default="prokaryote", show_default=True, help="database to download"
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
    generate_download_blastdb_pipeline -d prokaryote

    e.g.  generate_download_blastdb_pipeline -d prokaryote -m download_blastdb.mk -o /home/atks/downloads

    \b
    database prokaryote|fungi|invertebrate|
             mitochondrion|plant|plasmid|protozoa|plastid|
             viral|vertebrate_mammalian|vertebrate_other
    """
    # check on database validity
    if database not in [
        "prokaryote"
         ]:
        print("error : database not valid\n", file=sys.stderr)
        with click.Context(main) as ctx:
            click.echo(main.get_help(ctx))
            exit()

    if make_file is None:
        make_file = f"download_blastdb_{database}_db.mk"

    if output_directory is None:
        output_directory = os.getcwd()

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("database", database))
    print("\t{0:<20} :   {1:<10}".format("output_directory", output_directory))
    print("\n")
    print(
        "Please invoked pipeline with not more than 8 jobs running concurrently due to NCBI restrictions"
    )
    print("i.e. make -f download_blastdb.mk -j 8 -k")

    output_dir = f"{output_directory}/{database}"
    print(f"\nDatabase will be downloaded to {output_dir}")

    # try:
    #     os.makedirs(output_dir, exist_ok=True)
    # except OSError as error:
    #     print(f"Directory {output_dir} cannot be created")

    # get listing of files to download
    print("Getting directory listings")
    out = subprocess.run(
        ["curl", f"https://ftp.ncbi.nlm.nih.gov/blast/db/"],
        text=True,
        capture_output=True,
    )

    if database == "prokaryote":
        base_name = "nt_prok"

    files = []
    for line in out.stdout.splitlines():
        m = re.search(rf">({base_name}.\d+.tar.gz)<", line)
        if m is not None:
            print(m.group(1))
            files.append(m.group(1))

    print(f"Downloading {len(files)} files.")

    # generate make file
    print("Generating pipeline")
    pg = PipelineGenerator(make_file)

    # download each file
    for file_name in files:
        err = f"{output_dir}/{file_name}.err"
        tgt = f"{output_dir}/{file_name}.OK"
        dep = ""
        cmd = f"wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/{file_name} -P {output_dir} 2> {err}"
        pg.add(tgt, dep, cmd)




    # clean files
    cmd = f"rm {output_dir}/*.fna.gz  {output_dir}/*.OK {output_dir}/*.err"
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
