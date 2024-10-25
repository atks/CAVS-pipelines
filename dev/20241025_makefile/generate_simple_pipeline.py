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
import click
import sys
from shutil import copy2
from datetime import datetime


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="simple_stuff.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default= f"{os.getcwd()}/simple_stuff",
    show_default=True,
    help="working directory",
)
def main(make_file, working_dir):
    """
    Simple make file example to generate 10 files, concatenate them and clean up temporary files

    e.g. generate_simple_pipeline.py
    """

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working dir", working_dir))

    #version
    version = "1.0.0"

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    trace_dir = f"{working_dir}/trace"
    try:
        os.makedirs(working_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    input_files = ""
    input_files_ok = ""

    ######################
    #1. Generate 100 files
    ######################
    for i in range(1, 100):
        input_file = f" {working_dir}/{i}.txt"
        tgt = f"{working_dir}/{i}.OK"
        dep = ""
        cmd = f"echo {i} > {input_file}"
        pg.add(tgt, dep, cmd)
        input_files += f" {input_file}"
        input_files_ok += f" {working_dir}/{i}.OK"

    #########################
    #2. Concatenate 100 files
    #########################
    output_file = f"{working_dir}/all.txt"
    tgt = f"{working_dir}/concatenated.OK"
    dep = input_files_ok
    cmd = f"cat {input_files} > {output_file}"
    pg.add(tgt, dep, cmd)

    ###########################
    #3. Cleanup temporary files
    ###########################
    tgt = f"{working_dir}/cleaned.OK"
    dep = f"{working_dir}/concatenated.OK"
    cmd = f"rm {input_files}"
    pg.add(tgt, dep, cmd)

    # write make file
    print("Writing pipeline")
    pg.write()

    #copy files to trace
    copy2(__file__, trace_dir)
    copy2(make_file, trace_dir)

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

if __name__ == "__main__":
    main() # type: ignore
