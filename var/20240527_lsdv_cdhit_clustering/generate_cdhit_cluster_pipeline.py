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
from shutil import copy2

@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_cdhit_clustering.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option("-r", "--ref_fasta_file", required=True, help="sample file")
def main(make_file, working_dir, ref_fasta_file):
    """
    Cluster reference sequences over a set of thresholds for clustering

    e.g. generate_cdhit_cluster_pipeline
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("ref_fasta_file", ref_fasta_file))

    # create directories in destination folder directory
    log_dir = f"{working_dir}/log"
    try:
        os.makedirs(log_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    trace_dir = f"{working_dir}/trace"

    # initialize
    pg = PipelineGenerator(make_file)

    #programs
    cdhit =  "/usr/local/cd-hit-4.8.1/cd-hit"

    cutoffs = [0.80, 0.85, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    for i in range(1,10):
        cutoffs.append(0.99 + i*0.001)
    for i in range(1,10):
        cutoffs.append(0.999 + i*0.0001)
    for i in range(1,10):
        cutoffs.append(0.9999 + i*0.00001)
    for i in range(1,10):
        cutoffs.append(0.99999 + i*0.000001)
    cutoffs.append(1.00)

    # cluster sequences
    for cutoff in cutoffs:
        #cdhit -i nr -o nr100 -c 1.00 -n 5 -M 2000
        clustered_fasta_file = f"{working_dir}/clustered.{cutoff*100}.fasta"
        log_file = f"{log_dir}/{cutoff}.cutoff.log"
        dep = ""
        cmd = f"{cdhit} -i {ref_fasta_file} -o {clustered_fasta_file} -c {cutoff} -T 10 -M 2000 > {log_file}"
        tgt = f"{clustered_fasta_file}.OK"
        pg.add_srun(tgt, dep, cmd, 10)

    # write make file
    print("Writing pipeline")
    pg.write()

    # copy files to trace
    copy2(__file__, trace_dir)
    copy2(make_file, trace_dir)
    copy2(ref_fasta_file, trace_dir)

class PipelineGenerator(object):
    def __init__(self, make_file):
        self.make_file = make_file
        self.tgts = []
        self.deps = []
        self.cmds = []
        self.clean_cmd = ""

    def add_srun(self, tgt, dep, cmd, cpu):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(f"srun --mincpus {cpu} {cmd}")

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
