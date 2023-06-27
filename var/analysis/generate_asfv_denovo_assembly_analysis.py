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

import os
import click


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_asfv_bulan_analysis.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
def main(make_file, working_dir):
    """
    African Swine Fever Virus Denovi Assembly evaluation

    e.g. generate_asfv_denove_assembly_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    # initialize
    pg = PipelineGenerator(make_file)

    # programs
    efetch = "/usr/local/edirect-17.0/efetch"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    fastq_dir = f"{working_dir}/fastq"
    fasta_dir = f"{working_dir}/fasta"
    bam_dir = f"{working_dir}/bam"
    consensus_dir = f"{working_dir}/consensus"
    stats_dir = f"{working_dir}/stats"
    coverage_stats_dir = f"{working_dir}/stats/coverage"

    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(consensus_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(coverage_stats_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    #################
    # reference files
    #################
    id = "FR682468.2"
    output_fasta_file = f"{ref_dir}/{id}.fasta"
    dep = ""
    tgt = f"{output_fasta_file }.OK"
    cmd = f"{efetch} -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    ###############
    # sequence file
    ###############
    wb1_dir = "/home/atks/analysis/20230522_asfv_ilm51/M230255_wildboar"
    bulan_pig46_2_dir = "/home/atks/analysis//20230522_asfv_ilm51/M230446_pig2"
    bulan_pig48_4_dir = "/home/atks/analysis//20230522_asfv_ilm51/M230448_pig4"

    # extract paired sequences from original fastq file based on a bam file

    # denovo assembly

    # clean
    pg.add_clean(f"rm -fr {ref_dir} {fastq_dir} {fasta_dir} {bam_dir} {stats_dir} ")

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
                if "|" in self.cmds[i]:
                    f.write(f"\tset -o pipefail; {self.cmds[i]}\n")
                else:
                    f.write(f"\t{self.cmds[i]}\n")
                f.write(f"\ttouch {self.tgts[i]}\n\n")

            if self.clean_cmd != "":
                f.write(f"clean : \n")
                f.write(f"\t{self.clean_cmd}\n")


class Sample(object):
    def __init__(self):
        self.id = ""
        self.name = ""
        self.seq_tech = ""
        self.fastq1 = ""
        self.fastq2 = ""
        self.fastq1_OK = ""
        self.fastq2_OK = ""

    def __init__(self, id, name, seq_tech, fastq1, fastq2):
        self.id = id
        self.name = name
        self.seq_tech = seq_tech
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.fastq1_OK = ""
        self.fastq2_OK = ""

    def __init__(self, id, name, seq_tech, fastq1, fastq2, fastq1_OK, fastq2_OK):
        self.id = id
        self.name = name
        self.seq_tech = seq_tech
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.fastq1_OK = fastq1_OK
        self.fastq2_OK = fastq2_OK

    def print(self):
        print(f"id         : {self.id}")
        print(f"name       : {self.name}")
        print(f"seq tech   : {self.seq_tech}")
        print(f"fastq1     : {self.fastq1}")
        print(f"fastq2     : {self.fastq2}")
        print(f"fastq1_OK  : {self.fastq1_OK}")
        print(f"fastq2_OK  : {self.fastq2_OK}")


if __name__ == "__main__":
    main()
