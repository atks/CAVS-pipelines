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
    default="run_consensus_analysis.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(make_file, working_dir, sample_file):
    """
    Consensus and coverage calculations

    e.g. generate_coverage_consensus_analysis_pipeline
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # create directories in destination folder directory
    log_dir = f"{working_dir}/log"
    trace_dir = f"{working_dir}/trace"

    try:
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    align_and_consensus = (
        "/home/atks/programs/CAVS-pipelines/minipipes/align_and_consensus.py"
    )

    # initialize
    pg = PipelineGenerator(make_file)

    sequence = []
    # read from reference sample file
    ##acc-id	country	collection_year	submission_year	fasta_header
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                acc_id, country, collection_year, submission_year, fasta_header = (
                    line.rstrip().split("\t")
                )

                fasta_file = f"/net/singapura/var/projects/lsdv/lsdv_ref/{acc_id}.fasta"
                sequence.append(
                    Sequence(
                        acc_id,
                        country,
                        collection_year,
                        submission_year,
                        fasta_header,
                        fasta_file,
                    )
                )

    # combine data
    ##acc-id	country	collection_year	submission_year	fasta_header
    ilm_fastq1 = (
        "/net/singapura/var/hts/ilm53/53_1_M220323_LSDV_skinswab2_R1.fastq.gz,"
        "/net/singapura/var/hts/ilm53/53_2_M220323_LSDV_skinswab3_R1.fastq.gz,"
        "/net/singapura/var/hts/ilm53/53_3_unclassified_R1.fastq.gz"
    )

    ilm_fastq2 = (
        "/net/singapura/var/hts/ilm53/53_1_M220323_LSDV_skinswab2_R2.fastq.gz,"
        "/net/singapura/var/hts/ilm53/53_2_M220323_LSDV_skinswab3_R2.fastq.gz,"
        "/net/singapura/var/hts/ilm53/53_3_unclassified_R2.fastq.gz"
    )

    ont_fastq = (
        "/net/singapura/var/hts/ont22/22_1_M220338_LSDV_nasal2.fastq.gz,"
        "/net/singapura/var/hts/ont22/22_2_M220338_LSDV_nasal7.fastq.gz,"
        "/net/singapura/var/hts/ont22/22_3_M220338_LSDV_skin1.fastq.gz,"
        "/net/singapura/var/hts/ont22/22_4_M220338_LSDV_skin2.fastq.gz,"
        "/net/singapura/var/hts/ont22/22_5_negative_control.fastq.gz,"
        "/net/singapura/var/hts/ont22/22_6_unclassified.fastq.gz"
    )

    # align and consensus
    for s in sequence:
        log_file = f"{log_dir}/{s.acc_id}.log"
        output_dir = f"{working_dir}/{s.acc_id}"
        dep = ""
        cmd = f"{align_and_consensus} -1 {ilm_fastq1} -2 {ilm_fastq2} -n {ont_fastq} -r {s.fasta_file} -s {s.acc_id} -w {output_dir} > {log_file}"
        tgt = f"{log_dir}/{s.acc_id}.OK"
        pg.add(tgt, dep, cmd)

    # write make file
    print("Writing pipeline")
    pg.write()

    # copy files to trace
    copy2(__file__, trace_dir)
    copy2(make_file, trace_dir)
    copy2(sample_file, trace_dir)


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


class Sequence(object):
    def __init__(
        self,
        acc_id,
        country,
        collection_year,
        submission_year,
        fasta_header,
        fasta_file,
    ):
        self.acc_id = acc_id
        self.country = country
        self.collection_year = collection_year
        self.submission_year = submission_year
        self.fasta_header = fasta_header
        self.fasta_file = fasta_file

    def print(self):
        print(f"acc id            : {self.acc_id}")
        print(f"country           : {self.country}")
        print(f"collection year   : {self.collection_year}")
        print(f"submission year   : {self.submission_year}")
        print(f"fasta header      : {self.fasta_header}")
        print(f"fasta file        : {self.fasta_file}")


if __name__ == "__main__":
    main()  # type: ignore
