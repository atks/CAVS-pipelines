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
    default="run_alignment_coverage_analysis.mk",
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
@click.option("-r", "--reference_fasta_file", required=True, help="reference FASTA file")
def main(make_file, working_dir, sample_file, reference_fasta_file):
    """
    Generates alignments and coverage statistics.

    e.g. generate_alignment_coverage_analysis_pipeline
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("reference_fasta_file", reference_fasta_file))

    # create directories in destination folder directory
    alignments_dir = f"{working_dir}/alignments"
    log_dir = f"{working_dir}/log"
    trace_dir = f"{working_dir}/trace"
    stats_dir = f"{alignments_dir}/stats"
    try:
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    align_and_consensus = (
        "/home/atks/programs/CAVS-pipelines/var/20240626_asfv_ilm59_coverage/align_and_consensus.py"
    )
    multiqc = "docker run -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc"

    # initialize
    pg = PipelineGenerator(make_file)

    samples = []
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                sample_id, fastq1, fastq2 = line.rstrip().split("\t")
                samples.append(Sample(sample_id, fastq1, fastq2))

    samtools_multiqc_dep = ""

    # align and consensus
    for s in samples:
        log_file = f"{log_dir}/{s.id}.log"
        output_dir = f"{alignments_dir}/{s.id}"
        dep = ""
        cmd = f"{align_and_consensus} -1 {s.fastq1} -2 {s.fastq2} -r {reference_fasta_file} -s {s.id} -w {output_dir} > {log_file}"
        tgt = f"{log_dir}/{s.id}.OK"
        pg.add(tgt, dep, cmd)
        samtools_multiqc_dep += f"{tgt} "

    # plot samtools
    analysis = "samtools"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = samtools_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {alignments_dir}; {multiqc} . -m samtools -f -o {stats_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    pg.add_clean(f"rm -rf {alignments_dir} {log_dir}")

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


class Sample(object):
    def __init__(self, id, fastq1, fastq2):
        self.id = id
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"id        : {self.id}")
        print(f"fastq1    : {self.fastq1}")
        print(f"fastq2    : {self.fastq2}")

if __name__ == "__main__":
    main()  # type: ignore
