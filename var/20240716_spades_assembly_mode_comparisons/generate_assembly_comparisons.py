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
    default="run_assembly_comparisons_analysis.mk",
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
    Generates assembly comparisons statistics.

    e.g. generate_assembly_comparisons_pipeline.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # create directories in destination folder directory
    pairwise_alignments_dir = f"{working_dir}/pairwise_alignments"
    plot_dir = f"{working_dir}/plot"
    trace_dir = f"{working_dir}/trace"
    try:
        os.makedirs(pairwise_alignments_dir, exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    nucmer = "/usr/local/mummer-4.0.0/bin/nucmer"
    dnadiff = "/usr/local/mummer-4.0.0/bin/dnadiff"
    aggregate_nucmer_reports = "/home/atks/programs/CAVS-pipelines/var/20240716_spades_assembly_mode_comparisons/aggregate_nucmer_reports.py"

    # initialize
    pg = PipelineGenerator(make_file)

    samples = []
    idx = 0
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                idx += 1
                sample_id, contigs_file = line.rstrip().split("\t")
                samples.append(Sample(idx, sample_id, contigs_file))

    # pairwise alignment
    isolate_assembly_dir = f"{working_dir}/isolate"
    metaviral_assembly_dir = f"{working_dir}/metaviral"
    meta_assembly_dir = f"{working_dir}/meta"

    #####################################
    #isolate vs meta pairwise alignment
    #####################################
    all_nucmer_report_files_ok = ""
    input_nucmer_report_files = ""

    for s in samples:

        #nucmer pairwise alignment
        ref_fasta_file = f"{isolate_assembly_dir}/{s.contigs_file}"
        query_fasta_file = f"{meta_assembly_dir}/{s.contigs_file}"
        file_prefix = f"{pairwise_alignments_dir}/isolate_meta_{s.idx}"
        log_file = f"{pairwise_alignments_dir}/isolate_meta_{s.idx}.nucmer.log"
        dep = ""
        cmd = f"{nucmer} {ref_fasta_file} {query_fasta_file} -p {file_prefix} > {log_file}"
        tgt = f"{file_prefix}.delta.OK"
        pg.add(tgt, dep, cmd)

        #isolate vs meta dnadiff analysis
        input_delta_file = f"{pairwise_alignments_dir}/isolate_meta_{s.idx}.delta"
        log_file = f"{pairwise_alignments_dir}/isolate_meta_{s.idx}.dnadiff.log"
        dep = f"{file_prefix}.delta.OK"
        cmd = f"{dnadiff} -d {input_delta_file} -p {file_prefix} 2> {log_file}"
        tgt = f"{file_prefix}.report.OK"
        pg.add(tgt, dep, cmd)

        input_nucmer_report_files += f" {file_prefix}.report"
        all_nucmer_report_files_ok += f" {tgt}"

    # aggregrate reports
    tag = "isolate_meta"
    output_txt_file = f"{plot_dir}/{tag}.txt"
    dep = all_nucmer_report_files_ok
    cmd = f"{aggregate_nucmer_reports} {input_nucmer_report_files} -o {output_txt_file}"
    tgt = f"{output_txt_file}.OK"
    pg.add(tgt, dep, cmd)

    #####################################
    #metaviral vs meta pairwise alignment
    #####################################
    all_nucmer_report_files_ok = ""
    input_nucmer_report_files = ""

    for s in samples:

        #nucmer pairwise alignment
        ref_fasta_file = f"{metaviral_assembly_dir}/{s.contigs_file}"
        query_fasta_file = f"{meta_assembly_dir}/{s.contigs_file}"
        file_prefix = f"{pairwise_alignments_dir}/metaviral_meta_{s.idx}"
        log_file = f"{pairwise_alignments_dir}/metaviral_meta_{s.idx}.log"
        dep = ""
        cmd = f"{nucmer} {ref_fasta_file} {query_fasta_file} -p {file_prefix} > {log_file}"
        tgt = f"{file_prefix}.delta.OK"
        pg.add(tgt, dep, cmd)

        #metaviral vs meta dnadiff analysis
        input_delta_file = f"{pairwise_alignments_dir}/metaviral_meta_{s.idx}.delta"
        log_file = f"{pairwise_alignments_dir}/metaviral_meta_{s.idx}.dnadiff.log"
        dep = f"{file_prefix}.delta.OK"
        cmd = f"{dnadiff} -d {input_delta_file} -p {file_prefix} 2> {log_file}"
        tgt = f"{file_prefix}.report.OK"
        pg.add(tgt, dep, cmd)

        input_nucmer_report_files += f" {file_prefix}.report"
        all_nucmer_report_files_ok += f" {tgt}"

    # aggregrate reports
    tag = "metaviral_meta"
    output_txt_file = f"{plot_dir}/{tag}.txt"
    dep = all_nucmer_report_files_ok
    cmd = f"{aggregate_nucmer_reports} {input_nucmer_report_files} -o {output_txt_file}"
    tgt = f"{output_txt_file}.OK"
    pg.add(tgt, dep, cmd)

    #clean up
    pg.add_clean(f"rm -rf {pairwise_alignments_dir} {plot_dir}")

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
    def __init__(self, idx, id, contigs_file):
        self.idx = idx
        self.id = id
        self.contigs_file = contigs_file

    def print(self):
        print(f"idx           : {self.idx}")
        print(f"id            : {self.id}")
        print(f"contigs_file  : {self.contigs_file}")

if __name__ == "__main__":
    main()  # type: ignore
