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
import sys
import re
import sys
from shutil import copy2
from datetime import datetime


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="ilm_deploy_and_qc.mk",
    help="make file name",
)
@click.option("-r", "--run_id", required=True, help="Run ID")
@click.option("-i", "--illumina_dir", required=True, help="illumina directory")
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(make_file, run_id, illumina_dir, working_dir, sample_file):
    """
    Moves Illumina fastq files to a destination and performs QC

    e.g. var_generate_ilm_deploy_and_qc_pipeline -r ilm23 -i illumina_raw -i ilm23.sa
    """
    dest_dir = working_dir + "/" + run_id
    aux_dir = f"{working_dir}/aux"

    fastq_dir = ""
    illumina_dir = os.path.abspath(illumina_dir)
    for dirpath, dirnames, filenames in os.walk(illumina_dir):
        for dirname in dirnames:
            if dirname == "Fastq":
                fastq_dir = os.path.join(dirpath, dirname)

    sample_sheet_csv_file = ""
    for dirpath, dirnames, filenames in os.walk(illumina_dir):
        for filename in filenames:
            if filename == "SampleSheet.csv":
                sample_sheet_csv_file = os.path.join(dirpath, filename)

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<20} :   {1:<10}".format("illumina_dir", illumina_dir))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("dest_dir", dest_dir))
    print("\t{0:<20} :   {1:<10}".format("fastq_path", fastq_dir))

    # programs
    bcl2fastq = "/usr/local/bin/bcl2fastq"
    fastqc = "/usr/local/FastQC-0.12.1/fastqc"
    kraken2 = "/usr/local/kraken2-2.1.2/kraken2"
    kraken2_std_db = "/usr/local/ref/kraken2/20210908_standard"
    kt_import_taxonomy = "/usr/local/KronaTools-2.8.1/bin/ktImportTaxonomy"
    multiqc = "/usr/local/bin/multiqc"

    # read sample file
    run = Run(run_id)
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, fastq1, fastq2 = line.rstrip().split("\t")
                run.add_sample(index, sample_id, fastq1, fastq2)

    # create directories in destination folder directory
    log_dir = f"{working_dir}/log"
    untrimmed_fastq_dir = f"{working_dir}/untrimmed_fastq"
    untrimmed_dest_dir = f"{dest_dir}/untrimmed"
    analysis_dir = f"{dest_dir}/analysis"
    new_dir = ""
    try:
        new_dir = log_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = untrimmed_fastq_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = untrimmed_dest_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = aux_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = analysis_dir
        os.makedirs(new_dir, exist_ok=True)
        for sample in run.samples:
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/untrimmed_fastqc_result"
            os.makedirs(new_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)
    multiqc_dep = ""

    # create sample sheet
    make_sample_sheet = "make_ilm_sample_sheet_without_adaptor_trimming.py"
    input_csv_file = sample_sheet_csv_file
    output_csv_file = f"{untrimmed_fastq_dir}/no_trim_sample_sheet.csv"
    tgt = f"{output_csv_file}.OK"
    dep = ""
    cmd = f"{make_sample_sheet} {input_csv_file} -o {output_csv_file}"
    pg.add(tgt, dep, cmd)

    # generate untrimmed fastq
    # bcl2fastq -o untrimmed_fastq -R raw/20230309_FS10001208_21_BSB09417-2507 --sample-sheet notrim_SampleSheet.csv
    input_csv_file = f"{untrimmed_fastq_dir}/no_trim_sample_sheet.csv"
    illumina_project_dir = os.path.dirname(sample_sheet_csv_file)
    log = f"{untrimmed_fastq_dir}/trimmed.log"
    tgt = f"{untrimmed_fastq_dir}/trimmed.OK"
    dep = f"{input_csv_file}.OK"
    cmd = f"{bcl2fastq} -o {untrimmed_fastq_dir} -R {illumina_project_dir} --interop-dir {untrimmed_fastq_dir} --sample-sheet {input_csv_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    # analyze
    fastqc_multiqc_dep = ""
    fastqc_directories = ""

    # analyze untrimmed reads
    for idx, sample in enumerate(run.samples):
        # copy the files
        src_fastq1 = f"{untrimmed_fastq_dir}/{sample.fastq1}"
        dst_fastq1 = (
            f"{untrimmed_dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz"
        )
        tgt = f"{log_dir}/untrimmed_{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        dep = f"{untrimmed_fastq_dir}/trimmed.OK"
        cmd = f"cp {src_fastq1} {dst_fastq1}"
        pg.add(tgt, dep, cmd)

        src_fastq2 = f"{untrimmed_fastq_dir}/{sample.fastq2}"
        dst_fastq2 = (
            f"{untrimmed_dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz"
        )
        tgt = f"{log_dir}/untrimmed_{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        dep = f"{untrimmed_fastq_dir}/trimmed.OK"
        cmd = f"cp {src_fastq2} {dst_fastq2}"
        pg.add(tgt, dep, cmd)

        src_fastq1 = dst_fastq1
        src_fastq2 = dst_fastq2

        # fastqc
        fastqc = "/usr/local/FastQC-0.11.9/fastqc"
        fastqc_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/untrimmed_fastqc_result"
        fastqc_directories += f"{fastqc_dir}\n"

        log = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_fastqc1.log"
        err = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_fastqc1.err"
        tgt = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_fastqc1.OK"
        dep = f"{log_dir}/untrimmed_{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"{fastqc} {src_fastq1} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        log = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_fastqc2.log"
        err = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_fastqc2.err"
        tgt = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_fastqc2.OK"
        dep = f"{log_dir}/untrimmed_{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"{fastqc} {src_fastq2} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

    with open(f"{aux_dir}/untrimmed_fastqc_dir.txt", "w") as file:
        file.write(fastqc_directories)

    # plot fastqc multiqc results
    analysis = "untrimmed_fastqc"
    dir_list_file = f"{aux_dir}/untrimmed_fastqc_dir.txt"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} -l {dir_list_file} -f -m fastqc -o {output_dir} -n {analysis} -dd -1 --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    fastqc_multiqc_dep = ""
    fastqc_directories = ""
    kraken2_multiqc_dep = ""
    kraken2_reports = ""

    # analyze trimmed reads
    for idx, sample in enumerate(run.samples):
        # copy the files
        src_fastq1 = f"{fastq_dir}/{sample.fastq1}"
        dst_fastq1 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz"
        tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        dep = ""
        cmd = f"cp {src_fastq1} {dst_fastq1}"
        pg.add(tgt, dep, cmd)

        src_fastq2 = f"{fastq_dir}/{sample.fastq2}"
        dst_fastq2 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz"
        tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        dep = ""
        cmd = f"cp {src_fastq2} {dst_fastq2}"
        pg.add(tgt, dep, cmd)

        src_fastq1 = dst_fastq1
        src_fastq2 = dst_fastq2

        # fastqc
        fastqc_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"
        fastqc_directories += f"{fastqc_dir}\n"

        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.err"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"{fastqc} {src_fastq1} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.err"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"{fastqc} {src_fastq2} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        # kraken2
        input_fastq_file1 = f"{src_fastq1}"
        input_fastq_file2 = f"{src_fastq2}"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        report_file = f"{output_dir}/report.txt"
        log = f"{output_dir}/report.log"
        err = f"{output_dir}/run.log"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK {log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        if idx == 0:
            dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK {log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        else:
            dep = f"{log_dir}/{run.idx}_{run.samples[idx-1].idx}_{run.samples[idx-1].id}_R1.fastq.gz.OK {log_dir}/{run.idx}_{run.samples[idx-1].idx}_{run.samples[idx-1].id}_R2.fastq.gz.OK {log_dir}/{run.samples[idx-1].idx}_{run.samples[idx-1].id}.kraken2.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        kraken2_multiqc_dep += f" {tgt}"
        kraken2_reports += f" {report_file}"
        cmd = f"{kraken2} --db {kraken2_std_db} --paired {input_fastq_file1} {input_fastq_file2} --use-names --report {report_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # plot kronatools radial tree
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        input_txt_file = f"{output_dir}/report.txt"
        output_html_file = f"{output_dir}/krona_radial_tree.html"
        log = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.krona_radial_tree.OK"
        cmd = f"{kt_import_taxonomy} -m 3 -t 5 {input_txt_file} -o {output_html_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

    with open(f"{aux_dir}/trimmed_fastqc_dir.txt", "w") as file:
        file.write(fastqc_directories)

    # plot fastqc multiqc results
    analysis = "fastqc"
    dir_list_file = f"{aux_dir}/trimmed_fastqc_dir.txt"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} -l {dir_list_file} -f -m fastqc -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kraken2 multiqc results
    analysis = "kraken2"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = kraken2_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -f -m kraken -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kronatools radial tree
    kt_import_taxonomy = "/usr/local/KronaTools-2.8.1/bin/ktImportTaxonomy"
    analysis = "kraken2"
    output_dir = f"{analysis_dir}/all/{analysis}"
    input_txt_files = kraken2_reports
    output_html_file = f"{output_dir}/krona_radial_tree.html"
    log = f"{log_dir}/{analysis}.krona_radial_tree.log"
    err = f"{log_dir}/{analysis}.krona_radial_tree.err"
    dep = f"{kraken2_multiqc_dep}"
    tgt = f"{log_dir}/{analysis}.krona_radial_tree.OK"
    cmd = f"{kt_import_taxonomy} -m 3 -t 5 {input_txt_files} -o {output_html_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    pg.add_clean(f"rm -fr {log_dir} {untrimmed_fastq_dir} {dest_dir}")

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
                f.write(f"\t{self.cmds[i]}\n")
                f.write(f"\ttouch {self.tgts[i]}\n\n")

            if self.clean_cmd != "":
                f.write(f"clean : \n")
                f.write(f"\t{self.clean_cmd}\n")


class Sample(object):
    def __init__(self):
        self.idx = ""
        self.id = ""
        self.fastq1 = ""
        self.fastq2 = ""

    def __init__(self, idx, id, fastq1, fastq2):
        self.idx = idx
        self.id = id
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"index   : {self.idx}")
        print(f"id      : {self.id}")
        print(f"barcode : {self.fastq1}")
        print(f"barcode : {self.fastq2}")


class Run(object):
    def __init__(self, id):
        self.idx = id[4::]
        self.id = id
        self.samples = []

    def add_sample(self, idx, sample_id, fastq1, fastq2):
        self.samples.append(Sample(idx, sample_id, fastq1, fastq2))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id}")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")


if __name__ == "__main__":
    main()
