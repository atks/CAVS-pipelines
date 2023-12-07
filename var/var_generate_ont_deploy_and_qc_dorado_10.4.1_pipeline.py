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


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="ont_deploy_and_qc.mk",
    help="make file name",
)
@click.option("-r", "--run_id", required=True, help="Run ID")
@click.option("-i", "--nanopore_dir", required=True, help="nanopore directory")
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-q", "--qscore", help="Minimum passing read Q-score", required=False, default=7
)
@click.option(
    "-l", "--len", help="Minimum passing read length", required=False, default=20
)
@click.option("-x", "--memory", help="Memory for fastqc", required=False, default=2048)
@click.option("-k", "--kit", default="SQK-NBD114.24", show_default=True, help="Kit ID")
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(
    make_file,
    run_id,
    nanopore_dir,
    working_dir,
    qscore,
    len,
    memory,
    kit,
    sample_file,
):
    """
    Moves ONT fastq files to a destination and performs QC

    e.g. var_generate_ont_deploy_and_qc_pipeline  -r ont26 -i raw -s ont26.sa
         var_generate_ont_deploy_and_qc_pipeline  -r ont26 -i raw -s ont26.sa --no-base-call
         var_generate_ont_deploy_and_qc_pipeline.py  -r ont30 -i raw -s ont30.sa -b "EXP-NBD104 EXP-NBD114"
    """
    log_dir = f"{working_dir}/log"
    dest_dir = working_dir + "/" + run_id
    fastq_path = f"{working_dir}/demux"

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<20} :   {1:<10}".format("nanopore_dir", nanopore_dir))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("minumum qscore", qscore))
    print("\t{0:<20} :   {1:<10}".format("minumum length", len))
    print("\t{0:<20} :   {1:<10}".format("memory", memory))
    print("\t{0:<20} :   {1:<10}".format("kit", kit))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("dest_dir", dest_dir))
    print("\t{0:<20} :   {1:<10}".format("fastq_path", fastq_path))

    # version
    version = "20231004_ont_deploy_and_qc_dorado_10.4.1_pipeline"

    barcode_kit = ""
    dorado_basecall_model = ""
    if kit == "SQK-NBD114.24":
        barcode_kit = "EXP-NBD104 EXP-NBD114"
        dorado_basecall_model = (
            "/usr/local/dorado-0.5.0/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
        )
    # programs
    dorado_basecaller = "/usr/local/dorado-0.5.0/bin/dorado basecaller"
    guppy_barcoder = "/usr/local/ont-guppy-6.5.7/bin/guppy_barcoder"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    fastqc = f"/usr/local/FastQC-0.12.1/fastqc --adapters /usr/local/FastQC-0.12.1/Configuration/adapter_list.nanopore.txt --memory {memory}"
    multiqc = "/usr/local/bin/multiqc"
    nanoplot = "/usr/local/bin/NanoPlot"
    kraken2 = "/usr/local/kraken2-2.1.2/kraken2"
    kraken2_std_db = "/usr/local/ref/kraken2/20220816_standard"
    kt_import_taxonomy = "/usr/local/KronaTools-2.8.1/bin/ktImportTaxonomy"
    ft = "/usr/local/cavstools-0.0.1/ft"

    run = Run(run_id)
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, barcode = line.rstrip().split("\t")
                run.add_sample(index, sample_id, barcode)

    # create directories in destination folder directory
    analysis_dir = f"{dest_dir}/analysis"
    aux_dir = f"{working_dir}/aux"
    log_dir = f"{working_dir}/log"
    bam_dir = f"{working_dir}/bam"
    fastq_dir = f"{working_dir}/fastq"
    untrimmed_dir = f"{dest_dir}/raw/untrimmed_fastq"
    trimmed_dir = f"{dest_dir}/raw/trimmed_fastq"
    new_dir = ""
    try:
        new_dir = analysis_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = log_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = aux_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = bam_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = fastq_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = trimmed_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = untrimmed_dir
        os.makedirs(new_dir, exist_ok=True)
        for sample in run.samples:
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/untrimmed_fastqc_result"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
            os.makedirs(new_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # base call
    # dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.1.0 pod5s/ > calls.bam
    output_bam_file = f"{working_dir}/bam/basecalls.bam"
    pod5_files_dir = nanopore_dir
    log = f"{log_dir}/dorado_base_caller.log"
    err = f"{log_dir}/dorado_base_caller.err"
    tgt = f"{log_dir}/dorado_base_caller.OK"
    dep = ""
    cmd = f"{dorado_basecaller} -r {dorado_basecall_model} {pod5_files_dir} > {output_bam_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    # convert bam to fastq
    # dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.1.0 pod5s/ > calls.bam
    input_bam_file = f"{working_dir}/bam/basecalls.bam"
    output_fastq_file = f"{working_dir}/fastq/basecalls.fastq"
    tgt = f"{log_dir}/dorado_basecalled_fastq.OK"
    dep = f"{log_dir}/dorado_base_caller.OK"
    cmd = f"{samtools} bam2fq {input_bam_file} > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    # demux
    # guppy_barcoder -i {input_dir} -r -s {output_dir} --barcode_kits EXP-NBD104 -t 12 --compress_fastq
    input_dir = f"{working_dir}/fastq"
    output_dir = f"{working_dir}/demux/untrimmed"
    log = f"{log_dir}/demux.untrimmed.log"
    err = f"{log_dir}/demux.untrimmed.err"
    dep = f"{log_dir}/dorado_basecalled_fastq.OK"
    tgt = f"{log_dir}/demux.untrimmed.OK"
    cmd = f'{guppy_barcoder} -i "{input_dir}" -r -s {output_dir} --barcode_kits "{barcode_kit}" --compress_fastq -t 12 > {log} 2> {err}'
    pg.add(tgt, dep, cmd)

    # demux
    # guppy_barcoder -i {input_dir} -r -s {output_dir} --barcode_kits EXP-NBD104 -t 12 --compress_fastq --enable_trim_barcodes
    input_dir = f"{working_dir}/fastq"
    output_dir = f"{working_dir}/demux/trimmed"
    log = f"{log_dir}/demux.trimmed.log"
    err = f"{log_dir}/demux.trimmed.err"
    dep = f"{log_dir}/dorado_basecalled_fastq.OK"
    tgt = f"{log_dir}/demux.trimmed.OK"
    cmd = f'{guppy_barcoder} -i "{input_dir}" -r -s {output_dir} --barcode_kits "{barcode_kit}" --compress_fastq -t 12 --enable_trim_barcodes > {log} 2> {err}'
    pg.add(tgt, dep, cmd)

    untrimmed_fastqc_multiqc_dep = ""
    untrimmed_fastqc_directories = ""
    trimmed_filtered_fastqc_multiqc_dep = ""
    trimmed_filtered_fastqc_directories = ""

    kraken2_multiqc_dep = ""
    kraken2_reports = ""

    for idx, sample in enumerate(run.samples):
        # combine untrimmed gzip files
        input_dir = f"{working_dir}/demux/untrimmed/{sample.barcode}"
        output_file = f"{untrimmed_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.err"
        dep = f"{log_dir}/demux.untrimmed.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.OK"
        cmd = f"zcat {input_dir}/*.fastq.gz | gzip -c > {output_file}"
        pg.add(tgt, dep, cmd)

        # combine trimmed fastq gzip files
        input_dir = f"{working_dir}/demux/trimmed/{sample.barcode}"
        output_fastq_file = f"{trimmed_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.err"
        dep = f"{log_dir}/demux.trimmed.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.OK"
        cmd = f"zcat {input_dir}/*.fastq.gz | gzip -c > {output_fastq_file}"
        pg.add(tgt, dep, cmd)

        # filter trimmed fastq gzip files
        input_fastq_file = f"{trimmed_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.OK"
        cmd = f"{ft} filter -q {qscore} -l {len} {input_fastq_file}  -o {output_fastq_file}"
        pg.add(tgt, dep, cmd)

        # untrimmed fastqc
        input_file = f"{untrimmed_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/untrimmed_fastqc_result"
        untrimmed_fastqc_directories += f"{output_dir}\n"
        log = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.fastqc.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.fastqc.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.fastqc.OK"
        untrimmed_fastqc_multiqc_dep += f" {tgt}"
        cmd = f"{fastqc} {input_file} -o {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # trimmed filtered fastqc
        input_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"
        trimmed_filtered_fastqc_directories += f"{output_dir}\n"
        log = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.fastqc.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.fastqc.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.fastqc.OK"
        trimmed_filtered_fastqc_multiqc_dep += f" {tgt}"
        cmd = f"{fastqc} {input_file} -o {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # untrimmed nanoplot
        input_file = f"{untrimmed_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = (
            f"{analysis_dir}/{sample.idx}_{sample.id}/untrimmed_nanoplot_result"
        )
        log = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed_nanoplot.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed_nanoplot.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed.zcat.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.untrimmed_nanoplot.OK"
        cmd = f"{nanoplot} --fastq {input_file} -o {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # trimmed filtered nanoplot
        input_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/nanoplot_result"
        log = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.nanoplot.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.nanoplot.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.filtered.nanoplot.OK"
        cmd = f"{nanoplot} --fastq {input_file} -o {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # kraken2
        input_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        report_file = f"{output_dir}/report.txt"
        log = f"{output_dir}/report.log"
        err = f"{output_dir}/run.log"
        dep = (
            f"{log_dir}/{sample.idx}_{sample.id}.trimmed.zcat.OK"
            if idx == 0
            else f"{log_dir}/{run.samples[idx-1].idx}_{run.samples[idx-1].id}.kraken2.OK"
        )
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        kraken2_multiqc_dep += f" {tgt}"
        kraken2_reports += f" {report_file}"
        cmd = f"set -o pipefail; {kraken2} --db {kraken2_std_db} {input_fastq_file} --use-names --report {report_file} > {log} 2> {err}"
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

    with open(f"{aux_dir}/untrimmed_fastqc_dir.txt", "w") as file:
        file.write(untrimmed_fastqc_directories)

    # plot untrimmed fastqc multiqc results
    analysis = "untrimmed_fastqc"
    dir_list_file = f"{aux_dir}/untrimmed_fastqc_dir.txt"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = trimmed_filtered_fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} -l {dir_list_file} -f -m fastqc -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    with open(f"{aux_dir}/filtered_trimmed_fastqc_dir.txt", "w") as file:
        file.write(trimmed_filtered_fastqc_directories)

    # plot filtered trimmed fastqc multiqc results
    analysis = "fastqc"
    dir_list_file = f"{aux_dir}/filtered_trimmed_fastqc_dir.txt"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = trimmed_filtered_fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} -l {dir_list_file} -f -m fastqc -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kraken2 results
    analysis = "kraken2"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = kraken2_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -f -m kraken -o {output_dir} -n {analysis} -dd -1 --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kronatools radial tree
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
    def __init__(self):  # type: ignore
        self.idx = ""
        self.id = ""
        self.barcode = ""

    def __init__(self, idx, id, barcode):
        self.idx = idx
        self.id = id
        self.barcode = barcode

    def print(self):
        print(f"index   : {self.idx}")
        print(f"id      : {self.id}")
        print(f"barcode : {self.barcode}")


class Run(object):
    def __init__(self, id):
        self.idx = id[3::]
        self.id = id
        self.samples = []

    def add_sample(self, idx, sample_id, barcode):
        self.samples.append(Sample(idx, sample_id, barcode))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id}")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")


if __name__ == "__main__":
    main()  # type: ignore
