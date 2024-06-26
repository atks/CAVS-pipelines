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

import os
import click
import sys


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="ont_basecall_deploy_and_qc.mk",
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
    "-q", "--qscore", help="Minimum passing read Q-score", required=False, default=8
)
@click.option(
    "-f", "--flowcell", help="Flow cell ID", required=False, default="FLO-MIN106"
)
@click.option("-k", "--kit", default="SQK-LSK109", show_default=True, help="Kit ID")
@click.option(
    "-b",
    "--barcode_kit",
    default="EXP-NBD104",
    show_default=True,
    help="Barcode Kit ID",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
@click.option(
    "--base_call/--no-base-call",
    default=False,
    help="do or do not perform base calling",
)
def main(
    make_file,
    run_id,
    nanopore_dir,
    working_dir,
    qscore,
    flowcell,
    kit,
    barcode_kit,
    sample_file,
    base_call,
):
    """
    Moves ONT fastq files to a destination and performs QC

    e.g. var_generate_ont_deploy_and_qc_pipeline -r ont23 -m runthis.mk -w . -i nanopore_raw
    """
    log_dir = f"{working_dir}/log"
    dest_dir = working_dir + "/" + run_id
    fastq_path = ""
    if base_call:
        fastq_path = f"{working_dir}/demux"
    else:
        for dirpath, dirnames, filenames in os.walk(nanopore_dir):
            for dirname in dirnames:
                if dirname == "fastq_pass":
                    fastq_path = os.path.join(dirpath, dirname)

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<20} :   {1:<10}".format("nanopore_dir", nanopore_dir))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("qscore", qscore))
    print("\t{0:<20} :   {1:<10}".format("flowcell", flowcell))
    print("\t{0:<20} :   {1:<10}".format("kit", kit))
    print("\t{0:<20} :   {1:<10}".format("barcode_kit", barcode_kit))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("base_call", "True" if base_call else "False"))
    print("\t{0:<20} :   {1:<10}".format("dest_dir", dest_dir))
    print("\t{0:<20} :   {1:<10}".format("fastq_path", fastq_path))

    try:
        os.makedirs(f"{working_dir}/log", exist_ok=True)
    except OSError as error:
        print(f"Directory {working_dir} cannot be created")

    run = Run(run_id)
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, barcode = line.rstrip().split("\t")
                run.add_sample(index, sample_id, barcode)

    # initialize
    pg = PipelineGenerator(make_file)

    if base_call:
        guppy_base_caller = "/usr/local/ont-guppy-6.1.5/bin/guppy_basecaller"
        # base call
        # guppy_basecaller -i {input_dir} -r -s {output_dir} --flowcell FLO-MIN106  --kit SQK-LSK109 --min_qscore 8 -x auto --compress_fastq
        output_dir = f"{working_dir}/fastq/"
        log = f"{log_dir}/guppy_base_caller.log"
        err = f"{log_dir}/guppy_base_caller.err"
        tgt = f"{log_dir}/guppy_base_caller.OK"
        dep = ""
        cmd = f'{guppy_base_caller} -i "{nanopore_dir}" -r -s {output_dir} --flowcell {flowcell}  --kit {kit} --min_qscore {qscore}  -x auto --compress_fastq > {log} 2> {err}'
        pg.add(tgt, dep, cmd)

        # demux
        guppy_barcoder = "/usr/local/ont-guppy-6.1.5/bin/guppy_barcoder"
        # guppy_barcoder -i {input_dir} -r -s {output_dir} --barcode_kits EXP-NBD104 --trim_barcodes -t 12 --compress_fastq
        input_dir = f"{working_dir}/fastq/pass"
        output_dir = f"{working_dir}/demux"
        log = f"{log_dir}/demux.log"
        err = f"{log_dir}/demux.err"
        dep = f"{log_dir}/guppy_base_caller.OK"
        tgt = f"{log_dir}/demux.OK"
        cmd = f'{guppy_barcoder} -i "{input_dir}" -r -s {output_dir} --barcode_kits {barcode_kit} --trim_barcodes --compress_fastq -t 2 > {log} 2> {err}'
        pg.add(tgt, dep, cmd)

    # create directories in destination folder directory
    analysis_dir = f"{dest_dir}/analysis"
    new_dir = ""
    try:
        new_dir = analysis_dir
        os.makedirs(new_dir, exist_ok=True)
        for sample in run.samples:
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
            os.makedirs(new_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    nanoplot_multiqc_dep = ""
    kraken2_multiqc_dep = ""
    kraken2_reports_dep = ""
    kraken2_reports = ""

    for idx, sample in enumerate(run.samples):
        # combine gzip files
        input_dir = f"{fastq_path}/{sample.barcode}"
        output_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}.zcat.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.zcat.err"
        dep = f"{log_dir}/demux.OK" if base_call else ""
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.zcat.OK"
        cmd = f"zcat {input_dir}/*.fastq.gz | gzip -c > {output_file}"
        pg.add(tgt, dep, cmd)

        # nanoplot
        nanoplot = "/usr/local/bin/NanoPlot"
        input_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/nanoplot_result"
        log = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.zcat.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.OK"
        nanoplot_multiqc_dep += f" {tgt}"
        cmd = f"{nanoplot} --fastq {input_file} -o {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # kraken2
        kraken2 = "/usr/local/kraken2-2.1.2/kraken2"
        kraken2_std_db = "/usr/local/ref/kraken2/20220816_standard"
        input_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        report_file = f"{output_dir}/report.txt"
        log = f"{output_dir}/report.log"
        err = f"{output_dir}/run.log"
        dep = (
            f"{log_dir}/{sample.idx}_{sample.id}.zcat.OK"
            if idx == 0
            else f"{log_dir}/{run.samples[idx-1].idx}_{run.samples[idx-1].id}.kraken2.OK"
        )
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        kraken2_multiqc_dep += f" {tgt}"
        kraken2_reports += f" {log}"
        cmd = (
            f"{kraken2} --db {kraken2_std_db} {input_fastq_file} --use-names --report {report_file} "
            + '| perl -lane \'{@F=split("\\t"); $$F[2]=~/(.+) \(taxid (\d+)\)/; print "$$F[0]\\t$$F[1]\\t$$1\\t$$2\\t$$F[3]\\t$$F[4]"}\''
            + f"> {log} 2> {err}"
        )
        pg.add(tgt, dep, cmd)

        # plot kronatools radial tree
        kt_import_taxonomy = "/usr/local/KronaTools-2.8.1/bin/ktImportTaxonomy"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        input_txt_file = f"{output_dir}/report.log"
        output_html_file = f"{output_dir}/krona_radial_tree.html"
        log = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.krona_radial_tree.OK"
        cmd = f"{kt_import_taxonomy} -q 2 -t 4 {input_txt_file} -o {output_html_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

    multiqc = "/usr/local/bin/multiqc"

    # plot nanoplot results
    analysis = "nanoplot"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = nanoplot_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -m nanostat -o {output_dir} -n {analysis} -dd -1 --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kraken2 results
    analysis = "kraken2"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = kraken2_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -m kraken -o {output_dir} -n {analysis} -dd -1 --no-ansi > {log} 2> {err}"
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
    cmd = f"{kt_import_taxonomy} -q 2 -t 4 {input_txt_files} -o {output_html_file} > {log} 2> {err}"
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
    def __init__(self):
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
    main()
