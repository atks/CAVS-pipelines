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
import re
import sys
from shutil import copy2
from datetime import datetime

@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="novogene_ilm_deploy_and_qc.mk",
    help="make file name",
)
@click.option("-r", "--run_id", required=True, help="Run ID")
@click.option("-i", "--novogene_illumina_dir", required=True, help="Novogene Illumina directory")
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(make_file, run_id, novogene_illumina_dir, working_dir, sample_file):
    """
    Moves Novogene Illumina fastq files to a destination and performs QC

    e.g. vm_generate_novogene_ilm_deploy_and_qc_pipeline -r ilm1 -i raw -i ilm23.sa
    """
    log_dir = f"{working_dir}/log"
    dest_dir = working_dir + "/" + run_id
    illumina_dir = os.path.abspath(novogene_illumina_dir)
    fastq_dir = illumina_dir

    print("\t{0:<21} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<21} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<21} :   {1:<10}".format("novogene_illumina_dir", novogene_illumina_dir))
    print("\t{0:<21} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<21} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<21} :   {1:<10}".format("dest_dir", dest_dir))
    print("\t{0:<21} :   {1:<10}".format("fastq_path", fastq_dir))

    # read sample file
    ## novogene-sample-id fastq-infix
    ## Aq1_21	aquatic_1-21
    ## Aq46_21	aquatic_46-21
    ## Aq47_21	aquatic_47-21
    run = Run(run_id)
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                novogene_sample_id, sample_id = line.rstrip().split("\t")
                #search directory for fastq files
                fastq_dir = f"{novogene_illumina_dir}/01.RawData"
                sample_dir = f"{os.path.abspath(fastq_dir)}/{novogene_sample_id}"
                files = {}
                no_files = 0
                novogene_fastq1s = ""
                novogene_fastq2s = ""
                for file_name in os.listdir(sample_dir):
                    if file_name.endswith("fq.gz"):
                        m = re.match("^(.+)_[12].fq.gz", file_name)
                        if m is not None:
                            prefix_fastq_file_name = m.group(1)
                            if prefix_fastq_file_name in files:
                                files[prefix_fastq_file_name] += 1
                            else:
                                files[prefix_fastq_file_name] = 1

                no_files = len(files.keys())

                for prefix_file_name in sorted(list(files.keys())):
                    novogene_fastq1s += f"{sample_dir}/{prefix_file_name}_1.fq.gz "
                    novogene_fastq2s += f"{sample_dir}/{prefix_file_name}_2.fq.gz "

                run.add_sample(index, novogene_sample_id, sample_id, novogene_fastq1s, novogene_fastq2s, no_files)

    #run.print()

    # create directories in destination folder directory
    analysis_dir = f"{dest_dir}/analysis"
    contigs_dir = f"{dest_dir}/contigs"
    new_dir = ""
    try:
        os.makedirs(log_dir, exist_ok=True)
        new_dir = analysis_dir
        os.makedirs(new_dir, exist_ok=True)
        new_dir = contigs_dir
        os.makedirs(new_dir, exist_ok=True)
        for sample in run.samples:
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"
            os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/spades_result"
            os.makedirs(new_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)
    multiqc_dep = ""

    # analyze
    fastqc_multiqc_dep = ""
    kraken2_multiqc_dep = ""
    kraken2_reports = ""

    for idx, sample in enumerate(run.samples):
        # copy the files
        dst_fastq1 = ""
        dst_fastq2 = ""

        if sample.no_files == 1:
            src_fastq1 = f"{sample.novogene_fastq1s}"
            dst_fastq1 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz"
            tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
            dep = ""
            cmd = f"cp {src_fastq1} {dst_fastq1}"
            pg.add(tgt, dep, cmd)

            src_fastq2 = f"{sample.novogene_fastq2s}"
            dst_fastq2 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz"
            tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
            dep = ""
            cmd = f"cp {src_fastq2} {dst_fastq2}"
            pg.add(tgt, dep, cmd)

            sample.fastq1 = dst_fastq1
            sample.fastq2 = dst_fastq2
        else:
            src_fastq1 = f"{sample.novogene_fastq1s}"
            dst_fastq1 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz"
            tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
            dep = ""
            cmd = f"zcat {src_fastq1} | gzip > {dst_fastq1}"
            pg.add(tgt, dep, cmd)

            src_fastq2 = f"{sample.novogene_fastq2s}"
            dst_fastq2 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz"
            tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
            dep = ""
            cmd = f"zcat {src_fastq2} | gzip > {dst_fastq2}"
            pg.add(tgt, dep, cmd)

            sample.fastq1 = dst_fastq1
            sample.fastq2 = dst_fastq2

        # fastqc
        fastqc = "/usr/local/FastQC-0.11.9/fastqc"
        fastqc_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"

        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.err"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"{fastqc} {sample.fastq1} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.err"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        cmd = f"{fastqc} {sample.fastq2} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        # kraken2
        kraken2 = "/usr/local/kraken2-2.1.2/kraken2"
        kraken2_std_db = "/usr/local/ref/kraken2/20210908_standard"
        input_fastq_file1 = f"{sample.fastq1}"
        input_fastq_file2 = f"{sample.fastq2}"
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
        kraken2_reports += f" {log}"
        cmd = f"{kraken2} --db {kraken2_std_db} --paired {input_fastq_file1} {input_fastq_file2} --use-names --report {report_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # plot kronatools radial tree
        kt_import_taxonomy = "/usr/local/KronaTools-2.8.1/bin/ktImportTaxonomy"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        input_txt_file = f"{output_dir}/report.txt"
        output_html_file = f"{output_dir}/krona_radial_tree.html"
        log = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.krona_radial_tree.OK"
        cmd = f"{kt_import_taxonomy} -m 3 -t 5 {input_txt_file} -o {output_html_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # assemble
        # /usr/local/SPAdes-3.15.2/bin/spades.py -1 Siniae-1086-20_S3_L001_R1_001.fastq.gz -2 Siniae-1086-20_S3_L001_R2_001.fastq.gz -o 1086 --isolate
        spades = "/usr/local/SPAdes-3.15.4/bin/spades.py"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/spades_result/assembly"
        input_fastq_file1 = f"{sample.fastq1}"
        input_fastq_file2 = f"{sample.fastq2}"
        log = f"{log_dir}/{sample.idx}_{sample.id}.spades_assembly.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.spades_assembly.err"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK {log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.spades_assembly.OK"
        cmd = f"{spades} -1 {input_fastq_file1} -2 {input_fastq_file2} -o {output_dir} --isolate > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        #copy contigs to main directory
        src_fasta = f"{output_dir}/contigs.fasta"
        dst_fasta = f"{contigs_dir}/{run.idx}_{sample.idx}_{sample.id}.contigs.fasta"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.spades_assembly.OK"
        tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.contigs.fasta.OK"
        cmd = f"cp {src_fasta} {dst_fasta}"
        pg.add(tgt, dep, cmd)

        # evaluate assembly
        # Salmonella enterica Enteritidis
        # Salmonella enterica Typhimurium
        # Escherichia coli
        # Vibrio spp.
        # Streptococcus iniae
        # Streptococcus equi
        # Streptococcus suis
        # bp_genbank2gff --stdout --accession AM933172 > AM933172.gff
        # /usr/local/quast-5.2.0/quast.py -o quat_report -r /home/atks/analysis/20221101_vm_ngs/ref/e_coli/NC_000913.3.fasta contigs.fasta -g /home/atks/analysis/20221101_vm_ngs/ref/e_coli/U00096.gff
        # choose reference

    multiqc = "/usr/local/bin/multiqc"

    # plot fastqc multiqc results
    analysis = "fastqc"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -m fastqc -o {output_dir} -n {analysis} -dd -1 --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kraken2 multiqc results
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

    def __init__(self, idx, novogene_id, id, novogene_fastq1s, novogene_fastq2s, no_novogene_fastq_files):
        self.idx = idx
        self.novogene_id = novogene_id
        self.id = id
        self.novogene_fastq1s = novogene_fastq1s
        self.novogene_fastq2s = novogene_fastq2s
        self.no_files = no_novogene_fastq_files
        self.fastq1 = ""
        self.fastq2 = ""

    def print(self):
        print(f"index                    : {self.idx}")
        print(f"novogene id              : {self.novogene_id}")
        print(f"id                       : {self.id}")
        print(f"novogene_fastq1s         : {self.novogene_fastq1s}")
        print(f"novogene_fastq2s         : {self.novogene_fastq2s}")
        print(f"no_novogene_fastq_files  : {self.no_files}")


class Run(object):
    def __init__(self, id):
        m = re.match(r"\D+(\d+)", id)
        if m is not None:
            self.idx = m.group(1)
        else:
            self.idx = 999
        self.id = id
        self.samples = []

    def add_sample(self, idx, sample_id, novogene_fastq1s, novogene_fastq2s, no_novogene_files, fastq_infix):
        self.samples.append(Sample(idx, sample_id, novogene_fastq1s, novogene_fastq2s, no_novogene_files, fastq_infix))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id}")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")

if __name__ == "__main__":
    main() # type: ignore
