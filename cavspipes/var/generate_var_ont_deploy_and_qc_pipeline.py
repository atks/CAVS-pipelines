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
import re
import sys
from shutil import copy2
from datetime import datetime


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
@click.option("-k", "--kit", default="SQK-NBD114-24", show_default=True, help="Kit ID")
#dna_r10.4.1_e8.2_400bps_sup@v5.0.0
@click.option("-y", "--basecall_model", default="dna_r10.4.1_e8.2_400bps_hac@v5.0.0", show_default=True)
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
    basecall_model,
    sample_file,
):
    """
    Moves ONT fastq files to a destination and performs QC

    e.g. generate_var_ont_deploy_and_qc_pipeline -r ont44 -i raw -s ont44.sa
    """
    log_dir = f"{working_dir}/log"
    dest_dir = working_dir + "/" + run_id
    fastq_path = f"{working_dir}/demux"
    basecall_model = f"/usr/local/dorado-0.8.3/models/{basecall_model}"

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<20} :   {1:<10}".format("nanopore_dir", nanopore_dir))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("minumum qscore", qscore))
    print("\t{0:<20} :   {1:<10}".format("minumum length", len))
    print("\t{0:<20} :   {1:<10}".format("memory", memory))
    print("\t{0:<20} :   {1:<10}".format("kit", kit))
    print("\t{0:<20} :   {1:<10}".format("basecall model", basecall_model))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("dest_dir", dest_dir))
    print("\t{0:<20} :   {1:<10}".format("fastq_path", fastq_path))

    #get acquisition run ID
    #print("\t{0:<20} :   {1:<10}".format("acquisition run id", acquisition_run_id))
    acquisition_run_id = ""
    dest_dir = working_dir + "/" + run_id
    nanopore_dir = os.path.abspath(nanopore_dir)
    for dirpath, dirnames, filenames in os.walk(nanopore_dir):
        for filename in filenames:
            if filename.startswith("final_summary_"):
                final_summary_txt_file = os.path.join(dirpath, filename)
                with open(final_summary_txt_file, "r") as file:
                    for line in file:
                        if line.startswith("acquisition_run_id"):
                            acquisition_run_id = line.rstrip().split("=")[1]
                            break
    print("\t{0:<20} :   {1:<10}".format("acquisition_run_id", acquisition_run_id))


    # version
    version = "1.0.0"

    # programs
    dorado = "/usr/local/dorado-0.8.3/bin/dorado"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    fastqc = f"/usr/local/FastQC-0.12.1/fastqc --adapters /usr/local/FastQC-0.12.1/Configuration/adapter_list.nanopore.txt --memory {memory}"
    multiqc = "docker run  -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc "
    nanoplot = f"docker run -u \"$$(id -u):$$(id -g)\" -t -v `pwd`:`pwd` -w `pwd` staphb/nanoplot:1.42.0 NanoPlot "
    kraken2 = "/usr/local/kraken2-2.1.2/kraken2"
    kraken2_std_db = "/usr/local/ref/kraken2/20220816_standard"
    kt_import_taxonomy = "/usr/local/KronaTools-2.8.1/bin/ktImportTaxonomy"
    ft = "/usr/local/cavstools-0.0.1/ft"
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    plot_bamstats = "/usr/local/samtools-1.17/bin/plot-bamstats"

    virus_genomes = {
        "ASFV":"/usr/local/ref/var/FR682468.2.fasta",
        "ISKNV":"/usr/local/ref/var/NC_003494.1.fasta",
        "KHV":"/usr/local/ref/var/NC_009127.1.fasta"
                     }

    run = Run(run_id)
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, barcode, virus = line.rstrip().split("\t")
                run.add_sample(index, sample_id, barcode, virus)

    # create directories in destination folder directory
    analysis_dir = f"{dest_dir}/analysis"
    trace_dir = f"{dest_dir}/trace"
    aux_dir = f"{working_dir}/aux"
    log_dir = f"{working_dir}/log"
    bam_dir = f"{working_dir}/bam"
    fastq_dir = f"{working_dir}/fastq"
    try:
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
        os.makedirs(analysis_dir, exist_ok=True)
        os.makedirs(aux_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        for sample in run.samples:
            sample_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(sample_dir, exist_ok=True)
            os.makedirs(f"{sample_dir}/kraken2_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/fastqc_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/ref", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/general_stats", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/coverage_stats", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/flag_stats", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/idx_stats", exist_ok=True)
            os.makedirs(f"{sample_dir}/nanoplot_result", exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # base call
    # dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.2.0  pod5s/ > calls.bam
    output_bam_file = f"{working_dir}/bam/basecalls.bam"
    pod5_files_dir = nanopore_dir
    log = f"{log_dir}/dorado_base_caller.log"
    err = f"{log_dir}/dorado_base_caller.err"
    tgt = f"{log_dir}/dorado_base_caller.OK"
    dep = ""
    cmd = f"{dorado} duplex -r {basecall_model} {pod5_files_dir} > {output_bam_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    # demux
    # trim barcode (default option)
    # dorado demux --output-dir demux --kit-name SQK-NBD114-24  basecalls.bam
    input_bam_file = f"{working_dir}/bam/basecalls.bam"
    output_dir = f"{working_dir}/demux"
    log = f"{log_dir}/demux.log"
    err = f"{log_dir}/demux.err"
    dep = f"{log_dir}/dorado_base_caller.OK"
    tgt = f"{log_dir}/demux.OK"
    cmd = f'{dorado} demux --output-dir {output_dir} --kit-name {kit} {input_bam_file} > {log} 2> {err}'
    pg.add(tgt, dep, cmd)

    fastqc_multiqc_dep = ""
    nanoplot_multiqc_dep = ""
    kraken2_multiqc_dep = ""
    samtools_multiqc_dep = ""

    for sample in run.samples:

        # samtools bam2fq <acquisition_run_id>_SQK-NBD114-24_barcode01.bam  -T "*"
        input_bam_file = f"{working_dir}/demux/{acquisition_run_id}_{kit}_{sample.barcode}.bam"
        if sample.barcode=="unclassified":
            input_bam_file = f"{working_dir}/demux/{acquisition_run_id}_unclassified.bam"
        output_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        log = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz.log.OK"
        tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz.OK"
        dep = f"{log_dir}/demux.OK"
        cmd = f"{samtools} bam2fq -T \"*\" {input_bam_file} 2> {log} | gzip > {output_fastq_file} "
        pg.add(tgt, dep, cmd)

        # # filter
        # input_fastq_file = f"{trimmed_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        # output_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        # tgt = f"{log_dir}/{sample.idx}_{sample.id}.trimmed.filtered.fastq.gz.OK"
        # dep = f"{log_dir}/demux.trimmed.OK"
        # cmd = f"{ft} filter -q {qscore} -l {len} -o {output_fastq_file} {input_fastq_file}"
        # pg.add(tgt, dep, cmd)

        # symbolic link for fastqc
        fastqc_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"

        src_fastq = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        dst_fastq = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}.fastq.gz"
        tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastq.gz.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz.OK"
        cmd = f"ln -sf {src_fastq} {dst_fastq}"
        pg.add(tgt, dep, cmd)

        # fastqc
        input_file = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"
        log = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastqc.log"
        err = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastqc.err"
        dep = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastq.gz.OK"
        tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastqc.OK"
        fastqc_multiqc_dep += f" {tgt}"
        cmd = f"{fastqc} {input_file} -o {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # nanoplot
        input_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/nanoplot_result"
        prefix = f"{sample.padded_idx}_{sample.id}_"
        log = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.err"
        dep = f"{log_dir}/{sample.padded_idx}_{sample.id}.fastq.gz.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.run.OK"
        cmd = f"{nanoplot} --fastq {input_file} -p {prefix} -o {output_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        #rename nanoplot file
        src_txt_file = f"{analysis_dir}/{sample.idx}_{sample.id}/nanoplot_result/{sample.padded_idx}_{sample.id}_NanoStats.txt"
        dst_txt_file = f"{analysis_dir}/{sample.idx}_{sample.id}/nanoplot_result/{sample.padded_idx}_{sample.id}.txt"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.run.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.run.renamed.OK"
        nanoplot_multiqc_dep += f" {tgt}"
        cmd = f"mv {src_txt_file} {dst_txt_file}"
        pg.add(tgt, dep, cmd)

        # #remove post filtering files from nanoplot
        # txt_file = f"{analysis_dir}/{sample.idx}_{sample.id}/nanoplot_result/{sample.padded_idx}_{sample.id}_NanoStats_post_filtering.txt"
        # dep = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.run.renamed.OK"
        # tgt = f"{log_dir}/{sample.idx}_{sample.id}.nanoplot.OK"
        # nanoplot_multiqc_dep += f" {tgt}"
        # cmd = f"rm {txt_file}"
        # pg.add(tgt, dep, cmd)

        if sample.virus!="n/a":

            ###########################
            # align to reference genome
            ###########################
            align_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/align_result"
            ref_fasta_file_base_name = os.path.basename(virus_genomes[sample.virus])
            reference_fasta_file = virus_genomes[sample.virus]

            #copy reference
            dep = ""
            src_fasta_file = f"{virus_genomes[sample.virus]}"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.align.ref.OK"
            cmd = f"cp {src_fasta_file} {align_dir}"
            pg.add(tgt, dep, cmd)

            # construct reference
            reference_fasta_file = f"{align_dir}/{ref_fasta_file_base_name}"
            log = f"{log_dir}/{sample.idx}_{sample.id}.align.ref.mmi.log"
            dep = f"{log_dir}/{sample.idx}_{sample.id}.align.ref.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.align.ref.mmi.OK"
            cmd = f"{minimap2} -d {reference_fasta_file}.mmi {reference_fasta_file} 2> {log}"
            pg.add(tgt, dep, cmd)

            # align
            input_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
            output_bam_file = f"{align_dir}/{sample.idx}_{sample.id}.bam"
            log = f"{log_dir}/{sample.idx}_{sample.id}.align.log"
            sort_log = f"{log_dir}/{sample.idx}_{sample.id}.align.sort.log"
            dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz.OK {log_dir}/{sample.idx}_{sample.id}.align.ref.mmi.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.bam.OK"
            cmd = f"{minimap2} -ax map-ont {reference_fasta_file}.mmi {input_fastq_file} 2> {log} | {samtools} view -h | {samtools} sort -o {output_bam_file} 2> {sort_log}"
            pg.add(tgt, dep, cmd)

            #  index
            input_bam_file = f"{align_dir}/{sample.idx}_{sample.id}.bam"
            dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
            cmd = f"{samtools} index {input_bam_file}"
            pg.add(tgt, dep, cmd)

            # stats
            output_stats_file = f"{align_dir}/general_stats/{sample.padded_idx}_{sample.id}.txt"
            dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.stats.OK"
            cmd = f"{samtools} stats {input_bam_file} > {output_stats_file}"
            samtools_multiqc_dep += f" {tgt}"
            pg.add(tgt, dep, cmd)

            # coverage stats
            output_stats_file = f"{align_dir}/coverage_stats/{sample.padded_idx}_{sample.id}.txt"
            dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.coverage.stats.OK"
            cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
            samtools_multiqc_dep += f" {tgt}"
            pg.add(tgt, dep, cmd)

            # flag stats
            output_stats_file = f"{align_dir}/flag_stats/{sample.padded_idx}_{sample.id}.txt"
            dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.flag.stats.OK"
            cmd = f"{samtools} flagstat {input_bam_file} > {output_stats_file}"
            samtools_multiqc_dep += f" {tgt}"
            pg.add(tgt, dep, cmd)

            #  idx stats
            output_stats_file = f"{align_dir}/idx_stats/{sample.padded_idx}_{sample.id}.txt"
            dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.idx.stats.OK"
            cmd = f"{samtools} idxstats {input_bam_file} > {output_stats_file}"
            samtools_multiqc_dep += f" {tgt}"
            pg.add(tgt, dep, cmd)

            # plot samtools stats
            input_stats_file = f"{align_dir}/general_stats/{sample.padded_idx}_{sample.id}.txt"
            dep = f"{log_dir}/{sample.idx}_{sample.id}.stats.OK"
            tgt = f"{log_dir}/{sample.idx}_{sample.id}.plot_bamstats.OK"
            cmd = f"{plot_bamstats} -p  {align_dir}/plot_bamstats/plot {input_stats_file}"
            pg.add(tgt, dep, cmd)

        # kraken2
        input_fastq_file = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        report_file = f"{output_dir}/{sample.padded_idx}_{sample.id}.txt"
        log = f"{output_dir}/{sample.padded_idx}_{sample.id}.log"
        err = f"{output_dir}/{sample.padded_idx}_{sample.id}.err"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.fastq.gz.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        kraken2_multiqc_dep += f" {tgt}"
        cmd = f"{kraken2} --db {kraken2_std_db} {input_fastq_file} --use-names --report {report_file} > {log} 2> {err}"
        pg.add_srun(tgt, dep, cmd, 15)

        # plot kronatools radial tree
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        input_txt_file = f"{output_dir}/{sample.padded_idx}_{sample.id}.txt"
        output_html_file = f"{output_dir}/krona_radial_tree.html"
        log = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}.krona_radial_tree.err"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.krona_radial_tree.OK"
        cmd = f"{kt_import_taxonomy} -m 3 -t 5 {input_txt_file} -o {output_html_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

    # plot fastqc multiqc results
    analysis = "fastqc"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -f -m fastqc -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot nanoplot multiqc results
    analysis = "nanoplot"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = nanoplot_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -f -m nanostat -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kraken2 multiqc results
    analysis = "kraken2"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = kraken2_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -f -m kraken -o {output_dir} -n {analysis} -dd -1 --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot samtools
    analysis = "samtools"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = samtools_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -m samtools -f -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # write make file
    print("Writing pipeline")
    pg.write()

    #copy files to trace
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


class Sample(object):
    def __init__(self, idx, id, barcode, virus):
        self.idx = idx
        self.padded_idx = f"{idx:02}"
        self.id = id
        self.barcode = barcode
        self.fastq = ""
        self.virus = virus

    def print(self):
        print(f"index   : {self.idx}")
        print(f"id      : {self.id}")
        print(f"barcode : {self.barcode}")
        print(f"fastq   : {self.fastq}")
        print(f"virus   : {self.virus}")

class Run(object):
    def __init__(self, id):
        m = re.match(r"\D+(\d+)", id)
        if m is not None:
            self.idx = m.group(1)
        else:
            self.idx = 999
        self.id = id
        self.samples = []

    def add_sample(self, idx, sample_id, barcode, virus):
        self.samples.append(Sample(idx, sample_id, barcode, virus))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id}")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")

if __name__ == "__main__":
    main()  # type: ignore
