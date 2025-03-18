#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2025 Adrian Tan <adrian_tan@nparks.gov.sg>
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
    default="novogene_ilm_deploy_and_qc.mk",
    help="make file name",
)
@click.option("-r", "--run_id", required=True, help="Run ID")
@click.option("-i", "--novogene_illumina_dir", required=True, help="illumina directory")
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
    Moves Illumina fastq files to a destination and performs QC

    e.g. generate_var_novogene_ilm_deploy_and_qc_pipeline -r ilm23 -i raw -i ilm23.sa
    """
    dest_dir = working_dir + "/" + run_id
    illumina_dir = os.path.abspath(novogene_illumina_dir)
    fastq_dir = ""
    for dirpath, dirnames, filenames in os.walk(illumina_dir):
        for dirname in dirnames:
            if dirname == "01.RawData":
                fastq_dir = os.path.join(dirpath, dirname)

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<20} :   {1:<10}".format("illumina_dir", illumina_dir))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("dest_dir", dest_dir))
    print("\t{0:<20} :   {1:<10}".format("fastq_dir", fastq_dir))

    # read sample file
    ## novogene-sample-id fastq-infix
    ## ASFV_T1 1_A112_22-1_ASFV_spleen
    ## ASFV_T2 2_A113_22-15_ASFV_spleen
    ## ASFV_T3 3_A113_22-17_ASFV_spleen
    ## KHV_1   A050910_KHV_isolate
    run = Run(run_id)
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                novogene_sample_id, sample_id, virus = line.rstrip().split("\t")
                #search directory for fastq files
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

                run.add_sample(index, novogene_sample_id, sample_id, virus, novogene_fastq1s, novogene_fastq2s, no_files)

    # create directories in destination folder directory
    untrimmed_fastq_dir = f"{working_dir}/untrimmed_fastq"
    log_dir = f"{working_dir}/log"
    analysis_dir = f"{dest_dir}/analysis"
    trace_dir = f"{dest_dir}/trace"
    try:
        os.makedirs(untrimmed_fastq_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(analysis_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
        for sample in run.samples:
            sample_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(sample_dir, exist_ok=True)
            os.makedirs(f"{sample_dir}/kraken2_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/fastqc_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/spades_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/ref", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/general_stats", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/coverage_stats", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/flag_stats", exist_ok=True)
            os.makedirs(f"{sample_dir}/align_result/idx_stats", exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    #version
    version = "1.0.0"

    # programs
    fastqc = f"/usr/local/FastQC-0.12.1/fastqc --adapters /usr/local/FastQC-0.12.1/Configuration/adapter_list.illumina.txt"
    trimmomatic = "java -jar /usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar PE"
    trimmomatic_trimmer = "ILLUMINACLIP:/usr/local/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36"
    kraken2 = "/usr/local/kraken2-2.1.2/kraken2"
    kraken2_std_db = "/usr/local/ref/kraken2/20210908_standard"
    kt_import_taxonomy = "/usr/local/KronaTools-2.8.1/bin/ktImportTaxonomy"
    multiqc = "docker run  -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc "
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    plot_bamstats = "/usr/local/samtools-1.17/bin/plot-bamstats"

    virus_genomes = {
        "ASFV":"/usr/local/ref/var/FR682468.2.fasta",
        "ISKNV":"/usr/local/ref/var/NC_003494.1.fasta",
        "KHV":"/usr/local/ref/var/NC_009127.1.fasta",
        "NDV":"/usr/local/ref/var/NC_039223.1.fasta",
        "H1N1":"/usr/local/ref/var/H1N1.fasta",
        "H9N2":"/usr/local/ref/var/H9N2.fasta",
        "H5N1":"/usr/local/ref/var/H5N1.fasta",
        "H3N2":"/usr/local/ref/var/H3N2.fasta",
        "H2N2":"/usr/local/ref/var/H2N2.fasta",
        "H7N9":"/usr/local/ref/var/H7N9.fasta",
        "H1N1a":"/usr/local/ref/var/H1N1a.fasta",
        "H1N1b":"/usr/local/ref/var/H1N1b.fasta"
    }

    # initialize
    pg = PipelineGenerator(make_file)

    # multiqc dependencies
    fastqc_multiqc_dep = ""
    kraken2_multiqc_dep = ""
    samtools_multiqc_dep = ""
    quast_multiqc_dep = ""

    for sample in run.samples:

        # copy the files
        dst_fastq1 = ""
        dst_fastq2 = ""

        #combine files if necessary
        if sample.no_files == 1:
            src_fastq1 = f"{sample.novogene_fastq1s}"
            dst_fastq1 = f"{untrimmed_fastq_dir}/untrimmed_{sample.idx}_{sample.id}_R1.fastq.gz"
            tgt = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
            dep = ""
            cmd = f"cp {src_fastq1} {dst_fastq1}"
            pg.add(tgt, dep, cmd)

            src_fastq2 = f"{sample.novogene_fastq2s}"
            dst_fastq2 = f"{untrimmed_fastq_dir}/untrimmed_{sample.idx}_{sample.id}_R2.fastq.gz"
            tgt = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
            dep = ""
            cmd = f"cp {src_fastq2} {dst_fastq2}"
            pg.add(tgt, dep, cmd)
        else:
            src_fastq1 = f"{sample.novogene_fastq1s}"
            dst_fastq1 = f"{untrimmed_fastq_dir}/untrimmed_{sample.idx}_{sample.id}_R1.fastq.gz"
            tgt = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
            dep = ""
            cmd = f"zcat {src_fastq1} | gzip > {dst_fastq1}"
            pg.add(tgt, dep, cmd)

            src_fastq2 = f"{sample.novogene_fastq2s}"
            dst_fastq2 = f"{untrimmed_fastq_dir}/untrimmed_{sample.idx}_{sample.id}_R2.fastq.gz"
            tgt = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
            dep = ""
            cmd = f"zcat {src_fastq2} | gzip > {dst_fastq2}"
            pg.add(tgt, dep, cmd)

 #trim
        #java -jar /usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar PE
        #illu17/17_1_aquatic_160816_R1.fastq.gz
        #illu17/17_1_aquatic_160816_R2.fastq.gz
        #nextera_output_forward_paired.fq.gz
        #nextera_output_forward_unpaired.fq.gz
        #nextera_output_reverse_paired.fq.gz
        #nextera_output_reverse_unpaired.fq.gz
        #ILLUMINACLIP:/usr/local/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
        src_fastq1 = f"{untrimmed_fastq_dir}/untrimmed_{sample.idx}_{sample.id}_R1.fastq.gz"
        dst_fastq1 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz"
        dst_unpaired_fastq1 = f"{untrimmed_fastq_dir}/unpaired_{sample.idx}_{sample.id}_R1.fastq.gz"
        src_fastq2 = f"{untrimmed_fastq_dir}/untrimmed_{sample.idx}_{sample.id}_R2.fastq.gz"
        dst_fastq2 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz"
        dst_unpaired_fastq2 = f"{untrimmed_fastq_dir}/unpaired_{sample.idx}_{sample.id}_R2.fastq.gz"
        log = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.trimmomatic.log"
        err = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}.trimmomatic.err"
        tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        dep = f"{log_dir}/untrimmed_{sample.idx}_{sample.id}_R1.fastq.gz.OK {log_dir}/untrimmed_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        cmd = f"{trimmomatic} {src_fastq1} {src_fastq2} {dst_fastq1} {dst_unpaired_fastq1} {dst_fastq2} {dst_unpaired_fastq2} {trimmomatic_trimmer} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        #create second read target
        tgt = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"echo create read 2 target : {tgt}"
        pg.add(tgt, dep, cmd)

        sample.fastq1 = dst_fastq1
        sample.fastq2 = dst_fastq2

        # symbolic link for fastq files
        fastqc_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/fastqc_result"

        src_fastq1 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz"
        dst_fastq1 = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}_R1.fastq.gz"
        tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}_R1.fastq.gz.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"ln -sf {src_fastq1} {dst_fastq1}"
        pg.add(tgt, dep, cmd)

        src_fastq1 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz"
        dst_fastq1 = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}_R2.fastq.gz"
        tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}_R2.fastq.gz.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        cmd = f"ln -sf {src_fastq1} {dst_fastq1}"
        pg.add(tgt, dep, cmd)

        #fastqc
        input_fastq_file = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}_R1.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.err"
        tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}_fastqc1.OK"
        dep = f"{log_dir}/{sample.padded_idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"{fastqc} {input_fastq_file} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        input_fastq_file = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}_R2.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.err"
        tgt = f"{log_dir}/{sample.padded_idx}_{sample.id}_fastqc2.OK"
        dep = f"{log_dir}/{sample.padded_idx}_{sample.id}_R2.fastq.gz.OK"
        cmd = f"{fastqc} {input_fastq_file} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        # kraken2
        input_fastq_file1 = f"{src_fastq1}"
        input_fastq_file2 = f"{src_fastq2}"
        output_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/kraken2_result"
        report_file = f"{output_dir}/{sample.padded_idx}_{sample.id}.txt"
        log = f"{output_dir}/report.log"
        err = f"{output_dir}/run.log"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK {log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.kraken2.OK"
        kraken2_multiqc_dep += f" {tgt}"
        cmd = f"{kraken2} --db {kraken2_std_db} --paired {input_fastq_file1} {input_fastq_file2} --use-names --report {report_file} > {log} 2> {err}"
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

        ###################################
        # align to virus reference sequence
        ###################################
        align_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/align_result"
        ref_fasta_file_base_name = os.path.basename(virus_genomes[sample.virus])
        reference_fasta_file =  virus_genomes[sample.virus]

        #copy reference
        dep = ""
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.ref.fasta.OK"
        cmd = f"cp {reference_fasta_file} {align_dir}"
        pg.add(tgt, dep, cmd)

        # construct reference
        reference_fasta_file = f"{align_dir}/{ref_fasta_file_base_name}"
        log = f"{log_dir}/{sample.idx}_{sample.id}.ref.bwa_index.log"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.ref.fasta.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.ref.bwa_index.OK"
        cmd = f"{bwa} index -a bwtsw {reference_fasta_file} 2> {log}"
        pg.add(tgt, dep, cmd)

        # align
        output_bam_file = f"{align_dir}/{sample.idx}_{sample.id}.bam"
        log = f"{log_dir}/{sample.idx}_{sample.id}.align.log"
        sort_log = f"{log_dir}/{sample.idx}_{sample.id}.align.sort.log"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.ref.bwa_index.OK {log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK {log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.bam.OK"
        cmd = f"{bwa} mem -t 2 -M {reference_fasta_file} {sample.fastq1} {sample.fastq2} 2> {log} | {samtools} view -h | {samtools} sort -o {output_bam_file} 2> {sort_log}"
        pg.add(tgt, dep, cmd)

        #  index
        input_bam_file = f"{align_dir}/{sample.idx}_{sample.id}.bam"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
        cmd = f"{samtools} index {input_bam_file}"
        pg.add(tgt, dep, cmd)

        #  coverage
        input_bam_file = f"{align_dir}/{sample.idx}_{sample.id}.bam"
        output_stats_file = f"{align_dir}/coverage_stats/{sample.padded_idx}_{sample.id}.txt"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.coverage.stats.OK"
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # stats
        output_stats_file = f"{align_dir}/general_stats/{sample.padded_idx}_{sample.id}.txt"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.stats.OK"
        cmd = f"{samtools} stats {input_bam_file} > {output_stats_file}"
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


    # plot fastqc multiqc results
    analysis = "fastqc"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -m fastqc -f -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # plot kraken2 multiqc results
    analysis = "kraken2"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = kraken2_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -m kraken -f -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
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

    pg.add_clean(f"rm -fr {log_dir} {dest_dir}")

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

    def __init__(self, idx, novogene_id, id, virus, novogene_fastq1s, novogene_fastq2s, no_novogene_fastq_files):
        self.idx = idx
        self.padded_idx = f"{idx:02}"
        self.novogene_id = novogene_id
        self.id = id
        self.virus = virus
        self.novogene_fastq1s = novogene_fastq1s
        self.novogene_fastq2s = novogene_fastq2s
        self.no_files = no_novogene_fastq_files
        self.fastq1 = ""
        self.fastq2 = ""

    def print(self):
        print(f"index                    : {self.idx}")
        print(f"novogene id              : {self.novogene_id}")
        print(f"id                       : {self.id}")
        print(f"virus                    : {self.virus}")
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

    def add_sample(self, idx, sample_id, virus, novogene_fastq1s, novogene_fastq2s, no_novogene_files, fastq_infix):
        self.samples.append(Sample(idx, sample_id, virus, novogene_fastq1s, novogene_fastq2s, no_novogene_files, fastq_infix))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id}")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")


if __name__ == "__main__":
    main() # type: ignore
