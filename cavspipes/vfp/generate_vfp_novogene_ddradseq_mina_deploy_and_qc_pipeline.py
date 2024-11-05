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
    default="novogene_ddradseq_mina_deploy_and_qc.mk",
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
@click.option("-g", "--genome_fasta_file", required=True, help="genome fasta file")
def main(make_file, run_id, novogene_illumina_dir, working_dir, sample_file, genome_fasta_file):
    """
    Moves Novogene Illumina fastq files to a destination and performs QC

    e.g. generate_vfp_novogene_ddradseq_mina_deploy_and_qc_pipeline -r mina1 -i raw -s mina1.sa -g /usr/local/ref/vfp/ManJav1.0_HiC.fasta
    """
    dest_dir = working_dir + "/" + run_id
    novogene_illumina_dir = os.path.abspath(novogene_illumina_dir)

    # read sample file
    ## novogene-sample-id fastq-infix
    ## ILLU2024	pangolin_radseq
    run = DDRadseq_Run(run_id, sample_file, novogene_illumina_dir)
    #run.print()

    print("\t{0:<21} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<21} :   {1:<10}".format("run_dir", run_id))
    print("\t{0:<21} :   {1:<10}".format("novogene_illumina_dir", novogene_illumina_dir))
    print("\t{0:<21} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<21} :   {1:<10}".format("sample_file", sample_file))
    print("\t{0:<21} :   {1:<10}".format("genome_fasta_file", genome_fasta_file))
    print("\t{0:<21} :   {1:<10}".format("dest_dir", dest_dir))

    # create directories in destination folder directory
    log_dir = f"{working_dir}/log"
    analysis_dir = f"{dest_dir}/analysis"
    trace_dir = f"{dest_dir}/trace"
    aux_dir = f"{working_dir}/aux"
    try:
        os.makedirs(aux_dir, exist_ok=True)
        os.makedirs(analysis_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
        for sample in run.novogene_samples:
            sample_dir = f"{aux_dir}/{sample.id}"
            os.makedirs(sample_dir, exist_ok=True)
            os.makedirs(f"{sample_dir}/demux", exist_ok=True)

        for sample in run.samples:
            sample_dir = f"{analysis_dir}/{sample.idx}_{sample.id}"
            os.makedirs(sample_dir, exist_ok=True)
            os.makedirs(f"{sample_dir}/fastqc_result", exist_ok=True)
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

    #programs
    fastqc = "/usr/local/FastQC-0.11.9/fastqc"
    trimmomatic = "java -jar /usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar PE"
    trimmomatic_trimmer = "ILLUMINACLIP:/usr/local/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36"
    process_radtags = "/usr/local/stacks-2.68/bin/process_radtags"
    multiqc = "docker run  -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc "
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    plot_bamstats = "/usr/local/samtools-1.17/bin/plot-bamstats"
    compute_effective_coverage = "/home/atks/programs/cavspipes/vfp/compute_effective_coverage.py"
    extract_general_stats = "/home/atks/programs/cavspipes/vfp/extract_general_stats.py"

    # initialize
    pg = PipelineGenerator(make_file)

    #trim and demultiplex
    for sample in run.novogene_samples:
        sample_dir = f"{aux_dir}/{sample.id}"
        #trim adaptors
        src_fastq1 = f"{sample.fastq1}"
        dst_fastq1 = f"{sample_dir}/trimmed_{sample.id}_R1.fastq.gz"
        dst_unpaired_fastq1 = f"{sample_dir}/unpaired_{sample.id}_R1.fastq.gz"
        src_fastq2 = f"{sample.fastq2}"
        dst_fastq2 = f"{sample_dir}/trimmed_{sample.id}_R2.fastq.gz"
        dst_unpaired_fastq2 = f"{sample_dir}/unpaired_{sample.id}_R2.fastq.gz"
        log = f"{log_dir}/{sample.id}.trimmomatic.log"
        err = f"{log_dir}/{sample.id}.trimmomatic.err"
        tgt = f"{log_dir}/{sample.id}.trimmomatic.OK"
        dep = f""
        cmd = f"{trimmomatic} {src_fastq1} {src_fastq2} {dst_fastq1} {dst_unpaired_fastq1} {dst_fastq2} {dst_unpaired_fastq2} {trimmomatic_trimmer} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        # write the barcode files
        # GCATG	BIOS0045
        # AACCA	BIOS0203
        # CGATC	BIOS0220
        # TCGAT	BIOS0472
        demux_dir = f"{sample_dir}/demux"
        demux_file = f"{sample_dir}/demux/demux.txt"
        with open(demux_file, "w") as f:
            for sub_sample in sample.samples:
                f.write(f"{sub_sample.barcode}\t{sub_sample.id}\n")

        #demultiplex
        #process_radtags -b ILLU2024_01_Idx01_demux.txt  -1 mina5/5_1_pangolin_ddradseq_R1.fastq.gz -2 mina5/5_1_pangolin_ddradseq_R2.fastq.gz  --inline-null --renz_1 ecoRI --renz_2 mspI -o demux -c -q  -r
        src_fastq1 = f"{sample_dir}/trimmed_{sample.id}_R1.fastq.gz"
        src_fastq2 = f"{sample_dir}/trimmed_{sample.id}_R2.fastq.gz"
        log = f"{log_dir}/{sample.id}.demux.log"
        err = f"{log_dir}/{sample.id}.demux.err"
        tgt = f"{log_dir}/{sample.id}.demux.OK"
        dep = f"{log_dir}/{sample.id}.trimmomatic.OK"
        cmd = f"{process_radtags} -b {demux_file} -1 {src_fastq1} -2 {src_fastq2} --inline-null --renz_1 ecoRI --renz_2 mspI -o {demux_dir} -c -q -r  > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        #copy the files to the destination directory
        for sub_sample in sample.samples:
            src_fastq_file = f"{sample_dir}/demux/{sub_sample.id}.1.fq.gz"
            dst_fastq_file = f"{dest_dir}/{run.idx}_{sub_sample.idx}_{sub_sample.id}_R1.fastq.gz"
            tgt = f"{log_dir}/{run.idx}_{sub_sample.idx}_{sub_sample.id}_R1.fastq.gz.OK"
            dep = f"{log_dir}/{sample.id}.demux.OK"
            cmd = f"mv {src_fastq_file} {dst_fastq_file}"
            pg.add(tgt, dep, cmd)

            src_fastq_file = f"{sample_dir}/demux/{sub_sample.id}.2.fq.gz"
            dst_fastq_file = f"{dest_dir}/{run.idx}_{sub_sample.idx}_{sub_sample.id}_R2.fastq.gz"
            tgt = f"{log_dir}/{run.idx}_{sub_sample.idx}_{sub_sample.id}_R2.fastq.gz.OK"
            dep = f"{log_dir}/{sample.id}.demux.OK"
            cmd = f"mv {src_fastq_file} {dst_fastq_file}"
            pg.add(tgt, dep, cmd)

    # analyze
    fastqc_multiqc_dep = ""
    samtools_multiqc_dep = ""

    #core analysis
    for sample in run.samples:
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

        input_fastq1 = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}_R1.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.err"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}_fastqc1.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK"
        cmd = f"{fastqc} {input_fastq1} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"

        input_fastq2 = f"{fastqc_dir}/{sample.padded_idx}_{sample.id}_R2.fastq.gz"
        log = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.log"
        err = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.err"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}_fastqc2.OK"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        cmd = f"{fastqc} {input_fastq2} -o {fastqc_dir} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)
        fastqc_multiqc_dep += f" {tgt}"


        #############################
        # align to reference assembly
        #############################
        align_dir = f"{analysis_dir}/{sample.idx}_{sample.id}/align_result"
        reference_fasta_file = genome_fasta_file

        #  align
        src_fastq1 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz"
        src_fastq2 = f"{dest_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz"
        output_bam_file = f"{align_dir}/{sample.idx}_{sample.id}.bam"
        log = f"{log_dir}/{sample.idx}_{sample.id}.align.log"
        sort_log = f"{log_dir}/{sample.idx}_{sample.id}.align.sort.log"
        dep = f"{log_dir}/{run.idx}_{sample.idx}_{sample.id}_R1.fastq.gz.OK {log_dir}/{run.idx}_{sample.idx}_{sample.id}_R2.fastq.gz.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.bam.OK"
        cmd = f"{bwa} mem -t 2 -M {reference_fasta_file} {src_fastq1} {src_fastq2} 2> {log} | {samtools} view -h | {samtools} sort -o {output_bam_file} 2> {sort_log}"
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

        # compute effective coverage
        input_stats_file = f"{align_dir}/coverage_stats/{sample.padded_idx}_{sample.id}.txt"
        output_stats_file = f"{align_dir}/coverage_stats/{sample.padded_idx}_{sample.id}.effective.coverage.stats.txt"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.coverage.stats.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.effective.coverage.stats.OK"
        cmd = f"{compute_effective_coverage} {input_stats_file} -o {output_stats_file} -s {sample.idx}_{sample.id}"
        pg.add(tgt, dep, cmd)

        #  stats
        output_stats_file = f"{align_dir}/general_stats/{sample.padded_idx}_{sample.id}.txt"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.bam.bai.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.stats.OK"
        cmd = f"{samtools} stats {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # extract general stats
        input_stats_file = f"{align_dir}/general_stats/{sample.padded_idx}_{sample.id}.txt"
        output_stats_file = f"{align_dir}/general_stats/{sample.padded_idx}_{sample.id}.extracted.stats.txt"
        dep = f"{log_dir}/{sample.idx}_{sample.id}.stats.OK"
        tgt = f"{log_dir}/{sample.idx}_{sample.id}.extracted.stats.OK"
        cmd = f"{extract_general_stats} {input_stats_file} -o {output_stats_file} -s {sample.idx}_{sample.id}"
        pg.add(tgt, dep, cmd)

        #  flag stats
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

    #plot fastqc multiqc results
    analysis = "fastqc"
    output_dir = f"{analysis_dir}/all/{analysis}"
    log = f"{log_dir}/{analysis}.multiqc_report.log"
    err = f"{log_dir}/{analysis}.multiqc_report.err"
    dep = fastqc_multiqc_dep
    tgt = f"{log_dir}/{analysis}.multiqc_report.OK"
    cmd = f"cd {analysis_dir}; {multiqc} . -m fastqc -f -o {output_dir} -n {analysis} --no-ansi > {log} 2> {err}"
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

    def __init__(self, idx, id, barcode, fastq1, fastq2):
        self.idx = idx
        self.padded_idx = f"{idx:02}"
        self.id = id
        self.barcode = barcode
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def print(self):
        print(f"  index           : {self.idx}")
        print(f"  sample id       : {self.id}")
        print(f"  barcode         : {self.barcode}")

class Novogene_Sample(object):
    def __init__(self, id, prefix_fastq_file_path):
        self.id = id
        self.dir = prefix_fastq_file_path
        self.fastq1 = f"{prefix_fastq_file_path}_1.fq.gz"
        self.fastq2 = f"{prefix_fastq_file_path}_2.fq.gz"
        self.samples = []

    def add_sample(self, idx, id, barcode):
        self.samples.append(Sample(idx, id, barcode, "", ""))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"sample_id   : {self.id}")
        print(f"directory   : {self.dir}")
        print(f"fastq1      : {self.fastq1}")
        print(f"fastq2      : {self.fastq2}")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")

class DDRadseq_Run(object):
    def __init__(self, id, sample_file, illumina_dir):
        m = re.match(r"\D+(\d+)", id)
        if m is not None:
            self.idx = m.group(1)
        else:
            self.idx = 999
        self.id = id
        self.novogene_samples = [] # for demultiplexing
        self.samples = []

        fastq_dir = ""
        for dirpath, dirnames, filenames in os.walk(illumina_dir):
            for dirname in dirnames:
                if dirname == "01.RawData":
                    fastq_dir = os.path.join(dirpath, dirname)
        fastq_dir = os.path.abspath(fastq_dir)

        novogene_sample_ids = {}

        with open(sample_file, "r") as file:
            novogene_sample_index = 0
            index = 0
            for line in file:
                if not line.startswith("#"):
                    index += 1
                    sample_id, barcode, novogene_sample_id = line.rstrip().split("\t")
                    #search directory for fastq files
                    sample_dir = f"{fastq_dir}/{novogene_sample_id}"
                    files = {}
                    prefix_fastq_file_name = ""
                    for file_name in os.listdir(sample_dir):
                        if file_name.endswith("fq.gz"):
                            m = re.match("^(.+)_[12].fq.gz", file_name)
                            if m is not None:
                                prefix_fastq_file_name = m.group(1)
                                if prefix_fastq_file_name in files:
                                    files[prefix_fastq_file_name] += 1
                                else:
                                    files[prefix_fastq_file_name] = 1

                    assert len(files.keys()) == 1

                    if prefix_fastq_file_name not in novogene_sample_ids:
                        novogene_sample_ids[prefix_fastq_file_name] = novogene_sample_index
                        self.novogene_samples.append(Novogene_Sample(novogene_sample_id, f"{sample_dir}/{prefix_fastq_file_name}"))
                        novogene_sample_index += 1

                    self.novogene_samples[novogene_sample_ids[prefix_fastq_file_name]].add_sample(index, sample_id, barcode)
                    self.samples.append(Sample(index, sample_id, barcode, f"{self.idx}_{index}_{sample_id}_R1.fastq.gz", f"{self.idx}_{index}_{sample_id}_R2.fastq.gz"))


    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id} {self.idx}")
        print(f"++++++++++++++++++++")
        print(f"novogene ids   ")
        print(f"++++++++++++++++++++")
        for sample in self.novogene_samples:
            sample.print()
        print(f"++++++++++++++++++++")
        print(f"sample ids    ")
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")


if __name__ == "__main__":
    main() # type: ignore
