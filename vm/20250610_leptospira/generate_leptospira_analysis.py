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
import random

@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_leptospira_analysis.mk",
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
    Leptospira characterisation  

    e.g. generate_leptospira_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample file", sample_file))

    # read sample file
    samples = {}
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                id, fastq1, fastq2, contigs = line.rstrip().split("\t")
                samples[id] = Sample(id, fastq1, fastq2, contigs)
            
    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    log_dir = f"{working_dir}/log"
    stats_dir = f"{working_dir}/stats"
    plot_dir = f"{working_dir}/plot"
    bam_dir = f"{working_dir}/bam"
    try:
        os.makedirs(ref_dir, exist_ok=True)
        # os.makedirs(log_dir, exist_ok=True)
        # os.makedirs(stats_dir, exist_ok=True)
        # os.makedirs(f"{stats_dir}/coverage", exist_ok=True)
        # os.makedirs(f"{stats_dir}/general", exist_ok=True)
        # os.makedirs(f"{stats_dir}/flag", exist_ok=True)
        # os.makedirs(f"{stats_dir}/idx", exist_ok=True)
        # os.makedirs(plot_dir, exist_ok=True)
        # os.makedirs(bam_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    ##########
    # programs
    ##########
    multiqc = "docker run  -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc "
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    plot_bamstats = "/usr/local/samtools-1.17/bin/plot-bamstats"
    compute_effective_coverage = "/home/atks/programs/cavspipes/vfp/compute_effective_coverage.py"
    extract_general_stats = "/home/atks/programs/cavspipes/vfp/extract_general_stats.py"
    filter_ddradseq_vcf = "/home/atks/programs/CAVS-pipelines/vfp/20241210_pangolin_ddradseq/filter_ddradseq_vcf.py"
    qwikplot = "/home/atks/programs/CAVS-pipelines/vfp/20241210_pangolin_ddradseq/qwikplot"
    
    #################
    # reference files
    #################
    # >NZ_CP043884.1 Leptospira interrogans serovar Canicola strain 782 chromosome 1, complete sequence
    # >NZ_CP043885.1 Leptospira interrogans serovar Canicola strain 782 chromosome 2, complete sequence
    # >NZ_CP043886.1 Leptospira interrogans serovar Canicola strain 782 plasmid p1, complete sequence
    # >NZ_CP043887.1 Leptospira interrogans serovar Canicola strain 782 plasmid p2, complete sequence
    # >NZ_CP043888.1 Leptospira interrogans serovar Canicola strain 782 plasmid p4, complete sequence
    # >NZ_CP043889.1 Leptospira interrogans serovar Canicola strain 782 plasmid p5, complete sequence
    # >NZ_CP043890.1 Leptospira interrogans serovar Canicola strain 782 plasmid p3, complete sequence
    # >NZ_CP044513.1 Leptospira interrogans serovar Canicola strain 611 chromosome 1, complete sequence
    # >NZ_CP044514.1 Leptospira interrogans serovar Canicola strain 611 chromosome 2, complete sequence
    # >NZ_CP044515.1 Leptospira interrogans serovar Canicola strain 611 plasmid p1, complete sequence
    # >NZ_CP044516.1 Leptospira interrogans serovar Canicola strain 611 plasmid p2, complete sequence
    # >NZ_CP044509.1 Leptospira interrogans serovar Canicola strain LJ178 chromosome 1, complete sequence
    # >NZ_CP044510.1 Leptospira interrogans serovar Canicola strain LJ178 chromosome 2, complete sequence
    # >NZ_CP044511.1 Leptospira interrogans serovar Canicola strain LJ178 plasmid p1, complete sequence
    # >NZ_CP044512.1 Leptospira interrogans serovar Canicola strain LJ178 plasmid p2, complete sequence
    # >NC_025197.1 Leptospira interrogans serovar Canicola plasmid pGui2, complete sequence
    # >NC_025136.1 Leptospira interrogans serovar Canicola plasmid pGui1, complete sequence
    # reference genome assembly for sunda pangolin
    # https://www.ncbi.nlm.nih.gov/assembly/GCF_014570535.1

    ref_seqs_acc_ids = [
        "NZ_CP043884.1",  # chromosome 1
        "NZ_CP043885.1",  # chromosome 2
        "NZ_CP043886.1",  # plasmid p1
        "NZ_CP043887.1",  # plasmid p2
        "NZ_CP043888.1",  # plasmid p3
        "NZ_CP043889.1",  # plasmid p4
        "NZ_CP043890.1"   # plasmid p5
        ]

    fasta_files =  ""
    fasta_files_dep = ""
    for acc_id in ref_seqs_acc_ids:
        output_fasta_file = f"{ref_dir}/{acc_id}.fasta"
        tgt = f"{output_fasta_file}.OK"
        fasta_files += f"{output_fasta_file} "
        fasta_files_dep += f"{tgt} "
        dep = ""
        cmd = f"efetch -db nuccore -id {acc_id} -format fasta > {output_fasta_file}"
        pg.add(tgt, dep, cmd)

    output_fasta_file = f"{ref_dir}/leptospira_interrogans_canicola.fasta"
    tgt = f"{output_fasta_file}.OK"
    dep = fasta_files_dep
    cmd = f"cat {fasta_files} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    #index reference sequence
    ref_fasta_file = f"{ref_dir}/leptospira_interrogans_canicola.fasta"
    log = f"{ref_dir}/bwa_index.log"
    tgt = f"{ref_dir}/bwa_index.OK"
    dep = f"{ref_fasta_file}.OK"
    cmd = f"{bwa} index -a bwtsw {ref_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

 
    # #filter short reads
    # input_fastq1_file = f"{fastq_dir}/{sample.id}.1.fq.gz"
    # input_fastq2_file = f"{fastq_dir}/{sample.id}.2.fq.gz"
    # output_fastq2_file = f"{fastq_dir}/{sample.id}.2.fq.gz"
    # tgt = f"{output_fastq2_file}.OK"
    # fastq_files_OK += f"{tgt} "
    # dep = f"{fastq_dir}/{sample.id}.1.fq.gz.OK {fastq_dir}/{sample.id}.2.fq.gz.OK"
    # cmd = f"{process_shortreads} -1 fastq/BIOS0006.1.fq.gz -2 fastq/BIOS0006.2.fq.gz -o {fastq_ulen_dir} --len-limit 145 "
    # pg.add(tgt, dep, cmd)

    # clean
    pg.add_clean(f"rm -fr {ref_dir}")

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

    def __init__(self, id, fastq1, fastq2, contigs):
        self.id = id
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.contigs = contigs

    def print(self):
        print(f"id       : {self.id}")
        print(f"fastq1   : {self.fastq1}")
        print(f"fastq2   : {self.fastq2}")
        print(f"contigs  : {self.contigs}")


if __name__ == "__main__":
    main() # type: ignore
