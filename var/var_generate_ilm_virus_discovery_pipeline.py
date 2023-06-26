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
import re
import sys
from shutil import copy2
from datetime import datetime


@click.command()
@click.option("-1", "--read1", required=True, help="read 1 FASTQ file")
@click.option("-2", "--read2", required=True, help="read 2 FASTQ file")
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="ilm_virus_discovery.mk",
    help="make file",
)
def main(read1, read2, working_dir, make_file):
    """
    Assembles virus sequences and annotates.

    e.g. var_generate_ilm_virus_discovery_pipeline.py -1 left.fastq.gz -2 right.fastq.gz
    """
    read1_fastq_file = read1
    read2_fastq_file = read2

    log_dir = f"{working_dir}/log"

    print("\t{0:<20} :   {1:<10}".format("read1", read1))
    print("\t{0:<20} :   {1:<10}".format("read2", read2))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))

    # create directories in destination folder directory
    analysis_dir = f"{working_dir}"
    new_dir = ""
    try:
        os.makedirs(log_dir, exist_ok=True)
        new_dir = f"{analysis_dir}/virus_reads"
        os.makedirs(new_dir, exist_ok=True)
        new_dir = f"{analysis_dir}/assembly"
        os.makedirs(new_dir, exist_ok=True)
        new_dir = f"{analysis_dir}/mapping"
        os.makedirs(new_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    # databases
    refseq_db = "/usr/local/ref/refseq/213/viral/refseq.213.viral.fasta.gz"

    # initialize
    pg = PipelineGenerator(make_file)

    # make diamond database
    # diamond makedb --in {uniprot_database} -d {args.output_directory}/nr > {args.output_directory}/nr.log 2> {args.output_directory}/nr.err\n")

    # filter for virus sequences
    # sensitivity modes
    # do a more comprehensive one later
    # ilm38/38_1_A220228_NDV_house34_tracheal_R1.fastq.gz
    # fast     376
    # normal   436
    # mid      476
    # more     484
    # very     488
    # ultra    488
    # iterate  436 (appears to be the default)

    # diamond  blastx -d uniprot90_virus.dmnd -q R1.fastq.gz -f 6 qseqid stitle -o out.m8
    diamond = "/usr/local/diamond-2.0.15/diamond"
    diamond_db = f"/usr/local/ref/uniprot/2021_02/virus/diamond/uniprot90_virus.dmnd"

    input_fastq_file = f"{read1_fastq_file}"
    output_alignment_txt_file = f"{analysis_dir}/virus_reads/r1.txt"
    log = f"{log_dir}/diamond_r1.log"
    err = f"{log_dir}/diamond_r1.err"
    tgt = f"{log_dir}/diamond_r1.OK"
    dep = ""
    cmd = f"{diamond} blastx -d {diamond_db} -q {input_fastq_file} -o {output_alignment_txt_file} -f 6 qseqid stitle > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    input_fastq_file = f"{read2_fastq_file}"
    output_alignment_txt_file = f"{analysis_dir}/virus_reads/r2.txt"
    log = f"{log_dir}/diamond_r2.log"
    err = f"{log_dir}/diamond_r2.err"
    tgt = f"{log_dir}/diamond_r2.OK"
    dep = ""
    cmd = f"{diamond} blastx -d {diamond_db} -q {input_fastq_file} -o {output_alignment_txt_file} -f 6 qseqid stitle > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # get uniq read IDs
    input_alignment_txt_file1 = f"{analysis_dir}/virus_reads/r1.txt"
    input_alignment_txt_file2 = f"{analysis_dir}/virus_reads/r2.txt"
    output_read_id_txt_file = f"{analysis_dir}/virus_reads/uniq_virus_read_id.txt"
    tgt = f"{log_dir}/uniq_virus_read_id.OK"
    dep = f"{log_dir}/diamond_r1.OK {log_dir}/diamond_r2.OK"
    cmd = f"cat {input_alignment_txt_file1} {input_alignment_txt_file2} | cut -f1 | sort | uniq > {output_read_id_txt_file}"
    pg.add(tgt, dep, cmd)

    # extract sequences
    seqtk = "/usr/local/seqtk-1.3/seqtk"

    input_read_id_txt_file = f"{analysis_dir}/virus_reads/uniq_virus_read_id.txt"
    input_fastq_file = f"{read1_fastq_file}"
    output_virus_fastq_file = f"{analysis_dir}/virus_reads/r1.fastq.gz"
    tgt = f"{log_dir}/r1.fastq.gz.OK"
    dep = f"{log_dir}/uniq_virus_read_id.OK"
    cmd = f"{seqtk} subseq {input_fastq_file} {input_read_id_txt_file} | gzip >  {output_virus_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_file = f"{read2_fastq_file}"
    output_virus_fastq_file = f"{analysis_dir}/virus_reads/r2.fastq.gz"
    tgt = f"{log_dir}/r2.fastq.gz.OK"
    dep = f"{log_dir}/uniq_virus_read_id.OK"
    cmd = f"{seqtk} subseq {input_fastq_file} {input_read_id_txt_file} | gzip > {output_virus_fastq_file}"
    pg.add(tgt, dep, cmd)

    # assemble filtered viral reads
    # Trinity --seqType fq --left virus_reads/r1.fastq.gz  --right virus_reads/r2.fastq.gz  --CPU 2 --max_memory 40G --output trinity
    trinity = "/usr/local/trinityrnaseq-v2.14.0/Trinity"
    input_fastq_file1 = f"{analysis_dir}/virus_reads/r1.fastq.gz"
    input_fastq_file2 = f"{analysis_dir}/virus_reads/r2.fastq.gz"
    assembly_dir = f"{analysis_dir}/assembly"
    tgt = f"{log_dir}/trinity.Trinity.fasta.OK"
    dep = f"{log_dir}/r1.fastq.gz.OK {log_dir}/r2.fastq.gz.OK"
    cmd = f"cd {assembly_dir}; {trinity} --seqType fq --left {input_fastq_file1} --right {input_fastq_file2} --CPU 2 --max_memory 40G --output trinity"
    pg.add(tgt, dep, cmd)

    # blastdb construction (manual)
    # makeblastdb -in refseq.214.viral.fasta -dbtype nucl -parse_seqids -out refseq.214.viral
    # copy blastdb to /usr/local/ref/refseq/214/viral

    # blast contigs to refseq
    # export BLASTDB=/usr/local/ref/refseq/214/viral/blastdb ; blastn -db refseq.214.viral  -query assembly/trinity.Trinity.fasta -outfmt 6 -max_target_seqs 1
    blastn = "/usr/local/ncbi-blast-2.12.0+/bin/blastn"
    blastdb_dir = "/usr/local/ref/refseq/214/viral/blastdb"
    input_fasta_file = f"{analysis_dir}/assembly/trinity.Trinity.fasta"
    blast_output_file = f"{analysis_dir}/blast/blast.output.psl"
    tgt = f"{log_dir}/blast.output.psl.OK"
    dep = f"{log_dir}/trinity.Trinity.fasta.OK"
    cmd = f'export BLASTDB={blastdb_dir}; {blastn} --seqType fa -query {input_fasta_file} -outfmt "6 stitle" -max_target_seqs 1 -out {blast_output_file}'
    pg.add(tgt, dep, cmd)

    #

    # map reads back to refseq hits to estimate coverage

    # map reads to assembly

    # stitch and complete refseq sequences

    # map reads to assembly

    # collect stats

    # generate report

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
        self.idx = id[3::]
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
