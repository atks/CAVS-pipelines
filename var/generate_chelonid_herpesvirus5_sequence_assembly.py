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
    default="run_amba_vera_analysis.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default="/home/atks/analysis/20220725_amba_vera",
    show_default=True,
    help="working directory",
)
def main(make_file, working_dir):
    """
    Amba and Vera analysis

    e.g. generate_amba_vera_analysis
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    try:
        os.makedirs(f"{working_dir}/data", exist_ok=True)
        os.makedirs(f"{working_dir}/data", exist_ok=True)
        os.makedirs(f"{working_dir}/data", exist_ok=True)
        os.makedirs(f"{working_dir}/data", exist_ok=True)
    except OSError as error:
        print(f"Directory {working_dir} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    data_dir = f"{working_dir}/data"
    log_dir = f"{working_dir}/log"
    try:
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {data_dir} cannot be created")

    # data directory
    # Wuhan reference NC_045512.fasta
    # singapore reference EPI_ISL_6600690.fasta

    # setup pooled sequences
    # zcat ont25/25_1_flipper_tumour_ChHV5.fastq.gz ont25/25_2_neck_tumour_ChHV5.fastq.gz | gzip > tissue.fastq.gz
    # zcat ont25/25_3_pump_12_ChHV5.fastq.gz ont25/25_4_pump_13_ChHV5.fastq.gz ont25/25_5_sterivex_14_ChHV5.fastq.gz ont25/25_6_sterivex_15_ChHV5.fastq.gz | gzip > water.fastq.gz
    # minimap2 -d HQ878327.mmi HQ878327.fasta
    # minimap2 -ax map-ont HQ878327.fasta tissue.fastq.gz  | samtools sort | samtools view -o tissue.bam
    # minimap2 -ax map-ont HQ878327.fasta water.fastq.gz  | samtools sort | samtools view -o water.bam
    # samtools index tissue.bam
    # samtools index water.bam
    # samtools coverage tissue.bam | transpose
    # samtools coverage water.bam | transpose
    # minimap2 -ax map-ont AB012104.1.fasta tissue.fastq.gz | samtools sort | samtools view -o tissue.mito.bam
    # minimap2 -ax map-ont AB012104.1.fasta water.fastq.gz | samtools sort | samtools view -o water.mito.bam
    # samtools index tissue.mito.bam
    # samtools coverage tissue.mito.bam | transpose
    # samtools index water.mito.bam
    # samtools coverage water.mito.bam | transpose

    # ont26
    # minimap2 -ax map-ont AB012104.1.fasta ont26/26_8_turtle_2020.fastq.gz | samtools sort | samtools view -o turtle.2020.mito.bam
    # samtools coverage turtle.2020.mito.bam | transpose
    # rname	AB012104.1
    # startpos	1
    # endpos	16497
    # numreads	35
    # covbases	2727
    # coverage	16.5303
    # meandepth	1.03516
    # meanbaseq	21.9
    # meanmapq	14.1

    # minimap2 -ax map-ont AB012104.1.fasta ont26/26_9_turtle_2021.fastq.gz | samtools sort | samtools view -o turtle.20201.mito.bam
    # samtools coverage turtle.2021.mito.bam | transpose
    # rname	AB012104.1
    # startpos	1
    # endpos	16497
    # numreads	25
    # covbases	4116
    # coverage	24.95
    # meandepth	0.901194
    # meanbaseq	21.3
    # meanmapq	13.7

    # minimap2 -ax map-ont AB012104.1.fasta ont26/26_6_M221254_EIV_T268.fastq.gz | samtools sort | samtools view -o eiv.t268.mito.bam
    # samtools coverage eiv.t268.mito.bam | transpose
    # rname	AB012104.1
    # startpos	1
    # endpos	16497
    # numreads	28
    # covbases	3802
    # coverage	23.0466
    # meandepth	0.779657
    # meanbaseq	20.6
    # meanmapq	21.6

    # minimap2 -ax map-ont AB012104.1.fasta ont26/26_7_M221254_EIV_T269.fastq.gz | samtools sort | samtools view -o eiv.t269.mito.bam
    # samtools coverage eiv.t269.mito.bam | transpose
    # rname	AB012104.1
    # startpos	1
    # endpos	16497
    # numreads	40
    # covbases	4546
    # coverage	27.5565
    # meandepth	1.34364
    # meanbaseq	21.4
    # meanmapq	20.9

    # minimap2 -ax map-ont AB012104.1.fasta ont26/26_13_unclassified.fastq.gz | samtools sort | samtools view -o unclassified.mito.bam
    # samtools coverage unclassified.mito.bam | transpose
    # rname	AB012104.1
    # startpos	1
    # endpos	16497
    # numreads	119
    # covbases	7233
    # coverage	43.8443
    # meandepth	3.82706
    # meanbaseq	20.2
    # meanmapq	15.9

    # zcat ont26/*.gz | gzip > all.fastq.gz
    # minimap2 -ax map-ont AB012104.1.fasta all.fastq.gz | samtools sort | samtools view -o all.mito.bam
    # samtools coverage all.mito.bam | transpose
    # rname	AB012104.1
    # startpos	1
    # endpos	16497
    # numreads	716
    # covbases	8676
    # coverage	52.5914
    # meandepth	22.5721
    # meanbaseq	20.9
    # meanmapq	16.9

    # minimap2 -d AB012104.1.mmi AB012104.1.fasta

    # mini_align -i tissue.fastq.gz -r HQ878327.fasta -m -p tissue.chhv5.minialign -t 2
    # medaka consensus tissue.minialign.bam results/out.hdf --model r941_min_hac_g507
    # medaka stitch results/out.hdf HQ878327.fasta tissue.polished.fasta

    # mini_align -i tissue.fastq.gz -r AB012104.1.fasta -m -p tissue.mito.minialign -t 2
    # medaka consensus tissue.mito.minialign.bam results/mito.hdf --model r941_min_hac_g507
    # medaka stitch results/mito.hdf AB012104.1.fasta tissue.mito.polished.fasta

    # mini_align -i water.fastq.gz -r HQ878327.fasta -m -p water.chhv5.minialign -t 2
    # medaka consensus water.chhv5.minialign.bam results/water.chhv5.hdf --model r941_min_hac_g507
    # medaka stitch results/water.chhv5.hdf HQ878327.fasta water.chhv5.polished.fasta

    # combine all amba related sequences
    fastq_file1 = (
        f"/net/singapura/var/hts/ont15/15_1_M211147_SARSCoV2_v1_nasal.fastq.gz"
    )
    fastq_file2 = (
        f"/net/singapura/var/hts/ont15/15_2_M211147_SARSCoV2_v3_nasal.fastq.gz"
    )
    output_fastq_file = f"{working_dir}/data/amba.fastq.gz"
    err = f"{log_dir}/amba.fastq.gz.err"
    dep = ""
    tgt = f"{log_dir}/amba.fastq.gz.OK"
    cmd = f"zcat {fastq_file1} {fastq_file2} | gzip > {output_fastq_file} 2> err"
    pg.add(tgt, dep, cmd)

    # combine all amba related sequences
    fastq_file1 = (
        f"/net/singapura/var/hts/ont18/18_1_M211147_SARSCoV2_v1_nasal.fastq.gz"
    )
    fastq_file2 = (
        f"/net/singapura/var/hts/ont18/18_2_M211147_SARSCoV2_v3_nasal.fastq.gz"
    )
    fastq_file3 = f"/net/singapura/var/hts/ont18/18_3_negative_control.fastq.gz"
    fastq_file4 = f"/net/singapura/var/hts/ont18/18_4_unclassified.fastq.gz"
    output_fastq_file = f"{working_dir}/data/vera.fastq.gz"
    err = f"{log_dir}/vera.fastq.gz.err"
    dep = ""
    tgt = f"{log_dir}/vera.fastq.gz.OK"
    cmd = f"zcat {fastq_file1} {fastq_file2} {fastq_file3} {fastq_file4} | gzip > {output_fastq_file} 2> {err}"
    pg.add(tgt, dep, cmd)

    # mapping with minimap2
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    samtools = "/usr/local/samtools-1.15/samtools"

    # EPI_ISL_6600690 database
    database_file = f"{working_dir}/minimap2/EPI_ISL_6600690.mmi"
    input_fasta_file = f"{working_dir}/data/EPI_ISL_6600690.fasta"
    output_fasta_file = f"{working_dir}/database/EPI_ISL_6600690.fasta"
    dep = ""
    tgt = f"{log_dir}/EPI_ISL_6600690.fasta.OK"
    cmd = f"ln -fs {input_fasta_file} {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    database_file = f"{working_dir}/minimap2/EPI_ISL_6600690.mmi"
    input_fasta_file = f"{working_dir}/database/EPI_ISL_6600690.fasta"
    log = f"{log_dir}/EPI_ISL_6600690.mmi.log"
    err = f"{log_dir}/EPI_ISL_6600690.mmi.err"
    dep = f"{log_dir}/EPI_ISL_6600690.fasta.OK"
    tgt = f"{log_dir}/EPI_ISL_6600690.mmi.OK"
    cmd = f"{minimap2} -d {database_file} {input_fasta_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # map amba
    database_file = f"{working_dir}/database/EPI_ISL_6600690.mmi"
    ref_fasta_file = f"{working_dir}/database/EPI_ISL_6600690.fasta"
    input_fastq_file = f"{working_dir}/data/amba.fastq.gz"
    output_bam_file = f"{working_dir}/minimap2/amba.bam"
    log = f"{log_dir}/amba.minimap2.log"
    err = f"{log_dir}/amba.minimap2.err"
    dep = f"{log_dir}/EPI_ISL_6600690.mmi.OK"
    tgt = f"{log_dir}/amba.bam.OK"
    cmd = f"{minimap2} -ax map-ont {ref_fasta_file} {input_fastq_file} | {samtools} sort | {samtools} view -o {output_bam_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # index
    input_bam_file = f"{working_dir}/minimap2/amba.bam"
    dep = f"{log_dir}/amba.bam.OK"
    tgt = f"{log_dir}/amba.bam.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    # map vera
    database_file = f"{working_dir}/database/EPI_ISL_6600690.mmi"
    ref_fasta_file = f"{working_dir}/database/EPI_ISL_6600690.fasta"
    input_fastq_file = f"{working_dir}/data/vera.fastq.gz"
    output_bam_file = f"{working_dir}/minimap2/vera.bam"
    log = f"{log_dir}/vera.minimap2.log"
    err = f"{log_dir}/vera.minimap2.err"
    dep = f"{log_dir}/EPI_ISL_6600690.mmi.OK"
    tgt = f"{log_dir}/vera.bam.OK"
    cmd = f"{minimap2} -ax map-ont {ref_fasta_file} {input_fastq_file} | {samtools} sort | {samtools} view -o {output_bam_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # index
    input_bam_file = f"{working_dir}/minimap2/vera.bam"
    dep = f"{log_dir}/vera.bam.OK"
    tgt = f"{log_dir}/vera.bam.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    # blast
    makeblastdb = "/usr/local/ncbi-blast-2.12.0+/bin/makeblastdb"
    blastn = "/usr/local/ncbi-blast-2.12.0+/bin/blastn"

    # make blastdb
    database_dir = f"{working_dir}/database"
    blast_database_name = f"EPI_ISL_6600690"
    input_fasta_file = f"{working_dir}/database/EPI_ISL_6600690.fasta"
    log = f"{log_dir}/EPI_ISL_6600690.blastdb.log"
    err = f"{log_dir}/EPI_ISL_6600690.blastdb.err"
    dep = f"{log_dir}/EPI_ISL_6600690.fasta.OK"
    tgt = f"{log_dir}/EPI_ISL_6600690.blastdb.OK"
    cmd = f"cd {database_dir}; {makeblastdb} -in {input_fasta_file} -out {blast_database_name} -dbtype nucl -input_type fasta > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # amba contig alignment against reference
    blast_database_name = f"EPI_ISL_6600690"
    blast_database_dir = f"{working_dir}/database"
    input_fasta_file = f"{working_dir}/gd/amba_contigs.fasta"
    output_blast_file = f"{working_dir}/blast/amba_ref_alignment.txt"
    log = f"{log_dir}/amba.blast.log"
    err = f"{log_dir}/amba.blast.err"
    dep = f"{log_dir}/EPI_ISL_6600690.blastdb.OK"
    tgt = f"{log_dir}/amba.blast.OK"
    cmd = f'cd {database_dir}; {blastn} -db {blast_database_name} -query {input_fasta_file} -max_target_seqs 1 -outfmt "6 qacc sacc qlen slen sstart send score length pident stitle" -out {output_blast_file} > {log} 2> {err}'
    pg.add(tgt, dep, cmd)

    # vera contig alignment against reference
    blast_database_name = f"EPI_ISL_6600690"
    blast_database_dir = f"{working_dir}/database"
    input_fasta_file = f"{working_dir}/gd/vera_contigs.fasta"
    output_blast_file = f"{working_dir}/blast/vera_ref_alignment.txt"
    log = f"{log_dir}/vera.blast.log"
    err = f"{log_dir}/vera.blast.err"
    dep = f"{log_dir}/EPI_ISL_6600690.blastdb.OK"
    tgt = f"{log_dir}/vera.blast.OK"
    cmd = f'cd {database_dir}; {blastn} -db {blast_database_name} -query {input_fasta_file} -max_target_seqs 1 -outfmt "6 qacc sacc qlen slen sstart send score length pident stitle" -out {output_blast_file} > {log} 2> {err}'
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


if __name__ == "__main__":
    main()
