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


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_chhv5_chmy_mito_sequence_assembly_pipeline.mk",
    help="make file name",
)
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
def main(make_file, working_dir):
    """
    ChHV5 and Chelonia mydas mitchondria sequence assembly

    e.g. generate_chhv5_chmy_mito_sequence_assembly.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    fastq_dir = f"{working_dir}/seq"
    fasta_dir = f"{working_dir}/fasta"
    bam_dir = f"{working_dir}/bam"
    assembly_dir = f"{working_dir}/assembly"
    stats_dir = f"{working_dir}/stats"
    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(assembly_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    #################
    # reference files
    #################
    id = "HQ878327.2"
    output_genbank_file = f"{ref_dir}/{id}.genbank"
    dep = ""
    tgt = f"{output_genbank_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format genbank > {output_genbank_file}"
    pg.add(tgt, dep, cmd)

    output_fasta_file = f"{ref_dir}/{id}.fasta"
    dep = ""
    tgt = f"{output_fasta_file}.OK"
    cmd = f"efetch -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    id = "JQ034420.1"
    output_genbank_file = f"{ref_dir}/{id}.genbank"
    dep = ""
    tgt = f"{output_genbank_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format genbank > {output_genbank_file}"
    pg.add(tgt, dep, cmd)

    output_fasta_file = f"{ref_dir}/{id}.fasta"
    dep = ""
    tgt = f"{output_fasta_file}.OK"
    cmd = f"efetch -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    id = "AB012104.1"
    output_genbank_file = f"{ref_dir}/{id}.genbank"
    dep = ""
    tgt = f"{output_genbank_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format genbank > {output_genbank_file}"
    pg.add(tgt, dep, cmd)

    output_fasta_file = f"{ref_dir}/{id}.fasta"
    dep = ""
    tgt = f"{output_fasta_file}.OK"
    cmd = f"efetch -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    ###############
    # sequence file
    ###############

    # ***********
    # 2020 turtle
    # ***********

    # ilm31
    input_fastq_file = (
        f"/net/singapura/var/hts/ilm31/31_4_A20068_ChHV5_tumour_R1.fastq.gz"
    )
    output_fastq_file = f"{fastq_dir}/2020_turtle_ilm31_R1.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_file} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_file = (
        f"/net/singapura/var/hts/ilm31/31_4_A20068_ChHV5_tumour_R2.fastq.gz"
    )
    output_fastq_file = f"{fastq_dir}/2020_turtle_ilm31_R2.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_file} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    # ont10
    input_fastq_file = (
        f"/net/singapura/var/hts/ont10/10_6_A200608_ChHV5_tissue.fastq.gz"
    )
    output_fastq_file = f"{fastq_dir}/2020_turtle_ont10.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_file} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    # ont26
    input_fastq_file = f"/net/singapura/var/hts/ont26/26_8_turtle_2020.fastq.gz "
    output_fastq_file = f"{fastq_dir}/2020_turtle_ont26.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_file} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    # ***********
    # 2021 turtle
    # ***********

    # ilm35
    input_fastq_file = (
        f"/net/singapura/var/hts/ilm35/35_4_A210806_ChHV5_tumour_R1.fastq.gz"
    )
    output_fastq_file = f"{fastq_dir}/2021_turtle_ilm35_R1.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_file} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_file = (
        f"/net/singapura/var/hts/ilm35/35_4_A210806_ChHV5_tumour_R2.fastq.gz"
    )
    output_fastq_file = f"{fastq_dir}/2021_turtle_ilm35_R2.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_file} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    # ont26
    input_fastq_file = f"/net/singapura/var/hts/ont26/26_9_turtle_2021.fastq.gz "
    output_fastq_file = f"{fastq_dir}/2021_turtle_ont26.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_file} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    # ************************
    # 2020 and 2021 turtle mix
    # ************************
    input_fastq_file = f"/net/singapura/var/hts/ont26/*.fastq.gz "
    output_fastq_file = f"{fastq_dir}/ont26.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_file} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    # ***********
    # 2022 turtle
    # ***********

    # ont25
    input_fastq_files = f"/net/singapura/var/hts/ont25/25_1_flipper_tumour_ChHV5.fastq.gz /net/singapura/var/hts/ont25/25_2_neck_tumour_ChHV5.fastq.gz "
    output_fastq_file = f"{fastq_dir}/2022_turtle_ont25.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = f"/net/singapura/var/hts/ont25/25_3_pump_12_ChHV5.fastq.gz /net/singapura/var/hts/ont25/25_4_pump_13_ChHV5.fastq.gz /net/singapura/var/hts/ont25/25_5_sterivex_14_ChHV5.fastq.gz /net/singapura/var/hts/ont25/25_6_sterivex_15_ChHV5.fastq.gz"
    output_fastq_file = f"{fastq_dir}/2022_turtle_edna_ont25.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = f"/net/singapura/var/hts/ont25/25_3_pump_12_ChHV5.fastq.gz"
    output_fastq_file = f"{fastq_dir}/2022_turtle_edna_pump12_ont25.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_files} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = f"/net/singapura/var/hts/ont25/25_4_pump_13_ChHV5.fastq.gz"
    output_fastq_file = f"{fastq_dir}/2022_turtle_edna_pump13_ont25.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_files} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = f"/net/singapura/var/hts/ont25/25_5_sterivex_14_ChHV5.fastq.gz"
    output_fastq_file = f"{fastq_dir}/2022_turtle_edna_sterivex14_ont25.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_files} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = f"/net/singapura/var/hts/ont25/25_6_sterivex_15_ChHV5.fastq.gz"
    output_fastq_file = f"{fastq_dir}/2022_turtle_edna_sterivex15_ont25.fastq.gz"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"cp {input_fastq_files} {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    #########
    # mapping
    #########
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.16/bin/samtools"
    seqtk = "/usr/local/seqtk-1.3/seqtk"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

    # create bwa indices
    input_fasta_file = f"{ref_dir}/HQ878327.2.fasta"
    output_bwt_file = f"{input_fasta_file}.bwt"
    log = f"{output_bwt_file}.log"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_bwt_file}.OK"
    cmd = f"{bwa} index -a bwtsw {input_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{ref_dir}/JQ034420.1.fasta"
    output_bwt_file = f"{input_fasta_file}.bwt"
    log = f"{output_bwt_file}.log"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_bwt_file}.OK"
    cmd = f"{bwa} index -a bwtsw {input_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{ref_dir}/AB012104.1.fasta"
    output_bwt_file = f"{input_fasta_file}.bwt"
    log = f"{output_bwt_file}.log"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_bwt_file}.OK"
    cmd = f"{bwa} index -a bwtsw {input_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    # create minimap2 indices
    input_fasta_file = f"{ref_dir}/HQ878327.2.fasta"
    output_mmi_file = f"{input_fasta_file}.mmi"
    log = f"{output_mmi_file}.log"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_mmi_file}.OK"
    cmd = f"{minimap2} -d {output_mmi_file} {input_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{ref_dir}/JQ034420.1.fasta"
    output_mmi_file = f"{input_fasta_file}.mmi"
    log = f"{output_mmi_file}.log"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_mmi_file}.OK"
    cmd = f"{minimap2} -d {output_mmi_file} {input_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{ref_dir}/AB012104.1.fasta"
    output_mmi_file = f"{input_fasta_file}.mmi"
    log = f"{output_mmi_file}.log"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_mmi_file}.OK"
    cmd = f"{minimap2} -d {output_mmi_file} {input_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    # map all the relevant sequences

    # input_fastq_file1 = f'{fastq_dir}/2020_turtle_ilm31_R1.fastq.gz'
    # input_fastq_file2 = f'{fastq_dir}/2020_turtle_ilm31_R2.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2020_turtle_ont10.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2020_turtle_ont26.fastq.gz'
    # input_fastq_file1 = f'{fastq_dir}/2021_turtle_ilm35_R1.fastq.gz'
    # input_fastq_file2 = f'{fastq_dir}/2021_turtle_ilm35_R2.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2021_turtle_ont26.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/ont26.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2022_turtle_ont25.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2022_turtle_edna_ont25.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2022_turtle_edna_pump12_ont25.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2022_turtle_edna_pump13_ont25.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2022_turtle_edna_sterivex14_ont25.fastq.gz'
    # input_fastq_file = f'{fastq_dir}/2022_turtle_edna_sterivex15_ont25.fastq.gz'

    #############
    # 2020 Turtle
    #############

    # ilm31
    # HQ878327.2
    input_fastq_file1 = f"{fastq_dir}/2020_turtle_ilm31_R1.fastq.gz"
    input_fastq_file2 = f"{fastq_dir}/2020_turtle_ilm31_R2.fastq.gz"
    ref_fasta_file = f"{ref_dir}/HQ878327.2.fasta"
    output_bam_file = f"{bam_dir}/2020_turtle_ilm31.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_fasta_file}.bwt.OK {input_fastq_file1}.OK {input_fastq_file2}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {bwa} mem -t 2 -M {ref_fasta_file} {input_fastq_file1} {input_fastq_file2} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_turtle_ilm31.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file1 = f"{fastq_dir}/2020_turtle_ilm31_R1.fastq.gz"
    input_fastq_file2 = f"{fastq_dir}/2020_turtle_ilm31_R2.fastq.gz"
    ref_fasta_file = f"{ref_dir}/JQ034420.1.fasta"
    output_bam_file = f"{bam_dir}/2020_turtle_ilm31.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_fasta_file}.bwt.OK {input_fastq_file1}.OK {input_fastq_file2}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {bwa} mem -t 2 -M {ref_fasta_file} {input_fastq_file1} {input_fastq_file2} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_turtle_ilm31.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ont10
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2020_turtle_ont10.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2020_turtle_ont10.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2}  -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_turtle_ont10.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2020_turtle_ont10.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2020_turtle_ont10.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2}  -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_turtle_ont10.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ont26
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2020_turtle_ont26.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2020_turtle_ont26.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2}  -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_turtle_ont26.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2020_turtle_ont26.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2020_turtle_ont26.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2}  -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_turtle_ont26.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    #############
    # 2021 Turtle
    #############
    # ILM35
    # HQ878327.2
    input_fastq_file1 = f"{fastq_dir}/2021_turtle_ilm35_R1.fastq.gz"
    input_fastq_file2 = f"{fastq_dir}/2021_turtle_ilm35_R2.fastq.gz"
    ref_fasta_file = f"{ref_dir}/HQ878327.2.fasta"
    output_bam_file = f"{bam_dir}/2021_turtle_ilm35.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_fasta_file}.bwt.OK {input_fastq_file1}.OK {input_fastq_file2}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {bwa} mem -t 2 -M {ref_fasta_file} {input_fastq_file1} {input_fastq_file2} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2021_turtle_ilm35.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file1 = f"{fastq_dir}/2021_turtle_ilm35_R1.fastq.gz"
    input_fastq_file2 = f"{fastq_dir}/2021_turtle_ilm35_R2.fastq.gz"
    ref_fasta_file = f"{ref_dir}/JQ034420.1.fasta"
    output_bam_file = f"{bam_dir}/2021_turtle_ilm35.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_fasta_file}.bwt.OK {input_fastq_file1}.OK {input_fastq_file2}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {bwa} mem -t 2 -M {ref_fasta_file} {input_fastq_file1} {input_fastq_file2} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2021_turtle_ilm35.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ONT26
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2021_turtle_ont26.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2021_turtle_ont26.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2}  -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2021_turtle_ont26.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2021_turtle_ont26.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2021_turtle_ont26.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2021_turtle_ont26.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # 2020_2021 Turtle
    # ONT26
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/ont26.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2020_2021_turtle_ont26.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_2021_turtle_ont26.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/ont26.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2020_2021_turtle_ont26.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2020_2021_turtle_ont26.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    #############
    # 2022 Turtle
    #############
    # ONT25
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2022_turtle_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_ont25.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_ont25.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2022_turtle_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_ont25.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_ont25.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ONT25 EDNA
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_ont25.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_ont25.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_ont25.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_ont25.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ONT25 EDNA pump 12
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_pump12_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_pump12_ont25.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_pump12_ont25.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_pump12_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_pump12_ont25.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_pump12_ont25.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ONT25 EDNA pump 13
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_pump13_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_pump13_ont25.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_pump13_ont25.ChHV5.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_pump13_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_pump13_ont25.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_pump13_ont25.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ONT25 EDNA sterivex 14
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_sterivex14_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_sterivex14_ont25.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = (
        f"{stats_dir}/2022_turtle_edna_sterivex14_ont25.ChHV5.coverage.txt"
    )
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_sterivex14_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_sterivex14_ont25.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_sterivex14_ont25.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # ONT25 EDNA sterivex 15
    # HQ878327.2
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_sterivex15_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/HQ878327.2.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_sterivex15_ont25.ChHV5.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = (
        f"{stats_dir}/2022_turtle_edna_sterivex15_ont25.ChHV5.coverage.txt"
    )
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    # JQ034420.1
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_sterivex15_ont25.fastq.gz"
    ref_mmi_file = f"{ref_dir}/JQ034420.1.fasta.mmi"
    output_bam_file = f"{bam_dir}/2022_turtle_edna_sterivex15_ont25.mito.bam"
    log = f"{output_bam_file}.log"
    dep = f"{ref_mmi_file}.OK {input_fastq_file}.OK"
    tgt = f"{output_bam_file}.OK"
    cmd = f"set -o pipefail; {minimap2} -ax map-ont {ref_mmi_file} {input_fastq_file} 2> {log} | {samtools} sort - 2>> {log} | {samtools} view -o {output_bam_file}"
    pg.add(tgt, dep, cmd)

    input_bam_file = output_bam_file
    dep = f"{input_bam_file}.OK"
    tgt = f"{input_bam_file}.bai.OK"
    cmd = f"{samtools} index {input_bam_file}"
    pg.add(tgt, dep, cmd)

    output_txt_file = f"{stats_dir}/2022_turtle_edna_sterivex15_ont25.mito.coverage.txt"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
    pg.add(tgt, dep, cmd)

    ##########################
    # get consensus assemblies
    ##########################
    # have to install merdaka via pip3 install medaka
    mini_align = "/home/atks/.miniconda3/bin/mini_align"
    medaka = "/home/atks/.miniconda3/bin/medaka"

    # create directories in destination folder directory
    try:
        os.makedirs(f"{assembly_dir}/2020_turtle_ilm31_ChHV5", exist_ok=True)
        os.makedirs(f"{assembly_dir}/2020_turtle_ilm31_mito", exist_ok=True)
        os.makedirs(f"{assembly_dir}/2021_turtle_ilm35_ChHV5", exist_ok=True)
        os.makedirs(f"{assembly_dir}/2021_turtle_ilm35_mito", exist_ok=True)
        os.makedirs(f"{assembly_dir}/2022_turtle_ont25_ChHV5", exist_ok=True)
        os.makedirs(f"{assembly_dir}/2022_turtle_ont25_mito", exist_ok=True)
        os.makedirs(f"{assembly_dir}/2022_turtle_edna_ont25_ChHV5", exist_ok=True)
        os.makedirs(f"{assembly_dir}/2022_turtle_edna_ont25_mito", exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    #############
    # 2020 Turtle
    #############
    consensus_dir = f"{assembly_dir}/2020_turtle_ilm31_ChHV5"
    input_bam_file = f"{bam_dir}/2020_turtle_ilm31.ChHV5.bam"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2020_turtle_ilm31_ChHV5.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2020_Turtle_ChHV5 {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    consensus_dir = f"{assembly_dir}/2020_turtle_ilm31_mito"
    input_bam_file = f"{bam_dir}/2020_turtle_ilm31.mito.bam"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2020_turtle_ilm31_mito.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2020_Turtle_mito {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    #############
    # 2021 Turtle
    #############
    consensus_dir = f"{assembly_dir}/2021_turtle_ilm35_ChHV5"
    input_bam_file = f"{bam_dir}/2021_turtle_ilm35.ChHV5.bam"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2021_turtle_ilm35_ChHV5.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2021_Turtle_ChHV5 {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    consensus_dir = f"{assembly_dir}/2021_turtle_ilm35_mito"
    input_bam_file = f"{bam_dir}/2021_turtle_ilm35.mito.bam"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{samtools} consensus {input_bam_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2021_turtle_ilm35_mito.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2021_Turtle_mito {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    #############
    # 2022 Turtle
    #############
    # tissue
    # ONT25 ChHV5
    consensus_dir = f"{assembly_dir}/2022_turtle_ont25_ChHV5"
    ref_fasta_file = f"{ref_dir}/HQ878327.2.fasta"
    input_fastq_file = f"{fastq_dir}/2022_turtle_ont25.fastq.gz"

    bam_file_prefix = f"{consensus_dir}/reads"
    log = f"{consensus_dir}/align.log"
    dep = f"{input_fastq_file}.OK"
    tgt = f"{consensus_dir}/aligned.OK"
    cmd = f"{mini_align} -i {input_fastq_file} -r {ref_fasta_file} -m -p {bam_file_prefix} -t 2 2> {log}"
    pg.add(tgt, dep, cmd)

    input_bam_file = f"{consensus_dir}/reads.bam"
    output_hdf_file = f"{consensus_dir}/reads.hdf"
    log = f"{consensus_dir}/consensus.log"
    dep = f"{consensus_dir}/aligned.OK"
    tgt = f"{consensus_dir}/consensus.OK"
    cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model r941_min_hac_g507 2> {log}"
    pg.add(tgt, dep, cmd)

    input_hdf_file = f"{consensus_dir}/reads.hdf"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    log = f"{consensus_dir}/stitched.log"
    dep = f"{consensus_dir}/consensus.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{medaka} stitch {input_hdf_file} {ref_fasta_file} {output_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2022_turtle_ont25_ChHV5.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2022_Turtle_ChHV5 {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    # ONT25 mito
    consensus_dir = f"{assembly_dir}/2022_turtle_ont25_mito"
    ref_fasta_file = f"{ref_dir}/JQ034420.1.fasta"
    input_fastq_file = f"{fastq_dir}/2022_turtle_ont25.fastq.gz"

    bam_file_prefix = f"{consensus_dir}/reads"
    log = f"{consensus_dir}/align.log"
    dep = f"{input_fastq_file}.OK"
    tgt = f"{consensus_dir}/aligned.OK"
    cmd = f"{mini_align} -i {input_fastq_file} -r {ref_fasta_file} -m -p {bam_file_prefix} -t 2 2> {log}"
    pg.add(tgt, dep, cmd)

    input_bam_file = f"{consensus_dir}/reads.bam"
    output_hdf_file = f"{consensus_dir}/reads.hdf"
    log = f"{consensus_dir}/consensus.log"
    dep = f"{consensus_dir}/aligned.OK"
    tgt = f"{consensus_dir}/consensus.OK"
    cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model r941_min_hac_g507 2> {log}"
    pg.add(tgt, dep, cmd)

    input_hdf_file = f"{consensus_dir}/reads.hdf"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    log = f"{consensus_dir}/stitched.log"
    dep = f"{consensus_dir}/consensus.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{medaka} stitch {input_hdf_file} {ref_fasta_file} {output_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2022_turtle_ont25_mito.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2022_Turtle_mito {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    # edna
    # ONT25 ChHV5
    consensus_dir = f"{assembly_dir}/2022_turtle_edna_ont25_ChHV5"
    ref_fasta_file = f"{ref_dir}/HQ878327.2.fasta"
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_ont25.fastq.gz"

    bam_file_prefix = f"{consensus_dir}/reads"
    log = f"{consensus_dir}/align.log"
    dep = f"{input_fastq_file}.OK"
    tgt = f"{consensus_dir}/aligned.OK"
    cmd = f"{mini_align} -i {input_fastq_file} -r {ref_fasta_file} -m -p {bam_file_prefix} -t 2 2> {log}"
    pg.add(tgt, dep, cmd)

    input_bam_file = f"{consensus_dir}/reads.bam"
    output_hdf_file = f"{consensus_dir}/reads.hdf"
    log = f"{consensus_dir}/consensus.log"
    dep = f"{consensus_dir}/aligned.OK"
    tgt = f"{consensus_dir}/consensus.OK"
    cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model r941_min_hac_g507 2> {log}"
    pg.add(tgt, dep, cmd)

    input_hdf_file = f"{consensus_dir}/reads.hdf"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    log = f"{consensus_dir}/stitched.log"
    dep = f"{consensus_dir}/consensus.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{medaka} stitch --fill_char N {input_hdf_file} {ref_fasta_file} {output_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2022_turtle_edna_ont25_ChHV5.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2022_Turtle_eDNA_ChHV5 {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    # ONT25 mito
    consensus_dir = f"{assembly_dir}/2022_turtle_edna_ont25_mito"
    ref_fasta_file = f"{ref_dir}/JQ034420.1.fasta"
    input_fastq_file = f"{fastq_dir}/2022_turtle_edna_ont25.fastq.gz"

    bam_file_prefix = f"{consensus_dir}/reads"
    log = f"{consensus_dir}/align.log"
    dep = f"{input_fastq_file}.OK"
    tgt = f"{consensus_dir}/aligned.OK"
    cmd = f"{mini_align} -i {input_fastq_file} -r {ref_fasta_file} -m -p {bam_file_prefix} -t 2 2> {log}"
    pg.add(tgt, dep, cmd)

    input_bam_file = f"{consensus_dir}/reads.bam"
    output_hdf_file = f"{consensus_dir}/reads.hdf"
    log = f"{consensus_dir}/consensus.log"
    dep = f"{consensus_dir}/aligned.OK"
    tgt = f"{consensus_dir}/consensus.OK"
    cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model r941_min_hac_g507 2> {log}"
    pg.add(tgt, dep, cmd)

    input_hdf_file = f"{consensus_dir}/reads.hdf"
    output_fasta_file = f"{consensus_dir}/consensus.fasta"
    log = f"{consensus_dir}/stitched.log"
    dep = f"{consensus_dir}/consensus.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{medaka} stitch --fill_char N {input_hdf_file} {ref_fasta_file} {output_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/consensus.fasta"
    output_fasta_file = f"{fasta_dir}/2022_turtle_edna_ont25_mito.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r 2022_Turtle_eDNA_mito {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    # clean
    pg.add_clean(
        f"rm -fr {ref_dir} {fastq_dir} {fasta_dir} {bam_dir} {assembly_dir} {stats_dir}"
    )

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


if __name__ == "__main__":
    main()
