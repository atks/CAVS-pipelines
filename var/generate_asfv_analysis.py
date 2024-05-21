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
    default="run_asfv_analysis.mk",
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
    African Swine Fever Virus Assembly

    e.g. generate_asfv_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    fastq_dir = f"{working_dir}/fastq"
    fasta_dir = f"{working_dir}/fasta"
    bam_dir = f"{working_dir}/bam"
    consensus_dir = f"{working_dir}/consensus"
    stats_dir = f"{working_dir}/stats"
    coverage_stats_dir = f"{working_dir}/stats/coverage"

    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(consensus_dir, exist_ok=True)
        os.makedirs(f"{consensus_dir}/ilm", exist_ok=True)
        os.makedirs(f"{consensus_dir}/ont", exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(coverage_stats_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    #################
    # reference files
    #################
    id = "FR682468.2"
    output_fasta_file = f"{ref_dir}/{id}.fasta"
    dep = ""
    tgt = f"{output_fasta_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    id = "FR682468.2"
    output_genbank_file = f"{ref_dir}/{id}.genbank"
    dep = ""
    tgt = f"{output_genbank_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format genbank > {output_genbank_file}"
    pg.add(tgt, dep, cmd)

    ###############
    # sequence file
    ###############
    fastq_files = [
        "M230221_ASFV_liver",
        "M230221_ASFV_lung",
        "M230221_ASFV_heart",
        "M230221_ASFV_spleen",
        "M230221_ASFV_lymph_node",
        "M230221_ASFV_kidney",
        "M230221_ASFV_tonsil",
        "M230221_ASFV_abdominal_fluid",
        "M230221_ASFV_thoracic_fluid",
        "unclassified",
    ]

    all_ilm_fastq1_files = ""
    all_ilm_fastq2_files = ""
    all_ont_fastq_files = ""

    samples = []
    for i, root_name in enumerate(fastq_files):
        id = i + 1
        tech = "ilm"
        if root_name == "unclassified":
            name = f"41_{id+1}_{root_name}"
            fastq1 = f"/net/singapura/var/hts/ilm41/41_{id+1}_{root_name}_R1.fastq.gz"
            fastq2 = f"/net/singapura/var/hts/ilm41/41_{id+1}_{root_name}_R2.fastq.gz"
        else:
            name = f"41_{id}_{root_name}"
            fastq1 = f"/net/singapura/var/hts/ilm41/41_{id}_{root_name}_R1.fastq.gz"
            fastq2 = f"/net/singapura/var/hts/ilm41/41_{id}_{root_name}_R2.fastq.gz"
        samples.append(Sample(id, name, tech, fastq1, fastq2, "", ""))

        all_ilm_fastq1_files += f" {fastq1}"
        all_ilm_fastq2_files += f" {fastq2}"

        tech = "ont"
        name = f"27-28_{id}_{root_name}"
        fastq1 = f"{fastq_dir}/27-28_{id}_{root_name}.fastq.gz"
        fastq2 = f""
        samples.append(Sample(id, name, tech, fastq1, fastq2, "", ""))

        all_ont_fastq_files += f" {fastq1}"

        # merge fastq files for ont
        output_fastq_file = f"{fastq_dir}/27-28_{id}_{root_name}.fastq.gz"
        fastq_files = f"/net/singapura/var/hts/ont27/27_{id}_{root_name}.fastq.gz /net/singapura/var/hts/ont28/28_{id}_{root_name}.fastq.gz"
        dep = ""
        tgt = f"{output_fastq_file}.OK"
        cmd = f"zcat {fastq_files} | gzip > {output_fastq_file}"
        pg.add(tgt, dep, cmd)

    # generate aggregate fastq files
    output_fastq_file = f"{fastq_dir}/ilm41.all_R1.fastq.gz"
    fastq1 = output_fastq_file
    fastq1_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {all_ilm_fastq1_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    output_fastq_file = f"{fastq_dir}/ilm41.all_R2.fastq.gz"
    fastq2 = output_fastq_file
    fastq2_OK = f"{fastq2}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {all_ilm_fastq2_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    samples.append(Sample(11, "ilm_all", "ilm", fastq1, fastq2, fastq1_OK, fastq2_OK))

    output_fastq_file = f"{fastq_dir}/ont27-28.all.fastq.gz"
    fastq1 = output_fastq_file
    fastq1_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {all_ont_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    samples.append(Sample(11, "ont_all", "ont", fastq1, "", fastq1_OK, ""))

    #################
    # align sequences
    #################
    # mapping with minimap2
    minimap2 = "/usr/local/minimap2-2.24/minimap2"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.15/samtools"

    # bwa index -a bwtsw NC_039223.1.fa
    # minimap2 -d HQ878327.mmi HQ878327.fasta
    bwa_ref_fasta_file = f"{ref_dir}/FR682468.2.fasta"
    minimap2_ref_mmi_file = f"{ref_dir}/FR682468.2.mmi"

    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"
    ref_fasta_file = f"{ref_dir}/FR682468.2.fasta"

    # bwa index
    log = f"{ref_dir}/bwa_index.log"
    err = f"{ref_dir}/bwa_index.err"
    dep = f"{ref_fasta_file}.OK"
    tgt = f"{ref_dir}/bwa_index.OK"
    cmd = f"{bwa} index -a bwtsw {ref_fasta_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # minimap2 index
    log = f"{ref_dir}/minimap2_index.log"
    err = f"{ref_dir}/minimap2_index.err"
    dep = f"{ref_fasta_file}.OK"
    tgt = f"{ref_dir}/minimap2_index.OK"
    cmd = f"{minimap2} -d {ref_dir}/FR682468.2.mmi {ref_fasta_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    for sample in samples:
        if sample.seq_tech == "ont":
            output_bam_file = f"{bam_dir}/{sample.name}.bam"
            log = f"{output_bam_file}.log"
            err = f"{output_bam_file}.err"
            dep = f"{sample.fastq1_OK} {ref_dir}/minimap2_index.OK"
            tgt = f"{output_bam_file}.OK"
            cmd = f"{minimap2} -ax map-ont {minimap2_ref_mmi_file} {sample.fastq1} 2> {err}| {samtools} sort - | {samtools} view -o {output_bam_file} > {log} "
            pg.add(tgt, dep, cmd)

            input_bam_file = f"{bam_dir}/{sample.name}.bam"
            dep = f"{input_bam_file}.OK"
            tgt = f"{input_bam_file}.bai.OK"
            cmd = f"{samtools} index {input_bam_file}"
            pg.add(tgt, dep, cmd)

            input_bam_file = f"{bam_dir}/{sample.name}.bam"
            output_txt_file = f"{coverage_stats_dir}/{sample.name}.coverage.txt"
            dep = f"{input_bam_file}.OK"
            tgt = f"{output_txt_file}.OK"
            cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
            pg.add(tgt, dep, cmd)

        else:
            output_bam_file = f"{bam_dir}/{sample.name}.bam"
            log = f"{output_bam_file}.log"
            err = f"{output_bam_file}.log"
            dep = f"{sample.fastq1_OK} {sample.fastq2_OK} {ref_dir}/minimap2_index.OK"
            tgt = f"{output_bam_file}.OK"
            cmd = f"{bwa} mem -t 2 -M {bwa_ref_fasta_file} {sample.fastq1} {sample.fastq2} 2> {err}| {samtools} sort - | {samtools} view -o {output_bam_file} > {log} 2> {err}"
            pg.add(tgt, dep, cmd)

            input_bam_file = f"{bam_dir}/{sample.name}.bam"
            dep = f"{input_bam_file}.OK"
            tgt = f"{input_bam_file}.bai.OK"
            cmd = f"{samtools} index {input_bam_file}"
            pg.add(tgt, dep, cmd)

            input_bam_file = f"{bam_dir}/{sample.name}.bam"
            output_txt_file = f"{coverage_stats_dir}/{sample.name}.coverage.txt"
            dep = f"{input_bam_file}.OK"
            tgt = f"{output_txt_file}.OK"
            cmd = f"{samtools} coverage {input_bam_file} > {output_txt_file}"
            pg.add(tgt, dep, cmd)

    ##############################
    # generate consensus sequences
    ##############################

    # FR682468.2:104591-105047

    # ILM41
    input_bam_file = f"{bam_dir}/ilm_all.bam"
    output_fasta_file = f"{consensus_dir}/ilm/asfv.ilm.fasta"
    log = f"{output_fasta_file}.log"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f'{samtools} consensus {input_bam_file} | {seqkit} replace -p "^.+$\\" -r asfv_cavs_ilm_FR682468.2 > {output_fasta_file} 2> {log}'
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/ilm/asfv.ilm.fasta"
    output_fasta_file = f"{fasta_dir}/asfv.ilm.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"cp {input_fasta_file} {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    # ONT27/28
    # tools from medaka
    medaka = "/home/atks/.miniconda3/bin/medaka"
    ref_fasta_file = f"{ref_dir}/FR682468.2.fasta"

    input_bam_file = f"{bam_dir}/ont_all.bam"
    output_hdf_file = f"{consensus_dir}/ont/reads.hdf"
    log = f"{output_hdf_file}.log"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_hdf_file}.OK"
    cmd = f"{medaka} consensus {input_bam_file} {output_hdf_file} --model r941_min_hac_g507 2> {log}"
    pg.add(tgt, dep, cmd)

    input_hdf_file = f"{consensus_dir}/ont/reads.hdf"
    output_fasta_file = f"{consensus_dir}/ont/consensus.fasta"
    log = f"{output_fasta_file}.log"
    dep = f"{input_hdf_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{medaka} stitch --fill_char N {input_hdf_file} {ref_fasta_file} {output_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/ont/consensus.fasta"
    output_fasta_file = f"{consensus_dir}/ont/asfv.ont.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"{seqkit} replace -p '^.+$\\' -r asfv_cavs_ont_FR682468.2 {input_fasta_file} > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    input_fasta_file = f"{consensus_dir}/ont/asfv.ont.fasta"
    output_fasta_file = f"{fasta_dir}/asfv.ont.fasta"
    dep = f"{input_fasta_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f"cp {input_fasta_file} {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    # comparison
    stretcher = "/usr/local/emboss-6.6.0/bin/stretcher"
    input_fasta_file1 = f"{consensus_dir}/ilm/asfv.ilm.fasta"
    input_fasta_file2 = f"{consensus_dir}/ont/asfv.ont.fasta"
    output_txt_file = f"{consensus_dir}/ilm_ont.stretcher.txt"
    dep = f"{input_fasta_file1}.OK {input_fasta_file2}.OK"
    tgt = f"{output_txt_file}.OK"
    cmd = f"{stretcher} {input_fasta_file1} {input_fasta_file2} {output_txt_file}"
    pg.add(tgt, dep, cmd)

    ###################################
    # extract sequence for phylogenetic
    ###################################

    ####################
    # compute viral load
    ####################

    ######################
    # IGR analysis
    ######################

    # clean
    pg.add_clean(f"rm -fr {ref_dir} {fastq_dir} {fasta_dir} {bam_dir} {stats_dir} ")

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
                if "|" in self.cmds[i]:
                    f.write(f"\tset -o pipefail; {self.cmds[i]}\n")
                else:
                    f.write(f"\t{self.cmds[i]}\n")
                f.write(f"\ttouch {self.tgts[i]}\n\n")

            if self.clean_cmd != "":
                f.write(f"clean : \n")
                f.write(f"\t{self.clean_cmd}\n")


class Sample(object):
    def __init__(self):
        self.id = ""
        self.name = ""
        self.seq_tech = ""
        self.fastq1 = ""
        self.fastq2 = ""
        self.fastq1_OK = ""
        self.fastq2_OK = ""

    def __init__(self, id, name, seq_tech, fastq1, fastq2):
        self.id = id
        self.name = name
        self.seq_tech = seq_tech
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.fastq1_OK = ""
        self.fastq2_OK = ""

    def __init__(self, id, name, seq_tech, fastq1, fastq2, fastq1_OK, fastq2_OK):
        self.id = id
        self.name = name
        self.seq_tech = seq_tech
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.fastq1_OK = fastq1_OK
        self.fastq2_OK = fastq2_OK

    def print(self):
        print(f"id         : {self.id}")
        print(f"name       : {self.name}")
        print(f"seq tech   : {self.seq_tech}")
        print(f"fastq1     : {self.fastq1}")
        print(f"fastq2     : {self.fastq2}")
        print(f"fastq1_OK  : {self.fastq1_OK}")
        print(f"fastq2_OK  : {self.fastq2_OK}")


if __name__ == "__main__":
    main()
