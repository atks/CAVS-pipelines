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
    default="run_asfv_bulan_analysis.mk",
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

    e.g. generate_asfv_bulan_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    # initialize
    pg = PipelineGenerator(make_file)

    # programs
    efetch = "/usr/local/edirect-17.0/efetch"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    bwa = "/usr/local/bwa-0.7.17/bwa"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

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
    cmd = f"{efetch} -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    ###############
    # sequence file
    ###############
    fastq_files = [
        "M230446_ASFV_pig1_kidney",
        "M230446_ASFV_pig1_lymph1",
        "M230446_ASFV_pig1_lymph2",
        "M230446_ASFV_pig2_spleen",
        "M230446_ASFV_pig2_lymph1",
        "M230446_ASFV_pig2_lymph2",
        "unclassified",
    ]

    pig1_fastq1_files = ""
    pig1_fastq2_files = ""
    pig2_fastq1_files = ""
    pig2_fastq2_files = ""

    samples = []
    for i, root_name in enumerate(fastq_files):
        id = i + 1
        tech = "ilm"
        name = f"50_{id}_{root_name}"
        fastq1 = f"/net/singapura/var/hts/ilm50/50_{id}_{root_name}_R1.fastq.gz"
        fastq2 = f"/net/singapura/var/hts/ilm50/50_{id}_{root_name}_R2.fastq.gz"
        samples.append(Sample(id, name, tech, fastq1, fastq2, "", ""))

        if "pig1" in root_name:
            pig1_fastq1_files += f" {fastq1}"
            pig1_fastq2_files += f" {fastq2}"
        elif "pig2" in root_name:
            pig2_fastq1_files += f" {fastq1}"
            pig2_fastq2_files += f" {fastq2}"

    # generate aggregate fastq files for pig1 and pig2 and combo pig
    input_fastq_files = pig1_fastq1_files
    output_fastq_file = f"{fastq_dir}/pig1_R1.fastq.gz"
    fastq1 = output_fastq_file
    fastq1_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = pig1_fastq2_files
    output_fastq_file = f"{fastq_dir}/pig1_R2.fastq.gz"
    fastq2 = output_fastq_file
    fastq2_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    samples.append(Sample(8, "pig1", "ilm", fastq1, fastq2, fastq1_OK, fastq2_OK))

    # generate aggregate fastq files for pig2
    input_fastq_files = pig2_fastq1_files
    output_fastq_file = f"{fastq_dir}/pig2_R1.fastq.gz"
    fastq1 = output_fastq_file
    fastq1_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = pig2_fastq2_files
    output_fastq_file = f"{fastq_dir}/pig2_R2.fastq.gz"
    fastq2 = output_fastq_file
    fastq2_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    samples.append(Sample(9, "pig2", "ilm", fastq1, fastq2, fastq1_OK, fastq2_OK))

    # generate aggregate fastq files for combo pig
    input_fastq_files = f"{pig1_fastq1_files} {pig2_fastq1_files}"
    output_fastq_file = f"{fastq_dir}/combo_pig_R1.fastq.gz"
    fastq1 = output_fastq_file
    fastq1_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    input_fastq_files = f"{pig1_fastq2_files} {pig2_fastq2_files}"
    output_fastq_file = f"{fastq_dir}/combo_pig_R2.fastq.gz"
    fastq2 = output_fastq_file
    fastq2_OK = f"{fastq1}.OK"
    dep = ""
    tgt = f"{output_fastq_file}.OK"
    cmd = f"zcat {input_fastq_files} | gzip > {output_fastq_file}"
    pg.add(tgt, dep, cmd)

    samples.append(Sample(10, "combo_pig", "ilm", fastq1, fastq2, fastq1_OK, fastq2_OK))

    #################
    # align sequences
    #################
    # bwa index -a bwtsw NC_039223.1.fa
    bwa_ref_fasta_file = f"{ref_dir}/FR682468.2.fasta"
    ref_fasta_file = f"{ref_dir}/FR682468.2.fasta"

    # bwa index
    log = f"{ref_dir}/bwa_index.log"
    err = f"{ref_dir}/bwa_index.err"
    dep = f"{ref_fasta_file}.OK"
    tgt = f"{ref_dir}/bwa_index.OK"
    cmd = f"{bwa} index -a bwtsw {ref_fasta_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    for sample in samples:
        output_bam_file = f"{bam_dir}/{sample.name}.bam"
        log = f"{output_bam_file}.log"
        err = f"{output_bam_file}.log"
        dep = f"{sample.fastq1_OK} {sample.fastq2_OK} {ref_dir}/bwa_index.OK"
        tgt = f"{output_bam_file}.OK"
        cmd = f"{bwa} mem -t 2 -M {bwa_ref_fasta_file} {sample.fastq1} {sample.fastq2} 2> {err}| {samtools} view -h -F 4 | {samtools} sort - | {samtools} view -o {output_bam_file} > {log} 2> {err}"
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

    #################################################
    # generate consensus for pig1, pig2 and combo pig
    #################################################
    input_bam_file = f"{bam_dir}/pig1.bam"
    output_fasta_file = f"{consensus_dir}/pig1.asfv.fasta"
    log = f"{output_fasta_file}.log"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f'{samtools} consensus {input_bam_file} | {seqkit} replace -p "^.+$\\" -r CAVS_ASFV_PIG1 > {output_fasta_file} 2> {log}'
    pg.add(tgt, dep, cmd)

    input_bam_file = f"{bam_dir}/pig2.bam"
    output_fasta_file = f"{consensus_dir}/pig2.asfv.fasta"
    log = f"{output_fasta_file}.log"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f'{samtools} consensus {input_bam_file} | {seqkit} replace -p "^.+$\\" -r CAVS_ASFV_PIG2 > {output_fasta_file} 2> {log}'
    pg.add(tgt, dep, cmd)

    input_bam_file = f"{bam_dir}/combo_pig.bam"
    output_fasta_file = f"{consensus_dir}/combo_pig.asfv.fasta"
    log = f"{output_fasta_file}.log"
    dep = f"{input_bam_file}.OK"
    tgt = f"{output_fasta_file}.OK"
    cmd = f'{samtools} consensus {input_bam_file} | {seqkit} replace -p "^.+$\\" -r CAVS_ASFV_COMBO_PIG > {output_fasta_file} 2> {log}'
    pg.add(tgt, dep, cmd)

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
