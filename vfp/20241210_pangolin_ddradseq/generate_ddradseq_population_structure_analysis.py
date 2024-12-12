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


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_ddradseq_population_structure_analysis.mk",
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
@click.option("-g", "--genome_fasta_file", default="", required=False, help="genome FASTA file")
def main(make_file, working_dir, sample_file, genome_fasta_file):
    """
    Population structure of Pangolins

    e.g. generate_ddradseq_population_structure_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("genome fasta file", genome_fasta_file))

    # read sample file
    samples = {}
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                id, fastq1, fastq2 = line.rstrip().split("\t")
                if id not in samples:
                    samples[id] = Sample(id, fastq1, fastq2)
                else:
                    samples[id].add_fastq(fastq1, fastq2)

    # for id in samples.keys():
    #     samples[id].print()

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    bam_dir = f"{working_dir}/bam"
    fastq_dir = f"{working_dir}/fastq"
    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    #################
    # reference files
    #################
    # reference genome assembly for sunda pangolin
    # https://www.ncbi.nlm.nih.gov/assembly/GCF_014570535.1

    id = "AJ810453.1"
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

    #################
    # variant calling
    #################

    #################
    # variant QC
    # * HWE stats
    # * PCA
    # * * combine with other data sets
    # * *
    #################

    #########
    # mapping
    #########
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.16/bin/samtools"
    seqtk = "/usr/local/seqtk-1.3/seqtk"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"

    for id, sample in samples.items():
        #combine files
        if len(sample.fastq1s) == 1:
            input_fastq1_file = sample.fastq1s[0]
            output_fastq1_file = f"{fastq_dir}/{sample.id}.R1.fastq.gz"
            tgt = f"{output_fastq1_file}.OK"
            cmd = f"ln -s {input_fastq1_file} {output_fastq1_file}"
            pg.add(tgt, dep, cmd)

            input_fastq2_file = sample.fastq2s[0]
            output_fastq2_file = f"{fastq_dir}/{sample.id}.R2.fastq.gz"
            tgt = f"{output_fastq2_file}.OK"
            cmd = f"ln -s {input_fastq2_file} {output_fastq2_file}"
            pg.add(tgt, dep, cmd)

        elif len(sample.fastq1s) > 1:
            input_fastq1_files = " ".join(sample.fastq1s)
            output_fastq1_file = f"{fastq_dir}/{sample.id}.R1.fastq.gz"
            tgt = f"{output_fastq1_file}.OK"
            cmd = f"zcat {input_fastq1_files} | gzip > {output_fastq1_file}"
            pg.add(tgt, dep, cmd)

            input_fastq2_files = " ".join(sample.fastq2s)
            output_fastq2_file = f"{fastq_dir}/{sample.id}.R2.fastq.gz"
            tgt = f"{output_fastq2_file}.OK"
            cmd = f"zcat {input_fastq2_files} | gzip > {output_fastq2_file}"
            pg.add(tgt, dep, cmd)

    # clean
    pg.add_clean(f"rm -fr {ref_dir} {bam_dir} {fastq_dir}")

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

    def __init__(self, id, fastq1, fastq2):
        self.id = id
        self.fastq1s = []
        self.fastq1s.append(fastq1)
        self.fastq2s = []
        self.fastq2s.append(fastq2)

    def add_fastq(self, fastq1, fastq2):
        self.fastq1s.append(fastq1)
        self.fastq2s.append(fastq2)

    def print(self):
        print(f"id       : {self.id}")
        print(f"no files : {len(self.fastq1s)}")
        print(f"fastq1s  : {','.join(self.fastq1s)}")
        print(f"fastq2s  : {','.join(self.fastq2s)}")


if __name__ == "__main__":
    main() # type: ignore
