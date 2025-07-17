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
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import click
import python-docx

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
@click.option("-p", "--population_map_file", required=True, help="population map file")
@click.option("-g", "--genome_fasta_file", default="", required=False, help="genome FASTA file")
def main(make_file, working_dir, sample_file, population_map_file, genome_fasta_file):
    """
    Population structure of Pangolins

    e.g. generate_ddradseq_population_structure_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("population map file", population_map_file))
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

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    log_dir = f"{working_dir}/log"
    stats_dir = f"{working_dir}/stats"
    plot_dir = f"{working_dir}/plot"
    fastq_dir = f"{working_dir}/fastq"
    fastq_ulen_dir = f"{working_dir}/fastq_ulen"
    bam_dir = f"{working_dir}/bam"
    vcf_dir = f"{working_dir}/vcf"
    annotation_dir = f"{working_dir}/annotations"
    denovo_stacks_dir = f"{working_dir}/denovo_stacks"
    ref_stacks_dir = f"{working_dir}/ref_stacks"
    qc_dir = f"{working_dir}/qc"
    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(f"{stats_dir}/coverage", exist_ok=True)
        os.makedirs(f"{stats_dir}/general", exist_ok=True)
        os.makedirs(f"{stats_dir}/flag", exist_ok=True)
        os.makedirs(f"{stats_dir}/idx", exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)
        os.makedirs(fastq_dir, exist_ok=True)
        os.makedirs(fastq_ulen_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(vcf_dir, exist_ok=True)
        os.makedirs(denovo_stacks_dir, exist_ok=True)
        os.makedirs(ref_stacks_dir, exist_ok=True)
        os.makedirs(qc_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    
    

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
