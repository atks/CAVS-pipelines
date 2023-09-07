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


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="amr_resfinder_plasmid_finder.mk",
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
    Runs resFinder and plasmidFinder on isolate sequences

    e.g. vm_generate_amr_analysis_pipeline.py -r ilm1 -i raw -i ilm23.sa
    """
    log_dir = f"{working_dir}/log"

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, species, fastq1, fastq2, fasta = line.rstrip().split("\t")
                renamed_species = ""
                if species == "escherichia_coli":
                    renamed_species = "Escherichia coli"
                elif species in [
                    "salmonella_typhimurium",
                    "salmonella_enteritidis",
                    "salmonella",
                ]:
                    renamed_species = "Salmonella enterica"
                samples.append(
                    Sample(index, sample_id, renamed_species, fastq1, fastq2, fasta)
                )
                # print(f'{species} => {renamed_species}')
    # create directories in destination folder directory
    analysis_dir = f"{working_dir}/analysis"
    new_dir = ""
    try:
        os.makedirs(log_dir, exist_ok=True)

        for sample in samples:
            # new_dir = f'{working_dir}/resfinder/reads/{sample.id}'
            # os.makedirs(new_dir, exist_ok=True)
            # new_dir = f'{working_dir}/resfinder/contigs/{sample.id}'
            # os.makedirs(new_dir, exist_ok=True)
            # new_dir = f'{working_dir}/pointfinder/reads/{sample.id}'
            # os.makedirs(new_dir, exist_ok=True)
            # new_dir = f'{working_dir}/pointfinder/contigs/{sample.id}'
            # os.makedirs(new_dir, exist_ok=True)
            new_dir = f"{working_dir}/results/{sample.id}"
            os.makedirs(new_dir, exist_ok=True)

    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # analyze
    for idx, sample in enumerate(samples):
        # /usr/local/resfinder-4.1.5/run_resfinder.py
        #    -ifq salmonella_R1.fastq.gz salmonella_R2.fastq.gz
        #    -o salmonella  -s 'Salmonella' --min_cov 0.6 -t 0.8 --acquired
        #    -k /usr/local/kma-1.4.7/kma
        # /usr/local/resfinder-4.1.5/run_resfinder.py
        #    -ifq salmonella_R1.fastq.gz salmonella_R2.fastq.gz
        #    -o salmonella  -s 'Salmonella' --min_cov 0.6 -t 0.8 --acquired
        #    -k /usr/local/kma-1.4.7/kma
        # resfinder
        resfinder = "/usr/local/resfinder-4.1.5/run_resfinder.py"
        db_resfinder = "/usr/local/resfinder-4.1.5/db_resfinder"
        db_pointfinder = "/usr/local/resfinder-4.1.5/db_pointfinder"
        kma = "/usr/local/kma-1.4.7/kma"
        blastn = "/usr/local/ncbi-blast-2.13.0+/bin/blastn"
        species = sample.species
        fastqs = f"{sample.fastq1} {sample.fastq2}"
        fasta = sample.fasta

        # resfinder, pointfinder, contigs
        output_dir = f"{working_dir}/results/{sample.id}"
        log = f"{output_dir}/run.log"
        err = f"{output_dir}/run.err"
        tgt = f"{log_dir}/{sample.id}_resfinder_pointfinder_contigs.OK"
        dep = ""
        cmd = f'{resfinder} -o {output_dir} -ifa {fasta} -b {blastn} --min_cov 0.6 -t 0.8 --acquired -s "{species}" --point > {log} 2> {err}'
        pg.add(tgt, dep, cmd)

        # #resfinder, reads
        # output_dir = f'{working_dir}/resfinder/reads/{sample.id}'
        # log = f'{output_dir}/run.log'
        # err = f'{output_dir}/run.err'
        # tgt = f'{log_dir}/{sample.id}_resfinder_reads.OK'
        # dep = ''
        # cmd = f'{resfinder} -o {output_dir} -ifq {fastqs} -k {kma} --min_cov 0.6 -t 0.8 --acquired > {log} 2> {err}'
        # pg.add(tgt, dep, cmd)

        # #pointfinder, contigs
        # output_dir = f'{working_dir}/resfinder/contigs/{sample.id}'
        # log = f'{output_dir}/run.log'
        # err = f'{output_dir}/run.err'
        # tgt = f'{log_dir}/{sample.id}_resfinder_contigs.OK'
        # dep = ''
        # cmd = f'{resfinder} -o {output_dir} -ifa {fasta} -b {blastn} --min_cov 0.6 -t 0.8  --acquired > {log} 2> {err}'
        # pg.add(tgt, dep, cmd)

        # #pointfinder, reads
        # output_dir = f'{working_dir}/pointfinder/reads/{sample.id}'
        # log = f'{output_dir}/run.log'
        # err = f'{output_dir}/run.err'
        # tgt = f'{log_dir}/{sample.id}_pointfinder_reads.OK'
        # dep = ''
        # cmd = f'{resfinder} -o {output_dir} -ifq {fastqs} -k {kma} -s {species} --point > {log} 2> {err}'
        # pg.add(tgt, dep, cmd)

        # #pointfinder, contigs
        # output_dir = f'{working_dir}/pointfinder/contigs/{sample.id}'
        # log = f'{output_dir}/run.log'
        # err = f'{output_dir}/run.err'
        # tgt = f'{log_dir}/{sample.id}_pointfinder_contigs.OK'
        # dep = ''
        # cmd = f'{resfinder} -o {output_dir} -ifa {fasta} -b {blastn} -s {species} --point > {log} 2> {err}'
        # pg.add(tgt, dep, cmd)

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
        self.fasta = ""

    def __init__(self, idx, id, species, fastq1, fastq2, fasta):
        self.idx = idx
        self.id = id
        self.species = species
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.fasta = fasta

    def print(self):
        print(f"index   : {self.idx}")
        print(f"id      : {self.id}")
        print(f"fastq1  : {self.fastq1}")
        print(f"fastq2 : {self.fastq2}")


if __name__ == "__main__":
    main()
