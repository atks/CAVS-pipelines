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
from shutil import copy2, which

@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_vm_salmonella_typing_pipeline.mk",
    help="make file name",
)
@click.option(
    "-o",
    "--output_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(make_file, output_dir, sample_file):
    """
    Analyse run raw reads and assembled contigs on
        a. seroseq2 v1.3.1
        b. SISTR 1.1.3
        c. mlst 2.23.0

    e.g. generate_vm_salmonella_typing_pipeline.py -s salmonella_pt.sa
    """
    if not os.path.isabs(output_dir):
        cur_dir =  os.getcwd()
        output_dir = f"{cur_dir}/{output_dir}"

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("output_dir", output_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, fastq1, fastq2, contigs_fasta = line.rstrip().split("\t")
                samples.append(
                    Sample(index, sample_id, fastq1, fastq2, contigs_fasta)
                )

    # create directories in destination folder directory
    trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(trace_dir, exist_ok=True)
        for sample in samples:
            os.makedirs(f"{output_dir}/{sample.id}/seqsero2", exist_ok=True)
            os.makedirs(f"{output_dir}/{sample.id}/mlst", exist_ok=True)
            os.makedirs(f"{output_dir}/{sample.id}/sistr", exist_ok=True)

    except OSError as error:
        print(f"{error.filename} cannot be created")

    #version
    version = "1.3.0"

    # programs
    seqsero2 = "/usr/local/SeqSero2-v1.3.2/bin/SeqSero2_package.py"
    sistr = "/usr/local/sistr-1.1.3/bin/sistr"
    mlst = "/usr/local/mlst-2.23.0/bin/mlst"

    print("")
    print(f"  Checking for required programs")
    print(f"  ==============================")
    commands = ['blastn', 'any2fasta', 'samtools', 'spades.py', 'SalmID.py', 'seqkit', 'bedtools', 'mafft', 'aggregate_salmonella_typing_results.py']
    tools = ['mlst', 'mlst', 'seqsero2', 'seqsero2', 'seqsero2', 'seqsero2', 'seqsero2', 'sistr', 'vmst']

    for i,cmd in enumerate(commands):
        command_path = which(cmd)
        if command_path:
            print(f"  '{cmd}' found at {command_path} [{tools[i]}]")
        else:
            print(f"  '{cmd}' not found, please make sure command is accessible in $PATH.")
    print("")

    # initialize
    pg = PipelineGenerator(make_file)

    analyse_ok_files = ""

    # analyze
    for idx, sample in enumerate(samples):

        # serovar and antigen
        # /usr/local/SeqSero2/bin/SeqSero2_package.py -d results -n 23_1704 -t 2 -i r1.fastq.gz r2.fastq.gz
        out_dir = f"{output_dir}/{sample.id}/seqsero2"
        log = f"{out_dir}/run.log"
        err = f"{out_dir}/run.err"
        tgt = f"{out_dir}/seqsero2.OK"
        dep = ""
        cmd = f'{seqsero2} -d {out_dir} -n {sample.id} -t 2 -i {sample.fastq1} {sample.fastq2} > {log} 2> {err}'
        analyse_ok_files += f"{tgt} "
        pg.add(tgt, dep, cmd)

        # sequence typing
        out_dir = f"{output_dir}/{sample.id}/mlst"
        input_contig_fasta_file = f"{sample.contigs_fasta}"
        results_txt_file = f"{out_dir}/results.txt"
        log = f"{out_dir}/mlst.log"
        tgt = f"{out_dir}/mlst.OK"
        dep = f""
        cmd = f'{mlst} {input_contig_fasta_file} --json {out_dir}/typing.json --scheme senterica_achtman_2  --nopath > {results_txt_file} 2> {log}'
        analyse_ok_files += f"{tgt} "
        pg.add(tgt, dep, cmd)

        # SISTR
        out_dir = f"{output_dir}/{sample.id}/sistr"
        output_file = f"{out_dir}/sistr"
        input_contig_fasta_file = f"{sample.contigs_fasta}"
        log = f"{out_dir}/run.log"
        err = f"{out_dir}/run.err"
        tgt = f"{out_dir}/sistr.OK"
        dep = f""
        cmd = f'{sistr} {input_contig_fasta_file} -f tab -o {output_file} -K -T {out_dir} > {log} 2> {err}'
        analyse_ok_files += f"{tgt} "
        pg.add(tgt, dep, cmd)

    #aggregate files
    output_xlsx = f"{output_dir}/summary.xlsx"
    tgt = f"{output_dir}/summary.xlsx.OK"
    dep = analyse_ok_files
    cmd = f'aggregate_salmonella_typing_results.py -i {output_dir} -s {sample_file} -o {output_xlsx}'
    analyse_ok_files += f"{tgt} "
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
    def __init__(self, idx, id, fastq1, fastq2, contigs_fasta):
        self.idx = idx
        self.id = id
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.contigs_fasta = contigs_fasta

    def print(self):
        print(f"index           : {self.idx}")
        print(f"id              : {self.id}")
        print(f"fastq1          : {self.fastq1}")
        print(f"fastq2          : {self.fastq2}")
        print(f"contigs fasta   : {self.contigs_fasta}")


if __name__ == "__main__":
    main() # type: ignore