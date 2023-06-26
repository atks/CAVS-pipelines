#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2021 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import argparse
import textwrap
import gzip
import re
import sys
from shutil import copy2
from datetime import datetime


def main():
    # programs
    guppy_base_caller = "/usr/local/ont-guppy-6.1.5/bin/guppy_basecaller"
    guppy_barcoder = "/usr/local/ont-guppy-6.1.5/bin/guppy_barcoder"
    nanoplot = "/usr/local/bin/NanoPlot"

    # options
    parser = argparse.ArgumentParser(
        description="Moves ONT fastq files to a destination and performs QC",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
           usage: generate_ont_deploy_and_qc_pipeline -m <make_file> -d <dest_dir>

           """
        ),
    )
    parser.add_argument(
        "-m",
        "--make_file",
        help="make file name",
        type=str,
        default="ont_basecall_deploy_and_qc.mk",
    )
    parser.add_argument(
        "-d", "--dest_dir", help="destination directory", type=str, required=True
    )
    parser.add_argument(
        "-w",
        "--working_dir",
        help="working directory",
        type=str,
        required=False,
        default=os.getcwd(),
    )
    parser.add_argument(
        "-r", "--runs_file", help="sample runs file", type=str, required=True
    )
    parser.add_argument(
        "-s", "--sample_file", help="sample file", type=str, required=True
    )
    args = parser.parse_args()

    for arg in vars(args):
        print("\t{0:<20} :   {1:<10}".format(arg, getattr(args, arg) or ""))

    try:
        os.makedirs(f"{args.working_dir}/log", exist_ok=True)
    except OSError as error:
        print(f"Directory {args.working_dir} cannot be created")

    RUNS = dict()
    file = open(args.runs_file, "r")
    for line in file:
        if not line.startswith("#"):
            run_id, sub_run_id, dir = line.rstrip().split("\t")
            if run_id not in RUNS:
                RUNS[run_id] = Run(run_id)
            RUNS[run_id].add_sub_run(sub_run_id, dir)
    file.close()

    file = open(args.sample_file, "r")
    for line in file:
        if not line.startswith("#"):
            run_id, sample_id, barcode = line.rstrip().split("\t")
            if run_id not in RUNS:
                sys.exit("Run ID not in runs file")
            else:
                RUNS[run_id].add_sample(sample_id, barcode)
    file.close()

    #    for run in RUNS.values():
    #        run.print()
    #    sys.exit()

    # initialize
    pg = PipelineGenerator(args.make_file)

    working_dir = args.working_dir
    log_dir = f"{args.working_dir}/log"
    previous_run = ""
    for run in RUNS.values():
        for sub_run in run.sub_runs:
            # base call
            output_dir = f"{working_dir}/{run.id}/{sub_run.id}"
            log = f"{log_dir}/{run.id}_{sub_run.id}.log"
            err = f"{log_dir}/{run.id}_{sub_run.id}.err"
            tgt = f"{log_dir}/{run.id}_{sub_run.id}.OK"
            dep = previous_run
            cmd = f'{guppy_base_caller} -i "{sub_run.dir}" -r -s {output_dir} --flowcell FLO-MIN106  --kit SQK-LSK109 -x auto  > {log} 2> {err}'
            pg.add(tgt, dep, cmd)
            previous_run = f"{log_dir}/{run.id}_{sub_run.id}.OK"

            # demux
            input_dir = f"{working_dir}/{run.id}/{sub_run.id}/pass"
            output_dir = f"{working_dir}/{run.id}/{sub_run.id}/demux"
            log = f"{log_dir}/{run.id}_{sub_run.id}.demux.log"
            err = f"{log_dir}/{run.id}_{sub_run.id}.demux.err"
            dep = f"{log_dir}/{run.id}_{sub_run.id}.OK"
            tgt = f"{log_dir}/{run.id}_{sub_run.id}.demux.OK"
            cmd = f'{guppy_barcoder} -i "{input_dir}" -r -s {output_dir} --barcode_kits EXP-NBD104 --compress_fastq -t 2 > {log} 2> {err}'
            pg.add(tgt, dep, cmd)

    # create directories in destination folder directory
    for run in RUNS.values():
        if len(run.sub_runs) > 1:
            for sub_run in run.sub_runs:
                for sample in run.samples:
                    new_dir = f"{args.dest_dir}/{run.id}/{sub_run.id}/{sample.id}"
                    # print(f'Directory {new_dir} created')
                    try:
                        os.makedirs(f"{new_dir}", exist_ok=True)
                    except OSError as error:
                        print(f"Directory {new_dir} cannot be created")
        else:
            for sample in run.samples:
                new_dir = f"{args.dest_dir}/{run.id}/{sample.id}"
                # print(f'Directory {new_dir} created')
                try:
                    os.makedirs(f"{new_dir}", exist_ok=True)
                except OSError as error:
                    print(f"Directory {new_dir} cannot be created")

    # create directories in destination folder directory
    for run in RUNS.values():
        if len(run.sub_runs) > 1:
            for sub_run in run.sub_runs:
                for sample in run.samples:
                    input_dir = (
                        f"{working_dir}/{run.id}/{sub_run.id}/demux/{sample.barcode}"
                    )
                    output_file = f"{args.dest_dir}/{run.id}/{sub_run.id}/{sample.id}/{sample.id}.fastq.gz"
                    log = f"{log_dir}/{run.id}_{sub_run.id}.zcat.log"
                    err = f"{log_dir}/{run.id}_{sub_run.id}.zcat.err"
                    dep = f"{log_dir}/{run.id}_{sub_run.id}.demux.OK"
                    tgt = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.zcat.OK"
                    cmd = f"zcat {input_dir}/*.fastq.gz | gzip -c > {output_file}"
                    pg.add(tgt, dep, cmd)
                    # print(f"\twriting {output_file}")

                    input_file = f"{args.dest_dir}/{run.id}/{sub_run.id}/{sample.id}/{sample.id}.fastq.gz"
                    output_dir = f"{args.dest_dir}/{run.id}/{sub_run.id}/{sample.id}/nanoplot_result"
                    log = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.nanoplot.log"
                    err = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.nanoplot.err"
                    dep = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.zcat.OK"
                    tgt = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.nanoplot.OK"
                    cmd = f"{nanoplot} --fastq {input_file} -o {output_dir} > {log} 2> {err}"
                    pg.add(tgt, dep, cmd)

        else:
            for sample in run.samples:
                input_dir = f"{working_dir}/{run.id}/run1/demux/{sample.barcode}"
                output_file = (
                    f"{args.dest_dir}/{run.id}/{sample.id}/{sample.id}.fastq.gz"
                )
                log = f"{log_dir}/{run.id}_run1.zcat.log"
                err = f"{log_dir}/{run.id}_run1.zcat.err"
                dep = f"{log_dir}/{run.id}_run1.demux.OK"
                tgt = f"{log_dir}/{run.id}_run1_{sample.id}.zcat.OK"
                cmd = f"zcat {input_dir}/*.fastq.gz | gzip -c > {output_file}"
                pg.add(tgt, dep, cmd)
                # print(f"\twriting {output_file}")

                input_file = (
                    f"{args.dest_dir}/{run.id}/{sample.id}/{sample.id}.fastq.gz"
                )
                output_dir = f"{args.dest_dir}/{run.id}/{sample.id}/nanoplot_result"
                log = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.nanoplot.log"
                err = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.nanoplot.err"
                dep = f"{log_dir}/{run.id}_run1_{sample.id}.zcat.OK"
                tgt = f"{log_dir}/{run.id}_{sub_run.id}_{sample.id}.nanoplot.OK"
                cmd = (
                    f"{nanoplot} --fastq {input_file} -o {output_dir} > {log} 2> {err}"
                )
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

    def add(self, tgt, dep, cmd):
        self.tgts.append(tgt)
        self.deps.append(dep)
        self.cmds.append(cmd)

    def print(self):
        print(".DELETE_ON_ERROR:")
        for i in range(len(self.tgts)):
            print(f"{self.tgts[i]} : {self.deps[i]}")
            print(f"\t{self.cmds[i]}")
            print(f"\ttouch {self.tgts[i]}")

    def write(self):
        f = open(self.make_file, "w")
        f.write(".DELETE_ON_ERROR:\n\n")
        f.write("all : ")
        for i in range(len(self.tgts)):
            f.write(f"{self.tgts[i]} ")
        f.write("\n\n")
        for i in range(len(self.tgts)):
            f.write(f"{self.tgts[i]} : {self.deps[i]}\n")
            f.write(f"\t{self.cmds[i]}\n")
            f.write(f"\ttouch {self.tgts[i]}\n\n")
        f.close()


class Sample(object):
    def __init__(self):
        self.id = ""
        self.barcode = ""

    def __init__(self, id, barcode):
        self.id = id
        self.barcode = barcode

    def print(self):
        print(f"id      : {self.id}")
        print(f"barcode : {self.barcode}")


class Run(object):
    def __init__(self, id):
        self.id = id
        self.sub_runs = []
        self.samples = []

    def add_sub_run(self, sub_run_id, dir):
        self.sub_runs.append(SubRun(sub_run_id, dir))

    def add_sample(self, sample_id, barcode):
        self.samples.append(Sample(sample_id, barcode))

    def print(self):
        print(f"++++++++++++++++++++")
        print(f"run id       : {self.id}")
        print(f"++++++++++++++++++++")
        for sub_run in self.sub_runs:
            sub_run.print()
        print(f"++++++++++++++++++++")
        for sample in self.samples:
            sample.print()
        print(f"++++++++++++++++++++")


class SubRun(object):
    def __init__(self):
        self.id = ""
        self.dir = ""

    def __init__(self, id, dir):
        self.id = id
        self.dir = dir

    def print(self):
        print(f"sub run id   : {self.id}")
        print(f"dir          : {self.dir}")


main()
