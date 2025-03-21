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
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import click
import random

@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_structure_pca_analysis.mk",
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
@click.option("-i", "--input_vcf_file", required=True, help="VCF file")
@click.option("-d", "--dataset", required=True, help="dataset name")
def main(make_file, working_dir, sample_file, input_vcf_file, dataset):
    """
    Population structure of Pangolins

    e.g. generate_structure_pca_analysis_pipeline.py
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("input VCF file", input_vcf_file))
    print("\t{0:<20} :   {1:<10}".format("dataset", dataset))

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    working_dir = os.path.abspath(working_dir)
    ref_dir = f"{working_dir}/ref"
    log_dir = f"{working_dir}/log"
    stats_dir = f"{working_dir}/stats"
    plot_dir = f"{working_dir}/plot"
    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(stats_dir, exist_ok=True)
        os.makedirs(f"{stats_dir}/coverage", exist_ok=True)
        os.makedirs(f"{stats_dir}/general", exist_ok=True)
        os.makedirs(f"{stats_dir}/flag", exist_ok=True)
        os.makedirs(f"{stats_dir}/idx", exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    ##########
    # programs
    ##########
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    bcftools = "/usr/local/bcftools-1.17/bin/bcftools"
    structure = "/usr/local/structure-2.3.4/structure"
    fpca = "/usr/local/fratools-1.0/fpca"
    distruct = "/usr/local/distruct-1.1/distruct"

    script_dir = os.path.dirname(__file__)
    vcf_to_structure = f"{script_dir}/vcf_to_structure.py"
    vcf_to_tg = f"{script_dir}/vcf_to_tg.py"
    structure_to_clumpp_distruct = f"{script_dir}/structure_to_clumpp_distruct.py"
    structure_gis_to_sa = f"{script_dir}/structure_gis_to_sa.py"
    structure_pca_to_sa = f"{script_dir}/structure_pca_to_sa.py"
    plot_gis_structure = f"{script_dir}/plot_gis_structure.py"
    plot_pca_structure = f"{script_dir}/plot_pca_structure.py"

    # create directories in destination folder directory
    structure_dir = f"{working_dir}/structure"
    pca_dir = f"{working_dir}/pca"
    try:
        os.makedirs(structure_dir, exist_ok=True)
        os.makedirs(f"{structure_dir}/barplots", exist_ok=True)
        os.makedirs(f"{structure_dir}/gisplots", exist_ok=True)
        os.makedirs(pca_dir, exist_ok=True)
        os.makedirs(f"{pca_dir}/gisplots", exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    ##########
    #structure
    ##########
    #convert VCF file to structure format
    output_dir = structure_dir
    log = f"{output_dir}/structure_files.log"
    tgt = f"{output_dir}/structure_files.OK"
    dep = f""
    cmd = f"{vcf_to_structure} {input_vcf_file} -o {output_dir} > {log}"
    pg.add(tgt, dep, cmd)

    #run structure
    random.seed(3323)
    for k in range(2, 5):
        for rep in range(1, 4):
            seed = int(1000*random.random())
            input_structure_file = f"{structure_dir}/{dataset}.structure"
            output_structure_results_file = f"{structure_dir}/K{k}_R{rep}"
            mainparams = f"{structure_dir}/mainparams"
            extraparams = f"{structure_dir}/extraparams"
            log = f"{structure_dir}/K{k}_R{rep}.log"
            tgt = f"{structure_dir}/K{k}_R{rep}.OK"
            dep = f"{structure_dir}/structure_files.OK"
            cmd = f"{structure} -i {input_structure_file} -o {output_structure_results_file} -m {mainparams} -e {extraparams} -K {k} -D {seed} > {log}"
            pg.add(tgt, dep,  cmd)

        #prepare distruct files
        input_structure_results_files = ""
        log = f"{structure_dir}/K{k}.distruct.log"
        tgt = f"{structure_dir}/K{k}.distruct.OK"
        dep = ""
        for rep in range(1, 4):
            input_structure_results_files += f"{structure_dir}/K{k}_R{rep}_f "
            dep += f"{structure_dir}/K{k}_R{rep}.OK "
        cmd = f"{structure_to_clumpp_distruct} {input_structure_results_files} > {log}"
        pg.add(tgt, dep, cmd)

        #run distruct
        for rep in range(1, 4):
            input_drawparam = f"K{k}_R{rep}.drawparams"
            log = f"{structure_dir}/K{k}_R{rep}.distruct.log"
            tgt = f"{structure_dir}/K{k}_R{rep}.distruct.OK"
            dep = f"{structure_dir}/K{k}.distruct.OK "
            cmd = f"cd {structure_dir}; {distruct} -d {input_drawparam} > {log}; set $? 0"
            pg.add(tgt, dep, cmd)

            input_ps_file = f"{structure_dir}/K{k}_R{rep}.ps"
            output_pdf_file = f"{structure_dir}/barplots/K{k}_R{rep}.pdf"
            tgt = f"{structure_dir}/K{k}_R{rep}.pdf.ok"
            dep = f"{structure_dir}/K{k}_R{rep}.distruct.OK "
            cmd = f"ps2pdf {input_ps_file} {output_pdf_file}"
            pg.add(tgt, dep, cmd)

            #generate sample files for plotting GIS scatterplots
            input_structure_file = f"{structure_dir}/K{k}_R{rep}_f"
            output_sa_file = f"{structure_dir}/K{k}_R{rep}.sa"
            log = f"{structure_dir}/K{k}_R{rep}.sa.log"
            tgt = f"{structure_dir}/K{k}_R{rep}.sa.OK"
            dep = f"{structure_dir}/K{k}_R{rep}.OK "
            cmd = f"{structure_gis_to_sa} -g {sample_file} -s {input_structure_file} -o {output_sa_file} > {log}"
            pg.add(tgt, dep, cmd)

            #plot geospatial plot with structure pie charts
            input_sa_file = f"{structure_dir}/K{k}_R{rep}.sa"
            output_r_file = f"{structure_dir}/K{k}_R{rep}.r"
            output_pdf_file = f"{structure_dir}/gisplots/K{k}_R{rep}.pdf"
            log = f"{structure_dir}/K{k}_R{rep}_gis.log"
            tgt = f"{structure_dir}/K{k}_R{rep}_gis.pdf.OK"
            dep = f"{structure_dir}/K{k}_R{rep}.sa.OK"
            cmd = f"{plot_gis_structure} {input_sa_file} -o {output_r_file} -z {output_pdf_file} > {log}"
            pg.add(tgt, dep, cmd)

    ####
    #PCA
    ####
    #convert VCF file to tg format
    output_tg_file = f"{pca_dir}/{dataset}.tg"
    log = f"{pca_dir}/pca_files.log"
    tgt = f"{output_tg_file}.OK"
    dep = f""
    cmd = f"{vcf_to_tg} {input_vcf_file} -o {output_tg_file} > {log}"
    pg.add(tgt, dep, cmd)

    #pca
    input_tg_file = f"{pca_dir}/{dataset}.tg"
    log = f"{pca_dir}/pca.log"
    err = f"{pca_dir}/pca.err"
    tgt = f"{pca_dir}/pca.OK"
    dep = f"{input_tg_file}.OK"
    cmd = f"cd {pca_dir}; {fpca} -i {input_tg_file} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    #generate sample files for plotting PCA - Structure scatterplots
    for k in range(2, 5):
        for rep in range(1, 4):
            input_pca_file = f"{pca_dir}/{dataset}.pca"
            input_structure_file = f"{structure_dir}/K{k}_R{rep}_f"
            output_sa_file = f"{pca_dir}/K{k}_R{rep}.sa"
            log = f"{output_sa_file}.log"
            tgt = f"{output_sa_file}.OK"
            dep = f"{pca_dir}/pca.OK {structure_dir}/K{k}_R{rep}.OK"
            cmd = f"{structure_pca_to_sa} -s {input_structure_file} -p {input_pca_file} -o {output_sa_file} > {log}"
            pg.add(tgt, dep, cmd)

            #plot PCA with structure pie charts
            #plot geospatial plot with structure pie charts
            input_sa_file = f"{pca_dir}/K{k}_R{rep}.sa"
            output_pdf_file = f"{pca_dir}/gisplots/K{k}_R{rep}.pdf"
            log = f"{pca_dir}/K{k}_R{rep}_gis.log"
            tgt = f"{pca_dir}/K{k}_R{rep}.pdf.OK"
            dep = f"{input_sa_file}.OK"
            cmd = f"{plot_pca_structure} {input_sa_file} -o {pca_dir}/K{k}_R{rep}.r -z {output_pdf_file} > {log}"
            pg.add(tgt, dep, cmd)

    # clean
    pg.add_clean(f"rm -fr {ref_dir} {log_dir} {stats_dir} {plot_dir} {structure_dir} {pca_dir}")

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
