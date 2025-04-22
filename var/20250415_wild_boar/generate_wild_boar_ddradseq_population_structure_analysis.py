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
    default="run_wild_boar_ddradseq_population_structure_analysis.mk",
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
    Population structure of Wild Boars

    e.g. generate_wild_boar_ddradseq_population_structure_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample file", sample_file))
    print("\t{0:<20} :   {1:<10}".format("population map file", population_map_file))
    print("\t{0:<20} :   {1:<10}".format("genome fasta file", genome_fasta_file))

    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                id, bam_file = line.rstrip().split("\t")
                samples.append(Sample(id, bam_file))

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    log_dir = f"{working_dir}/log"
    plot_dir = f"{working_dir}/plot"
    bam_dir = f"{working_dir}/bam"
    vcf_dir = f"{working_dir}/vcf"
    annotation_dir = f"{working_dir}/annotations"
    ref_stacks_dir = f"{working_dir}/ref_stacks"
    qc_dir = f"{working_dir}/qc"
    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)
        os.makedirs(vcf_dir, exist_ok=True)
        os.makedirs(ref_stacks_dir, exist_ok=True)
        os.makedirs(qc_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    ##########
    # programs
    ##########
    ref_stacks = "/usr/local/stacks-2.68/bin/ref_map.pl"
    structure = "/usr/local/structure-2.3.4/structure"
    fpca = "/usr/local/fratools-1.0/fpca"
    distruct = "/usr/local/distruct-1.1/distruct"
    bcftools = "/usr/local/bcftools-1.17/bin/bcftools"
    script_dir = "/home/atks/programs/CAVS-pipelines/var/20250415_wild_boar"
    filter_ddradseq_vcf = f"{script_dir}/filter_ddradseq_vcf.py"
    qwikplot = f"{script_dir}/qwikplot"
    vcf_to_structure = f"{script_dir}/vcf_to_structure.py"
    vcf_to_plink = f"{script_dir}/vcf_to_plink.py"
    vcf_to_tg = f"{script_dir}/vcf_to_tg.py"
    structure_to_clumpp_distruct = f"{script_dir}/structure_to_clumpp_distruct.py"
    structure_pca_to_sa = f"{script_dir}/structure_pca_to_sa.py"
    plot_pca_structure = f"{script_dir}/plot_pca_structure.py"

    ####################
    # bam files staging
    ####################
    for sample in samples:
        target = sample.bam_file
        link =f"{bam_dir}/{sample.id}.bam"
        if os.path.islink(link) or os.path.exists(link):
            os.unlink(link)
        os.symlink(target, link)
        target = f"{sample.bam_file}.bai"
        link =f"{bam_dir}/{sample.id}.bam.bai"
        if os.path.islink(link) or os.path.exists(link):
            os.unlink(link)
        os.symlink(target, link)

    #################
    # reference files
    #################
    # reference genome assembly for sunda pangolin
    # https://www.ncbi.nlm.nih.gov/assembly/GCF_014570535.1

    #copy reference fasta file to reference directory
    ref_fasta_file = f"{ref_dir}/{os.path.basename(genome_fasta_file)}"
    tgt = f"{ref_fasta_file}.OK"
    dep = ""
    cmd = f"cp {genome_fasta_file} {ref_fasta_file}"
    pg.add(tgt, dep, cmd)

    ################################
    # Variant Calling
    ################################

    #ref stacks
    log = f"{log_dir}/ref_stacks.log"
    tgt = f"{log_dir}/ref_stacks.OK"
    dep = ""
    cmd = f"{ref_stacks} -T 30 -X \"populations: --vcf\" -o {ref_stacks_dir} --popmap {population_map_file} --samples {bam_dir} 2> {log}"
    pg.add(tgt, dep, cmd)

    ################################
    # Quality Checking and Filtering
    ################################

    #filter SNPs
    input_vcf_file = f"{ref_stacks_dir}/populations.snps.vcf"
    log = f"{log_dir}/qc.log"
    tgt = f"{log_dir}/qc.OK"
    dep = f"{log_dir}/ref_stacks.OK"
    cmd = f"{filter_ddradseq_vcf} {input_vcf_file} -o {qc_dir} -s 0.1 > {log}"
    pg.add(tgt, dep, cmd)

    #draw qc plots of call rates and mafs
    for i in range(1, 4):
        input_txt_file = f"{qc_dir}/sample_call_rate_iter_{i}.txt"
        output_pdf_file = f"{qc_dir}/sample_call_rate_iter_{i}.pdf"
        tgt = f"{output_pdf_file}.OK"
        dep = f"{log_dir}/qc.OK"
        cmd = f"{qwikplot} {input_txt_file} --ylim \"c(0,1)\" -y sample_call_rate -t \"Iteration {i} Sample Call Rate\" -z {output_pdf_file} -p 20 -c green >/dev/null"
        pg.add(tgt, dep, cmd)

        input_txt_file = f"{qc_dir}/snp_call_rate_iter_{i}.txt"
        output_pdf_file = f"{qc_dir}/snp_call_rate_iter_{i}.pdf"
        tgt = f"{output_pdf_file}.OK"
        dep = f"{log_dir}/qc.OK"
        cmd = f"{qwikplot} {input_txt_file} --ylim \"c(0,1)\" -y variant_call_rate -t \"Iteration {i} SNP Call Rate\" -z {output_pdf_file} -p 20 -c green -q 0.1 >/dev/null"
        pg.add(tgt, dep, cmd)

        input_txt_file = f"{qc_dir}/maf_iter_{i}.txt"
        output_pdf_file = f"{qc_dir}/maf_iter_{i}.pdf"
        tgt = f"{output_pdf_file}.OK"
        dep = f"{log_dir}/qc.OK"
        cmd = f"{qwikplot} {input_txt_file} --ylim \"c(0,0.5)\" -y variant_maf -t \"Iteration {i} MAF\" -z {output_pdf_file} -p 20 -c green -q 0.1 >/dev/null"
        pg.add(tgt, dep, cmd)

    ##################
    #Relative analysis
    ##################

    #data sets

    #compute heterozygosity

    #compute pairwise IBS stats

    #plot heterozygosity plots

    #plot pairwise IBS plots

    ##################
    #prepare data sets
    ##################

    #copy over 42 samples, 48514 variants data set
    input_vcf_file = f"{qc_dir}/populations.snps.filtered.vcf"
    output_vcf_file = f"{vcf_dir}/42samples_48514snps_wild_boar.vcf"
    tgt = f"{output_vcf_file}.OK"
    dep = f"{log_dir}/qc.OK"
    cmd = f"cp {input_vcf_file} {output_vcf_file}"
    pg.add(tgt, dep, cmd)

    for dataset in ["42samples_48514snps"]:
        # create directories in destination folder directory
        structure_dir = f"{working_dir}/{dataset}/structure"
        pca_dir = f"{working_dir}/{dataset}/pca"
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
        input_vcf_file = f"{vcf_dir}/{dataset}_wild_boar.vcf"
        output_dir = f"{working_dir}/{dataset}/structure"
        tgt = f"{output_dir}/structure_files.OK"
        dep = f"{input_vcf_file}.OK"
        cmd = f"{vcf_to_structure} {input_vcf_file} -o {output_dir}"
        pg.add(tgt, dep, cmd)

        #run structure
        random.seed(3323)
        for k in range(2, 5):
            for rep in range(1, 4):
                seed = int(1000*random.random())
                output_dir = f"{working_dir}/{dataset}/structure"
                input_structure_file = f"{output_dir}/{dataset}_wild_boar.structure"
                output_structure_results_file = f"{output_dir}/K{k}_R{rep}"
                mainparams = f"{output_dir}/mainparams"
                extraparams = f"{output_dir}/extraparams"
                log = f"{output_dir}/K{k}_R{rep}.log"
                tgt = f"{output_dir}/K{k}_R{rep}.OK"
                dep = f"{output_dir}/structure_files.OK"
                cmd = f"{structure} -i {input_structure_file} -o {output_structure_results_file} -m {mainparams} -e {extraparams} -K {k} -D {seed} > {log}"
                pg.add(tgt, dep, cmd)

        #prepare distruct files
        for k in range(2, 5):
            output_dir = f"{working_dir}/{dataset}/structure"
            input_structure_results_files = ""
            log = f"{output_dir}/K{k}.distruct.log"
            tgt = f"{output_dir}/K{k}.distruct.OK"
            dep = ""
            for rep in range(1, 4):
                input_structure_results_files += f"{output_dir}/K{k}_R{rep}_f "
                dep += f"{output_dir}/K{k}_R{rep}.OK "
            cmd = f"{structure_to_clumpp_distruct} {input_structure_results_files} > {log}"
            pg.add(tgt, dep, cmd)


        for k in range(2, 5):
            for rep in range(1, 4):
                #run distruct
                input_drawparam = f"K{k}_R{rep}.drawparams"
                log = f"{output_dir}/K{k}_R{rep}.distruct.log"
                tgt = f"{output_dir}/K{k}_R{rep}.distruct.OK"
                dep = f"{output_dir}/K{k}.distruct.OK "
                cmd = f"cd {output_dir}; {distruct} -d {input_drawparam} > {log}; set $? 0"
                pg.add(tgt, dep, cmd)

                input_ps_file = f"{output_dir}/K{k}_R{rep}.ps"
                output_pdf_file = f"{output_dir}/K{k}_R{rep}.pdf"
                tgt = f"{output_dir}/K{k}_R{rep}.pdf.ok"
                dep = f"{output_dir}/K{k}_R{rep}.distruct.OK "
                cmd = f"ps2pdf {input_ps_file} {output_pdf_file}"
                pg.add(tgt, dep, cmd)

        ####
        #PCA
        ####
        #convert VCF file to tg format
        input_vcf_file = f"{vcf_dir}/{dataset}_wild_boar.vcf"
        output_dir = f"{working_dir}/{dataset}/pca"
        tgt = f"{output_dir}/pca_files.OK"
        dep = f"{input_vcf_file}.OK"
        cmd = f"{vcf_to_tg} {input_vcf_file} -o {output_dir}"
        pg.add(tgt, dep, cmd)

        #pca
        output_dir = f"{working_dir}/{dataset}/pca"
        input_tg_file = f"{output_dir}/{dataset}_wild_boar.tg"
        log = f"{output_dir}/pca.log"
        err = f"{output_dir}/pca.err"
        tgt = f"{output_dir}/pca.OK"
        dep = f"{input_vcf_file}.OK"
        cmd = f"cd {output_dir}; {fpca} -i {input_tg_file} > {log} 2> {err}"
        pg.add(tgt, dep, cmd)

        #generate sample files for plotting PCA - Structure scatterplots
        for k in range(2, 5):
            output_dir = f"{working_dir}/{dataset}/pca"
            input_pca_file = f"{output_dir}/{dataset}_wild_boar.pca"
            input_structure_file = f"{structure_dir}/K{k}_R1_f"
            output_sa_file = f"{output_dir}/gisplots/K{k}.sa"
            log = f"{output_sa_file}.log"
            tgt = f"{output_sa_file}.OK"
            dep = f"{output_dir}/pca.OK {structure_dir}/K{k}_R1.OK"
            cmd = f"{structure_pca_to_sa} -s {input_structure_file} -p {input_pca_file} -o {output_sa_file} > {log}"
            pg.add(tgt, dep, cmd)

            #plot PCA with structure pie charts
            input_sa_file = f"{output_dir}/gisplots/K{k}.sa"
            output_pdf_file = f"{output_dir}/gisplots/K{k}.pdf"
            tgt = f"{output_pdf_file}.OK"
            dep = f"{output_dir}/gisplots/K{k}.sa"
            cmd = f"{plot_pca_structure} {input_sa_file} -p K{k} -o {output_dir}/gisplots -z {output_pdf_file}"
            pg.add(tgt, dep, cmd)

    # clean
    pg.add_clean(f"rm -fr {ref_dir}")

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

    def __init__(self, id, bam_file):
        self.id = id
        self.bam_file = bam_file

    def print(self):
        print(f"id       : {self.id}")
        print(f"bam files : {len(self.bam_file)}")


if __name__ == "__main__":
    main() # type: ignore
