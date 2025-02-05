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
    denovo_stacks_dir = f"{working_dir}/denovo_stacks"
    ref_stacks_dir = f"{working_dir}/ref_stacks"
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
        os.makedirs(denovo_stacks_dir, exist_ok=True)
        os.makedirs(ref_stacks_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    ##########
    # programs
    ##########
    multiqc = "docker run  -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc "
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.17/bin/samtools"
    plot_bamstats = "/usr/local/samtools-1.17/bin/plot-bamstats"
    compute_effective_coverage = "/home/atks/programs/cavspipes/vfp/compute_effective_coverage.py"
    extract_general_stats = "/home/atks/programs/cavspipes/vfp/extract_general_stats.py"
    process_shortreads = "/usr/local/stacks-2.68/bin/process_shortreads"
    denovo_stacks = "/usr/local/stacks-2.68/bin/denovo_map.pl"
    ref_stacks = "/usr/local/stacks-2.68/bin/ref_map.pl"

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

    #index reference sequence
    log = f"{ref_dir}/bwa_index.log"
    tgt = f"{ref_dir}/bwa_index.OK"
    dep = f"{ref_fasta_file}.OK"
    cmd = f"{bwa} index -a bwtsw {ref_fasta_file} 2> {log}"
    pg.add(tgt, dep, cmd)

    fastq_files_OK = ""
    bam_files_OK = ""
    samtools_multiqc_dep = ""

    for id, sample in samples.items():

        #combine files
        if len(sample.fastq1s) == 1:
            input_fastq1_file = sample.fastq1s[0]
            output_fastq1_file = f"{fastq_dir}/{sample.id}.1.fq.gz"
            tgt = f"{output_fastq1_file}.OK"
            fastq_files_OK += f"{tgt} "
            dep = ""
            cmd = f"ln -s {input_fastq1_file} {output_fastq1_file}"
            pg.add(tgt, dep, cmd)

            input_fastq2_file = sample.fastq2s[0]
            output_fastq2_file = f"{fastq_dir}/{sample.id}.2.fq.gz"
            tgt = f"{output_fastq2_file}.OK"
            fastq_files_OK += f"{tgt} "
            dep = ""
            cmd = f"ln -s {input_fastq2_file} {output_fastq2_file}"
            pg.add(tgt, dep, cmd)

        elif len(sample.fastq1s) > 1:
            input_fastq1_files = " ".join(sample.fastq1s)
            output_fastq1_file = f"{fastq_dir}/{sample.id}.1.fq.gz"
            tgt = f"{output_fastq1_file}.OK"
            fastq_files_OK += f"{tgt} "
            dep = ""
            cmd = f"zcat {input_fastq1_files} | gzip > {output_fastq1_file}"
            pg.add(tgt, dep, cmd)

            input_fastq2_files = " ".join(sample.fastq2s)
            output_fastq2_file = f"{fastq_dir}/{sample.id}.2.fq.gz"
            tgt = f"{output_fastq2_file}.OK"
            fastq_files_OK += f"{tgt} "
            dep = ""
            cmd = f"zcat {input_fastq2_files} | gzip > {output_fastq2_file}"
            pg.add(tgt, dep, cmd)

        # #filter short reads
        # input_fastq1_file = f"{fastq_dir}/{sample.id}.1.fq.gz"
        # input_fastq2_file = f"{fastq_dir}/{sample.id}.2.fq.gz"
        # output_fastq2_file = f"{fastq_dir}/{sample.id}.2.fq.gz"
        # tgt = f"{output_fastq2_file}.OK"
        # fastq_files_OK += f"{tgt} "
        # dep = f"{fastq_dir}/{sample.id}.1.fq.gz.OK {fastq_dir}/{sample.id}.2.fq.gz.OK"
        # cmd = f"{process_shortreads} -1 fastq/BIOS0006.1.fq.gz -2 fastq/BIOS0006.2.fq.gz -o {fastq_ulen_dir} --len-limit 145 "
        # pg.add(tgt, dep, cmd)

        # align
        src_fastq1 = f"{fastq_dir}/{sample.id}.1.fq.gz"
        src_fastq2 = f"{fastq_dir}/{sample.id}.2.fq.gz"
        output_bam_file = f"{bam_dir}/{sample.id}.bam"
        log = f"{log_dir}/{sample.id}.align.log"
        sort_log = f"{log_dir}/{sample.id}.align.sort.log"
        dep = f"{src_fastq1}.OK {src_fastq2}.OK {ref_dir}/bwa_index.OK"
        tgt = f"{output_bam_file}.OK"
        cmd = f"{bwa} mem -t 2 -M {ref_fasta_file} {src_fastq1} {src_fastq2} 2> {log} | {samtools} view -h | {samtools} sort -o {output_bam_file} 2> {sort_log}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # index
        input_bam_file = f"{bam_dir}/{sample.id}.bam"
        dep = f"{bam_dir}/{sample.id}.bam.OK"
        tgt = f"{bam_dir}/{sample.id}.bam.bai.OK"
        cmd = f"{samtools} index {input_bam_file}"
        pg.add(tgt, dep, cmd)

        #  coverage
        output_stats_file = f"{stats_dir}/coverage/{sample.id}.txt"
        dep = f"{bam_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{output_stats_file}.OK"
        cmd = f"{samtools} coverage {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # compute effective coverage
        input_stats_file = f"{stats_dir}/coverage/{sample.id}.txt"
        output_stats_file = f"{stats_dir}/coverage/{sample.id}.ecov.txt"
        dep = f"{input_stats_file}.OK"
        tgt = f"{output_stats_file}.OK"
        cmd = f"{compute_effective_coverage} {input_stats_file} -o {output_stats_file} -s {sample.id}"
        pg.add(tgt, dep, cmd)

        # stats
        output_stats_file = f"{stats_dir}/general/{sample.id}.txt"
        dep = f"{bam_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{output_stats_file}.OK"
        cmd = f"{samtools} stats {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # extract general stats
        input_stats_file = f"{stats_dir}/general/{sample.id}.txt"
        output_stats_file = f"{stats_dir}/general/{sample.id}.extracted.txt"
        dep = f"{input_stats_file}.OK"
        tgt = f"{output_stats_file}.OK"
        cmd = f"{extract_general_stats} {input_stats_file} -o {output_stats_file} -s {sample.id}"
        pg.add(tgt, dep, cmd)

        # flag stats
        output_stats_file = f"{stats_dir}/flag/{sample.id}.txt"
        dep = f"{bam_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{output_stats_file}.OK"
        cmd = f"{samtools} flagstat {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # idx stats
        output_stats_file = f"{stats_dir}/idx/{sample.id}.txt"
        dep = f"{bam_dir}/{sample.id}.bam.bai.OK"
        tgt = f"{output_stats_file}.OK"
        cmd = f"{samtools} idxstats {input_bam_file} > {output_stats_file}"
        samtools_multiqc_dep += f" {tgt}"
        pg.add(tgt, dep, cmd)

        # plot samtools stats
        input_stats_file = f"{stats_dir}/general/{sample.id}.txt"
        dep = f"{input_stats_file}.OK"
        tgt = f"{stats_dir}/general/{sample.id}.plot_bamstats.OK"
        cmd = f"{plot_bamstats} -p {stats_dir}/plot {input_stats_file}"
        pg.add(tgt, dep, cmd)


    # plot samtools
    output_dir = f"{plot_dir}"
    log = f"{plot_dir}/samtools.multiqc_report.log"
    err = f"{plot_dir}/samtools.multiqc_report.err"
    dep = samtools_multiqc_dep
    tgt = f"{plot_dir}/samtools.multiqc_report.OK"
    cmd = f"cd {stats_dir}; {multiqc} . -m samtools -f -o {output_dir} -n samtools --no-ansi > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # denovo stacks
    # log = f"{working_dir}/denovo_stacks.log"
    # tgt = f"{working_dir}/denovo_stacks.OK"
    # dep = fastq_files_OK
    # cmd = f"{denovo_stacks} -T 45 -M 7 -o {denovo_stacks_dir} --popmap {population_map_file} --samples {fastq_dir} --paired -X \"ustacks: --force-diff-len\" 2> {log}"
    # pg.add(tgt, dep, cmd)

    #ref stacks
    log = f"{working_dir}/ref_stacks.log"
    tgt = f"{working_dir}/ref_stacks.OK"
    dep = bam_files_OK
    cmd = f"{ref_stacks} -T 30 -o {ref_stacks_dir} --popmap {population_map_file} --samples {bam_dir} 2> {log}"
    pg.add(tgt, dep, cmd)

    #denovo_map.pl -T 45 -M 7 -o ./stacks --popmap ./population.map --samples ./fastq --paired -X "ustacks: --force-diff-len"
    #populations -P . --vcf -O . --no-hap-exports --structure --genepop --vcf-all  --phylip --write-random-snp -t 30

    #filter SNPs

    #compute heterozygosity

    #compute pairwise IBS stats

    #plot heterozygosity plots

    #plot pairwise IBS plots

    #filter individuals

    #plot SNP/Sample call rate at each

    #run structure

    #run PCA

    #plot geospatial plot with structure pie charts

    #plot PCA with structure pie charts





    # clean
    pg.add_clean(f"rm -fr {ref_dir} {denovo_stacks_dir} {fastq_dir}")

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
