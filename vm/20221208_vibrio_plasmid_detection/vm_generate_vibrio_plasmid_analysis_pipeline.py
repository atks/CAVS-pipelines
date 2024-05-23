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
    default="vibrio_plasmid_detection.mk",
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
    Look for pirAB homologs on isolate sequences

    e.g. vm_generate_vibrio_plasmid_analysis_pipeline.py
    """
    log_dir = f"{working_dir}/log"
    ref_dir = f"{working_dir}/ref"
    scripts_dir = f"{working_dir}/scripts"
    blast_dir = f"{working_dir}/blast"
    vibrio_contigs_db_dir = f"{blast_dir}/vibrio_contigs_db"
    bwa_dir = f"{working_dir}/bwa"
    bwa_db_dir = f"{bwa_dir}/db"

    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("sample_file", sample_file))

    # Reference files
    # GenBank (accession no. KM067908)
    # PirA cds position: 17198..17533
    # PirB cds position: 17546..18862

    # read sample file
    samples = []
    with open(sample_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                index += 1
                sample_id, read1_fastq, read2_fastq, contig_fasta = line.rstrip().split(
                    "\t"
                )
                species = ""
                samples.append(
                    Sample(
                        index,
                        sample_id,
                        species,
                        read1_fastq,
                        read2_fastq,
                        contig_fasta,
                    )
                )

    # create directories in destination folder directory
    try:
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(vibrio_contigs_db_dir, exist_ok=True)
        os.makedirs(blast_dir, exist_ok=True)
        os.makedirs(bwa_db_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory {new_dir} cannot be created")

    # initialize
    pg = PipelineGenerator(make_file)

    # download sequences
    output_fasta = f"{ref_dir}/KM067908.fasta"
    tgt = f"{log_dir}/KM067908.fasta.OK"
    dep = ""
    cmd = f"efetch -db nuccore -id KM067908 -format fasta > {output_fasta}"
    pg.add(tgt, dep, cmd)

    output_gb = f"{ref_dir}/KM067908.gb"
    tgt = f"{log_dir}/KM067908.gb.OK"
    dep = ""
    cmd = f"efetch -db nuccore -id KM067908 -format genbank > {output_gb}"
    pg.add(tgt, dep, cmd)

    # annotate species
    annotate_isolate_species = f"{scripts_dir}/annotate_isolate_species.py"
    input_sa = f"{sample_file}"
    output_sa = f"{ref_dir}/vibrio.annotated.sa"
    tgt = f"{log_dir}/vibrio.annotated.sa.OK"
    dep = ""
    cmd = f"{annotate_isolate_species} -s {input_sa} -o {output_sa}"
    pg.add(tgt, dep, cmd)

    # rename contigs to include sample ID and species annotation
    rename_contigs = f"{scripts_dir}/rename_contigs.py"
    input_sa = f"{ref_dir}/vibrio.annotated.sa"
    output_fasta = f"{ref_dir}/vibrio.contigs.fasta"
    log = f"{log_dir}/vibrio.contigs.fasta.log"
    tgt = f"{log_dir}/vibrio.contigs.fasta.OK"
    dep = f"{log_dir}/vibrio.annotated.sa.OK"
    cmd = f"{rename_contigs} -s {input_sa} -o {output_fasta} > {log}"
    pg.add(tgt, dep, cmd)

    # construct DB from vibrio contigs
    makeblastdb = "/usr/local/ncbi-blast-2.13.0+/bin/makeblastdb"
    blastdb_dir = vibrio_contigs_db_dir
    input_fasta = f"{ref_dir}/vibrio.contigs.fasta"
    log = f"{log_dir}/vibrio_contigs_db.log"
    err = f"{log_dir}/vibrio_contigs_db.err"
    tgt = f"{log_dir}/vibrio_contigs_db.OK"
    dep = f"{log_dir}/vibrio.contigs.fasta.OK"
    cmd = f"cd {blastdb_dir}; {makeblastdb} -in {input_fasta} -out vibrio_contigs_db -dbtype nucl -parse_seqids > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # make bed file
    # GenBank (accession no. KM067908)
    # PirA cds position: 17198..17533
    # PirB cds position: 17546..18862
    bed_line = "KM067908.1 17197 17533-KM067908.1 17545 18862"
    output_bed = f"{ref_dir}/pirab.bed"
    tgt = f"{log_dir}/pirab.bed.OK"
    dep = f"{log_dir}/KM067908.fasta.OK"
    cmd = f'echo {bed_line} | tr " " "\\t" | tr "-" "\\n" > {output_bed}'
    pg.add(tgt, dep, cmd)

    # prepare query fasta
    seqtk = "/usr/local/seqtk-1.3/seqtk"
    input_fasta = f"{ref_dir}/KM067908.fasta"
    input_bed = f"{ref_dir}/pirab.bed"
    output_fasta = f"{ref_dir}/pirab.fasta"
    tgt = f"{log_dir}/pirab.fasta.OK"
    dep = f"{log_dir}/pirab.bed.OK"
    cmd = f"seqtk subseq {input_fasta} {input_bed} > {output_fasta}"
    pg.add(tgt, dep, cmd)

    # blastn
    blastn = "/usr/local/ncbi-blast-2.13.0+/bin/blastn"
    blastdb_dir = vibrio_contigs_db_dir
    input_fasta = f"{ref_dir}/pirab.fasta"
    output_txt = f"{blast_dir}/pirab.vibrio.blastn.txt"
    tgt = f"{log_dir}/pirab.vibrio.blastn.txt.OK"
    dep = f"{log_dir}/pirab.fasta.OK {log_dir}/vibrio_contigs_db.OK"
    cmd = f"cd {blastdb_dir}; {blastn} -query {input_fasta} -db vibrio_contigs_db -out {output_txt}"
    pg.add(tgt, dep, cmd)

    # convert nucleotides to amini acids
    nt_to_aa = f"{scripts_dir}/nt_to_aa.py"
    input_fasta = f"{ref_dir}/pirab.fasta"
    output_fasta = f"{ref_dir}/pirab.aa.fasta"
    tgt = f"{log_dir}/pirab.aa.fasta.OK"
    dep = f"{log_dir}/pirab.fasta.OK"
    cmd = f"{nt_to_aa} -i {input_fasta} -o {output_fasta}"
    pg.add(tgt, dep, cmd)

    # tblastn
    tblastn = "/usr/local/ncbi-blast-2.13.0+/bin/tblastn"
    blastdb_dir = vibrio_contigs_db_dir
    input_fasta = f"{ref_dir}/pirab.aa.fasta"
    output_txt = f"{blast_dir}/pirab.vibrio.tblastn.txt"
    tgt = f"{log_dir}/pirab.vibrio.tblastn.txt.OK"
    dep = f"{log_dir}/pirab.aa.fasta.OK {log_dir}/vibrio_contigs_db.OK"
    cmd = f"cd {blastdb_dir}; {tblastn} -query {input_fasta} -db vibrio_contigs_db -out {output_txt}"
    pg.add(tgt, dep, cmd)

    # build bwa index
    bwa = "/usr/local/bwa-0.7.17/bwa"
    samtools = "/usr/local/samtools-1.16/bin/samtools"
    input_fasta = f"{bwa_db_dir}/pirab.fasta"
    log = f"{log_dir}/pirab_bwa_db.log"
    err = f"{log_dir}/pirab_bwa_db.err"
    tgt = f"{log_dir}/pirab_bwa_db.OK"
    dep = f"{log_dir}/pirab.fasta.OK"
    cmd = f"cp {ref_dir}/pirab.fasta {bwa_db_dir}; {bwa} index -a bwtsw {input_fasta} > {log} 2> {err}"
    pg.add(tgt, dep, cmd)

    # bwa map
    for idx, sample in enumerate(samples):
        input_fasta = f"{bwa_db_dir}/pirab.fasta"
        output_bam = f"{bwa_dir}/{sample.id}.bam"
        log = f"{log_dir}/bwa.{sample.id}.log"
        err = f"{log_dir}/bwa.{sample.id}.err"
        tgt = f"{log_dir}/{sample.id}.bam.OK"
        dep = f"{log_dir}/pirab_bwa_db.OK"
        cmd = f"{bwa} mem -t 2 -M {input_fasta} {sample.fastq1} {sample.fastq2} 2> {err} | {samtools} view -F 4 -o {output_bam}  > {log} "
        pg.add(tgt, dep, cmd)

    # to check output from bam
    # cd {working_dir}/bwa
    # ls *.bam | perl -lane '{print "samtools coverage $_ | cut -f4"}' | sh | sort  | uniq -c
    # 130 0
    #  65 numreads

    # clean
    cmd = f"rm -fr {log_dir} {blast_dir} {bwa_dir}"
    pg.add_clean(cmd)

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


def append_file_suffix(file_name, suffix):
    file_name = os.path.basename(file_name)
    tokens = file_name.split(".")
    tokens.insert(-1, suffix)
    return ".".join(tokens)


if __name__ == "__main__":
    main()
