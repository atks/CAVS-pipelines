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


@click.command()
@click.option(
    "-m",
    "--make_file",
    show_default=True,
    default="run_cdv_analysis.mk",
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
    Canine Distemper Virus Analysis Pipeline

    e.g. generate_cdv_analysis.py
    """
    print("\t{0:<20} :   {1:<10}".format("make_file", make_file))
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))

    # initialize
    pg = PipelineGenerator(make_file)

    # create directories in destination folder directory
    ref_dir = f"{working_dir}/ref"
    assembly_dir = f"{working_dir}/assembly"
    blast_dir = f"{working_dir}/blast"
    pairwise_alignment_dir = f"{working_dir}/pairwise_alignment"
    alignment_dir = f"{working_dir}/alignment"
    quast_dir = f"{working_dir}/quast"
    annotation_dir = f"{working_dir}/annotation"
    log_dir = f"{working_dir}/log"

    try:
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(assembly_dir, exist_ok=True)
        os.makedirs(blast_dir, exist_ok=True)
        os.makedirs(pairwise_alignment_dir, exist_ok=True)
        os.makedirs(alignment_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
    except OSError as error:
        print(f"Directory cannot be created")

    #########
    #programs
    #########
    spades = "/usr/local/SPAdes-4.0.0/bin/spades.py"
    blastn = "/usr/local/ncbi-blast-2.16.0+/bin/blastn "
    blastdb_nt = "/db/blast/nt/nt"
    blastdb_tx = "/db/blast/nt"
    seqkit = "/usr/local/seqkit-2.10.1/seqkit"
    seqtk = "/usr/local/seqtk-1.4/seqtk"
    needle = "/usr/local/emboss-6.6.0/bin/needle"
    align_and_consense = "/home/atks/programs/CAVS-pipelines/minipipes/align_and_consense.py"
    quast = f"docker run -t -v  `pwd`:`pwd` -w `pwd` fischuu/quast quast.py"
    prokka = "docker run -t -v  `pwd`:`pwd` -w `pwd` staphb/prokka:1.14.6 prokka "
    

    #################
    # reference files
    #################
    id = "AF014953.1"
    output_fasta_file = f"{ref_dir}/{id}.fasta"
    dep = ""
    tgt = f"{output_fasta_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format fasta > {output_fasta_file}"
    pg.add(tgt, dep, cmd)

    id = "AF014953.1"
    output_genbank_file = f"{ref_dir}/{id}.genbank"
    dep = ""
    tgt = f"{output_genbank_file }.OK"
    cmd = f"efetch -db nuccore -id {id} -format genbank > {output_genbank_file}"
    pg.add(tgt, dep, cmd)

    ###############
    # sequence file
    ###############
    samples = []

    name = f"M250740_CDV_lung"
    fastq1 = f"/net/singapura/var/hts/ilm74/74_2_M250740_CDV_lung_R1.fastq.gz"
    fastq2 = f"/net/singapura/var/hts/ilm74/74_2_M250740_CDV_lung_R2.fastq.gz"
    samples.append(Sample(1, name, fastq1, fastq2, "", ""))

    name = f"M250740_CDV_bladder"
    fastq1 = f"/net/singapura/var/hts/ilm74/74_3_M250740_CDV_bladder_R1.fastq.gz"
    fastq2 = f"/net/singapura/var/hts/ilm74/74_3_M250740_CDV_bladder_R2.fastq.gz"
    samples.append(Sample(1, name, fastq1, fastq2, "", ""))

    name = f"M250740_CDV_brain"
    fastq1 = f"/net/singapura/var/hts/ilm74/74_4_M250740_CDV_brain_R1.fastq.gz"
    fastq2 = f"/net/singapura/var/hts/ilm74/74_4_M250740_CDV_brain_R2.fastq.gz"
    samples.append(Sample(1, name, fastq1, fastq2, "", ""))

    ######################
    # Assemble
    ######################
    for sample in samples:
        
        #metaviral assembly
        output_dir = f"{assembly_dir}/{sample.name}/metaviral_assembly"
        dep = ""
        tgt = f"{log_dir}/{sample.name}_metaviral_assembly.OK"
        log = f"{log_dir}/{sample.name}_metaviral_assembly.OK"
        err = f"{log_dir}/{sample.name}_metaviral_assembly.err"
        cmd = f"{spades} -1 {sample.fastq1} -2 {sample.fastq2} --threads 10 --metaviral -o {output_dir} > {log} 2> {err}"
        pg.add_srun(tgt, dep, cmd, 10)

        #rnaviral assembly
        output_dir = f"{assembly_dir}/{sample.name}/rnaviral_assembly"
        dep = ""
        tgt = f"{log_dir}/{sample.name}_rnaviral_assembly.OK"
        log = f"{log_dir}/{sample.name}_rnaviral_assembly.OK"
        err = f"{log_dir}/{sample.name}_rnaviral_assembly.err"
        cmd = f"{spades} -1 {sample.fastq1} -2 {sample.fastq2} --threads 10 --rnaviral -o {output_dir} > {log} 2> {err}"
        pg.add_srun(tgt, dep, cmd, 10)

    ########
    # blast  
    ########
    src_fasta_file = f"{assembly_dir}/M250740_CDV_bladder/metaviral_assembly/contigs.fasta"
    output_txt_file = f"{blast_dir}/blast.results.txt"
    log = f"{log_dir}/blast.log"
    tgt = f"{log_dir}/blast.OK"
    dep = f"{log_dir}/M250740_CDV_bladder_metaviral_assembly.OK"
    cmd = f"export BLASTDB={blastdb_tx}/; {blastn} -db {blastdb_nt} -query {src_fasta_file} -outfmt \"6 qacc sacc qlen slen score length pident stitle staxids sscinames scomnames sskingdoms\" -max_target_seqs 20 -evalue 1e-5 -task megablast -out {output_txt_file} > {log}"
    pg.add(tgt, dep, cmd)


    ####################
    # Pairwise alignment  
    ####################
    src_fasta_file = f"{assembly_dir}/M250740_CDV_bladder/metaviral_assembly/contigs.fasta"
    dst_fasta_file = f"{pairwise_alignment_dir}/metaviral.fasta"
    tgt = f"{log_dir}/metaviral.fasta.OK"
    dep = f"{log_dir}/M250740_CDV_bladder_metaviral_assembly.OK"
    cmd = f"{seqkit} replace {src_fasta_file} -p \"^.*$$\" -r metaviral | {seqtk} seq > {dst_fasta_file}"
    pg.add(tgt, dep, cmd)

    src_fasta_file = f"{assembly_dir}/M250740_CDV_bladder/rnaviral_assembly/contigs.fasta"
    dst_fasta_file = f"{pairwise_alignment_dir}/rnaviral.fasta"
    tgt = f"{log_dir}/rnaviral.fasta.OK"
    dep = f"{log_dir}/M250740_CDV_bladder_rnaviral_assembly.OK"
    cmd = f"echo NODE_1_length_15672_cov_101.763379 | {seqtk} subseq {src_fasta_file} - | {seqkit} replace -p \"^.*$$\"  -r rnaviral | {seqtk} seq -r > {dst_fasta_file}"
    pg.add(tgt, dep, cmd)

    fasta1_file = f"{pairwise_alignment_dir}/metaviral.fasta"
    fasta2_file = f"{pairwise_alignment_dir}/rnaviral.fasta"
    output_txt_file = f"{pairwise_alignment_dir}/metaviral_rnaviral_alignment.txt"
    log = f"{log_dir}/metaviral_rnaviral_needle_alignment.log"
    tgt = f"{log_dir}/metaviral_rnaviral_needle_alignment.OK"
    dep = f"{log_dir}/metaviral.fasta.OK {log_dir}/rnaviral.fasta.OK"
    cmd = f"{needle} {fasta1_file} {fasta2_file} -gapopen 10 -gapextend 0.5 -outfile {output_txt_file} > {log} 2>&1"
    pg.add(tgt, dep, cmd)

    ################
    # Read alignment  
    ################
    fastq1_file = f"/net/singapura/var/hts/ilm74/74_3_M250740_CDV_bladder_R1.fastq.gz"
    fastq2_file = f"/net/singapura/var/hts/ilm74/74_3_M250740_CDV_bladder_R2.fastq.gz"

    #against metaviral assembly
    ref_fasta_file = f"{pairwise_alignment_dir}/metaviral.fasta"
    output_dir = f"{alignment_dir}/metaviral_read_alignment"
    log = f"{log_dir}/metaviral_read_alignment.log"
    tgt = f"{log_dir}/metaviral_read_alignment.OK"
    dep = f"{log_dir}/metaviral_rnaviral_needle_alignment.OK"
    cmd = f"{align_and_consense} -r {ref_fasta_file} -1 {fastq1_file} -2 {fastq2_file} -o {output_dir} > {log} 2>&1"
    pg.add(tgt, dep, cmd)

    #against rnaviral assembly
    ref_fasta_file = f"{pairwise_alignment_dir}/rnaviral.fasta"
    output_dir = f"{alignment_dir}/rnaviral_read_alignment"
    log = f"{log_dir}/rnaviral_read_alignment.log"
    tgt = f"{log_dir}/rnaviral_read_alignment.OK"
    dep = f"{log_dir}/metaviral_rnaviral_needle_alignment.OK"
    cmd = f"{align_and_consense} -r {ref_fasta_file} -1 {fastq1_file} -2 {fastq2_file} -o {output_dir} > {log} 2>&1"
    pg.add(tgt, dep, cmd)

    #########################################
    # Trim and copy chosen assembled sequence  
    #########################################
    #trim last there nucleotides off RNAViral assembly due to ambiguity observed in alignments
    #seqkit subseq pairwise_alignment/rnaviral.fasta -r 1:15669 | seqkit replace -p "^.*$" -r "M250740"  -o new.fasta
    ref_fasta_file = f"{pairwise_alignment_dir}/rnaviral.fasta"
    output_fasta_file = f"{ref_dir}/M250740.fasta"
    log = f"{log_dir}/reference_seq_trimming_renaming.log"
    tgt = f"{log_dir}/M250740.fasta.OK"
    dep = f"{log_dir}/metaviral_read_alignment.OK {log_dir}/rnaviral_read_alignment.OK" 
    cmd = f"{seqkit} subseq {ref_fasta_file} -r 1:15669 | {seqkit} replace -p \"^.*$$\" -r M250740 -o {output_fasta_file} > {log} 2>&1"
    pg.add(tgt, dep, cmd)

    ################################
    # Read alignment against M250740
    ################################
    fastq1_file = f"/net/singapura/var/hts/ilm74/74_3_M250740_CDV_bladder_R1.fastq.gz"
    fastq2_file = f"/net/singapura/var/hts/ilm74/74_3_M250740_CDV_bladder_R2.fastq.gz"

    ref_fasta_file = f"{ref_dir}/M250740.fasta"
    output_dir = f"{alignment_dir}/M250740_read_alignment"
    log = f"{log_dir}/M250740_read_alignment.log"
    tgt = f"{log_dir}/M250740_read_alignment.OK"
    dep = f"{log_dir}/M250740.fasta.OK"
    cmd = f"{align_and_consense} -r {ref_fasta_file} -1 {fastq1_file} -2 {fastq2_file} -o {output_dir} > {log} 2>&1"
    pg.add(tgt, dep, cmd)

    #######
    # Quast
    #######   
    ref_fasta_file = f"{ref_dir}/M250740.fasta"
    input_bam_file = f"{alignment_dir}/M250740_read_alignment/bam/ilm.bam"
    log = f"{log_dir}/M250740_quast.log"
    tgt = f"{log_dir}/M250740_quast.OK"
    dep = f"{log_dir}/M250740_read_alignment.OK"
    cmd = f"{quast} {ref_fasta_file} --bam {input_bam_file} -o {quast_dir} > {log} 2>&1"
    pg.add(tgt, dep, cmd) 
    
    ##########
    # Annotate
    ##########
    fasta_file = f"{ref_dir}/M250740.fasta"
    genbank_file = f"{ref_dir}/AF014953.1.genbank"
    output_dir = f"{annotation_dir}/M250740"
    log = f"{log_dir}/M250740_prokka.log"
    tgt = f"{log_dir}/M250740_prokka.OK"
    dep = f"{log_dir}/M250740.fasta.OK"
    cmd = f"{prokka} --kingdom Viruses {fasta_file} --proteins {genbank_file} --force --outdir {output_dir} --prefix M250740 > {log} 2>&1"
    pg.add(tgt, dep, cmd) 

#docker run  -t -v  `pwd`:`pwd` -w `pwd` staphb/prokka:1.14.6 prokka  --kingdom Viruses /home/atks/analysis/var/20251024_cdv/ref/M250740.fasta  --proteins ref/AF014953.1.genbank  --force --outdir M250740 --prefix M250740
#docker run  -t -v  `pwd`:`pwd` -w `pwd` staphb/prokka:1.14.6 prokka  
#docker run -u \"root:root\" -t -v  `pwd`:`pwd` -w `pwd` staphb/prokka: prokka  --kingdom Viruses ../fasta/contigs.fasta --proteins ../1915_genbank_cdv/NC_001921.1.genbank  --force --outdir M250740 --prefix M250740
#docker run   staphb/prokka:1.14.6 prokka  --kingdom Viruses /home/atks//fasta/contigs.fasta --proteins ../1915_genbank_cdv/NC_001921.1.genbank  --force --outdir M250740 --prefix M250740
#scripts/filter_prokka_tbl.py  M220338/M220338.tbl -r M220338.fasta -o M220338.filtered.tbl

    #######################
    # prepare whole genomes
    #######################



    #####
    # phylogenetics
    #####




    # clean
    pg.add_clean(f"rm -fr {ref_dir} {assembly_dir} {log_dir}")

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
        self.fastq1 = ""
        self.fastq2 = ""
        self.fastq1_OK = ""
        self.fastq2_OK = ""

    def __init__(self, id, name, fastq1, fastq2):
        self.id = id
        self.name = name
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.fastq1_OK = ""
        self.fastq2_OK = ""

    def __init__(self, id, name, fastq1, fastq2, fastq1_OK, fastq2_OK):
        self.id = id
        self.name = name
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.fastq1_OK = fastq1_OK
        self.fastq2_OK = fastq2_OK

    def print(self):
        print(f"id         : {self.id}")
        print(f"name       : {self.name}")
        print(f"fastq1     : {self.fastq1}")
        print(f"fastq2     : {self.fastq2}")
        print(f"fastq1_OK  : {self.fastq1_OK}")
        print(f"fastq2_OK  : {self.fastq2_OK}")


if __name__ == "__main__":
    main()
