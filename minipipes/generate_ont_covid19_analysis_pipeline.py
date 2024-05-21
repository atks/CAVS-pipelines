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

# Example run:
# ./generate_illumina_virus_detection_pipeline.py -o /home/melody/Out -s /home/melody/ilm30.sa
# ./generate_illumina_virus_detection_pipeline.py -v -o /home/melody/Out_test -s /home/melody/test.sa
# ./generate_illumina_virus_detection_pipeline.py -s /home/melody/testdir/test2.sa


def main():
    user_directory = os.getcwd()
    home_directory = os.path.expanduser("~")
    multiqc = home_directory + "/.local/bin/multiqc"
    fastqc = "/usr/local/FastQC-0.11.9/fastqc"
    uniprot_database = home_directory + "/db/viral_proteins/uniref90.fasta.gz"
    # uniprot_database = home_directory + "/out.fa"
    diamond = "/usr/local/diamond-2.0.11/diamond"
    seqtk = "/usr/local/seqtk-1.3/seqtk"
    trinity = "/usr/local/trinityrnaseq-v2.13.1/Trinity"
    refseq = home_directory + "/db/refseq/blastdb/refseq.virus.fasta"

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--output_directory",
        help="Takes in path of directory that you want to create to store the output",
        type=str,
    )
    parser.add_argument(
        "-v", "--version", help="Tells you the version of script", action="store_true"
    )
    parser.add_argument(
        "-s",
        "--sample_file",
        help="Takes in the data that you wish to use to run FastQC on",
        type=str,
    )
    args = parser.parse_args()

    # Reading the file
    try:
        data_file = open(args.sample_file, "r")
        samples = data_file.readlines()
        data_file.close()
    except:
        print("You did not input a valid argument for -s, please see -h for more info.")
        return None

    if args.output_directory != None:
        if args.output_directory[0] != "/":
            args.output_directory = user_directory + "/" + args.output_directory
        print(
            f"Your directory is {args.output_directory}, your sample file is {args.sample_file}."
        )
    else:
        print(
            f"Default directory is QC_Out which has been created in your current directory, sample file is {args.sample_file}"
        )
        args.output_directory = user_directory + "/QC_out"

    if args.version:
        print("version 1.0")

    # Create the directory
    # 'Out' in
    # '/home/melody'
    if os.path.exists(args.output_directory):
        print(f"Warning: {args.output_directory} is an existing directory.")
    else:
        os.mkdir(args.output_directory)
        print(f"Directory {args.output_directory} created")

    # creates nested list for each sample
    sample_list = []
    for sample in samples:
        nested_sample_ls = []
        ls_samples = sample.split("\t")
        for sam in ls_samples:
            clean_sam = sam.rstrip("\n")
            nested_sample_ls.append(clean_sam)
        sample_list.append(nested_sample_ls)
        print(sample_list)

    # Creates a dictionary of samples
    d1 = {}
    for nested_sample_ls in sample_list:
        sample_name = nested_sample_ls[0]
        sample_R1 = nested_sample_ls[1]
        sample_R2 = nested_sample_ls[2]
        s = Sample(sample_name, sample_R1, sample_R2)
        d1[s.name] = s
        print(d1)

    # change directory
    os.chdir(f"{args.output_directory}")
    cwd = os.getcwd()

    # creates directories for each sample
    for key in d1:
        sample_path = os.path.join(cwd, key)
        os.mkdir(sample_path)
        print(f"Directory {key} created")

    # write makefile to user_directory
    os.chdir(f"{user_directory}")
    cwd = os.getcwd()
    print(f"Makefile will be created in {cwd}")

    WRITE_FLAG = "w"
    is_exist = os.path.isfile("./illumina_virus_detection_pipeline.mk")
    if not is_exist:
        # Create a Makefile in the same location
        print("Created makefile")
        f = open(f"./illumina_virus_detection_pipeline.mk", WRITE_FLAG)

        f.close()

    else:
        print("Make file exists, overwritten")

    f = open(f"./illumina_virus_detection_pipeline.mk", WRITE_FLAG)

    f.write("all:\t")

    f.write(f"{args.output_directory}/generate_db.OK ")

    for key in d1:
        f.write(f"{args.output_directory}/{key}/fastqc1.OK ")
        f.write(f"{args.output_directory}/{key}/diamond_fastq1.OK ")
        f.write(f"{args.output_directory}/{key}/subseq1.OK ")

        f.write(f"{args.output_directory}/{key}/fastqc2.OK ")
        f.write(f"{args.output_directory}/{key}/diamond_fastq2.OK ")
        f.write(f"{args.output_directory}/{key}/subseq2.OK ")

        f.write(f"{args.output_directory}/{key}/trinity.OK ")
        f.write(f"{args.output_directory}/{key}/refseqblast.OK ")

    f.write(f"{args.output_directory}/multiqc_output.OK")

    f.write("\n\n")

    # generates diamond reference database for indexing
    f.write(f"{args.output_directory}/generate_db.OK:\n")
    f.write(
        f"\t{diamond} makedb --in {uniprot_database} -d {args.output_directory}/nr > {args.output_directory}/nr.log 2> {args.output_directory}/nr.err\n"
    )
    f.write(f"\ttouch {args.output_directory}/generate_db.OK\n\n")

    for key in d1:
        sample_object = d1[key]

        # run fastqc on fastq1
        #####
        f.write(f"{args.output_directory}/{key}/fastqc1.OK:\n")
        f.write(
            f"\t{fastqc} -o {args.output_directory}/{key} --extract {sample_object.fastq1} > {args.output_directory}/{key}/fastqc1.log 2> {args.output_directory}/{key}/fastqc1.err\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/fastqc1.OK\n\n")

        # run diamond on fastq1
        #####
        f.write(
            f"{args.output_directory}/{key}/diamond_fastq1.OK: {args.output_directory}/generate_db.OK\n"
        )
        f.write(
            f"\tzcat {sample_object.fastq1} | seqtk seq -A > {args.output_directory}/{key}/input{key}1.fasta\n"
        )
        f.write(
            f"\twc -l {args.output_directory}/{key}/input{key}2.fasta > count1.txt\n"
        )
        f.write(
            f"\t{diamond} blastx -d {args.output_directory}/nr.dmnd -q {args.output_directory}/{key}/input{key}1.fasta -o {args.output_directory}/{key}/stitle{key}1.m8 -f 6 qseqid stitle 2> {args.output_directory}/{key}/diamond1.log 2> {args.output_directory}/{key}/diamond1.err\n"
        )
        f.write(
            f"\twc -l {args.output_directory}/{key}/stitle{key}1.m8 > {args.output_directory}/{key}/length_{key}1_diamond.txt\n"
        )
        f.write(
            f"\tcat {args.output_directory}/{key}/stitle{key}1.m8 | cut -f2 | sort | uniq -c > {args.output_directory}/{key}/grouped1.txt\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/diamond_fastq1.OK\n\n")

        # subseq fastq1
        #####
        f.write(
            f"{args.output_directory}/{key}/subseq1.OK: {args.output_directory}/{key}/diamond_fastq1.OK\n"
        )
        f.write(
            f"\t{seqtk} subseq {sample_object.fastq1} {args.output_directory}/{key}/stitle{key}1.m8 > {args.output_directory}/{key}/{key}_1.fq\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/subseq1.OK\n\n")

        # run fastqc on fastq2
        #####
        f.write(f"{args.output_directory}/{key}/fastqc2.OK:\n")
        f.write(
            f"\t{fastqc} -o {args.output_directory}/{key} --extract {sample_object.fastq2} > {args.output_directory}/{key}/fastqc2.log 2> {args.output_directory}/{key}/fastqc2.err\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/fastqc2.OK\n\n")

        # run diamond on fastq2
        #####
        f.write(
            f"{args.output_directory}/{key}/diamond_fastq2.OK: {args.output_directory}/generate_db.OK\n"
        )
        f.write(
            f"\tzcat {sample_object.fastq2} | seqtk seq -A > {args.output_directory}/{key}/input{key}2.fasta\n"
        )
        f.write(
            f"\twc -l {args.output_directory}/{key}/input{key}2.fasta > count2.txt\n"
        )
        f.write(
            f"\t{diamond} blastx -d {args.output_directory}/nr.dmnd -q {args.output_directory}/{key}/input{key}2.fasta -o {args.output_directory}/{key}/stitle{key}2.m8 -f 6 qseqid stitle 2> {args.output_directory}/{key}/diamond2.log 2> {args.output_directory}/{key}/diamond2.err\n"
        )
        f.write(
            f"\twc -l {args.output_directory}/{key}/stitle{key}2.m8 > {args.output_directory}/{key}/length_{key}2_diamond.txt\n"
        )
        f.write(
            f"\tcat {args.output_directory}/{key}/stitle{key}2.m8 | cut -f2 | sort | uniq -c > {args.output_directory}/{key}/grouped2.txt\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/diamond_fastq2.OK\n\n")

        # subseq fastq2
        #####
        f.write(
            f"{args.output_directory}/{key}/subseq2.OK: {args.output_directory}/{key}/diamond_fastq2.OK\n"
        )
        f.write(
            f"\t{seqtk} subseq {sample_object.fastq2} {args.output_directory}/{key}/stitle{key}2.m8 > {args.output_directory}/{key}/{key}_2.fq\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/subseq2.OK\n\n")

        # trinity
        #####
        f.write(
            f"{args.output_directory}/{key}/trinity.OK: {args.output_directory}/{key}/subseq1.OK {args.output_directory}/{key}/subseq2.OK\n"
        )
        f.write(
            f"\t{trinity} --seqType fq --left {args.output_directory}/{key}/{key}_1.fq --right {args.output_directory}/{key}/{key}_2.fq --CPU 2 --max_memory 40G --output {args.output_directory}/{key}/trinity_assembly > {args.output_directory}/{key}/trinity_assembly.txt.log 2> {args.output_directory}/{key}/trinity_assembly.txt.err\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/trinity.OK\n\n")

        # blastn against refseq
        # {args.output_directory}/{key}/trinity_assembly
        #####
        f.write(
            f"{args.output_directory}/{key}/refseqblast.OK: {args.output_directory}/{key}/trinity.OK\n"
        )
        f.write(
            f"\tblastn -db {refseq} -query {args.output_directory}/{key}/trinity_assembly.Trinity.fasta  -outfmt "
            + '"6 stitle pident"'
            + f" -out {args.output_directory}/{key}/{key}_blast.psl -max_target_seqs 1 2> {args.output_directory}/{key}/blast_results{key}.txt.log > {args.output_directory}/{key}/blast_results{key}.txt.err\n"
        )
        f.write(
            f"\tcut -f1 {args.output_directory}/{key}/{key}_blast.psl | sort | uniq -c > {args.output_directory}/{key}/blast_{key}.txt\n"
        )
        f.write(f"\ttouch {args.output_directory}/{key}/refseqblast.OK\n\n")

    # multiqc
    f.write(f"{args.output_directory}/multiqc_output.OK:\n")
    f.write(
        f"\t{multiqc} {args.output_directory} -o {args.output_directory}/multiqc_output/\n"
    )
    f.write(f"\ttouch {args.output_directory}/multiqc_output.OK\n\n")

    f.close()


# creates sample class
class Sample(object):
    def __init__(self, name, fastq1, fastq2):
        self.name = name
        self.fastq1 = fastq1
        self.fastq2 = fastq2


main()
