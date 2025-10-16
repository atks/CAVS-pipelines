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

import sys
import os
import click
import subprocess
from shutil import copy2

@click.command()
@click.option(
    "-i",
    "--contigs_fasta_file",
    show_default=True,
    help="Contigs fasta filerectory",
)
@click.option(
    "-c",
    "--coverage_txt_file",
    required=True,
    help="Coverage stats file generate using samtools coverage",
)
@click.option(
    "-o",
    "--output_dir",
    default="mitoseq",
    show_default=True,
    help="Output Directory",
)
@click.option(
    "--min_len",
    default=15000,
    help="Minimum Sequence length",
)
@click.option(
    "--max_len",
    default=17000,
    help="Maximum Sequence length",
)
@click.option(
    "--min_cov",
    default=30,
    help="Minimum Coverage",
)
def main(
    contigs_fasta_file,
    coverage_txt_file,
    output_dir,
    min_len,
    max_len,
    min_cov,
):
    """
    Extracts candidate mitochondria sequences and blasts them from a genome skimming sequence run.

    e.g. search_mitoseq.py -i contigs.fasta -c coverage.txt -o blast_mito.txt
    """
    
    print("\t{0:<20} :   {1:<10}".format("contigs FASTA file", contigs_fasta_file))
    print("\t{0:<20} :   {1:<10}".format("coverage TXT file", coverage_txt_file))
    print("\t{0:<20} :   {1:<10}".format("output BLAST TXT file", output_dir))
    print("\t{0:<20} :   {1:<10}".format("mininum contig length", min_len))
    print("\t{0:<20} :   {1:<10}".format("maximum contig length", max_len))
    print("\t{0:<20} :   {1:<10}".format("minimum coverage", min_cov))

    # version
    version = "1.0.0"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/search_mitoseq.log")

    # programs
    seqtk = "/usr/local/seqtk-1.4/seqtk"
    blastn = "/usr/local/ncbi-blast-2.16.0+/bin/blastn"

    # make directories
    output_dir = os.path.abspath(output_dir)
    trace_dir = os.path.join(output_dir, "trace")
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # select candidate sequences
    # read coverage file
    candidate_seq_ids = []
    with open(coverage_txt_file, "r") as file:
        index = 0
        for line in file:
            if not line.startswith("#"):
                rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq = line.rstrip().split("\t")
                
                startpos = int(startpos)
                endpos = int(endpos)
                meandepth = float(meandepth)

                # if endpos > 15000:
                #     print(rname, startpos, endpos, meandepth)

                if endpos >= min_len and endpos <= max_len and meandepth >= min_cov:
                    #print(rname)
                    candidate_seq_ids.append(rname)

    seq_no = len(candidate_seq_ids)
    if seq_no == 0:
        exit("No candidate sequences")      

    #write out 
    id_txt_file = os.path.join(output_dir, "id.txt")
    with open(id_txt_file, "w") as out:
        for id in candidate_seq_ids:
            out.write(f"{id}\n")

    mpm.log(f"{seq_no} selected for blasting")

    #extract candidate seq
    input_fasta_file = contigs_fasta_file
    output_fasta_file = os.path.join(output_dir, f"candidate_mitoseq.fasta")
    cmd = f"{seqtk} subseq {input_fasta_file} {id_txt_file} > {output_fasta_file}"
    tgt = f"{output_fasta_file}.OK"
    desc = f"extract candidate sequences"
    mpm.run(cmd, tgt, desc)

    #blast candidate sequences
    input_fasta_file = os.path.join(output_dir, f"candidate_mitoseq.fasta")
    output_blast_txt_file = os.path.join(output_dir, f"candidate_blast.txt")
    log = os.path.join(output_dir, f"blast.log")
    blastdb = "/db/blast/nt/nt"
    cmd = f"{blastn} -db {blastdb} -query {input_fasta_file} -outfmt \"6 qacc sacc qlen slen score length pident stitle staxids sscinames scomnames sskingdoms\" -max_target_seqs 10 -evalue 1e-5 -out {output_blast_txt_file} > {log}"
    tgt = f"{output_blast_txt_file}.OK"
    desc = f"blast candidate sequences"
    mpm.run(cmd, tgt, desc)






    # write log file
    mpm.print_log()

    # copy files to trace
    copy2(__file__, trace_dir)
    subprocess.run(f'echo {" ".join(sys.argv)} > {trace_dir}/cmd.txt', shell=True, check=True)

class MiniPipeManager(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log_msg = []

    def run(self, cmd, tgt, desc):
        try:
            if os.path.exists(tgt):
                self.log(f"{desc} -  already executed")
                self.log(cmd)
                return
            else:
                self.log(f"{desc}")
                self.log(cmd)
                subprocess.run(cmd, shell=True, check=True)
                subprocess.run(f"touch {tgt}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            self.log(f" - failed")
            exit(1)

    def log(self, msg):
        print(msg)
        self.log_msg.append(msg)

    def print_log(self):
        self.log(f"\nlogs written to {self.log_file}")
        with open(self.log_file, "w") as f:
            f.write("\n".join(self.log_msg))


if __name__ == "__main__":
    main() # type: ignore
