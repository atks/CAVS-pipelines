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

import sys
import os
import click
import subprocess
import re
from shutil import copy2


@click.command()
@click.option(
    "-w",
    "--working_dir",
    default=f"{os.getcwd()}/phylo_tree",
    show_default=True,
    help="working directory",
)
@click.option(
    "-f",
    "--fasta_file",
    required=False,
    show_default=True,
    help="sequence panel fasta file",
)
@click.option(
    "-r",
    "--ref_fasta_file",
    required=False,
    show_default=True,
    help="reference fasta file",
)
@click.option(
    "-m",
    "--ref_msa_file",
    required=False,
    show_default=True,
    help="reference multiple sequence alignment file",
)
@click.option(
    "-s",
    "--sample_file",
    required=False,
    show_default=True,
    help="for naming the nodes of the tree",
)
@click.option(
    "-p",
    "--prefix",
    default="genus_species",
    help="for RAXML file naming",
)
def main(working_dir, fasta_file, ref_fasta_file, ref_msa_file, sample_file, prefix):
    """
    Generates a phylogenetic tree from a panel of sequences and reference sequences

     #combines both fasta files, multiple align, build tree \n
         make_phylogenetic_tree.py -f lsdv.fasta -r lsdv_ref.fasta

     #multiple align reference fasta file, build tree \n
         make_phylogenetic_tree.py -r lsdv_ref.fasta

     #multiple align sequences in fasta file to existing multiple alignment, build tree \n
         make_phylogenetic_tree.py -f lsdv.fasta -m lsdv_ref.msa

     #build tree from the multiple sequence alignment \n
         make_phylogenetic_tree.py -m lsdv_ref.msa
    """

    #make sure at least one of the two files is provided
    if (ref_fasta_file is None and ref_msa_file is None) or (ref_fasta_file is not None and ref_msa_file is not None):
        print("Please provide either a reference fasta file or a reference multiple sequence alignment file")
        sys.exit(1)

    combined_fasta_file = True
    multiple_align = "all"

    if fasta_file is not None and ref_fasta_file is None:
        combined_fasta_file = False

    if ref_msa_file is None:
        multiple_align = "all"
    elif fasta_file is not None and ref_msa_file is not None:
        multiple_align = "add"
    else:
        multiple_align = "none"

    print(f"msa status = {multiple_align}")

    # make directories
    output_dir = f"{os.getcwd()}/{prefix}_phylo"
    if working_dir != "":
        output_dir = os.path.abspath(working_dir)
        trace_dir = f"{output_dir}/trace"
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(trace_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    # version
    version = "1.0.0"

    # initialize
    mpm = MiniPipeManager(f"{output_dir}/make_phylogenetic_tree.log")

    # programs
    mafft = "/usr/local/mafft-7.490/bin/mafft"
    raxml = "/usr/local/raxml-ng-1.1.0/raxml-ng"
    seqkit = "/usr/local/seqkit-2.1.0/bin/seqkit"
    gotree = "/usr/local/gotree-0.4.5/gotree"

    # on uluaotu
    # mafft = "/usr/local/mafft-7.525/bin/mafft"
    # raxml = "/usr/local/raxml-ng-1.2.2/raxml-ng"
    # seqkit = "/usr/local/seqkit-2.8.2/seqkit"
    # gotree = "/usr/local/gotree-0.4.5/gotree"

    # log text
    log_text = ""

    # if fasta_file not empty, compare sequences to reference
    desc = (
        f"Generating combined FASTA file with clean IDs for multiple sequence alignment"
    )
    if fasta_file is None:
        fasta_file = ""
        desc = f"Generating FASTA file with clean IDs from reference FASTA only for multiple sequence alignment"
    if ref_fasta_file is None:
        ref_fasta_file = ""
        desc = f"Generating FASTA file with clean IDs from reference multiple sequence alignment only for multiple sequence alignment"
    combined_fasta_file = f"{output_dir}/combined.fasta"
    cmd = fr'cat {fasta_file} {ref_fasta_file} | {seqkit} replace -p "[\s;:,\(\)\']" -r "_"  > {combined_fasta_file}'
    tgt = f"{combined_fasta_file}.OK"
    mpm.run(cmd, tgt, desc)

    if multiple_align == "all":
        # perform multiple sequence alignment
        input_fasta_file = f"{output_dir}/combined.fasta"
        output_msa_file = f"{output_dir}/{prefix}.msa"
        log = f"{output_dir}/msa.log"
        cmd = f"{mafft} --thread -1 {input_fasta_file} > {output_msa_file} 2>{log}"
        tgt = f"{output_msa_file}.OK"
        desc = f"Multiple sequence alignment"
        mpm.run(cmd, tgt, desc)
    elif multiple_align == "add":
        # perform add on to multiple sequence alignment
        input_fasta_file = f"{combined_fasta_file}"
        output_msa_file = f"{output_dir}/{prefix}.msa"
        log = f"{output_dir}/msa.log"
        cmd = f"{mafft} --add {input_fasta_file} {ref_msa_file} > {output_msa_file} 2>{log}"
        tgt = f"{output_msa_file}.OK"
        desc = f"Additive multiple sequence alignment"
        mpm.run(cmd, tgt, desc)
    else:
        pass

    # construct phylogenetic tree with bootstrap
    input_msa_file = f"{output_dir}/{prefix}.msa"
    log = f"{output_dir}/construct_trees.log"
    cmd = f"cd {output_dir}; {raxml} --threads 10 --msa {input_msa_file} --model GTR+G --redo --prefix {prefix} --bootstrap > {log}"
    tgt = f"{output_dir}/construct_trees.OK"
    desc = f"Constructing phylogenetic tree"
    mpm.run(cmd, tgt, desc)

    # construct consensus tree
    log = f"{output_dir}/consensus_tree.log"
    cmd = f"cd {output_dir}; {raxml} --consense MRE --tree {prefix}.raxml.bootstraps --redo --prefix {prefix} > {log}"
    tgt = f"{output_dir}/consensus_tree.OK"
    desc = f"Constructing consensus tree"
    mpm.run(cmd, tgt, desc)

    if sample_file is not None:
        #prepare renaming file
        fasta_hdr_idx = 0
        tree_node_name_idx = 0
        rename_file = f"{output_dir}/rename_tree.txt"
        with open(rename_file, "w") as out:
            with open(sample_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        header_names = line.lstrip('#').strip().split("\t")
                        fasta_hdr_idx = header_names.index("fasta_header")
                        tree_node_name_idx = header_names.index("tree_node_name")
                        out.write(f"#old_name\tnew_name\n")
                    else:
                        values = line.strip().split("\t")
                        old_name = re.sub(r"[\s;:,()']", "_", values[fasta_hdr_idx].lstrip(">"))
                        out.write(f"{old_name}\t{values[tree_node_name_idx]}\n")

        # rename bootstrap trees
        input_tree_file = f"{output_dir}/{prefix}.raxml.bootstraps"
        output_tree_file = f"{output_dir}/{prefix}.raxml.bootstraps.renamed.tree"
        log = f"{output_dir}/rename_bootstrap_trees.log"
        cmd = f"{gotree} rename -i {input_tree_file} -o {output_tree_file} -m {rename_file} > {log}"
        tgt = f"{output_dir}/rename_bootstrap_trees.OK"
        desc = f"Renaming bootstrap trees"
        mpm.run(cmd, tgt, desc)

        # rename consensus tree
        input_tree_file = f"{output_dir}/{prefix}.raxml.consensusTreeMRE"
        output_tree_file = f"{output_dir}/{prefix}.raxml.consensus.renamed.tree"
        log = f"{output_dir}/rename_consensus_tree.log"
        cmd = f"{gotree} rename -i {input_tree_file} -o {output_tree_file} -m {rename_file} > {log}"
        tgt = f"{output_dir}/rename_consensus_tree.OK"
        desc = f"Renaming consensus tree"
        mpm.run(cmd, tgt, desc)

    # copy files to trace
    copy2(__file__, trace_dir)

    # write log file
    mpm.print_log()

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
    main() # type: ignore[arg-type]
