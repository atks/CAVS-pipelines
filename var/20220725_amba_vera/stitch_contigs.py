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

from cmath import pi
import os
import click
import sys
import pysam

lion = "vera"


@click.command()
@click.option(
    "-r",
    "--reference_fasta_file",
    default="EPI_ISL_6600690.fasta",
    show_default=True,
    help="working directory",
)
@click.option(
    "-cr",
    "--contigs_ref_aln_file",
    default=f"{lion}_ref_alignment.txt",
    show_default=True,
    help="working directory",
)
@click.option(
    "-c",
    "--contigs_file",
    default=f"{lion}_contigs.fasta",
    show_default=True,
    help="working directory",
)
@click.option(
    "-b",
    "--bam_file",
    default=f"{lion}.bam",
    show_default=True,
    help="working directory",
)
@click.option(
    "-s",
    "--stitched_fasta_file",
    default=f"{lion}_stitched.fasta",
    show_default=True,
    help="working directory",
)
@click.option(
    "-g",
    "--gaps_contigs_fasta_file",
    default=f"{lion}_gaps_contigs.fasta",
    show_default=True,
    help="working directory",
)
@click.option(
    "-m", "--min_base_quality", default=10, show_default=True, help="working directory"
)
@click.option(
    "-p",
    "--proportion_consensus_cutoff",
    default=0.8,
    show_default=True,
    help="working directory",
)
def main(
    reference_fasta_file,
    contigs_ref_aln_file,
    contigs_file,
    bam_file,
    stitched_fasta_file,
    gaps_contigs_fasta_file,
    min_base_quality,
    proportion_consensus_cutoff,
):
    """
    Stitch Contigs

    e.g.
    """
    print("\t{0:<20} :   {1:<10}".format("reference fasta file", reference_fasta_file))
    print("\t{0:<20} :   {1:<10}".format("contigs alignment", contigs_ref_aln_file))
    print("\t{0:<20} :   {1:<10}".format("contigs", contigs_file))
    print("\t{0:<20} :   {1:<10}".format("bam", bam_file))
    print("\t{0:<20} :   {1:<10}".format("stitched fasta file", stitched_fasta_file))
    print(
        "\t{0:<20} :   {1:<10}".format(
            "gaps and contigs fasta file", gaps_contigs_fasta_file
        )
    )
    print("\t{0:<20} :   {1:<10}".format("minmum base quality score", min_base_quality))
    print(
        "\t{0:<20} :   {1:<10}".format(
            "proportion consensus cutoff", proportion_consensus_cutoff
        )
    )

    print("=================")
    print("reading reference")
    print("=================")
    REF = []
    name = ""
    with open(reference_fasta_file, "r") as file:
        index = 0
        for line in file:
            line = line.rstrip()
            if line.startswith(">"):
                name = line.lstrip(">")
            else:
                for c in line:
                    REF.append(c)
    CHROM = Chromosome(name, len(REF))
    print(f"Reference length = {len(REF)}")

    print("==========================")
    print("reading contigs alignments")
    print("==========================")
    # qacc sacc qlen slen sstart send score length pident stitle
    # qacc - Severe_acute_respiratory_syndrome-related_coronavirus_contig_1
    # sacc - hCoV-19/Singapore/8152/2021|EPI_ISL_6600690|2021-11-07
    # qlen - 1313
    # slen - 29858
    # sstart - 47
    # send - 1359
    # score - 1313
    # length - 1313
    # pidebt - 100.000
    # stitle - hCoV-19/Singapore/8152/2021|EPI_ISL_6600690|2021-11-07
    ALIGNMENTS = []
    CONTIG2ALIGNMENT = {}
    with open(contigs_ref_aln_file, "r") as file:
        index = 0
        for line in file:
            (
                qacc,
                sacc,
                qlen,
                slen,
                sstart,
                send,
                score,
                length,
                pident,
                stitle,
            ) = line.rstrip().split("\t")
            ALIGNMENTS.append(
                Contig(
                    qacc, sacc, qlen, slen, sstart, send, score, length, pident, stitle
                )
            )
            CONTIG2ALIGNMENT[qacc] = len(ALIGNMENTS) - 1
            ALIGNMENTS[-1].print()

    print("===============")
    print("reading contigs")
    print("===============")
    CONSENSUS = []
    with open(contigs_file, "r") as file:
        seq = ""
        contig_index = -1
        for line in file:
            line = line.rstrip()
            if line.startswith(">"):
                if contig_index != -1:
                    ALIGNMENTS[contig_index].set_seq(seq)
                    seq = ""
                line = line.lstrip(">")
                # >_Severe_acute_respiratory_syndrome-related_coronavirus_contig_1
                # Depth-of-coverage:2603.59
                # Taxon:694009
                # refseq:NC_045512.2
                qacc, depth, taxon, ref = line.split(" ")
                contig_index = CONTIG2ALIGNMENT[qacc]
            else:
                seq = seq + line
        ALIGNMENTS[contig_index].set_seq(seq)

    print("========")
    print("add gaps")
    print("========")
    SEGMENTS = []
    last_end = 0
    gap_no = 1
    for contig in ALIGNMENTS:
        if contig.sstart - 1 != last_end:
            SEGMENTS.append(Gap(f"gap_{gap_no}", last_end + 1, contig.sstart - 1))
            gap_no += 1
        SEGMENTS.append(contig)
        last_end = contig.send
    if last_end != CHROM.len:
        SEGMENTS.append(Gap(f"gap_{gap_no}", last_end + 1, CHROM.len))

    for segment in SEGMENTS:
        segment.print()

    print("===========")
    print("reading bam")
    print("===========")
    sam_file = pysam.AlignmentFile(bam_file, "rb")
    # for each gap
    t_no_cons = 0
    t_no_nons = 0
    t_no_total = 0
    t_depth = 0
    for gap in SEGMENTS:
        if type(gap).__name__ == "Gap":
            gap.print()
            consensus = ""
            no_dels = 0
            no_cons = 0
            no_nons = 0
            no_total = 0
            depth = 0
            for pileupcolumn in sam_file.pileup(CHROM.name, gap.sstart - 1, gap.send):
                if (
                    gap.sstart - 1 <= pileupcolumn.pos
                    and pileupcolumn.pos <= gap.send - 1
                ):
                    # print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                    pileupcolumn.set_min_base_quality(min_base_quality)
                    # print(pileupcolumn.get_mapping_qualities())
                    n = 0
                    B = {}
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del:
                            base = pileupread.alignment.query_sequence[
                                pileupread.query_position
                            ]
                            if base not in B:
                                B[base] = 1
                            else:
                                B[base] += 1
                            n += 1
                        else:
                            if "D" not in B:
                                B["D"] = 1
                            else:
                                B["D"] += 1
                            n += 1
                        depth += 1
                    # decision making
                    print(f"{pileupcolumn.pos+1} ", end="")

                    max_b = ""
                    max_n = 0
                    tie = False
                    for b in B:
                        if B[b] > max_n:
                            max_b = b
                            max_n = B[b]
                            tie = False
                        elif B[b] == max_n:
                            # this case does not exist for amba/vera
                            tie = True
                        print(f":{B[b]}{b}", end="")
                    if tie:
                        print(f" => TIE N")
                        consensus += "N"
                        no_nons += 1
                        no_total += 1
                    else:
                        print(f" / {n}", end="")
                        # choose consensus
                        if n != 0:
                            p = max_n / n
                            con_b = ""
                            if p >= proportion_consensus_cutoff:
                                if max_b != "D":
                                    consensus += max_b
                                    con_b = max_b
                                    no_cons += 1
                                else:
                                    print("Deletion ", end="")
                                    no_dels += 1
                                    con_b = "N"
                                    no_nons += 1
                                no_total += 1
                            else:
                                consensus += "N"
                                con_b = "N"
                                no_nons += 1
                                no_total += 1
                            print(f" => {p:0.2f} {con_b}")
                        else:
                            print("   N = 0!!!")
            gap.set_seq(consensus)
            if len(consensus) != (gap.send - gap.sstart + 1):
                print(
                    f"consensus sequence length not expected: {len(consensus)} vs {gap.send-gap.sstart+1}"
                )
            else:
                print(
                    f"consensus sequence length expected: {len(consensus)} vs {gap.send-gap.sstart+1}"
                )
            print(f"consensus# : {no_cons}")
            print(f"n#         : {no_nons}")
            print(f"total#     : {no_total}")
            print(f"depth      : {depth}")
            print(f"mean depth : {depth/no_total:.2f}")

            t_no_cons += no_cons
            t_no_nons += no_nons
            t_no_total += no_total
            t_depth += depth

    print(f"consensus# : {t_no_cons}")
    print(f"n#         : {t_no_nons}")
    print(f"total#     : {t_no_total}")
    print(f"depth      : {t_depth}")
    print(f"mean depth : {t_depth/t_no_total:.2f}")

    print("======================")
    print("print gaps and contigs")
    print("======================")
    with open(gaps_contigs_fasta_file, "w") as f:
        for contig in SEGMENTS:
            if type(contig).__name__ == "Gap":
                f.write(f">{contig.sacc}\n")
                f.write(f"{contig.seq}\n")
            else:
                f.write(f">{contig.sacc}\n")
                f.write(f"{contig.seq}\n")

    print("==============")
    print("print stitched")
    print("==============")
    stitched_seq = ""
    for contig in SEGMENTS:
        stitched_seq += contig.seq
    with open(stitched_fasta_file, "w") as f:
        f.write(f">{lion} stitched sequence len={len(stitched_seq)}\n")
        f.write(f"{stitched_seq}\n")


class Chromosome(object):
    def __init__(self):
        self.name = ""
        self.len = 0

    def __init__(self, name, len):
        self.name = name
        self.len = len

    def print(self):
        print(f"+++++++++++++++++++")
        print(f"chrom  : {self.name}")
        print(f"len    : {self.len}")


class Contig(object):
    def __init__(self):
        self.qacc = ""
        self.sacc = ""
        self.qlen = 0
        self.slen = 0
        self.sstart = 0
        self.send = 0
        self.score = 0
        self.length = 0
        self.pident = ""
        self.stitle = ""
        self.seq = ""

    def __init__(
        self, qacc, sacc, qlen, slen, sstart, send, score, length, pident, stitle
    ):
        self.qacc = qacc
        self.sacc = sacc
        self.qlen = int(qlen)
        self.slen = int(slen)
        self.sstart = int(sstart)
        self.send = int(send)
        self.score = float(score)
        self.length = int(length)
        self.pident = float(pident)
        self.stitle = stitle

    def set_seq(self, seq):
        if len(seq) == self.qlen:
            self.seq = seq
            print(f"{self.qacc} successfully populated")
        else:
            self.print()
            print(f"Sequence not consistent with reported length : {len(seq)}")

    def have_indels(self):
        return self.send - self.sstart + 1 != self.qlen

    # positive - insertions wrt ref
    # negative - deletions wrt ref
    def get_delta(self):
        return self.qlen - (self.send - self.sstart + 1)

    def print(self):
        print(f"+++++++++++++++++++")
        print(f"qacc    : {self.qacc}")
        print(f"slen    : {self.slen}")
        print(f"qlen    : {self.qlen}")
        print(f"sstart  : {self.sstart}")
        print(f"send    : {self.send}")
        print(f"indels  : {self.have_indels()}")
        print(f"delta   : {self.get_delta()}")


class Gap(object):
    def __init__(self):
        self.sacc = ""
        self.sstart = 0
        self.send = 0
        self.seq = ""

    def __init__(self, sacc, sstart, send):
        self.sacc = sacc
        self.sstart = int(sstart)
        self.send = int(send)

    def set_seq(self, seq):
        self.seq = seq

    def print(self):
        print(f"+++++++++++++++++++")
        print(f"sacc    : {self.sacc}")
        print(f"sstart  : {self.sstart}")
        print(f"send    : {self.send}")
        print(f"len     : {self.send-self.sstart+1}")


if __name__ == "__main__":
    main()
