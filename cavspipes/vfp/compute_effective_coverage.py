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

import click

@click.command()
@click.argument("coverage_stats_txt", required=True)
@click.option("-o", "--output_global_coverage_stats_file", required=True, help="output global coverage statistics file")
@click.option("-s", "--sample_id", required=True, help="sample ID")
def main(coverage_stats_txt, output_global_coverage_stats_file, sample_id):
    """
    Computes global coverage

    e.g. generate_vfp_novogene_ddradseq_mina_deploy_and_qc_pipeline -r mina1 -i raw -s mina1.sa -g /usr/local/ref/vfp/ManJav1.0_HiC.fasta
    """
    # #rname              startpos  endpos     numreads  covbases  coverage   meandepth   meanbaseq  meanmapq
    # HiC_scaffold_1      1         62011232   83151     1988296   3.20635    0.175215    35.9       53.3
    # HiC_scaffold_2      1         24701697   39132     949688    3.84463    0.208846    36         56
    # HiC_scaffold_3      1         111765584  152409    3490039   3.12264    0.178298    35.9       52.4
    # HiC_scaffold_4      1         91963190   134393    2751959   2.99246    0.17853     35.1       47.4
    # HiC_scaffold_5      1         88098640   136733    2917699   3.31185    0.200028    35.7       48.8

    with open(coverage_stats_txt, "r") as file:
        total_covbases = 0
        total_bases = 0
        total_reads = 0
        for line in file:
            if not line.startswith("#"):
                rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq = line.rstrip().split("\t")
                total_bases += float(meandepth)*int(endpos)
                total_covbases += int(covbases)
                total_reads += int(numreads)

    effective_coverage = total_bases / total_covbases

    with open(output_global_coverage_stats_file, "w") as file:
       file.write(f"#id\ttotbases\ttotreads\tcovbases\teffcov\n")
       file.write(f"{sample_id}\t{int(total_bases)}\t{int(total_reads)}\t{total_covbases}\t{effective_coverage:.2f}\n")

if __name__ == "__main__":
    main() # type: ignore
