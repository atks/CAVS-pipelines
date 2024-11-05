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
import re

@click.command()
@click.argument("general_stats_txt", required=True)
@click.option("-o", "--output_extracted_general_stats_file", required=True, help="output global coverage statistics file")
@click.option("-s", "--sample_id", required=True, help="sample ID")
def main(general_stats_txt, output_extracted_general_stats_file, sample_id):
    """
    Extracts general stats

    e.g. generate_vfp_novogene_ddradseq_mina_deploy_and_qc_pipeline -r mina1 -i raw -s mina1.sa -g /usr/local/ref/vfp/ManJav1.0_HiC.fasta
    """
    total_reads = 0
    mapped_reads = 0
    mapping_rate = 0
    insert_mean = 0
    insert_sd = 0
    readlen_mean = 0

    with open(general_stats_txt, "r") as file:
        for line in file:
            if line.startswith("SN"):
                m = re.search(r"SN\s+raw total sequences:\s+(\d+)", line)
                if m is not None:
                    total_reads = int(m.group(1))
                    continue
                m = re.search(r"SN\s+reads mapped:\s+(\d+)", line)
                if m is not None:
                    mapped_reads = int(m.group(1))
                    continue
                m = re.search(r"SN\s+insert size average:\s+(\d+)", line)
                if m is not None:
                    insert_mean = int(m.group(1))
                    continue
                m = re.search(r"SN\s+insert size standard deviation:\s+(\d+)", line)
                if m is not None:
                    insert_sd = float(m.group(1))
                    continue
                m = re.search(r"SN\s+average length:\s+(\d+)", line)
                if m is not None:
                    readlen_mean = m.group(1)
                    continue

    mapping_rate = mapped_reads / total_reads

    with open(output_extracted_general_stats_file, "w") as file:
       file.write(f"#id\ttotrawreads\tmapped_reads\tmapping_rate\tinsertsize_mean\tinsert_sd\treadlen_mean\n")

       file.write(f"{sample_id}\t{total_reads}\t{mapped_reads}\t{mapping_rate:.2f}\t{insert_mean}\t{insert_sd:.2f}\t{readlen_mean}\n")

if __name__ == "__main__":
    main() # type: ignore
