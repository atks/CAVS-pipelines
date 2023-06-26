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
import re
import fnmatch


@click.command()
@click.argument("ilm_dir", required=True)
def main(ilm_dir):
    """
    Takes in a ilm directory and generates a fastq file list

    e.g. make_ilm_fastq_list.py ilm38
    """
    fastq_path = ""

    for dirpath, dirnames, filenames in os.walk(ilm_dir):
        for dirname in dirnames:
            if dirname == "Fastq":
                fastq_path = os.path.join(dirpath, dirname)

    ilm_dir = os.path.abspath(fastq_path)

    no_fastq_files = int(len(fnmatch.filter(os.listdir(ilm_dir), "*.fastq.gz")) / 2)
    # print(f'#fastq : {no_fastq_files}')

    samples = [""] * no_fastq_files
    run_id = ""
    for file_name in os.listdir(ilm_dir):
        if file_name.endswith("fastq.gz"):
            # print(file_name)
            m = re.match("^(.+_S(\d+)_L\d+)_R[12].+.fastq.gz", file_name)
            sample_id = int(m.group(2))
            name = m.group(1)
            # print(f'{sample_id} {name}')
            samples[sample_id] = name

    # shift undeterminate sample to last
    temp = samples.pop(0)
    samples.append(temp)

    for id, name in enumerate(samples):
        print(f"{id+1}\t{name}_R1_001.fastq.gz\t{name}_R2_001.fastq.gz")


if __name__ == "__main__":
    main()  # type: ignore
