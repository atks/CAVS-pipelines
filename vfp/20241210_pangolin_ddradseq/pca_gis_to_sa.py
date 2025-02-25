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
import re

@click.command()
@click.option("-p", "--pca_sa_file", required=True, help="pca sample file")
@click.option("-g", "--gis_sa_file", required=True, help="gis sample file")
@click.option("-o", "--output_sa_file", required=True, help="output sample file")
def main(pca_sa_file, gis_sa_file, output_sa_file):
    """
    Combined pca PCs with gis coordinates.

    e.g. pca_gis_to_sa.py
    """
    print("\t{0:<20} :   {1:<10}".format("pca sample file", pca_sa_file))
    print("\t{0:<20} :   {1:<10}".format("gis_sample_file", gis_sa_file))
    print("\t{0:<20} :   {1:<10}".format("output_sa_file", output_sa_file))

    print(f"processing {pca_sa_file}")

    #sample-id	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15	PC16	PC17	PC18	PC19	PC20
    print(f"processing {gis_sa_file}")

    samples = {}

    with open(pca_sa_file, "r") as f:
        for line in f.readlines():
            if line.startswith("sample-id"):
                continue
            sample_id, *pc = line.split()
            samples[sample_id] = SampleRecord(sample_id, *pc)


    # sample-id	latitude	longitude	sex
    # BIOS0005	Missing	Missing	Unknown
    # BIOS0006	1.33555	103.7986	F
    # BIOS0007	1.376321	103.712874	M
    # BIOS0008	1.382438	103.775458	F
    # BIOS0009	1.410534	103.827636	M
    # BIOS0010	1.368223	103.778963	F
    # BIOS0011	1.35404	103.682833	M
    # BIOS0012	1.371556	103.778181	Unknown
    # BIOS0013	1.36449	103.77235	M

    with open(output_sa_file, "w") as o:
        o.write("sample_id\tlatitude\tlongitude\tsex")
        for i in range(int(20)):
            o.write(f"\tPC{i+1}")
        o.write("\n")
        with open(gis_sa_file, "r") as f:
            for line in f.readlines():
                sample_id, latitude, longitude, sex = line.split()
                if latitude!="Missing" and sample_id in samples:
                    l = f"{sample_id}\t{latitude}\t{longitude}\t{sex}"
                    l += f"\t{samples[sample_id].print_to_str()}"
                    o.write(f"{l}\n")

class SampleRecord(object):

    def __init__(self, sample_id, *pcs):
        self.sample_id = sample_id
        self.pcs = pcs

    def print_to_str(self):
        line = f"{self.sample_id}"
        for pc in self.pcs:
            line += f"\t{float(pc):.3f}"
        return line

if __name__ == "__main__":
    main() # type: ignore
