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
@click.option("-s", "--structure_file", required=True, help="structure file")
@click.option("-p", "--pca_sa_file", required=True, help="pca sample file")
@click.option("-o", "--output_sa_file", required=True, help="output sample file")
def main(structure_file, pca_sa_file, output_sa_file):
    """
    Combined structure coefficients with pca coordinates.

    e.g. structure_pca_to_sa.py
    """
    print("\t{0:<20} :   {1:<10}".format("structure files", structure_file))
    print("\t{0:<20} :   {1:<10}".format("pca_sa_file", pca_sa_file))
    print("\t{0:<20} :   {1:<10}".format("output_sa_file", output_sa_file))

    structure_file = os.path.abspath(structure_file)
    print(f"processing {structure_file}")
    no_clusters = 0

    with open(structure_file, "r") as f:
        no_individuals = -1
        in_individuals_inferred_cluster = False
        no_individuals_read = 0
        k = 0
        rep = {}
        lines = f.readlines()
        for line in lines:
            if no_individuals == -1:
                m = re.match(r"\s+(\d+) individuals", line)
                if m is not None:
                    no_individuals = int(m.group(1))
                    print(f"\tno_individuals: {no_individuals}")
                    continue
            else:
                if no_individuals_read < no_individuals:
                    if in_individuals_inferred_cluster:
                        no_individuals_read += 1
                        line = line.rstrip()
                        rec = IndividualInferredClusterRecord(line)
                        no_clusters = len(rec.cluster_memberships)
                        rep[rec.sample_id] = rec
                    else:
                        if line.endswith("Inferred clusters\n"):
                            in_individuals_inferred_cluster = True
                            continue
                else:
                    break

    print(f"processing {pca_sa_file}")

    #sample-id	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15	PC16	PC17	PC18	PC19	PC20

    with open(output_sa_file, "w") as o:
        o.write("sample_id")
        for i in range(20):
            o.write(f"\tPC{i+1}")
        for i in range(int(no_clusters)):
            o.write(f"\tC{i+1}")
        o.write("\n")
        with open(pca_sa_file, "r") as f:
            for line in f.readlines():
                sample_id, *pcs = line.split()
                if sample_id in rep:
                    l = f"{sample_id}"
                    for pc in pcs:
                        l += f"\t{float(pc):.3f}"
                    for m in rep[sample_id].cluster_memberships:
                        l += f"\t{m:.3f}"
                    o.write(f"{l}\n")
class IndividualInferredClusterRecord(object):

    def __init__(self, line):
        tokens = line.split()
        self.sample_idx = int(tokens[0])
        self.sample_id = tokens[1]
        self.missing = int(tokens[2].strip("()"))
        self.cluster_memberships = []
        for i in range(4, len(tokens)):
            self.cluster_memberships.append(float(tokens[i]))
            #print(f"adding cluster membership {i-3}")
        self.line = line

    def print_to_str(self):
        line = f"{self.sample_id}"
        for i in range(len(self.cluster_memberships)):
            line += f"\t{self.cluster_memberships[i]:.3f}"
        return line

    def print(self):
        print(f"idx: {self.sample_idx}")
        print(f"id: {self.sample_id}")
        print(f"missing: {self.missing}")
        print(f"cluster memberships: {self.cluster_memberships}")
        print(f"line: {self.line}")

if __name__ == "__main__":
    main() # type: ignore
