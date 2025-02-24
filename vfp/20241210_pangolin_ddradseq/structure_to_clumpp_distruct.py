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
@click.argument("structure_files", nargs=-1)
def main(structure_files):
    """
    Convert structure output files to clump and distruct input files.

    e.g. structure_to_clump_distruct.py
    """
    print("\t{0:<20} :   {1:<10}".format("structure files", ",".join(structure_files)))


    K = 0
    reps = []

    for file in structure_files:
        with open(file, "r") as f:
            no_line = 1
            no_individuals = -1
            k = 0
            no_individuals_read = 0
            in_individuals_inferred_cluster = False
            lines = f.readlines()
            for line in lines:
                if no_line == 15:
                    m = re.search(r"(\d+) individuals\n", line)
                    if m is not None:
                        no_individuals = int(m.group(1))
                    else:
                        raise ValueError("Could not determine number of individuals from line: {0}".format(line))
                        exit(1)

                if no_line == 17:
                    m = re.search(r"(\d+) populations assumed\n", line)
                    if m is not None:
                        k = int(m.group(1))
                        if K == 0:
                            K = k
                        else:
                            assert K==k, "Number of populations is not consistent across files"
                    else:
                        raise ValueError("Could not determine assumed populations from line: {0}".format(line))
                        exit(1)

                if in_individuals_inferred_cluster:
                    if no_individuals_read < no_individuals:
                        no_individuals_read += 1
                        continue
                    else:
                        break





                if line == "Inferred ancestry of individuals:":
                    in_individuals_inferred = True
                    continue

                no_line += 1


class StructureRep(object):

    def __init__(self):
        self.indiv_records = []

    def add_individual_record(self, record):
        self.indiv_records.append(record)

    def sort_individual_records(self):
        self.indiv_records.sort(key=lambda x: x.individual_id)
class IndividualInferredClusterLine(object):

    def __init__(self, cluster1, record):
        self.cluster1 = cluster1
        self.record = record


if __name__ == "__main__":
    main() # type: ignore
