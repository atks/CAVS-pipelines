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
        file = os.path.abspath(file)
        print(f"processing {file}")

        #read file
        with open(file, "r") as f:
            no_individuals = -1
            in_individuals_inferred_cluster = False
            no_individuals_read = 0
            k = 0
            reps = []
            lines = f.readlines()
            for line in lines:
                if no_individuals == -1:
                    m = re.match(r"\s+(\d+) individuals", line)
                    if m is not None:
                        no_individuals = int(m.group(1))
                        print(f"\tno_individuals: {no_individuals}")
                        continue
                elif k == 0:
                    m = re.match(r"\s+(\d+) populations assumed", line)
                    if m is not None:
                        k = int(m.group(1))
                        reps.append(StructureRep(file, k))
                        print(f"\tk: {k}")
                        continue
                else:
                    if no_individuals_read < no_individuals:
                        if in_individuals_inferred_cluster:
                            no_individuals_read += 1
                            line = line.rstrip()
                            reps[-1].add_individual_record_line(line)
                            #print(line)
                        else:
                            if line.endswith("Inferred clusters\n"):
                                in_individuals_inferred_cluster = True
                                #print("\tin cluster")
                                continue
                    else:
                        #print("\tout of cluster")
                        break

        print("writing out distruct files")
        #write distruct files
        for rep in reps:
            rep.sort_individual_records()
            indivq_file = rep.file.replace("_f", ".indivq")
            with open(indivq_file, "w") as f:
                f.write(rep.print_indiv_str())
            popq_file = rep.file.replace("_f", ".popq")
            with open(popq_file, "w") as f:
                f.write(rep.print_popq_str())

            ps_file = rep.file.replace("_f", ".ps")
            drawparams = f"""
PARAMETERS FOR THE PROGRAM distruct.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.

"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean
        (1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes)
"(d)"   means that this is a double (a real number).

Data settings

#define INFILE_POPQ        {os.path.basename(popq_file)}      // (str) input file of population q's
#define INFILE_INDIVQ      {os.path.basename(indivq_file)}    // (str) input file of individual q's
#define INFILE_LABEL_BELOW input.names     // (str) input file of labels for below figure
#define INFILE_LABEL_ATOP  input.languages // (str) input file of labels for atop figure
#define INFILE_CLUST_PERM  input.perm     // (str) input file of permutation of clusters to print
#define OUTFILE            {os.path.basename(ps_file)}       //(str) name of output file

#define K	{rep.k}    // (int) number of clusters
#define NUMPOPS 1    // (int) number of pre-defined populations
#define NUMINDS {rep.n}  // (int) number of individuals

Main usage options

#define PRINT_INDIVS      1  // (B) 1 if indiv q's are to be printed, 0 if only population q's
#define PRINT_LABEL_ATOP  0  // (B) print labels above figure
#define PRINT_LABEL_BELOW 0  // (B) print labels below figure
#define PRINT_SEP         0  // (B) print lines to separate populations

Figure appearance

#define FONTHEIGHT 6	// (d) size of font
#define DIST_ABOVE 5	// (d) distance above plot to place text
#define DIST_BELOW -7	// (d) distance below plot to place text
#define BOXHEIGHT  36	// (d) height of the figure
#define INDIVWIDTH 1.5	// (d) width of an individual


Extra options

#define ORIENTATION 0	     // (int) 0 for horizontal orientation (default)
			     //       1 for vertical orientation
			     //	      2 for reverse horizontal orientation
                             //       3 for reverse vertical orientation
#define XORIGIN 72		// (d) lower-left x-coordinate of figure
#define YORIGIN 288		// (d) lower-left y-coordinate of figure
#define XSCALE 1		// (d) scale for x direction
#define YSCALE 1		// (d) scale for y direction
#define ANGLE_LABEL_ATOP 60	// (d) angle for labels atop figure (in [0,180])
#define ANGLE_LABEL_BELOW 60    // (d) angle for labels below figure (in [0,180])
#define LINEWIDTH_RIM  3	// (d) width of "pen" for rim of box
#define LINEWIDTH_SEP 0.3	// (d) width of "pen" for separators between pops and for tics
#define LINEWIDTH_IND 0.3	// (d) width of "pen" used for individuals
#define GRAYSCALE 0	        // (B) use grayscale instead of colors
#define ECHO_DATA 1             // (B) print some of the data to the screen
#define REPRINT_DATA 1          // (B) print the data as a comment in the ps file
#define PRINT_INFILE_NAME 0     // (B) print the name of INFILE_POPQ above the figure
                                //     this option is meant for use only with ORIENTATION=0
#define PRINT_COLOR_BREWER 1    // (B) print ColorBrewer settings in the output file
                                //     this option adds 1689 lines and 104656 bytes to the output
                                //     and is required if using ColorBrewer colors


Command line options:

-d drawparams
-K K
-M NUMPOPS
-N NUMINDS
-p input file (population q's)
-i input file (individual q's)
-a input file (labels atop figure)
-b input file (labels below figure)
-c input file (cluster permutation)
-o output file
"""
            drawparam_file = rep.file.replace("_f", ".drawparams")
            with open(drawparam_file, "w") as f:
                f.write(drawparams)






            #rep.print()

    #output clumpp files

class StructureRep(object):

    def __init__(self, file, k):
        self.file = file
        self.k = k
        self.n = 0
        self.indiv_records = []

    def add_individual_record_line(self, line):
        self.indiv_records.append(IndividualInferredClusterRecord(line))

    def sort_individual_records(self):
        #detect highest cluster representation
        membership_sums = [0]*self.k
        for rec in self.indiv_records:
            for i in range(self.k):
                membership_sums[i] += rec.cluster_memberships[i]

        self.n = len(self.indiv_records)

        #check highest cluster representation
        highest_ksum = 0
        highest_k = 0
        for i in range(self.k):
            if membership_sums[i] > highest_ksum:
                highest_ksum = membership_sums[i]
                highest_k = i

        self.indiv_records.sort(key=lambda x: x.cluster_memberships[highest_k], reverse=True)

    def print_popq_str(self):
        total = 0.0
        membership_proportion = [0.0]*self.k
        for rec in self.indiv_records:
            for i in range(self.k):
                total += rec.cluster_memberships[i]
                membership_proportion[i] += rec.cluster_memberships[i]

        #normalise
        for i in range(self.k):
             membership_proportion[i] = membership_proportion[i]/total

        line = f"1:"
        for i in range(len(membership_proportion)):
            line += f" {membership_proportion[i]:.3f}"
        line += f" {len(self.indiv_records)}\n"

        return line


    def print_indiv_str(self):
        str = ""
        for rec in self.indiv_records:
            str += rec.print_to_str() + "\n"
        return str

    def print(self):
        print(f"file: {self.file}")
        print(f"k: {self.k}")
        for rec in self.indiv_records:
            print(rec.print_to_str())

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
        # 47 BIOS1470    (0)   :  0.493 0.507
        line = f"{self.sample_idx:6d}  {self.sample_id}    ({self.missing})   1 : "
        for i in range(len(self.cluster_memberships)):
            line += f" {self.cluster_memberships[i]:.3f}"
        return line

    def print(self):
        print(f"idx: {self.sample_idx}")
        print(f"id: {self.sample_id}")
        print(f"missing: {self.missing}")
        print(f"cluster memberships: {self.cluster_memberships}")
        print(f"line: {self.line}")

if __name__ == "__main__":
    main() # type: ignore
