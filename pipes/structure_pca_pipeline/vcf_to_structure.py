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

@click.command()
@click.argument("vcf_file")
@click.option("-o", "--output_dir", required=True, help="output directory")
@click.option("-d", "--dataset", required=True, help="dataset")
def main(vcf_file, output_dir, dataset):
    """
    Convert VCF file to structure format.  Generates the extraparams and mainparams files

    e.g. vcf_to_structure.py
    """
    print("\t{0:<20} :   {1:<10}".format("vcf file", vcf_file))

    # read VCF file, obtain master matrix of data
    data = []
    samples = []
    vcf_hdr = ""
    no_variants = 0
    no_samples = 0
    with open(vcf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.rstrip().split("\t")[9:]
                    no_samples = len(samples)
                else:
                    vcf_hdr += line
            else:
                chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = line.rstrip().split("\t")
                data.append(Variant(id, chrom, pos, ref, alt, genotypes))
                no_variants +=1

    out_structure_file = os.path.join(output_dir, f"{dataset}.structure")
    with open(out_structure_file, "w") as file:
        #write locus header
        for i in range(no_variants):
            file.write(f"\tsnp_{i}")
        file.write(f"\n")
        #write sample data
        for j in range(no_samples):
            sample_line1 = samples[j]
            sample_line2 = samples[j]
            for i in range(no_variants):
                gt = data[i].genotypes[j].gt
                if gt == -1:
                    sample_line1 += "\t-9"
                    sample_line2 += "\t-9"
                elif gt == 0:
                    sample_line1 += "\t0"
                    sample_line2 += "\t0"
                elif gt == 1:
                    sample_line1 += "\t0"
                    sample_line2 += "\t1"
                elif gt == 2:
                    sample_line1 += "\t1"
                    sample_line2 += "\t1"
            file.write(f"{sample_line1}\n")
            file.write(f"{sample_line2}\n")

    out_mainparams_file = os.path.join(output_dir, "mainparams")
    with open(out_mainparams_file, "w") as file:
        mainparams = f"""
        KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
        IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
        FILE extraparams.


        "(int)" means that this takes an integer value.
        "(B)"   means that this variable is Boolean
                (ie insert 1 for True, and 0 for False)
        "(str)" means that this is a string (but not enclosed in quotes!)


        Basic Program Parameters

        #define MAXPOPS    2      // (int) number of populations assumed
        #define BURNIN    10000   // (int) length of burnin period
        #define NUMREPS   20000   // (int) number of MCMC reps after burnin

        Input/Output files

        #define INFILE   infile   // (str) name of input data file
        #define OUTFILE  outfile  //(str) name of output data file

        Data file format

        #define NUMINDS   {no_samples}  // (int) number of diploid individuals in data file
        #define NUMLOCI   {no_variants} // (int) number of loci in data file
        #define PLOIDY       2    // (int) ploidy of data
        #define MISSING     -9    // (int) value given to missing genotype data
        #define ONEROWPERIND 0    // (B) store data for individuals in a single line


        #define LABEL     1     // (B) Input file contains individual labels
        #define POPDATA   0     // (B) Input file contains a population identifier
        #define POPFLAG   0     // (B) Input file contains a flag which says
                                    whether to use popinfo when USEPOPINFO==1
        #define LOCDATA   0     // (B) Input file contains a location identifier

        #define PHENOTYPE 0     // (B) Input file contains phenotype information
        #define EXTRACOLS 0     // (int) Number of additional columns of data
                                    before the genotype data start.

        #define MARKERNAMES      1  // (B) data file contains row of marker names
        #define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                                    // and a row to indicate which alleles are recessive
        #define MAPDISTANCES     0  // (B) data file contains row of map distances
                                    // between loci


        Advanced data file options

        #define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
        #define PHASEINFO        0 // (B) the data for each individual contains a line
                                        indicating phase (linkage model)
        #define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
        #define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data

        Command line options:

        -m mainparams
        -e extraparams
        -s stratparams
        -K MAXPOPS
        -L NUMLOCI
        -N NUMINDS
        -i input file
        -o output file
        -D SEED
        """
        file.write(mainparams)

    out_extraparams_file = os.path.join(output_dir, "extraparams")
    with open(out_extraparams_file, "w") as file:
        extraparams = f"""
            EXTRA PARAMS FOR THE PROGRAM structure.  THESE PARAMETERS CONTROL HOW THE
            PROGRAM RUNS.  ATTRIBUTES OF THE DATAFILE AS WELL AS K AND RUNLENGTH ARE
            SPECIFIED IN mainparams.

            "(int)" means that this takes an integer value.
            "(d)"   means that this is a double (ie, a Real number such as 3.14).
            "(B)"   means that this variable is Boolean
                    (ie insert 1 for True, and 0 for False).

            PROGRAM OPTIONS

            #define NOADMIX     0 // (B) Use no admixture model (0=admixture model, 1=no-admix)
            #define LINKAGE     0 // (B) Use the linkage model model
            #define USEPOPINFO  0 // (B) Use prior population information to pre-assign individuals
                                        to clusters
            #define LOCPRIOR    0 //(B)  Use location information to improve weak data

            #define FREQSCORR   1 // (B) allele frequencies are correlated among pops
            #define ONEFST      0 // (B) assume same value of Fst for all subpopulations.

            #define INFERALPHA  1 // (B) Infer ALPHA (the admixture parameter)
            #define POPALPHAS   0 // (B) Individual alpha for each population
            #define ALPHA     1.0 // (d) Dirichlet parameter for degree of admixture
                                        (this is the initial value if INFERALPHA==1).

            #define INFERLAMBDA 0 // (B) Infer LAMBDA (the allele frequencies parameter)
            #define POPSPECIFICLAMBDA 0 //(B) infer a separate lambda for each pop
                                (only if INFERLAMBDA=1).
            #define LAMBDA    1.0 // (d) Dirichlet parameter for allele frequencies

            PRIORS

            #define FPRIORMEAN 0.01 // (d) Prior mean and SD of Fst for pops.
            #define FPRIORSD   0.05  // (d) The prior is a Gamma distribution with these parameters

            #define UNIFPRIORALPHA 1 // (B) use a uniform prior for alpha;
                                            otherwise gamma prior
            #define ALPHAMAX     10.0 // (d) max value of alpha if uniform prior
            #define ALPHAPRIORA   1.0 // (only if UNIFPRIORALPHA==0): alpha has a gamma
                                        prior with mean A*B, and
            #define ALPHAPRIORB   2.0 // variance A*B^2.

            #define LOG10RMIN     -4.0   //(d) Log10 of minimum allowed value of r under linkage model
            #define LOG10RMAX      1.0   //(d) Log10 of maximum allowed value of r
            #define LOG10RPROPSD   0.1   //(d) standard deviation of log r in update
            #define LOG10RSTART   -2.0   //(d) initial value of log10 r

            USING PRIOR POPULATION INFO (USEPOPINFO)

            #define GENSBACK    2  //(int) For use when inferring whether an indiv-
                                    idual is an immigrant, or has an immigrant an-
                                    cestor in the past GENSBACK generations.  eg, if
                                    GENSBACK==2, it tests for immigrant ancestry
                                    back to grandparents.
            #define MIGRPRIOR 0.01 //(d) prior prob that an individual is a migrant
                                        (used only when USEPOPINFO==1).  This should
                                        be small, eg 0.01 or 0.1.
            #define PFROMPOPFLAGONLY 0 // (B) only use individuals with POPFLAG=1 to update	P.
                                            This is to enable use of a reference set of
                                            individuals for clustering additional "test"
                                            individuals.

            LOCPRIOR MODEL FOR USING LOCATION INFORMATION

            #define LOCISPOP      0    //(B) use POPDATA for location information
            #define LOCPRIORINIT  1.0  //(d) initial value for r, the location prior
            #define MAXLOCPRIOR  20.0  //(d) max allowed value for r

            OUTPUT OPTIONS

            #define PRINTNET     1 // (B) Print the "net nucleotide distance" to screen during the run
            #define PRINTLAMBDA  1 // (B) Print current value(s) of lambda to screen
            #define PRINTQSUM    1 // (B) Print summary of current population membership to screen

            #define SITEBYSITE   0  // (B) whether or not to print site by site results.
                                (Linkage model only) This is a large file!
            #define PRINTQHAT    0  // (B) Q-hat printed to a separate file.  Turn this
                                    on before using STRAT.
            #define UPDATEFREQ   100  // (int) frequency of printing update on the screen.
                                            Set automatically if this is 0.
            #define PRINTLIKES   0  // (B) print current likelihood to screen every rep
            #define INTERMEDSAVE 0  // (int) number of saves to file during run

            #define ECHODATA     1  // (B) Print some of data file to screen to check
                                        that the data entry is correct.
            (NEXT 3 ARE FOR COLLECTING DISTRIBUTION OF Q:)
            #define ANCESTDIST   0  // (B) collect data about the distribution of an-
                                        cestry coefficients (Q) for each individual
            #define NUMBOXES   1000 // (int) the distribution of Q values is stored as
                                        a histogram with this number of boxes.
            #define ANCESTPINT 0.90 // (d) the size of the displayed probability
                                        interval on Q (values between 0.0--1.0)

            MISCELLANEOUS

            #define COMPUTEPROB 1     // (B) Estimate the probability of the Data under
                                        the model.  This is used when choosing the
                                        best number of subpopulations.
            #define ADMBURNIN  500    // (int) [only relevant for linkage model]:
                                        Initial period of burnin with admixture model (see Readme)
            #define ALPHAPROPSD 0.025 // (d) SD of proposal for updating alpha
            #define STARTATPOPINFO 0  // Use given populations as the initial condition
                                        for population origins.  (Need POPDATA==1).  It
                                        is assumed that the PopData in the input file
                                        are between 1 and k where k<=MAXPOPS.
            #define RANDOMIZE      1  // (B) use new random seed for each run
            #define SEED        2245  // (int) seed value for random number generator
                                    (must set RANDOMIZE=0)
            #define METROFREQ    10   // (int) Frequency of using Metropolis step to update
                                        Q under admixture model (ie use the metr. move every
                                        i steps).  If this is set to 0, it is never used.
                                        (Proposal for each q^(i) sampled from prior.  The
                                        goal is to improve mixing for small alpha.)
            #define REPORTHITRATE 0 //   (B) report hit rate if using METROFREQ
        """
        file.write(extraparams)

class Variant(object):
    def __init__(self, id, chrom, pos, ref, alt, genotypes):
        self.id = id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.genotypes = []
        for g in genotypes:
            self.add_genotype(g)
        self.ts = self.is_ts(ref,alt)
        #print(f"tv: {self.tv} | {ref} | {alt}")

    def is_ts(self, ref, alt):
        return (ref == "A" and alt == "G") or \
        (ref == "G" and alt == "A") or \
        (ref == "C" and alt == "T") or \
        (ref == "T" and alt == "C")

    def add_genotype(self, genotype):
        if genotype == "./.":
            self.genotypes.append(Genotype(-1, -1, -1, -1, "-1,-1,-1"))
        else:
            gt, dp, ad, gq, gl = genotype.split(":")
            if gt == "./.":
                #print(genotype)
                self.genotypes.append(Genotype(-1, -1, -1, -1, "-1,-1,-1"))
            else:
                if gt == "0/0":
                    gt = 0
                elif gt == "0/1":
                    gt = 1
                elif gt == "1/1":
                    gt = 2
                self.genotypes.append(Genotype(gt, int(dp), ad, int(gq), gl))

class Genotype(object):
    def __init__(self, gt, dp, ad, gq, gl):
        self.gt = gt
        self.dp = dp
        if ad == -1:
            self.ad = -1
        else:
            self.ad = int(ad.split(",")[1])
        self.gq = gq
        self.gl = gl.split(",")

if __name__ == "__main__":
    main() # type: ignore
