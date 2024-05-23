#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2022 Adrian Tan <adrian_tan@nparks.gov.sg>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permsit persons to whom the Software is
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
@click.option('-w', '--results_dir', default=os.getcwd(), show_default=True, help='results directory')
@click.option('-s', '--sample_file', required=True, help='sample file')
@click.option('-a', '--amr_ast_results_file', required=True, help='AMR AST results file')
def main(results_dir, sample_file, amr_ast_results_file):
    """
    Collate AMR results

    e.g. collate_amr_results.py ilm38
   """
    print(f"results directory: {results_dir}")
    print(f"sample  file: {sample_file}")
    print(f"AMR AST results file: {amr_ast_results_file}")

    #read AST results file
    samples = []
    headers = []
    with open(amr_ast_results_file, 'r') as file:
        index = 0
        for line in file:
            if line.startswith('#'):
                headers = line.rstrip().split('\t')
            else:
                fields = line.rstrip().split('\t')
                sample_id = fields[0]
                species = fields[1]
                sample = Sample(sample_id, species)
                for i in range(2,len(fields)):
                    antimicrobial = Antimicrobial(headers[i])
                    antimicrobial.ast_pheno = fields[i]
                    sample.add_am(antimicrobial)
                samples.append(sample)

    #go through results directories
    for sample in samples:
        resfinder_results_file = f'{results_dir}/{sample.id}/pheno_table.txt'
        to_process = False
        with open(resfinder_results_file, 'r') as file:
            for line in file:
                if not to_process:
                    if line.startswith('# Antimicrobial'):
                        to_process = True
                    else:
                        pass
                else:
                    if line != '\n':
                        fields = line.rstrip().split('\t')
                        antimicrobial = fields[0]
                        am_class = fields[1]
                        pred_pheno = fields[2]
                        if pred_pheno == 'Resistant':
                            pred_pheno = 'R'
                        elif pred_pheno == 'No resistance':
                            pred_pheno = 'S'
                        else:
                            print(f"Unexpected prediction {pred_pheno}")
                        match = fields[3]
                        genetic_background = fields[4] if match == 0 else ''
                        if antimicrobial in ('trimethoprim', 'sulfamethoxazole'):
                            combined_antimicrobial = 'trimethoprim/sulfamethoxazole'
                            if antimicrobial == 'trimethoprim':
                                sample.antimicrobials[combined_antimicrobial].ast_pred = pred_pheno + "/" + sample.antimicrobials[combined_antimicrobial ].ast_pred
                            else:
                                sample.antimicrobials[combined_antimicrobial].ast_pred = pred_pheno
                        else:
                            if antimicrobial in sample.antimicrobials:
                                sample.antimicrobials[antimicrobial].ast_pred = pred_pheno
                    else:
                        break

    sorted_antimicrobials = sorted(samples[0].antimicrobials)
    with open('summary.txt', 'w') as file:
        file.write('sample-id')
        for key in sorted_antimicrobials:
            file.write(f'\t{key}\t')
        file.write('\n')

        for key in sorted_antimicrobials:
            file.write(f'\tAST\tresfinder')
        file.write('\n')

        for sample in samples:
            file.write(sample.id)
            for i, antimicrobial in enumerate(sorted_antimicrobials):
                am = sample.antimicrobials[antimicrobial]
                file.write(f'\t{am.ast_pheno}\t{am.ast_pred}')
            file.write('\n')

            #sample.print()

    sample_stats = {}
    antimicrobial_stats = {}
    for sample in samples:
        sample_stats[sample.id] = Stats()
    for i, antimicrobial in enumerate(sorted_antimicrobials):
        antimicrobial_stats[antimicrobial] = Stats()

    #collect stats
    for sample in samples:
        for i, antimicrobial in enumerate(sorted_antimicrobials):
            am = sample.antimicrobials[antimicrobial]
            pheno = am.ast_pheno
            pred = am.ast_pred
            if '/' in pred:
                if pred == 'R/R':
                    pred = 'R'
                elif pred == 'R/S':
                    pred = 'S'
                elif pred == 'S/R':
                    pred = 'S'
                elif pred == 'S/S':
                    pred = 'S'
            if pheno != 'NI' and pred != '':
                sample_stats[sample.id].n += 1
                antimicrobial_stats[antimicrobial].n += 1
                if pheno == 'R' and pred == 'R':
                    sample_stats[sample.id].n_r += 1
                    sample_stats[sample.id].n_rr += 1
                    antimicrobial_stats[antimicrobial].n_rr += 1
                    antimicrobial_stats[antimicrobial].n_r += 1
                elif pheno == 'I' and pred == 'R':
                    sample_stats[sample.id].n_i += 1
                    sample_stats[sample.id].n_ir += 1
                    antimicrobial_stats[antimicrobial].n_ir += 1
                    antimicrobial_stats[antimicrobial].n_i += 1
                elif pheno == 'S' and pred == 'R':
                    sample_stats[sample.id].n_s += 1
                    sample_stats[sample.id].n_sr += 1
                    antimicrobial_stats[antimicrobial].n_sr += 1
                    antimicrobial_stats[antimicrobial].n_s += 1
                elif pheno == 'R' and pred == 'S':
                    sample_stats[sample.id].n_r += 1
                    sample_stats[sample.id].n_rs += 1
                    antimicrobial_stats[antimicrobial].n_rs += 1
                    antimicrobial_stats[antimicrobial].n_r += 1
                elif pheno == 'I' and pred == 'S':
                    sample_stats[sample.id].n_i += 1
                    sample_stats[sample.id].n_is += 1
                    antimicrobial_stats[antimicrobial].n_is += 1
                    antimicrobial_stats[antimicrobial].n_i += 1
                elif pheno == 'S' and pred == 'S':
                    sample_stats[sample.id].n_s += 1
                    sample_stats[sample.id].n_ss += 1
                    antimicrobial_stats[antimicrobial].n_ss += 1
                    antimicrobial_stats[antimicrobial].n_s += 1
                else:
                    print(f'not caught {antimicrobial}: "{pheno}" \t "{pred}"')
            else:
                #adhoc!!!!!
                if pheno == 'NI':
                    sample_stats[sample.id].n_m += 1
                    antimicrobial_stats[antimicrobial].n_m += 1
                elif pheno == 'S':
                    sample_stats[sample.id].n += 1
                    sample_stats[sample.id].n_s += 1
                    antimicrobial_stats[antimicrobial].n += 1
                    antimicrobial_stats[antimicrobial].n_s += 1


    print(f'sample-id\tn\tr\ti\ts\tni\trr\tir\tsr\trs\tis\tss')
    for sample in samples:
        stats = sample_stats[sample.id]
        print(f'{sample.id}\t{stats.n}\t{stats.n_r}\t{stats.n_i}\t{stats.n_s}\t{stats.n_m}\t{stats.n_rr}\t{stats.n_ir}\t{stats.n_sr}\t{stats.n_rs}\t{stats.n_is}\t{stats.n_ss}')

    print(f'antimicrobial\tn\tr\ti\ts\tni\trr\tir\tsr\trs\tis\tss')
    for i, antimicrobial in enumerate(sorted_antimicrobials):
        stats = antimicrobial_stats[antimicrobial]
        print(f'{antimicrobial}\t{stats.n}\t{stats.n_r}\t{stats.n_i}\t{stats.n_s}\t{stats.n_m}\t{stats.n_rr}\t{stats.n_ir}\t{stats.n_sr}\t{stats.n_rs}\t{stats.n_is}\t{stats.n_ss}' , end="")
        if antimicrobial != 'streptomycin':
            print(f'\t{stats.n_ir+stats.n_sr+stats.n_rs}/{stats.n}\t{stats.n_is+stats.n_sr+stats.n_rs}/{stats.n}')
        else:
            print(f'\tNA')

#concordance values
#  *sample
#  *antimicrobial
#
# S,I - no resistance
# R - resistant
#
# collapse S,I
# collapse I,R
#
# N - observations
#
# novel counts
# pheno - Resistant
# pred - not resistant
#
# NI - missing data
class Stats(object):
    def __init__(self):
        self.n_r = 0
        self.n_i = 0
        self.n_s = 0
        self.n_m = 0
        self.n_rr = 0
        self.n_ir = 0
        self.n_sr = 0
        self.n_rs = 0
        self.n_is = 0
        self.n_ss = 0
        self.n = 0

class Sample(object):
    def __init__(self):
        self.id = ''
        self.species = ''
        self.antimicrobials = {}

    def __init__(self, id, species):
        self.id = id
        self.species = species
        self.antimicrobials = {}

    def add_am(self, antimicrobial):
        self.antimicrobials[antimicrobial.id] = antimicrobial

    def print(self):
        print(f'id      : {self.id}')
        print(f'species : {self.species}')
        for key in sorted(self.antimicrobials):
            self.antimicrobials[key].print()

# Antimicrobial	Class	WGS-predicted phenotype	Match	Genetic background
# amikacin	aminoglycoside	Resistant	3	aac(6')-Iaa (aac(6')-Iaa_NC_003197)
# colistin	polymyxin	No resistance	0
# ampicillin+clavulanic acid	beta-lactam	No resistance	0
# pipercallin+tazobactam	NA	NA	NA	Not in database
# ampicillin	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# cefoxitin	beta-lactam	No resistance	0
# ertapenem	beta-lactam	No resistance	0
# tigecycline	tetracycline	No resistance	0
#*chloramphenicol	amphenicol	No resistance	0
#*tetracycline	tetracycline	No resistance	0
# tobramycin	aminoglycoside	Resistant	3	aac(6')-Iaa (aac(6')-Iaa_NC_003197)
#*azithromycin	macrolide	No resistance	0
# ciprofloxacin	quinolone	Resistant	3	qnrB19 (qnrB19_EU432277)
#*nalidixic acid	quinolone	No resistance	0
# ceftazidime	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# imipenem	beta-lactam	No resistance	0
# temocillin	beta-lactam	No resistance	0
#*cefotaxime	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# cefepime	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
#*sulfamethoxazole	folate pathway antagonist	No resistance	0
#*meropenem	beta-lactam	No resistance	0
#*gentamicin	aminoglycoside	No resistance	0
#*trimethoprim	folate pathway antagonist	No resistance	0

# Antimicrobial	Class	WGS-predicted phenotype	Match	Genetic background
# gentamicin	aminoglycoside	No resistance	0
# tobramycin	aminoglycoside	Resistant	3	aac(6')-Iaa (aac(6')-Iaa_NC_003197)
# streptomycin	aminoglycoside	No resistance	0
# amikacin	aminoglycoside	Resistant	3	aac(6')-Iaa (aac(6')-Iaa_NC_003197)
# isepamicin	aminoglycoside	No resistance	0
# dibekacin	aminoglycoside	No resistance	0
# kanamycin	aminoglycoside	No resistance	0
# neomycin	aminoglycoside	No resistance	0
# lividomycin	aminoglycoside	No resistance	0
# paromomycin	aminoglycoside	No resistance	0
# ribostamycin	aminoglycoside	No resistance	0
# unknown aminoglycoside	aminoglycoside	No resistance	0
# butiromycin	aminoglycoside	No resistance	0
# butirosin	aminoglycoside	No resistance	0
# hygromycin	aminoglycoside	No resistance	0
# netilmicin	aminoglycoside	No resistance	0
# apramycin	aminoglycoside	No resistance	0
# sisomicin	aminoglycoside	No resistance	0
# arbekacin	aminoglycoside	No resistance	0
# kasugamycin	aminoglycoside	No resistance	0
# astromicin	aminoglycoside	No resistance	0
# fortimicin	aminoglycoside	No resistance	0
# spectinomycin	aminocyclitol	No resistance	0
# fluoroquinolone	quinolone	No resistance	0
# ciprofloxacin	quinolone	Resistant	3	qnrB19 (qnrB19_EU432277)
# unknown quinolone	quinolone	No resistance	0
# nalidixic acid	quinolone	No resistance	0
# amoxicillin	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# amoxicillin+clavulanic acid	beta-lactam	No resistance	0
# ampicillin	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# ampicillin+clavulanic acid	beta-lactam	No resistance	0
# cefepime	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# cefixime	beta-lactam	No resistance	0
# cefotaxime	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# cefoxitin	beta-lactam	No resistance	0
# ceftazidime	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# ertapenem	beta-lactam	No resistance	0
# imipenem	beta-lactam	No resistance	0
# meropenem	beta-lactam	No resistance	0
# piperacillin	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# piperacillin+tazobactam	beta-lactam	No resistance	0
# unknown beta-lactam	beta-lactam	No resistance	0
# aztreonam	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# cefotaxime+clavulanic acid	beta-lactam	No resistance	0
# temocillin	beta-lactam	No resistance	0
# ticarcillin	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# ceftazidime+avibactam	beta-lactam	No resistance	0
# penicillin	beta-lactam	No resistance	0
# ceftriaxone	beta-lactam	Resistant	3	blaCTX-M-55 (blaCTX-M-55_DQ810789)
# ticarcillin+clavulanic acid	beta-lactam	No resistance	0
# cephalothin	beta-lactam	No resistance	0
# cephalotin	beta-lactam	No resistance	0
# piperacillin+clavulanic acid	beta-lactam	No resistance	0
# ceftiofur	under_development	No resistance	0
# sulfamethoxazole	folate pathway antagonist	No resistance	0
# trimethoprim	folate pathway antagonist	No resistance	0
# fosfomycin	fosfomycin	No resistance	0
# vancomycin	glycopeptide	No resistance	0
# teicoplanin	glycopeptide	No resistance	0
# bleomycin	glycopeptide	No resistance	0
# lincomycin	lincosamide	No resistance	0
# clindamycin	lincosamide	No resistance	0
# dalfopristin	streptogramin a	No resistance	0
# pristinamycin iia	streptogramin a	No resistance	0
# virginiamycin m	streptogramin a	No resistance	0
# quinupristin+dalfopristin	streptogramin a	No resistance	0
# tiamulin	pleuromutilin	No resistance	0
# carbomycin	macrolide	No resistance	0
# erythromycin	macrolide	No resistance	0
# azithromycin	macrolide	No resistance	0
# oleandomycin	macrolide	No resistance	0
# spiramycin	macrolide	No resistance	0
# tylosin	macrolide	No resistance	0
# telithromycin	macrolide	No resistance	0
# tetracycline	tetracycline	No resistance	0
# doxycycline	tetracycline	No resistance	0
# minocycline	tetracycline	No resistance	0
# tigecycline	tetracycline	No resistance	0
# quinupristin	streptogramin b	No resistance	0
# pristinamycin ia	streptogramin b	No resistance	0
# virginiamycin s	streptogramin b	No resistance	0
# linezolid	oxazolidinone	No resistance	0
# chloramphenicol	amphenicol	No resistance	0
# florfenicol	amphenicol	No resistance	0
# colistin	polymyxin	No resistance	0
# fusidic acid	steroid antibacterial	No resistance	0
# mupirocin	pseudomonic acid	No resistance	0
# rifampicin	rifamycin	No resistance	0
# metronidazole	nitroimidazole	No resistance	0
# narasin	ionophores	No resistance	0
# salinomycin	ionophores	No resistance	0
# maduramicin	ionophores	No resistance	0

#specific = {True, False}
#ast_pheno = {R, I, S, NI}
#ast_pred = {R, S}
class Antimicrobial(object):
    def __init__(self):
        self.id = ''
        self.type = ''
        self.specific = ''
        self.ast_pheno = ''
        self.ast_pred = ''

    def __init__(self, id):
        self.id = id
        self.type = ''
        self.specific = False
        self.ast_pheno = ''
        self.ast_pred = ''

    def print(self):
        print(f'\t===================')
        print(f'\tid        : {self.id}')
        print(f'\ttype      : {self.type}')
        print(f'\tast_pheno : {self.ast_pheno}')
        print(f'\tast_pred  : {self.ast_pred}')

if __name__ == '__main__':
    main()
