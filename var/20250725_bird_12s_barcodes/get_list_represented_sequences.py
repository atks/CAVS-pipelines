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
import random

@click.command()
@click.option("-s", "--singapore_bird_list_file", required=True, help="singapore bird list file")
@click.option("-g", "--genbank_sequence_list_file", required=True, help="genbank sequence file")
@click.option("-o", "--output_bird_sequence_list_file", required=True, help="sequence list file")
def main(singapore_bird_list_file, genbank_sequence_list_file, output_bird_sequence_list_file): 
    """
    Extract 12Sand COI sequences from the Singapore Bird List and GenBank sequence list files.   

    e.g. get_list_represented_sequences.py -s singapore_bird_list.txt -g genbank_sequence_list.txt
    """
    print("\t{0:<20} :   {1:<10}".format("singapore bird list", singapore_bird_list_file))
    print("\t{0:<20} :   {1:<10}".format("genbank sequencelist", genbank_sequence_list_file))
    
    # read sample file
    bird_species = {}
    with open(singapore_bird_list_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                idx, common_name, species_name, rarity, status = line.rstrip().split("\t")
                bird_species[species_name] = BirdSpecies(idx, common_name, species_name, rarity, status)

    print(f"Total bird species found: {len(bird_species)}")
   
    no_seq = 0
    no_bird_seq = 0
    with open(output_bird_sequence_list_file, "w") as ofile:
        with open(genbank_sequence_list_file, "r") as file:
            for line in file:
                if not line.startswith("#"):
                    acc, title, length, tax_id, kingdom, species, common_name, common_scientific_namesc = line.rstrip().split("\t")
                    if species in bird_species:
                        ofile.write(line)
                        no_bird_seq += 1
                    no_seq += 1
                    if no_seq % 100000 == 0:
                        print(f"Processed {no_seq} sequences so far with {no_bird_seq} matching bird species...")
class BirdSpecies(object):
    def __init__(self, idx, common_name, species_name, rarity, status):
        self.idx = idx
        self.common_name = common_name
        self.species_name = species_name
        self.rarity = rarity
        self.status = status



if __name__ == "__main__":
    main() # type: ignore
