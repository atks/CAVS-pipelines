#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2023 Adrian Tan <adrian_tan@nparks.gov.sg>
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

import os.path
import click
import re


@click.command()
@click.argument("input_fasta_files", nargs=-1)
@click.option(
    "-w",
    "--write",
    required=False,
    default=False,
    show_default=True,
    help="reference fasta file",
)
def main(input_fasta_files, write):
    """
    Rename FASTA header

    e.g. rename_seq_header.py target.fa
    """
    n = 0
    n_country_annotated = 0
    n_date_annotated = 0

    # @HD	VN:1.6	SO:unknown
    # @PG	ID:basecaller	PN:dorado	VN:0.3.4+5f5cd02	CL:dorado basecaller /usr/local/dorado-0.3.4/models/dna_r9.4.1_e8_hac@v3.3 pod5
    # @PG	ID:samtools	PN:samtools	PP:basecaller	VN:1.16.1	CL:samtools view -h out.sam
    # @RG	ID:244e09e8425500585d8d9ededfb108745cd80b65_dna_r9.4.1_e8_hac@v3.3	PU:FAW17565	PM:Virus-Nanopore-Laptop	DT:2023-06-08T06:42:22.538+00:00	PL:ONT	DS:basecall_model=dna_r9.4.1_e8_hac@v3.3 runid=244e09e8425500585d8d9ededfb108745cd80b65	LB:20230608_fishfarm_edna	SM:20230608_fishfarm_edna

    # QNAME  1	e8749f2b-39c3-4939-bf0f-8bc08f00485c
    # FLAG   2	4
    # RNAME  3	*
    # POS    4	0
    # MAPQ   5	0
    # CIGAR  6	*
    # RNEXT  7	*
    # PNEXT  8	0
    # TLEN   9	0
    # SEQ    10	TATACTTCGTTCAGTTACGTATTGCTAAGGTTAAAAGGATTCATTCCACGGTAACACCAGCACCTCAAAGAACTGCTCCATAGTTCCTGGGCCCTTGCCAGCGCCTCGGCATCAAGCCTTCAGGGATTGCCGTTCTGGAGCCTTTCGCTTGTCAGCATCTGCGCTGCGCGCAAAGCGATCGGGAAGACCTGGTACAAAACCGCACCTGATGCCCACTTCGGAACTTTGTTGTAAGTGTTTATCTTGAAATCATCTCGGTCTGGCTGGTCCAGCTCGTGCAGCCCCGTGGCGTTGGGCCAAAAACTCTCGCCGGATCTCAGTTCGATAAAGAACCGGTAGTACCTCTGGAGTTGTACATGACAATTGGAGCCTCATACCAGCTCCAGCCGTGACAGCGAGGCAAACACCTTTAGCCCGGGCGTGAAGAAGCCCTCACCAGTCGCTCTGCAGATCCAGCACCTTAGCGATTGGGCCAAAGATGAGTGGATCAGGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAACTGA
    # QUAL   11	%&'''()=;:;<55;;11.,+***+))43699?A<;@1***,==<<:)289;<=@@==>8999;::;@:?@B?;00243)&&&)-0358;762(((&&&&$$$%'(34456:9''*01121*)))''&&(&$%&&,,/034468;<333822221333=5476675621235:5***.34))(,(),,..6:;>>>>AIC><;;<>>AAC?=5'&()*811:=;<<<::;:==<9;<<?>@>=008-,54D===<@8888<<<>>?:791-18;:((''))))524667=>>?6&&&&'&&*.9:9;<=@855/:<::*):=>?@?A@880.++244545---49..3('&'9;=?@:::<@@@@7:=4444=<?@@?=?<<;;>9:99<<<@877,+(*1341100.0/-.,-,+(&%%'''%)+))++-632%$$$$%-(()),.112444(%'()+*('')-000=>=<;;=@@F<;;;:99<=8=)).789?@8787<7::;=<6567.,,,/08:<A<?@@==>>=:<=@CA;;;885+,
    # QSCORE 12	qs:i:12
    # DURATIONOFREADINSECONDS  13	du:f:1.22475
    # NOSAMPLESINSIGNALPRIOTOTRIMMING           14	ns:i:4899
    # NOSAMPLESTRIMMEDFROMSTARTOFSIGNAL           15	ts:i:130
    # READMUX    16	mx:i:3
    # CHANNEL    17	ch:i:324
    # STARTTIME  18	st:Z:2023-06-08T06:49:11.381+00:00
    # READNUMBER 19	rn:i:713
    # FILENAME   20	fn:Z:FAW17565_fc68debf_244e09e8_2.pod5
    # scaling midpoint   21	sm:f:75.2258
    # scaling dispersion 22	sd:f:13.4649
    # scaling version    23	sv:Z:med_mad
    # DUPLEX TAG 24	dx:i:0
    ### 1 for duplex reads
    ### 0 simplex reads with no duplex offsprings
    ### -1 simplex reads that have duplex offsprings
    # READGROUP  25	RG:Z:244e09e8425500585d8d9ededfb108745cd80b65_dna_r9.4.1_e8_hac@v3.3

    # https://github.com/nanoporetech/bonito/blob/master/documentation/SAM.md

    for input_fasta_file in input_fasta_files:
        n += 1
        acc = ""
        header = ""
        with open(input_fasta_file) as file:
            for line in file:
                if line.startswith(">"):
                    header = line.strip()
                    acc = header.lstrip(">").split(" ")[0]
                    break
        genbank_file = input_fasta_file.replace("fasta", "seq")
        country = "no country"
        date = "no date"
        year = "no year"
        if os.path.isfile(genbank_file):
            with open(genbank_file) as file:
                for line in file:
                    if line.startswith("ORIGIN"):
                        break
                    else:
                        line = line.strip()
                        m = re.match(r'/country="(.+)"', line)
                        if m:
                            n_country_annotated += 1
                            country = m.group(1)
                            continue
                        m = re.match(r"JOURNAL   Submitted \((.+)\).+", line)
                        if m:
                            n_date_annotated += 1
                            date = m.group(1)
                            year = date.split("-")[2]
                            continue

        print(f"{acc}\t{country}\t{year}\t{header}")

    print(f"no fasta files       : {n}")
    print(f"no country annotated : {n_country_annotated}")
    print(f"no date annotated    : {n_date_annotated}")

    #  JOURNAL   Submitted (15-MAY-2007) Chapman D.A., Microbiology, Institute Of
    #             Animal Health, Pirbright Laboratory, Ash Road, Pirbright, Woking,
    #             Surrey, GU24 0NF, UNITED KINGDOM

    # FEATURES             Location/Qualifiers
    #      source          1..182284
    #                      /organism="African swine fever virus Benin 97/1"
    #                      /mol_type="genomic DNA"
    #                      /isolate="Benin 97/1"
    #                      /db_xref="taxon:443876"
    #                      /country="Benin"


if __name__ == "__main__":
    main()  # type: ignore
