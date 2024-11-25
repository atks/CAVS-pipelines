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

import os
import click
import random
import math
import re
import sys


@click.command()
@click.argument("file_name")
def main(file_name):
    """
    summarise aphc xmas 2024 formSG results

    e.g. summarise_aphc_xmas_formsg.py
    """
    # read phrase file
    phrases = []
    line_no = 0
    no_participating = 0
    no_skipping = 0
    section_no = {
        "CAR": 0,
        "CWR": 0,
        "VAR": 0,
        "VM": 0,
        "VFP": 0,
        "PQA": 0,
        "PBPS": 0,
        "PHLA": 0,
    }

    participants = []
    non_participants = []

    with open(file_name, "r") as file:
        for line in file:
            line_no += 1
            if line_no > 6:
                line = line.rstrip()
                # print(line)
                tokens = line.split(",")
                l = len(tokens)
                attendance = tokens[4]
                name = tokens[3]
                section = tokens[-1]
                if attendance == "Yes":
                    if not (
                        name == "Adrian" or (name == "Adrian Tan" and section == "PHLA")
                    ):
                        wish = ",".join(tokens[5:-2])
                        wish = wish.rstrip('"').lstrip('"')
                        section_no[section] += 1
                        email = tokens[-2]
                        participants.append(Person(name, section, email, wish))
                        # print(f"YES\t{name}\t{email}\t{section}\t{wish}")
                        no_participating += 1
                else:
                    if name == "Ng Oi Wing":
                        section = "VAR"
                    if name == "Amy Chan":
                        section = "VFP"
                    if not (name == "Tara" or name == "Kum Chew"):
                        no_skipping += 1
                        non_participants.append(Person(name, section, email, wish))
                        # print(f"NO\t{name}\tn/a\t{section}\tn/a")

    print("Randomising ", end="", file=sys.stderr)
    # randomise assignments

    no_participants = len(participants)
    print(f"{no_participants} assigments ... ", file=sys.stderr)

    sample = []

    while True:
        sample = random.sample(range(no_participants), no_participants)
        # check self assignment
        collision_detected = False
        for idx, i in enumerate(sample):
            if idx == i:
                collision_detected = True
                print("\tCollision detected, resampling", file=sys.stderr)
                break
        if not collision_detected:
            break

    for idx in range(no_participants):
        santa = participants[sample[idx]]
        participants[idx].set_santa(santa.name, santa.section, santa.email)
    print("\tsuccessful. ", file=sys.stderr)

    # move committee members (santas) to end of list
    print("Move committee members (santas) to the back\n", file=sys.stderr)
    old_participants = participants.copy()
    participants.clear()
    for person in old_participants:
        if person.santa in (
            "Adrian Tan",
            "Jasmine",
            "Noemi Gesmundo",
            "Loke Li Yan",
            "Brina",
            "Li Pei",
            "Shahreza Darwis",
            "Reg",
        ):
            person.set_committee()
            participants.append(person)
        else:
            participants.insert(0, person)

    # tweaks
    # affixed: kum chew randomly sampled mdm tay
    print(f"======Swap KC's santee======", file=sys.stderr)

    idx_santa_kc = 0
    idx_mdm_tay = 0
    for idx, person in enumerate(participants):
        if person.santa == "Kum Chew":
            idx_santa_kc = idx
        if person.name == "Tay Yih Hong":
            idx_mdm_tay = idx

    # copy madam tay's santa to idx_santa_kc
    # participants[idx].set_santa(santa.name, santa.section, santa.email)

    participants[idx_mdm_tay].print()
    print(f"+++++++++++++++++++++++", file=sys.stderr)
    participants[idx_santa_kc].print()

    print(f"=======After swap======", file=sys.stderr)

    participants[idx_santa_kc].set_santa(
        participants[idx_mdm_tay].santa,
        participants[idx_mdm_tay].santa_section,
        participants[idx_mdm_tay].santa_email,
    )
    participants[idx_mdm_tay].set_santa(
        "Kum Chew", "VFP", "hiong_kum_chew@nparks.gov.sg"
    )
    participants[idx_mdm_tay].print()
    print(f"+++++++++++++++++++++++", file=sys.stderr)
    participants[idx_santa_kc].print()

    output_file_name = re.sub("\.csv$", ".assigned.txt", file_name)
    print(f"Write out to {output_file_name}", file=sys.stderr)
    with open(output_file_name, "w") as f:
        f.write("#santa\t#santa_email\t#santee\t#santee_section\t#santee_wish\n")
        for person in participants:
            f.write(person.print_str_to_santa() + "\n")

    print(f"\n\n", file=sys.stderr)
    print(f"Summary", file=sys.stderr)
    print(f"No respondants:   {no_participating+no_skipping}", file=sys.stderr)
    print(f"No participating: {no_participating}", file=sys.stderr)
    print(f"No skipping:      {no_skipping}", file=sys.stderr)

    print(f"car:      {section_no['CAR']}", file=sys.stderr)
    print(f"cwr:      {section_no['CWR']}", file=sys.stderr)
    print(f"var:      {section_no['VAR']}", file=sys.stderr)
    print(f"vm:       {section_no['VM']}", file=sys.stderr)
    print(f"vfp:      {section_no['VFP']}", file=sys.stderr)
    print(f"pqa:      {section_no['PQA']}", file=sys.stderr)
    print(f"pbps:     {section_no['PBPS']}", file=sys.stderr)
    print(f"phla:     {section_no['PHLA']}", file=sys.stderr)


class Person(object):
    def __init__(self, name, section, email, wish):
        self.name = name
        self.section = section
        self.email = email
        self.wish = wish
        self.santa = ""
        self.santa_email = ""
        self.santa_section = ""
        self.committee = "aphc"

    def set_santa(self, name, section, email):
        self.santa = name
        self.santa_section = section
        self.santa_email = email

    def set_committee(self):
        self.committee = "committee"

    def print(self):
        print(f"name    : {self.name}")
        print(f"section : {self.section}")
        print(f"email   : {self.email}")
        print(f"wish    : {self.wish}")
        print(f"santa_name    : {self.santa}")
        print(f"santa_section : {self.santa_section}")
        print(f"santa_email   : {self.santa_email}")
        print(f"committee     : {self.committee}")

    def print_str(self):
        return f"{self.name}\t{self.section}\t{self.email}\t{self.wish}\t{self.santa}\t{self.santa_section}\t{self.santa_email}"

    def print_str_to_santa(self):
        return f"{self.committee}\t{self.santa}\t{self.santa_section}\t{self.santa_email}\t{self.name}\t{self.section}\t{self.email}\t{self.wish}"


if __name__ == "__main__":
    main()  # type: ignore
