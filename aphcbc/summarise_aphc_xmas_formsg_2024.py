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
        "Senior Management": 0,
    }

    participants = []
    non_participants = []

    with open(file_name, "r") as file:
        for line in file:
            line_no += 1
            if line_no > 6:
                line = line.rstrip()
                # print(line)
                tokens = line.split("\t")
                l = len(tokens)
                attendance = tokens[4]
                name = tokens[3]
                section = tokens[5]

                # print(f"0: {tokens[0]}")
                # print(f"1: {tokens[1]}")
                # print(f"2: {tokens[2]}")
                # print(f"3: {tokens[3]}")
                # print(f"4: {tokens[4]}")
                # print(f"5: {tokens[5]}")


                if attendance == "Yes":
                    # print(f"6: {tokens[6]}")
                    # print(f"7: {tokens[7]}")

                    if not (
                        name == "Adrian" or (name == "Adrian Tan" and section == "PHLA")
                    ):
                        wish = tokens[6]
                        section_no[section] += 1
                        email = tokens[7]
                        #print(f"YES\t{name}\t{email}\t{section}\t{wish}")
                        participants.append(Person(name, section, email, wish))
                        no_participating += 1
                else:
                        email = "n/a"
                        wish = "n/a"
                        no_skipping += 1
                        #print(f"NO\t{name}\tn/a\t{section}\tn/a")
                        non_participants.append(Person(name, section, email, wish))


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
        angel = participants[sample[idx]]
        participants[idx].set_angel(angel.name, angel.section, angel.email)
    print("\tsuccessful. ", file=sys.stderr)

    # move committee members (santas) to end of list
    print("Move committee members (santas) to the back\n", file=sys.stderr)
    old_participants = participants.copy()
    participants.clear()
    for person in old_participants:
        if person.angel in (
            "Adrian Tan",
            "Jasmine",
            "Noemi Gesmundo",
            "Loke Li Yan",
            "Brina Chan",
            "Lee Li Pei",
            "Shahreza Darwis",
            "Reg",
        ):
            person.set_committee()
            participants.append(person)
        else:
            participants.insert(0, person)

    output_file_name = re.sub("\.txt$", ".assigned.txt", file_name)
    print(f"Write out to {output_file_name}", file=sys.stderr)
    with open(output_file_name, "w") as f:
        f.write("#angel\t#angel_email\t#mortal\t#mortal_section\t#mortal_wish\n")
        for person in participants:
            f.write(person.print_str_to_angel() + "\n")

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
    print(f"snr mgmt: {section_no['Senior Management']}", file=sys.stderr)


class Person(object):
    def __init__(self, name, section, email, wish):
        self.name = name
        self.section = section
        self.email = email
        self.wish = wish
        self.angel = ""
        self.angel_email = ""
        self.angel_section = ""
        self.committee = "aphc"

    def set_angel(self, name, section, email):
        self.angel = name
        self.angel_section = section
        self.angel_email = email

    def set_committee(self):
        self.committee = "committee"

    def print(self):
        print(f"name    : {self.name}")
        print(f"section : {self.section}")
        print(f"email   : {self.email}")
        print(f"wish    : {self.wish}")
        print(f"angel_name    : {self.angel}")
        print(f"angel_section : {self.angel_section}")
        print(f"angel_email   : {self.angel_email}")
        print(f"committee     : {self.committee}")

    def print_str(self):
        return f"{self.name}\t{self.section}\t{self.email}\t{self.wish}\t{self.angel}\t{self.angel_section}\t{self.angel_email}"

    def print_str_to_angel(self):
        return f"{self.committee}\t{self.angel}\t{self.angel_section}\t{self.angel_email}\t{self.name}\t{self.section}\t{self.email}\t{self.wish}"


if __name__ == "__main__":
    main()  # type: ignore
