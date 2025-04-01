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

import csv
import click
from shutil import copy2
from openpyxl import Workbook
from openpyxl import load_workbook

@click.command()
@click.option(
    "-t",
    "--txt_file",
    show_default=True,
    help="txt file",
)
@click.option(
    "-x",
    "--xlsx_file",
    show_default=True,
    help="xlsx file",
)
def main(txt_file, xlsx_file):
    """
    Convert excel file into text file

    e.g. xlsx_to_txt.py -t poe.txt -x poe_results.xlsx
    """
    print("\t{0:<20} :   {1:<10}".format("text file", txt_file))
    print("\t{0:<20} :   {1:<10}".format("excel file", xlsx_file))

    #read in text file
    with open(txt_file, mode='w') as file:
        file.write(f"investigator\tgroup\tdog_id\tvideo_id\tframe_id\tmethod\tla\tao\n")
        wb = load_workbook(xlsx_file, data_only=True)
        ws = wb['For calculation']
        for row in ws.iter_rows():
            if row[ws["A1"].column - 1].value == "Investigator":
                continue

            investigator = row[ws["A1"].column - 1].value
            group = row[ws["B1"].column - 1].value
            if group == "Self-studied vet student":
                group = "Self-study vet student"

            dog_vid = row[ws["C1"].column - 1].value
            if dog_vid == "MU-653-10-08":
                dog_vid = "MU-653-01-08"
            dog_id = "-".join(dog_vid.split("-")[0:3])
            video_id = dog_vid.split("-")[3]
            r1_la = row[ws["D1"].column - 1].value
            r1_ao = row[ws["E1"].column - 1].value
            h1_la = row[ws["F1"].column - 1].value
            h1_ao = row[ws["G1"].column - 1].value
            r2_la = row[ws["H1"].column - 1].value
            r2_ao = row[ws["I1"].column - 1].value
            h2_la = row[ws["J1"].column - 1].value
            h2_ao = row[ws["K1"].column - 1].value
            r3_la = row[ws["L1"].column - 1].value
            r3_ao = row[ws["M1"].column - 1].value
            h3_la = row[ws["N1"].column - 1].value
            h3_ao = row[ws["O1"].column - 1].value

            if r1_la != h1_la or r2_la != h2_la or r3_la != h3_la:
                print(f"Warning: {investigator} {group} {dog_id} {video_id} has different la values")
                print(f"\t{r1_la} {h1_la}")
                print(f"\t{r2_la} {h2_la}")
                print(f"\t{r3_la} {h3_la}")

            file.write(f"{investigator}\t{group}\t{dog_id}\tv{video_id}\tf1\trishniw\t{r1_la}\t{r1_ao}\n")
            file.write(f"{investigator}\t{group}\t{dog_id}\tv{video_id}\tf1\thansson\t{h1_la}\t{h1_ao}\n")
            file.write(f"{investigator}\t{group}\t{dog_id}\tv{video_id}\tf2\trishniw\t{r2_la}\t{r2_ao}\n")
            file.write(f"{investigator}\t{group}\t{dog_id}\tv{video_id}\tf2\thansson\t{h2_la}\t{h2_ao}\n")
            file.write(f"{investigator}\t{group}\t{dog_id}\tv{video_id}\tf3\trishniw\t{r3_la}\t{r3_ao}\n")
            file.write(f"{investigator}\t{group}\t{dog_id}\tv{video_id}\tf3\thansson\t{h3_la}\t{h3_ao}\n")

if __name__ == "__main__":
    main() # type: ignore
