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
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import click
import docx

@click.command()
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
@click.option(
    "-d",
    "--data_dir",
    required = True,
    show_default=True,
    help="data directory",
)
@click.option("-s", "--sample_file", required=True, help="sample file")
def main(working_dir, data_dir):
    """
    Combine data from 3 BCP files and generate various reports.

    e.g. generate_bcp_reports.py
    """
    reports_dir = f"{working_dir}/reports"
    print("\t{0:<20} :   {1:<10}".format("working dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("data dir", data_dir))
    print("\t{0:<20} :   {1:<10}".format("reports dir", reports_dir))

    # create directories
    try:
        os.makedirs(reports_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")





class BCP(object):
    def __init__(self, bcp_dir):
        self.bcp_dir = bcp_dir

    # read bcp files
    # BCP_CAVS Laboratory Tests_24 Feb 2025.xlsx
    # CAVS_Lab_combined.xlsx
    # bcp_files/CAVS\ FormSG\ sample\ submission\ spreadsheets
    def initialise_bcp_files(self):
        print(f"initializing bcp files")
 
    def generate_reports(self):
        print(f"Generating reports for BCP:")




if __name__ == "__main__":
    main() # type: ignore
