#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2021 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import argparse
import textwrap
import re
import openpyxl


def main():
    cwd = os.getcwd()

    parser = argparse.ArgumentParser(
        description="Generate sample file for dog panel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
           usage: generate_sa_from_dir.py <input_directory> -s <sample_file>
           """
        ),
    )
    parser.add_argument("input_directory", type=str, help="input directory")
    parser.add_argument(
        "-s",
        "--sample_file",
        type=str,
        default="dogs_fsa_csv_files.sa",
        help="sample file to write to",
    )
    args = parser.parse_args()

    for arg in vars(args):
        print("\t{0:<20} :   {1:<10}".format(arg, getattr(args, arg) or ""))

    #######################
    # aggregate sample files
    #######################
    dogs = {}
    dog_files = {}

    for subdir, dirs, files in os.walk(args.input_directory):
        # print(f"SUBDIR\t{subdir}")
        for file in files:
            if re.search("\.fsa\.csv$", file) != None:
                # process one whole STR panel for each dog
                id = file.split("_")[1]
                # print(f"ID: {id}")
                m = re.match("([a-zA-Z]+)\s*(\d+)", id)
                id_norm = m.group(1) + m.group(2).lstrip("0").rjust(4, "0")
                # print(f"ID: {id} => {id_norm}")
                if id_norm in dogs:
                    dogs[id_norm] += 1
                    id_norm += f"_{dogs[id_norm]}"
                else:
                    dogs[id_norm] = 1
                dog_files[id_norm] = os.path.abspath(subdir + os.sep + file)

    sa_file = open(args.sample_file, "w")
    sa_file.write("#sample_id\tfile\n")
    for key in sorted(dog_files.keys()):
        sa_file.write(f'{key}\t"{dog_files[key]}"\n')
    sa_file.close()

    no_files = len(dog_files)
    print(f"No of files: {no_files}")


main()
