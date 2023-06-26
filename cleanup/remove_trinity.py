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
from os.path import exists
import subprocess
import argparse
import textwrap
import re
import sys
import filecmp
from pathlib import Path
from shutil import copy2
from shutil import rmtree

def main():

    # options
    parser = CAVSParser(
        description='Removes trinity folders from backups.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
           usage: remove_trinity <dir>

           '''))
    parser.add_argument('dir', help='directory to scan for trinity subfolder folders')
    parser.add_argument('-l', '--level', default=0, type=int, required=False,  help='Level of action, 0 for dry run, 1 for copy, 2 for copy and delete')
    args = parser.parse_args()

    for arg in vars(args):
        print('\t{0:<20} :   {1:<10}'.format(arg, getattr(args, arg)))

    level = args.level

    no_read_partitions = 0
    if not exists('read_partitions_path.txt'):
        #scan all directories
        print("Searching for read_partitions")
        print(f"find {args.dir} -type d -name \"read_partitions\"")
        dirs = subprocess.run(f"find {args.dir} -type d -name \"read_partitions\"", text=True, shell=True, capture_output=True)
        print(dirs.stdout)
        no_read_partitions = 0
        f = open('read_partitions_path.txt', 'w')
        dir_list = dirs.stdout.split('\n')
        for i, path in enumerate(dir_list):
            if i < len(dir_list) - 1:
                f.write(f'{path}\n')
            else:
                f.write(path)
            no_read_partitions += 1
        f.close()
    else:
        no_read_partitions = subprocess.run(f"wc -l read_partitions_path.txt", text=True, shell=True, capture_output=True).stdout

    print(f'Found {no_read_partitions} read partitions')

    no_dirs = 0
    no_dirs_with_fasta = 0
    no_dirs_without_fasta = 0
    total_fasta_size = 0
    total_dir_size = 0
    total_fasta_files_copied = 0
    total_fasta_files_exist = 0
    no_fasta_files_same = 0
    no_dirs_deleted = 0

    file = open('read_partitions_path.txt')
    for line in file:
        if line.endswith('read_partitions\n'):
            dir_path = Path(line.rstrip()[:-15:])
            print(dir_path)
        else:
            print(f'UNEXPECTED DIR ENCOUNTERED : {line}')
            break

        #look for contig file
        fasta_file = None
        for file in os.listdir(dir_path):
            if re.search('\.fasta$', file) is not None:
                #print('\t' + file)
                fasta_file = file

        #move fasta file up one directory
        file_has_been_copied =  False
        no_fasta_file =  False
        if fasta_file is not None:
            no_dirs_with_fasta += 1
            parent_dir_path = dir_path.parent.absolute()
            fasta_size = os.stat(f'{dir_path}/{fasta_file}').st_size
            total_fasta_size += fasta_size

            src_file_path = dir_path.joinpath(fasta_file)
            dst_file_path = parent_dir_path.joinpath(fasta_file)
            dst_dir_path = parent_dir_path
            if os.path.exists(dst_file_path):
                total_fasta_files_exist += 1
                if filecmp.cmp(src_file_path, dst_file_path):
                    no_fasta_files_same += 1
                    file_has_been_copied = True
            else:
                if level == 1:
                    print(f'copy {src_file_path} to {dst_dir_path}')
                    copy2(src_file_path, dst_dir_path)
                    file_has_been_copied = True

                total_fasta_files_copied += 1
        else:
            print(f'NO FASTA: {dir_path}')
            no_fasta_file = True
            no_dirs_without_fasta += 1

        total_dir_size += get_dir_size(dir_path)

        try:
            if level == 2:
                if file_has_been_copied or no_fasta_file:
                    print(f'delete {dir_path}')
                    rmtree(dir_path)
                    no_dirs_deleted += 1
        except OSError as ex:
            print(ex)

    total_fasta_size /= 1000000000
    total_dir_size /= 1000000000

    print(f'===========================')
    print(f'statistics')
    print(f'===========================')
    print(f'#dirs:              {no_dirs}')
    print(f'   #with fasta:     {no_dirs_with_fasta}')
    print(f'   #w/o  fasta:     {no_dirs_without_fasta}')
    print(f'total fasta size:   {total_fasta_size}GB')
    print(f'total dir size:     {total_dir_size}GB')
    print(f'total fasta copied: {total_fasta_files_copied}')
    print(f'total fasta exist:  {total_fasta_files_exist}')
    print(f'   #no_same:        {no_fasta_files_same}')
    print(f'no dir deleted:     {no_dirs_deleted}')

    # for line in dirs.stdout.splitlines():
    #     print(line)

def get_dir_size(path='.'):
    total = 0
    with os.scandir(path) as it:
        for entry in it:
            if entry.is_file():
                total += entry.stat().st_size
            elif entry.is_dir():
                total += get_dir_size(entry.path)
    return total

class CAVSParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

main()