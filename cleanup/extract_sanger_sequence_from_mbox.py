#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2023 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import click
import filecmp
from pathlib import Path
from shutil import copy2
from shutil import rmtree


@click.command()
@click.option('-s', '--src_dir', default=os.getcwd(), show_default=True, help='working directory')
@click.option('-d', '--dst_dir', default='/data/var/backup/categories/sanger', show_default=True, help='working directory')
def main(src_dir, dst_dir):
    """
    Looks for all .ab1 and .seq files in source directory and copies it to destination directory

    e.g. copy_sanger_sequences
    """

    ab1_files = dict()
    seq_files = dict()

    no_ab1 = 0
    no_seq = 0

    no_ab1_x = 0
    no_ab1_o = 0
    no_seq_x = 0
    no_seq_o = 0

    for dirpath, dirnames, filenames in os.walk(src_dir):
        for file in filenames:
            # print(file)
            # print(f'\t{dirpath}/{file}')
            if file.endswith('.ab1'):
                no_ab1 += 1
                ab1 = os.path.basename(file)
                if ab1 in ab1_files.keys():
                    if ab1_files[ab1].update(f'{dirpath}/{file}'):
                        no_ab1_o += 1
                    else:
                        no_ab1_x += 1
                else:
                    ab1_files[ab1] = Record(f'{dirpath}/{file}')
            elif file.endswith('.seq'):
                no_seq += 1
                seq = os.path.basename(file)
                if seq in seq_files.keys():
                    if seq_files[seq].update(f'{dirpath}/{file}'):
                        no_seq_o += 1
                    else:
                        no_seq_x += 1
                else:
                    seq_files[seq] = Record(f'{dirpath}/{file}')

    ab1_collisions = [0 for i in range(25)]
    seq_collisions = [0 for i in range(25)]
    for key in ab1_files.keys():
        ab1_collisions[len(ab1_files[key].collisions)] += 1
    for key in seq_files.keys():
        seq_collisions[len(seq_files[key].collisions)] += 1

    for key in ab1_files.keys():
        src = ab1_files[key].file_path
        dst = f'{dst_dir}/{ab1_files[key].base}'
        # copy2(src, dst)
        ab1_files[key].print_tab()
        for i, file in enumerate(ab1_files[key].collisions):
            src = file
            dst = f'{dst_dir}/{ab1_files[key].base}.{i+1}'
            # copy2(src, dst)

    for key in seq_files.keys():
        src = seq_files[key].file_path
        dst = f'{dst_dir}/{seq_files[key].base}'
        # copy2(src, dst)
        seq_files[key].print_tab()
        for i, file in enumerate(seq_files[key].collisions):
            src = file
            dst = f'{dst_dir}/{seq_files[key].base}.{i+1}'
            # copy2(src, dst)

    print(f'ab1 files observed  : {no_ab1}')
    print(f'unique ab1 files    : {len(ab1_files.keys())}')
    print(f'seq files observed  : {no_seq}')
    print(f'unique seq files    : {len(seq_files.keys())}')
    print(f'no_ab1_o            : {no_ab1_o}')
    print(f'no_ab1_x            : {no_ab1_x}')
    print(f'no_seq_o            : {no_seq_o}')
    print(f'no_seq_x            : {no_seq_x}')
    print(f'ab1 collisions      : {ab1_collisions}')
    print(f'seq collisions      : {seq_collisions}')


class Record(object):
    def __init__(self, file_path):
        self.base = os.path.basename(file_path)
        # self.base = self.base[:-4]
        self.n = 1
        self.file_path = file_path
        self.collisions = []

    def update(self, file_path):
        if filecmp.cmp(self.file_path, file_path):
            self.n += 1
            return True
        else:
            # print(f'{self.file_path} vs {file_path} not the same')
            self.collisions.append(file_path)
            return False

    def print_tab(self):
        print(f'{self.base}\t{self.n}\t{self.file_path}{":".join(self.collisions)}')

    def print(self):
        print(f'++++++++++++++++++++')
        print(f'base       : {self.base}')
        print(f'n          : {self.n}')
        print(f'path       : {self.file_path}')
        print(f'++++++++++++++++++++')


main()
