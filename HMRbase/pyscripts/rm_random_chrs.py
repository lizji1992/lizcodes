#!/usr/usc/python/2.7.8/bin/python2.7

# tabSubHMR.py
#
# Copyright (C) 2013-2017 University of Southern California and Liz Ji
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import re
import argparse
import pwd
import shutil
import glob

CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', \
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
        'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chrX', 'chrY', 'chrM']

#=============================================

def main():
    
    #=============================================================
    parser = argparse.ArgumentParser(description = 'del all random chrs')
    parser.add_argument('-i', '--input', dest='infile', help='input file')
    parser.add_argument('-o', '--output', dest='outfile', help='output file')

    args = parser.parse_args()

    #=============================================================

    fpin = os.path.join(args.infile)
    fpout = os.path.join(args.outfile)
    IN = open(fpin)
    OUT = open(fpout, 'w')
    for line in IN:
        items = line.rstrip().split('\t')
        if items[0] in CHRS:
            print >>OUT, line.rstrip()

if __name__ == "__main__":
    main()
