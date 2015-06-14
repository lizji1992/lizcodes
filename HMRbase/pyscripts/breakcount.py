# bed2tab.py
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


#=============================================

def getSampleVec(IDmapfile):
    samples = []
    IN = open(IDmapfile)
    for line in IN:
        items = line.rstrip().split('\t')
        fnm = items[1]
        if (fnm not in samples):
            samples.append(fnm)
    return samples


def checkExistence(occurs, samples):
    vec = []
    for s in samples:    
        if s in occurs: 
           vec.append('1')
        else:
           vec.append('0')
    return vec


def getSubHMRtab(bedfile, samples, cutoff):
    IN = open(bedfile)
    subhmrID = 1
    subtab = {}    
    for line in IN:
        items = line.rstrip().split('\t')
        freq = float(items[4])
        if freq >= cutoff:
            setID = 'subHMR_%d' % (subhmrID)
            subtab[setID] = {'freq':str(freq)}
            subtab[setID]['vec'] = checkExistence(items[3].split(','), samples)
            subhmrID = subhmrID + 1
    return subtab

def printHMRtab(outfile, samples, subtab):
    OUT = open(outfile, 'w')
    header = ' '.join(['id', 'freq'] + samples)
    print >>OUT, header

    for subhmr in subtab.keys():
        print >>OUT, ' '.join([subhmr, subtab[subhmr]['freq']] +\
                subtab[subhmr]['vec'])    

#=============================================
def main():
    
    #=============================================================
    parser = argparse.ArgumentParser(description = 'bed to table')
    parser.add_argument('-i', '--input', dest='bedfile', required=True, \
                        help='BED file')
    parser.add_argument('-i', '--input', dest='bedfile', required=True, \
                          help='BED file')
    parser.add_argument('-o', '--output', dest='outfile', \
                        default='/panfs/cmb-panasas2/xiaojinj/subhmr/subhmr.tab', \
                        help='the output file (default: %(default))')   
    parser.add_argument('-c', '--cutoff', dest='cutoff', \
                        default=0.5, help='cutoff (default: %(default))')    

    args = parser.parse_args()

    #=============================================================
    
    samples = getSampleVec(IDmapfile)

    subtab =  getSubHMRtab(args.bedfile, samples, args.cutoff)
    
    printHMRtab(args.outfile, samples, subtab)


if __name__ == "__main__":
    main()
