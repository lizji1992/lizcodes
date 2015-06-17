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
import argparse

################################################################################
#-------------------------------------------------------------------------------
def print_large_intervals(fin, fout, size):
    IN = open(fin)
    OUT = open(fout, 'w')

    size_int = int(size)

    for line in IN:
        line = line.rstrip()
        items = line.split('\t')
	start = items[1]
	end = items[2]
        if (int(end) - int(start)) > size_int:
            print >>OUT, line


#==============================================================================
def main():
    parser = argparse.ArgumentParser(description = 'BED filter')
    parser.add_argument('-m', '--mode', dest='mode', nargs='+', \
                        default=['size_greater'], \
                        help='the modes you want to run with BED filter')
    parser.add_argument('-s', '--size', dest='size', required=False, \
                        help="size threshold")
    parser.add_argument('-i', '--input', dest='fin', required=True, \
                        help="input file")
    parser.add_argument('-o', '--output', dest='fout', required=True, \
                        help="output file")
 #   parser.add_argument('--chrsize2', dest='chrsize_lift', nargs='*',\
 #           help="chromsize file for ucsc trackbuilding - liftover assembly")
 #   parser.add_argument('-c', '--chain', dest='chain', \
 #           help="ucsc chain file for liftover")

   
    args = parser.parse_args()
    
    if 'size_greater' in args.mode:
        print_large_intervals(args.fin, args.fout, args.size)
    

#==============================================================================
if __name__ == '__main__':
  main()
