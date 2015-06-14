# getHMRs.py
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

#!/usr/bin/python

import os
import sys
import re
import argparse
import yaml
import pwd
import shutil
import glob


#=============================================

def getHMRpaths(self, rootdir, config):
    IN = open(config)
    x = yaml.load(IN)
    fpaths = {}
    for celltype in x['hg19']:
        for proj in x['hg19'][celltype]:
            projdir = os.path.join(rootdir, proj)
            for sample in x['hg19'][celltype][proj]:
                files = glob.glob(os.path.join(projdir, sample, '*/results_hg19', \
                        sample+'*.hmr'))
                if len(files) <= 1:
                    files = glob.glob(os.path.join(projdir, sample, 'results_hg19', \
                            sample+'*.hmr'))
                fnm = '_'.join([proj, sample.split('_')[1]])
                fpaths[fnm] = files
    return fpaths

#=============================================
def main():
    
    #=============================================================
    parser = argparse.ArgumentParser(description = 'getHMRs included in cell type tracks')
    parser.add_argument('-i', '--input', dest='orgfile', \
                        default='/home/xiaojing/codes/HMRbase/organization.txt', \
                        help='the organization file (default: %(default))') 
    parser.add_argument('-o', '--outdir', dest='outdir', \
                        default='/home/xiaojing/hmr', \
                        help='the output dir (default: %(default))')   
    parser.add_argument('-r', '--rootdir', dest='rootdir', \
                        default='/labdata/methylome/public', \
                        help='the output  (default: %(default))')    

    args = parser.parse_args()

    #=============================================================

    fpaths = getHMRpaths(args.rootdir, args.rootdir, args.orgfile)
    for dataset in fpaths:
        for reppath in fpaths[dataset]:
            basenm = os.path.basename(reppath)
            basenm = basenm.split('.')[0]
            rep = basenm.split('_')[-1]
            if rep[0] == 'R':
                outpath = os.path.join(args.outdir, dataset+'_'+rep+'.hmr')
            else:
                outpath = os.path.join(args.outdir, dataset+'.hmr')
            shutil.copyfile(reppath, outpath)


if __name__ == "__main__":
    main()
