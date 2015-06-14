#! /bin/bash

#//////////////////////////////////////////////////////////////////////////////
#compress=".all.meth .mr"

#//////////////////////////////////////////////////////////////////////////////

/usr/usc/python/2.7.8/bin/python2.7 AutoMeth.py --mode ${mode} -p ${projdir} \
	-m ${mrdir} -F ${projconfig} -s ${skip} -M ${methpipe} -U ${ucsc} \
	-a ${assembly} --chrsize ${chromsize} --compress ${compress}
