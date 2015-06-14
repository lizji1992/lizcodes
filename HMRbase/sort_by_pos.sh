#!/bin/bash


#-----------I/O----------:
nm=$1
pref=~/panfs/subhmr/subhmr/$nm


##########################################################


##########################################################

LC_ALL=C; sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o ${pref}.subhmr.sorted ${pref}.subhmr
echo Sort: Done!

