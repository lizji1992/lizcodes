#!/bin/bash

bin=~/bin/bedtools2/bin
promo_ori=~/data/Annotation/2kb_promo_ori.bed
promo=~/data/Annotation/2kb_promo.bed
TSS=~/data/Annotation/TSS_ori.bed
dir=~/data/subHMR/state_stat/0520
group=$1
data=$2

#-----------I/O----------:
pref=$dir/$group/$2
##########################################################

#awk '{print $6}’ $pref > $pref.id

#awk 'BEGIN{t=0; OFS="\t"}{print $1, $2, $3, t, $5, $6; t+=1}' $pref > $pref.num
#awk 'BEGIN{t=0; OFS="\t"}{print $1, int(($2+$3)/2), ints(($2+$3)/2), t, $5, $6; t+=1}' $pref > $pref.mid
echo "File: " $pref

echo "Remove random chrs..."
python rm_random_chrs.py -i $pref -o ${pref}.tmp
echo "Sorting..."
LC_ALL=C; sort -k 1,1 -k 2,2n -k 3,3n -o ${pref}.nordm ${pref}.tmp
rm ${pref}.tmp

echo "Getting promoter intervals ..."
$bin/intersectBed -wb -a $promo -b ${pref}.nordm > $pref.p
echo "Getting non-promoter intervals ..."
$bin/intersectBed -wb -b $promo -a ${pref}.nordm > $pref.np -v

echo "Getting closest TSS ..."
$bin/closestBed -a ${pref}.nordm -b $TSS > $pref.tss
#$bin/closestBed -d -a $pref.mid -b $TSS > $pref.tssmid

