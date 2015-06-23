#!/bin/bash

bin=~/codes/bedtools2/bin
HMRbase=~/codes/lizcodes/HMRbase
promo_ori=~/data/Annotation/2kb_promo_ori.bed
promo=~/data/Annotation/2kb_promo.bed
TSS=~/data/Annotation/TSS_ori.bed

#-----------I/O----------:
pref=$1
#state2id=$2
##########################################################

#if [ $state2id == "y" ]
#then
#awk 'BEGIN{t=0; OFS="\t"}{print $1, int(($2+$3)/2), ints(($2+$3)/2), t, $5, $6; t+=1}' $pref > $pref.mid
#fi

echo "Assign IDs to intervals...: " $pref
#awk 'BEGIN{t=0; OFS="\t"}{print $1, $2, $3, t, $5; t+=1}' $pref > $pref.withid
awk 'BEGIN{t=0; OFS="\t"}{print t, $4; t+=1}' $pref > $pref.id
echo "Remove random chrs..."
#python $HMRbase/pyscripts/rm_random_chrs.py -i $pref.withid -o ${pref}.tmp
#rm ${pref}.withid
#echo "Sorting..."
#LC_ALL=C; sort -k 1,1 -k 2,2n -k 3,3n -o ${pref}.nordm ${pref}.tmp
#rm ${pref}.tmp

#echo "Getting promoter intervals ..."
#$bin/intersectBed -wb -a $promo -b ${pref}.nordm > $pref.p
#echo "Getting non-promoter intervals ..."
#$bin/intersectBed -wb -b $promo -a ${pref}.nordm > $pref.np -v
#
#echo "Getting closest TSS ..."
#$bin/closestBed -a ${pref}.nordm -b $TSS > $pref.tss
##$bin/closestBed -d -a $pref.mid -b $TSS > $pref.tssmid
