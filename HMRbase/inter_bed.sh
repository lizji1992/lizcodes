f=$1


intersectBed=~/panfs/bin/bedtools2/bin/intersectBed
annot=~/data/Annotation/2kb_promo.bed
in=~/panfs/subhmr/subhmr/$f


$intersectBed -wb -a $annot -b $in > $in.p
$intersectBed -wb -a $annot -b $in > $in.tss
