#!/bin/sh

# set -xe

FA=$1
VCF=$2
CONTIGS=$3
WD=$4

mkdir -p $WD

python3 $(dirname $0)/scripts/vcf2regions.py $VCF 500 > $WD/regions.txt

rm -f $WD/hap1.fa $WD/hap2.fa
grep "first" $WD/regions.txt | cut -f1 -d" " | while read region
do
    samtools faidx $FA $region | bcftools consensus $VCF -H 1
done > $WD/haps.fa 2> $WD/hap1.log

grep "second" $WD/regions.txt | cut -f1 -d" " | while read region
do
    samtools faidx $FA $region | bcftools consensus $VCF -H 2
done >> $WD/haps.fa 2> $WD/hap2.log

~/software/minimap2/minimap2 -t4 -c $CONTIGS $WD/haps.fa > $WD/haps.paf
