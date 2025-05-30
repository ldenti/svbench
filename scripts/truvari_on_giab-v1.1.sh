#!/bin/sh

set -e

WD=$1
TVCF=$2
BED=$3

mkdir -p $WD/giab-v1.1

for vcf in $(ls $WD/callsets/*.vcf.gz)
do
    tool=$(basename $vcf .vcf.gz)
    echo $tool
    truvari bench -c $vcf -b $TVCF -o $WD/giab-v1.1/$tool --passonly --dup-to-ins --includebed $BED
done
