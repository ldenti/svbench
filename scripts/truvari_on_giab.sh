#!/bin/sh

set -e

WD=$1
TVCF=$2
BED=$3
MODE=$4

OD=$WD/giab-v1.1
if [[ MODE -eq "06" ]]
then
    OD=$WD/giab-v0.6
fi
mkdir -p $OD

for vcf in $(ls $WD/callsets/*.vcf.gz)
do
    tool=$(basename $vcf .vcf.gz)
    echo $tool
    truvari bench -c $vcf -b $TVCF -o $OD/$tool --passonly --dup-to-ins --includebed $BED
done
