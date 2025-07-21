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
mkdir -p $OD/truvari-def
mkdir -p $OD/truvari-bed
mkdir -p $OD/truvari-nosim


for vcf in $(ls $WD/callsets/*.vcf.gz)
do
    tool=$(basename $vcf .vcf.gz)
    echo $tool
    truvari bench -c $vcf -b $TVCF -o $OD/truvari-def/$tool --passonly --dup-to-ins &
    truvari bench -c $vcf -b $TVCF -o $OD/truvari-bed/$tool --passonly --dup-to-ins --includebed $BED &
    truvari bench -c $vcf -b $TVCF -o $OD/truvari-nosim/$tool --passonly --dup-to-ins --pctseq 0 &
    wait
done
