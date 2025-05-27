#!/bin/sh

TD=$1  # truth directory

mkdir -p $TD/comparison-def

for truth1 in $TD/*.vcf.gz
do
    name1=$(basename $truth1 .vcf.gz)
    for truth2 in $TD/*.vcf.gz
    do
        name2=$(basename $truth2 .vcf.gz)
        truvari bench -s 50 -S 50 -c $truth2 -b $truth1 -o $TD/comparison-def/$name2-against-$name1 --passonly &
    done
done

wait

mkdir -p $TD/comparison-bed

for truth1 in $TD/*.vcf.gz
do
    name1=$(basename $truth1 .vcf.gz)
    for truth2 in $TD/*.vcf.gz
    do
        name2=$(basename $truth2 .vcf.gz)
        truvari bench -s 50 -S 50 -c $truth2 -b $truth1 -o $TD/comparison-bed/$name2-against-$name1 --includebed $TD/dipcall.bed &
    done
done

wait
