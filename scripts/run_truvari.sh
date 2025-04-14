#!/bin/sh

TD=$1

mkdir -p $TD/comparison

for truth1 in $TD/*.vcf.gz
do
    name1=$(basename $truth1 .vcf.gz)
    for truth2 in $TD/*.vcf.gz
    do
        name2=$(basename $truth2 .vcf.gz)
        truvari bench -s 50 -S 50 -c $truth2 -b $truth1 -o $TD/comparison/$name2-against-$name1 --passonly &
    done
done

wait