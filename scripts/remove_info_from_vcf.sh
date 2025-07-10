#!/bin/sh

VCF=$1

bcftools view --header-only $VCF

bcftools view --no-header $VCF | cut -f1-7,9- | while read a b c d e f g h i
do
    echo -e "$a\t$b\t$c\t$d\t$e\t$f\t$g\t.\t$h\t$i"
done


