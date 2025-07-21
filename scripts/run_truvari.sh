#!/bin/sh

WD=$1
WD1=$2
WD2=$3
giab11=$4
giab06=$5

mkdir -p $WD

dipcall1=$WD1/truths/dipcall.vcf.gz
dipcall2=$WD2/truths/dipcall.vcf.gz

hapdiff1=$WD1/truths/hapdiff.vcf.gz
hapdiff2=$WD2/truths/hapdiff.vcf.gz

svimasm1=$WD1/truths/svim-asm.vcf.gz
svimasm2=$WD2/truths/svim-asm.vcf.gz

if [ "$giab06" = "." ]
then
    names=( "dipcall" "dipcallgiab" "hapdiff" "hapdiffgiab" "svimasm" "svimasmgiab" "giab11" )
    vcfs=( $dipcall1 $dipcall2 $hapdiff1 $hapdiff2 $svimasm1 $svimasm2 $giab11 )
else
    names=( "dipcall" "dipcallgiab" "hapdiff" "hapdiffgiab" "svimasm" "svimasmgiab" "giab11" "giab06" )
    vcfs=( $dipcall1 $dipcall2 $hapdiff1 $hapdiff2 $svimasm1 $svimasm2 $giab11 $giab06 )
fi
n=${#vcfs[@]}

for i1 in $(seq 0 $((n-1)) )
do
    name1=${names[i1]}
    vcf1=${vcfs[i1]}
    for i2 in $(seq $i1 $((n-1)) )
    do
	name2=${names[i2]}
	vcf2=${vcfs[i2]}
	echo $name1 $name2
	echo truvari bench -s 50 -S 50 -c $vcf2 -b $vcf1 -o $WD/$name2-against-$name1 --passonly
	truvari bench -s 50 -S 50 -c $vcf2 -b $vcf1 -o $WD/$name2-against-$name1 --passonly
    done
done
