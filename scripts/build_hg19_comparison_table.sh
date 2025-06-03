#!/bin/sh

WD=$1

for truth in dipcall svim-asm severus-paper
do
    csv=$WD/$truth.truvari-nosim.csv

    echo $truth,SVDSS,$(grep "SVDSS-w4" $csv | cut -f5,6,7 -d","),$(grep "svdss-severus" $csv | cut -f5,6,7 -d",")
    echo $truth,cuteSV,$(grep "cutesv-w4" $csv | cut -f5,6,7 -d","),$(grep "cutesv-severus" $csv | cut -f5,6,7 -d",")
    echo $truth,debreak,$(grep "debreak," $csv | cut -f5,6,7 -d","),$(grep "debreak-severus" $csv | cut -f5,6,7 -d",")
    echo $truth,sniffles,$(grep "sniffles," $csv | cut -f5,6,7 -d","),$(grep "sniffles-severus" $csv | cut -f5,6,7 -d",")
    echo $truth,severus,$(grep "severus-w4," $csv | cut -f5,6,7 -d","),$(grep "severus-severus" $csv | cut -f5,6,7 -d",")
    echo $truth,svision-pro,$(grep "svisionpro-w4" $csv | cut -f5,6,7 -d","),$(grep "svisionpro-severus" $csv | cut -f5,6,7 -d",")
done
