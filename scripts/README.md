``` sh

python3 plot_truvari.py rank ~/data/severus-comparison/results/hg38 hg38

---

bedtools intersect -a ~/data/severus-comparison/truths/t2t/comparison-def/dipcall-against-hapdiff/fn.vcf.gz -b ~/data/severus-comparison/truths/t2t/dipcall.bed -c | cut -f11 | sort -n | uniq -c

---

htsbox pileup -q20 -evcf ~/data/palss/chm13v2.0.fa ~/data/severus-comparison/truths/t2t/hapdiff_pat.bam ~/data/severus-comparison/truths/t2t/hapdiff_mat.bam | dipcall-aux.js vcfpair - | bcftools view -v indels -i '(ILEN <= -50 || ILEN >= 50)' | bcftools norm -Oz --multiallelics - > ~/data/severus-comparison/truths/t2t/hapdiff.htsbox.vcf.gz
tabix -p vcf ~/data/severus-comparison/truths/t2t/hapdiff.htsbox.vcf.gz

truvari bench -s 50 -S 50 -c dipcall.vcf.gz -b hapdiff.vcf.gz -o comparison-def/dipcall-against-hapdiff-all
```
