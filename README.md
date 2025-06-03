# svbench

Modular and extensible framework to evaluate SV callers (from long reads) against SV truthsets created from diploid genomes.

The goal of this repository is to simplify SV calling performance evaluation while enabling a fair and standardised evaluation.

This framework is provided as a Snakemake pipeline. Additional python scripts can be used to summarize and plot the results.

### Prerequisites
* conda/mamba
* snakemake
* seaborn (if you want plots)

``` sh
mamba create -c bioconda -c conda-forge -n svbench snakemake
conda activate svbench
# edit config/config.yml
snakemake -c 16 --use-conda --configfile config/config.yml -p [-n]

WD=$(grep "wd:" ./config/config.yml | cut -f2 -d" ")
ls $WD/*.csv $WD/*.png
```

### Other analyses
For now, these scripts assume `truvari` and `seaborn` to be installed.
``` sh
WD=$(grep "wd:" ./config/config.yml | cut -f2 -d" ")

# plot INS/DEL and size distribution of assembly-based callsets
python3 scripts/plot_truth.py distr $t2t_wd/truths $hg38_wd/truths $hg19_wd/truths
# plot GT and neighboring from assembly-based callsets
python3 scripts/plot_truth.py gt $t2t_wd/truths $hg38_wd/truths $hg19_wd/truths

# plot statistics from truth callsets (from single reference run)
python3 scripts/plot_truth_single.py $WD/truths > dipcall-170bp.bed
# intersect dipcall peak around 170bp with Centromeric Satellite Annotation
# https://genome.cse.ucsc.edu/cgi-bin/hgTrackUi?db=hub_3671779_hs1&c=chr12&g=hub_3671779_censat
bedtools intersect -a dipcall-170bp.bed -b censat.bed -wb | cut -f 7 | sort | uniq -c

bash run_truvari.sh $WD/truths

# single genome heatmap
python3 scripts/plot_comparison.py $WD/truths
# 3 genomes heatmap (jaccard similarity)
python3 scripts/plot_comparison.py avga $t2t_wd/truths $hg38_wd/truths $hg19_wd/truths
# 3 genomes heatmap (precision)
python3 scripts/plot_comparison.py avgp $t2t_wd/truths $hg38_wd/truths $hg19_wd/truths

# summarize F1 in a single plot
python3 scripts/plot_truvari.py all $t2t_wd $hg38_wd $hg19_wd
# rank map depending on F1
python3 scripts/plot_truvari.py rank $t2t_wd t2t

# F1 results on GIAB stratification
# uncomment lines 28:32 in scripts/plot_truvari.py , then run:
python3 scripts/plot_truvari.py all $t2t_wd $hg38_wd $hg19_wd

# Compare results on hg19 between "our" callsets and severus callsets
bash ./scripts/build_hg19_comparison_table.sh $hg19_wd
---

# Check how many calls from hapdiff are in dipcall "non-confident" regions (0s in the histogram)
bedtools intersect -a $WD/comparison-def/dipcall-against-hapdiff/fn.vcf.gz -b $WD/truths/dipcall.bed -c | cut -f11 | sort -n | uniq -c

---

# Run last dipcall step on hapdiff .bam files and check
htsbox pileup -q20 -evcf $ref $WD/truths/hapdiff_pat.bam $WD/truths/hapdiff_mat.bam | dipcall-aux.js vcfpair - | bcftools view -v indels -i '(ILEN <= -50 || ILEN >= 50)' | bcftools norm -Oz --multiallelics - > $WD/truths/hapdiff.htsbox.vcf.gz
tabix -p vcf $WD/truths/hapdiff.htsbox.vcf.gz
truvari bench -s 50 -S 50 -c $WD/truths/hapdiff.htsbox.vcf.gz -b $WD/truths/hapdiff.vcf.gz -o OUT
```

### Evaluation against GIAB v1.1
```
# Get *_stvar.vcf.gz and *_stvar.benchmark.bed from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/
bcftools view -Oz -v indels -i '(ILEN <= -50 || ILEN >= 50)' CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz > CHM13v2.0-HG2.sv.vcf.gz
tabix -p vcf CHM13v2.0-HG2.sv.vcf.gz
bash ./scripts/truvari_on_giab-v1.1.sh $WD/ CHM13v2.0-HG2.sv.vcf.gz CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed
python3 ./scripts/format_truvari.py $WD/giab-v1.1/ > $WD/giab-v1.1.csv

python3 ./scripts/plot_pr.py giab-v1.1.t2t.csv giab-v1.1.hg38.csv giab-v1.1.hg19.csv
```

### Supported tools
Read alignment:
* minimap2 (v2.28)

Small variants caller:
* deepvariant (v1.8.0)

Small variants phasing:
* whatshap (v2.5)

SV calling from long reads:
* SVDSS (v1.0.5, v2.1.0)
* sniffles (v2.3)
* cuteSV (v2.1.2)
* debreak (v1.3)
* SVision-pro (v2.4)
* Severus (v1.4.0)
* sawfish (v2.0.0)

SV calling from diploid assemblies:
* dipcall (v0.3)
* svim-asm (v1.0.3)
* hapdiff (commit e0abbb9a8095c70a0de23c49408a530901361b12)

Benchmarkers:
* truvari (v5.2.0)
* minda (commit 47d0fb5484b2b15865a94a9ba81436beaf52cf16)
  * no analysis on minda results
