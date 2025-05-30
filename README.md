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

# summarize F1 in a single plot (assuming max 3 truth callsets)
python3 plot_truvari.py matrix $WD

# plot statistics from truth callsets
python3 scripts/plot_truth.py $WD/truths
bash run_truvari.sh $WD/truths
python3 scripts/plot_comparison.py $WD/truths
```

### Evaluation against GIAB v1.1
```
# Get *_stvar.vcf.gz and *_stvar.benchmark.bed from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/
bcftools view -Oz -v indels -i '(ILEN <= -50 || ILEN >= 50)' CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz > CHM13v2.0-HG2.sv.vcf.gz
tabix -p vcf CHM13v2.0-HG2.sv.vcf.gz
bash ./scripts/truvari_on_giab-v1.1.sh $WD/ CHM13v2.0-HG2.sv.vcf.gz CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed
python3 ./scripts/format_truvari.py $WD/giab-v1.1/ > $WD/giab-v1.1.csv

paste -d"," <(cut -f1,5,6,7 -d"," hg19.csv) <(cut -f5,6,7 -d"," hg38.csv) <(cut -f5,6,7 -d"," t2t.csv)

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

SV calling from diploid assemblies:
* dipcall (v0.3)
* svim-asm (v1.0.3)
* hapdiff (commit e0abbb9a8095c70a0de23c49408a530901361b12)

Benchmarkers:
* truvari (v5.2.0)
* minda (commit 47d0fb5484b2b15865a94a9ba81436beaf52cf16)
  * no analysis on minda results
