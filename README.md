# svbench

Modular and extensible framework to evaluate SV callers (from long reads) against SV truthsets created from diploid genomes.

The goal of this repository is to simplify SV calling performance evaluation while enabling a fair and standardised evaluation.

This framework is provided as a Snakemake pipeline. Additional python scripts can be used to summarize and plot the results.

### Prerequisites
* conda/mamba
* snakemake
* seaborn (if you want plots)

``` sh
mamba create -c bioconda -c conda-forge -n svbench snakemake-minimal seaborn
conda activate svbench
# edit config/config.yml
snakemake -c 16 --use-conda --configfile config/config.yml -p [-n]

WD=$(grep "wd:" ./config/config.yml | cut -f2 -d" ")
ls $WD/*.csv $WD/*.png
```

### Analyses
These scripts assume `truvari` and `seaborn` to be installed. Moreover, everything is hardcoded to 3 reference sequences (so please run the snakemake pipeline three times).

##### Assembly-based
``` sh
# Figure on assembly-based callsets
python3 scripts/plot_asm.py $t2t_wd $hg38_wd $hg19_wd

# Three genomes heatmap (jaccard similarity/accuracy)
python3 scripts/plot_comparison_asm.py $t2t_wd/truths $hg38_wd/truths $hg19_wd/truths
```

##### Read-based vs assembly-based
``` sh
# summarize F1 in a single plot
python3 scripts/plot_truvari.py all $t2t_wd $hg38_wd $hg19_wd
# rank map depending on F1
python3 scripts/plot_truvari.py rank $t2t_wd t2t

# F1 results on GIAB stratification
python3 scripts/plot_giabstrat.py $t2t_wd $hg38_wd $hg19_wd
```

##### Evaluation against GIAB (v1.1 and v0.6)
``` sh
# Get GIAB v1.1
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/{CHM13v2.0,GRCh37,GRCh38}_HG2-T2TQ100-V1.1_stvar.{vcf.gz,vcf.gz.tbi,benchmark.bed}
# run this for all three reference sequences:
bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz > CHM13v2.0-HG2.sv.vcf.gz
tabix -p vcf CHM13v2.0-HG2.sv.vcf.gz
bash ./scripts/truvari_on_giab.sh $t2t_wd/ CHM13v2.0-HG2.sv.vcf.gz CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed 11
python3 ./scripts/format_truvari.py $t2t_wd/giab-v1.1/ > $t2t_wd/giab-v1.1.csv

# Get GIAB v0.6
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.{vcf.gz,vcf.gz.tbi,bed}
bash ./scripts/truvari_on_giab.sh $hg19_wd/ HG002_SVs_Tier1_v0.6.vcf.gz HG002_SVs_Tier1_v0.6.bed 06
python3 ./scripts/format_truvari.py $WD/giab-v0.6/ > $hg_19wd/giab-v0.6.csv

python3 ./scripts/plot_giab.py giab-v1.1.t2t.csv giab-v1.1.hg38.csv giab-v1.1.hg19.csv giab-v0.6.csv
```

### Data
- Reference sequences and annotations:
```
# hg19
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
wget https://github.com/fritzsedlazeck/Sniffles/blob/master/annotations/human_hs37d5.trf.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh37@all/Union/GRCh37_alldifficultregions.bed.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh37@all/Union/GRCh37_notinalldifficultregions.bed.gz

# hg38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/refs/heads/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/Union/GRCh38_notinalldifficultregions.bed.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz

# t2t
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/CHM13@all/Union/CHM13_alldifficultregions.bed.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz
```

- HG002 assembly (from HPRC):
```
curl -o HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc
agc getset ../../HPRC-yr1.agc HG002.1 > HG002.1.fa
agc getset ../../HPRC-yr1.agc HG002.2 > HG002.2.fa
```

- HG002 HiFi sample: `m64012_190920_173625`, `m64012_190921_234837`, `m64015_190920_185703`, and `m64015_190922_010918` (from [HPRC](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/15kb/))


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
