# Analyses

This folder contains all the scripts needed to reproduce the results presented in the manuscript. All the scripts assume `truvari` and `seaborn` to be installed (e.g., via conda).

We ran the Snakemake pipeline three times, one per reference genome (i.e., hg19, hg38, t2t). For each reference genome, update the [configuration `config/config.yml`](https://github.com/ldenti/svbench/blob/main/config/config.yml) and run the Snakemake pipeline. We assume the variables `{hg19,hg38,t2t}_smk_wd` to be the path to the output directories of the three Snakemake runs (one per reference sequence).

#### Assembly-based analysis
The folder `truth` in the snakemake output directory contains the assembly-based truth sets. They can be analyzed using the following commands:
``` sh
# Figure on assembly-based callsets
python3 ./plot_asm.py $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd

# Three genomes heatmap (jaccard similarity/accuracy)
python3 ./plot_comparison_asm.py $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd
```

#### Read-based vs assembly-based
``` sh
# summarize F1 in a single plot
python3 ./plot_truvari.py all $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd
# rank map depending on F1
python3 ./plot_truvari.py rank $t2t_smk_wd T2T

# F1 results on GIAB stratification
python3 ./plot_giabstrat.py $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd
```

#### Evaluation against GIAB (v1.1 and v0.6)
``` sh
python3 ./plot_giab.py $t2t_smk_wd/giab-v1.1-def.csv $hg38_smk_wd/giab-v1.1-def.csv $hg19_smk_wd/giab-v1.1-def.csv $hg19_smk_wd/giab-v0.6-def.csv
# python3 ./plot_giab.py $t2t_smk_wd/giab-v1.1-bed.csv $hg38_smk_wd/giab-v1.1-bed.csv $hg19_smk_wd/giab-v1.1-bed.csv $hg19_smk_wd/giab-v0.6-bed.csv
```

<!--
## Double assembly analyses
These scripts analyze the results obtained from both assemblies:
1. run the snakemake on 3 references using HPRC contigs
2. run the snakemake on 3 references using GIAB haplotypes (you can symlink the callsets directory)
3. run the GIAB v1.1 scripts (see above) from one of the two runs (since this is independent from the assembly used)

#### Pairwise comparison of truth sets
This will compare all assembly-based truth sets and GIAB v1.1 and v0.6 (on hg19).
``` sh
bash ./run_truvari.sh t2t.output_directory t2t.smk_workdir_on_hprc t2t.smk_workdir_on_giab hg38.giab-v11.vcf.gz .
bash ./run_truvari.sh hg38.output_directory hg38.smk_workdir_on_hprc hg38.smk_workdir_on_giab hg38.giab-v11.vcf.gz .
bash ./run_truvari.sh hg19.output_directory hg19.smk_workdir_on_hprc hg19.smk_workdir_on_giab hg19.giab-v11.vcf.gz hg19.giab-v06.vcf.gz

# plot the accuracy heatmap
python3 ./plot_comparison_asm_full.py t2t.output_directory hg38.output_directory hg19.output_directory
```

#### Full rankmap
This will produce a rankmap containing both assemblies and the GIAB v1.1 truth set.
``` sh
python3 scripts/plot_rankmap.py hg38.smk_workdir_on_hprc hg38.smk_workdir_on_giab PlotTitle
# ^ adapt for other references
```
-->

## Data
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

- HG002 assembly (from GIAB):
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz
# split fasta in maternal and paternal using samtools faidx
```

- HG002 HiFi sample: `m64012_190920_173625`, `m64012_190921_234837`, `m64015_190920_185703`, and `m64015_190922_010918` (from [HPRC](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/15kb/))

- GIAB curated SV callsets:
```
# Get GIAB v0.6
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.{vcf.gz,vcf.gz.tbi,bed}

# Get GIAB v1.1
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/{CHM13v2.0,GRCh37,GRCh38}_HG2-T2TQ100-V1.1_stvar.{vcf.gz,vcf.gz.tbi,benchmark.bed}

# run this for all three reference sequences:
bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz > CHM13v2.0-HG2.sv.vcf.gz
tabix -p vcf CHM13v2.0-HG2.sv.vcf.gz
```
