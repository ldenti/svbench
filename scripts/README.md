# Analyses

This folder contains all the scripts needed to reproduce the results presented in the manuscript. We ran the snakemake pipeline three times, one per reference genome.

All the scripts assume `truvari` and `seaborn` to be installed.

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
# Get GIAB v1.1
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/{CHM13v2.0,GRCh37,GRCh38}_HG2-T2TQ100-V1.1_stvar.{vcf.gz,vcf.gz.tbi,benchmark.bed}

# run this for all three reference sequences:
bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz > CHM13v2.0-HG2.sv.vcf.gz
tabix -p vcf CHM13v2.0-HG2.sv.vcf.gz

bash ./truvari_on_giab.sh $t2t_wd/ CHM13v2.0-HG2.sv.vcf.gz CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed 11
python3 ./format_truvari.py $t2t_wd/giab-v1.1/truvari-def > $t2t_wd/truvari-def.giab-v11.csv
python3 ./format_truvari.py $t2t_wd/giab-v1.1/truvari-bed > $t2t_wd/truvari-bed.giab-v11.csv
python3 ./format_truvari.py $t2t_wd/giab-v1.1/truvari-nosim > $t2t_wd/truvari-nosim.giab-v11.csv

# Get GIAB v0.6
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.{vcf.gz,vcf.gz.tbi,bed}
bash ./truvari_on_giab.sh $hg19_wd/ HG002_SVs_Tier1_v0.6.vcf.gz HG002_SVs_Tier1_v0.6.bed 06
python3 ./format_truvari.py $WD/giab-v0.6/truvari-def > $hg_19wd/truvari-def.giab-v0.6.csv

python3 ./plot_giab.py truvari-def.giab-v1.1.t2t.csv truvari-def.giab-v1.1.hg38.csv truvari-def.giab-v1.1.hg19.csv truvari-def.giab-v0.6.csv
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
