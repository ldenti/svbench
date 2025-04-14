# svbench

Modular and extensible framework to evaluate SV callers (from long reads) against SV truthsets created from diploid genomes.

The goal of this repository is to
* simplify SV calling performance evaluation
* enable fair and standardised benchmarks of SV callers

This framework is provided as a Snakemake pipeline and several python scripts to analyze its output (and create plots).

### Supported tools
SV calling from long reads:
* SVDSS (v1.0.5)
* SVDSS (v2.1.0)
* sniffles (v2.3)
* cuteSV (v2.1.2)
* debreak (v1.3)
* SVision-pro (v2.4)
* Severus (v1.4.0)

SV calling from diploid assemblies:
* dipcall (v0.3)
* svim-asm (v1.0.3)
* hapdiff (commit: e0abbb9a8095c70a0de23c49408a530901361b12)
