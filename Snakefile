# from snakemake.utils import min_version
# min_version("6.4.1")
from os.path import join as pjoin


##### config file #####
# configfile: "config/config.yml"

SAMPLE_NAME = config["name"]
REF = config["fa"]
FQ = config["fq"]
HAP1 = config["hap-pat"]
HAP2 = config["hap-mat"]
WD = config["wd"]
TRF = config["trf"]

EASYBED = config["strat"]["easy"]
HARDBED = config["strat"]["hard"]

GIAB11 = config["giab11"]["vcf"]
GIAB11_BED = config["giab11"]["bed"]
GIAB06 = config["giab06"]["vcf"]
GIAB06_BED = config["giab06"]["bed"]

TRUTHS = config["truths"]
CALLERS = config["callers"]


# callers from assembly
include: "rules/callers-asm.smk"
include: "rules/callers-asm-post.smk"
# alignment and phasing
include: "rules/preprocess.smk"
# callers from BAM
include: "rules/callers.smk"
# callsets from severus paper
include: "rules/callers-severus-hg19.smk"
# SVDSS
include: "rules/svdss2.smk"
include: "rules/svdss2-ht.smk"
include: "rules/svdss.smk"
# giab
include: "rules/giab.smk"


# postprocessing
truvari_options = {
    "def": "--passonly --dup-to-ins",
    "nosim": "--passonly --dup-to-ins -p 0",
    "bed": "--passonly --dup-to-ins --includebed " + DIPBED,
    "wbed": "--passonly --dup-to-ins --includebed ",
    "easybed": "--passonly --dup-to-ins --includebed " + EASYBED,
    "hardbed": "--passonly --dup-to-ins --includebed " + HARDBED,
    "sev": "--passonly --typeignore --dup-to-ins -p 0 -s 30 -S 0",
}


include: "rules/truvari.smk"
include: "rules/minda.smk"


rule all:
    input:
        # from callers-asm.smk
        expand(pjoin(WD, "truths", "{truth}.vcf.gz"), truth=TRUTHS),
        # from callers-asm-post.smk
        expand(pjoin(WD, "truths", "{truth}.haps-w500.paf"), truth=TRUTHS),
        expand(
            pjoin(WD, "truths", "comparison-{mode}", "{truth2}-against-{truth1}"),
            mode=["def", "bed"],
            truth1=TRUTHS,
            truth2=TRUTHS,
        ),
        # from truvari.smk
        expand(
            pjoin(WD, "{truth}.truvari-{opt}.csv.png"),
            truth=TRUTHS,
            opt=["def", "nosim", "bed", "easybed", "hardbed"],
        ),
        # from giab.smk
        pjoin(WD, "giab-v1.1-def.csv"),
        pjoin(WD, "giab-v1.1-wbed.csv"),
        pjoin(WD, "giab-v0.6-def.csv") if GIAB06 != "." else [],
        pjoin(WD, "giab-v0.6-wbed.csv") if GIAB06 != "." else [],


rule copy_truth:
    input:
        "data/HG002_grch37.vcf.gz",
    output:
        pjoin(WD, "truths", "severus-paper.vcf.gz"),
    shell:
        """
        cp {input} {output}
        cp {input}.tbi {output}.tbi
        """


# rule plot_truth:
#     input:
#         expand(pjoin(WD, "truths", "{truth}.vcf.gz"), truth=TRUTHS),
#     output:
#         pjoin(WD, "truths", "stats.png"),
#     params:
#         tdir=pjoin(WD, "truths"),
#     conda:
#         "./envs/seaborn.yml"
#     shell:
#         """
#         python3 ./scripts/plot_truth.py {params.tdir}
#         """


rule plot_truvari_all:
    input:
        expand(
            pjoin(WD, "{truth}.truvari-{opt}.csv"), truth=TRUTHS, opt=truvari_options
        ),
    output:
        pjoin(WD, "truvari-all.f1.png"),
    conda:
        "./envs/seaborn.yml"
    shell:
        """
        python3 ./scripts/plot_truvari.py all {WD}
        """


# rule vcf2gz:
#     input:
#         "{fname}.vcf",
#     output:
#         "{fname}.vcf.gz",
#     params:
#         tmp_prefix="{fname}.bcftools-sort-tmp",
#     conda:
#         "./envs/bcftools.yml"
#     shell:
#         """
#         mkdir -p {params.tmp_prefix}
#         bcftools sort -T {params.tmp_prefix} -Oz {input} > {output}
#         tabix -p vcf {output}
#         """
