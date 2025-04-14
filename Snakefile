# from snakemake.utils import min_version
# min_version("6.4.1")
from os.path import join as pjoin


##### config file #####
# configfile: "config/config.yml"

SAMPLE_NAME = config["name"]
REF = config["fa"]
FQ = config["fq"]
HAP1 = config["hap1"]
HAP2 = config["hap2"]
WD = config["wd"]
TRF = config["trf"]

TRUTHS = config["truths"]
CALLERS = config["callers"]


# callers from assembly
include: "rules/callers-asm.smk"
# alignment and phasing
include: "rules/preprocess.smk"
# callers from BAM
include: "rules/callers.smk"
# SVDSS
include: "rules/svdss2.smk"
include: "rules/svdss2-ht.smk"
include: "rules/svdss.smk"


# postprocessing
truvari_options = {
    "def": "--passonly --dup-to-ins",
    "nopass": "--dup-to-ins",
    "nosim": "--passonly --dup-to-ins -p 0",
    # "cus":"--passonly --dup-to-ins -s 50 -S 50",
    "bed": "--passonly --dup-to-ins --includebed "
    + DIPBED,  # XXX: DIPBED is actually not an input for truvari rule
    "sev": "--passonly --typeignore --dup-to-ins -p 0 -s 30 -S 0",
}


include: "rules/truvari.smk"
include: "rules/minda.smk"


rule all:
    input:
        # pjoin(WD, "truths", "hapdiff"),
        # from truvari.smk
        expand(
            pjoin(WD, "{truth}.truvari-{opt}.csv"), truth=TRUTHS, opt=truvari_options
        ),
        expand(
            pjoin(WD, "{truth}.truvari-{opt}.csv.png"),
            truth=TRUTHS,
            opt=truvari_options,
        ),
        # expand(pjoin(WD, "{truth}.truvari41-{opt}.csv"), truth=TRUTHS, opt=truvari_options),
        # expand(pjoin(WD, "minda-{truth}", "{caller}"), truth=TRUTHS, caller=CALLERS),


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
