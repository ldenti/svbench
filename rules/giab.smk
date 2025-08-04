wildcard_constraints:
    giab="(1.1)|(0.6)",
    option="(def)|(wbed)",


rule truvari_on_giab_11:
    input:
        vcf=pjoin(WD, "callsets", "{caller}.vcf.gz"),
        truth=GIAB11,
        bed=GIAB11_BED,
    output:
        directory(pjoin(WD, "giab-v1.1", "truvari-{option}", "{caller}")),
    params:
        opt=lambda wildcards: truvari_options[wildcards.option]
        + (GIAB11_BED if wildcards.option == "wbed" else ""),
    conda:
        "../envs/truvari.yml"
    shell:
        """
        truvari bench -c {input.vcf} -b {input.truth} -o {output} {params.opt}
        """


rule truvari_on_giab_06:
    input:
        vcf=pjoin(WD, "callsets", "{caller}.vcf.gz"),
        truth=GIAB06,
        bed=GIAB06_BED,
    output:
        directory(pjoin(WD, "giab-v0.6", "truvari-{option}", "{caller}")),
    params:
        opt=lambda wildcards: truvari_options[wildcards.option]
        + (GIAB06_BED if wildcards.option == "wbed" else ""),
    conda:
        "../envs/truvari.yml"
    shell:
        """
        truvari bench -c {input.vcf} -b {input.truth} -o {output} {params.opt}
        """


rule format_truvari_on_giab:
    input:
        expand(
            pjoin(WD, "giab-v{{giab}}", "truvari-{{option}}", "{caller}"),
            caller=CALLERS,
        ),
    output:
        pjoin(WD, "giab-v{giab}-{option}.csv"),
    params:
        bd=pjoin(WD, "giab-v{giab}", "truvari-{option}"),
    shell:
        """
        python3 ./scripts/format_truvari.py {params.bd} > {output}
        """
