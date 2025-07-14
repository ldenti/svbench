rule truvari:
    input:
        vcf=pjoin(WD, "callsets", "{caller}.vcf.gz"),
        truth=pjoin(WD, "truths", "{truth}.vcf.gz"),
        bed=DIPBED,
    output:
        directory(pjoin(WD, "truvari-{truth}-{option}", "{caller}")),
    params:
        opt=lambda wildcards: truvari_options[wildcards.option],
    conda:
        "../envs/truvari.yml"
    shell:
        """
        truvari bench -c {input.vcf} -b {input.truth} -o {output} {params.opt}
        """


rule format_truvari:
    input:
        expand(pjoin(WD, "truvari-{{truth}}-{{option}}", "{caller}"), caller=CALLERS),
    output:
        pjoin(WD, "{truth}.truvari-{option}.csv"),
    params:
        bd=pjoin(WD, "truvari-{truth}-{option}"),
    shell:
        """
        python3 ./scripts/format_truvari.py {params.bd} > {output}
        """


rule plot_truvari:
    input:
        pjoin(WD, "{truth}.truvari-{option}.csv"),
    output:
        pjoin(WD, "{truth}.truvari-{option}.csv.png"),
    conda:
        "../envs/seaborn.yml"
    shell:
        """
        python3 ./scripts/plot_truvari.py single {input}
        """
