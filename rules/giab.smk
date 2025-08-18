wildcard_constraints:
    giab="(1.1)|(0.6)",
    giabopt="(def)|(wbed)",


rule truvari_on_giab_11:
    input:
        fa=REF,
        vcf=pjoin(WD, "callsets", "{caller}.vcf.gz"),
        truth=GIAB11,
        bed=GIAB11_BED,
    output:
        d=directory(pjoin(WD, "giab-v1.1", "truvari-{giabopt}", "{caller}")),
        dd=directory(pjoin(WD, "giab-v1.1", "truvari-{giabopt}", "{caller}", "phab_bench")),
        ddd=pjoin(WD, "giab-v1.1", "truvari-{giabopt}", "{caller}", "ga4gh_with_refine.summary.json"),
    params:
        bed=lambda wildcards: "--includebed "+ GIAB11_BED if wildcards.giabopt == "wbed" else "",
        aligner="mafft" # lambda wildcards: "poa" if wildcards.giabopt == "wbed" else "mafft",
    conda:
        "../envs/truvari.yml"
    threads: workflow.cores / 2
    shell:
        """
        rm -rf {output.d}
        truvari bench --reference {input.fa} {params.bed} --base {input.truth} --comp {input.vcf} --output {output.d} --passonly --pick ac --dup-to-ins
        truvari refine --reference {input.fa} --regions {output.d}/candidate.refine.bed --coords R --use-original-vcfs --threads {threads} --align {params.aligner} {output.d}
        truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine # --with-refine
        """

rule truvari_on_giab_06:
    input:
        fa=REF,
        vcf=pjoin(WD, "callsets", "{caller}.vcf.gz"),
        truth=GIAB06,
        bed=GIAB06_BED,
    output:
        d=directory(pjoin(WD, "giab-v0.6", "truvari-{giabopt}", "{caller}")),
        dd=directory(pjoin(WD, "giab-v0.6", "truvari-{giabopt}", "{caller}", "phab_bench")),
        ddd=pjoin(WD, "giab-v0.6", "truvari-{giabopt}", "{caller}", "ga4gh_with_refine.summary.json"),
    params:
        bed=lambda wildcards: "--includebed "+ GIAB06_BED if wildcards.giabopt == "wbed" else "",
        aligner="mafft" # lambda wildcards: "poa" if wildcards.giabopt == "wbed" else "mafft",
    conda:
        "../envs/truvari.yml"
    shell:
        """
        rm -rf {output.d}
        truvari bench --reference {input.fa} {params.bed} --base {input.truth} --comp {input.vcf} --output {output.d} --passonly --pick ac --dup-to-ins
        truvari refine --reference {input.fa} --regions {output.d}/candidate.refine.bed --coords R --use-original-vcfs --threads {threads} --align {params.aligner} {output.d}
        truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine # --with-refine
        """


rule format_truvari_on_giab:
    input:
        expand(
            pjoin(WD, "giab-v{{giab}}", "truvari-{{giabopt}}", "{caller}"),
            caller=CALLERS,
        ),
    output:
        csv=pjoin(WD, "giab-v{giab}-{giabopt}.csv"),
        refcsv=pjoin(WD, "giab-v{giab}-{giabopt}.refine.csv"),
    params:
        bd=pjoin(WD, "giab-v{giab}", "truvari-{giabopt}"),
    shell:
        """
        python3 ./scripts/format_truvari.py {params.bd} > {output.csv}
        python3 ./scripts/format_truvari.py --refine {params.bd} > {output.refcsv}
        """
