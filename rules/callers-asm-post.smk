wildcard_constraints:
    h="(1)|(2)",
    mode="(def)|(wbed)",
    reg="()|(-confident)"


# === Pairwise comparison
# ===========================
rule assemblybased_pairwise_truvari:
    input:
        fa=REF,
        vcf1=pjoin(WD, "truths", "{truth1}.vcf.gz"),
        vcf2=pjoin(WD, "truths", "{truth2}.vcf.gz"),
        bed=lambda wildcards: BED if wildcards.truth1 == "hapdiff" or wildcards.truth2 == "hapdiff" else DIPBED,
    output:
        d=directory(pjoin(WD, "truths", "comparison-{mode}", "{truth2}-against-{truth1}")),
        dd=directory(pjoin(WD, "truths", "comparison-{mode}", "{truth2}-against-{truth1}", "phab_bench")),
    params:
        opt=lambda wildcards, input: truvari_options[wildcards.mode] + (" " + input.bed if wildcards.mode == "wbed" else ""),
        aligner="mafft" # lambda wildcards: "poa" if wildcards.mode == "bed" else "mafft",
    conda:
        "../envs/truvari.yml"
    threads: workflow.cores / 2
    shell:
        """
        rm -rf {output.d}
        truvari bench -s 50 -S 50 {params.opt} --reference {input.fa} --base {input.vcf1} --comp {input.vcf2} --output {output.d}
        truvari refine --reference {input.fa} --regions {output.d}/candidate.refine.bed --coords R --use-original-vcfs --threads {threads} --align {params.aligner} {output.d}
        truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine # --with-refine
        """


# === TTmars-like analysis
# ============================


rule remove_info:
    input:
        vcf=pjoin(WD, "truths{reg}", "{asmcaller}.vcf.gz"),
    output:
        vcf=pjoin(WD, "truths{reg}", "{asmcaller}.noinfo.vcf.gz"),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        bash ./scripts/remove_info_from_vcf.sh {input.vcf} | $CONDA_PREFIX/bin/bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule vcf2regions:
    input:
        vcf=rules.remove_info.output.vcf,
    output:
        txt=pjoin(WD, "truths{reg}", "{asmcaller}.noinfo.regions-w{w}.txt"),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        python3 ./scripts/vcf2regions.py {input.vcf} {wildcards.w} > {output.txt}
        """


rule get_contigs:
    input:
        fa=REF,
        vcf=rules.remove_info.output.vcf,
        txt=rules.vcf2regions.output.txt,
    output:
        fa=pjoin(WD, "truths{reg}", "{asmcaller}.hap{h}-w{w}.fa"),
    params:
        tag=lambda wildcards: "first" if wildcards.h == "1" else "second",
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        grep {params.tag} {input.txt} | cut -f1 -d" " | while read region ; do samtools faidx {input.fa} $region | $CONDA_PREFIX/bin/bcftools consensus {input.vcf} -H {wildcards.h} ; done > {output.fa} 2> {output.fa}.log
        """


rule cat_real_contigs:
    input:
        fa1=HAP1,
        fa2=HAP2,
    output:
        fa=pjoin(WD, "real-contigs.fa"),
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        """


rule cat_alternative_contigs:
    input:
        fa1=pjoin(WD, "truths{reg}", "{asmcaller}.hap1-w{w}.fa"),
        fa2=pjoin(WD, "truths{reg}", "{asmcaller}.hap2-w{w}.fa"),
    output:
        fa=pjoin(WD, "truths{reg}", "{asmcaller}.haps-w{w}.fa"),
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        """


rule align_alternative_contigs:
    input:
        tfa=pjoin(WD, "real-contigs.fa"),
        qfa=pjoin(WD, "truths{reg}", "{asmcaller}.haps-w{w}.fa"),
    output:
        paf=pjoin(WD, "truths{reg}", "{asmcaller}.haps-w{w}.paf"),
    conda:
        "../envs/minimap2.yml"
    threads: workflow.cores
    shell:
        """
        minimap2 -t{threads} -c {input.tfa} {input.qfa} > {output.paf}
        """
