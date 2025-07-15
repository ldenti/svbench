DIPBED = pjoin(WD, "truths", "dipcall", "prefix.dip.bed")


rule dipcall:
    input:
        fa=REF,
        # bed=PARBED,
        hap1=HAP1,  # paternal
        hap2=HAP2,  # maternal
    output:
        vcf=pjoin(WD, "truths", "dipcall", "prefix.dip.vcf.gz"),
        bed=pjoin(WD, "truths", "dipcall", "prefix.dip.bed"),  # DIPBED
        bam1=pjoin(WD, "truths", "dipcall", "prefix.hap1.bam"),
        bam2=pjoin(WD, "truths", "dipcall", "prefix.hap2.bam"),
    params:
        wdir=pjoin(WD, "truths", "dipcall"),
        prefix=pjoin(WD, "truths", "dipcall", "prefix"),
        mak=pjoin(WD, "truths", "dipcall.mak"),
    threads: workflow.cores
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        run-dipcall -t {threads} {params.prefix} {input.fa} {input.hap1} {input.hap2} > {params.mak}
        mkdir -p {params.wdir}
        make -j 2 -f {params.mak}
        tabix -p vcf {output.vcf}
        samtools index {output.bam1}
        samtools index {output.bam2}
        """


rule clean_dipcall:
    input:
        fai=REF + ".fai",
        bed=pjoin(WD, "truths", "dipcall", "prefix.dip.bed"),
        vcf=pjoin(WD, "truths", "dipcall", "prefix.dip.vcf.gz"),
    output:
        vcf=pjoin(WD, "truths", "dipcall.vcf.gz"),
        bed=pjoin(WD, "truths", "dipcall.bed"),
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools reheader --fai {input.fai} {input.vcf} | bcftools norm --multiallelics - | bcftools view -v indels -i '(ILEN <= -30 || ILEN >= 30)' -Oz > {output.vcf}
        tabix -p vcf {output.vcf}
        cp {input.bed} {output.bed}
        """


######################################################################
######################################################################
######################################################################

# rule merge_bam:
#     input:
#         bam1=pjoin(WD, "dipcall", "prefix.hap1.bam"),
#         bam2=pjoin(WD, "dipcall", "prefix.hap2.bam"),
#     output:
#         bam=pjoin(WD, "dipcall", "prefix.dip.bam"),
#     conda:
#         "../envs/samtools.yml"
#     shell:
#         """
#         samtools merge {output.bam} {input.bam1} {input.bam2}
#         samtools index {output.bam}
#         """


# rule cutesv_asm:
#     input:
#         fa=REF,
#         bam=pjoin(WD, "dipcall", "prefix.dip.bam"),
#     output:
#         vcf=pjoin(WD, "cutesv-asm.vcf"),
#     params:
#         tmp=pjoin(WD, "cutesv-asm.tmp"),
#     conda:
#         "../envs/cutesv.yml"
#     shell:
#         """
#         mkdir -p {params.tmp}
#         cuteSV {input.bam} {input.fa} {output.vcf}.pre {params.tmp} -s 1 --genotype --report_readid -p -1 -mi 500 -md 500 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
#         python3 $CONDA_PREFIX/lib/python3.6/site-packages/cuteSV/diploid_calling.py {output.vcf}.pre {output.vcf}
#         """

######################################################################
######################################################################
######################################################################


rule svim_asm:
    input:
        fa=REF,
        bam1=pjoin(WD, "truths", "dipcall", "prefix.hap1.bam"),
        bam2=pjoin(WD, "truths", "dipcall", "prefix.hap2.bam"),
    output:
        vcf=pjoin(WD, "truths", "svim-asm", "variants.vcf"),
    params:
        wd=pjoin(WD, "truths", "svim-asm"),
    conda:
        "../envs/svimasm.yml"
    shell:
        """
        mkdir -p {params.wd}
        svim-asm diploid {params.wd} {input.bam1} {input.bam2} {input.fa}
        """


rule clean_svim_asm:
    input:
        vcf=pjoin(WD, "truths", "svim-asm", "variants.vcf"),
    output:
        vcf=pjoin(WD, "truths", "svim-asm.vcf.gz"),
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


######################################################################
######################################################################
######################################################################


rule get_hapdiff:
    output:
        exe=pjoin(WD, "software", "hapdiff", "hapdiff.py"),
        repo=directory(pjoin(WD, "software", "hapdiff")),
    shell:
        """
        git clone https://github.com/KolmogorovLab/hapdiff {output.repo}
        cd {output.repo}
        git checkout e0abbb9a8095c70a0de23c49408a530901361b12
        git submodule update --init
        make
        """


rule hapdiff:
    input:
        exe=pjoin(WD, "software", "hapdiff", "hapdiff.py"),
        fa=REF,
        hap1=HAP1,
        hap2=HAP2,
    output:
        vcf=pjoin(WD, "truths", "hapdiff", "hapdiff_phased.vcf.gz"),
    params:
        outd=pjoin(WD, "truths", "hapdiff"),
    conda:
        "../envs/hapdiff.yml"
    threads: workflow.cores
    shell:
        """
        {input.exe} --reference {input.fa} --pat {input.hap1} --mat {input.hap2} --out-dir {params.outd} -t {threads}
        """


rule hapdiff_post:
    input:
        vcf=pjoin(WD, "truths", "hapdiff", "hapdiff_phased.vcf.gz"),
    output:
        vcf=pjoin(WD, "truths", "hapdiff.vcf.gz"),
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """
