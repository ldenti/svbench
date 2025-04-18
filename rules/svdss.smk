rule svdss_smooth:
    input:
        fa=REF,
        bam=BAM,
    output:
        bam=pjoin(WD, "SVDSS", "smoothed.selective.bam"),
    params:
        wd=pjoin(WD, "SVDSS"),
    threads: workflow.cores
    log:
        pjoin(WD, "times", "svdss-smooth.time"),
    conda:
        "../envs/svdss.yml"
    shell:
        """
        SVDSS smooth --reference {input.fa} --bam {input.bam} --workdir {params.wd} --threads {threads}
        samtools index {output.bam}
        """


rule svdss_search:
    input:
        fmd=REF + ".fmd",
        bam=pjoin(WD, "SVDSS", "smoothed.selective.bam"),
    output:
        sfs=pjoin(WD, "SVDSS", "solution_batch_0.assembled.sfs"),
    params:
        wd=pjoin(WD, "SVDSS"),
    threads: workflow.cores
    log:
        pjoin(WD, "times", "svdss-search.time"),
    conda:
        "../envs/svdss.yml"
    shell:
        """
        SVDSS search --index {input.fmd} --bam {input.bam} --threads {threads} --workdir {params.wd} --assemble
        """


rule svdss_call:
    input:
        fa=REF,
        bam=pjoin(WD, "SVDSS", "smoothed.selective.bam"),
        sfs=pjoin(WD, "SVDSS", "solution_batch_0.assembled.sfs"),
    output:
        vcf=pjoin(WD, "callsets", "SVDSS-w{w}.vcf.gz"),
    params:
        wd=pjoin(WD, "SVDSS"),
    threads: workflow.cores
    log:
        pjoin(WD, "times", "svdss-call-w{w}.time"),
    conda:
        "../envs/svdss.yml"
    shell:
        """
        n=$(ls {params.wd}/solution_batch_*.assembled.sfs | wc -l)
        SVDSS call --reference {input.fa} --bam {input.bam} --threads {threads} --workdir {params.wd} --batches ${{n}} --min-cluster-weight {wildcards.w}
        bgzip -c {params.wd}/svs_poa.vcf > {output.vcf}
        tabix -p vcf {output.vcf}
        """
