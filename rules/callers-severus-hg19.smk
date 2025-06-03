rule svdss_severus:
    input:
        vcf=pjoin(os.getcwd(), "data", "svdss.vcf.gz"),
    output:
        vcf=pjoin(WD, "callsets", "svdss-severus.vcf.gz"),
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.vcf}.tbi {output.vcf}.tbi
        """


rule cutesv_severus:
    input:
        vcf=pjoin(os.getcwd(), "data", "cuteSV.vcf.gz"),
    output:
        vcf=pjoin(WD, "callsets", "cutesv-severus.vcf.gz"),
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.vcf}.tbi {output.vcf}.tbi
        """


rule debreak_severus:
    input:
        vcf=pjoin(os.getcwd(), "data", "debreak.vcf.gz"),
    output:
        vcf=pjoin(WD, "callsets", "debreak-severus.vcf.gz"),
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.vcf}.tbi {output.vcf}.tbi
        """


rule severus_severus:
    input:
        vcf=pjoin(os.getcwd(), "data", "severus.vcf.gz"),
    output:
        vcf=pjoin(WD, "callsets", "severus-severus.vcf.gz"),
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.vcf}.tbi {output.vcf}.tbi
        """


rule sniffles_severus:
    input:
        vcf=pjoin(os.getcwd(), "data", "sniffles2.vcf.gz"),
    output:
        vcf=pjoin(WD, "callsets", "sniffles-severus.vcf.gz"),
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.vcf}.tbi {output.vcf}.tbi
        """


rule svisionpro_severus:
    input:
        vcf=pjoin(os.getcwd(), "data", "svision.vcf.gz"),
    output:
        vcf=pjoin(WD, "callsets", "svisionpro-severus.vcf.gz"),
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.vcf}.tbi {output.vcf}.tbi
        """
