rule sniffles:
    input:
        fa=REF,
        bam=pjoin(WD, "{aligner}", "whatshap-haplotagged.bam"),
        bed=TRF,
    output:
        vcf=pjoin(WD, "{aligner}", "sniffles", "variations.vcf"),
    params:
        supp=lambda wildcards: "2" if wildcards.support == "2" else "auto",
    conda:
        "../envs/sniffles.yml"
    threads: THREADS
    shell:
        """
        sniffles --phase --threads {THREADS} --reference {input.fa} --input {input.bam} --vcf {output.vcf} --tandem-repeats {input.bed} --mapq 20 --minsvlen 50 --minsupport {params.support}
        """
