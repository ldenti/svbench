
rule remove_info:
    input:
        vcf=pjoin(WD, "truths", "{asm}.vcf.gz"),
    output:
        vcf=pjoin(WD, "truths", "{asm}.noinfo.vcf.gz"),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        python3 ../scripts/remove_info_from_vcf.py {input.vcf} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule vcf2regions:
    input:
        vcf=rules.remove_info.output.vcf,
    output:
        txt=pjoin(WD, "truths", "{asm}.noinfo.regions-{w}.txt"),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        python3 ../scripts/vcf2regions.py {input.vcf} {wildcards.w} > {output.txt}
        """


rule get_contigs:
    input:
        fa=REF,
        vcf=pjoin(WD, "truths", "{asm}.vcf.gz"),
        txt=rules.vcf2regions.output.txt,
    output:
        fa=pjoin(WD, "truths", "{asm}.hap{h}-{w}.fa"),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        grep "first" {input.txt} | cut -f1 -d" " | while read region do ; samtools faidx {input.fa} $region | bcftools consensus {input.vcf} -H {wildcards.h} ; done > {output.fa} 2> {output.fa}.log
        """


rule minimap2:
    input:
        hap=HAP1,
        fa=pjoin(WD, "truths", "{asm}.hap{h}-{w}.fa"),
    output:
        paf=pjoin(WD, "truths", "{asm}.hap{h}-{w}.paf"),
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -t2 -c {input.hap} {input.fa} > {output.paf}
        """


rule cat_contigs:
    input:
        fa1=pjoin(WD, "truths", "{asm}.hap1-{w}.paf"),
        fa2=pjoin(WD, "truths", "{asm}.hap2-{w}.paf"),
    output:
        fa=pjoin(WD, "truths", "{asm}.haps}.paf"),
    shell:
        """
        cat {input} > {output}
        """
