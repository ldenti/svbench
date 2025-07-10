# TTmars-like analysis


wildcard_constraints:
    h="(1)|(2)",


rule remove_info:
    input:
        vcf=pjoin(WD, "truths", "{asm}.vcf.gz"),
    output:
        vcf=pjoin(WD, "truths", "{asm}.noinfo.vcf.gz"),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        bash ./scripts/remove_info_from_vcf.sh {input.vcf} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule vcf2regions:
    input:
        vcf=rules.remove_info.output.vcf,
    output:
        txt=pjoin(WD, "truths", "{asm}.noinfo.regions-w{w}.txt"),
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
        fa=pjoin(WD, "truths", "{asm}.hap{h}-w{w}.fa"),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        grep "first" {input.txt} | cut -f1 -d" " | while read region do ; samtools faidx {input.fa} $region | bcftools consensus {input.vcf} -H {wildcards.h} ; done > {output.fa} 2> {output.fa}.log
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


rule cat_artificial_contigs:
    input:
        fa1=pjoin(WD, "truths", "{asm}.hap1-w{w}.fa"),
        fa2=pjoin(WD, "truths", "{asm}.hap2-w{w}.fa"),
    output:
        fa=pjoin(WD, "truths", "{asm}.haps-w{w}.fa"),
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        """


rule minimap2:
    input:
        tfa=pjoin(WD, "real-contigs.fa"),
        qfa=pjoin(WD, "truths", "{asm}.haps-w{w}.fa"),
    output:
        paf=pjoin(WD, "truths", "{asm}.haps-w{w}.paf"),
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -t2 -c {input.tfa} {input.qfa} > {output.paf}
        """
