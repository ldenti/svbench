"""
* Note 1: I tried to run all tools with their default parameters or
setting the same value for all (e.g., min mapq, min sv length...).
Minimum support is tested when the tool does not provide a way
to compute it automatically.
* Note 2: I tried to use output VCFs as they are. But sometimes,
truvari was complaining about BND records, so I decided to remove them.
"""


rule cutesv:
    input:
        fa=REF,
        bam=BAM,
        bed=TRF,
    output:
        vcf=pjoin(WD, "cutesv-{w}", "cutesv-w{w}.full.vcf"),
    params:
        tmp=pjoin(WD, "cutesv-tmp-w{w}"),
    log:
        time=pjoin(WD, "times", "cutesv-w{w}.time"),
    conda:
        "../envs/cutesv.yml"
    threads: workflow.cores
    shell:
        """
        mkdir -p {params.tmp}
        /usr/bin/time -vo {log.time} cuteSV --threads {threads} --min_support {wildcards.w} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 {input.bam} {input.fa} {output.vcf} {params.tmp} --min_size 50
        # --genotype --min_mapq 20
        """


rule cutesv_post:
    input:
        vcf=pjoin(WD, "cutesv-{w}", "cutesv-w{w}.full.vcf"),
    output:
        vcf=pjoin(WD, "callsets", "cutesv-w{w}.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        # bgzip -c {input.vcf} > {output.vcf}
        # bcftools view -Oz --exclude "INFO/SVTYPE='BND' | INFO/SVTYPE='INS' | INFO/SVTYPE='DUP'" {input.vcf} > {output.vcf}
        bcftools view -Oz --exclude "INFO/SVTYPE='BND'" {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


######################################################################
######################################################################
######################################################################


rule sniffles:
    input:
        fa=REF,
        bam=BAM,
        bed=TRF,
    output:
        vcf=pjoin(WD, "callsets", "sniffles.full.vcf"),
    log:
        time=pjoin(WD, "times", "sniffles.time"),
    conda:
        "../envs/sniffles.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} sniffles --threads {threads} --reference {input.fa} --input {input.bam} --vcf {output.vcf} --tandem-repeats {input.bed}
        # --mapq 20 --minsupport auto --phase --minsvlen 50
        """


rule sniffles_post:
    input:
        vcf=pjoin(WD, "callsets", "sniffles.full.vcf"),
    output:
        vcf=pjoin(WD, "callsets", "sniffles.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        # bgzip -c {input.vcf} > {output.vcf}
        # bcftools view -Oz --include "INFO/SVTYPE='DEL' | INFO/SVTYPE='INS' | INFO/SVTYPE='DUP'" {input.vcf} > {output.vcf}
        bcftools view -Oz --exclude "INFO/SVTYPE='BND'" {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


######################################################################
######################################################################
######################################################################


rule severus:
    input:
        fa=REF,
        bam=BAM,
        bed=TRF,
        vcf=pjoin(WD, "deepvariant-phased.vcf.gz"),
    output:
        vcf=pjoin(WD, "severus-w{w}", "all_SVs", "severus_all.vcf"),
    params:
        tmp=pjoin(WD, "severus-w{w}"),
    log:
        time=pjoin(WD, "times", "severus-w{w}.time"),
    conda:
        "../envs/severus.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} severus --min-mapq 20 --min-support {wildcards.w} --target-bam {input.bam} --vntr-bed {input.bed} --out-dir {params.tmp} -t {threads} --phasing-vcf {input.vcf}
        # --min-sv-size 50
        """


rule severus_post:
    input:
        vcf=pjoin(WD, "severus-w{w}", "all_SVs", "severus_all.vcf"),
    output:
        vcf=pjoin(WD, "callsets", "severus-w{w}.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        # bgzip -c {input.vcf} > {output.vcf}
        # bcftools view -Oz --include "INFO/SVTYPE='DEL' | INFO/SVTYPE='INS' | INFO/SVTYPE='DUP'" {input.vcf} > {output.vcf}
        bcftools view -Oz --exclude "INFO/SVTYPE='BND'" {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


######################################################################
######################################################################
######################################################################


rule debreak:
    input:
        bam=BAM,
    output:
        vcf=pjoin(WD, "debreak", "debreak.vcf"),
    params:
        odir=pjoin(WD, "debreak"),
    log:
        time=pjoin(WD, "times", "debreak.time"),
    conda:
        "../envs/debreak.yml"
    threads: workflow.cores
    shell:
        """
        debreak -t {threads} --bam {input.bam} -o {params.odir} --min_size 50 --min_quality 20 --aligner minimap2
        # --min_support 2 # it is estimated from depth, so I will let him decide
        """


rule debreak_post:
    input:
        vcf=pjoin(WD, "debreak", "debreak.vcf"),
    output:
        vcf=pjoin(WD, "callsets", "debreak.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf}
        # bcftools view -Oz --include "INFO/SVTYPE='DEL' | INFO/SVTYPE='INS' | INFO/SVTYPE='DUP'" {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


######################################################################
######################################################################
######################################################################

# svision-pro v2.4 from conda wasn't working
# I had some conflicts between modules compiled using NumPy 1.x that cannot be run in NumPy 2.0.2


rule get_svisionpro:
    output:
        exe=pjoin(WD, "software", "svision-pro", "SVision-pro"),
        repo=directory(pjoin(WD, "software", "svision-pro")),
        # env=pjoin(WD, "software", "svision-pro", "environment.fixed.yml"),
        model=pjoin(
            WD,
            "software",
            "svision-pro",
            "src",
            "pre_process",
            "model_liteunet_256_8_16_32_32_32.pth",
        ),
    shell:
        """
        rm -rf {output.repo}
        git clone https://github.com/songbowang125/SVision-pro.git {output.repo}
        cd {output.repo}
        git checkout f19eb78bb8430e2cc9df49cde276c59e67560c35
        chmod +x {output.exe}
        """


rule svisionpro:
    input:
        exe=pjoin(WD, "software", "svision-pro", "SVision-pro"),
        # model=pjoin(WD, "software", "svision-pro", "src", "pre_process", "model_liteunet_256_8_16_32_32_32.pth"),
        # ph=pjoin(WD, "software", "svision-pro", "smk-alldone"),
        fa=REF,
        bam=BAM,
        model=pjoin(
            WD,
            "software",
            "svision-pro",
            "src",
            "pre_process",
            "model_liteunet_256_8_16_32_32_32.pth",
        ),
    output:
        vcf=pjoin(WD, "svisionpro-{w}", f"{SAMPLE_NAME}.svision_pro_v2.4.s{{w}}.vcf"),
    params:
        odir=pjoin(WD, "svisionpro-{w}"),
    log:
        time=pjoin(WD, "times", "svisionpro-{w}.time"),
    conda:
        "../envs/svisionpro.yml"  # pjoin(WD, "software", "svision-pro", "environment.fixed.yml"),
    threads: workflow.cores
    shell:
        """
        {input.exe} --target_path {input.bam} --genome_path {input.fa} --model_path {input.model} --out_path {params.odir} --sample_name {SAMPLE_NAME} --detect_mode germline --process_num {threads} --min_supp {wildcards.w} --preset hifi
        # --min_mapq 20 --min_sv_size 20
        """


rule svisionpro_post:
    input:
        vcf=pjoin(WD, "svisionpro-{w}", f"{SAMPLE_NAME}.svision_pro_v2.4.s{{w}}.vcf"),
    output:
        vcf=pjoin(WD, "callsets", "svisionpro-w{w}.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        # bgzip -c {input.vcf} > {output.vcf}
        # bcftools view -Oz --include "INFO/SVTYPE='DEL' | INFO/SVTYPE='INS' | INFO/SVTYPE='DUP'" {input.vcf} > {output.vcf}
        bcftools sort {input.vcf} | bcftools view -Oz --exclude "INFO/SVTYPE='BND'" > {output.vcf}
        tabix -p vcf {output.vcf}
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


######################################################################
######################################################################
######################################################################


rule sawfish:
    input:
        fa=REF,
        bam=BAM,
    output:
        odir=directory(pjoin(WD, "sawfish")),
    log:
        time=pjoin(WD, "times", "sawfish-discover.time"),
    conda:
        "../envs/sawfish.yml"
    threads: workflow.cores
    shell:
        """
        sawfish discover --min-indel-size 50 --min-sv-mapq 20 --threads 16 --ref {input.fa} --bam {input.bam} --output-dir {output.odir}
        # --min-sv-mapq 10
        """


# --expected-cn ${DISTRO_ROOT_DIR}/data/expected_cn/expected_cn.hg38.XX.bed \
# --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz


rule sawfish_jointcall:
    input:
        odir=pjoin(WD, "sawfish"),
    output:
        vcf=pjoin(WD, "sawfish-final", "genotyped.sv.vcf.gz"),
    params:
        odir=pjoin(WD, "sawfish-final"),
    log:
        time=pjoin(WD, "times", "sawfish-jointcall.time"),
    conda:
        "../envs/sawfish.yml"
    threads: workflow.cores
    shell:
        """
        sawfish joint-call --threads 16 --sample {input.odir} --output-dir {params.odir}
        """


rule sawfish_post:
    input:
        vcf=pjoin(WD, "sawfish-final", "genotyped.sv.vcf.gz"),
    output:
        vcf=pjoin(WD, "callsets", "sawfish.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view -Oz --include "INFO/SVTYPE='DEL' | INFO/SVTYPE='INS' | INFO/SVTYPE='DUP'" {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """
