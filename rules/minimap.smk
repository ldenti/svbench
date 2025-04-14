rule minimap:
    input:
        fa=REF,
        fq=READS,
    output:
        bam=pjoin(WD, "minimap2", "alignments.bam"),
    threads: workflow.cores
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -ax map-hifi --MD --eqx -Y -R '@RG\\tID:{SAMPLE_NAME}\\tSM:{SAMPLE_NAME}' -t {threads} {input.fa} {input.fq} | samtools view -bS | samtools sort -@ $(({threads}-1)) -T {output.bam} > {output.bam}
        samtools index {output.bam}
        """
