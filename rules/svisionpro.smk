SVPRO_REPO = pjoin(WD, "software", "SVision-pro")


rule get_minda:
    output:
        exe=pjoin(WD, "software", "minda", "minda.py"),
    shell:
        """
        git clone https://github.com/songbowang125/SVision-pro.git

        git clone https://github.com/KolmogorovLab/minda.git {MINDA_REPO}
        """


rule minda:
    input:
        exe=pjoin(WD, "software", "minda", "minda.py"),
        vcf=pjoin(WD, "callsets", "{caller}.vcf.gz"),
        truth=pjoin(WD, "truths", "{truth}.vcf.gz"),
    output:
        directory(pjoin(WD, "minda-{truth}", "{caller}")),
    conda:
        pjoin(MINDA_REPO, "environment.yml")
    shell:
        """
        {input.exe} truthset --out_dir {output} --base {input.truth} --vcfs {input.vcf}
        """
