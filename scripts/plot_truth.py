import sys
import glob
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# sns.set(font_scale=2)
sns.set(style="whitegrid")

# REFSEQS = ["GRCh37", "GRCh38", "T2T"]
REFSEQS = ["hg19", "hg38", "t2t"]

def parse_vcf(vcf_fn, refseq="", name=""):
    data = []
    for record in VariantFile(vcf_fn):
        l = len(record.alts[0]) - len(record.ref)
        gt1 = record.samples[0]["GT"][0]
        gt2 = gt1
        if len(record.samples[0]["GT"]) == 2:
            gt1, gt2 = record.samples[0]["GT"]
            gt1 = gt1 if gt1 != None else 0
            gt2 = gt2 if gt2 != None else 0

        if abs(l) >= 50 and (gt1 > 0 or gt2 > 0):
            data.append(
                [
                    refseq,
                    name,
                    record.contig,
                    # record.pos,
                    l,
                    "INS" if l > 0 else "DEL",
                ]
            )
    return data


# hg19:dipcall:/Users/ld/code/svbench/data/hg19/truths/dipcall.vcf.gz hg19:hapdiff:/Users/ld/code/svbench/data/hg19/truths/hapdiff.vcf.gz hg19:svim-asm:/Users/ld/code/svbench/data/hg19/truths/svim-asm.vcf.gz hg38:dipcall:/Users/ld/code/svbench/data/hg38/truths/dipcall.vcf.gz hg38:hapdiff:/Users/ld/code/svbench/data/hg38/truths/hapdiff.vcf.gz hg38:svim-asm:/Users/ld/code/svbench/data/hg38/truths/svim-asm.vcf.gz t2t:dipcall:/Users/ld/code/svbench/data/t2t/truths/dipcall.vcf.gz t2t:hapdiff:/Users/ld/code/svbench/data/t2t/truths/hapdiff.vcf.gz t2t:svim-asm:/Users/ld/code/svbench/data/t2t/truths/svim-asm.vcf.gz hg19:dipcall.giab:/Users/ld/code/svbench/data/hg19-giab/truths/dipcall.vcf.gz hg19:hapdiff.giab:/Users/ld/code/svbench/data/hg19-giab/truths/hapdiff.vcf.gz hg19:svim-asm.giab:/Users/ld/code/svbench/data/hg19-giab/truths/svim-asm.vcf.gz hg38:dipcall.giab:/Users/ld/code/svbench/data/hg38-giab/truths/dipcall.vcf.gz hg38:hapdiff.giab:/Users/ld/code/svbench/data/hg38-giab/truths/hapdiff.vcf.gz hg38:svim-asm.giab:/Users/ld/code/svbench/data/hg38-giab/truths/svim-asm.vcf.gz t2t:dipcall.giab:/Users/ld/code/svbench/data/t2t-giab/truths/dipcall.vcf.gz t2t:hapdiff.giab:/Users/ld/code/svbench/data/t2t-giab/truths/hapdiff.vcf.gz t2t:svim-asm.giab:/Users/ld/code/svbench/data/t2t-giab/truths/svim-asm.vcf.gz hg19:giab11:/Users/ld/code/svbench/data/giab-callsets/v11/hg19.vcf.gz t2t:giab11:/Users/ld/code/svbench/data/giab-callsets/v11/t2t.vcf.gz hg38:giab11:/Users/ld/code/svbench/data/giab-callsets/v11/hg38.vcf.gz

def main():
    df = []
    for arg in sys.argv[1:]:
        ref, name, vcf = arg.split(":")
        df += parse_vcf(vcf, ref, name)

    df = pd.DataFrame(
        df, columns=["RefSeq", "Truth", "Chrom", "Len", "Type"]
    )
    print(df)
    TRUTHS = df["Truth"].unique()
    print(TRUTHS)
    
    fig, axes = plt.subplots(2, 3, figsize=(9, 11))

    for i, refseq in enumerate(REFSEQS):
        subdf = df[df["RefSeq"] == refseq]
        df2 = subdf.groupby(["Truth", "Type"]).count()
        sns.barplot(
            data=df2,
            x="Truth",
            order=TRUTHS,
            y="Chrom",
            hue="Type",
            legend=True if i == 2 else None,
            palette="Set2",
            ax=axes[0][i],
        )
        axes[0][i].set_xlabel("")  # Truth
        axes[0][i].tick_params(axis="x", labelrotation=30)
        # axes[0][i].set_ylim(0, 19000)  # Count
        axes[0][i].set_ylabel("")  # Count
        if i == 0:
            axes[0][i].set_ylabel("(a)")
        if i != 0:
            # remove y-ticks
            axes[0][i].set_yticklabels([])
        if i == 2:
            # move legends
            sns.move_legend(axes[0][i], "center left", bbox_to_anchor=(1, 0.5))

        axes[0][i].set_title(refseq)

        # variant length distribution per truthset
        sns.histplot(
            data=subdf[abs(subdf["Len"]) <= 500],
            x="Len",
            hue="Truth",
            hue_order=TRUTHS,
            element="poly",
            fill=False,
            legend=True if i == 2 else None,
            ax=axes[1][i],
        )
        axes[1][i].set_xlabel("Length")
        axes[1][i].set_xticks([-500, -250, 0, 250, 500])
        axes[1][i].set_ylabel("")  # Count
        # axes[1][i].set_ylim(0, 2000)
        if i == 0:
            axes[1][i].set_ylabel("(b)")
        if i != 0:
            # remove y-ticks
            axes[1][i].set_yticklabels([])
        if i == 2:
            sns.move_legend(axes[1][i], "center left", bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
