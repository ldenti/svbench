import sys
import glob
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# sns.set(font_scale=2)
sns.set(style="whitegrid")

TRUTHS = ["dipcall", "svim-asm", "hapdiff"]


def parse_dir(ddir, refseq=""):
    data = []
    for vcf_fn in glob.glob(f"{ddir}/*.vcf.gz"):
        name = vcf_fn.split("/")[-1].split(".")[0]
        if name not in TRUTHS:
            continue
        for record in VariantFile(vcf_fn):
            l = len(record.alts[0]) - len(record.ref)
            # filters=[x.name for x in record.filter.values()]
            # if name == "dipcall" and len(filters) > 0 and ("GAP1" in filters or "GAP2" in filters):
            #     continue
            # elif name != "dipcall" and "PASS" not in filters:
            #     continue

            gt1, gt2 = record.samples[0]["GT"]
            gt1 = gt1 if gt1 != None else 0
            gt2 = gt2 if gt2 != None else 0

            if abs(l) >= 50 and (gt1 > 0 or gt2 > 0):
                data.append(
                    [
                        refseq,
                        name,
                        record.contig,
                        record.pos,
                        l,
                        "INS" if l > 0 else "DEL",
                        gt1 == gt2,
                    ]
                )
    return data


def main():
    t2t_ddir = sys.argv[1]
    hg38_ddir = sys.argv[2]
    hg19_ddir = sys.argv[3]

    df = []
    df += parse_dir(t2t_ddir, "t2t")
    df += parse_dir(hg38_ddir, "hg38")
    df += parse_dir(hg19_ddir, "hg19")

    REFSEQS = ["t2t", "hg38", "hg19"]

    df = pd.DataFrame(
        df, columns=["RefSeq", "Truth", "Chrom", "Pos", "Len", "Type", "GT"]
    )

    fig, axes = plt.subplots(2, 3, figsize=(10, 5))

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
        axes[0][i].set_ylim(0, 19000)  # Count
        axes[0][i].set_ylabel("")  # Count
        if i != 0:
            # remove y-ticks
            axes[0][i].set_yticklabels([])
        if i == 2:
            # move legends
            sns.move_legend(axes[0][i], "center left", bbox_to_anchor=(1, 0.5))

        # ax1.bar_label(ax1.containers[0])
        # ax1.bar_label(ax1.containers[1])
        axes[0][i].set_title(refseq)

        # variant length distribution per truthset
        sns.histplot(
            data=subdf[abs(subdf["Len"]) <= 500],
            x="Len",
            hue="Truth",
            element="poly",
            fill=False,
            legend=True if i == 2 else None,
            ax=axes[1][i],
        )
        axes[1][i].set_xlabel("Length")
        axes[1][i].set_ylabel("")  # Count
        axes[1][i].set_ylim(0, 2000)
        if i != 0:
            # remove y-ticks
            axes[1][i].set_yticklabels([])
        if i == 2:
            sns.move_legend(axes[1][i], "center left", bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
