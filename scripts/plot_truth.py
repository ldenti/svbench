import sys
import glob
from collections import Counter

from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=2)


sns.set(style="whitegrid")

TRUTHS = ["dipcall", "svim-asm", "hapdiff"]


def main():
    vcf_dir = sys.argv[1]
    nb_dist = 50

    # we may not need this dict and do everything on df but it's more convenient to me
    truths = {}
    df = []
    for vcf_fn in glob.glob(f"{vcf_dir}/*.vcf.gz"):
        name = vcf_fn.split("/")[-1].split(".")[0]
        if name not in TRUTHS:
            continue
        truths[name] = {}
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
                if abs(l) >= 167 and abs(l) <= 171:
                    print(
                        record.contig,
                        record.pos,
                        record.pos + abs(l) if l < 0 else record.pos + 1,
                        sep="\t",
                    )
                df.append(
                    [
                        name,
                        record.contig,
                        record.pos,
                        l,
                        "INS" if l > 0 else "DEL",
                        gt1 == gt2,
                    ]
                )
                if record.contig not in truths[name]:
                    truths[name][record.contig] = []
                truths[name][record.contig].append(record.pos)

    df = pd.DataFrame(df, columns=["Truth", "Chrom", "Pos", "Len", "Type", "GT"])

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 9))

    # number of variants per truthset
    df2 = df.groupby(["Truth", "Type"]).count()
    sns.barplot(
        data=df2,
        x="Truth",
        order=TRUTHS,
        y="Chrom",
        hue="Type",
        palette="Set2",
        ax=ax1,
    )
    ax1.set_ylabel("")  # Count
    ax1.legend(loc=4)
    ax1.bar_label(ax1.containers[0])
    ax1.bar_label(ax1.containers[1])
    ax1.set_title("(a)")

    # variant length distribution per truthset
    sns.histplot(
        data=df[abs(df["Len"]) <= 500],
        x="Len",
        hue="Truth",
        element="poly",
        fill=False,
        ax=ax2,
    )
    ax2.set_xlabel("Length")
    ax2.set_ylabel("")  # Count
    ax2.set_title("(b)")

    # GT distribution
    gtdf = []
    for truth in truths:
        for b in [True, False]:
            gt = "1/1" if b else "0/1"
            gtdf.append([truth, gt, len(df[(df["Truth"] == truth) & (df["GT"] == b)])])
    gtdf = pd.DataFrame(gtdf, columns=["Truth", "GT", "Count"])

    sns.barplot(
        gtdf, x="Truth", order=TRUTHS, y="Count", hue="GT", palette="Set2", ax=ax3
    )
    ax3.set_ylabel("")  # Count
    ax3.legend(loc=4)
    ax3.set_title("(c)")

    # neighbor distribution per truthset
    df2 = []
    for truth in truths:
        neighbors = []
        for chrom in truths[truth]:
            last_p = truths[truth][chrom][0]
            neighbors.append(0)
            for p in truths[truth][chrom][1:]:
                if p - last_p > nb_dist:
                    neighbors.append(0)
                else:
                    neighbors[-1] += 1
                last_p = p
        d = {}
        for x in neighbors:
            k = str(x)
            if x >= 2:
                k = "2+"
            d[k] = d[k] + 1 if k in d else 1
        for k, v in d.items():
            df2.append([truth, k, v])

    df2 = pd.DataFrame(df2, columns=["Truth", f"#Neighbors-{nb_dist}bp", "Count"])
    sns.barplot(
        data=df2,
        x=f"#Neighbors-{nb_dist}bp",
        order=["0", "1", "2+"],
        y="Count",
        hue="Truth",
        ax=ax4,
    )
    ax4.set_ylabel("")  # Count
    ax4.set_title("(d)")

    plt.tight_layout()
    plt.show()
    # plt.savefig(vcf_dir + "/stats.png")


if __name__ == "__main__":
    main()
