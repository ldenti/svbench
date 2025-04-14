import sys
import glob
from collections import Counter

from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")


def main():
    vcf_dir = sys.argv[1]
    nb_dist = 50

    # we may not need this dict and do everything on df but it's more convenient to me
    truths = {}
    df = []
    for vcf_fn in glob.glob(f"{vcf_dir}/*.vcf.gz"):
        name = vcf_fn.split("/")[-1].split(".")[0]
        truths[name] = {}
        for record in VariantFile(vcf_fn):
            l = len(record.alts[0]) - len(record.ref)

            # filters=[x.name for x in record.filter.values()]
            # if name == "dipcall" and len(filters) > 0 and ("GAP1" in filters or "GAP2" in filters):
            #     continue
            # elif name != "dipcall" and "PASS" not in filters:
            #     continue

            if abs(l) >= 50:
                df.append(
                    [name, record.contig, record.pos, l, "INS" if l > 0 else "DEL"]
                )
                if record.contig not in truths[name]:
                    truths[name][record.contig] = []
                truths[name][record.contig].append(record.pos)

    df = pd.DataFrame(df, columns=["Truth", "Chrom", "Pos", "Len", "Type"])

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    # number of variants per truthset
    df2 = df.groupby(["Truth", "Type"]).count()
    sns.barplot(
        data=df2,
        x="Truth",
        y="Chrom",
        hue="Type",
        ax=ax1,
    )
    ax1.set_ylabel("Count")
    ax1.bar_label(ax1.containers[0])
    ax1.bar_label(ax1.containers[1])

    # variant length distribution per truthset
    sns.histplot(
        data=df[abs(df["Len"]) <= 500],
        x="Len",
        hue="Truth",
        element="poly",
        fill=False,
        ax=ax2,
    )

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
            if x >= 3:
                k = "3+"
            d[k] = d[k] + 1 if k in d else 1
        for k, v in d.items():
            df2.append([truth, k, v])

    df2 = pd.DataFrame(df2, columns=["Truth", f"#Neighbors-{nb_dist}bp", "Count"])
    sns.barplot(
        data=df2,
        x=f"#Neighbors-{nb_dist}bp",
        order=["0", "1", "2", "3+"],
        y="Count",
        hue="Truth",
        ax=ax3,
    )

    plt.show()


if __name__ == "__main__":
    main()
