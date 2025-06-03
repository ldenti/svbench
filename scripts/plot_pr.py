import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")

TOOLS = [
    "cutesv-w4",
    "debreak",
    "sawfish",
    "severus-w4",
    "sniffles",
    "svisionpro-w4",
    "SVDSS-w4",
    "SVDSS2ht-q0.98-w4",
]


def parse_csv(fn, refseq):
    data = []
    for line in open(fn):
        if line.startswith("Tool"):
            continue
        tool, TP, FP, FN, P, R, F1 = line.strip("\n").split(",")
        if tool not in TOOLS:
            continue
        if "-" in tool:
            tool = tool.split("-")[0]
        if float(F1) == 0.0:
            continue
        data.append([refseq, tool, float(P), float(R), float(F1)])
    return data


def main():
    sns.set(font_scale=0.75)
    t2t_fn = sys.argv[1]
    hg38_fn = sys.argv[2]
    hg19_fn = sys.argv[3]

    data = []
    data += parse_csv(t2t_fn, "t2t")
    data += parse_csv(hg38_fn, "hg38")
    data += parse_csv(hg19_fn, "hg19")

    df = pd.DataFrame(data, columns=["Ref", "Tool", "P", "R", "F1"])

    fig, axes = plt.subplots(2, 3)  # , figsize=(9, 5))
    for col, ref in enumerate(["t2t", "hg38", "hg19"]):
        sns.scatterplot(
            data=df[df["Ref"] == ref],
            x="P",
            y="R",
            hue="Tool",
            ax=axes[0][col],
            legend=None,
        )

        sns.barplot(
            data=df[df["Ref"] == ref],
            x="Tool",
            y="F1",
            ax=axes[1][col],
            hue="Tool",
            # alpha=0.75,
        )

        axes[0][col].set_title(ref)
        axes[0][col].set_xlim([87, 93])
        axes[0][col].set_ylim([45, 85])

        axes[1][col].set_ylim([60, 90])
        # axes[1][col].bar_label(axes[1][col].containers[0])
        axes[1][col].tick_params(axis="x", labelrotation=45)

    # for index, row in df.iterrows():
    #     if row["P"] >= 50 and row["R"] >= 50:
    #         ax1.text(row["P"] + 0.2, row["R"] + 0.2, row["Tool"])
    # ax1.set_xlim(50, 100)
    # ax1.set_ylim(50, 100)

    # sns.barplot(data=df, x="Tool", y="F1", ax=ax2, color="steelblue", alpha=0.75)
    # ax2.tick_params(axis="x", labelrotation=90)
    # ax2.bar_label(ax2.containers[0])
    # ax2.set_ylim(50, 100)
    plt.tight_layout()
    plt.show()
    # plt.savefig(csv_fn + ".png", dpi=300)


if __name__ == "__main__":
    main()
