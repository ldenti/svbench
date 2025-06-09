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
    # "svisionpro-w4",
    "SVDSS-w4",
    "SVDSS2ht-q0.98-w4",
]

TOOLS2 = [
    "svdss-severus",
    "cutesv-severus",
    "debreak-severus",
    "severus-severus",
    "sniffles-severus",
    # "svisionpro-severus",
]

x_order = [
    "SVDSS",
    "SVDSS2ht",
    "cuteSV",
    "debreak",
    "sawfish",
    "severus",
    "sniffles",
    # "svisionpro",
]


def parse_csv(fn, refseq, tools=TOOLS):
    data = []
    for line in open(fn):
        if line.startswith("Tool"):
            continue
        tool, TP, FP, FN, P, R, F1 = line.strip("\n").split(",")
        if tool not in tools:
            continue
        if "-" in tool:
            tool = tool.split("-")[0]
        if tool == "cutesv":
            tool = "cuteSV"
        elif tool == "svdss":
            tool = "SVDSS"
        if float(F1) == 0.0:
            continue
        data.append([refseq, tool, float(P), float(R), float(F1)])
    return data


def main():
    sns.set(font_scale=0.75)
    t2t_fn = sys.argv[1]
    hg38_fn = sys.argv[2]
    hg19_fn = sys.argv[3]
    hg19_v06_fn = sys.argv[4]

    fig, axes = plt.subplots(2, 4, figsize=(10, 5))

    # GIAB v0.6
    col = 0
    ref = "hg19"
    data = parse_csv(hg19_v06_fn, "hg19", TOOLS)
    df = pd.DataFrame(data, columns=["Ref", "Tool", "P", "R", "F1"])
    sns.scatterplot(
        data=df[(df["Ref"] == ref) & (df["F1"] > 0)],
        x="P",
        y="R",
        hue="Tool",
        hue_order=x_order,
        ax=axes[0][col],
        legend=None,
    )
    sns.barplot(
        data=df[df["Ref"] == ref],
        x="Tool",
        y="F1",
        order=x_order,
        hue_order=x_order,
        ax=axes[1][col],
        hue="Tool",
    )
    axes[0][col].set_title("(a)\n" + ref + " (GIAB v0.6, Tier 1)")
    axes[0][col].set_xlim([75, 100])
    axes[0][col].set_ylim([45, 100])
    axes[1][col].set_ylim([0, 100])
    # axes[1][col].bar_label(axes[1][col].containers[0])
    axes[1][col].tick_params(axis="x", labelrotation=90)

    data = []
    data += parse_csv(t2t_fn, "t2t")
    data += parse_csv(hg38_fn, "hg38")
    data += parse_csv(hg19_fn, "hg19")

    df = pd.DataFrame(data, columns=["Ref", "Tool", "P", "R", "F1"])
    print(df)

    for col, ref in enumerate(["hg19", "hg38", "t2t"], 1):
        sns.scatterplot(
            data=df[df["Ref"] == ref],
            x="P",
            y="R",
            hue="Tool",
            hue_order=x_order,
            ax=axes[0][col],
            legend=None,
        )

        sns.barplot(
            data=df[df["Ref"] == ref],
            x="Tool",
            y="F1",
            order=x_order,
            hue_order=x_order,
            ax=axes[1][col],
            hue="Tool",
            # alpha=0.75,
        )

        axes[0][col].set_title(f"({('abcd'[col])})\n" + ref + " (GIAB v1.1, w/ BED)")
        axes[0][col].set_xlim([75, 100])
        axes[0][col].set_ylim([45, 100])

        axes[1][col].set_ylim([0, 100])
        # axes[1][col].bar_label(axes[1][col].containers[0])
        axes[1][col].tick_params(axis="x", labelrotation=90)

        if col != 0:
            axes[0][col].set_ylabel("")
            axes[1][col].set_ylabel("")

    plt.tight_layout()
    plt.show()
    # plt.savefig(csv_fn + ".png", dpi=300)


if __name__ == "__main__":
    main()
