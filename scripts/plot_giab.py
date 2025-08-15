import sys
import os
import argparse
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
    "SVDSS2",
    "cuteSV",
    "debreak",
    "sawfish",
    "severus",
    "sniffles",
    "SVision-pro",
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
        elif tool == "SVDSS2ht":
            tool = "SVDSS2"
        elif tool == "svisionpro":
            tool = "SVision-pro"
        if float(F1) == 0.0:
            continue
        data.append([refseq, tool, float(P), float(R), float(F1)])
    data.sort(key=lambda x: x[-1], reverse=True)
    for i, row in enumerate(data, 1):
        row.append(i)
    return data


def main():
    sns.set(font_scale=0.75)

    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", action='store_true')
    parser.add_argument("t2t")
    parser.add_argument("hg38")
    parser.add_argument("hg19")
    parser.add_argument("hg19old")
    args = parser.parse_args()

    fig, axes = plt.subplots(2, 4, figsize=(10, 5))

    # GIAB v0.6
    col = 0
    ref = "GRCh37"
    data = parse_csv(args.hg19old, "GRCh37", TOOLS)
    df = pd.DataFrame(data, columns=["Ref", "Tool", "P", "R", "F1", "rank"])
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
    axes[0][col].set_title("(a)\n" + ref + " (GIAB v0.6" + (")" if not args.bed else ", Tier 1)"))
    axes[1][col].set_ylim([0, 100])
    for i, container in enumerate(axes[1][col].containers):
        axes[1][col].bar_label(container, labels=df[df["Tool"] == x_order[i]]["rank"])
    axes[1][col].tick_params(axis="x", labelrotation=90)

    data = []
    data += parse_csv(args.t2t, "T2T")
    data += parse_csv(args.hg38, "GRCh38")
    data += parse_csv(args.hg19, "GRCh37")

    df = pd.DataFrame(data, columns=["Ref", "Tool", "P", "R", "F1", "rank"])
    print(df)

    for col, ref in enumerate(["GRCh37", "GRCh38", "T2T"], 1):
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

        axes[0][col].set_title(f"({('abcd'[col])})\n" + ref + " (GIAB v1.1" + (")" if not args.bed else ", w/ BED)"))
        axes[1][col].set_ylim([0, 100])

        for i, container in enumerate(axes[1][col].containers):
            axes[1][col].bar_label(
                container,
                labels=df[(df["Ref"] == ref) & (df["Tool"] == x_order[i])]["rank"],
            )

        axes[1][col].tick_params(axis="x", labelrotation=90)

        if col != 0:
            axes[0][col].set_ylabel("")
            axes[1][col].set_ylabel("")

    xlim = 40 # 80 if args.bed else 40
    ylim = 30 # 40 if args.bed else 30
    for col in range(4):
        axes[0][col].set_xlim([xlim, 100])
        axes[0][col].set_ylim([ylim, 100])
    plt.tight_layout()
    plt.show()
    # plt.savefig(csv_fn + ".png", dpi=300)


if __name__ == "__main__":
    main()
