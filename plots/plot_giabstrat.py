import sys
import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sns.set(style="whitegrid")

TRUTHS = ["dipcall", "svim-asm", "hapdiff"]

TOOLS = {
    "cutesv-w4": "cuteSV",
    "debreak": "debreak",
    "sawfish": "sawfish",
    "severus-w4": "severus",
    "sniffles": "sniffles",
    "svisionpro-w4": "SVision-pro",
    "SVDSS-w4": "SVDSS",
    "SVDSS2ht-q0.98-w4": "SVDSS2",
}

BENCHS = [
    "truvari-easybed",
    "truvari-hardbed",
]


def parse_ddir(ddir, refine, refseq=""):
    data = []
    for csv_fp in glob.glob(f"{ddir}/*.csv"):
        truth, bench, *rest = csv_fp.split("/")[-1].split(".")

        if truth not in TRUTHS:
            continue
        if truth == "svim-asm":
            truth = "SVIM-asm"

        if bench not in BENCHS:
            continue
        if "easy" in bench:
            bench = "Easy"
        else:
            bench = "Hard"

        if refine and "refine" not in rest:
            continue

        for line in open(csv_fp):
            if line.startswith("Tool"):
                continue
            line = line.strip("\n").split(",")
            tool = line[0]
            if tool not in TOOLS:
                continue
            tool = TOOLS[tool]
            f1 = float(line[-1])
            if f1 == 0:
                continue
            data.append([refseq, truth, bench, refine, tool, f1])
    return data


def main():
    sns.set(font_scale=0.8)

    parser = argparse.ArgumentParser()
    parser.add_argument("--refine", action="store_true")
    parser.add_argument("t2t")
    parser.add_argument("hg38")
    parser.add_argument("hg19")
    args = parser.parse_args()

    data = []
    data += parse_ddir(args.t2t, args.refine, "T2T")
    data += parse_ddir(args.hg38, args.refine, "GRCh38")
    data += parse_ddir(args.hg19, args.refine, "GRCh37")

    df = pd.DataFrame(
        data, columns=["RefSeq", "Truth", "Bench", "Refine", "Tool", "F1"]
    )
    print(df)

    tools = df["Tool"].unique()
    tools.sort()
    xticks = {tool: i for i, tool in enumerate(tools, 1)}

    xoff = 0.25
    x_offsets = {"dipcall": (-xoff, "r"), "SVIM-asm": (0, "g"), "hapdiff": (xoff, "b")}

    markers = {"Easy": "o", "Hard": "x"}

    refseqs = ["GRCh37", "GRCh38", "T2T"]

    fig, axes = plt.subplots(1, 3, figsize=(10, 4), sharey=True)

    for i, refseq in enumerate(refseqs):
        xticksvalues = []
        xtickslabels = []
        for tool, x in xticks.items():
            for truth, (offset, color) in x_offsets.items():
                for bench, marker in markers.items():
                    f1 = df[
                        (df["RefSeq"] == refseq)
                        & (df["Truth"] == truth)
                        & (df["Bench"] == bench)
                        & (df["Tool"] == tool)
                    ]["F1"].iloc[0]
                    print(tool, x, truth, offset, bench, f1)

                    x1 = x + offset
                    if offset == 0:
                        xticksvalues.append(x)
                        xtickslabels.append(tool)
                    axes[i].plot(x1, f1, f"{color}{marker}")
        axes[i].set_xticks(xticksvalues, xtickslabels, rotation=30)
        axes[i].set_title(refseq)
        if i == 0:
            axes[i].set_ylabel("F1")

    legend_elements = [
        Patch(facecolor="r", edgecolor="r", label="dipcall"),
        Patch(facecolor="g", edgecolor="g", label="SVIM-asm"),
        Patch(facecolor="b", edgecolor="b", label="hapdiff"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            label="Easy",
            linewidth=0,
            markerfacecolor="black",
            # markersize=15,
        ),
        Line2D(
            [0],
            [0],
            marker="x",
            color="black",
            label="Hard",
            linewidth=0,
            markerfacecolor="black",
            # markersize=15,
        ),
    ]
    axes[0].legend(handles=legend_elements)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
