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

BENCHS = [
    "truvari-def",
    "truvari-bed",
    "truvari-nosim",
    # "truvari-sev",
]

# BENCHS = [
#     "truvari-def",
#     "truvari-easybed",
#     "truvari-hardbed",
# ]


def parse_csv(csv_fp):
    data = []
    for line in open(csv_fp):
        if line.startswith("Tool"):
            continue
        line = line.strip("\n").split(",")
        tool = line[0]
        if tool not in TOOLS:
            continue
        if "-" in tool:
            tool = tool.split("-")[0]
        if tool == "SVDSS2ht":
            tool = "SVDSS2"
        f1 = float(line[-1])
        if f1 == 0:
            continue
        data.append([tool, f1])
    return data
            
def parse_ddir(ddir, refseq=""):

    return data


def main():
    CMAP = "Purples_r"

    hprc_dir = sys.argv[1]
    giab_dir = sys.argv[2]
    refseq = sys.argv[3]

    data = []
    for csv_fp in glob.glob(f"{hprc_dir}/*.csv"):
        print(csv_fp)
        truth, bench, _ = csv_fp.split("/")[-1].split(".")
        if truth == "severus-paper":
            continue  # truth = "hapdiff"
        if bench not in BENCHS:
            continue
        if "giab" not in truth:
            truth = truth + ".hprc"
        data += [[refseq, truth, bench, tool, f1] for tool, f1 in parse_csv(csv_fp)]
    for csv_fp in glob.glob(f"{giab_dir}/*.csv"):
        print(csv_fp)
        truth, bench, _ = csv_fp.split("/")[-1].split(".")
        if truth == "severus-paper":
            continue  # truth = "hapdiff"
        if bench not in BENCHS:
            continue
        data += [[refseq, truth + ".giab", bench, tool, f1] for tool, f1 in parse_csv(csv_fp)]

    df = pd.DataFrame(data, columns=["RefSeq", "Truth", "Bench", "Tool", "F1"])
    print(df.to_string())

    truths = list(df["Truth"].unique())
    truths.sort()
    tools = list(df["Tool"].unique())
    tools.sort()

    fig, axes = plt.subplots(1, 3)
    for i, bench in enumerate(BENCHS):
        M = [[len(tools) for _ in truths] for _ in tools]
        M_annot = [["-" for _ in truths] for _ in tools]
        for col, truth in enumerate(truths):
            df2 = df[(df["Bench"] == bench) & (df["Truth"] == truth)]
            print(df2)
            f1s = [(tool, f1) for tool, f1 in zip(df2["Tool"], df2["F1"])]
            f1s.sort(key=lambda x: x[1], reverse=True)
            print(f1s)
            for rank, (tool, _) in enumerate(f1s, 1):
                M[tools.index(tool)][col] = rank
                M_annot[tools.index(tool)][col] = str(rank)
        sns.heatmap(
            M,
            ax=axes[i],
            annot=M_annot,
            fmt="",
            xticklabels=truths,
            yticklabels=tools if i == 0 else False,
            cbar=False,
            cmap=CMAP,
        )
        axes[i].set_title(bench)

    plt.suptitle(refseq)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
