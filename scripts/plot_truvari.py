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


def parse_ddir(ddir, refseq=""):
    data = []
    for csv_fp in glob.glob(f"{ddir}/*.csv"):
        truth, bench, _ = csv_fp.split("/")[-1].split(".")
        if truth == "severus-paper":
            truth = "hapdiff"
        for line in open(csv_fp):
            if line.startswith("Tool"):
                continue
            line = line.strip("\n").split(",")
            tool = line[0]
            if tool not in TOOLS:
                continue
            if bench not in BENCHS:
                continue
            if "-" in tool:
                tool = tool.split("-")[0]
            f1 = float(line[-1])
            if f1 == 0:
                continue
            data.append([refseq, truth, bench, tool, f1])
    return data


def main_rankmap():
    CMAP = "Purples_r"

    ddir = sys.argv[1]
    refseq = sys.argv[2]

    data = parse_ddir(ddir, refseq)
    df = pd.DataFrame(data, columns=["RefSeq", "Truth", "Bench", "Tool", "F1"])

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


def main_matrix():
    # sns.set(font_scale=0.7)
    t2t_ddir = sys.argv[1]
    hg38_ddir = sys.argv[2]
    hg19_ddir = sys.argv[3]

    data = []
    data += parse_ddir(t2t_ddir, "T2T")
    data += parse_ddir(hg38_ddir, "hg38")
    data += parse_ddir(hg19_ddir, "hg19")

    df = pd.DataFrame(data, columns=["RefSeq", "Truth", "Bench", "Tool", "F1"])
    print(df)

    tools = df["Tool"].unique()
    tools.sort()
    ro = df["Truth"].unique()
    ro.sort()
    co = BENCHS

    g = sns.catplot(
        data=df,
        x="Tool",
        y="F1",
        hue="RefSeq",
        col="Bench",
        row="Truth",
        kind="bar",
        order=tools,
        row_order=ro,
        col_order=co,
        sharex=True,
        sharey=True,
        height=2.5,
        aspect=1.7,
        margin_titles=True,
        legend_out=False,
    )

    g.tick_params(axis="x", labelrotation=45)  # set_xticklabels(rotation=90)
    g.set(ylim=(0, 100))

    # for row in g.axes:
    #     for ax in row:
    #         ax.bar_label(ax.containers[0], rotation=90)
    #         ax.bar_label(ax.containers[1], rotation=90)
    #         ax.bar_label(ax.containers[2], rotation=90)

    sns.move_legend(g, "upper left", bbox_to_anchor=(0.065, 0.95), title="", ncol=3)
    plt.tight_layout()
    plt.show()
    # plt.savefig(ddir + "/truvari-all.f1.png", dpi=300)
    plt.close()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5), sharey=True)

    sns.stripplot(
        df,
        x="F1",
        y="Tool",
        hue="RefSeq",
        alpha=0.6,
        s=6,
        order=tools,
        linewidth=1,
        ax=ax1,
    )
    ax1.set_xlim(35, 100)
    ax1.set_title("(a)")

    sns.stripplot(
        df[(df["RefSeq"] == "hg38") & (df["Bench"] == "truvari-nosim")],
        x="F1",
        y="Tool",
        hue="Truth",
        alpha=0.6,
        s=6,
        order=tools,
        linewidth=1,
        palette="Oranges_d",
        ax=ax2,
    )
    ax2.set_xlim(35, 100)
    ax2.set_title("(b) hg38 (truvari-nosim)")

    plt.tight_layout()
    plt.show()


def main_bar():
    sns.set(font_scale=0.75)
    csv_fn = sys.argv[1]
    df = pd.read_csv(csv_fn)
    df = df.sort_values(by=["Tool"], ascending=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
    sns.scatterplot(data=df, x="P", y="R", hue="Tool", ax=ax1, legend=None)

    for index, row in df.iterrows():
        if row["P"] >= 50 and row["R"] >= 50:
            ax1.text(row["P"] + 0.2, row["R"] + 0.2, row["Tool"])
    ax1.set_xlim(50, 100)
    ax1.set_ylim(50, 100)

    sns.barplot(data=df, x="Tool", y="F1", ax=ax2, color="steelblue", alpha=0.75)
    ax2.tick_params(axis="x", labelrotation=90)
    ax2.bar_label(ax2.containers[0])
    ax2.set_ylim(50, 100)
    plt.tight_layout()
    # plt.show()
    plt.savefig(csv_fn + ".png", dpi=300)


if __name__ == "__main__":
    if sys.argv[1] == "all1":
        sys.argv.pop(0)
        main_matrix1()
    elif sys.argv[1] == "all":
        sys.argv.pop(0)
        main_matrix()
    elif sys.argv[1] == "rank":
        sys.argv.pop(0)
        main_rankmap()
    elif sys.argv[1] == "single":
        sys.argv.pop(0)
        main_bar()
