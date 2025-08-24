import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    # "SVDSS2ht-q0.98-w4": "SVDSS2",
}

BENCHS = {
    "truvari-def": "Full",
    "truvari-wbed": "Conf",
}

# BENCHS = [
#     "truvari-def",
#     "truvari-easybed",
#     "truvari-hardbed",
# ]


def parse_ddir(ddir, refseq=""):
    data = []
    for csv_fp in glob.glob(f"{ddir}/*.csv"):
        truth, bench, *rest = csv_fp.split("/")[-1].split(".")
        if truth not in TRUTHS:
            continue
        if truth == "svim-asm":
            truth = "SVIM-asm"
        if bench not in BENCHS:
            continue
        bench = BENCHS[bench]
        refine = False
        if "refine" in rest:
            refine = True
        for line in open(csv_fp):
            if line.startswith("Tool"):
                continue
            line = line.strip("\n").split(",")
            tool = line[0]
            if tool not in TOOLS:
                continue
            tool = TOOLS[tool]
            p = round(float(line[-3]), 2)
            r = round(float(line[-2]), 2)
            f1 = round(float(line[-1]), 2)
            if f1 == 0:
                continue
            data.append([refseq, truth, bench, refine, tool, p, r, f1])
    return data


def main_rankmap():
    CMAP = "Purples_r"

    ddir = sys.argv[1]
    refseq = sys.argv[2]  # plot title

    data = parse_ddir(ddir, refseq)
    df = pd.DataFrame(
        data, columns=["RefSeq", "Truth", "Bench", "Refine", "Tool", "P", "R", "F1"]
    )
    print(df)

    # truths = list(df["Truth"].unique())
    # truths.sort()
    truths = ["dipcall", "SVIM-asm", "hapdiff"]
    tools = list(df["Tool"].unique())
    tools.sort()

    fig, axes = plt.subplots(2, 2, figsize=(5, 7))
    for col, bench in enumerate(["Full", "Conf"]):  # HARDCODED
        for row, refine in enumerate([False, True]):
            M = [[len(tools) for _ in truths] for _ in tools]
            M_annot = [["-" for _ in truths] for _ in tools]
            for hm_col, truth in enumerate(truths):
                df2 = df[
                    (df["Bench"] == bench)
                    & (df["Truth"] == truth)
                    & (df["Refine"] == refine)
                ]
                print(df2)
                f1s = [(tool, f1) for tool, f1 in zip(df2["Tool"], df2["F1"])]
                f1s.sort(key=lambda x: x[1], reverse=True)
                print(f1s)
                for rank, (tool, _) in enumerate(f1s, 1):
                    M[tools.index(tool)][hm_col] = rank
                    M_annot[tools.index(tool)][hm_col] = str(rank)
            g = sns.heatmap(
                M,
                ax=axes[row][col],
                annot=M_annot,
                fmt="",
                xticklabels=truths if row == 1 else False,
                yticklabels=tools if col == 0 else False,
                cbar=False,
                cmap=CMAP,
            )
            g.set_xticklabels(g.get_xticklabels(), rotation=30)
            if col == 0:
                g.set_ylabel("w/" + ("" if refine else "o") + " refine")
            else:
                g.set_ylabel("")
            if row == 0:
                title = "Full genome" if bench == "Full" else "Confident regions"
                axes[row][col].set_title(title)

    plt.suptitle(refseq)
    plt.tight_layout()
    plt.show()


def main_all():
    # sns.set(font_scale=0.7)
    t2t_ddir = sys.argv[1]
    hg38_ddir = sys.argv[2]
    hg19_ddir = sys.argv[3]

    data = []
    data += parse_ddir(t2t_ddir, "T2T-CHM13")
    data += parse_ddir(hg38_ddir, "GRCh38")
    data += parse_ddir(hg19_ddir, "GRCh37")

    df = pd.DataFrame(
        data, columns=["RefSeq", "Truth", "Bench", "Refine", "Tool", "P", "R", "F1"]
    )

    # Supplementary Table 1 (full table)
    df.sort_values(by=["RefSeq", "Truth", "Bench", "Refine", "Tool"])
    df.to_csv(sys.stdout, index=False)

    # Numbers for manuscript: avg f1
    for refseq in df["RefSeq"].unique():
        for bench in df["Bench"].unique():
            for refine in df["Refine"].unique():
                avg_f1 = df[
                    (df["RefSeq"] == refseq)
                    & (df["Bench"] == bench)
                    & (df["Refine"] == refine)
                ]["F1"].mean()
                print(refseq, bench, "ref" if refine else "noref", avg_f1, sep="\t")

    # Supplementary Table 2 (delta(refine,norefine)
    print(
        "RefSeq", "Truth", "Bench", "Tool", "F1-refine", "F1-norefine", "delta", sep=","
    )
    for refseq in df["RefSeq"].unique():
        for truth in df["Truth"].unique():
            for bench in df["Bench"].unique():
                for tool in df["Tool"].unique():
                    avg_f1_refine = df[
                        (df["RefSeq"] == refseq)
                        & (df["Truth"] == truth)
                        & (df["Bench"] == bench)
                        & (df["Tool"] == tool)
                        & (df["Refine"] == True)
                    ]["F1"].iloc[0]
                    avg_f1_norefine = df[
                        (df["RefSeq"] == refseq)
                        & (df["Truth"] == truth)
                        & (df["Bench"] == bench)
                        & (df["Tool"] == tool)
                        & (df["Refine"] == False)
                    ]["F1"].iloc[0]
                    print(
                        refseq,
                        truth,
                        bench,
                        tool,
                        avg_f1_refine,
                        avg_f1_norefine,
                        round(avg_f1_refine - avg_f1_norefine, 2),
                        sep=",",
                    )

    fig, ax1 = plt.subplots(1, 1, figsize=(6, 6), sharey=True)
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
    sns.move_legend(ax1, "upper left")

    plt.tight_layout()
    plt.show()


def main():
    import itertools

    t2t_ddir = sys.argv[1]
    hg38_ddir = sys.argv[2]
    hg19_ddir = sys.argv[3]

    data = []
    data += parse_ddir(t2t_ddir, "T2T-CHM13")
    data += parse_ddir(hg38_ddir, "GRCh38")
    data += parse_ddir(hg19_ddir, "GRCh37")

    df = pd.DataFrame(
        data, columns=["RefSeq", "Truth", "Bench", "Refine", "Tool", "P", "R", "F1"]
    )
    df.loc[df["Refine"] == True, "Refine"] = " (w/ ref)"
    df.loc[df["Refine"] == False, "Refine"] = " (w/o ref)"

    df["x"] = df["Bench"] + df["Refine"]  # .astype(str)
    df["Rank"] = df.groupby(["RefSeq", "Truth", "Bench", "Refine"])["F1"].rank(
        ascending=True
    )
    print(df)

    nrows = len(df["RefSeq"].unique())
    ncols = len(df["Truth"].unique())

    fig, axes = plt.subplots(nrows, ncols, figsize=(9, 6), sharex=True, sharey=True)

    legend = None
    for row, refseq in enumerate(sorted(df["RefSeq"].unique())):
        for col, truth in enumerate(sorted(df["Truth"].unique())):
            sub_df = df[(df["RefSeq"] == refseq) & (df["Truth"] == truth)]
            sns.lineplot(
                sub_df,
                x="x",
                y="Rank",
                hue="Tool",
                ax=axes[row][col],
                legend=(row == nrows - 1 and col == ncols - 1),
            )
            axes[row][col].set_ylim(0, 9)
            if row == 0:
                axes[row][col].set_title(truth)
            if row == 2:
                if col == 1:
                    axes[row][col].set_xlabel("Truvari setting")
                else:
                    axes[row][col].set_xlabel("")
                axes[row][col].set_xticklabels(
                    axes[row][col].get_xticklabels(), rotation=90
                )
            if col == 0:
                axes[row][col].set_ylabel(refseq)
            if row == nrows - 1 and col == ncols - 1:
                legend = axes[row][col].get_legend()
                axes[row][col].get_legend().remove()
                # sns.move_legend(
                #     axes[row][col],
                #     "upper left",
                #     bbox_to_anchor=(1, 1),
                #     title="",
                #     ncol=3,
                #     handletextpad=0.2,
                #     columnspacing=1,
                # )
    plt.legend(bbox_to_anchor=(1, 1))
    # fig.set_legend(legend)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # main()
    if sys.argv[1] == "all":
        sys.argv.pop(0)
        main_all()
    elif sys.argv[1] == "rank":
        sys.argv.pop(0)
        main_rankmap()
    # elif sys.argv[1] == "single":
    #     sys.argv.pop(0)
    #     main_bar()
