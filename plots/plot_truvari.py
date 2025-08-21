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
    "SVDSS2ht-q0.98-w4": "SVDSS2",
}

BENCHS = [
    "truvari-def",
    "truvari-wbed",
]

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
            if bench not in BENCHS:
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
    # force svision-pro to be the last one
    # tools = tools[:2] + tools[3:] + [tools[2]]

    fig, axes = plt.subplots(2, 2, figsize=(5, 7))
    for col, bench in enumerate(BENCHS):
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
            axes[row][col].set_title(bench + (" (refine)" if refine else ""))

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
    data += parse_ddir(hg38_ddir, "GRCh38")
    data += parse_ddir(hg19_ddir, "GRCh37")

    df = pd.DataFrame(
        data, columns=["RefSeq", "Truth", "Bench", "Refine", "Tool", "P", "R", "F1"]
    )
    print(df)
    df.sort_values(by=["RefSeq", "Truth", "Bench", "Refine", "Tool"])
    df.to_csv(sys.stdout, index=False)

    # df_latex = []
    # for refseq in ["T2T", "GRCh38", "GRCh37"]:
    #     for truth in ["dipcall", "SVIM-asm", "hapdiff"]:
    #         for tool in sorted(df["Tool"].unique()):
    #             row = [refseq, truth, tool]
    #             for bench in BENCHS:
    #                 for refine in [False, True]:
    #                     df_row = df[(df["RefSeq"] == refseq) & (df["Truth"] == truth) & (df["Bench"] == bench) & (df["Refine"] == refine) & (df["Tool"] == tool)]
    #                     row.append(float(df_row["P"].iloc[0]))
    #                     row.append(float(df_row["R"].iloc[0]))
    #                     row.append(float(df_row["F1"].iloc[0]))
    #             df_latex.append(row)
    # # first three: full genome (no refine), full genome (refine), confident regions (no refine), confident regions (refine)
    # df_latex = pd.DataFrame(df_latex, columns = ["RefSeq", "Truth", "Caller", "P", "R", "F1", "P", "R", "F1", "P", "R", "F1", "P", "R", "F1"])
    # # df_latex.to_latex(sys.stdout, index=False)
    # df_latex.to_csv(sys.stdout, index=False)

    # avg f1
    for refseq in ["T2T", "GRCh38", "GRCh37"]:
        for bench in BENCHS:
            for refine in [False, True]:
                avg_f1 = df[
                    (df["RefSeq"] == refseq)
                    & (df["Bench"] == bench)
                    & (df["Refine"] == refine)
                ]["F1"].mean()
                print(refseq, bench, "ref" if refine else "noref", avg_f1, sep="\t")

    tools = df["Tool"].unique()
    tools.sort()
    # ro = df["Truth"].unique()
    # ro.sort()
    ro = ["dipcall", "SVIM-asm", "hapdiff"]
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
        height=2.3,
        aspect=1.6,
        margin_titles=True,
        legend_out=False,
    )

    g.tick_params(axis="x", labelrotation=60)  # set_xticklabels(rotation=90)
    g.set(ylim=(0, 100))

    # for row in g.axes:
    #     for ax in row:
    #         ax.bar_label(ax.containers[0], rotation=90)
    #         ax.bar_label(ax.containers[1], rotation=90)
    #         ax.bar_label(ax.containers[2], rotation=90)

    sns.move_legend(
        g,
        "upper left",
        bbox_to_anchor=(0.07, 0.95),
        title="",
        ncol=3,
        handletextpad=0.2,
        columnspacing=1,
    )
    plt.tight_layout()
    plt.show()
    # plt.savefig(ddir + "/truvari-all.f1.png", dpi=300)
    plt.close()

    # original: all points on left, hg38 only on right plot
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
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
    # ax1.set_title("(a)")
    # sns.stripplot(
    #     df[(df["RefSeq"] == "hg38") & (df["Bench"] == "truvari-nosim")],
    #     x="F1",
    #     y="Tool",
    #     hue="Truth",
    #     alpha=0.6,
    #     s=6,
    #     order=tools,
    #     linewidth=1,
    #     palette="Oranges_d",
    #     ax=ax2,
    # )
    # ax2.set_xlim(35, 100)
    # ax2.set_title("(b) hg38 (truvari-nosim)")

    plt.tight_layout()
    plt.show()


# def main_bar():
#     sns.set(font_scale=0.75)
#     csv_fn = sys.argv[1]
#     df = pd.read_csv(csv_fn)
#     df = df.sort_values(by=["Tool"], ascending=True)

#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
#     sns.scatterplot(data=df, x="P", y="R", hue="Tool", ax=ax1, legend=None)

#     for index, row in df.iterrows():
#         if row["P"] >= 50 and row["R"] >= 50:
#             ax1.text(row["P"] + 0.2, row["R"] + 0.2, row["Tool"])
#     ax1.set_xlim(50, 100)
#     ax1.set_ylim(50, 100)

#     sns.barplot(data=df, x="Tool", y="F1", ax=ax2, color="steelblue", alpha=0.75)
#     ax2.tick_params(axis="x", labelrotation=90)
#     ax2.bar_label(ax2.containers[0])
#     ax2.set_ylim(50, 100)
#     plt.tight_layout()
#     # plt.show()
#     plt.savefig(csv_fn + ".png", dpi=300)


if __name__ == "__main__":
    if sys.argv[1] == "all":
        sys.argv.pop(0)
        main_matrix()
    elif sys.argv[1] == "rank":
        sys.argv.pop(0)
        main_rankmap()
    # elif sys.argv[1] == "single":
    #     sys.argv.pop(0)
    #     main_bar()
