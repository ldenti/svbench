import sys
import os
import argparse
import glob
import pandas as pd
import numpy as np
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
    # "SVDSS2ht-q0.98-w4",
]

x_order = [
    "SVDSS",
    # "SVDSS2",
    "cuteSV",
    "debreak",
    "sawfish",
    "severus",
    "sniffles",
    "SVision-pro",
]


def parse_csv(fn, refseq, giabv, mode, refine):
    data = []
    for line in open(fn):
        if line.startswith("Tool"):
            continue
        tool, TP, FP, FN, P, R, F1 = line.strip("\n").split(",")
        if tool not in TOOLS:
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
        data.append([refseq, giabv, mode, refine, tool, float(P), float(R), float(F1)])
    data.sort(key=lambda x: x[-1], reverse=True)
    for i, row in enumerate(data, 1):
        row.append(i)
    return data


def main():
    sns.set(font_scale=0.75)

    parser = argparse.ArgumentParser()
    parser.add_argument("t2t")
    parser.add_argument("hg38")
    parser.add_argument("hg19")
    args = parser.parse_args()

    data = []
    for mode in ["def", "wbed"]:
        for refine in [".refine", ""]:
            m = "Full" if mode == "def" else "Conf"
            r = refine != ""
            data += parse_csv(
                args.hg19 + f"/giab-v0.6-{mode}{refine}.csv", "GRCh37", "v0.6", m, r
            )
            data += parse_csv(
                args.t2t + f"/giab-v1.1-{mode}{refine}.csv", "T2T", "v1.1", m, r
            )
            data += parse_csv(
                args.hg38 + f"/giab-v1.1-{mode}{refine}.csv", "GRCh38", "v1.1", m, r
            )
            data += parse_csv(
                args.hg19 + f"/giab-v1.1-{mode}{refine}.csv", "GRCh37", "v1.1", m, r
            )

    df = pd.DataFrame(
        data,
        columns=["RefSeq", "Truth", "Bench", "Refine", "Tool", "P", "R", "F1", "rank"],
    )

    # Supplementary Table 5 (full table)
    df.sort_values(["RefSeq", "Truth", "Bench", "Refine"]).to_csv(
        sys.stdout, index=False
    )

    # Supplementary Table 6 (delta(refine,norefine)
    print(
        "RefSeq", "Truth", "Bench", "Tool", "F1-refine", "F1-norefine", "delta", sep=","
    )
    for refseq in df["RefSeq"].unique():
        for truth in df["Truth"].unique():
            if truth == "v0.6" and refseq != "GRCh37":
                continue
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

    for bench in ["Full", "Conf"]:
        for refine in [True, False]:
            print(f"### {bench} /", "Refine" if refine else "NoRefine", "###")
            old_f1 = df[
                (df["Bench"] == bench)
                & (df["Refine"] == refine)
                & (df["Truth"] == "v0.6")
            ].sort_values(["Tool"])[["Tool", "F1"]]
            new_f1 = df[
                (df["Bench"] == bench)
                & (df["Refine"] == refine)
                & (df["Truth"] == "v1.1")
                & (df["RefSeq"] == "GRCh37")
            ].sort_values(["Tool"])[["Tool", "F1"]]
            print(old_f1)
            print(new_f1)
            corr = np.corrcoef(old_f1["F1"], new_f1["F1"])[0, 1]
            print(corr)

            fig, axes = plt.subplots(2, 4, figsize=(10, 5))

            # GIAB v0.6
            sub_df = df[
                (df["Bench"] == bench)
                & (df["Refine"] == refine)
                & (df["Truth"] == "v0.6")
            ]
            ref = sub_df["RefSeq"].unique()[0]
            col = 0

            sns.scatterplot(
                data=sub_df,
                x="P",
                y="R",
                hue="Tool",
                hue_order=x_order,
                ax=axes[0][col],
                legend=None,
            )
            sns.barplot(
                data=sub_df,
                x="Tool",
                y="F1",
                order=x_order,
                hue_order=x_order,
                ax=axes[1][col],
                hue="Tool",
            )

            axes[0][col].set_title(
                f"(a)\n"
                + ref
                + " (GIAB v0.6"
                + (")" if bench == "Full" else ", Tier 1)")
            )

            axes[1][col].set_ylim([0, 100])
            for i, container in enumerate(axes[1][col].containers):
                axes[1][col].bar_label(
                    container, labels=sub_df[sub_df["Tool"] == x_order[i]]["rank"]
                )
            axes[1][col].tick_params(axis="x", labelrotation=90)

            # GIAB v1.1
            for col, ref in enumerate(["GRCh37", "GRCh38", "T2T"], 1):
                sub_df = df[
                    (df["Bench"] == bench)
                    & (df["Refine"] == refine)
                    & (df["RefSeq"] == ref)
                    & (df["Truth"] == "v1.1")
                ]

                sns.scatterplot(
                    data=sub_df,
                    x="P",
                    y="R",
                    hue="Tool",
                    hue_order=x_order,
                    ax=axes[0][col],
                    legend=None,
                )

                sns.barplot(
                    data=sub_df,
                    x="Tool",
                    y="F1",
                    order=x_order,
                    hue_order=x_order,
                    ax=axes[1][col],
                    hue="Tool",
                    # alpha=0.75,
                )

                axes[0][col].set_title(
                    f"({('abcd'[col])})\n"
                    + ref
                    + " (GIAB v1.1"
                    + (")" if bench == "Full" else ", w/ BED)")
                )
                axes[1][col].set_ylim([0, 100])

                for i, container in enumerate(axes[1][col].containers):
                    axes[1][col].bar_label(
                        container,
                        labels=sub_df[(sub_df["Tool"] == x_order[i])]["rank"],
                    )

                axes[1][col].tick_params(axis="x", labelrotation=90)

                if col != 0:
                    axes[0][col].set_ylabel("")
                    axes[1][col].set_ylabel("")

            xlim = 40  # 80 if args.bed else 40
            ylim = 30  # 40 if args.bed else 30
            for col in range(4):
                axes[0][col].set_xlim([xlim, 100])
                axes[0][col].set_ylim([ylim, 100])
            plt.tight_layout()
            plt.show()
            plt.close()


if __name__ == "__main__":
    main()
