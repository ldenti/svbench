import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")


def main_matrix():
    sns.set(font_scale=0.9)
    ddir = sys.argv[1]

    data = []
    for csv_fp in glob.glob(f"{ddir}/*.csv"):
        truth, bench, _ = csv_fp.split("/")[-1].split(".")
        for line in open(csv_fp):
            if line.startswith("Tool"):
                continue
            line = line.strip("\n").split(",")
            tool = line[0]
            f1 = float(line[-1])
            data.append([truth, bench, tool, f1])
    df = pd.DataFrame(data, columns=["Truth", "Bench", "Tool", "F1"])
    print(df)

    tools = df["Tool"].unique()
    tools.sort()

    if len(df["Truth"].unique()) == 2:
        new_data = []
        for row in data:
            if row[0] == "dipcall":
                new_data.append(["Severus", row[1], row[2], 0.0])
        new_df = pd.DataFrame(new_data, columns=["Truth", "Bench", "Tool", "F1"])
        df = pd.concat([df, new_df], ignore_index=True)

    g = sns.FacetGrid(df, row="Truth", col="Bench", margin_titles=True)
    g.map(
        sns.barplot,
        "Tool",
        "F1",
        color="forestgreen",
        alpha=0.75,
        order=tools,
    )
    g.tick_params(axis="x", labelrotation=90)  # set_xticklabels(rotation=90)
    g.set(ylim=(40, 100))

    for row in g.axes:
        for ax in row:
            ax.bar_label(ax.containers[0])

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
        ax1.text(row["P"] + 0.2, row["R"] + 0.2, row["Tool"])
    # ax1.set_xlim(50, 100)
    # ax1.set_ylim(50, 100)

    sns.barplot(data=df, x="Tool", y="F1", ax=ax2, color="steelblue", alpha=0.75)
    ax2.tick_params(axis="x", labelrotation=90)
    ax2.bar_label(ax2.containers[0])

    plt.tight_layout()
    # plt.show()
    plt.savefig(csv_fn + ".png", dpi=300)


if __name__ == "__main__":
    if sys.argv[1] == "all":
        sys.argv.pop(0)
        main_matrix()
    elif sys.argv[1] == "single":
        sys.argv.pop(0)
        main_bar()
