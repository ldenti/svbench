import sys
import os
import glob
import itertools

import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")

CMAP = "bone_r"

TRUTHS = ["dipcall", "dipcallgiab", "hapdiff", "hapdiffgiab", "svimasm", "svimasmgiab", "giab11", "giab06"]

def parse_summary(fpath):
    P, R, F = 0, 0, 0
    TP, FP, FN = 0, 0, 0
    for line in open(fpath):
        line = line.strip('\n \t"')
        if line.startswith('precision"'):
            P = round(float(line.split(" ")[1][:-1]), 2)
        elif line.startswith('recall"'):
            R = round(float(line.split(" ")[1][:-1]), 2)
        elif line.startswith('f1"'):
            F = round(float(line.split(" ")[1][:-1]), 2)
        elif line.startswith('TP-base"'):
            TP = int(line.split(" ")[1][:-1])
        elif line.startswith('FP"'):
            FP = int(line.split(" ")[1][:-1])
        elif line.startswith('FN"'):
            FN = int(line.split(" ")[1][:-1])
    return P, R, F, TP, FP, FN

def build_matrix(wd):
    data = {t:{} for t in TRUTHS}
    for folder in glob.glob(os.path.join(wd, "*")):
        print(folder)
        fname = folder.split("/")[-1]
        t1, _, t2 = fname.split("-")
        fn = f"{folder}/summary.json"
        P, R, F, TP, FP, FN = parse_summary(fn)
        ACC = round(TP / (TP + FP + FN), 2)
        if t1 not in data:
            data[t1] = {}
        data[t1][t2] = ACC
    print(data)
    M = [[0 for _ in TRUTHS] for _ in TRUTHS]
    M_annot = [[0 for _ in TRUTHS] for _ in TRUTHS]
    for i1,t1 in enumerate(TRUTHS):
        for i2,t2 in enumerate(TRUTHS):
            M[i2][i1] = data[t1][t2] if t2 in data[t1] else 0
            M_annot[i2][i1] = str(M[i2][i1]) if M[i2][i1] != 0 else ""
    return M, M_annot

def main():
    t2t_wd=sys.argv[1]
    hg38_wd=sys.argv[2]
    hg19_wd=sys.argv[3]

    fig, axes = plt.subplots(1,3, figsize=(11,5))

    titles = ["GRCh37", "GRCh38", "T2T"]

    for i, wd in enumerate([hg19_wd, hg38_wd, t2t_wd]):
        M, M_annot = build_matrix(wd)

        sns.heatmap(
            M,
            square=True,
            annot=M_annot,
            fmt="",
            xticklabels=TRUTHS,
            cmap=CMAP,
            yticklabels=TRUTHS if i == 0 else False,
            cbar=False,
            vmin=0.3,
            vmax=1.0,
            ax=axes[i],
        )
        axes[i].set_title(titles[i])
    #         
    #         for ax, title in zip(axes, x_titles):
    #             ax.set_title(title)
    #             ax.tick_params(axis="x", labelrotation=0)
    #             ax.tick_params(axis="y", labelrotation=90)
    axes[1].set_xlabel("How much of this...")
    axes[0].set_ylabel("...matches with this")
    plt.tight_layout()
    # plt.savefig("x.png")
    plt.show()


if __name__ == "__main__":
    main()
