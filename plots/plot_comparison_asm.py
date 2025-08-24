import sys
import os
import argparse
import glob
import itertools

import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")

CMAP = "bone_r"


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


def main():
    print("This script assumes <= 3 truth callsets!", file=sys.stderr)

    parser = argparse.ArgumentParser()
    parser.add_argument("--refine", action="store_true")
    parser.add_argument("t2t")
    parser.add_argument("hg38")
    parser.add_argument("hg19")
    args = parser.parse_args()

    labels = ["dipcall", "SVIM-asm", "hapdiff"]
    indexes = {"dipcall": 0, "svim-asm": 1, "hapdiff": 2}

    INDIRS = [args.t2t, args.hg38, args.hg19]

    for bench in ["def", "wbed"]:
        fig, axes = plt.subplots(1, 3, figsize=(9, 4))
        for i, ddir in enumerate(INDIRS):
            # assuming this order: t2t, hg38, and hg19, but we want reverse order
            i = 2 - i

            M = [[0 for _ in labels] for _ in labels]
            M_annot = [[0 for _ in labels] for _ in labels]
            for t1, t2 in itertools.product(labels, labels):
                t1 = t1.lower()
                t2 = t2.lower()
                ACC = 0
                if t1 == t2:
                    ACC = 1.0
                else:
                    n = "ga4gh_with_refine.summary" if args.refine else "summary"
                    fn = f"{ddir}/truths/comparison-{bench}/{t1}-against-{t2}/{n}.json"
                    if os.path.exists(fn):
                        P, R, F, TP, FP, FN = parse_summary(fn)
                        ACC = round(TP / (TP + FP + FN), 2)
                        print(fn, ACC)
                    else:
                        ACC = 0.0
                M[indexes[t2]][indexes[t1]] = ACC
                M_annot[indexes[t2]][indexes[t1]] = str(ACC)
            # Force diagonal
            M[1][0] = 0
            M_annot[1][0] = ""
            M[2][0] = 0
            M_annot[2][0] = ""
            M[2][1] = 0
            M_annot[2][1] = ""

            sns.heatmap(
                M,
                square=True,
                annot=M_annot,
                fmt="",
                xticklabels=labels,
                cmap=CMAP,
                yticklabels=labels if i == 0 else False,
                cbar=False,
                vmin=0.3,
                vmax=1.0,
                ax=axes[i],
            )

        x_titles = ["GRCh37", "GRCh38", "T2T-CHM13"]
        for ax, title in zip(axes, x_titles):
            ax.set_title(title)
            ax.tick_params(axis="x", labelrotation=0)
            ax.tick_params(axis="y", labelrotation=90)
        axes[1].set_xlabel("How much of this...")
        axes[0].set_ylabel("...matches with this")
        plt.tight_layout()
        plt.show()
        plt.close()


if __name__ == "__main__":
    main()
