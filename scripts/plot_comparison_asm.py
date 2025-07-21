import sys
import os
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

    labels = ["dipcall", "SVIM-asm", "hapdiff"]
    indexes = {"dipcall": 0, "svim-asm": 1, "hapdiff": 2}

    for bench in ["def", "bed"]:
        fig, axes = plt.subplots(1, 3, figsize=(9, 4))
        for i, ddir in enumerate(sys.argv[1:]):
            # assuming this order: t2t, hg38, and hg19, but we want reverse order
            i = 2 - i

            M = [[0 for _ in labels] for _ in labels]
            M_annot = [[0 for _ in labels] for _ in labels]
            for t1, t2 in itertools.product(labels, labels):
                t1 = t1.lower()
                t2 = t2.lower()
                fn = f"{ddir}/truths/comparison-{bench}/{t1}-against-{t2}/summary.json"
                P, R, F, TP, FP, FN = parse_summary(fn)

                ACC = round(TP / (TP + FP + FN), 2)

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

        x_titles = ["GRCh37", "GRCh38", "T2T"]
        for ax, title in zip(axes, x_titles):
            ax.set_title(title)
            ax.tick_params(axis="x", labelrotation=0)
            ax.tick_params(axis="y", labelrotation=90)
        axes[1].set_xlabel("How much of this...")
        axes[0].set_ylabel("...matches with this")
        plt.tight_layout()
        plt.show()
        plt.close()


def main_all():
    print("This script assumes <= 3 truth callsets!", file=sys.stderr)

    labels = ["dipcall", "svim-asm", "hapdiff"]
    indexes = {"dipcall": 0, "svim-asm": 1, "hapdiff": 2}
    pairs = [
        ("dipcall", "svim-asm"),
        ("dipcall", "hapdiff"),
        ("svim-asm", "hapdiff"),
    ]

    fig, axes = plt.subplots(3, 3)
    # axes = axes.reshape(-1)
    for i, ddir in enumerate(sys.argv[1:]):
        # assuming this order: t2t, hg38, and hg19
        M1 = [[0 for _ in labels] for _ in labels]
        M2 = [[0 for _ in labels] for _ in labels]
        M3 = [[0 for _ in labels] for _ in labels]
        M3_labels = [["" for _ in labels] for _ in labels]
        for t1, t2 in itertools.product(labels, labels):
            # for t1, t2 in pairs:
            fn = f"{ddir}/truths/comparison-def/{t1}-against-{t2}/summary.json"
            P1, R1, F1, TP1, FP1, FN1 = parse_summary(fn)
            fn = f"{ddir}/truths/comparison-bed/{t1}-against-{t2}/summary.json"
            P2, R2, F2, TP2, FP2, FN2 = parse_summary(fn)

            M1[indexes[t2]][indexes[t1]] = P1
            # M1[indexes[t1]][indexes[t2]] = R1
            # M1[indexes[t1]][indexes[t1]] = 1
            # M1[indexes[t2]][indexes[t2]] = 1

            M2[indexes[t2]][indexes[t1]] = P2
            # M2[indexes[t1]][indexes[t2]] = R2
            # M2[indexes[t1]][indexes[t1]] = 1
            # M2[indexes[t2]][indexes[t2]] = 1

            M3[indexes[t2]][indexes[t1]] = (P1 + P2) / 2
            # M3[indexes[t1]][indexes[t2]] = (R1 + R2) / 2
            # M3[indexes[t1]][indexes[t1]] = 1
            # M3[indexes[t2]][indexes[t2]] = 1

            M3_labels[indexes[t2]][indexes[t1]] = f"{P1}:{P2}"
            # M3_labels[indexes[t1]][indexes[t2]] = f"{R1}:{R2}"
            # M3_labels[indexes[t1]][indexes[t1]] = "1.0:1.0"
            # M3_labels[indexes[t2]][indexes[t2]] = "1.0:1.0"

        sns.heatmap(
            M1,
            square=True,
            annot=True,
            xticklabels=labels if i == 2 else False,
            yticklabels=labels,
            cbar=False,
            cmap=CMAP,
            vmin=0.5,
            vmax=1.0,
            ax=axes[i][0],
        )

        sns.heatmap(
            M2,
            square=True,
            annot=True,
            xticklabels=labels if i == 2 else False,
            yticklabels=False,
            cbar=False,
            cmap=CMAP,
            vmin=0.5,
            vmax=1.0,
            ax=axes[i][1],
        )

        sns.heatmap(
            M3,
            square=True,
            annot=True,
            xticklabels=labels if i == 2 else False,
            yticklabels=False,
            cbar=False,
            cmap=CMAP,
            vmin=0.5,
            vmax=1.0,
            ax=axes[i][2],
        )

    x_titles = ["truvari-def", "truvari-bed", "avg(def, bed)"]
    for ax, title in zip(axes[0], x_titles):
        ax.set_title(title)
    y_titles = ["T2T", "hg38", "hg19"]
    for ax, title in zip(axes[:, 2], y_titles):
        ax.set_ylabel(title, rotation=270, labelpad=15)
        ax.yaxis.set_label_position("right")

    axes[1][0].set_ylabel("...match with this")
    axes[2][1].set_xlabel("How much of this...")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
