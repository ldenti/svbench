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


def main_avg(acc=False):
    print("This script assumes <= 3 truth callsets!", file=sys.stderr)

    labels = ["dipcall", "svim-asm", "hapdiff"]
    indexes = {"dipcall": 0, "svim-asm": 1, "hapdiff": 2}
    pairs = [
        ("dipcall", "svim-asm"),
        ("dipcall", "hapdiff"),
        ("svim-asm", "hapdiff"),
    ]

    fig, axes = plt.subplots(1, 3)
    for i, ddir in enumerate(sys.argv[1:]):
        # assuming this order: t2t, hg38, and hg19
        M3 = [[0 for _ in labels] for _ in labels]
        M3_labels = [["" for _ in labels] for _ in labels]
        for t1, t2 in itertools.product(labels, labels):
            fn = f"{ddir}/comparison-def/{t1}-against-{t2}/summary.json"
            P1, R1, F1, TP1, FP1, FN1 = parse_summary(fn)
            fn = f"{ddir}/comparison-bed/{t1}-against-{t2}/summary.json"
            P2, R2, F2, TP2, FP2, FN2 = parse_summary(fn)

            ACC1 = round(TP1 / (TP1 + FP1 + FN1), 2)
            ACC2 = round(TP2 / (TP2 + FP2 + FN2), 2)

            if acc:
                M3[indexes[t2]][indexes[t1]] = (ACC1 + ACC2) / 2
                M3_labels[indexes[t2]][indexes[t1]] = f"def:bed\n{ACC1}:{ACC2}"
            else:
                M3[indexes[t2]][indexes[t1]] = (P1 + P2) / 2
                M3_labels[indexes[t2]][indexes[t1]] = f"def:bed\n{P1}:{P2}"

        sns.heatmap(
            M3,
            square=True,
            annot=M3_labels,
            fmt="",
            xticklabels=labels,
            cmap=CMAP,
            yticklabels=labels if i == 0 else False,
            cbar=False,
            vmin=0.5,
            vmax=1.0,
            ax=axes[i],
        )

    x_titles = ["T2T", "hg38", "hg19"]
    for ax, title in zip(axes, x_titles):
        ax.set_title(title)
        ax.set_xlabel("How much of this...")
    axes[0].set_ylabel("...matches with this")
    plt.tight_layout()
    plt.show()


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
            fn = f"{ddir}/comparison-def/{t1}-against-{t2}/summary.json"
            P1, R1, F1, TP1, FP1, FN1 = parse_summary(fn)
            fn = f"{ddir}/comparison-bed/{t1}-against-{t2}/summary.json"
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


def main():
    print("This script assumes <= 3 truth callsets!", file=sys.stderr)
    in_dir = sys.argv[1]
    run = "bed"  # sys.argv[2]

    labels = []  # ["dipcall", "svim-asm", "severus-paper"]
    indexes = {}  # {"dipcall": 0, "svim-asm": 1, "severus-paper": 2}
    for i, f in enumerate(glob.glob(f"{in_dir}/*.vcf.gz")):
        fname = f.split("/")[-1][:-7]
        labels.append(fname)
        indexes[fname] = i

    pairs = [
        ("dipcall", "svim-asm"),
        ("dipcall", "hapdiff"),
        ("svim-asm", "hapdiff"),
    ]
    M = [[0 for _ in labels] for _ in labels]
    two = False
    for t1, t2 in itertools.product(labels, labels):
        fn = f"{in_dir}/comparison-{run}/{t1}-against-{t2}/summary.json"
        if not os.path.exists(fn):
            print(f"Skipping {t1}/{t2}")
            two = True
            continue
        P, R, F, TP, FP, FN = parse_summary(fn)
        M[indexes[t1]][indexes[t2]] = P
        M[indexes[t2]][indexes[t1]] = R
        M[indexes[t1]][indexes[t1]] = 1
        M[indexes[t2]][indexes[t2]] = 1
    if two:
        M = [x[:2] for x in M[:2]]
        labels = labels[:2]
    sns.heatmap(
        M,
        square=True,
        annot=True,
        xticklabels=labels,
        yticklabels=labels,
        # vmin=0.7,
        # vmax=1.0,
    )
    plt.xlabel("How much of this...")
    plt.ylabel("...match with this")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if sys.argv[1] == "all":
        sys.argv.pop(0)
        main_all()
    elif sys.argv[1] == "avgp":
        sys.argv.pop(0)
        main_avg(acc=False)
    elif sys.argv[1] == "avga":
        sys.argv.pop(0)
        main_avg(acc=True)
    else:
        sys.argv.pop(0)
        main()
