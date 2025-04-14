import sys
import os
import glob

import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")


def parse_summary(fpath):
    P, R, F = 0, 0, 0
    for line in open(fpath):
        line = line.strip('\n \t"')
        if line.startswith('precision"'):
            P = round(float(line.split(" ")[1][:-1]), 2)
        elif line.startswith('recall"'):
            R = round(float(line.split(" ")[1][:-1]), 2)
        elif line.startswith('f1"'):
            F = round(float(line.split(" ")[1][:-1]), 2)
    return P, R, F


def main():
    print("This script assumes <= 3 truth callsets!", file=sys.stderr)
    in_dir = sys.argv[1]

    labels = []  # ["dipcall", "svim-asm", "severus-paper"]
    indexes = {}  # {"dipcall": 0, "svim-asm": 1, "severus-paper": 2}
    for i, f in enumerate(glob.glob(f"{in_dir}/*.vcf.gz")):
        fname = f.split("/")[-1][:-7]
        labels.append(fname)
        indexes[fname] = i

    pairs = [
        ("dipcall", "svim-asm"),
        ("dipcall", "severus-paper" if "severus-paper" in labels else "hapdiff"),
        ("svim-asm", "severus-paper" if "severus-paper" in labels else "hapdiff"),
    ]
    M = [[0 for _ in labels] for _ in labels]
    two = False
    for t1, t2 in pairs:
        fn = f"{in_dir}/comparison/{t1}-against-{t2}/summary.json"
        if not os.path.exists(fn):
            print(f"Skipping {t1}/{t2}")
            two = True
            continue
        P, R, F = parse_summary(fn)
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
    main()
