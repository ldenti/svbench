import sys
import pysam
import pandas as pd

# import seaborn as sns
# import matplotlib.pyplot as plt


def main():
    mode = sys.argv[1]
    in_fn = sys.argv[2]

    data = []
    if mode == "fai":
        for line in open(in_fn):
            l = int(line.split("\t")[1])
            data.append([l, True])
    else:
        bam = pysam.AlignmentFile(in_fn, "rb")
        for aln in bam:
            if aln.is_unmapped:
                continue
            primary = True
            if aln.is_supplementary or aln.is_secondary:
                primary = False
            l = aln.query_length
            s = aln.get_cigar_stats()[0][4]
            data.append([l - s, primary])

    print(len(data))
    for primary in [1, 0]:
        short = sum([l < 100000 for l, p in data if p == primary])
        total = len([l for l, p in data if p == primary])
        print(primary, short, total, short / total if total > 0 else 0)

    df = pd.DataFrame(data, columns=["Length", "Primary"])
    print(df.describe())

    # sns.boxplot(df, x="Length", hue="Primary")
    # plt.show()


if __name__ == "__main__":
    main()
