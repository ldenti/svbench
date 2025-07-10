import sys
import glob
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# sns.set(font_scale=2)
sns.set(style="whitegrid")

REFSEQS = ["hg19", "hg38", "T2T"]
TRUTHS = ["dipcall", "svim-asm", "hapdiff"]


def parse_dir(ddir, refseq=""):
    # we may not need this dict and do everything on df but it's more convenient to me
    truths = {}
    data = []
    for vcf_fn in glob.glob(f"{ddir}/*.vcf.gz"):
        if "noinfo" in vcf_fn:
            continue
        name = vcf_fn.split("/")[-1].split(".")[0]
        if name not in TRUTHS:
            continue
        truths[name] = {}
        for record in VariantFile(vcf_fn):
            l = len(record.alts[0]) - len(record.ref)
            # filters=[x.name for x in record.filter.values()]
            # if name == "dipcall" and len(filters) > 0 and ("GAP1" in filters or "GAP2" in filters):
            #     continue
            # elif name != "dipcall" and "PASS" not in filters:
            #     continue

            gt1, gt2 = record.samples[0]["GT"]
            gt1 = gt1 if gt1 != None else 0
            gt2 = gt2 if gt2 != None else 0

            if abs(l) >= 50 and (gt1 > 0 or gt2 > 0):
                data.append(
                    [
                        refseq,
                        name,
                        record.contig,
                        record.pos,
                        l,
                        "INS" if l > 0 else "DEL",
                        gt1 == gt2,
                    ]
                )
                if record.contig not in truths[name]:
                    truths[name][record.contig] = []
                truths[name][record.contig].append(record.pos)
    return data, truths


def parse_pafs(ddir, refseq=""):
    data = []
    for paf_fn in glob.glob(f"{ddir}/ttmars-like/*/haps.paf"):
        name = paf_fn.split("/")[-2]
        if name not in TRUTHS:
            continue
        for line in open(paf_fn):
            line = line.strip("\n").split("\t")
            nm = int(line[12].split(":")[-1])
            tp = line[16][-1]
            de = float(line[20].split(":")[-1])
            if tp == "P":
                data.append([refseq, name, nm, de])
    return data


def main_gt():
    nb_dist = 50

    fig, axes = plt.subplots(2, 3, figsize=(10, 5))
    for i, vcf_dir in enumerate(sys.argv[1:4]):
        # we may not need this dict and do everything on df but it's more convenient to me
        df, truths = parse_dir(vcf_dir)
        df = pd.DataFrame(
            df, columns=["x", "Truth", "Chrom", "Pos", "Len", "Type", "GT"]
        )

        # GT distribution
        gtdf = []
        for truth in truths:
            for b in [True, False]:
                gt = "1/1" if b else "0/1"
                gtdf.append(
                    [truth, gt, len(df[(df["Truth"] == truth) & (df["GT"] == b)])]
                )
        gtdf = pd.DataFrame(gtdf, columns=["Truth", "GT", "Count"])

        sns.barplot(
            gtdf,
            x="Truth",
            order=TRUTHS,
            y="Count",
            hue="GT",
            palette="Set2",
            ax=axes[0][i],
            legend=True if i == 2 else None,
        )
        axes[0][i].set_ylim(0, 31000)
        axes[0][i].set_xlabel("")
        axes[0][i].set_ylabel("")  # Count
        if i == 0:
            axes[0][i].set_ylabel("(a)")
        # if i != 0:
        #     # remove y-ticks
        #     axes[0][i].set_yticklabels([])
        if i == 2:
            # move legends
            sns.move_legend(axes[0][i], "center left", bbox_to_anchor=(1, 0.5))

        # ax1.bar_label(ax1.containers[0])
        # ax1.bar_label(ax1.containers[1])
        axes[0][i].set_title(REFSEQS[i])

        # neighbor distribution per truthset
        df2 = []
        for truth in truths:
            neighbors = []
            for chrom in truths[truth]:
                last_p = truths[truth][chrom][0]
                neighbors.append(0)
                for p in truths[truth][chrom][1:]:
                    if p - last_p > nb_dist:
                        neighbors.append(0)
                    else:
                        neighbors[-1] += 1
                    last_p = p
            d = {}
            for x in neighbors:
                k = str(x)
                if x >= 2:
                    k = "2+"
                d[k] = d[k] + 1 if k in d else 1
            for k, v in d.items():
                df2.append([truth, k, v])

        df2 = pd.DataFrame(df2, columns=["Truth", f"#Neighbors-{nb_dist}bp", "Count"])
        sns.barplot(
            data=df2,
            x=f"#Neighbors-{nb_dist}bp",
            order=["0", "1", "2+"],
            y="Count",
            hue="Truth",
            ax=axes[1][i],
            legend=True if i == 2 else None,
        )
        axes[1][i].set_ylim(0, 35000)
        axes[1][i].set_ylabel("")  # Count
        if i == 0:
            axes[1][i].set_ylabel("(b)")
        if i == 2:
            # move legends
            sns.move_legend(axes[1][i], "center left", bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.show()
    # plt.savefig(vcf_dir + "/stats.png")


def main_distr():
    t2t_ddir = sys.argv[1]
    hg38_ddir = sys.argv[2]
    hg19_ddir = sys.argv[3]

    df = []
    df += parse_dir(t2t_ddir, "t2t")[0]
    df += parse_dir(hg38_ddir, "hg38")[0]
    df += parse_dir(hg19_ddir, "hg19")[0]

    df = pd.DataFrame(
        df, columns=["RefSeq", "Truth", "Chrom", "Pos", "Len", "Type", "GT"]
    )

    fig, axes = plt.subplots(2, 3, figsize=(10, 5))

    for i, refseq in enumerate(REFSEQS):
        subdf = df[df["RefSeq"] == refseq]
        df2 = subdf.groupby(["Truth", "Type"]).count()
        sns.barplot(
            data=df2,
            x="Truth",
            order=TRUTHS,
            y="Chrom",
            hue="Type",
            legend=True if i == 2 else None,
            palette="Set2",
            ax=axes[0][i],
        )
        axes[0][i].set_xlabel("")  # Truth
        axes[0][i].set_ylim(0, 19000)  # Count
        axes[0][i].set_ylabel("")  # Count
        if i == 0:
            axes[0][i].set_ylabel("(a)")
        if i != 0:
            # remove y-ticks
            axes[0][i].set_yticklabels([])
        if i == 2:
            # move legends
            sns.move_legend(axes[0][i], "center left", bbox_to_anchor=(1, 0.5))

        # ax1.bar_label(ax1.containers[0])
        # ax1.bar_label(ax1.containers[1])
        axes[0][i].set_title(refseq)

        # variant length distribution per truthset
        sns.histplot(
            data=subdf[abs(subdf["Len"]) <= 500],
            x="Len",
            hue="Truth",
            element="poly",
            fill=False,
            legend=True if i == 2 else None,
            ax=axes[1][i],
        )
        axes[1][i].set_xlabel("Length")
        axes[1][i].set_ylabel("")  # Count
        axes[1][i].set_ylim(0, 2000)
        if i == 0:
            axes[1][i].set_ylabel("(b)")
        if i != 0:
            # remove y-ticks
            axes[1][i].set_yticklabels([])
        if i == 2:
            sns.move_legend(axes[1][i], "center left", bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.show()


def main():
    t2t_ddir = sys.argv[1]
    hg38_ddir = sys.argv[2]
    hg19_ddir = sys.argv[3]
    nb_dist = 500

    df = []
    d, t2t_truth = parse_dir(t2t_ddir, "T2T")
    df += d
    d, hg38_truth = parse_dir(hg38_ddir, "hg38")
    df += d
    d, hg19_truth = parse_dir(hg19_ddir, "hg19")
    df += d

    df = pd.DataFrame(
        df, columns=["RefSeq", "Truth", "Chrom", "Pos", "Len", "Type", "GT"]
    )

    fig, axes = plt.subplots(4, 3, figsize=(9, 11))

    for i, refseq in enumerate(REFSEQS):
        subdf = df[df["RefSeq"] == refseq]
        df2 = subdf.groupby(["Truth", "Type"]).count()
        sns.barplot(
            data=df2,
            x="Truth",
            order=TRUTHS,
            y="Chrom",
            hue="Type",
            legend=True if i == 2 else None,
            palette="Set2",
            ax=axes[0][i],
        )
        axes[0][i].set_xlabel("")  # Truth
        axes[0][i].tick_params(axis="x", labelrotation=20)
        axes[0][i].set_ylim(0, 19000)  # Count
        axes[0][i].set_ylabel("")  # Count
        if i == 0:
            axes[0][i].set_ylabel("(a)")
        if i != 0:
            # remove y-ticks
            axes[0][i].set_yticklabels([])
        if i == 2:
            # move legends
            sns.move_legend(axes[0][i], "center left", bbox_to_anchor=(1, 0.5))

        # ax1.bar_label(ax1.containers[0])
        # ax1.bar_label(ax1.containers[1])
        axes[0][i].set_title(refseq)

        # variant length distribution per truthset
        sns.histplot(
            data=subdf[abs(subdf["Len"]) <= 500],
            x="Len",
            hue="Truth",
            hue_order=TRUTHS,
            element="poly",
            fill=False,
            legend=True if i == 2 else None,
            ax=axes[1][i],
        )
        axes[1][i].set_xlabel("Length")
        axes[1][i].set_xticks([-500, -250, 0, 250, 500])
        axes[1][i].set_ylabel("")  # Count
        axes[1][i].set_ylim(0, 2000)
        if i == 0:
            axes[1][i].set_ylabel("(b)")
        if i != 0:
            # remove y-ticks
            axes[1][i].set_yticklabels([])
        if i == 2:
            sns.move_legend(axes[1][i], "center left", bbox_to_anchor=(1, 0.5))

    # neighbor distribution per truthset
    for i, truths in enumerate([hg19_truth, hg38_truth, t2t_truth]):
        df2 = []
        for truth in truths:
            neighbors = []
            for chrom in truths[truth]:
                last_p = truths[truth][chrom][0]
                neighbors.append(0)
                for p in truths[truth][chrom][1:]:
                    if p - last_p > nb_dist:
                        neighbors.append(0)
                    else:
                        neighbors[-1] += 1
                    last_p = p
            d = {}
            for x in neighbors:
                k = str(x)
                if x >= 2:
                    k = "2+"
                d[k] = d[k] + 1 if k in d else 1
            for k, v in d.items():
                df2.append([truth, k, v])

        df2 = pd.DataFrame(df2, columns=["Truth", f"#Neighbors-{nb_dist}bp", "Count"])
        sns.barplot(
            data=df2,
            x=f"#Neighbors-{nb_dist}bp",
            order=["0", "1", "2+"],
            y="Count",
            hue="Truth",
            hue_order=TRUTHS,
            ax=axes[2][i],
            legend=True if i == 2 else None,
        )
        axes[2][i].set_ylim(0, 30000)
        axes[2][i].set_ylabel("")  # Count
        if i == 0:
            axes[2][i].set_ylabel("(c)")
        if i == 2:
            # move legends
            sns.move_legend(axes[2][i], "center left", bbox_to_anchor=(1, 0.5))

    # --- NM ttmars-like
    df = []
    df += parse_pafs(t2t_ddir, "T2T")
    df += parse_pafs(hg38_ddir, "hg38")
    df += parse_pafs(hg19_ddir, "hg19")

    df = pd.DataFrame(df, columns=["RefSeq", "Truth", "NM", "de"])
    for i, refseq in enumerate(REFSEQS):
        subdf = df[df["RefSeq"] == refseq]
        sns.histplot(
            data=subdf,
            x="NM",
            hue="Truth",
            hue_order=TRUTHS,
            binrange=[0, 50],
            discrete=True,
            element="step",
            legend=True if i == 2 else None,
            ax=axes[3][i],
        )
        axes[3][i].set_ylim(0, 4000)
        if i == 0:
            axes[3][i].set_ylabel("(d)")
        else:
            axes[3][i].set_ylabel("")
        if i == 2:
            # move legends
            sns.move_legend(axes[3][i], "center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()
    # plt.savefig("x.pdf")


if __name__ == "__main__":
    # mode = sys.argv.pop(1)
    # if mode == "distr":
    #     main_distr()
    # elif mode == "gt":
    #     main_gt()
    # else:
    main()
