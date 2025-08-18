import sys
import os
import argparse
import glob


def parse_summary(fpath):
    TP, FP, FN = 0, 0, 0
    for line in open(fpath):
        line = line.strip('\n \t"')
        if line.startswith('TP-base"'):
            if TP == 0:
                TP = int(line.split(" ")[1][:-1])
        elif line.startswith('FP"'):
            if FP == 0:
                FP = int(line.split(" ")[1][:-1])
        elif line.startswith('FN"'):
            if FN == 0:
                FN = int(line.split(" ")[1][:-1])
    return TP, FP, FN


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refine", action='store_true')
    parser.add_argument("indir")
    args = parser.parse_args()

    fname = "summary.json"
    if args.refine:
        fname = "ga4gh_with_refine.summary.json"

    print("Tool,TP,FP,FN,P,R,F1")
    for f in glob.glob(os.path.join(args.indir, "*", fname)):
        tool = f.split("/")[-2]
        tp, fp, fn = parse_summary(f)
        P = round(tp / (tp + fp) * 100, 1) if tp+fp > 0 else 0.0
        R = round(tp / (tp + fn) * 100, 1) if tp+fn > 0 else 0.0
        F = round(2 * (P * R) / (P + R) if P + R != 0 else 0.0, 1)
        print(tool, tp, fp, fn, P, R, F, sep=",")


if __name__ == "__main__":
    main()
