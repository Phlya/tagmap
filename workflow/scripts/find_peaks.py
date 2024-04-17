from pathlib import Path

import pandas as pd
import bioframe
from skimage.filters import threshold_otsu

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", "-i", type=str)
argparser.add_argument("--min-peak-reads", type=int, default=10)
argparser.add_argument("--min-peak-frac", type=float, default=0.01)
argparser.add_argument("--min-peak-width", type=int, default=50)
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()

coverage = pd.read_csv(
    args.input,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "coverage", "fraction"],
)

# coverage_non_singleton = coverage[coverage["coverage"] > 1]
# if coverage_non_singleton.shape[0] == 0:
# Path(args.output).touch()
# exit()

if coverage.shape[0] == 0:
    Path(args.output).touch()
    exit()

merged = (
    bioframe.cluster(coverage, min_dist=args.min_peak_width * 10)
    .groupby("cluster")
    .agg({"chrom": lambda x: x.iloc[0], "start": "min", "end": "max"})
)

merged = (
    bioframe.overlap(merged, coverage, how="left")[
        ["chrom", "start", "end", "coverage_", "fraction_"]
    ]
    .groupby(["chrom", "start", "end"])
    .agg({"coverage_": "sum", "fraction_": "sum"})
    .reset_index()
    .rename(columns={"coverage_": "coverage", "fraction_": "fraction"})
)

merged = merged[merged["end"] - merged["start"] >= args.min_peak_width]
merged[
    (merged["fraction"] >= args.min_peak_frac)
    & (merged["coverage"] >= args.min_peak_reads)
].to_csv(args.output, sep="\t", header=False, index=False)
