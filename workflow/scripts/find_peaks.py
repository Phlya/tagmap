import pandas as pd
import bioframe
from skimage.filters import threshold_otsu

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", "-i", type=str)
argparser.add_argument("--min-peak-width", type=int, default=50)
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()

coverage = pd.read_csv(
    args.input,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "coverage", "fraction"],
)
coverage_non_singleton = coverage[coverage["coverage"] > 1]
merged = (
    bioframe.cluster(coverage_non_singleton, min_dist=500)
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

threshold = 0.001  # threshold_otsu(merged["coverage"].to_numpy())

merged = merged[merged["end"] - merged["start"] >= args.min_peak_width]
merged[merged["fraction"] >= threshold].to_csv(
    args.output, sep="\t", header=False, index=False
)
