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
    args.input, sep="\t", header=None, names=["chrom", "start", "end", "coverage"]
)
coverage_non_singleton = coverage[coverage["coverage"] > 1]
merged = (
    bioframe.cluster(coverage_non_singleton, min_dist=50)
    .groupby("cluster")
    .agg(
        {"chrom": lambda x: x.iloc[0], "start": "min", "end": "max", "coverage": "sum"}
    )
)
threshold = threshold_otsu(merged["coverage"].to_numpy())

merged = merged[merged["end"] - merged["start"] >= args.min_peak_width]
merged[merged["coverage"] >= threshold].to_csv(
    args.output, sep="\t", header=False, index=False
)
