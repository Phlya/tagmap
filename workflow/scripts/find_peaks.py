import pandas as pd
import bioframe
from skimage.filters import threshold_li

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", "-i", type=str)
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()

coverage = pd.read_csv(
    args.input, sep="\t", header=None, names=["chrom", "start", "end", "coverage"]
)
coverage = coverage[coverage["coverage"] > 0]
merged = (
    bioframe.cluster(coverage)
    .groupby("cluster")
    .agg(
        {"chrom": lambda x: x.iloc[0], "start": "min", "end": "max", "coverage": "sum"}
    )
)
threshold = threshold_li(merged["coverage"].to_numpy())
merged[merged["coverage"] >= threshold].to_csv(
    args.output, sep="\t", header=False, index=False
)
