from pathlib import Path

import pandas as pd
import bioframe
from skimage.filters import threshold_li

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", "-i", type=str)
argparser.add_argument(
    "--auto-li",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Use Li's method to determine threshold of real peaks",
)
argparser.add_argument("--cluster", action=argparse.BooleanOptionalAction, default=True)
argparser.add_argument("--min-peak-reads", type=int, default=10)
argparser.add_argument("--min-peak-frac", type=float, default=0.01)
argparser.add_argument("--min-peak-width", type=int, default=50)
argparser.add_argument("--min-peak-positions", type=int, default=1)
argparser.add_argument("--min-peak-dist", type=int, default=1000)
argparser.add_argument("--ignore-chrom", type=str, default=None)
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

if args.ignore_chrom:
    coverage = coverage[coverage["chrom"] != args.ignore_chrom]

if args.cluster:
    merged = bioframe.merge(coverage, min_dist=args.min_peak_dist)
    coverage = (
        bioframe.overlap(merged, coverage, how="left")[
            ["chrom", "start", "end", "coverage_", "fraction_", "n_intervals"]
        ]
        .groupby(["chrom", "start", "end"])
        .agg({"coverage_": "sum", "fraction_": "sum", "n_intervals": "first"})
        .reset_index()
        .rename(
            columns={
                "fraction_": "fraction",
                "n_intervals": "n_positions",
            }
        )
    )
coverage = coverage.rename(columns={"coverage": "n_reads"})

if args.auto_li:
    threshold = threshold_li(coverage["n_reads"].to_numpy())
    coverage = coverage[coverage["n_reads"] >= threshold]
    coverage = bioframe.cluster(
        coverage, min_dist=1, return_cluster_ids=True, return_cluster_intervals=False
    )
    coverage = (
        coverage.sort_values(["chrom", "start", "end"])
        .groupby("cluster")
        .agg("first")
        .reset_index(drop=True)
    )

coverage = coverage[coverage["end"] - coverage["start"] >= args.min_peak_width]
if args.cluster:
    coverage = coverage[
        (coverage["fraction"] >= args.min_peak_frac)
        & (coverage["n_reads"] >= args.min_peak_reads)
        & (coverage["n_positions"] >= args.min_peak_positions)
    ]
else:
    coverage = coverage[
        (coverage["fraction"] >= args.min_peak_frac)
        & (coverage["n_reads"] >= args.min_peak_reads)
    ]

coverage.to_csv(args.output, sep="\t", header=False, index=False)
