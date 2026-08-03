"""Merge the forward and reverse peaks of one sample into a single file."""

import argparse

import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--fwd", type=str)
argparser.add_argument("--rev", type=str)
argparser.add_argument("--blacklist", type=str, default=None)
argparser.add_argument("--sample-name", type=str)
argparser.add_argument("--output", "-o", type=str)
argparser.add_argument("--output-genome-browser", type=str)
args = argparser.parse_args()

PEAK_COLUMNS = [
    "chrom",
    "start",
    "end",
    "counts",
    "fraction",
    "n_positions",
    "orientation",
]

fwd = tagmaplib.read_peaks(args.fwd, columns=PEAK_COLUMNS)
fwd["side"] = "+"

rev = tagmaplib.read_peaks(args.rev, columns=PEAK_COLUMNS)
rev["side"] = "-"

merged = pd.concat([fwd, rev]).sort_values(["chrom", "start", "end"])
merged = tagmaplib.apply_blacklist(merged, args.blacklist)
merged["sample"] = args.sample_name

merged[
    [
        "chrom",
        "start",
        "end",
        "sample",
        "counts",
        "side",
        "fraction",
        "n_positions",
        "orientation",
    ]
].sort_values(["chrom", "start", "end"]).to_csv(
    args.output, sep="\t", header=False, index=False
)


def norm_counts(x):
    x["counts"] = (x["counts"] / x["counts"].sum() * 1000).astype(int)
    return x


merged[["chrom", "start", "end", "sample", "side", "counts"]].groupby("side").apply(
    norm_counts, include_groups=False
).reset_index()[
    ["chrom", "start", "end", "sample", "side", "counts"]
].sort_values(["chrom", "start", "end"]).to_csv(
    args.output_genome_browser, sep="\t", header=False, index=False
)
