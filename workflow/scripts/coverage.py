import numpy as np

import bioframe
import pandas as pd
from pairtools.lib import fileio, headerops
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", "-i", type=str)
argparser.add_argument("--side", "-s", type=int, choices=[1, 2], default=1)
argparser.add_argument("--output", "-o", type=str)
argparser.add_argument("--output-bigwig", type=str, default=None, required=False)


def read_pairs(pairs, chunksize=None):
    pairs_stream = (
        fileio.auto_open(
            pairs,
            mode="r",
            nproc=10,
        )
        if isinstance(pairs, str)
        else pairs
    )

    header, pairs_body = headerops.get_header(pairs_stream)

    cols = headerops.extract_column_names(header)

    chromsizes = headerops.extract_chromsizes(header)

    pairs_df = pd.read_csv(
        pairs_body,
        header=None,
        names=cols,
        chunksize=None,
        sep="\t",
        dtype={"chrom1": str, "chrom2": str},
    )

    return pairs_df, chromsizes


def pairs_to_pairs_merged(pairs_df):
    pairs_bf1 = pairs_df[[c for c in pairs_df.columns if c.endswith("1")]]
    pairs_bf1.columns = [c[:-1] for c in pairs_df.columns if c.endswith("1")]

    pairs_bf2 = pairs_df[[c for c in pairs_df.columns if c.endswith("2")]]
    pairs_bf2.columns = [c[:-1] for c in pairs_df.columns if c.endswith("2")]

    pairs_df = pd.concat([pairs_bf1, pairs_bf2])

    return pairs_df


def intervals_to_increments(df):
    inc_df = pd.concat(
        [
            df[["chrom", "start"]].rename(columns={"start": "pos"}).eval("inc=1"),
            df[["chrom", "end"]].rename(columns={"end": "pos"}).eval("inc=-1"),
        ]
    )
    return inc_df


def aggregate_increments(inc_df):
    inc_df = (
        inc_df.groupby(["chrom", "pos"])
        .sum()
        .reset_index()
        .sort_values(["chrom", "pos"], ignore_index=True)
    )
    return inc_df


def coverage_single_chrom(chrom_df, chromsize):
    count_df = aggregate_increments(intervals_to_increments(chrom_df))

    count_df.insert(2, "count", count_df["inc"].cumsum())
    count_df.insert(2, "end", count_df["pos"].shift(-1, fill_value=chromsize))
    count_df.rename(columns={"pos": "start"}, inplace=True)

    # insert clean-up step here: make sure it starts with 0, doesn't contain zero-length intervals at the end
    if count_df.iloc[-1].start == chromsize:
        count_df.drop(count_df.index[-1], inplace=True)
    return count_df


args = argparser.parse_args()

pairs, chromsizes = read_pairs(args.input)

if pairs.shape[0] == 0:
    open(args.output, "w").close()
    exit()

s = args.side
pairs["start"] = pairs[f"pos3{s}"]
pairs["end"] = pairs[f"pos3{s}"] + np.sign(
    (pairs[f"strand{s}"] == "+").astype(int) - 0.5
).astype(int)
pairs[["start", "end"]] = np.sort(pairs[["start", "end"]], axis=1)
pairs["chrom"] = pairs[f"chrom{s}"]

pairs = pairs[["chrom", "start", "end"]]
pairs.sort_values(["chrom", "start", "end"], inplace=True)
pairs.reset_index(drop=True, inplace=True)

coverage_df = pd.concat(
    [
        coverage_single_chrom(chrom_reads, chromsizes[chrom])
        for chrom, chrom_reads in pairs.groupby("chrom")
    ]
).reset_index(drop=True)[["chrom", "start", "end", "count"]]
coverage_df = coverage_df[coverage_df["count"] > 0]

coverage_df.to_csv(args.output, sep="\t", index=False, header=False)

if args.output_bigwig is not None:
    bioframe.to_bigwig(coverage_df, chromsizes, args.output_bigwig)
