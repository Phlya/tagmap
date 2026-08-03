from pathlib import Path
import time
import numpy as np

import bioframe
import pandas as pd
import argparse

from tagmaplib import cassette_orientation, normalise_junction_pos, read_pairs

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", "-i", type=str)
argparser.add_argument("--side", "-s", type=int, choices=[1, 2], default=1)
argparser.add_argument(
    "--itr-side",
    choices=["forward", "reverse"],
    default=None,
    help="Which ITR primer these pairs are anchored at. Given it, each position "
    "also gets the orientation its reads imply, which find_peaks.py carries "
    "through so that find_insertion_sites.py does not have to guess one from "
    "the order of the two sides.",
)
argparser.add_argument(
    "--end", choices=[5, 3, 0], default=0
)  # 5' end, 3' end, or whatever is reported as pos1 and pos2 in pairs [default: 0]
argparser.add_argument("--threads", "-t", type=int, default=1)
argparser.add_argument("--output", "-o", type=str)
argparser.add_argument("--output-bigwig", type=str, default=None, required=False)


# def pairs_to_pairs_merged(pairs_df):
#     pairs_bf1 = pairs_df[[c for c in pairs_df.columns if c.endswith("1")]]
#     pairs_bf1.columns = [c[:-1] for c in pairs_df.columns if c.endswith("1")]

#     pairs_bf2 = pairs_df[[c for c in pairs_df.columns if c.endswith("2")]]
#     pairs_bf2.columns = [c[:-1] for c in pairs_df.columns if c.endswith("2")]

#     pairs_df = pd.concat([pairs_bf1, pairs_bf2])

#     return pairs_df


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


if __name__ == "__main__":
    args = argparser.parse_args()

    pairs, chromsizes = read_pairs(args.input, threads=args.threads)

    if pairs.shape[0] == 0:
        Path(args.output).touch()
        if args.output_bigwig is not None:
            Path(args.output_bigwig).touch()
        exit()

    s = args.side
    e = args.end if args.end != 0 else ""
    position = pairs[f"pos{e}{s}"]
    if args.end != 3:
        # pos and pos5 both name the alignment's 5' end, which is the junction
        # the reads were sequenced outwards from - and which pairtools places
        # one base apart on the two strands.
        position = normalise_junction_pos(position, pairs[f"strand{s}"])
    pairs["start"] = position - 1
    pairs["end"] = position
    pairs[["start", "end"]] = np.sort(pairs[["start", "end"]], axis=1)
    pairs["chrom"] = pairs[f"chrom{s}"]
    # Reads run outwards from the cassette, from the junction towards the
    # tagmentation cut, so which way they run says which way the insertion
    # faces. Kept per position, because both ITR sides of one integration sit
    # on the same base and so cannot be told apart by their order.
    if args.itr_side is not None:
        pairs["rightwards"] = pairs[f"pos3{s}"] > pairs[f"pos5{s}"]
        facing = (
            pairs.groupby(["chrom", "start"])["rightwards"].mean().rename("facing")
        )

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
    coverage_df["fraction"] = coverage_df["count"] / coverage_df["count"].sum()

    if args.output_bigwig is not None:
        bioframe.to_bigwig(coverage_df, chromsizes, args.output_bigwig)

    if args.itr_side is not None:
        merged = coverage_df.merge(facing, on=["chrom", "start"], how="left")
        coverage_df["orientation"] = cassette_orientation(
            args.itr_side, merged["facing"].fillna(0) > 0.5
        )
    else:
        coverage_df["orientation"] = "."
    coverage_df.to_csv(args.output, sep="\t", index=False, header=False)
