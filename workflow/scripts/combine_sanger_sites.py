"""Collapse per-read Sanger integration sites into one site per clone.

Several reads of the same clone - different wells, or the same well sequenced
from both ITR primers - should report the same integration site. Reads within
--max-dist of each other are grouped, and how much they agree is reported so
that a site seen by one read from one end can be told apart from one confirmed
from both.
"""

import argparse

import bioframe
import numpy as np
import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--sites", nargs="+", required=True)
argparser.add_argument("--max-dist", type=int, default=100)
argparser.add_argument("--blacklist", default=None)
argparser.add_argument("--output", "-o", required=True)
argparser.add_argument("--output-for-ucsc", required=True)

READ_SITE_COLUMNS = [
    "chrom",
    "start",
    "end",
    "sample_name",
    "mapq",
    "strand",
    "clone",
    "readname",
    "direction",
]

OUTPUT_COLUMNS = [
    "chrom",
    "start",
    "end",
    "sample_name",
    "score",
    "strand",
    "clone",
    "n_reads",
    "n_forward",
    "n_reverse",
    "both_directions",
    "strands_agree",
]


def consensus(group):
    """One site from the reads supporting it."""
    strands = group["strand"].unique()
    directions = group["direction"].value_counts()
    # Reads should agree exactly once snapped onto the insertion motif; take the
    # most common position so a single stray read cannot move the site.
    start = group["start"].mode().min()
    return pd.Series(
        {
            "chrom": group["chrom"].iloc[0],
            "start": int(start),
            "end": int(group.loc[group["start"] == start, "end"].iloc[0]),
            "sample_name": group["sample_name"].iloc[0],
            "clone": group["clone"].iloc[0],
            "strand": strands[0] if len(strands) == 1 else ".",
            "strands_agree": len(strands) == 1,
            "n_reads": group.shape[0],
            "n_forward": int(directions.get("forward", 0)),
            "n_reverse": int(directions.get("reverse", 0)),
        }
    )


if __name__ == "__main__":
    args = argparser.parse_args()

    sites = []
    for path in args.sites:
        sites.append(tagmaplib.read_peaks(path, columns=READ_SITE_COLUMNS))
    sites = pd.concat(sites, ignore_index=True)

    if sites.shape[0] == 0:
        print("No passing Sanger reads to combine")
        empty = pd.DataFrame({c: pd.Series(dtype=object) for c in OUTPUT_COLUMNS})
        empty.to_csv(args.output, sep="\t", index=False)
        empty[tagmaplib.SITE_COLUMNS].to_csv(
            args.output_for_ucsc, sep="\t", index=False, header=False
        )
        raise SystemExit(0)

    sites[["start", "end"]] = sites[["start", "end"]].astype(int)

    if args.blacklist is not None:
        before = sites.shape[0]
        sites = tagmaplib.apply_blacklist(sites, args.blacklist)
        print(f"Blacklist removed {before - sites.shape[0]} reads")

    clustered = bioframe.cluster(
        sites.sort_values(["sample_name", "clone", "chrom", "start"]),
        min_dist=args.max_dist,
        on=["sample_name", "clone"],
    )
    combined = (
        clustered.groupby("cluster", sort=True)
        .apply(consensus, include_groups=False)
        .reset_index(drop=True)
    )

    combined["both_directions"] = (combined["n_forward"] > 0) & (
        combined["n_reverse"] > 0
    )
    # BED scores are capped at 1000; a site with ten or more reads is as
    # convincing as it is going to get.
    combined["score"] = np.minimum(1000, combined["n_reads"] * 100)

    combined = combined[OUTPUT_COLUMNS].sort_values(
        ["sample_name", "clone", "chrom", "start"]
    )
    combined.to_csv(args.output, sep="\t", index=False)

    combined[tagmaplib.SITE_COLUMNS].sort_values(["chrom", "start", "end"]).to_csv(
        args.output_for_ucsc, sep="\t", index=False, header=False
    )

    print(
        f"{combined.shape[0]} sites in {combined['clone'].nunique()} clones, "
        f"{int(combined['both_directions'].sum())} seen from both ends"
    )
