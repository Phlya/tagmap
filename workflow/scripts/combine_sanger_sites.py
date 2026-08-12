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
argparser.add_argument(
    "--chromsizes",
    default=None,
    help="Chrom sizes of the genome alone (chrom_sizes_path_no_cassette). "
    "Sites outside it - i.e. on the cassette/landing-pad contigs, which "
    "genome browsers don't know about - are dropped from --output-for-ucsc. "
    "Left unfiltered if not given.",
)

UCSC_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]

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
        empty["name"] = empty["sample_name"]
        tagmaplib.write_bed(empty, args.output_for_ucsc, columns=UCSC_COLUMNS)
        raise SystemExit(0)

    sites[["start", "end"]] = sites[["start", "end"]].astype(int)

    # The two ITR primers read outwards in opposite directions, so a forward-
    # primer read and a reverse-primer read of the same real insertion land
    # on opposite raw genomic strands (see sanger_sites.py, which reports
    # each read's actual mapped strand, uncorrected). Forward is the
    # pipeline's reference; reverse is flipped here so both sides of one
    # clone report the same final strand and strands_agree means something.
    # Matches tagmaplib.cassette_orientation's own convention on the NGS
    # side, so both branches call the same insertion the same way.
    sites["strand"] = np.where(
        sites["direction"] == "reverse",
        tagmaplib.flip_strand(sites["strand"]),
        sites["strand"],
    )

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

    # sample_name alone repeats across every clone on a plate; plate+clone
    # actually identifies which well a browser hit came from.
    combined["name"] = combined["sample_name"] + "_" + combined["clone"]
    for_ucsc = combined[UCSC_COLUMNS].sort_values(["chrom", "start", "end"])
    if args.chromsizes is not None:
        chromsizes = bioframe.read_chromsizes(args.chromsizes)
        for_ucsc = bioframe.trim(for_ucsc, chromsizes).dropna()
        for_ucsc[["start", "end"]] = for_ucsc[["start", "end"]].astype(int)
    track_name = (
        "_".join(sorted(combined["sample_name"].unique())) if combined.shape[0] else None
    )
    tagmaplib.write_bed(
        for_ucsc, args.output_for_ucsc, track_name=track_name, color_by_strand=True
    )

    print(
        f"{combined.shape[0]} sites in {combined['clone'].nunique()} clones, "
        f"{int(combined['both_directions'].sum())} seen from both ends"
    )
