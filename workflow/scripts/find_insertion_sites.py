"""Turn peaks either side of an insertion into single-base integration sites.

A real integration is read from both ITR primers, so it shows up as a pair of
nearby peaks pointing at each other. Peaks within --max-dist are combined, the
pair's orientation gives the strand, and the site is then pinpointed on the
motif the transposon integrates into.
"""

import argparse

import numpy as np
import pandas as pd
import bioframe

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--peaks", required=True)
argparser.add_argument("--max-dist", type=int, default=100)
argparser.add_argument("--genome", "-g", required=True)
argparser.add_argument(
    "--genome-index",
    default=None,
    help="pyfastx index of --genome. Built next to the FASTA if not given.",
)
argparser.add_argument("--insertion-seq", default="TA")
argparser.add_argument("--output", "-o", required=True)
argparser.add_argument("--output-for-ucsc", required=True)


def determine_direction(series):
    """Orientation of an insertion from the order of its two peaks.

    A last resort, for peaks from a caller that could not say which way the
    reads ran. It only works while the two sides land on different bases, which
    they need not: both ITRs of one integration sit at the same motif, so a
    caller that positions on the sequenced junction puts them on the very same
    base and the order here becomes a coin toss.
    """
    if series.shape[0] != 2:
        return "."
    if np.all(series.to_numpy() == np.asarray(["+", "-"])):
        return "-"
    elif np.all(series.to_numpy() == np.asarray(["-", "+"])):
        return "+"
    else:
        return "."


def cluster_orientation(group):
    """Orientation of one cluster of peaks, preferring what the reads said."""
    stated = {o for o in group.get("orientation", []) if o in ("+", "-")}
    if len(stated) == 1:
        return stated.pop()
    if len(stated) > 1:
        # The two sides of one insertion cannot face opposite ways, so this is
        # two insertions merged, or noise. Better to say nothing.
        return "."
    return determine_direction(group["side"])


def determine_sides(series):
    """Which side(s) of the cassette a site's peaks came from.

    Independent of determine_direction: a cluster can have peaks from both
    sides without their order supporting a confident strand call, and that is
    still "both", not one-sided.
    """
    sides = set(series)
    if sides == {"+", "-"}:
        return "both"
    elif sides == {"+"}:
        return "forward_only"
    elif sides == {"-"}:
        return "reverse_only"
    else:
        return "unknown"


if __name__ == "__main__":
    args = argparser.parse_args()
    ins_seq = args.insertion_seq

    peaks = tagmaplib.read_peaks(args.peaks)
    if peaks.shape[0] == 0:
        print(f"No peaks in {args.peaks}, so there are no insertion sites")
        columns = tagmaplib.SITE_COLUMNS + [f"{ins_seq}_found", "site_sides"]
        empty = pd.DataFrame({c: pd.Series(dtype=object) for c in columns})
        empty.to_csv(args.output, sep="\t", index=False, header=True)
        empty[tagmaplib.SITE_COLUMNS].to_csv(
            args.output_for_ucsc, sep="\t", index=False, header=False
        )
        raise SystemExit(0)

    peaks = (
        peaks.astype({"chrom": str, "start": int, "end": int, "count": int})
        .sort_values(["sample_name", "chrom", "start", "end"])
        .reset_index(drop=True)
    )
    if peaks["n_positions"].isnull().all():
        peaks["n_positions"] = 1

    peaks = bioframe.cluster(
        peaks,
        min_dist=args.max_dist,
        on=["sample_name"],
        return_cluster_ids=True,
    )

    orientations = peaks.groupby("cluster").apply(
        cluster_orientation, include_groups=False
    )
    peaks["strand"] = peaks["cluster"].map(orientations)
    peaks["site_sides"] = peaks.groupby(["cluster"])["side"].transform(determine_sides)

    peaks["start"] = peaks["cluster_start"]
    peaks["end"] = peaks["cluster_end"]

    peaks = peaks[
        ["chrom", "start", "end", "sample_name", "fraction", "strand", "site_sides"]
    ]
    peaks = (
        peaks.groupby(
            ["chrom", "start", "end", "sample_name", "strand", "site_sides"]
        )["fraction"]
        .mean()
        .reset_index()
    )
    peaks["score"] = (peaks["fraction"] * 1000).round().astype(int)
    peaks = peaks.drop_duplicates().reset_index(drop=True)
    peaks = bioframe.expand(peaks, len(ins_seq))

    motif_start = tagmaplib.find_insertion_seq(
        peaks, args.genome, ins_seq, window=0, mode="first", index_file=args.genome_index
    )
    found = motif_start >= 0
    pinpointed = peaks.copy()
    pinpointed[f"{ins_seq}_found"] = found
    pinpointed.loc[found, "start"] = motif_start[found]
    pinpointed.loc[found, "end"] = motif_start[found] + 1
    pinpointed = pinpointed[tagmaplib.SITE_COLUMNS + [f"{ins_seq}_found", "site_sides"]]

    pinpointed.sort_values(["chrom", "start", "end", "sample_name"]).to_csv(
        args.output, sep="\t", index=False, header=True
    )

    pinpointed[tagmaplib.SITE_COLUMNS].sort_values(
        ["sample_name", "chrom", "start", "end"]
    ).to_csv(args.output_for_ucsc, sep="\t", index=False, header=False)
