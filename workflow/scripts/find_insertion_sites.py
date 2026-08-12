"""Turn peaks either side of an insertion into single-base integration sites.

A real integration is read from both ITR primers, so it shows up as a pair of
nearby peaks pointing at each other. Each peak is first snapped onto its
nearest motif (see --snap-window), so both sides of one real integration land
on the same base before --max-dist combines them; the pair's orientation
gives the strand, and the combined site is then re-pinpointed on the motif
for a final position and TA_found check.

Also splits out the sites seen from both sides with a resolved strand -
confirmed by NGS alone, with no reference to Sanger - into their own BED,
plus a copy of that with the construct/landing-pad contigs dropped.
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
argparser.add_argument(
    "--snap-window",
    type=int,
    default=5,
    help="Before combining forward/reverse peaks, snap each peak's own "
    "position to the nearest insertion_seq motif within this many bases, so "
    "peaks from one real integration land on the exact same base regardless "
    "of a few bp of slop (the +/- strand pos5 offset, indel wobble, "
    "Sleeping Beauty's own target-site duplication) - letting --max-dist "
    "combine sides on precise positions instead of tolerating that slop "
    "with a wide distance. Checked against real data: a peak's own raw "
    "position can be off by up to ~4bp even when correct. 0 disables "
    "snapping and combines on the raw peak positions instead.",
)
argparser.add_argument(
    "--min-orientation-support",
    type=int,
    default=0,
    help="Peaks with fewer than this many distinct molecules don't get a "
    "say in their cluster's orientation - at low support a peak's own "
    "stated orientation is little better than a coin flip, so trusting it "
    "risks a confident, wrong strand call. 0 disables the check.",
)
argparser.add_argument(
    "--chromsizes",
    default=None,
    help="Chrom sizes of the genome alone (chrom_sizes_path_no_cassette). "
    "Sites outside it - i.e. on the cassette/landing-pad contigs, which "
    "genome browsers don't know about - are dropped from "
    "--output-confirmed-no-cassette. Left unfiltered if not given.",
)
argparser.add_argument("--output", "-o", required=True)
argparser.add_argument("--output-for-ucsc", required=True)
argparser.add_argument(
    "--output-confirmed",
    required=True,
    help="Sites seen from both sides (site_sides == 'both') with a resolved "
    "strand (strand != '.') - confirmed by NGS alone, without reference to "
    "Sanger. Both-sided sites the two peaks couldn't agree a direction for "
    "are left out, not just unconfirmed ones.",
)
argparser.add_argument(
    "--output-confirmed-no-cassette",
    required=True,
    help="--output-confirmed with the construct/landing-pad contigs dropped "
    "(see --chromsizes).",
)


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


def cluster_orientation(group, min_support=0):
    """Orientation of one cluster of peaks, preferring what the reads said.

    A single peak's own molecule count is how much its stated orientation is
    trusted alone - checked against Sanger-confirmed sites, a thin peak's
    call is barely better than a coin flip. But two peaks that independently
    state the *same* orientation corroborate each other, so support is
    pooled per stated orientation before checking min_support, rather than
    requiring one single peak to individually clear the bar - one peak of
    32 molecules and another of 45 agreeing are together as trustworthy as
    one peak of 77.

    determine_direction is a *different*, weaker fallback, only for peak
    callers (e.g. "coverage") that never compute an orientation at all - it
    must not catch peaks that did state one but whose pooled support was
    still too thin to trust, since its own position-order guess has no
    support check of its own and can contradict what the (rejected-as-too-
    thin) peaks agreed on.
    """
    all_orientations = group.get("orientation", pd.Series(dtype=object))
    if not all_orientations.isin(["+", "-"]).any():
        return determine_direction(group["side"])
    stated = group[group["orientation"].isin(["+", "-"])]
    pooled = stated.groupby("orientation")["count"].sum()
    confident = pooled[pooled >= min_support] if min_support else pooled
    if len(confident) == 1:
        return confident.index[0]
    # Either multiple orientations each pooled enough support to be
    # confident - a real conflict - or none did - either way, better to say
    # nothing than fall back to a guess the peaks themselves didn't clear
    # the bar for.
    return "."


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
        tagmaplib.write_bed(empty, args.output_confirmed, columns=tagmaplib.SITE_COLUMNS)
        tagmaplib.write_bed(
            empty, args.output_confirmed_no_cassette, columns=tagmaplib.SITE_COLUMNS
        )
        raise SystemExit(0)

    peaks = (
        peaks.astype({"chrom": str, "start": int, "end": int, "count": int})
        .sort_values(["sample_name", "chrom", "start", "end"])
        .reset_index(drop=True)
    )
    if peaks["n_positions"].isnull().all():
        peaks["n_positions"] = 1

    if args.snap_window:
        # Every peak already has its own orientation ("+"/"-", never "." -
        # that only happens once sides are combined below), which is all
        # find_insertion_seq needs to not skip a row; insertion_seq (TA,
        # TTAA) is its own reverse complement for every transposon this
        # pipeline supports, so which strand is named here doesn't change
        # what gets matched.
        motif_start = tagmaplib.find_insertion_seq(
            peaks.rename(columns={"orientation": "strand"}),
            args.genome,
            ins_seq,
            window=args.snap_window,
            mode="nearest",
            index_file=args.genome_index,
        )
        found = motif_start >= 0
        site = tagmaplib.insertion_site_from_motif(motif_start, ins_seq)
        peaks.loc[found, "start"] = site[found]
        peaks.loc[found, "end"] = site[found] + 1

    peaks = bioframe.cluster(
        peaks,
        min_dist=args.max_dist,
        on=["sample_name"],
        return_cluster_ids=True,
    )

    orientations = peaks.groupby("cluster").apply(
        cluster_orientation,
        min_support=args.min_orientation_support,
        include_groups=False,
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
    site = tagmaplib.insertion_site_from_motif(motif_start, ins_seq)
    pinpointed = peaks.copy()
    pinpointed[f"{ins_seq}_found"] = found
    pinpointed.loc[found, "start"] = site[found]
    pinpointed.loc[found, "end"] = site[found] + 1
    pinpointed = pinpointed[tagmaplib.SITE_COLUMNS + [f"{ins_seq}_found", "site_sides"]]

    pinpointed.sort_values(["chrom", "start", "end", "sample_name"]).to_csv(
        args.output, sep="\t", index=False, header=True
    )

    pinpointed[tagmaplib.SITE_COLUMNS].sort_values(
        ["sample_name", "chrom", "start", "end"]
    ).to_csv(args.output_for_ucsc, sep="\t", index=False, header=False)

    confirmed = pinpointed.loc[
        (pinpointed["site_sides"] == "both") & (pinpointed["strand"] != "."),
        tagmaplib.SITE_COLUMNS,
    ].sort_values(["chrom", "start", "end", "sample_name"])
    tagmaplib.write_bed(confirmed, args.output_confirmed, columns=tagmaplib.SITE_COLUMNS)

    no_cassette = confirmed
    if args.chromsizes is not None:
        chromsizes = bioframe.read_chromsizes(args.chromsizes)
        no_cassette = bioframe.trim(confirmed, chromsizes).dropna()
        no_cassette[["start", "end"]] = no_cassette[["start", "end"]].astype(int)
    tagmaplib.write_bed(
        no_cassette, args.output_confirmed_no_cassette, columns=tagmaplib.SITE_COLUMNS
    )
