"""Call peaks from pairs, positioning on junctions and counting on shear ends.

A TagMap fragment runs

    [tagmentation site] --genome-- [insertion site] --cassette-- [ITR primer]

and sequencing it yields two kinds of evidence that find_peaks.py, working off
a coverage track of pos1, cannot tell apart:

  * A pair whose read walked through the cassette/genome junction reports that
    junction as the genomic 5' end. Checked against Sanger-confirmed 2G6
    clones, this is exact to the base - median offset 1-3bp - so one such read
    already places an insertion, however rare it is in the pool.
  * A pair that only reaches the ITR primer never sees the junction. Its
    genomic 5' end is wherever the read happened to stop, a median ~130bp short
    of the true site with a long tail, so it says an insertion is somewhere
    within a fragment length but not where. Treating those as positions invents
    one peak per molecule: on the short-read 2G6 run it turned ~800 real sites
    into 30,053.

So junctions decide *where*, and need no corroboration; primer-anchored pairs
can only corroborate a junction, or - when no junction was sequenced at all -
form a low-confidence peak that still has to clear the usual thresholds. That
keeps rare, sequence-perfect singletons while staying usable on reads too short
to span the junction.

Abundance is the number of distinct tagmentation positions, not of reads. Every
molecule of one insertion shares the same junction but is cut at its own random
tagmentation site, so distinct far ends count independent captures whereas reads
count PCR - the same reasoning as the shear-site abundance estimators used for
retroviral integration sites.
"""

import argparse
from pathlib import Path

import bioframe
import numpy as np
import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--input", "-i", type=str, help="A per-side .pairs file")
argparser.add_argument(
    "--side",
    choices=["forward", "reverse"],
    required=True,
    help="Which ITR primer the pairs in --input are anchored at",
)
argparser.add_argument("--output", "-o", type=str)
argparser.add_argument(
    "--output-detail",
    type=str,
    default=None,
    help="Optional TSV with the per-peak evidence behind the call",
)
argparser.add_argument(
    "--junction-jitter",
    type=int,
    default=2,
    help="Junction positions this close are the same site (mapping wobble)",
)
argparser.add_argument(
    "--max-fragment",
    type=int,
    default=1000,
    help="Longest junction-to-tagmentation distance an anchored pair may span",
)
argparser.add_argument(
    "--min-junction-frags",
    type=int,
    default=1,
    help="Distinct tagmentation positions among junction reads to call a site",
)
argparser.add_argument(
    "--min-junction-mapq",
    type=int,
    default=0,
    help="Drop junction pairs below this MAPQ instead of trusting them alone",
)
argparser.add_argument(
    "--motif-max-support",
    type=int,
    default=0,
    help="Thinly supported junction peaks must carry the insertion motif",
)
argparser.add_argument("--insertion-seq", default="TA")
argparser.add_argument("--genome", default=None)
argparser.add_argument("--genome-index", default=None)
argparser.add_argument("--min-peak-reads", type=int, default=10)
argparser.add_argument("--min-peak-frac", type=float, default=0.01)
argparser.add_argument("--min-peak-width", type=int, default=1)
argparser.add_argument("--min-peak-positions", type=int, default=1)
argparser.add_argument("--min-peak-dist", type=int, default=1000)
argparser.add_argument(
    "--ignore-chrom",
    type=str,
    nargs="*",
    default=[],
    help="Contigs that are part of the construct, so never insertion sites",
)
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
EVIDENCE_COLUMNS = [
    "n_junction_reads",
    "n_junction_frags",
    "n_anchor_reads",
    "n_anchor_frags",
]
DETAIL_COLUMNS = PEAK_COLUMNS + ["tier"] + EVIDENCE_COLUMNS


def write_empty():
    Path(args.output).touch()
    if args.output_detail is not None:
        pd.DataFrame(columns=DETAIL_COLUMNS).to_csv(
            args.output_detail, sep="\t", index=False
        )


def assign_to_junctions(anchors, junction_positions):
    """Junction position each anchored pair belongs to, NaN where none fits.

    An anchored read stops short of the junction, so the site has to lie beyond
    its pos5 on the side away from the tagmentation end, and the fragment as a
    whole - junction to tagmentation site - has to stay within --max-fragment.
    Of the junctions satisfying both, the nearest is the only sensible choice,
    so this walks outwards from pos5 in the reading direction.
    """
    assigned = pd.Series(np.nan, index=anchors.index)
    for chrom, rows in anchors.groupby("chrom", sort=False):
        positions = junction_positions.get(chrom)
        if positions is None or len(positions) == 0:
            continue
        pos5 = rows["pos51"].to_numpy()
        pos3 = rows["pos31"].to_numpy()
        # Reading direction: from the tagmentation end towards the junction.
        direction = np.sign(pos5 - pos3)
        for step in (1, -1):
            facing = direction == step
            if not facing.any():
                continue
            if step == 1:
                index = np.searchsorted(positions, pos5[facing] - args.junction_jitter)
            else:
                index = (
                    np.searchsorted(
                        positions, pos5[facing] + args.junction_jitter, side="right"
                    )
                    - 1
                )
            inside = (index >= 0) & (index < len(positions))
            candidate = positions[np.clip(index, 0, len(positions) - 1)]
            fits = inside & (np.abs(candidate - pos3[facing]) <= args.max_fragment)
            assigned.loc[rows.index[facing][fits]] = candidate[fits]
    return assigned


pairs, chromsizes = tagmaplib.read_pairs(args.input)
if pairs.shape[0] == 0:
    write_empty()
    exit()

pairs = pairs.rename(columns={"chrom1": "chrom"})
pairs["pos51"] = tagmaplib.normalise_junction_pos(pairs["pos51"], pairs["strand1"])
# Reads run outwards from the cassette, from the junction end of the genomic
# segment towards the tagmentation cut, so this is which way the insertion
# faces - the only honest source of a strand once both ITR sides are pinned to
# the same base.
pairs["rightwards"] = pairs["pos31"] > pairs["pos51"]
if args.ignore_chrom:
    pairs = pairs[~pairs["chrom"].isin(args.ignore_chrom)]
if pairs.shape[0] == 0:
    write_empty()
    exit()

is_junction = pairs["walk_pair_type"].isin(tagmaplib.JUNCTION_WALK_TYPES)
if args.min_junction_mapq:
    # A singleton is only as good as its alignment, so a junction read that
    # could sit anywhere is dropped rather than demoted to an anchor: its pos5
    # is a real junction, which it would be a lie to reuse as a read that
    # stopped short of one.
    pairs = pairs[~is_junction | (pairs["mapq1"] >= args.min_junction_mapq)]
    is_junction = pairs["walk_pair_type"].isin(tagmaplib.JUNCTION_WALK_TYPES)

junctions = pairs[is_junction]
anchors = pairs[~is_junction].copy()

# --- tier 1: sites placed by a sequenced junction -------------------------
if junctions.shape[0]:
    sites = junctions[["chrom", "pos51", "pos31", "rightwards"]].copy()
    sites["start"] = sites["pos51"]
    sites["end"] = sites["pos51"] + 1
    # Collapse the couple of bases of mapping wobble around one real junction,
    # and nothing more - genuine neighbouring insertions can be tens of bases
    # apart in a local re-mobilisation assay.
    sites = bioframe.cluster(
        sites,
        min_dist=args.junction_jitter,
        return_cluster_ids=True,
        return_cluster_intervals=False,
    )
    # The best-supported position in a cluster represents it; ties go leftwards
    # so the result never depends on row order.
    representative = (
        sites.groupby(["cluster", "pos51"])
        .size()
        .reset_index(name="n")
        .sort_values(["cluster", "n", "pos51"], ascending=[True, False, True])
        .groupby("cluster")["pos51"]
        .first()
    )
    sites["pos"] = sites["cluster"].map(representative)
    tier1 = (
        sites.groupby(["chrom", "pos"])
        .agg(
            n_junction_reads=("pos31", "size"),
            n_junction_frags=("pos31", "nunique"),
            n_positions=("pos51", "nunique"),
            rightwards=("rightwards", "mean"),
        )
        .reset_index()
    )
else:
    tier1 = pd.DataFrame(
        columns=[
            "chrom",
            "pos",
            "n_junction_reads",
            "n_junction_frags",
            "n_positions",
            "rightwards",
        ]
    )
tier1["pos"] = tier1["pos"].astype(float)

# --- anchored pairs: corroborate a junction, or stand on their own --------
junction_positions = {
    chrom: np.sort(rows["pos"].to_numpy())
    for chrom, rows in tier1.groupby("chrom", sort=False)
}
anchors["pos"] = (
    assign_to_junctions(anchors, junction_positions)
    if anchors.shape[0] and junction_positions
    else np.nan
)

support = (
    anchors[anchors["pos"].notna()]
    .groupby(["chrom", "pos"])
    .agg(n_anchor_reads=("pos31", "size"), n_anchor_frags=("pos31", "nunique"))
    .reset_index()
)
tier1 = tier1.merge(support, on=["chrom", "pos"], how="left")
tier1["width"] = 1

# --- tier 2: nothing but anchored pairs -----------------------------------
orphans = anchors[anchors["pos"].isna()].copy()
if orphans.shape[0]:
    # Without a junction the site is only known to lie beyond pos5, so the read
    # that got closest to it is the best estimate the data supports. Reads
    # facing opposite ways belong to different insertions, so they are
    # clustered apart.
    orphans["direction"] = np.sign(orphans["pos51"] - orphans["pos31"])
    orphans["start"] = orphans["pos51"]
    orphans["end"] = orphans["pos51"] + 1
    clustered = []
    for direction, rows in orphans.groupby("direction", sort=False):
        rows = bioframe.cluster(
            rows,
            min_dist=args.min_peak_dist,
            return_cluster_ids=True,
            return_cluster_intervals=False,
        )
        rows["cluster"] = f"{direction}_" + rows["cluster"].astype(str)
        clustered.append(rows)
    orphans = pd.concat(clustered)
    closest = orphans.groupby("cluster").apply(
        lambda g: g["pos51"].max() if g["direction"].iloc[0] > 0 else g["pos51"].min(),
        include_groups=False,
    )
    orphans["pos"] = orphans["cluster"].map(closest).astype(float)
    tier2 = (
        orphans.groupby(["chrom", "cluster"])
        .agg(
            pos=("pos", "first"),
            n_anchor_reads=("pos31", "size"),
            n_anchor_frags=("pos31", "nunique"),
            n_positions=("pos51", "nunique"),
            rightwards=("rightwards", "mean"),
            lowest=("pos51", "min"),
            highest=("pos51", "max"),
        )
        .reset_index()
    )
    tier2["width"] = tier2["highest"] - tier2["lowest"] + 1
    tier2["n_junction_reads"] = 0
    tier2["n_junction_frags"] = 0
    tier2 = tier2.drop(columns=["cluster", "lowest", "highest"])
else:
    tier2 = pd.DataFrame(columns=tier1.columns)

peaks = pd.concat([tier1, tier2], ignore_index=True)
if peaks.shape[0] == 0:
    write_empty()
    exit()

peaks[EVIDENCE_COLUMNS] = peaks[EVIDENCE_COLUMNS].fillna(0)
peaks["tier"] = np.where(peaks["n_junction_frags"] > 0, "junction", "anchor")
peaks["orientation"] = tagmaplib.cassette_orientation(
    args.side, peaks["rightwards"] > 0.5
)
# Independent molecules, not reads: PCR copies of one fragment share both ends,
# so only distinct tagmentation positions say how often a site was captured.
peaks["counts"] = peaks["n_junction_frags"] + peaks["n_anchor_frags"]
peaks["fraction"] = peaks["counts"] / peaks["counts"].sum()
peaks["pos"] = peaks["pos"].astype(int)
# Same half-open, 1bp convention as the coverage track find_peaks.py reads.
peaks["start"] = peaks["pos"] - 1
peaks["end"] = peaks["pos"]

# A read that walked through a junction really did walk through one, but at
# this depth a lone one is usually a chimeric PCR product or a mismapping
# rather than a rare integration: on the 2G6 long-read library, singleton
# junctions carry a TA at the junction only 48% of the time against a 23%
# background, so barely a third of them can be genuine, while junctions seen on
# four or more molecules are at 90%. Asking a thinly supported peak to look
# like an integration - the motif the transposon actually inserts into, right
# at the junction - is the cheap way to keep the rare ones without the noise,
# and needs no labelled data to calibrate.
has_motif = pd.Series(True, index=peaks.index)
if args.motif_max_support and args.genome:
    thin = peaks["n_junction_frags"] <= args.motif_max_support
    candidates = peaks.loc[thin].assign(strand="+")
    if candidates.shape[0]:
        motif_start = tagmaplib.find_insertion_seq(
            candidates,
            args.genome,
            args.insertion_seq,
            window=args.junction_jitter,
            mode="nearest",
            index_file=args.genome_index,
        )
        has_motif.loc[thin] = (motif_start >= 0).to_numpy()

# A sequenced junction is its own evidence, so it only has to clear the
# junction threshold. Anything resting on anchored pairs alone is a guess at a
# position, and has to earn its place the way find_peaks.py makes every peak.
keep = np.where(
    peaks["tier"] == "junction",
    (peaks["n_junction_frags"] >= args.min_junction_frags) & has_motif,
    (peaks["counts"] >= args.min_peak_reads)
    & (peaks["fraction"] >= args.min_peak_frac)
    & (peaks["n_positions"] >= args.min_peak_positions)
    & (peaks["width"] >= args.min_peak_width),
)
peaks = peaks[keep].sort_values(["chrom", "start", "end"])
for column in ["counts", "n_positions"] + EVIDENCE_COLUMNS:
    peaks[column] = peaks[column].astype(int)

peaks[PEAK_COLUMNS].to_csv(args.output, sep="\t", header=False, index=False)
if args.output_detail is not None:
    peaks[DETAIL_COLUMNS].to_csv(args.output_detail, sep="\t", index=False)
