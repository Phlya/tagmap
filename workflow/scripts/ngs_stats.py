"""Summarise raw NGS usability and insertion-site sidedness per library.

One row per NGS library, combining pre-deduplication pair classes and support
for retained peaks/sites with how many resulting insertion sites were seen
from both sides of the cassette versus only one.
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
import yaml

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--stats-yml", nargs="*", default=[])
argparser.add_argument(
    "--pairs",
    nargs="*",
    default=[],
    help="{sample}_forward.pairs/{sample}_reverse.pairs files",
)
argparser.add_argument("--raw-pairs", nargs="*", default=[])
argparser.add_argument(
    "--dedup-pairs",
    nargs="*",
    default=[],
    help="{sample}_dupmarked.pairs (pair_type DD for duplicates, otherwise "
    "unchanged) if the workflow deduped, else the same {sample}_sorted.pairs "
    "as --raw-pairs - denominator for frac_mobilized.",
)
argparser.add_argument("--peaks", required=True, help="all_peaks.bed")
argparser.add_argument("--sites", required=True, help="all_sites.bed")
argparser.add_argument("--primer-positions", required=True)
argparser.add_argument("--cassette-name", required=True)
argparser.add_argument("--cassette-length", type=int, required=True)
argparser.add_argument("--only-junctions", action="store_true")
argparser.add_argument("--output", "-o", required=True)
args = argparser.parse_args()

# Column groups shared with report.py/report_pdf.py, which split this
# script's output back into a raw-pairs table and a deduplicated/sidedness
# table - see tagmaplib.NGS_QC_RAW_COLUMNS for why.
RAW_QC_COLUMNS = tagmaplib.NGS_QC_RAW_COLUMNS
DEDUPLICATED_COLUMNS = tagmaplib.NGS_QC_DEDUPLICATED_COLUMNS
SIDEDNESS_COLUMNS = tagmaplib.NGS_QC_SIDEDNESS_COLUMNS

QC_COLUMNS = RAW_QC_COLUMNS + DEDUPLICATED_COLUMNS + SIDEDNESS_COLUMNS[1:]


def sample_name_from_path(path):
    return os.path.basename(path).rsplit("_stats.yml", 1)[0]


def sample_and_side_from_pairs_path(path):
    name = os.path.basename(path)
    for side in ("forward", "reverse"):
        suffix = f"_{side}.pairs"
        if name.endswith(suffix):
            return name[: -len(suffix)], side
    raise ValueError(f"Can't tell sample/side from {path!r}")


def count_pairs(path):
    """Number of data lines in a .pairs file - pairtools dedup's own filtered
    stats don't reliably populate a "total" for named filters (unlike
    no_filter's), so forward/reverse counts are taken directly from the
    already-filtered per-side pairs files instead of the stats yml.
    """
    with open(path) as f:
        return sum(1 for line in f if not line.startswith("#"))


def sample_name_from_dedup_pairs_path(path):
    name = os.path.basename(path)
    for suffix in ("_dupmarked.pairs", "_sorted.pairs"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"Can't tell sample from {path!r}")


def count_deduplicated_junction_pairs(path):
    """Non-duplicate pairs whose junction was actually sequenced - the pool
    frac_mobilized draws mobilized_pairs from, so that a library's raw
    sequencing depth (mostly PCR duplicates and off-target pairs, not
    informative either way) doesn't drown out its real mobilization rate.

    pairtools dedup --mark-dups overwrites a duplicate's own pair_type with
    "DD" but leaves everything else, including walk_pair_type, untouched -
    so this works whether the pipeline deduped (a _dupmarked.pairs, mixing
    "DD" and kept pairs) or not (a plain _sorted.pairs, with no "DD" pairs
    at all - every row already counts as "deduplicated" the same as this
    library's own raw_junction_pairs).
    """
    pairs, _chromsizes = tagmaplib.read_pairs(path)
    if pairs.shape[0] == 0:
        return 0
    kept = pairs["pair_type"] != "DD"
    junction = pairs["walk_pair_type"].isin(tagmaplib.JUNCTION_WALK_TYPES)
    return int((kept & junction).sum())


def raw_pair_classes(
    path, primer_positions, cassette_name, cassette_length, only_junctions
):
    """Classify every pre-dedup pair using the workflow's ITR-end predicate."""
    pairs, _chromsizes = tagmaplib.read_pairs(path)
    if pairs.shape[0] == 0:
        return (
            pairs,
            pd.Series(False, index=pairs.index),
            pd.Series(np.nan, index=pairs.index, dtype=object),
        )

    anchored = (
        pairs["pair_type"].isin(["UR", "UU", "RU"])
        & (pairs["chrom1"] != pairs["chrom2"])
        & (pairs["chrom2"] == cassette_name)
    )
    forward_end = pairs["pos2"].sub(cassette_length).abs().le(2) | pairs["pos2"].sub(
        primer_positions["forward_ITR_primer_position"] + pairs["read_len2"]
    ).abs().le(2)
    reverse_end = pairs["pos2"].sub(
        primer_positions["reverse_ITR_primer_position"] - pairs["read_len2"]
    ).clip(lower=0).abs().le(2) | pairs["pos2"].le(2)
    # The pairtools expression uses the forward/reverse predicate separately;
    # retain both here so the side can be assigned for peak/site support.
    forward = anchored & forward_end
    reverse = anchored & reverse_end
    if only_junctions:
        junction = pairs["walk_pair_type"].isin(tagmaplib.JUNCTION_WALK_TYPES)
        forward &= junction
        reverse &= junction
    side = pd.Series(np.nan, index=pairs.index, dtype=object)
    side.loc[forward] = "forward"
    side.loc[reverse] = "reverse"
    anchored = side.notna()
    return pairs, anchored, side


def feature_support(pairs, anchored, side, features, sample_name):
    """Count raw anchored pairs whose junction coordinate is in a final call."""
    if features is None or features.shape[0] == 0:
        return 0
    features = features[features["sample_name"] == sample_name]
    if features.shape[0] == 0:
        return 0
    points = pairs.loc[anchored, ["chrom1", "pos51"]].copy()
    points["side"] = side.loc[anchored].to_numpy()
    points["pos51"] = tagmaplib.normalise_junction_pos(
        points["pos51"], pairs.loc[anchored, "strand1"]
    )
    supported = 0
    for (chrom, point_side), point_group in points.groupby(["chrom1", "side"]):
        candidates = features[features["chrom"] == chrom]
        if "side" in features:
            candidates = candidates[
                candidates["side"] == ("+" if point_side == "forward" else "-")
            ]
        if candidates.shape[0] == 0:
            continue
        candidates = candidates.sort_values("start")
        starts = candidates["start"].to_numpy()
        ends = candidates["end"].to_numpy()
        for position in point_group["pos51"]:
            index = np.searchsorted(starts, position, side="right") - 1
            if index >= 0 and position < ends[index]:
                supported += 1
    return supported


pair_counts = {}
for path in args.pairs:
    sample_name, side = sample_and_side_from_pairs_path(path)
    pair_counts[(sample_name, side)] = count_pairs(path)

raw_pair_paths = {
    os.path.basename(path).rsplit("_sorted.pairs", 1)[0]: path
    for path in args.raw_pairs
}
dedup_pair_paths = {
    sample_name_from_dedup_pairs_path(path): path for path in args.dedup_pairs
}
with open(args.primer_positions) as f:
    primer_positions = json.load(f)
peaks = tagmaplib.read_peaks(args.peaks)
sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str})

rows = []
for path in args.stats_yml:
    with open(path) as f:
        stat = yaml.safe_load(f) or {}
    sample_name = sample_name_from_path(path)
    total = stat.get("no_filter", {}).get("total", 0)
    mapped = stat.get("no_filter", {}).get("total_mapped", 0)
    forward = pair_counts.get((sample_name, "forward"), 0)
    reverse = pair_counts.get((sample_name, "reverse"), 0)
    mobilized = forward + reverse
    dedup_path = dedup_pair_paths.get(sample_name)
    total_dedup_junction = (
        count_deduplicated_junction_pairs(dedup_path) if dedup_path else 0
    )
    raw_path = raw_pair_paths.get(sample_name)
    if raw_path:
        raw_pairs, anchored, raw_side = raw_pair_classes(
            raw_path,
            primer_positions,
            args.cassette_name,
            args.cassette_length,
            args.only_junctions,
        )
        raw_junction = raw_pairs["walk_pair_type"].isin(tagmaplib.JUNCTION_WALK_TYPES)
        raw_anchored = int(anchored.sum())
        raw_forward_anchored = int((raw_side == "forward").sum())
        raw_reverse_anchored = int((raw_side == "reverse").sum())
        raw_peak_support = feature_support(
            raw_pairs, anchored, raw_side, peaks, sample_name
        )
        raw_site_support = feature_support(
            raw_pairs, anchored, raw_side, sites, sample_name
        )
        raw_junction_count = int((anchored & raw_junction).sum())
        raw_forward_junction_count = int(
            ((raw_side == "forward") & raw_junction).sum()
        )
        raw_reverse_junction_count = int(
            ((raw_side == "reverse") & raw_junction).sum()
        )
        raw_non_junction_count = raw_anchored - raw_junction_count
        raw_unusable = int(raw_pairs.shape[0] - raw_anchored)
    else:
        raw_anchored = raw_junction_count = raw_non_junction_count = 0
        raw_forward_anchored = raw_reverse_anchored = 0
        raw_forward_junction_count = raw_reverse_junction_count = 0
        raw_peak_support = raw_site_support = raw_unusable = 0
    rows.append(
        {
            "sample_name": sample_name,
            "total_pairs": total,
            "mapped_pairs": mapped,
            "frac_mapped": mapped / total if total else float("nan"),
            "forward_junction_pairs": forward,
            "reverse_junction_pairs": reverse,
            "mobilized_pairs": mobilized,
            "total_deduplicated_junction_pairs": total_dedup_junction,
            "frac_mobilized": (
                mobilized / total_dedup_junction if total_dedup_junction else float("nan")
            ),
            "raw_itr_anchored_pairs": raw_anchored,
            "raw_forward_itr_anchored_pairs": raw_forward_anchored,
            "raw_reverse_itr_anchored_pairs": raw_reverse_anchored,
            "raw_junction_pairs": raw_junction_count,
            "raw_forward_junction_pairs": raw_forward_junction_count,
            "raw_reverse_junction_pairs": raw_reverse_junction_count,
            "raw_non_junction_pairs": raw_non_junction_count,
            "raw_peak_support_pairs": raw_peak_support,
            "raw_site_support_pairs": raw_site_support,
            "raw_unusable_pairs": raw_unusable,
        }
    )

mapping = pd.DataFrame(rows, columns=RAW_QC_COLUMNS + DEDUPLICATED_COLUMNS).sort_values(
    "sample_name"
)

int_columns = [
    "n_sites",
    "n_two_sided",
    "n_one_sided",
    "n_forward_only",
    "n_reverse_only",
]

if sites.shape[0] == 0:
    sidedness = pd.DataFrame(columns=SIDEDNESS_COLUMNS)
else:

    def summarize(group):
        n = group.shape[0]
        forward_only = (group["site_sides"] == "forward_only").sum()
        reverse_only = (group["site_sides"] == "reverse_only").sum()
        two_sided = (group["site_sides"] == "both").sum()
        return pd.Series(
            {
                "n_sites": n,
                "n_two_sided": two_sided,
                "frac_two_sided": two_sided / n if n else float("nan"),
                "n_one_sided": forward_only + reverse_only,
                "n_forward_only": forward_only,
                "n_reverse_only": reverse_only,
            }
        )

    sidedness = (
        sites.groupby("sample_name")
        .apply(summarize, include_groups=False)
        .reset_index()[SIDEDNESS_COLUMNS]
    )
    # groupby.apply returns one Series per group; mixing ints and the float
    # frac_two_sided in that Series forces the whole thing to float64.
    sidedness[int_columns] = sidedness[int_columns].astype(int)

# A library with no insertion sites has no row in sidedness at all, rather
# than a row of zeros - fill that in instead of leaving it blank.
qc = mapping.merge(sidedness, on="sample_name", how="left")
qc[int_columns] = qc[int_columns].fillna(0).astype(int)
qc = qc[QC_COLUMNS]
qc.to_csv(args.output, sep="\t", index=False)

print(
    f"{qc.shape[0]} NGS samples: "
    f"{int(qc['mobilized_pairs'].sum()) if qc.shape[0] else 0} mobilized pairs; "
    f"{int(qc['n_sites'].sum()) if qc.shape[0] else 0} insertion sites, "
    f"{int(qc['n_two_sided'].sum()) if qc.shape[0] else 0} seen from both sides"
)
