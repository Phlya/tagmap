"""Judge each NGS library as a clone: one dominant insertion, nothing else.

A clonal line carries a single insertion, so its library should hold one site
far above any other. The site calling is the pool pipeline's, unchanged -
sites are positioned to the base from the sequenced junction, and abundance is
counted in independent molecules (distinct tagmentation positions) rather than
reads, so PCR copies cannot inflate the dominant site or hide a contaminant.
What this adds is the verdict per library:

  clean          one site, from enough molecules, and nothing else at or above
                 clone_max_contamination_frac
  contaminated   one dominant site, but another holds that much or more (index
                 hopping, barcode cross-talk, carry-over from a neighbouring
                 clone - see shared_with)
  multiple       another site holds clone_multi_site_frac or more: the clone
                 really carries several insertions, or is a mix of clones
  unmobilized    clean, but the site is the original, pre-mobilization locus
  weak           too few molecules behind the dominant site (or, if asked for,
                 it was only seen from one ITR side)
  no_insertion   reads, but no site was called
  too_few_reads  not enough mobilized read pairs to judge

shared_with names the other clones whose dominant site shows up as a secondary
site here, which is what tells cross-talk from a second real insertion.
"""

import argparse
import json

import numpy as np
import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--support", required=True, help="all_sites_support.tsv")
argparser.add_argument("--ngs-qc", required=True, help="ngs_qc_stats.tsv")
argparser.add_argument(
    "--samples",
    nargs="+",
    required=True,
    help="sample=library:plate for every NGS sample (plate may be empty)",
)
argparser.add_argument("--original-site", default=None, help="original_site.json")
argparser.add_argument("--original-max-dist", type=int, default=0)
argparser.add_argument("--min-reads", type=int, default=20)
argparser.add_argument("--min-molecules", type=int, default=3)
argparser.add_argument("--max-contamination-frac", type=float, default=0.02)
argparser.add_argument("--multi-site-frac", type=float, default=0.2)
argparser.add_argument(
    "--min-dominant-frac",
    type=float,
    default=0.0,
    help="The dominant site must hold at least this share of the clone's "
    "molecules. Unlike the contamination tests, this counts every molecule, "
    "including the single-molecule scatter no individual site is judged on - "
    "so it catches a library that is mostly noise even though no one "
    "secondary site stands out. 0 disables it.",
)
argparser.add_argument("--require-both-sides", action="store_true")
argparser.add_argument(
    "--contaminant-min-molecules",
    type=int,
    default=2,
    help="A secondary site needs at least this many distinct molecules before "
    "it counts against the clone at all. At a few dozen molecules per clone a "
    "single stray one already clears a low percentage, so the fraction alone "
    "cannot separate contamination from mismapping.",
)
argparser.add_argument(
    "--contaminant-require-both-sides",
    action="store_true",
    help="A secondary site also has to have been seen from both ITR primers. "
    "Strict: one-sided calls are overwhelmingly noise, but a genuine second "
    "insertion whose other side was not captured is then missed too.",
)
argparser.add_argument(
    "--merge-dist",
    type=int,
    default=0,
    help="Sites within this distance of a better-supported one are folded "
    "into it, as one insertion reported more than once rather than "
    "neighbouring integrations. 0 (default) never merges.",
)
argparser.add_argument(
    "--shared-dist",
    type=int,
    default=2,
    help="A secondary site this close to another clone's dominant site is "
    "that site",
)
argparser.add_argument("--output", "-o", required=True)

ONE_SIDED = ("forward_only", "reverse_only")
SIDE_COLUMNS = {"n_forward": "forward_only", "n_reverse": "reverse_only"}


def merge_nearby(sites, distance):
    """Fold sites within `distance` of a better-supported one into it.

    One clone carries one insertion, so several called sites within a few tens
    of bases of each other are one insertion reported more than once, not
    neighbouring integrations: the two ITR sides can be positioned a few bases
    apart (the strand pos5 offset, the target-site duplication, indel wobble),
    and reads that stopped short of the junction scatter further still. Left
    unmerged they show up as a "second site" holding a sizeable share of the
    molecules, which would read as contamination.

    Sites are taken strongest first, so the position and strand of the merged
    site are the best-supported one's, and the molecule counts are pooled -
    which can turn two one-sided calls into a single two-sided one.

    The cost of setting this too wide is in the other direction: in a local
    re-mobilization assay two genuinely independent insertions can sit under
    100bp apart, and merging would hide one of them. Keep it well below the
    spacing of insertions the experiment is meant to resolve.
    """
    if distance <= 0 or sites.shape[0] < 2:
        return sites
    sites = sites.sort_values(
        ["n_molecules", "start"], ascending=[False, True]
    ).reset_index(drop=True)
    kept = []
    for site in sites.to_dict("records"):
        target = None
        for other in kept:
            if (
                other["chrom"] == site["chrom"]
                and abs(other["start"] - site["start"]) <= distance
            ):
                target = other
                break
        if target is None:
            kept.append(site)
            continue
        for column in ("n_forward", "n_reverse", "n_molecules"):
            target[column] += site[column]
        # Pooling can complete a side the stronger call never saw.
        seen = {
            label for column, label in SIDE_COLUMNS.items() if target[column] > 0
        }
        target["site_sides"] = seen.pop() if len(seen) == 1 else "both"
    return pd.DataFrame(kept, columns=sites.columns)


def near(chrom, start, others, distance):
    """Rows of `others` on `chrom` within `distance` of `start`."""
    return others[(others["chrom"] == chrom) & ((others["start"] - start).abs() <= distance)]


def position(row):
    return f"{row['chrom']}:{row['start']}"


def dominant_owners(row, dominant, distance):
    """Clones whose own dominant site sits at this clone's second_position."""
    where = getattr(row, "second_position", None)
    if not where or pd.isna(where):
        return []
    chrom, start = str(where).rsplit(":", 1)
    hits = near(chrom, int(start), dominant, distance)
    return sorted(set(hits.loc[hits["sample_name"] != row.sample_name, "sample_name"]))


def at_original(sites, original, max_dist):
    """Which of `sites` sit at the original, pre-mobilization locus."""
    if original is None or sites.shape[0] == 0:
        return pd.Series(False, index=sites.index)
    return (sites["chrom"] == original["chrom"]) & (
        (sites["start"] - original["pos"]).abs() <= max_dist
    )


def counts_as_contaminant(sites, total, args):
    """Which secondary sites are real enough to hold against a clone.

    A called site is not automatically a competing insertion. At this depth
    most secondary calls are a single molecule seen from one ITR only, which
    is what mismapping and chimeric PCR products look like; a fraction test
    alone cannot tell them apart, because one stray molecule out of the few
    dozen a clone yields already clears a low percentage. So a contaminant
    has to look like an insertion in its own right - enough molecules, and
    optionally both ITR sides - before its share is even considered.
    """
    real = sites["n_molecules"] >= args.contaminant_min_molecules
    if args.contaminant_require_both_sides:
        real &= sites["site_sides"] == "both"
    return real & (sites["n_molecules"] / total >= args.max_contamination_frac)


def judge(sites, n_pairs, original, args):
    """The summary row (without sample/plate/well/shared_with) for one clone."""
    row = {
        "n_pairs": n_pairs,
        "n_sites": sites.shape[0],
        "n_molecules": int(sites["n_molecules"].sum()) if sites.shape[0] else 0,
        "n_secondary_sites": 0,
        "at_original_site": False,
        "founder_molecules": 0,
        "founder_frac": 0.0,
        "reason": "",
    }
    # Too thin to judge, but still described below: what little was seen is
    # worth reading even when it decides nothing, so the row carries the same
    # columns as any other and only the status says not to trust them.
    too_few_reads = n_pairs < args.min_reads
    if sites.shape[0] == 0:
        row["status"] = (
            tagmaplib.CLONE_STATUS_TOO_FEW_READS
            if too_few_reads
            else tagmaplib.CLONE_STATUS_NO_INSERTION
        )
        row["reason"] = (
            f"{n_pairs} mobilized pairs, need {args.min_reads}"
            if too_few_reads
            else "no insertion site called"
        )
        return row, None

    sites = sites.sort_values(["n_molecules", "start"], ascending=[False, True])
    # A clone that mobilized has left the donor locus, so signal still sitting
    # there is not expected residue - it is cassette from cells that never
    # mobilized, i.e. contamination, and counts as such. It is only tracked
    # separately (founder_molecules/founder_frac) because knowing a
    # contaminant is unmobilized carry-over rather than another clone's
    # insertion says something different about where it came from.
    founder = at_original(sites, original, args.original_max_dist)
    if founder.any() and row["n_molecules"]:
        row["founder_molecules"] = int(sites.loc[founder, "n_molecules"].sum())
        row["founder_frac"] = row["founder_molecules"] / row["n_molecules"]

    top = sites.iloc[0]
    row["at_original_site"] = bool(founder.iloc[0])
    pool = sites
    total = pool["n_molecules"].sum()
    row.update(
        chrom=top["chrom"],
        start=int(top["start"]),
        end=int(top["end"]),
        strand=top["strand"],
        sides=top["site_sides"],
        dominant_molecules=int(top["n_molecules"]),
        dominant_frac=top["n_molecules"] / total,
    )

    others = pool.iloc[1:]
    second_is_founder = False
    real = counts_as_contaminant(others, total, args) if others.shape[0] else None
    row["n_secondary_sites"] = int(real.sum()) if real is not None else 0
    if real is not None and real.any():
        second = others[real].iloc[0]
        second_is_founder = bool(founder.loc[second.name])
        row.update(
            second_position=position(second),
            second_molecules=int(second["n_molecules"]),
            second_frac=second["n_molecules"] / total,
        )
    else:
        row.update(second_position="", second_molecules=0, second_frac=0.0)

    where = (
        f"second site {row['second_position']}"
        + (" (the unmobilized donor locus)" if second_is_founder else "")
        + f" holds {row['second_frac']:.0%} ({row['second_molecules']} molecules)"
    )
    if too_few_reads:
        status = tagmaplib.CLONE_STATUS_TOO_FEW_READS
        reason = (
            f"{n_pairs} mobilized pairs, need {args.min_reads} - the site "
            "below is what was seen, not a call"
        )
    elif row["at_original_site"]:
        # Mostly donor locus, but with a real integration underneath it: a
        # mixed well, not a clone that simply never mobilized.
        status = (
            tagmaplib.CLONE_STATUS_CONTAMINATED
            if row["n_secondary_sites"]
            else tagmaplib.CLONE_STATUS_UNMOBILIZED
        )
        reason = (
            f"dominant site is the original insertion site, but {where}"
            if row["n_secondary_sites"]
            else "dominant site is the original insertion site"
        )
    elif top["n_molecules"] < args.min_molecules:
        status = tagmaplib.CLONE_STATUS_WEAK
        reason = f"dominant site has {int(top['n_molecules'])} molecules, need {args.min_molecules}"
    elif row["second_frac"] >= args.multi_site_frac:
        # A heavily contaminating unmobilized population is still contamination,
        # not the clone carrying two insertions.
        status = (
            tagmaplib.CLONE_STATUS_CONTAMINATED
            if second_is_founder
            else tagmaplib.CLONE_STATUS_MULTIPLE
        )
        reason = where
    elif row["n_secondary_sites"]:
        status = tagmaplib.CLONE_STATUS_CONTAMINATED
        reason = where
    elif row["n_molecules"] and (
        top["n_molecules"] / row["n_molecules"] < args.min_dominant_frac
    ):
        status = tagmaplib.CLONE_STATUS_WEAK
        reason = (
            f"dominant site holds only "
            f"{top['n_molecules'] / row['n_molecules']:.0%} of all molecules, "
            f"scattered over {row['n_sites']} called sites"
        )
    elif args.require_both_sides and top["site_sides"] != "both":
        status = tagmaplib.CLONE_STATUS_WEAK
        reason = f"dominant site only seen as {top['site_sides']}"
    else:
        status = tagmaplib.CLONE_STATUS_CLEAN
        reason = ""
    row["status"] = status
    row["reason"] = reason
    return row, sites


if __name__ == "__main__":
    args = argparser.parse_args()

    support = pd.read_csv(args.support, sep="\t", dtype={"chrom": str})
    qc = pd.read_csv(args.ngs_qc, sep="\t").set_index("sample_name")
    pairs = qc["mobilized_pairs"] if "mobilized_pairs" in qc.columns else pd.Series(dtype=float)
    original = None
    if args.original_site:
        with open(args.original_site) as f:
            original = json.load(f)

    sample_map = {}
    for entry in args.samples:
        sample, rest = entry.split("=", 1)
        library, plate = rest.split(":", 1)
        sample_map[sample] = (library, plate)

    rows = []
    sites_by_sample = {}
    for sample, (library, plate) in sample_map.items():
        sites = merge_nearby(
            support[support["sample_name"] == sample], args.merge_dist
        )
        n_pairs = pairs.get(sample, 0)
        n_pairs = 0 if pd.isna(n_pairs) else int(n_pairs)
        row, ranked = judge(sites, n_pairs, original, args)
        row.update(sample_name=sample, plate=plate, well=library)
        rows.append(row)
        if ranked is not None:
            sites_by_sample[sample] = ranked

    summary = pd.DataFrame(rows)
    for column in tagmaplib.CLONE_SUMMARY_COLUMNS:
        if column not in summary.columns:
            summary[column] = np.nan

    # Which other clones' dominant site turns up as a secondary site here.
    dominant = summary.loc[
        summary["chrom"].notna(), ["sample_name", "chrom", "start"]
    ].astype({"start": int})
    shared = {}
    for sample, ranked in sites_by_sample.items():
        sources = set()
        for _, secondary in ranked.iloc[1:].iterrows():
            hits = near(secondary["chrom"], secondary["start"], dominant, args.shared_dist)
            sources.update(hits.loc[hits["sample_name"] != sample, "sample_name"])
        shared[sample] = sorted(sources)
    summary["shared_with"] = summary["sample_name"].map(
        lambda sample: ", ".join(shared.get(sample, [])[:5])
        + (" ..." if len(shared.get(sample, [])) > 5 else "")
    )

    # Where the contamination this clone was actually called on came from: the
    # clone whose own insertion sits at the site reported in second_position.
    # Narrower than shared_with on purpose - that one lists every secondary
    # site traceable to another clone, however thin, including the ones that
    # never cleared the thresholds, so it still hints at carry-over in a clone
    # called clean. This column names only the source of the contamination the
    # status rests on, and is empty when nothing qualified.
    summary["contaminated_by"] = [
        ", ".join(dominant_owners(row, dominant, args.shared_dist))
        for row in summary.itertuples()
    ]

    # In the order --samples listed them, i.e. the sample sheet's own order, so
    # the table lines up with how the plate was laid out rather than with how
    # the names happen to sort.
    summary = summary.set_index("sample_name").loc[list(sample_map)].reset_index()
    summary = summary[tagmaplib.CLONE_SUMMARY_COLUMNS]
    # Nullable integers, so that a clone with no site (NaN) does not turn
    # every position in the column into a float.
    for column in ("n_pairs", "n_sites", "n_molecules", "start", "end",
                   "dominant_molecules", "second_molecules", "n_secondary_sites"):
        summary[column] = summary[column].astype("Int64")
    summary.to_csv(args.output, sep="\t", index=False)

    counts = summary["status"].value_counts().to_dict()
    print("Clone calls:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
