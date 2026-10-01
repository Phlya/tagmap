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
argparser.add_argument("--require-both-sides", action="store_true")
argparser.add_argument(
    "--merge-dist",
    type=int,
    default=0,
    help="A one-sided site within this distance of one seen only from the "
    "opposite ITR is folded into it, as the same insertion whose two sides "
    "were positioned a few bases apart. 0 (default) never merges.",
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


def merge_complementary(sites, distance):
    """Fold one-sided sites into nearby ones seen from the opposite side."""
    if distance <= 0 or sites.shape[0] < 2:
        return sites
    sites = sites.sort_values("n_molecules", ascending=False).reset_index(drop=True)
    kept = []
    for site in sites.to_dict("records"):
        target = None
        if site["site_sides"] in ONE_SIDED:
            for other in kept:
                if (
                    other["chrom"] == site["chrom"]
                    and other["site_sides"] in ONE_SIDED
                    and other["site_sides"] != site["site_sides"]
                    and abs(other["start"] - site["start"]) <= distance
                ):
                    target = other
                    break
        if target is None:
            kept.append(site)
        else:
            for column in ("n_forward", "n_reverse", "n_molecules"):
                target[column] += site[column]
            target["site_sides"] = "both"
    return pd.DataFrame(kept, columns=sites.columns)


def near(chrom, start, others, distance):
    """Rows of `others` on `chrom` within `distance` of `start`."""
    return others[(others["chrom"] == chrom) & ((others["start"] - start).abs() <= distance)]


def position(row):
    return f"{row['chrom']}:{row['start']}"


def judge(sites, n_pairs, original, args):
    """The summary row (without sample/plate/well/shared_with) for one clone."""
    row = {
        "n_pairs": n_pairs,
        "n_sites": sites.shape[0],
        "n_molecules": int(sites["n_molecules"].sum()),
        "n_secondary_sites": 0,
        "at_original_site": False,
        "reason": "",
    }
    if n_pairs < args.min_reads:
        row["status"] = tagmaplib.CLONE_STATUS_TOO_FEW_READS
        row["reason"] = f"{n_pairs} mobilized pairs, need {args.min_reads}"
        return row, None
    if sites.shape[0] == 0:
        row["status"] = tagmaplib.CLONE_STATUS_NO_INSERTION
        row["reason"] = "no insertion site called"
        return row, None

    sites = sites.sort_values(["n_molecules", "start"], ascending=[False, True])
    total = sites["n_molecules"].sum()
    top = sites.iloc[0]
    row.update(
        chrom=top["chrom"],
        start=int(top["start"]),
        end=int(top["end"]),
        strand=top["strand"],
        sides=top["site_sides"],
        dominant_molecules=int(top["n_molecules"]),
        dominant_frac=top["n_molecules"] / total,
    )
    others = sites.iloc[1:]
    if others.shape[0]:
        second = others.iloc[0]
        row.update(
            second_position=position(second),
            second_molecules=int(second["n_molecules"]),
            second_frac=second["n_molecules"] / total,
        )
        row["n_secondary_sites"] = int(
            (others["n_molecules"] / total >= args.max_contamination_frac).sum()
        )
    else:
        row.update(second_position="", second_molecules=0, second_frac=0.0)

    if original is not None:
        row["at_original_site"] = bool(
            top["chrom"] == original["chrom"]
            and abs(top["start"] - original["pos"]) <= args.original_max_dist
        )

    if top["n_molecules"] < args.min_molecules:
        status = tagmaplib.CLONE_STATUS_WEAK
        reason = f"dominant site has {int(top['n_molecules'])} molecules, need {args.min_molecules}"
    elif row["second_frac"] >= args.multi_site_frac:
        status = tagmaplib.CLONE_STATUS_MULTIPLE
        reason = f"second site {row['second_position']} holds {row['second_frac']:.0%}"
    elif row["second_frac"] >= args.max_contamination_frac:
        status = tagmaplib.CLONE_STATUS_CONTAMINATED
        reason = f"second site {row['second_position']} holds {row['second_frac']:.0%}"
    elif row["at_original_site"]:
        status = tagmaplib.CLONE_STATUS_UNMOBILIZED
        reason = "dominant site is the original insertion site"
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
        sites = merge_complementary(
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

    summary = summary[tagmaplib.CLONE_SUMMARY_COLUMNS].sort_values(
        ["plate", "well"], kind="stable"
    )
    # Nullable integers, so that a clone with no site (NaN) does not turn
    # every position in the column into a float.
    for column in ("n_pairs", "n_sites", "n_molecules", "start", "end",
                   "dominant_molecules", "second_molecules", "n_secondary_sites"):
        summary[column] = summary[column].astype("Int64")
    summary.to_csv(args.output, sep="\t", index=False)

    counts = summary["status"].value_counts().to_dict()
    print("Clone calls:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
