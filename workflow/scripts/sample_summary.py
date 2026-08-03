"""Per-sample summary of one NGS run.

For each sample: how many sites were called, how many of those came from a
sequenced junction versus an ITR-anchored guess alone, how often the site
landed on the transposon's motif, and how much of the library is still on the
founder/construct contig rather than mobilised into the genome.
"""

import argparse
import os

import bioframe
import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--sites", required=True, help="all_sites.bed")
argparser.add_argument("--peaks", required=True, help="all_peaks.bed")
argparser.add_argument(
    "--evidence",
    nargs="*",
    default=[],
    help="{sample}_{side}_evidence.tsv files. Only the junction_tiered caller "
    "writes these; without them there is no junction/anchor split to report.",
)
argparser.add_argument(
    "--coverage",
    nargs="+",
    required=True,
    help="{sample}_{side}_coverage.bedgraph files",
)
argparser.add_argument(
    "--construct-contigs",
    nargs="*",
    default=[],
    help="Contigs that are founder/vector sequence, not genome",
)
argparser.add_argument(
    "--max-dist",
    type=int,
    required=True,
    help="Same value find_insertion_sites.py clustered peaks with, so peaks "
    "here recluster into the same sites it already reported",
)
argparser.add_argument("--output", "-o", required=True)


def sample_and_side(path, suffix):
    base = os.path.basename(path)[: -len(suffix)]
    sample, side = base.rsplit("_", 1)
    return sample, side


def read_coverage(paths):
    """Per-position coverage, plus the full set of configured sample names -
    every sample gets a coverage file, even an empty one, so this is also
    where the complete sample list for the summary comes from."""
    frames = []
    samples = set()
    for path in paths:
        sample, _side = sample_and_side(path, "_coverage.bedgraph")
        samples.add(sample)
        try:
            df = pd.read_csv(
                path,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "count", "fraction", "orientation"],
                dtype={"chrom": str},
            )
        except pd.errors.EmptyDataError:
            continue
        if df.shape[0] == 0:
            continue
        df["sample_name"] = sample
        frames.append(df)
    coverage = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["chrom", "count", "sample_name"]
    )
    return coverage, samples


def mobilization_efficiency(coverage, construct_contigs):
    """Fraction of ITR-anchored read coverage that is not on a construct
    contig, i.e. actually mobilised into the genome rather than still sitting
    at the founder/donor site. Counts raw pairs, not the distinct
    tagmentation positions peaks are counted on - the question here is what
    fraction of the library moved, not how many independent insertions it
    represents.
    """
    by_sample = coverage.groupby("sample_name")["count"].sum()
    on_construct = (
        coverage[coverage["chrom"].isin(construct_contigs)]
        .groupby("sample_name")["count"]
        .sum()
    )
    return (
        1 - on_construct.reindex(by_sample.index, fill_value=0) / by_sample
    ).rename("mobilization_efficiency")


def tier_counts(peaks, evidence_paths, max_dist):
    """Junction-tier vs anchor-tier site counts per sample.

    Peaks are re-clustered with the same --max-dist find_insertion_sites.py
    used, so each cluster here is the same site it reported; a cluster with
    any junction-tier peak behind it is a junction-tier site. Callers other
    than junction_tiered never produce evidence files, so they correctly skip
    this - they have no per-peak tier to report.
    """
    if not evidence_paths or peaks.shape[0] == 0:
        return None

    evidence_frames = []
    for path in evidence_paths:
        sample, side = sample_and_side(path, "_evidence.tsv")
        try:
            ev = pd.read_csv(path, sep="\t", dtype={"chrom": str})
        except pd.errors.EmptyDataError:
            continue
        if ev.shape[0] == 0:
            continue
        ev = ev[["chrom", "start", "end", "tier"]].copy()
        ev["sample_name"] = sample
        # combine_peaks.py encodes forward/reverse as +/- in all_peaks.bed's
        # "side" column, not the "forward"/"reverse" the filenames use.
        ev["side"] = "+" if side == "forward" else "-"
        evidence_frames.append(ev)
    if not evidence_frames:
        return None
    evidence = pd.concat(evidence_frames, ignore_index=True)

    tagged = peaks.merge(
        evidence, on=["chrom", "start", "end", "sample_name", "side"], how="left"
    )
    tagged = bioframe.cluster(
        tagged, min_dist=max_dist, on=["sample_name"], return_cluster_ids=True
    )
    cluster_tier = tagged.groupby(["sample_name", "cluster"])["tier"].apply(
        lambda s: "junction" if (s == "junction").any() else "anchor"
    )
    return (
        cluster_tier.reset_index()
        .groupby(["sample_name", "tier"])
        .size()
        .unstack("tier", fill_value=0)
        .reindex(columns=["junction", "anchor"], fill_value=0)
        .rename(columns={"junction": "n_junction_tier", "anchor": "n_anchor_tier"})
    )


if __name__ == "__main__":
    args = argparser.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str})
    coverage, cov_samples = read_coverage(args.coverage)
    peaks = tagmaplib.read_peaks(args.peaks).astype(
        {"chrom": str, "start": int, "end": int}
    )

    samples = sorted(cov_samples | set(sites.get("sample_name", pd.Series(dtype=str))))
    summary = pd.DataFrame(index=pd.Index(samples, name="sample_name"))

    summary["n_sites"] = (
        sites.groupby("sample_name").size().reindex(summary.index, fill_value=0)
    )

    tiers = tier_counts(peaks, args.evidence, args.max_dist)
    if tiers is not None:
        summary = summary.join(tiers.reindex(summary.index, fill_value=0))

    ta_columns = [c for c in sites.columns if c.endswith("_found")]
    if ta_columns:
        summary["ta_rate"] = (
            sites.groupby("sample_name")[ta_columns[0]]
            .mean()
            .reindex(summary.index)
            .round(4)
        )

    summary = summary.join(
        mobilization_efficiency(coverage, args.construct_contigs).round(4)
    )

    summary.reset_index().sort_values("sample_name").to_csv(
        args.output, sep="\t", index=False
    )
