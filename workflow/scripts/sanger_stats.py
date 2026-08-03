"""Summarise Sanger sequencing QC and per-clone coverage.

Reads that pass QC already end up in the per-clone consensus sites
(combine_sanger_sites.py), but a read that failed, or a primer direction that
was never sequenced for a clone, simply has no site to show for it - so
"missing" and "failed" look the same downstream. This works from the raw
per-read tables instead, which keep both, and reports for every clone whether
each ITR primer's direction passed, and if not, why: never sequenced, or the
QC failure reason itself. When NGS data for the same material is available,
also reports whether that clone's site was independently confirmed there
too.
"""

import argparse
import json

import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument(
    "--reads", nargs="*", default=[], help="{sample}_reads.tsv files"
)
argparser.add_argument(
    "--sites",
    default=None,
    help="all_sanger_sites.bed, to check whether a clone's forward and "
    "reverse reads agree on position rather than just each passing QC "
    "independently",
)
argparser.add_argument(
    "--validation", default=None, help="sanger_vs_ngs.tsv, if NGS data is available"
)
argparser.add_argument(
    "--original-site", default=None, help="original_site.json, if configured"
)
argparser.add_argument(
    "--original-max-dist",
    type=int,
    default=25,
    help="A clustered site within this many bp of --original-site is unmobilized",
)
argparser.add_argument("--output-qc", required=True)
argparser.add_argument("--output-clone-summary", required=True)
args = argparser.parse_args()

QC_COLUMNS = [
    "sample_name",
    "direction",
    "n_reads",
    "n_pass",
    "n_fail",
    "frac_pass",
    "fail_reasons",
]

BASE_CLONE_COLUMNS = [
    "sample_name",
    "clone",
    "n_forward_reads",
    "n_reverse_reads",
    "forward_status",
    "reverse_status",
]

ORIGINAL_SITE_COLUMN = "unmobilized"
CLONE_VALIDATION_COLUMNS = ["ngs_validated", "ngs_side"]

original_site = None
if args.original_site is not None:
    with open(args.original_site) as f:
        original_site = json.load(f)

CLONE_COLUMNS = (
    BASE_CLONE_COLUMNS
    + ["position", "positions_agree", "both_sides_confirmed"]
    + ([ORIGINAL_SITE_COLUMN] if original_site is not None else [])
    + ["summary"]
)


def at_original_site(sites, original_site, max_dist):
    return (sites["chrom"] == original_site["chrom"]) & (
        (sites["start"] - original_site["pos"]).abs() <= max_dist
    )


def classify(direction_reads):
    """no data / failed / confirmed for one clone's reads from one primer."""
    if direction_reads.shape[0] == 0:
        return "no data"
    if direction_reads["pass"].astype(bool).any():
        return "confirmed"
    return "failed"


def describe(forward_status, reverse_status, positions_agree):
    """A one-line summary of a clone's forward/reverse coverage.

    forward_status/reverse_status are "PASSED", "not sequenced", or the QC
    failure reason itself - so a passing side is identified by that exact
    string, not by absence of a separate reason column.

    Both sides can independently pass QC while pointing at different loci -
    e.g. one read stuck on a residual copy of the donor construct - so
    "PASSED" from both primers isn't enough; their sites also have to
    cluster together (positions_agree, from the already-clustered
    all_sanger_sites.bed) for the clone to really be confirmed.
    """
    forward_ok = forward_status == "PASSED"
    reverse_ok = reverse_status == "PASSED"
    if forward_ok and reverse_ok:
        if not positions_agree:
            return "forward and reverse confirmed but disagree on position"
        return "both sides confirmed"
    if forward_ok or reverse_ok:
        if forward_ok:
            confirmed_side, other_side, other_status = (
                "forward",
                "reverse",
                reverse_status,
            )
        else:
            confirmed_side, other_side, other_status = (
                "reverse",
                "forward",
                forward_status,
            )
        return f"{confirmed_side} only confirmed ({other_side}: {other_status})"
    if forward_status == "not sequenced" and reverse_status == "not sequenced":
        return "no Sanger data"
    return f"neither confirmed (forward: {forward_status}, reverse: {reverse_status})"


def fail_reason_summary(group):
    """'reason (n); reason (n)', most common first, for one failed group."""
    counts = group["reason"].value_counts()
    return "; ".join(f"{reason} ({n})" for reason, n in counts.items())


def clone_fail_reasons(direction_reads):
    """'reason (n); reason (n)' for the failed reads of one clone's primer
    direction, or "" if none failed."""
    failed = direction_reads[~direction_reads["pass"].astype(bool)].copy()
    if failed.empty:
        return ""
    failed["fail_reason"] = failed["fail_reason"].fillna("unknown").astype(str)
    failed.loc[failed["fail_reason"] == "", "fail_reason"] = "unknown"
    reasons = failed["fail_reason"].str.split("; ").explode()
    return fail_reason_summary(pd.DataFrame({"reason": reasons}))


def side_status(status, direction_reads):
    """The status of one clone's reads from one primer: "PASSED", "not
    sequenced", or the QC failure reason itself."""
    if status == "confirmed":
        return "PASSED"
    if status == "no data":
        return "not sequenced"
    return clone_fail_reasons(direction_reads)


def format_site(row):
    """One site cluster as 'chrom:start-end(strand, n reads)'."""
    return f"{row.chrom}:{row.start}-{row.end}({row.strand}, {row.n_reads} reads)"


def summarize_clone_sites(group, original_site=None, max_dist=0):
    """A clone's site(s), whether they agree on position, and - when the
    original, pre-mobilization locus is configured - whether any of them
    sits there rather than at a new locus.

    Usually one site backed by both directions, but a clone whose forward
    and reverse reads don't cluster together (see combine_sanger_sites.py)
    ends up with more than one single-direction site here - listing all of
    them shows where the disagreement actually is, rather than just that
    there is one.
    """
    sites = group.sort_values(["chrom", "start"])
    result = {
        "position": "; ".join(format_site(row) for row in sites.itertuples()),
        "positions_agree": bool(group["both_directions"].any()),
    }
    if original_site is not None:
        result[ORIGINAL_SITE_COLUMN] = bool(
            at_original_site(sites, original_site, max_dist).any()
        )
    return pd.Series(result)


def summarize_clone_validation(group):
    """Whether a clone's Sanger site was independently confirmed by NGS, and
    which side(s) of the NGS insertion site did the confirming."""
    confirmed = group[group["confirmed_by_ngs"]]
    sides = sorted(confirmed["ngs_site_sides"].dropna().unique())
    return pd.Series(
        {
            "ngs_validated": bool(confirmed.shape[0]),
            "ngs_side": ", ".join(sides) if sides else pd.NA,
        }
    )


if not args.reads:
    qc = pd.DataFrame(columns=QC_COLUMNS)
    clone_summary = pd.DataFrame(columns=BASE_CLONE_COLUMNS)
else:
    reads = pd.concat(
        (pd.read_csv(path, sep="\t", dtype={"chrom": str}) for path in args.reads),
        ignore_index=True,
    )
    reads["pass"] = tagmaplib.to_bool(reads["pass"])

    if reads.shape[0] == 0:
        qc = pd.DataFrame(columns=QC_COLUMNS)
        clone_summary = pd.DataFrame(columns=BASE_CLONE_COLUMNS)
    else:

        def summarize_reads(group):
            n = group.shape[0]
            n_pass = int(group["pass"].sum())
            return pd.Series(
                {
                    "n_reads": n,
                    "n_pass": n_pass,
                    "n_fail": n - n_pass,
                    "frac_pass": n_pass / n if n else float("nan"),
                }
            )

        qc = (
            reads.groupby(["sample_name", "direction"])
            .apply(summarize_reads, include_groups=False)
            .reset_index()
        )
        # groupby.apply returns one Series per group; mixing ints and the
        # float frac_pass in that Series forces the whole thing to float64.
        qc[["n_reads", "n_pass", "n_fail"]] = qc[
            ["n_reads", "n_pass", "n_fail"]
        ].astype(int)

        failed = reads[~reads["pass"]].copy()
        failed["fail_reason"] = failed["fail_reason"].fillna("unknown").astype(str)
        failed.loc[failed["fail_reason"] == "", "fail_reason"] = "unknown"
        exploded = failed.assign(reason=failed["fail_reason"].str.split("; ")).explode(
            "reason"
        )
        if exploded.shape[0]:
            fail_reasons = (
                exploded.groupby(["sample_name", "direction"])
                .apply(fail_reason_summary, include_groups=False)
                .reset_index(name="fail_reasons")
            )
        else:
            fail_reasons = pd.DataFrame(
                columns=["sample_name", "direction", "fail_reasons"]
            )

        qc = qc.merge(fail_reasons, on=["sample_name", "direction"], how="left")
        qc["fail_reasons"] = qc["fail_reasons"].fillna("")
        qc = qc[QC_COLUMNS].sort_values(["sample_name", "direction"])

        def summarize_clone(group):
            forward = group[group["direction"] == "forward"]
            reverse = group[group["direction"] == "reverse"]
            return pd.Series(
                {
                    "n_forward_reads": forward.shape[0],
                    "n_reverse_reads": reverse.shape[0],
                    "forward_status": side_status(classify(forward), forward),
                    "reverse_status": side_status(classify(reverse), reverse),
                }
            )

        clone_summary = (
            reads.groupby(["sample_name", "clone"])
            .apply(summarize_clone, include_groups=False)
            .reset_index()[BASE_CLONE_COLUMNS]
            .sort_values(["sample_name", "clone"])
        )
        clone_summary[["n_forward_reads", "n_reverse_reads"]] = clone_summary[
            ["n_forward_reads", "n_reverse_reads"]
        ].astype(int)

# Both primers can independently pass QC while pointing at different loci, so
# "confirmed" from each side isn't enough on its own - cross-check against the
# already-clustered sites, which group a clone's passing reads by position.
if args.sites is not None:
    site_clusters = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str})
    if site_clusters.shape[0] and "both_directions" in site_clusters.columns:
        site_clusters["both_directions"] = tagmaplib.to_bool(
            site_clusters["both_directions"]
        )
        site_summary = (
            site_clusters.groupby(["sample_name", "clone"])
            .apply(
                summarize_clone_sites,
                original_site=original_site,
                max_dist=args.original_max_dist,
                include_groups=False,
            )
            .reset_index()
        )
    else:
        site_summary = pd.DataFrame(
            columns=["sample_name", "clone", "position", "positions_agree"]
            + ([ORIGINAL_SITE_COLUMN] if original_site is not None else [])
        )
    clone_summary = clone_summary.merge(
        site_summary, on=["sample_name", "clone"], how="left"
    )
    clone_summary["position"] = clone_summary["position"].fillna("")
    clone_summary["positions_agree"] = (
        clone_summary["positions_agree"].fillna(False).astype(bool)
    )
    if original_site is not None:
        clone_summary[ORIGINAL_SITE_COLUMN] = (
            clone_summary[ORIGINAL_SITE_COLUMN].fillna(False).astype(bool)
        )
else:
    clone_summary["position"] = ""
    clone_summary["positions_agree"] = False
    if original_site is not None:
        clone_summary[ORIGINAL_SITE_COLUMN] = False

clone_summary["both_sides_confirmed"] = (
    (clone_summary["forward_status"] == "PASSED")
    & (clone_summary["reverse_status"] == "PASSED")
    & clone_summary["positions_agree"]
)
clone_summary["summary"] = [
    describe(f, r, p)
    for f, r, p in zip(
        clone_summary["forward_status"],
        clone_summary["reverse_status"],
        clone_summary["positions_agree"],
    )
]
clone_summary = clone_summary[CLONE_COLUMNS]

if args.validation is not None:
    validation = pd.read_csv(args.validation, sep="\t", dtype={"chrom": str})
    if validation.shape[0] and "confirmed_by_ngs" in validation.columns:
        validation["confirmed_by_ngs"] = tagmaplib.to_bool(
            validation["confirmed_by_ngs"]
        )
        clone_validation = (
            validation.groupby(["sample_name", "clone"])
            .apply(summarize_clone_validation, include_groups=False)
            .reset_index()
        )
    else:
        clone_validation = pd.DataFrame(
            columns=["sample_name", "clone"] + CLONE_VALIDATION_COLUMNS
        )
    clone_summary = clone_summary.merge(
        clone_validation, on=["sample_name", "clone"], how="left"
    )
    clone_summary["ngs_validated"] = (
        clone_summary["ngs_validated"].fillna(False).astype(bool)
    )
    clone_summary = clone_summary[CLONE_COLUMNS + CLONE_VALIDATION_COLUMNS]

qc.to_csv(args.output_qc, sep="\t", index=False)
clone_summary.to_csv(args.output_clone_summary, sep="\t", index=False)

print(
    f"{clone_summary.shape[0]} clones: "
    f"{int(clone_summary['both_sides_confirmed'].sum()) if clone_summary.shape[0] else 0} "
    "confirmed from both primers"
)
