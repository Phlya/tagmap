"""Summarise Sanger sequencing QC and per-clone coverage.

Reads that pass QC already end up in the per-clone consensus sites
(combine_sanger_sites.py), but a read that failed, or a primer direction that
was never sequenced for a clone, simply has no site to show for it - so
"missing" and "failed" look the same downstream. This works from the raw
per-read tables instead, which keep both, and reports for every clone whether
each ITR primer's direction is confirmed, failed QC, or was never attempted -
and, when NGS data for the same material is available, whether that clone's
site was independently confirmed there too.
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
    "--validation", default=None, help="sanger_vs_ngs.tsv, if NGS data is available"
)
argparser.add_argument(
    "--original-site", default=None, help="original_site.json, if configured"
)
argparser.add_argument(
    "--original-max-dist",
    type=int,
    default=25,
    help="A site within this many bp of --original-site is unmobilized",
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

CLONE_BASE_COLUMNS = [
    "clone",
    "samples",
    "n_forward_reads",
    "n_reverse_reads",
    "forward_status",
    "reverse_status",
    "both_sides_confirmed",
]
ORIGINAL_SITE_COLUMN = "unmobilized"
CLONE_VALIDATION_COLUMNS = ["ngs_validated", "ngs_side"]

original_site = None
if args.original_site is not None:
    with open(args.original_site) as f:
        original_site = json.load(f)

CLONE_COLUMNS = (
    CLONE_BASE_COLUMNS
    + ([ORIGINAL_SITE_COLUMN] if original_site is not None else [])
    + ["summary"]
)


def at_original_site(sites, original_site, max_dist):
    return (sites["chrom"] == original_site["chrom"]) & (
        (sites["start"] - original_site["pos"]).abs() <= max_dist
    )


def classify(direction_reads, original_site=None, max_dist=0):
    """no data / failed / confirmed / unmobilized for a clone's one primer."""
    if direction_reads.shape[0] == 0:
        return "no data"
    passing = direction_reads[direction_reads["pass"].astype(bool)]
    if passing.empty:
        return "failed"
    if (
        original_site is not None
        and at_original_site(passing, original_site, max_dist).any()
    ):
        return "unmobilized"
    return "confirmed"


def describe(forward_status, reverse_status):
    """A one-line summary of a clone's forward/reverse coverage."""
    good = {"confirmed", "unmobilized"}
    label = {
        "confirmed": "confirmed",
        "unmobilized": "at the original site (unmobilized)",
    }
    if forward_status in good and reverse_status in good:
        if forward_status == reverse_status:
            return f"both sides {label[forward_status]}"
        return f"forward {forward_status}, reverse {reverse_status} (mixed)"
    if forward_status in good or reverse_status in good:
        if forward_status in good:
            side, status, other_side, other_status = (
                "forward",
                forward_status,
                "reverse",
                reverse_status,
            )
        else:
            side, status, other_side, other_status = (
                "reverse",
                reverse_status,
                "forward",
                forward_status,
            )
        return f"{side} {label[status]} ({other_side}: {other_status})"
    if forward_status == "no data" and reverse_status == "no data":
        return "no Sanger data"
    if forward_status == "failed" and reverse_status == "failed":
        return "both sides failed QC"
    return f"neither confirmed (forward: {forward_status}, reverse: {reverse_status})"


def fail_reason_summary(group):
    """'reason (n); reason (n)', most common first, for one failed group."""
    counts = group["reason"].value_counts()
    return "; ".join(f"{reason} ({n})" for reason, n in counts.items())


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
    clone_summary = pd.DataFrame(columns=CLONE_COLUMNS)
else:
    reads = pd.concat(
        (pd.read_csv(path, sep="\t", dtype={"chrom": str}) for path in args.reads),
        ignore_index=True,
    )
    reads["pass"] = tagmaplib.to_bool(reads["pass"])

    if reads.shape[0] == 0:
        qc = pd.DataFrame(columns=QC_COLUMNS)
        clone_summary = pd.DataFrame(columns=CLONE_COLUMNS)
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
            forward_status = classify(forward, original_site, args.original_max_dist)
            reverse_status = classify(reverse, original_site, args.original_max_dist)
            result = {
                "samples": ", ".join(sorted(group["sample_name"].unique())),
                "n_forward_reads": forward.shape[0],
                "n_reverse_reads": reverse.shape[0],
                "forward_status": forward_status,
                "reverse_status": reverse_status,
                "both_sides_confirmed": forward_status == "confirmed"
                and reverse_status == "confirmed",
            }
            if original_site is not None:
                statuses = {forward_status, reverse_status}
                result[ORIGINAL_SITE_COLUMN] = (
                    "confirmed" not in statuses and "unmobilized" in statuses
                )
            result["summary"] = describe(forward_status, reverse_status)
            return pd.Series(result)

        clone_summary = (
            reads.groupby("clone")
            .apply(summarize_clone, include_groups=False)
            .reset_index()[CLONE_COLUMNS]
            .sort_values("clone")
        )
        clone_summary[["n_forward_reads", "n_reverse_reads"]] = clone_summary[
            ["n_forward_reads", "n_reverse_reads"]
        ].astype(int)
        clone_summary["both_sides_confirmed"] = clone_summary[
            "both_sides_confirmed"
        ].astype(bool)
        if original_site is not None:
            clone_summary[ORIGINAL_SITE_COLUMN] = clone_summary[
                ORIGINAL_SITE_COLUMN
            ].astype(bool)

if args.validation is not None:
    validation = pd.read_csv(args.validation, sep="\t", dtype={"chrom": str})
    if validation.shape[0] and "confirmed_by_ngs" in validation.columns:
        validation["confirmed_by_ngs"] = tagmaplib.to_bool(
            validation["confirmed_by_ngs"]
        )
        clone_validation = (
            validation.groupby("clone")
            .apply(summarize_clone_validation, include_groups=False)
            .reset_index()
        )
    else:
        clone_validation = pd.DataFrame(columns=["clone"] + CLONE_VALIDATION_COLUMNS)
    clone_summary = clone_summary.merge(clone_validation, on="clone", how="left")
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
