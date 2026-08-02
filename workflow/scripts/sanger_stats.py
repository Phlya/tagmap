"""Summarise Sanger sequencing coverage and QC outcomes per clone.

Reads that pass QC already end up in the per-clone consensus sites
(combine_sanger_sites.py), but a read that failed, or a primer direction that
was never sequenced for a clone, simply has no site to show for it - so
"missing" and "failed" look the same downstream. This works from the raw
per-read tables instead, which keep both, and reports for every clone whether
each ITR primer's direction is confirmed, failed QC, or was never attempted.
"""

import argparse

import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--reads", nargs="*", default=[], help="{sample}_reads.tsv files")
argparser.add_argument("--output-read-summary", required=True)
argparser.add_argument("--output-fail-reasons", required=True)
argparser.add_argument("--output-clone-summary", required=True)
args = argparser.parse_args()

READ_SUMMARY_COLUMNS = [
    "sample_name",
    "direction",
    "n_reads",
    "n_pass",
    "n_fail",
    "frac_pass",
]

FAIL_REASON_COLUMNS = ["direction", "reason", "n_reads"]

CLONE_COLUMNS = [
    "clone",
    "samples",
    "n_forward_reads",
    "n_reverse_reads",
    "forward_status",
    "reverse_status",
    "both_sides_confirmed",
    "summary",
]


def classify(direction_reads):
    """no data / failed / confirmed for one clone's reads from one primer."""
    if direction_reads.shape[0] == 0:
        return "no data"
    if direction_reads["pass"].astype(bool).any():
        return "confirmed"
    return "failed"


def describe(forward_status, reverse_status):
    """A one-line summary of a clone's forward/reverse coverage."""
    if forward_status == "confirmed" and reverse_status == "confirmed":
        return "both sides confirmed"
    if forward_status == "confirmed" or reverse_status == "confirmed":
        if forward_status == "confirmed":
            confirmed_side, other_side, other_status = "forward", "reverse", reverse_status
        else:
            confirmed_side, other_side, other_status = "reverse", "forward", forward_status
        return f"{confirmed_side} only confirmed ({other_side}: {other_status})"
    if forward_status == "no data" and reverse_status == "no data":
        return "no Sanger data"
    if forward_status == "failed" and reverse_status == "failed":
        return "both sides failed QC"
    return f"neither confirmed (forward: {forward_status}, reverse: {reverse_status})"


if not args.reads:
    read_summary = pd.DataFrame(columns=READ_SUMMARY_COLUMNS)
    fail_reasons = pd.DataFrame(columns=FAIL_REASON_COLUMNS)
    clone_summary = pd.DataFrame(columns=CLONE_COLUMNS)
else:
    reads = pd.concat(
        (pd.read_csv(path, sep="\t", dtype={"chrom": str}) for path in args.reads),
        ignore_index=True,
    )
    reads["pass"] = tagmaplib.to_bool(reads["pass"])

    if reads.shape[0] == 0:
        read_summary = pd.DataFrame(columns=READ_SUMMARY_COLUMNS)
        fail_reasons = pd.DataFrame(columns=FAIL_REASON_COLUMNS)
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

        read_summary = (
            reads.groupby(["sample_name", "direction"])
            .apply(summarize_reads, include_groups=False)
            .reset_index()[READ_SUMMARY_COLUMNS]
            .sort_values(["sample_name", "direction"])
        )
        # groupby.apply returns one Series per group; mixing ints and the
        # float frac_pass in that Series forces the whole thing to float64.
        read_summary[["n_reads", "n_pass", "n_fail"]] = read_summary[
            ["n_reads", "n_pass", "n_fail"]
        ].astype(int)

        failed = reads[~reads["pass"]].copy()
        failed["fail_reason"] = failed["fail_reason"].fillna("unknown").astype(str)
        failed.loc[failed["fail_reason"] == "", "fail_reason"] = "unknown"
        exploded = failed.assign(reason=failed["fail_reason"].str.split("; ")).explode(
            "reason"
        )
        fail_reasons = (
            exploded.groupby(["direction", "reason"])
            .size()
            .reset_index(name="n_reads")[FAIL_REASON_COLUMNS]
            .sort_values(["direction", "n_reads"], ascending=[True, False])
        )

        def summarize_clone(group):
            forward = group[group["direction"] == "forward"]
            reverse = group[group["direction"] == "reverse"]
            forward_status = classify(forward)
            reverse_status = classify(reverse)
            return pd.Series(
                {
                    "samples": ", ".join(sorted(group["sample_name"].unique())),
                    "n_forward_reads": forward.shape[0],
                    "n_reverse_reads": reverse.shape[0],
                    "forward_status": forward_status,
                    "reverse_status": reverse_status,
                    "both_sides_confirmed": forward_status == "confirmed"
                    and reverse_status == "confirmed",
                    "summary": describe(forward_status, reverse_status),
                }
            )

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

read_summary.to_csv(args.output_read_summary, sep="\t", index=False)
fail_reasons.to_csv(args.output_fail_reasons, sep="\t", index=False)
clone_summary.to_csv(args.output_clone_summary, sep="\t", index=False)

print(
    f"{clone_summary.shape[0]} clones: "
    f"{int(clone_summary['both_sides_confirmed'].sum()) if clone_summary.shape[0] else 0} "
    "confirmed from both primers"
)
