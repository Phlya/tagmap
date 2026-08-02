"""Summarise the Sanger-vs-NGS cross-validation table into per-run counts.

compare_sanger_ngs.py reports one row per Sanger site with the NGS evidence
found for it; this collapses that into how many sites each Sanger run
confirmed, plus a totals row.
"""

import argparse

import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--validation", required=True, help="sanger_vs_ngs.tsv")
argparser.add_argument("--output", "-o", required=True)
args = argparser.parse_args()

COLUMNS = [
    "sample_name",
    "n_sanger_sites",
    "n_matched_ngs_site",
    "n_confirmed",
    "frac_confirmed",
    "n_confirmed_two_sided",
]

validation = pd.read_csv(args.validation, sep="\t", dtype={"chrom": str})

if validation.shape[0] == 0:
    summary = pd.DataFrame(columns=COLUMNS)
else:
    validation["confirmed_by_ngs"] = tagmaplib.to_bool(validation["confirmed_by_ngs"])
    validation["ngs_two_sided"] = tagmaplib.to_bool(validation["ngs_two_sided"])

    def summarize(group):
        n = group.shape[0]
        confirmed = int(group["confirmed_by_ngs"].sum())
        return pd.Series(
            {
                "n_sanger_sites": n,
                "n_matched_ngs_site": int(group["ngs_site_dist"].notna().sum()),
                "n_confirmed": confirmed,
                "frac_confirmed": confirmed / n if n else float("nan"),
                "n_confirmed_two_sided": int(group["ngs_two_sided"].sum()),
            }
        )

    summary = (
        validation.groupby("sample_name")
        .apply(summarize, include_groups=False)
        .reset_index()[COLUMNS]
        .sort_values("sample_name")
    )
    totals = summarize(validation)
    totals["sample_name"] = "all"
    summary = pd.concat([summary, pd.DataFrame([totals])[COLUMNS]], ignore_index=True)
    # summarize returns one Series mixing ints and a float, which forces the
    # whole Series - and so every row built from it - to float64.
    int_columns = ["n_sanger_sites", "n_matched_ngs_site", "n_confirmed", "n_confirmed_two_sided"]
    summary[int_columns] = summary[int_columns].astype(int)

summary.to_csv(args.output, sep="\t", index=False)

print(
    f"{validation.shape[0]} Sanger sites, "
    f"{int(validation['confirmed_by_ngs'].sum()) if validation.shape[0] else 0} confirmed by NGS"
)
