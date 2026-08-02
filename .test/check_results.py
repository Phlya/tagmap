"""Check that the workflow found the insertions make_test_data.py planted.

Run from the .test directory after the workflow:  python3 check_results.py
Exits non-zero with a description of what is wrong, so CI fails loudly rather
than on a missing file.
"""

import os
import sys

import pandas as pd

EXPECTED = pd.read_csv(
    "resources/expected_sites.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "planted_strand"],
)

problems = []


def check(condition, message):
    if not condition:
        problems.append(message)


ngs = pd.read_csv("results/insertion_sites/all_sites.bed", sep="\t")
sanger = pd.read_csv("results/sanger/all_sanger_sites.bed", sep="\t")
validation = pd.read_csv("results/validation/sanger_vs_ngs.tsv", sep="\t")

check(ngs.shape[0] == EXPECTED.shape[0], f"expected {EXPECTED.shape[0]} NGS sites, got {ngs.shape[0]}")
check(sanger.shape[0] == EXPECTED.shape[0], f"expected {EXPECTED.shape[0]} Sanger sites, got {sanger.shape[0]}")

for _, site in EXPECTED.iterrows():
    ngs_hit = ngs[(ngs["chrom"] == site["chrom"]) & (ngs["start"] == site["start"])]
    check(ngs_hit.shape[0] == 1, f"no NGS site at {site['chrom']}:{site['start']}")
    if ngs_hit.shape[0] == 1:
        check(
            bool(ngs_hit["TA_found"].iloc[0]),
            f"NGS site at {site['chrom']}:{site['start']} did not land on a TA",
        )

    sanger_hit = sanger[
        (sanger["chrom"] == site["chrom"]) & (sanger["start"] == site["start"])
    ]
    check(sanger_hit.shape[0] == 1, f"no Sanger site at {site['chrom']}:{site['start']}")
    if sanger_hit.shape[0] == 1 and ngs_hit.shape[0] == 1:
        check(
            sanger_hit["strand"].iloc[0] == ngs_hit["strand"].iloc[0],
            f"branches disagree on the strand at {site['chrom']}:{site['start']}: "
            f"Sanger {sanger_hit['strand'].iloc[0]}, NGS {ngs_hit['strand'].iloc[0]}",
        )

# Clone A01 was sequenced from both primers, B01 from only the forward one -
# which is the case the NGS library has to vouch for.
both = sanger[sanger["clone"] == "A01"]
check(both.shape[0] == 1 and bool(both["both_directions"].iloc[0]),
      "clone A01 should have been seen from both ITR primers")
one_sided = sanger[sanger["clone"] == "B01"]
check(one_sided.shape[0] == 1 and not bool(one_sided["both_directions"].iloc[0]),
      "clone B01 should have been seen from one ITR primer only")

check(
    bool(validation["confirmed_by_ngs"].all()),
    "every Sanger site should have been confirmed by the NGS library",
)
check(
    bool(validation["ngs_two_sided"].all()),
    "every matching NGS site should have been two-sided",
)

# The new stats/ tables should agree with what the checks above already
# established directly from the site files.
ngs_sidedness = pd.read_csv("results/stats/ngs_site_sidedness.tsv", sep="\t")
clone1_sidedness = ngs_sidedness[ngs_sidedness["sample_name"] == "clone1"]
check(
    clone1_sidedness.shape[0] == 1
    and clone1_sidedness["n_sites"].iloc[0] == EXPECTED.shape[0]
    and clone1_sidedness["n_two_sided"].iloc[0] == EXPECTED.shape[0]
    and clone1_sidedness["n_one_sided"].iloc[0] == 0,
    "ngs_site_sidedness.tsv should show both sites as two-sided for clone1",
)

ngs_mapping = pd.read_csv("results/stats/ngs_mapping_stats.tsv", sep="\t")
clone1_mapping = ngs_mapping[ngs_mapping["sample_name"] == "clone1"]
check(
    clone1_mapping.shape[0] == 1 and clone1_mapping["mobilized_pairs"].iloc[0] > 0,
    "ngs_mapping_stats.tsv should show mobilized pairs for clone1",
)

# A01 was sequenced from both primers, B01 only ever had a forward read - the
# reverse primer was never attempted for it, rather than attempted and failed.
sanger_clones = pd.read_csv("results/stats/sanger_clone_summary.tsv", sep="\t")
a01 = sanger_clones[sanger_clones["clone"] == "A01"]
check(
    a01.shape[0] == 1
    and bool(a01["both_sides_confirmed"].iloc[0])
    and a01["forward_status"].iloc[0] == "confirmed"
    and a01["reverse_status"].iloc[0] == "confirmed",
    "sanger_clone_summary.tsv should show clone A01 confirmed from both primers",
)
b01 = sanger_clones[sanger_clones["clone"] == "B01"]
check(
    b01.shape[0] == 1
    and not bool(b01["both_sides_confirmed"].iloc[0])
    and b01["forward_status"].iloc[0] == "confirmed"
    and b01["reverse_status"].iloc[0] == "no data",
    "sanger_clone_summary.tsv should show clone B01 confirmed forward-only, "
    "with the reverse primer never attempted",
)

validation_summary = pd.read_csv("results/stats/validation_summary.tsv", sep="\t")
totals = validation_summary[validation_summary["sample_name"] == "all"]
check(
    totals.shape[0] == 1 and totals["n_confirmed"].iloc[0] == EXPECTED.shape[0],
    "validation_summary.tsv should show every Sanger site confirmed",
)

check(os.path.getsize("results/stats/report.md") > 0, "report.md should not be empty")

if problems:
    print("FAILED:")
    for problem in problems:
        print(f"  - {problem}")
    sys.exit(1)
print(f"OK: {EXPECTED.shape[0]} insertions found by both branches and cross-confirmed")
