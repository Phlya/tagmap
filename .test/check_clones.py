"""Check the barcoded, per-clone run (config/config_clones.yaml).

Run from the .test directory after the workflow:  python3 check_clones.py
Expectations come from CLONE_LIBRARIES in make_test_data.py.
"""

import os
import sys

import pandas as pd

problems = []


def check(condition, message):
    if not condition:
        problems.append(message)


RESULTS = "results_clones"

demux = pd.read_csv(f"{RESULTS}/stats/demux_stats.tsv", sep="\t")
clones = pd.read_csv(f"{RESULTS}/stats/ngs_clone_summary.tsv", sep="\t").set_index(
    "sample_name"
)

# Every (plate, well) is a sample of its own, named "{plate}_{name}"
expected_status = {
    "plate1_A01": "clean",
    "plate2_A01": "unmobilized",
    "plate1_B01": "contaminated",
    "plate2_B01": "multiple",
    "plate1_C01": "too_few_reads",
    "plate2_C01": "too_few_reads",
}
check(
    set(clones.index) == set(expected_status),
    f"clone samples should be {sorted(expected_status)}, got {sorted(clones.index)}",
)
for sample, status in expected_status.items():
    if sample in clones.index:
        check(
            clones.loc[sample, "status"] == status,
            f"{sample} should be {status}, got {clones.loc[sample, 'status']} "
            f"({clones.loc[sample, 'reason']})",
        )

# Dominant sites: insertion 0 is chr1:8001, insertion 1 chr1:17001
for sample, start in {
    "plate1_A01": 8001,
    "plate2_A01": 17001,
    "plate1_B01": 8001,
}.items():
    if sample in clones.index:
        check(
            clones.loc[sample, "chrom"] == "chr1"
            and clones.loc[sample, "start"] == start,
            f"{sample} dominant site should be chr1:{start}, got "
            f"{clones.loc[sample, 'chrom']}:{clones.loc[sample, 'start']}",
        )

# The contaminant is 4 of 84 molecules - the runner-up, but well under half
if "plate1_B01" in clones.index:
    second = clones.loc["plate1_B01", "second_frac"]
    check(
        0.02 <= second < 0.2,
        f"plate1_B01's runner-up should hold a few percent, got {second}",
    )
    check(
        clones.loc["plate1_B01", "second_position"] == "chr1:17001",
        "plate1_B01's runner-up should be chr1:17001",
    )
if "plate2_B01" in clones.index:
    check(
        clones.loc["plate2_B01", "n_sites"] == 2,
        "plate2_B01 should carry two sites",
    )
    check(
        clones.loc["plate2_B01", "second_frac"] >= 0.2,
        "plate2_B01's two integrations should be comparable in strength",
    )

# plate1_B01's contaminant is the founder locus, so the reason should say so -
# unmobilized carry-over, not another clone's integration.
if "plate1_B01" in clones.index:
    check(
        "donor locus" in str(clones.loc["plate1_B01", "reason"]),
        "plate1_B01's contaminant is the donor locus and the reason should say "
        f"so - got {clones.loc['plate1_B01', 'reason']!r}",
    )
    check(
        clones.loc["plate1_B01", "founder_frac"] > 0,
        "plate1_B01 should report a non-zero founder_frac",
    )
# A clean clone has nothing else
if "plate1_A01" in clones.index:
    check(clones.loc["plate1_A01", "n_sites"] == 1, "plate1_A01 should carry one site")
    check(
        clones.loc["plate1_A01", "sides"] == "both",
        "plate1_A01's site should be seen from both ITR sides",
    )
# plate1_B01's contaminant is plate2_A01's dominant site
if "plate1_B01" in clones.index:
    check(
        "plate2_A01" in str(clones.loc["plate1_B01", "shared_with"]),
        "plate1_B01's contaminant should be traced to plate2_A01, whose "
        f"dominant site it is - got shared_with={clones.loc['plate1_B01', 'shared_with']!r}",
    )

# Demultiplexing: reads were assigned to the right plates, none to the wrong
# ones, and the barcode-less strays (3 per library) were left unassigned
assigned = demux[demux["category"] == "assigned"]
check(
    set(assigned["plate"]) == {"plate1", "plate2"},
    "demux should have assigned reads to plate1 and plate2",
)
for source in ("A01", "B01", "C01"):
    stray = demux[(demux["source"] == source) & (demux["category"] == "no_barcode")]
    check(
        stray["n_pairs"].sum() == 3,
        f"library {source} should have 3 barcode-less pairs, got "
        f"{stray['n_pairs'].sum()}",
    )
# A01 plate1: 40 molecules on each of two sides
a01_p1 = assigned[(assigned["source"] == "A01") & (assigned["plate"] == "plate1")]
check(a01_p1["n_pairs"].sum() == 80, f"A01/plate1 should have 80 pairs, got {a01_p1['n_pairs'].sum()}")
check(
    set(a01_p1.loc[a01_p1["n_pairs"] > 0, "side"]) == {"forward", "reverse"},
    "A01/plate1 should have pairs from both barcode sides",
)
# One planted barcode mismatch per (library, plate) with reads
check(
    assigned[assigned["mismatches"] > 0]["n_pairs"].sum() >= 1,
    "at least one pair should have matched a barcode with a mismatch",
)

# Rows should follow the sample sheet, not an alphabetical re-sort.
sheet = pd.read_csv("config/samples_clones.tsv", sep="\t", comment="#")
expected_order = [
    f"{plate}_{name}"
    for name in dict.fromkeys(sheet["name"])
    for plate in ("plate1", "plate2")
]
check(
    list(clones.index) == [s for s in expected_order if s in set(clones.index)],
    f"clone rows should follow the sample sheet's order, got {list(clones.index)}",
)

for name in ("report.md", "report.pdf"):
    path = f"{RESULTS}/stats/{name}"
    check(os.path.exists(path) and os.path.getsize(path) > 0, f"{name} should not be empty")
if os.path.exists(f"{RESULTS}/stats/report.md"):
    report = open(f"{RESULTS}/stats/report.md").read()
    check("## NGS clones" in report and "## Demultiplexing" in report, "report.md should have the demux and clone sections")

if problems:
    print("FAILED:")
    for problem in problems:
        print(f"  - {problem}")
    sys.exit(1)
print(f"OK: {len(expected_status)} clones called as expected, demultiplexing correct")
