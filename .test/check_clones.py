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
    # contaminated_by is the narrower column: the source of the contamination
    # the call itself rests on, not every secondary site traceable elsewhere.
    check(
        str(clones.loc["plate1_B01", "contaminated_by"]) == "plate2_A01",
        "plate1_B01 should name plate2_A01 as what contaminated it - got "
        f"{clones.loc['plate1_B01', 'contaminated_by']!r}",
    )
    check(
        pd.isna(clones.loc["plate1_A01", "contaminated_by"]),
        "a clean clone should have no contaminated_by - got "
        f"{clones.loc['plate1_A01', 'contaminated_by']!r}",
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
# A01 plate1: 40 molecules on each of two sides, each sequenced COPIES times,
# plus the planted single-read noise molecules.
COPIES, N_NOISE = 3, 6
expected_pairs = 40 * 2 * COPIES + N_NOISE
a01_p1 = assigned[(assigned["source"] == "A01") & (assigned["plate"] == "plate1")]
check(
    a01_p1["n_pairs"].sum() == expected_pairs,
    f"A01/plate1 should have {expected_pairs} pairs, got {a01_p1['n_pairs'].sum()}",
)
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

# min_reads_per_molecule: the planted noise is one read per molecule and every
# real molecule was sequenced several times, so none of the noise sites should
# have survived into the calls. NOISE_START/NOISE_SPACING in make_test_data.py.
NOISE_START, NOISE_SPACING = 4000, 1500
noise_positions = {NOISE_START + i * NOISE_SPACING for i in range(N_NOISE)}
support = pd.read_csv(
    f"{RESULTS}/insertion_sites/all_sites_support.tsv", sep="\t", dtype={"chrom": str}
)
survived = [
    (r.chrom, r.start)
    for r in support.itertuples()
    if r.chrom == "chr1" and any(abs(r.start - p) <= 50 for p in noise_positions)
]
check(
    not survived,
    f"single-read noise molecules should have been dropped by "
    f"min_reads_per_molecule, but sites remain at {survived}",
)
# ... and the real molecules, sequenced several times over, should not have been
check(
    clones.loc["plate1_A01", "dominant_molecules"] == 80,
    "plate1_A01 should keep all 80 of its molecules (40 per ITR side) through "
    f"the read-support filter, got {clones.loc['plate1_A01', 'dominant_molecules']}",
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
