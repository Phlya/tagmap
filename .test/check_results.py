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

# Clone C01 (make_test_data.py) crosses the whole of this short standalone
# contig and keeps going into unrelated real genome - it has no NGS coverage
# by design (this failure mode is Sanger-only), so it sits outside
# EXPECTED/the NGS-confirmed counts everywhere below, but should still show
# up as one real, passing, one-sided Sanger site.
LAUNCHPAD_CONTIG = "test_launchpad"

problems = []


def check(condition, message):
    if not condition:
        problems.append(message)


ngs = pd.read_csv("results/insertion_sites/all_sites.bed", sep="\t")
sanger = pd.read_csv("results/sanger/all_sanger_sites.bed", sep="\t")
# A plain, format-compliant BED - no header, no extra columns, name is
# sample_name + "_" + clone - so it must be read positionally.
confirmed_sanger = pd.read_csv(
    "results/sanger/confirmed_sanger_sites.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "name", "score", "strand"],
)
validation = pd.read_csv("results/validation/sanger_vs_ngs.tsv", sep="\t")

# filter_confirmed_sites.py: A01 (both primers agree) and B01 (one-sided in
# Sanger, but NGS backs it) both count as confirmed; C01's launchpad-escape
# site is one-sided with no NGS data at all, so it's the one site dropped.
check(
    set(zip(confirmed_sanger["name"], confirmed_sanger["chrom"]))
    == {("plate1_A01", "chr1"), ("plate1_B01", "chr1")},
    "confirmed_sanger_sites.bed should keep A01 (both primers) and B01 "
    f"(NGS-backed) but drop C01 (neither) - got names {sorted(confirmed_sanger['name'])}",
)

check(ngs.shape[0] == EXPECTED.shape[0], f"expected {EXPECTED.shape[0]} NGS sites, got {ngs.shape[0]}")
# +1: clone C01's test_launchpad site is a real, passing Sanger site, just
# not a "new integration" NGS also confirms, so it's outside EXPECTED.
check(
    sanger.shape[0] == EXPECTED.shape[0] + 1,
    f"expected {EXPECTED.shape[0] + 1} Sanger sites (including C01's "
    f"test_launchpad site), got {sanger.shape[0]}",
)

# Direct regression check for sanger_sites.py's runs_off_contig_end: C01's
# read crosses the whole of the short test_launchpad contig and keeps going
# into unrelated real genome, which used to make analyse_read() see a second
# qualifying genomic segment and fail the read as "multiple genomic
# alignments" even though nothing ambiguous happened - the contig just ran
# out. It should pass QC instead.
plate1_reads = pd.read_csv("results/sanger/plate1_reads.tsv", sep="\t")
c01_reads = plate1_reads[plate1_reads["clone"] == "C01"]
check(
    c01_reads.shape[0] == 1
    and str(c01_reads["pass"].iloc[0]) == "True"
    and "multiple genomic alignments" not in str(c01_reads["fail_reason"].iloc[0]),
    "clone C01's read should pass QC rather than being flagged as multiple "
    "genomic alignments for outrunning the short test_launchpad contig",
)

# sanger_read_qc.tsv's at_original_site - the well-level source for the "u"
# marker on the read-QC plate grids: only A01 (both sides confirmed at the
# configured original insertion site, chr1:8001) should be True; B01 and
# C01 are both mobilized elsewhere, and B01/C01's reverse was never
# sequenced at all.
read_qc = pd.read_csv("results/stats/sanger_read_qc.tsv", sep="\t")
at_site = read_qc.set_index(["well", "direction"])["at_original_site"]
check(
    bool(at_site[("A01", "forward")])
    and bool(at_site[("A01", "reverse")])
    and not at_site[("B01", "forward")]
    and not at_site[("B01", "reverse")]
    and not at_site[("C01", "forward")]
    and not at_site[("C01", "reverse")],
    "sanger_read_qc.tsv should mark only A01 (forward and reverse) as "
    "at_original_site",
)

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

# Excludes C01's test_launchpad site, which has no NGS coverage by design.
confirmable = validation[validation["chrom"] != LAUNCHPAD_CONTIG]
check(
    bool(confirmable["confirmed_by_ngs"].all()),
    "every Sanger site (other than C01's launchpad-escape site) should have "
    "been confirmed by the NGS library",
)
check(
    bool(confirmable["ngs_two_sided"].all()),
    "every matching NGS site should have been two-sided",
)

# The new stats/ tables should agree with what the checks above already
# established directly from the site files. NGS mobilization and sidedness
# are now combined into one file, keyed by sample_name.
ngs_qc = pd.read_csv("results/stats/ngs_qc_stats.tsv", sep="\t")
clone1_qc = ngs_qc[ngs_qc["sample_name"] == "clone1"]
check(
    clone1_qc.shape[0] == 1
    and clone1_qc["mobilized_pairs"].iloc[0] > 0
    and clone1_qc["n_sites"].iloc[0] == EXPECTED.shape[0]
    and clone1_qc["n_two_sided"].iloc[0] == EXPECTED.shape[0]
    and clone1_qc["n_one_sided"].iloc[0] == 0,
    "ngs_qc_stats.tsv should show mobilized pairs and both sites two-sided for clone1",
)

# sanger_qc_stats.tsv's n_clones_pass counts clones, not reads: all three
# clones (A01, B01, C01) have a passing forward read, but only A01 was ever
# sequenced from the reverse primer at all (B01 and C01 are forward-only by
# design), so reverse should show 1 of 3 clones passed.
sanger_qc = pd.read_csv("results/stats/sanger_qc_stats.tsv", sep="\t")
forward_qc = sanger_qc[sanger_qc["direction"] == "forward"]
reverse_qc = sanger_qc[sanger_qc["direction"] == "reverse"]
check(
    forward_qc.shape[0] == 1
    and forward_qc["n_clones"].iloc[0] == 3
    and forward_qc["n_clones_pass"].iloc[0] == 3
    and reverse_qc.shape[0] == 1
    and reverse_qc["n_clones"].iloc[0] == 3
    and reverse_qc["n_clones_pass"].iloc[0] == 1,
    "sanger_qc_stats.tsv should show 3 of 3 clones passed forward and 1 of 3 "
    "passed reverse (only A01 was ever sequenced from the reverse primer)",
)

# A01 was sequenced from both primers, B01 only ever had a forward read - the
# reverse primer was never attempted for it, rather than attempted and failed.
# Both clones' sites matched a two-sided NGS site (checked above). A01 is
# already "both" from Sanger alone; B01's ngs_verification should ALSO read
# "both" - not because Sanger saw its reverse side, but because NGS fills in
# exactly the side Sanger's own primers left unconfirmed (see
# tagmaplib.ngs_verification_label). The test config also points
# original_insertion_site at A01's own planted insertion (chr1:8000): A01
# still has a real, confirmed site there (from both primers, agreeing on
# position), it just happens to be the founder's own locus, so it is
# both_sides_confirmed and its summary reads "unmobilized" instead of "both
# sides confirmed" (there is no separate unmobilized column - both primers
# have to agree on the site for the label to apply at all).
sanger_clones = pd.read_csv("results/stats/sanger_clone_summary.tsv", sep="\t")
a01 = sanger_clones[sanger_clones["clone"] == "A01"]
check(
    a01.shape[0] == 1
    and bool(a01["both_sides_confirmed"].iloc[0])
    and a01["forward_status"].iloc[0] == "PASSED"
    and a01["reverse_status"].iloc[0] == "PASSED"
    and bool(a01["positions_agree"].iloc[0])
    and a01["summary"].iloc[0] == "unmobilized"
    and a01["ngs_verification"].iloc[0] == "both",
    "sanger_clone_summary.tsv should show clone A01 confirmed from both "
    "primers, agreeing on position, at the original insertion site "
    "(summary: unmobilized), and validated by a two-sided NGS site",
)
b01 = sanger_clones[sanger_clones["clone"] == "B01"]
check(
    b01.shape[0] == 1
    and not bool(b01["both_sides_confirmed"].iloc[0])
    and b01["forward_status"].iloc[0] == "PASSED"
    and b01["reverse_status"].iloc[0] == "not sequenced"
    and b01["summary"].iloc[0] != "unmobilized"
    and b01["ngs_verification"].iloc[0] == "both",
    "sanger_clone_summary.tsv should show clone B01 confirmed forward-only at "
    "a different (mobilized) locus, with the reverse primer never attempted, "
    "but still validated by a two-sided NGS site",
)

# A01 and B01 sit at two different loci (checked above), both sequenced on
# the one plate configured here, so plate1 and the "all" total should agree
# and both equal the number of planted insertions.
sanger_positions = pd.read_csv("results/stats/sanger_positions.tsv", sep="\t")
plate1_positions = sanger_positions[sanger_positions["sample_name"] == "plate1"]
all_positions = sanger_positions[sanger_positions["sample_name"] == "all"]
# +1: C01's test_launchpad locus, a third distinct position on plate1 (and
# so also in the "all" total) alongside A01's and B01's.
check(
    plate1_positions.shape[0] == 1
    and plate1_positions["n_distinct_positions"].iloc[0] == EXPECTED.shape[0] + 1
    and all_positions.shape[0] == 1
    and all_positions["n_distinct_positions"].iloc[0] == EXPECTED.shape[0] + 1,
    "sanger_positions.tsv should show "
    f"{EXPECTED.shape[0] + 1} distinct positions for plate1 and for the total",
)
# n_one_sided: B01 (forward only) and C01 (forward only) - A01 is confirmed
# from both primers, so it's the only one of the three not one-sided.
# n_ngs_confirmed: A01 and B01 are both confirmed by the NGS library; C01's
# launchpad-escape site has no NGS coverage by design (see LAUNCHPAD_CONTIG).
check(
    plate1_positions["n_one_sided"].iloc[0] == 2
    and all_positions["n_one_sided"].iloc[0] == 2
    and plate1_positions["n_ngs_confirmed"].iloc[0] == 2
    and all_positions["n_ngs_confirmed"].iloc[0] == 2,
    "sanger_positions.tsv should show 2 one-sided positions (B01, C01) and "
    "2 NGS-confirmed positions (A01, B01) for plate1 and for the total",
)

# A01 is unmobilized (checked below), so it's folded into one combined
# "unmobilized" row instead of appearing as its own chr1:8001 position -
# only B01's real, mobilized site shows up in the per-locus breakdown here.
# Plus the two running totals: 2 usable clones (A01+B01), 1 mobilized
# position (B01's). A01 was seen from both primers (n_one_sided=0 for
# unmobilized), B01 from forward only (n_one_sided=1 for chr1's row).
# Labels must match sanger_stats.py's constants.
UNMOBILIZED_LABEL = "unmobilized: at the original insertion site"
TOTAL_CLONES_LABEL = "clones with at least one good read"
TOTAL_MOBILIZED_POSITIONS_LABEL = "total mobilized positions"
# C01's test_launchpad locus adds a third clone and a second real, distinct
# mobilized position - one-sided in Sanger, and (unlike B01) not completed
# by NGS since it has none, so its own ngs_verification reads "forward only"
# rather than "both". It is not a founder/unmobilized site, so it does not
# touch the A01/UNMOBILIZED_LABEL row.
position_counts = pd.read_csv("results/stats/sanger_position_counts.tsv", sep="\t")
by_label = position_counts.set_index("chrom")[
    ["n_times_found", "n_one_sided", "ngs_verification"]
]
check(
    position_counts.shape[0] == 5
    and tuple(by_label.loc["chr1"]) == (1, 1, "both")
    and tuple(by_label.loc[LAUNCHPAD_CONTIG]) == (1, 1, "forward only")
    and tuple(by_label.loc[UNMOBILIZED_LABEL]) == (1, 0, "both")
    and by_label.loc[TOTAL_CLONES_LABEL, "n_times_found"] == 3
    and by_label.loc[TOTAL_MOBILIZED_POSITIONS_LABEL, "n_times_found"] == 2,
    "sanger_position_counts.tsv should list B01's one-sided-in-Sanger-but-"
    "NGS-completed-to-both site, C01's one-sided-and-NGS-unconfirmed "
    "launchpad-escape site, one combined two-sided unmobilized row for A01, "
    "and the 3-clones/2-positions totals",
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
