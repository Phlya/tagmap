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

import bioframe
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
argparser.add_argument(
    "--region",
    default=None,
    help="validated_clones_region (chrom:start-end), if configured - flags "
    "each clone's in_region column when its own site(s) overlap it, for "
    "report_pdf.py to draw a bold border around validated clones there.",
)
argparser.add_argument("--output-qc", required=True)
argparser.add_argument("--output-clone-summary", required=True)
argparser.add_argument("--output-positions", required=True)
argparser.add_argument("--output-position-counts", required=True)
argparser.add_argument("--output-read-qc", required=True)
args = argparser.parse_args()

QC_COLUMNS = [
    "sample_name",
    "direction",
    "n_reads",
    "n_pass",
    "n_fail",
    "frac_pass",
    "n_clones",
    "n_clones_pass",
    "frac_clones_pass",
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
CLONE_VALIDATION_COLUMNS = ["ngs_verification"]
REGION_COLUMNS = ["in_region"]

# Not one of the regex groups sanger_sites.py treats specially (well, clone,
# direction) - just a convention some configs' sanger_name_regex use to name
# the plate/submission id embedded in Sanger trace filenames (see e.g.
# .test/config/config.yaml). Kept optional throughout: most projects' clone
# ids don't come with one at all.
PLATE_RUN_COLUMN = "run"
plate_runs = None

original_site = None
if args.original_site is not None:
    with open(args.original_site) as f:
        original_site = json.load(f)

region = bioframe.parse_region(args.region, check_bounds=False) if args.region else None

CLONE_COLUMNS = BASE_CLONE_COLUMNS + [
    "position",
    "orientation",
    "positions_agree",
    "both_sides_confirmed",
    "summary",
]


def at_original_site(sites, original_site, max_dist):
    return (sites["chrom"] == original_site["chrom"]) & (
        (sites["start"] - original_site["pos"]).abs() <= max_dist
    )


def in_region(sites, region):
    chrom, start, end = region
    return (sites["chrom"] == chrom) & (sites["start"] < end) & (sites["end"] > start)


def classify(direction_reads):
    """no data / failed / confirmed for one clone's reads from one primer."""
    if direction_reads.shape[0] == 0:
        return "no data"
    if direction_reads["pass"].astype(bool).any():
        return "confirmed"
    return "failed"


def describe(forward_status, reverse_status, positions_agree, unmobilized=False):
    """A one-line summary of a clone's forward/reverse coverage.

    forward_status/reverse_status are "PASSED", "not sequenced", or the QC
    failure reason itself - so a passing side is identified by that exact
    string, not by absence of a separate reason column.

    Both sides can independently pass QC while pointing at different loci -
    e.g. one read stuck on a residual copy of the donor construct - so
    "PASSED" from both primers isn't enough; their sites also have to
    cluster together (positions_agree, from the already-clustered
    all_sanger_sites.bed) for the clone to really be confirmed.

    unmobilized (only ever True alongside positions_agree - see
    summarize_clone_sites) replaces "both sides confirmed" with
    "unmobilized" rather than adding a separate column: it is just the
    special case of a confirmed site that happens to sit at the
    pre-mobilization founder locus.
    """
    forward_ok = forward_status == "PASSED"
    reverse_ok = reverse_status == "PASSED"
    if forward_ok and reverse_ok:
        if not positions_agree:
            return "forward and reverse confirmed but disagree on position"
        if unmobilized:
            return "unmobilized"
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


CLONE_SIDE_COLUMNS = ["sample_name", "direction", "n_clones", "n_clones_pass"]


def summarize_clone_sides(clone_summary):
    """How many clones had each ITR primer's side pass QC, per sample - by
    clone, not by read: several reads for the same clone and direction (a
    rerun, or - for a project whose sanger_name_regex captures a primer pair
    as its own clone - a different primer pair tried on the same well) were
    already folded into one forward_status/reverse_status per clone
    (side_status), so a clone counts as passed here as soon as any of its
    reads on that side did, same as the per-clone table itself.
    """
    if clone_summary.shape[0] == 0:
        return pd.DataFrame(columns=CLONE_SIDE_COLUMNS)
    n_clones = clone_summary.groupby("sample_name").size().reset_index(name="n_clones")
    counts = []
    for direction, status_col in (
        ("forward", "forward_status"),
        ("reverse", "reverse_status"),
    ):
        passed = (
            clone_summary[status_col]
            .eq("PASSED")
            .groupby(clone_summary["sample_name"])
            .sum()
            .reset_index(name="n_clones_pass")
        )
        passed["direction"] = direction
        counts.append(passed.merge(n_clones, on="sample_name"))
    result = pd.concat(counts, ignore_index=True)
    result["n_clones_pass"] = result["n_clones_pass"].astype(int)
    return result[CLONE_SIDE_COLUMNS]


def format_site(row):
    """One site cluster as 'chrom:start-end(strand, Nf forward + Nr reverse
    reads)'. Spelling out the per-primer breakdown, rather than just the
    total, is what lets a reader tell which side actually supports this
    particular locus - the distinction that matters when a clone's sites
    disagree (see summarize_clone_sites) and position lists more than one.
    """
    return (
        f"{row.chrom}:{row.start}-{row.end}({row.strand}, "
        f"{row.n_forward} forward + {row.n_reverse} reverse reads)"
    )


def summarize_clone_sites(group, original_site=None, max_dist=0, region=None):
    """A clone's site(s), whether they agree on position, and - when the
    original, pre-mobilization locus is configured - whether the
    both-primers-agree site sits there rather than at a new locus.

    Usually one site backed by both directions, but a clone whose forward
    and reverse reads don't cluster together (see combine_sanger_sites.py)
    ends up with more than one single-direction site here - listing all of
    them shows where the disagreement actually is, rather than just that
    there is one.

    A single-direction site landing near the original locus doesn't count:
    only a site both primers agree on (both_directions) can call a clone
    unmobilized, otherwise a stray read on a residual donor copy could mislabel
    a genuinely mobilized clone whose other primer confirms a real, different
    site.

    in_region checks every site listed in position, not just an
    agreeing/both_directions one - a clone validated via a single primer plus
    NGS (see tagmaplib.clone_color) never has an agreeing site at all, but
    its one listed site is still its real evidence.
    """
    sites = group.sort_values(["chrom", "start"])
    result = {
        "position": "; ".join(format_site(row) for row in sites.itertuples()),
        # Already "." (see combine_sanger_sites.py) rather than a real strand
        # on the rare cluster where forward/reverse reads land on the same
        # site but disagree on which way the insertion faces.
        "orientation": "; ".join(sites["strand"]),
        "positions_agree": bool(group["both_directions"].any()),
    }
    if original_site is not None:
        agreeing_sites = sites[sites["both_directions"]]
        result[ORIGINAL_SITE_COLUMN] = bool(
            at_original_site(agreeing_sites, original_site, max_dist).any()
        )
    if region is not None:
        result["in_region"] = bool(in_region(sites, region).any())
    return pd.Series(result)


def summarize_clone_validation(group, original_site=None, max_dist=0):
    """A clone's combined Sanger+NGS verification: "both", "forward only",
    "reverse only", or "not verified" - the same per-row label
    compare_sanger_ngs.py already computes (ngs_verification), which uses
    NGS only to fill in a side Sanger's own primers left unconfirmed rather
    than override a side Sanger already confirmed.

    A clone whose Sanger reads disagree on locus (more than one row here)
    has no single locus to report a combined status for - disagreement
    alone already fails the clone (see tagmaplib.clone_color), regardless
    of what NGS shows. Purely as extra information, not something that can
    rescue it, each disagreeing site is checked against NGS's own side-
    calling there: do the two agree on which single side is present? The
    original, pre-mobilization locus (if configured) is exempted from this
    check - NGS showing both sides there is the ordinary, expected result
    for a real, long-established genomic feature, not a red flag.
    """
    if group.shape[0] == 1:
        return pd.Series({"ngs_verification": group["ngs_verification"].iloc[0]})
    consistent = True
    for row in group.itertuples():
        is_original = (
            original_site is not None
            and row.chrom == original_site["chrom"]
            and abs(row.start - original_site["pos"]) <= max_dist
        )
        if is_original:
            continue
        sanger_sides = set()
        if row.n_forward > 0:
            sanger_sides.add("forward")
        if row.n_reverse > 0:
            sanger_sides.add("reverse")
        ngs_sides = (
            tagmaplib.NGS_SIDE_SETS.get(row.ngs_site_sides, set())
            if row.confirmed_by_ngs
            else set()
        )
        if sanger_sides != ngs_sides:
            consistent = False
    label = "disagreeing (NGS consistent)" if consistent else "disagreeing (NGS inconsistent)"
    return pd.Series({"ngs_verification": label})


READ_QC_COLUMNS = [
    "sample_name",
    "well",
    "direction",
    "status",
    "max_read_length",
    "n_reads",
    "at_original_site",
]
READ_QC_STATUS_RANK = {
    tagmaplib.READ_QC_FAILED: 0,
    tagmaplib.READ_QC_SHORT: 1,
    tagmaplib.READ_QC_WORKED: 2,
}


def summarize_read_qc(reads, original_site=None, max_dist=0):
    """One row per (sample_name, well, direction): the best sequencing
    outcome (tagmaplib.read_qc_status) across every read that landed in that
    physical well - some configs' sanger_name_regex captures a finer-grained
    "clone" than the well it came from (e.g. one well tried with several ITR
    primer pairs), in which case several reads share one well/direction here
    and the best of them wins, since any one of them producing usable
    sequence means that well/side did sequence successfully. Falls back to
    "clone" as the well id when a project's regex has no separate well group
    at all (both then mean the same thing - see sanger_sites.py).

    A well/direction with no reads at all in reads.tsv gets
    READ_QC_NOT_SEQUENCED - distinct from READ_QC_FAILED, which is a trace
    that was actually run but came back under the clean-bases floor (see
    tagmaplib.read_qc_status). A trace that basecalled to literally nothing
    (see ab1_to_fastq.py) never gets a row at all, same as one that was
    never attempted - both count as NOT_SEQUENCED here, since neither left
    anything behind to tell them apart.

    at_original_site is True if any of a well/direction's *passing* reads
    (the ones that actually produced a reported site - see sanger_sites.py)
    lands within max_dist of the configured --original-site - the same
    "still at the founder locus" check clone_summary's own describe()
    applies per clone/both-sides, just per well/direction here instead.
    Always False when original_site is None (not configured).
    """
    if reads.shape[0] == 0 or "well" not in reads.columns:
        return pd.DataFrame(columns=READ_QC_COLUMNS)
    reads = reads.copy()
    reads["well"] = reads["well"].where(reads["well"].notna(), reads["clone"])
    reads = reads.dropna(subset=["well"])
    if reads.shape[0] == 0:
        return pd.DataFrame(columns=READ_QC_COLUMNS)
    reads["read_length"] = pd.to_numeric(reads["read_length"], errors="coerce")
    reads["status"] = reads["read_length"].map(tagmaplib.read_qc_status)
    if original_site is not None:
        reads["at_original_site"] = reads["pass"].astype(bool) & at_original_site(
            reads, original_site, max_dist
        )
    else:
        reads["at_original_site"] = False

    rows = []
    for sample_name, sample_reads in reads.groupby("sample_name"):
        wells = sorted(sample_reads["well"].astype(str).unique())
        for direction in ("forward", "reverse"):
            by_well = dict(
                tuple(sample_reads[sample_reads["direction"] == direction].groupby("well"))
            )
            for well in wells:
                group = by_well.get(well)
                if group is None:
                    status, max_length, n_reads, at_site = (
                        tagmaplib.READ_QC_NOT_SEQUENCED,
                        pd.NA,
                        0,
                        False,
                    )
                else:
                    status = max(group["status"], key=READ_QC_STATUS_RANK.get)
                    max_length = group["read_length"].max()
                    n_reads = group.shape[0]
                    at_site = bool(group["at_original_site"].any())
                rows.append(
                    {
                        "sample_name": sample_name,
                        "well": well,
                        "direction": direction,
                        "status": status,
                        "max_read_length": max_length,
                        "n_reads": n_reads,
                        "at_original_site": at_site,
                    }
                )
    result = pd.DataFrame(rows, columns=READ_QC_COLUMNS)
    result["n_reads"] = result["n_reads"].astype(int)
    return result


POSITION_BASE_COLUMNS = ["sample_name", "n_distinct_positions", "n_one_sided"]
POSITION_NGS_COLUMN = "n_ngs_confirmed"


def position_columns(has_validation):
    return POSITION_BASE_COLUMNS + ([POSITION_NGS_COLUMN] if has_validation else [])


def position_flags(site_clusters, group_cols):
    """One row per distinct position (group_cols is ["sample_name", "chrom",
    "start"] for the per-sample counts below, or just ["chrom", "start"] for
    the "all" total): one_sided is True if no clone confirming that position
    did so from both ITR primers; ngs_confirmed (only present when NGS
    validation data was merged into site_clusters) is True if any clone
    there was independently confirmed by NGS.
    """
    flags = (
        site_clusters.groupby(group_cols)["both_directions"]
        .any()
        .reset_index(name="any_both_directions")
    )
    flags["one_sided"] = ~flags["any_both_directions"]
    if "confirmed_by_ngs" in site_clusters.columns:
        ngs = (
            site_clusters.groupby(group_cols)["confirmed_by_ngs"]
            .any()
            .reset_index(name="ngs_confirmed")
        )
        flags = flags.merge(ngs, on=group_cols)
    return flags


def summarize_positions(site_clusters):
    """Distinct (chrom, start) positions per plate, plus an "all" totals row.

    Counted from the already-clustered sites rather than per-read, so a clone
    seen from both primers only contributes its one agreed-on position - and
    two clones that land on the very same locus (e.g. an unmobilized founder
    control sequenced more than once) count as one position, not two.

    n_one_sided counts how many of those positions have no clone confirming
    them from both ITR primers at all - every clone that landed there did so
    from a single primer, so the position rests on weaker evidence than a
    clone (or two agreeing clones) confirmed from both sides would give it.
    Where NGS data for the same material is available, n_ngs_confirmed
    separately counts how many positions were independently backed by NGS -
    reported whether or not the position is also one-sided, since NGS
    confirmation is informative either way, not just as a rescue for weak
    Sanger evidence.
    """
    has_validation = "confirmed_by_ngs" in site_clusters.columns
    columns = position_columns(has_validation)
    if site_clusters.shape[0] == 0:
        return pd.DataFrame(columns=columns)

    per_sample_flags = position_flags(site_clusters, ["sample_name", "chrom", "start"])
    per_sample = (
        per_sample_flags.groupby("sample_name")
        .agg(
            n_distinct_positions=("one_sided", "size"),
            n_one_sided=("one_sided", "sum"),
        )
        .reset_index()
    )
    if has_validation:
        ngs_per_sample = (
            per_sample_flags.groupby("sample_name")["ngs_confirmed"]
            .sum()
            .reset_index(name=POSITION_NGS_COLUMN)
        )
        per_sample = per_sample.merge(ngs_per_sample, on="sample_name")
    per_sample = per_sample.sort_values("sample_name")

    total_flags = position_flags(site_clusters, ["chrom", "start"])
    total = {
        "sample_name": "all",
        "n_distinct_positions": total_flags.shape[0],
        "n_one_sided": int(total_flags["one_sided"].sum()),
    }
    if has_validation:
        total[POSITION_NGS_COLUMN] = int(total_flags["ngs_confirmed"].sum())
    total = pd.DataFrame([total])

    result = pd.concat([per_sample, total], ignore_index=True)[columns]
    result["n_one_sided"] = result["n_one_sided"].astype(int)
    if has_validation:
        result[POSITION_NGS_COLUMN] = result[POSITION_NGS_COLUMN].astype(int)
    return result


POSITION_COUNT_BASE_COLUMNS = ["chrom", "start", "end", "n_times_found", "n_one_sided"]
NGS_VERIFICATION_COLUMN = "ngs_verification"
DISAGREEING_CLONES_LABEL = "invalid: forward/reverse reads disagree on position"
UNMOBILIZED_LABEL = "unmobilized: at the original insertion site"
TOTAL_CLONES_LABEL = "clones with at least one good read"
TOTAL_MOBILIZED_POSITIONS_LABEL = "total mobilized positions"


def locus_ngs_verification(group):
    """A locus's combined Sanger+NGS verification, pooling every clone that
    landed there (see position_counts) rather than one clone at a time -
    the same "both"/"forward only"/"reverse only"/"not verified" label as
    compare_sanger_ngs.py's per-row ngs_verification, just evaluated across
    the whole group's Sanger reads and NGS matches instead of one row's.
    """
    sanger_sides = []
    if (group["n_forward"] > 0).any():
        sanger_sides.append("forward")
    if (group["n_reverse"] > 0).any():
        sanger_sides.append("reverse")
    ngs_sides = []
    for sides_value in group.loc[group["confirmed_by_ngs"] == True, "ngs_site_sides"]:
        ngs_sides.extend(tagmaplib.NGS_SIDE_SETS.get(sides_value, []))
    return tagmaplib.ngs_verification_label(sanger_sides, ngs_sides)


def position_counts(site_clusters, original_site=None, max_dist=0):
    """How many clones land on each distinct (chrom, start, end).

    One count per clone, not per site row: a clone whose forward and reverse
    reads cluster into two different loci (see combine_sanger_sites.py) -
    e.g. one primer stuck on a residual copy of the donor construct while
    the other reaches the real integration site - doesn't confirm either
    locus, so counting it at both (or picking one arbitrarily) would
    misrepresent what the clone actually shows. Such clones are tallied
    separately instead, as one DISAGREEING_CLONES_LABEL row - never split
    across the several loci they touch.

    A clone whose one agreed-on site (already disagreement-free, or there
    would be no single site to speak of) is within max_dist of
    original_site is not a new integration either. Rather than let
    founder-locus clones fragment across several near-identical coordinates
    (mapping noise means they rarely land on the exact same base), they are
    pulled out of the per-locus breakdown and tallied together as a single
    UNMOBILIZED_LABEL row.

    Both that row and every real position row count clones confirmed by
    only one primer (the other side was never sequenced, or was sequenced
    but failed QC) the same as clones both primers agree on - counting one
    is not counting the other twice - but n_one_sided says how many of a
    row's clones that was, since one-sided evidence is weaker: the
    DISAGREEING_CLONES_LABEL clones above show that a primer stuck on the
    founder locus while the other primer, when it is usable at all, often
    turns out to disagree and land somewhere else - so a one-sided clone
    parked on the founder locus could equally be a mobilized clone whose
    other primer simply failed rather than a true non-mover. n_one_sided is
    NA for DISAGREEING_CLONES_LABEL (both primers there were confirmed, on
    different sites - the opposite situation) and for the two totals below.

    Two totals close out the table: TOTAL_CLONES_LABEL (every clone with at
    least one QC-passing, motif-snapped site - i.e. every clone represented
    somewhere above, whether in a real position row or one of the two label
    rows) and TOTAL_MOBILIZED_POSITIONS_LABEL (how many distinct new-
    integration loci that clone pool actually produced - the count of real
    position rows, excluding both label rows).

    When NGS validation is available (site_clusters carries confirmed_by_ngs
    and ngs_site_sides - see compare_sanger_ngs.py), an NGS_VERIFICATION_COLUMN
    is added: the combined Sanger+NGS verification for every clone pooled
    into that row (see locus_ngs_verification) - "both", "forward only",
    "reverse only", or "not verified". NA for DISAGREEING_CLONES_LABEL (both
    primers there were confirmed, just on different sites - not a case of
    one side being unverified) and for the two totals, neither of which is
    one locus.

    Sorted most-found first, with the label/total rows always last.
    """
    has_validation = "confirmed_by_ngs" in site_clusters.columns
    columns = POSITION_COUNT_BASE_COLUMNS + (
        [NGS_VERIFICATION_COLUMN] if has_validation else []
    )
    if site_clusters.shape[0] == 0:
        return pd.DataFrame(columns=columns)
    clone_sizes = site_clusters.groupby(["sample_name", "clone"]).size()
    n_total_clones = int(clone_sizes.shape[0])
    n_disagreeing = int((clone_sizes > 1).sum())
    agreeing = site_clusters.merge(
        clone_sizes.reset_index(name="n_sites"), on=["sample_name", "clone"]
    )
    agreeing = agreeing[agreeing["n_sites"] == 1].copy()
    agreeing["one_sided"] = ~agreeing["both_directions"]

    n_unmobilized = 0
    n_unmobilized_one_sided = 0
    unmobilized_ngs_verification = pd.NA
    if original_site is not None and agreeing.shape[0]:
        unmobilized = at_original_site(agreeing, original_site, max_dist)
        n_unmobilized = int(unmobilized.sum())
        n_unmobilized_one_sided = int((unmobilized & agreeing["one_sided"]).sum())
        if has_validation and n_unmobilized:
            unmobilized_ngs_verification = locus_ngs_verification(
                agreeing[unmobilized]
            )
        agreeing = agreeing[~unmobilized]

    counts = (
        agreeing.groupby(["chrom", "start", "end"])
        .agg(n_times_found=("clone", "count"), n_one_sided=("one_sided", "sum"))
        .reset_index()
        .sort_values(
            ["n_times_found", "chrom", "start"], ascending=[False, True, True]
        )
        .reset_index(drop=True)
    )
    counts["start"] = counts["start"].astype("Int64")
    counts["end"] = counts["end"].astype("Int64")
    counts["n_one_sided"] = counts["n_one_sided"].astype("Int64")
    n_mobilized_positions = int(counts.shape[0])
    if has_validation:
        if counts.shape[0]:
            ngs_col = (
                agreeing.groupby(["chrom", "start", "end"])
                .apply(locus_ngs_verification, include_groups=False)
                .reset_index(name=NGS_VERIFICATION_COLUMN)
            )
            counts = counts.merge(ngs_col, on=["chrom", "start", "end"], how="left")
        else:
            counts[NGS_VERIFICATION_COLUMN] = pd.Series(dtype=object)

    rows = []
    if n_unmobilized:
        rows.append(
            (UNMOBILIZED_LABEL, n_unmobilized, n_unmobilized_one_sided, unmobilized_ngs_verification)
        )
    if n_disagreeing:
        rows.append((DISAGREEING_CLONES_LABEL, n_disagreeing, pd.NA, pd.NA))
    rows.append((TOTAL_CLONES_LABEL, n_total_clones, pd.NA, pd.NA))
    rows.append((TOTAL_MOBILIZED_POSITIONS_LABEL, n_mobilized_positions, pd.NA, pd.NA))
    extra = pd.DataFrame(
        {
            "chrom": [row[0] for row in rows],
            "start": pd.array([pd.NA] * len(rows), dtype="Int64"),
            "end": pd.array([pd.NA] * len(rows), dtype="Int64"),
            "n_times_found": [row[1] for row in rows],
            "n_one_sided": pd.array([row[2] for row in rows], dtype="Int64"),
        }
    )
    if has_validation:
        extra[NGS_VERIFICATION_COLUMN] = [row[3] for row in rows]
    counts = pd.concat([counts, extra], ignore_index=True)
    return counts[columns]


if not args.reads:
    qc = pd.DataFrame(columns=QC_COLUMNS)
    clone_summary = pd.DataFrame(columns=BASE_CLONE_COLUMNS)
    read_qc = pd.DataFrame(columns=READ_QC_COLUMNS)
else:
    reads = pd.concat(
        (pd.read_csv(path, sep="\t", dtype={"chrom": str}) for path in args.reads),
        ignore_index=True,
    )
    reads["pass"] = tagmaplib.to_bool(reads["pass"])
    read_qc = summarize_read_qc(
        reads, original_site=original_site, max_dist=args.original_max_dist
    )

    if PLATE_RUN_COLUMN in reads.columns:
        plate_runs = (
            reads.groupby("sample_name")[PLATE_RUN_COLUMN]
            .agg(lambda s: "/".join(sorted(s.dropna().astype(str).unique())))
            .reset_index()
        )

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

        # clone_sides is the complete (sample_name, direction) grid - every
        # sample gets a forward and a reverse row regardless of whether any
        # read of that direction was ever attempted - so it, not qc (which
        # only has rows for directions reads.tsv actually has), anchors this
        # merge: a direction nobody ever sequenced still shows n_clones_pass=0
        # of the sample's full clone count, rather than silently vanishing.
        clone_sides = summarize_clone_sides(clone_summary)
        clone_sides["frac_clones_pass"] = (
            clone_sides["n_clones_pass"] / clone_sides["n_clones"]
        )
        qc = clone_sides.merge(qc, on=["sample_name", "direction"], how="left")
        qc[["n_reads", "n_pass", "n_fail"]] = (
            qc[["n_reads", "n_pass", "n_fail"]].fillna(0).astype(int)
        )
        qc["fail_reasons"] = qc["fail_reasons"].fillna("")
        qc = qc[QC_COLUMNS].sort_values(["sample_name", "direction"])

# Read once, up front, so both position_counts (pooled across clones sharing
# a locus) and the per-clone summary further down can use the same NGS
# verification without re-reading the file.
validation = None
if args.validation is not None:
    validation = pd.read_csv(args.validation, sep="\t", dtype={"chrom": str})
    if validation.shape[0] and "confirmed_by_ngs" in validation.columns:
        validation["confirmed_by_ngs"] = tagmaplib.to_bool(validation["confirmed_by_ngs"])
    else:
        validation = None

# Both primers can independently pass QC while pointing at different loci, so
# "confirmed" from each side isn't enough on its own - cross-check against the
# already-clustered sites, which group a clone's passing reads by position.
if args.sites is not None:
    site_clusters = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str})
    if site_clusters.shape[0] and "both_directions" in site_clusters.columns:
        site_clusters["both_directions"] = tagmaplib.to_bool(
            site_clusters["both_directions"]
        )
    if validation is not None:
        site_clusters = site_clusters.merge(
            validation[
                ["sample_name", "clone", "chrom", "start", "end"]
                + ["confirmed_by_ngs", "ngs_site_sides"]
            ],
            on=["sample_name", "clone", "chrom", "start", "end"],
            how="left",
        )
    positions = summarize_positions(site_clusters)
    counts = position_counts(
        site_clusters, original_site=original_site, max_dist=args.original_max_dist
    )
    if site_clusters.shape[0] and "both_directions" in site_clusters.columns:
        site_summary = (
            site_clusters.groupby(["sample_name", "clone"])
            .apply(
                summarize_clone_sites,
                original_site=original_site,
                max_dist=args.original_max_dist,
                region=region,
                include_groups=False,
            )
            .reset_index()
        )
    else:
        site_summary = pd.DataFrame(
            columns=["sample_name", "clone", "position", "orientation", "positions_agree"]
            + ([ORIGINAL_SITE_COLUMN] if original_site is not None else [])
            + (REGION_COLUMNS if region is not None else [])
        )
    clone_summary = clone_summary.merge(
        site_summary, on=["sample_name", "clone"], how="left"
    )
    clone_summary["position"] = clone_summary["position"].fillna("")
    clone_summary["orientation"] = clone_summary["orientation"].fillna("")
    clone_summary["positions_agree"] = (
        clone_summary["positions_agree"].fillna(False).astype(bool)
    )
    if original_site is not None:
        clone_summary[ORIGINAL_SITE_COLUMN] = (
            clone_summary[ORIGINAL_SITE_COLUMN].fillna(False).astype(bool)
        )
    if region is not None:
        clone_summary["in_region"] = clone_summary["in_region"].fillna(False).astype(bool)
else:
    clone_summary["position"] = ""
    clone_summary["orientation"] = ""
    clone_summary["positions_agree"] = False
    if original_site is not None:
        clone_summary[ORIGINAL_SITE_COLUMN] = False
    if region is not None:
        clone_summary["in_region"] = False
    positions = pd.DataFrame(columns=POSITION_COLUMNS)
    counts = pd.DataFrame(columns=POSITION_COUNT_BASE_COLUMNS)

clone_summary["both_sides_confirmed"] = (
    (clone_summary["forward_status"] == "PASSED")
    & (clone_summary["reverse_status"] == "PASSED")
    & clone_summary["positions_agree"]
)
unmobilized_values = (
    clone_summary[ORIGINAL_SITE_COLUMN]
    if original_site is not None
    else [False] * clone_summary.shape[0]
)
clone_summary["summary"] = [
    describe(f, r, p, u)
    for f, r, p, u in zip(
        clone_summary["forward_status"],
        clone_summary["reverse_status"],
        clone_summary["positions_agree"],
        unmobilized_values,
    )
]
output_columns = CLONE_COLUMNS + (REGION_COLUMNS if region is not None else [])
clone_summary = clone_summary[output_columns]

if validation is not None:
    clone_validation = (
        validation.groupby(["sample_name", "clone"])
        .apply(
            summarize_clone_validation,
            original_site=original_site,
            max_dist=args.original_max_dist,
            include_groups=False,
        )
        .reset_index()
    )
    clone_summary = clone_summary.merge(
        clone_validation, on=["sample_name", "clone"], how="left"
    )
    clone_summary["ngs_verification"] = clone_summary["ngs_verification"].fillna(
        "not verified"
    )
    clone_summary = clone_summary[output_columns + CLONE_VALIDATION_COLUMNS]

if plate_runs is not None:
    clone_summary = clone_summary.merge(plate_runs, on="sample_name", how="left")
    clone_summary[PLATE_RUN_COLUMN] = clone_summary[PLATE_RUN_COLUMN].fillna("")

qc.to_csv(args.output_qc, sep="\t", index=False)
clone_summary.to_csv(args.output_clone_summary, sep="\t", index=False)
positions.to_csv(args.output_positions, sep="\t", index=False)
counts.to_csv(args.output_position_counts, sep="\t", index=False)
read_qc.to_csv(args.output_read_qc, sep="\t", index=False)

print(
    f"{clone_summary.shape[0]} clones: "
    f"{int(clone_summary['both_sides_confirmed'].sum()) if clone_summary.shape[0] else 0} "
    "confirmed from both primers"
)
