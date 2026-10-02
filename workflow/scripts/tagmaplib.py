"""Helpers shared between the NGS and Sanger branches of the workflow.

The scripts in this folder are meant to stay runnable by hand, so they keep
their own argparse interfaces and simply ``import tagmaplib`` - Python puts a
script's own directory on ``sys.path``, so no installation is needed.
"""

import re

import numpy as np
import pandas as pd

# Columns of the per-sample peak files written by combine_peaks.py and
# concatenated into all_peaks.bed.
PEAK_COLUMNS = [
    "chrom",
    "start",
    "end",
    "sample_name",
    "count",
    "side",
    "fraction",
    "n_positions",
    "orientation",
]

# Columns of the final insertion site files, shared by the NGS and Sanger
# branches so that the two can be compared and loaded the same way.
SITE_COLUMNS = ["chrom", "start", "end", "sample_name", "score", "strand"]

# Columns of the bedgraph files written by coverage.py.
COVERAGE_COLUMNS = ["chrom", "start", "end", "count", "fraction", "orientation"]

# Column groups of ngs_qc_stats.tsv (written by ngs_stats.py) - split out
# here, rather than left as one flat column list, so report.py/report_pdf.py
# can render the raw pre-deduplication counts and the deduplicated/sidedness
# counts as two separate tables: together they're too wide for one table's
# columns to stay readable even on a landscape page.
NGS_QC_RAW_COLUMNS = [
    "sample_name",
    "total_pairs",
    "mapped_pairs",
    "frac_mapped",
    "raw_unusable_pairs",
    "raw_itr_anchored_pairs",
    "raw_forward_itr_anchored_pairs",
    "raw_reverse_itr_anchored_pairs",
    "raw_junction_pairs",
    "raw_forward_junction_pairs",
    "raw_reverse_junction_pairs",
    "raw_non_junction_pairs",
    "raw_peak_support_pairs",
    "raw_site_support_pairs",
]

NGS_QC_DEDUPLICATED_COLUMNS = [
    "forward_junction_pairs",
    "reverse_junction_pairs",
    "mobilized_pairs",
    "total_deduplicated_junction_pairs",
    "frac_mobilized",
]

NGS_QC_SIDEDNESS_COLUMNS = [
    "sample_name",
    "n_sites",
    "n_two_sided",
    "frac_two_sided",
    "n_one_sided",
    "n_forward_only",
    "n_reverse_only",
]


# Walk types of pairtools parse2 in which the cassette/genome junction was
# actually sequenced, because both sides of the pair come from within one read
# (or from reads that overlap). "R1-2" is the other case: the two sides come
# from the two different reads, so whatever lies between them was never read.
JUNCTION_WALK_TYPES = ("R1", "R2", "R1&2")


def cassette_orientation(side, runs_rightwards):
    """Orientation of an integration, from which way its reads leave it.

    Reads run outwards from the cassette, so which end of it a read came from
    plus which way along the genome it then ran fixes how the cassette sits -
    no need to compare the positions of the two sides, which say nothing once
    both are pinned to the same base. The flip is the pipeline's convention,
    set by combine_sanger_sites.py flipping its reverse reads (forward is the
    reference, left unflipped), so that both branches call the same
    insertion the same way.
    """
    rightwards = np.asarray(runs_rightwards, dtype=bool)
    from_reverse = np.asarray(side) == "reverse"
    return np.where(from_reverse == rightwards, "-", "+")


def normalise_junction_pos(pos, strand):
    """Make both strands name the same base as the junction.

    pairtools reports the 5' end of a plus-strand alignment one base further
    right than that of a minus-strand one, so a single insertion read from
    either direction comes out as two positions 1bp apart. Left in, that splits
    every site in two, and makes the order of the two ITR sides - which is all
    find_insertion_sites has to go on for the strand - a coin flip.
    """
    return np.asarray(pos) - (np.asarray(strand) == "+")


def read_pairs(pairs, threads=1):
    """Read a .pairs file into a frame, plus the chromsizes from its header."""
    from pairtools.lib import fileio, headerops

    pairs_stream = (
        fileio.auto_open(pairs, mode="r", nproc=threads)
        if isinstance(pairs, str)
        else pairs
    )
    header, pairs_body = headerops.get_header(pairs_stream)
    cols = headerops.extract_column_names(header)
    chromsizes = headerops.extract_chromsizes(header)
    pairs_df = pd.read_csv(
        pairs_body,
        header=None,
        names=cols,
        chunksize=None,
        sep="\t",
        dtype={"chrom1": str, "chrom2": str},
    )
    return pairs_df, chromsizes


def junction_position(start, end, strand):
    """0-based genomic coordinate of an alignment's cassette/genome junction.

    Reads are sequenced outwards from the cassette, so the junction is at the
    alignment's 5' end in read orientation: the start of a forward alignment,
    the last base of a reverse one. This is the definition of an integration
    site for both branches of the workflow. Works on scalars and on arrays.
    """
    pos = np.where(np.asarray(strand) == "+", np.asarray(start), np.asarray(end) - 1)
    return pos.item() if pos.ndim == 0 else pos


def flip_strand(strand):
    """Swap + and -, leaving anything else (e.g. '.') alone."""
    strand = pd.Series(strand, dtype=object)
    return strand.map({"+": "-", "-": "+"}).fillna(strand)


def sample_and_side(path, suffix):
    """Pull {sample}/{side} back out of a "{sample}_{side}{suffix}" filename -
    the naming convention find_peaks_junctions.py's per-side outputs and
    coverage.py's bedgraphs both use.
    """
    import os

    base = os.path.basename(path)[: -len(suffix)]
    sample, side = base.rsplit("_", 1)
    return sample, side


def read_peaks(path, columns=None):
    """Read a headerless peak/site BED, tolerating an empty file."""
    columns = PEAK_COLUMNS if columns is None else columns
    try:
        df = pd.read_csv(
            path, sep="\t", header=None, names=columns, dtype={"chrom": str}
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in columns})
    return df


def apply_blacklist(df, path):
    """Drop intervals overlapping a blacklist BED. No-op if path is None."""
    import bioframe

    if path is None:
        return df
    blacklist = pd.read_csv(
        path, sep="\t", header=None, comment="#", names=["chrom", "start", "end"]
    )
    return bioframe.setdiff(df, blacklist)


def find_insertion_seq(df, genome, ins_seq, window=0, mode="first", index_file=None):
    """Locate the insertion motif around each interval.

    Transposons integrate into a fixed short motif (TA for Sleeping Beauty,
    TTAA for PiggyBac), so the exact integration site can be recovered by
    searching around the mapped junction. Uses pysam's faidx-based random
    access rather than streaming the whole genome.

    Searches ``[start - window, end + window)`` and returns the 0-based genomic
    start of the chosen occurrence per row, or -1 where there is none. ``mode``
    is ``"first"`` for the leftmost occurrence, or ``"nearest"`` for the one
    closest to ``start`` - the latter suits single-junction intervals, where
    the motif should be right at the mapped position. Rows on strand '.' are
    skipped, since without a strand there is no junction to refine.

    Callers decide what to do with the result, because the NGS and Sanger
    branches report insertion sites with different widths.
    """
    import pysam

    ins_seq = ins_seq.upper()
    result = pd.Series(-1, index=df.index, dtype=int)
    if df.shape[0] == 0:
        return result

    # pyfastx's random-access reader segfaults on genomes in the multi-GB
    # range (reproduces even without calling .fetch() at all, so it is not
    # data-dependent) - pysam.FastaFile is the well-exercised equivalent.
    fasta = pysam.FastaFile(genome, filepath_index=index_file)
    chrom_lengths = dict(zip(fasta.references, fasta.lengths))

    for i, row in df.iterrows():
        if row.get("strand", ".") == "." or row["chrom"] not in chrom_lengths:
            continue
        lo = max(0, int(row["start"]) - window)
        hi = min(chrom_lengths[row["chrom"]], int(row["end"]) + window)
        if hi - lo < len(ins_seq):
            continue
        # pysam.FastaFile.fetch takes a 0-based, half-open interval.
        seq = fasta.fetch(row["chrom"], lo, hi).upper()
        hits = []
        j = seq.find(ins_seq)
        while j != -1:
            hits.append(lo + j)
            if mode == "first":
                break
            j = seq.find(ins_seq, j + 1)
        if not hits:
            continue
        if mode == "first":
            result[i] = hits[0]
        else:
            # Ties go leftwards, so the result does not depend on search order.
            result[i] = min(hits, key=lambda p: (abs(p - int(row["start"])), p))
    return result


def insertion_site_from_motif(motif_start, ins_seq):
    """The transposon's actual cut site, given the motif's own start.

    Target-site duplication is symmetric around the middle of the motif (TA
    for Sleeping Beauty, TTAA for PiggyBac) - the transposon inserts there,
    and the flanking copy grows to either side as a result. In 0-based
    half-open coordinates that midpoint has an exact integer representation
    (no half-base rounding needed): motif_start + len(ins_seq) // 2 is the
    coordinate of the gap between the motif's two halves, e.g. for TA at
    genomic position p (p is the T, p+1 is the A), the insertion site is
    p + 1 - the same coordinate as the A's own 0-based position.
    """
    return motif_start + len(ins_seq) // 2


NGS_SIDE_SETS = {
    "both": {"forward", "reverse"},
    "forward_only": {"forward"},
    "reverse_only": {"reverse"},
}

# all_peaks.bed's own "side" column (see combine_peaks.py) encodes which ITR
# primer a peak is anchored at as "+"/"-", the same symbols - but not the same
# meaning - as a site's genomic strand. Shared so compare_sanger_ngs.py and
# sanger_stats.py's exact-peak-match fallback (see ngs_verification_label)
# read it the same way.
PEAK_SIDE_DIRECTIONS = {"+": "forward", "-": "reverse"}


# Which sides each ngs_verification_label result stands for; the labels that
# say nothing about sides are simply absent.
NGS_LABEL_SIDES = {
    "both": {"forward", "reverse"},
    "forward only": {"forward"},
    "reverse only": {"reverse"},
}


def ngs_verification_label(ngs_sides, ngs_checked=True, ambiguous=False):
    """What NGS alone says about a Sanger site: "both" / "forward only" /
    "reverse only" / "ambiguous" / "no NGS site nearby" / "not verified".

    Deliberately says nothing about which primer(s) Sanger itself confirmed -
    that is already in the Sanger columns (n_forward/n_reverse,
    forward_status/reverse_status). ngs_sides is an iterable of
    "forward"/"reverse" (or empty) - see compare_sanger_ngs.py (per Sanger
    site) and sanger_stats.py's position_counts (pooled across every clone
    sharing a locus). It covers two tiers the caller has already merged: a
    side from a fully called, orientation-resolved NGS *site* nearby, or -
    lacking that - a side from a single NGS *peak* landing at the exact same
    base (see PEAK_SIDE_DIRECTIONS), weaker evidence than a called site but
    still an exact positional match rather than "somewhere nearby".

    ngs_checked distinguishes "NGS was matched to a site or exact-position
    peak here but gave no usable side (e.g. its strand contradicts Sanger's,
    so "not verified")" from "there was nothing from NGS within the matching
    distance at all" ("no NGS site nearby").

    ambiguous says the NGS sites around the Sanger position disagree with
    each other (orientation or which sides they show - see
    compare_sanger_ngs.py), so which of them the Sanger read really belongs
    to is unknown and none of their sides can be attributed to it.
    """
    if ambiguous:
        return "ambiguous"
    ngs_sides = set(ngs_sides)
    if ngs_sides >= NGS_LABEL_SIDES["both"]:
        return "both"
    if ngs_sides == {"forward"}:
        return "forward only"
    if ngs_sides == {"reverse"}:
        return "reverse only"
    return "not verified" if ngs_checked else "no NGS site nearby"


# UCSC's colorByStrand track attribute: '+' strand features render red,
# '-' strand render blue, without needing a per-row itemRgb column.
STRAND_COLORS = "255,0,0 0,0,255"


def write_bed(
    df, path, columns=None, track_name=None, color_by_strand=False, header=False
):
    """Write a BED/bedGraph, optionally with a UCSC track line.

    color_by_strand adds colorByStrand=STRAND_COLORS to the track line, so
    +/- features render in different colors in the browser. Only meaningful
    alongside track_name, for a BED with a strand column.
    """
    columns = list(df.columns) if columns is None else columns
    mode = "w"
    if track_name is not None:
        track_line = f"track name={track_name}"
        if color_by_strand:
            track_line += f" colorByStrand='{STRAND_COLORS}'"
        with open(path, "w") as f:
            f.write(track_line + "\n")
        mode = "a"
    df[columns].to_csv(path, sep="\t", index=False, header=header, mode=mode)


def to_bool(series):
    """Parse a True/False column that survived a round trip through a TSV.

    to_csv writes booleans as the strings 'True'/'False' and pandas.NA as an
    empty field, none of which astype(bool) handles correctly on reload
    (every non-empty string, including 'False', is truthy). Unresolved/empty
    values are treated as False, since here they mean "nothing to confirm
    it" rather than an unknown true value.
    """
    return (
        series.map({True: True, False: False, "True": True, "False": False})
        .fillna(False)
        .astype(bool)
    )


# Shared between report.py's markdown table and report_pdf.py's row fill, so
# a clone's status maps to the same colour in both.
CLONE_GREEN = (198, 239, 206)
CLONE_YELLOW = (255, 235, 156)
CLONE_RED = (255, 199, 206)


def clone_color(
    forward_status, reverse_status, both_sides_confirmed, ngs_verification=None
):
    """Which of CLONE_GREEN/YELLOW/RED a clone's status maps to.

    - green: both sides confirmed and agree on position; or only one side
      passed and NGS shows both sides at that locus (see
      ngs_verification_label) - regardless of whether the other primer was
      never sequenced or sequenced and failed, since neither case is
      positive confirmation of the missing side on its own.
    - yellow: only one side passed, and NGS did not confirm both sides
      there (no NGS data, nothing found, ambiguous, or NGS itself only
      one-sided on the same side Sanger saw).
    - red: neither side passed, or both passed individually but disagree on
      position - positive evidence of a conflict, which NGS cannot rescue
      (see ngs_verification_label's "disagreeing" clones - checking NGS
      there is informational only, not a path back to green); or only one
      side passed and NGS shows only the opposite side there, so the side
      Sanger saw is absent from NGS at its own site - also a conflict.
    """
    if both_sides_confirmed:
        return CLONE_GREEN
    forward_pass = forward_status == "PASSED"
    reverse_pass = reverse_status == "PASSED"
    if forward_pass and reverse_pass:
        return CLONE_RED  # both passed, but both_sides_confirmed is False: positions disagree
    if forward_pass or reverse_pass:
        passed_side = "forward" if forward_pass else "reverse"
        ngs_sides = NGS_LABEL_SIDES.get(ngs_verification, set())
        if ngs_sides and passed_side not in ngs_sides:
            return CLONE_RED  # NGS sees only the side Sanger did not
        return CLONE_GREEN if ngs_sides >= NGS_LABEL_SIDES["both"] else CLONE_YELLOW
    return CLONE_RED


# Thresholds for read_qc_status, on a read's length *after* quality trimming
# (ab1_to_fastq.py's Mott trimming/tracy - i.e. its count of clean, usable
# bases, not the raw trace length).
READ_QC_MIN_CLEAN_BASES = 20
READ_QC_SHORT_MAX_LENGTH = 100

READ_QC_WORKED = "worked"
READ_QC_SHORT = "very short"
READ_QC_FAILED = "failed"
# Distinct from READ_QC_FAILED: there is no trace file for this well/side at
# all - see summarize_read_qc in sanger_stats.py, which is the only place
# this is ever assigned (a real, existing read is always one of the three
# statuses above, never this one).
READ_QC_NOT_SEQUENCED = "not sequenced"


def read_qc_status(read_length):
    """READ_QC_WORKED/SHORT/FAILED for one *existing* read's clean (trimmed)
    length - i.e. a trace that was actually run, whether or not it produced
    much (or any) usable sequence. read_length is NaN if a trace basecalled
    to literally nothing (see ab1_to_fastq.py); that's treated the same as a
    read that came back a few garbage bases, not the same as a well/side
    with no trace file at all (READ_QC_NOT_SEQUENCED, which this function is
    never called for). Between the clean-bases floor and
    READ_QC_SHORT_MAX_LENGTH there is real sequence, just not much of it,
    hence the separate SHORT bucket.
    """
    if pd.isna(read_length) or read_length < READ_QC_MIN_CLEAN_BASES:
        return READ_QC_FAILED
    if read_length < READ_QC_SHORT_MAX_LENGTH:
        return READ_QC_SHORT
    return READ_QC_WORKED


# Shared between report.py's markdown grid and report_pdf.py's drawn grid, so
# a well's read QC status maps to the same colour (and the same colours as
# clone_color's, plus grey for "no trace file at all") in both.
READ_QC_GREY = (220, 220, 220)
READ_QC_COLORS = {
    READ_QC_WORKED: CLONE_GREEN,
    READ_QC_SHORT: CLONE_YELLOW,
    READ_QC_FAILED: CLONE_RED,
    READ_QC_NOT_SEQUENCED: READ_QC_GREY,
}


# ANSI/SLAS 1-2004: every SBS-format microplate (96, 384, ...) shares a 9mm
# well-to-well pitch, and a 96-well plate's round wells are ~6.4mm across -
# shared so report.py's HTML grid and report_pdf.py's drawn grid both come
# out to scale with a real plate, not just a same-shaped table.
WELL_PITCH_MM = 9
WELL_DIAMETER_MM = 6.4

WELL_RE = re.compile(r"^(?:\d+)?([A-Za-z])(\d{1,2})$")


def parse_well(clone):
    """(row, col), both 0-based, if `clone` looks like a plate well (e.g.
    "A01"/"A1"), optionally with a leading plate number as rename_clone_id
    adds (e.g. "1A01"), else None.

    Only meaningful for projects where a Sanger clone id is the well it was
    picked from - not universal, so callers should treat a project where no
    clone parses as "no plate to draw" rather than an error.
    """
    match = WELL_RE.match(str(clone).strip())
    if match is None:
        return None
    row, col = match.groups()
    return ord(row.upper()) - ord("A"), int(col) - 1


CLONE_ID_RE = re.compile(r"^(?:(batch2)_)?plate(\d+)_(.+)$")


def rename_clone_id(sample_name, clone):
    """"plateN_POOL" (or "batch2_plateN_POOL") sample name and its raw
    well/clone id, in this lab's own shorthand - "POOLplN" (or "b2POOLplN"
    for batch 2) for the sample, and the plate number folded onto the front
    of the well/clone id (e.g. "plate1_R"/"A02" -> "Rpl1"/"1A02",
    "batch2_plate1_R"/"A02" -> "b2Rpl1"/"1A02") so a clone id alone still
    says which plate it came from once the sample name is dropped (e.g. in a
    dedup list's "(+N others)" name, or a well grid's own label). Left
    untouched if sample_name isn't in either form this project's sample
    sheets use - e.g. NGS sample names, which never are.
    """
    match = CLONE_ID_RE.match(str(sample_name))
    if match is None:
        return sample_name, clone
    batch, plate_num, pool = match.groups()
    prefix = "b2" if batch else ""
    return f"{prefix}{pool}pl{plate_num}", f"{plate_num}{clone}"


def plate_layout(positions):
    """(n_rows, n_cols) of the smallest standard plate - at least 96-well,
    8x12 - that fits every (row, col) in `positions`."""
    positions = list(positions)
    if not positions:
        return 8, 12
    max_row = max(row for row, _ in positions)
    max_col = max(col for _, col in positions)
    return max(8, max_row + 1), max(12, max_col + 1)


def plate_label(sample_name, group):
    """Heading for a plate grid: the sample name, plus the plate/run id(s)
    embedded in the Sanger trace filenames when sanger_name_regex captured
    one (see sanger_stats.py's PLATE_RUN_COLUMN) - absent for most projects,
    whose clone ids carry no such id.
    """
    label = f"Plate: {sample_name}"
    if "run" in group.columns:
        runs = sorted({r for r in group["run"].dropna() if r})
        if runs:
            label += f" (run {', '.join(runs)})"
    return label


def tidy_numeric_dtypes(df):
    """Undo the int->float64 coercion a CSV round trip causes whenever an
    otherwise-integer column has a blank cell (e.g. a summary row that
    leaves position columns empty).

    pandas can't write "this column is really an integer" to a plain CSV, so
    read_csv falls back to float64 the moment any cell is blank - a whole
    coordinate like 15013799 then either grows a spurious ".0" or, worse,
    gets rounded to 2 decimal places by format_cell. Whole-number float
    columns are cast back to nullable Int64 so they render as plain integers
    again; genuinely fractional columns (e.g. frac_pass) are left alone.
    """
    df = df.copy()
    for col in df.columns:
        series = df[col]
        if pd.api.types.is_float_dtype(series):
            non_null = series.dropna()
            if non_null.shape[0] and (non_null == non_null.round()).all():
                df[col] = series.astype("Int64")
    return df


def to_markdown(df):
    """A minimal GitHub-flavoured markdown table, without a tabulate dependency."""
    columns = list(df.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, separator]
    for _, row in df.iterrows():
        cells = [format_cell(value) for value in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def format_cell(value):
    """One table cell as display text - shared by report.py's markdown
    table and report_pdf.py's, so a value renders the same in both. A float
    is a fraction (tidy_numeric_dtypes already casts whole-number float
    columns back to Int64), rounded to 2 decimal places rather than shown at
    full precision.
    """
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


# --- Barcode demultiplexing (demux_barcodes.py) ------------------------------

BARCODE_COLUMNS = ["plate", "side", "barcode"]
BARCODE_SIDES = ("forward", "reverse")


def prefix_distance(a, b):
    """Hamming distance between two barcodes over their shared prefix.

    Staggered barcodes differ in length, but a read is only ever compared to a
    barcode over that barcode's own length, so two of them are confusable to
    the extent that their common prefix is - never more.
    """
    n = min(len(a), len(b))
    return sum(x != y for x, y in zip(a[:n], b[:n]))


def read_barcodes(path, max_mismatch=1):
    """Read and validate a barcode sheet (columns plate, side, barcode).

    Several barcodes may share a plate and side - a plate pooled over a mix of
    staggered-length barcodes for Illumina's sake - and all of them label that
    plate. Checked here, at parse time, rather than left to misassign reads
    quietly later: every barcode has to stay distinguishable from every other
    after up to max_mismatch errors, i.e. their prefix distance has to exceed
    2*max_mismatch. That also rules out one barcode being a prefix of another,
    which would make the length to trim ambiguous.
    """
    barcodes = pd.read_table(path, comment="#", dtype=str).dropna(how="all")
    missing = [c for c in BARCODE_COLUMNS if c not in barcodes.columns]
    if missing:
        raise ValueError(
            f"Barcode sheet {path} is missing column(s) {missing}; it needs "
            f"{BARCODE_COLUMNS}."
        )
    barcodes = barcodes[BARCODE_COLUMNS].apply(lambda s: s.str.strip())
    barcodes["barcode"] = barcodes["barcode"].str.upper()
    barcodes["side"] = barcodes["side"].str.lower()
    bad_side = barcodes[~barcodes["side"].isin(BARCODE_SIDES)]
    if bad_side.shape[0]:
        raise ValueError(
            f"Barcode sheet {path}: side must be one of {BARCODE_SIDES}, got "
            f"{sorted(bad_side['side'].unique())}."
        )
    bad_seq = barcodes[~barcodes["barcode"].str.fullmatch(r"[ACGT]+")]
    if bad_seq.shape[0]:
        raise ValueError(
            f"Barcode sheet {path}: barcodes may only contain A, C, G, T; got "
            f"{sorted(bad_seq['barcode'].unique())}."
        )
    if barcodes.shape[0] == 0:
        raise ValueError(f"Barcode sheet {path} has no barcodes.")
    if barcodes["barcode"].duplicated().any():
        raise ValueError(
            f"Barcode sheet {path}: barcode(s) listed more than once: "
            f"{sorted(barcodes.loc[barcodes['barcode'].duplicated(), 'barcode'])}."
        )
    rows = list(barcodes.itertuples(index=False))
    for i, first in enumerate(rows):
        for second in rows[i + 1 :]:
            distance = prefix_distance(first.barcode, second.barcode)
            if distance <= 2 * max_mismatch:
                raise ValueError(
                    f"Barcodes {first.barcode} ({first.plate}, {first.side}) and "
                    f"{second.barcode} ({second.plate}, {second.side}) differ at "
                    f"only {distance} position(s) over their shared prefix, so "
                    f"reads can't be told apart with up to {max_mismatch} "
                    f"mismatch(es) allowed (needs more than {2 * max_mismatch})."
                )
    return barcodes.reset_index(drop=True)


# --- Per-clone NGS calls (ngs_clone_stats.py) --------------------------------

CLONE_STATUS_CLEAN = "clean"
CLONE_STATUS_CONTAMINATED = "contaminated"
CLONE_STATUS_MULTIPLE = "multiple"
CLONE_STATUS_UNMOBILIZED = "unmobilized"
CLONE_STATUS_WEAK = "weak"
CLONE_STATUS_NO_INSERTION = "no_insertion"
CLONE_STATUS_TOO_FEW_READS = "too_few_reads"

# Same palette as the Sanger clone grid, so a plate reads the same way
# whichever assay it came from: green is a usable single site, yellow is
# usable but worth a look, red is not usable, grey is not enough data to say.
CLONE_STATUS_COLORS = {
    CLONE_STATUS_CLEAN: CLONE_GREEN,
    CLONE_STATUS_UNMOBILIZED: CLONE_GREEN,
    CLONE_STATUS_CONTAMINATED: CLONE_YELLOW,
    CLONE_STATUS_WEAK: CLONE_YELLOW,
    CLONE_STATUS_MULTIPLE: CLONE_RED,
    CLONE_STATUS_NO_INSERTION: CLONE_RED,
    CLONE_STATUS_TOO_FEW_READS: READ_QC_GREY,
}

CLONE_SUMMARY_COLUMNS = [
    "sample_name",
    "plate",
    "well",
    "status",
    "n_pairs",
    "n_sites",
    "n_molecules",
    "chrom",
    "start",
    "end",
    "strand",
    "sides",
    "dominant_molecules",
    "dominant_frac",
    "second_position",
    "second_molecules",
    "second_frac",
    "n_secondary_sites",
    "contaminated_by",
    "at_original_site",
    "founder_molecules",
    "founder_frac",
    "shared_with",
    "reason",
]

# Columns of ngs_clone_summary.tsv shown in the reports - the rest (end, plate,
# well, the raw molecule counts of the runner-up) are in the TSV for whoever
# wants them, but make the table too wide to read.
CLONE_REPORT_COLUMNS = [
    "sample_name",
    "status",
    "n_pairs",
    "n_sites",
    "chrom",
    "start",
    "strand",
    "sides",
    "dominant_molecules",
    "dominant_frac",
    "second_position",
    "second_frac",
    "contaminated_by",
    "shared_with",
    "reason",
]


def clone_plate_wells(group):
    """{(row, col): (title, colour, marker)} for one plate of ngs_clone_summary
    rows, for drawing its well grid. Wells whose library name does not parse
    as a well (see parse_well) are left out, since not every project names its
    libraries after wells. Unmobilized clones get a "u", like the Sanger grid.
    """
    wells = {}
    for row in group.itertuples():
        position = parse_well(row.well)
        if position is None:
            continue
        wells[position] = (
            f"{row.well}: {row.status}",
            CLONE_STATUS_COLORS.get(row.status, READ_QC_GREY),
            "u" if row.status == CLONE_STATUS_UNMOBILIZED else None,
        )
    return wells


def clone_plate_label(plate):
    return f"Plate: {plate}" if isinstance(plate, str) and plate else "Plate"


def demux_tables(stats):
    """Per-library and per-plate summaries of demux_stats.tsv.

    The raw file has a row per library, barcode and outcome - a useful audit
    trail, but too long to read. Returns (libraries, plates): how each
    library's pairs divided into assigned/unassigned, and how the assigned
    ones divided over plates, with how many of the plate's barcodes turned up
    at all (a barcode that never does is the first sign of a pooling mistake)
    and how many matched only with mismatches.
    """
    outcomes = ["assigned", "no_barcode", "conflict", "primer_mismatch"]
    totals = stats.pivot_table(
        index="source", columns="category", values="n_pairs", aggfunc="sum", fill_value=0
    ).reindex(columns=outcomes, fill_value=0)
    totals["total_pairs"] = totals.sum(axis=1)
    totals["frac_assigned"] = totals["assigned"] / totals["total_pairs"].where(
        totals["total_pairs"] > 0
    )
    totals.columns.name = None
    libraries = (
        totals.rename(columns={"assigned": "assigned_pairs"})
        .reset_index()
        .rename(columns={"source": "library"})[
            [
                "library",
                "total_pairs",
                "assigned_pairs",
                "frac_assigned",
                "no_barcode",
                "conflict",
                "primer_mismatch",
            ]
        ]
    )

    assigned = stats[stats["category"] == "assigned"].copy()
    assigned["with_mismatch"] = assigned["n_pairs"].where(assigned["mismatches"] > 0, 0)
    plates = assigned.groupby(["source", "plate"]).agg(
        assigned_pairs=("n_pairs", "sum"), with_mismatch=("with_mismatch", "sum")
    )
    # A barcode that matched both exactly and with mismatches has a row for
    # each, so seen/listed are counted per distinct barcode.
    per_barcode = assigned.groupby(["source", "plate", "barcode"])["n_pairs"].sum()
    plates["barcodes"] = (
        (per_barcode > 0).groupby(level=["source", "plate"]).sum().astype(int).astype(str)
        + "/"
        + per_barcode.groupby(level=["source", "plate"]).size().astype(str)
    )
    plates["frac_of_library"] = plates["assigned_pairs"] / plates.index.get_level_values(
        "source"
    ).map(totals["total_pairs"]).to_numpy()
    plates["frac_with_mismatch"] = plates["with_mismatch"] / plates["assigned_pairs"].where(
        plates["assigned_pairs"] > 0
    )
    plates = (
        plates.reset_index()
        .rename(columns={"source": "library"})[
            [
                "library",
                "plate",
                "assigned_pairs",
                "frac_of_library",
                "barcodes",
                "frac_with_mismatch",
            ]
        ]
    )
    return libraries, plates
