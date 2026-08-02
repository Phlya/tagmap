"""Helpers shared between the NGS and Sanger branches of the workflow.

The scripts in this folder are meant to stay runnable by hand, so they keep
their own argparse interfaces and simply ``import tagmaplib`` - Python puts a
script's own directory on ``sys.path``, so no installation is needed.
"""

import numpy as np
import pandas as pd
import pysam

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
]

# Columns of the final insertion site files, shared by the NGS and Sanger
# branches so that the two can be compared and loaded the same way.
SITE_COLUMNS = ["chrom", "start", "end", "sample_name", "score", "strand"]

# Columns of the bedgraph files written by coverage.py.
COVERAGE_COLUMNS = ["chrom", "start", "end", "count", "fraction"]


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


def write_bed(df, path, columns=None, track_name=None, header=False):
    """Write a BED/bedGraph, optionally with a UCSC track line."""
    columns = list(df.columns) if columns is None else columns
    mode = "w"
    if track_name is not None:
        with open(path, "w") as f:
            f.write(f"track name={track_name}\n")
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


def to_markdown(df):
    """A minimal GitHub-flavoured markdown table, without a tabulate dependency."""
    columns = list(df.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, separator]
    for _, row in df.iterrows():
        cells = [_format_cell(value) for value in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def _format_cell(value):
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.3g}"
    return str(value)
