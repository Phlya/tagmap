"""Drop molecules that too few reads support.

After deduplication each retained pair stands for one molecule - one physical
fragment that went into the PCR - and its PCR copies sit in the same file,
marked as duplicates. In an over-amplified library a real fragment is copied
many times over, so the number of reads sharing a molecule's two ends says how
much amplification it saw, which a molecule that only ever appeared once did
not.

Measured on a real per-clone run, that separates signal from scatter better
than anything downstream: molecules at a clone's dominant insertion site had a
median of 5 reads behind them and 34% were single-read, while molecules
anywhere else had a median of 1 read and 73% were single-read. Keeping only
molecules with at least 3 reads retained 57% of the dominant site's molecules
and 14% of the rest.

It is a trade, not a free filter: it discards about half the real signal too,
so a genuinely rare insertion captured only once is lost with the noise, and
the molecule counts everything downstream is judged on roughly halve.
Thresholds set against unfiltered counts need lowering to match.

Reads are grouped by the two ends deduplication itself keys on (the junction
pos5 and the tagmentation pos3 of each side), rather than by the parent each
duplicate was marked against: with pairtools' cython backend duplication is
not transitive, so those parent links form chains across neighbouring - and
genuinely distinct - molecules, and following them merges fragments that are
not copies of each other at all. Grouping on exact positions instead can only
split a molecule whose copies differ by the base or two of mapping wobble that
dedup's own --max-mismatch forgives, which costs a little support on a real
molecule rather than crediting support to a spurious one.
"""

import argparse
import sys
from collections import Counter

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument(
    "--dupmarked",
    required=True,
    help="The deduplication output, holding the retained pairs and the "
    "duplicates marked against them (pairtools dedup --mark-dups).",
)
argparser.add_argument(
    "--min-reads",
    type=int,
    required=True,
    help="Reads a molecule needs, counting the retained read itself",
)
argparser.add_argument(
    "--input", "-i", default=None, help="Pairs to filter (default: stdin)"
)
argparser.add_argument("--output", "-o", default=None, help="Default: stdout")

# The ends pairtools dedup was pointed at in this workflow (--p1 pos51
# --p2 pos31), plus the chromosomes they are on.
KEY_COLUMNS = ["chrom1", "pos51", "chrom2", "pos31"]


def header_columns(line):
    return line.split(":", 1)[1].split()


def key_indices(columns, path):
    missing = [c for c in KEY_COLUMNS if c not in columns]
    if missing:
        raise SystemExit(
            f"{path} is missing the column(s) {missing} that identify a "
            f"molecule; it has {', '.join(columns)}."
        )
    return [columns.index(c) for c in KEY_COLUMNS]


def read_selected(source):
    """The pairs to filter, kept as lines plus their molecule keys.

    Buffered rather than streamed so the duplicate counting below only has to
    hold counts for the molecules actually in play, instead of one entry per
    molecule in the whole library.
    """
    header, lines, keys = [], [], []
    indices = None
    for line in source:
        if line.startswith("#"):
            header.append(line)
            if line.startswith("#columns:"):
                indices = key_indices(header_columns(line.rstrip("\n")), "the input")
            continue
        if indices is None:
            raise SystemExit("the input has no #columns header line")
        fields = line.rstrip("\n").split("\t")
        lines.append(line)
        keys.append(tuple(fields[i] for i in indices))
    return header, lines, keys


def count_reads(path, wanted):
    """Reads per molecule, counting only the molecules in `wanted`."""
    counts = Counter()
    indices = None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("#columns:"):
                    indices = key_indices(header_columns(line.rstrip("\n")), path)
                continue
            if indices is None:
                raise SystemExit(f"{path} has no #columns header line")
            fields = line.rstrip("\n").split("\t")
            key = tuple(fields[i] for i in indices)
            if key in wanted:
                counts[key] += 1
    return counts


if __name__ == "__main__":
    args = argparser.parse_args()

    source = open(args.input) if args.input else sys.stdin
    try:
        header, lines, keys = read_selected(source)
    finally:
        if args.input:
            source.close()

    counts = count_reads(args.dupmarked, set(keys))

    sink = open(args.output, "w") if args.output else sys.stdout
    kept = dropped = 0
    try:
        for line in header:
            sink.write(line)
        for line, key in zip(lines, keys):
            if counts.get(key, 0) >= args.min_reads:
                sink.write(line)
                kept += 1
            else:
                dropped += 1
    finally:
        if args.output:
            sink.close()

    total = kept + dropped
    share = f"{dropped / total:.1%}" if total else "n/a"
    print(
        f"kept {kept} of {total} molecules supported by >= {args.min_reads} "
        f"read(s); dropped {dropped} ({share})",
        file=sys.stderr,
    )
