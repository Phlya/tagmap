"""Split one paired-end FASTQ by the plate barcode at the start of a read.

Barcoded ITR primers put a short barcode in front of the ITR sequence, so the
read that starts at the cassette begins

    [barcode][ITR primer][cassette ...]

A plate may use one barcode per side or a mix of several - staggered in length
so that the first bases of the library are not all the same, which Illumina
dislikes - and each read carries a random one of them. Forward-set barcodes sit
on forward-ITR reads and reverse-set ones on reverse-ITR reads, so the set a
barcode belongs to also says which ITR primer the read should start with.

Every barcode listed for a plate labels that plate, whichever it is. The
matched barcode is trimmed off, so what follows is read exactly as in a library
without barcodes - the read starts at the ITR primer.

Pairs with no barcode, with barcodes that disagree about the plate, or (with
--require-primer) whose barcode is not followed by its ITR primer are not
written anywhere, only counted.
"""

import argparse
import gzip
from collections import Counter

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--r1", required=True)
argparser.add_argument("--r2", required=True)
argparser.add_argument("--barcodes", required=True, help="TSV: plate, side, barcode")
argparser.add_argument("--source", required=True, help="Name of the input library")
argparser.add_argument(
    "--plates", nargs="+", required=True, help="Plates to write, in output order"
)
argparser.add_argument("--out-r1", nargs="+", required=True)
argparser.add_argument("--out-r2", nargs="+", required=True)
argparser.add_argument("--stats", required=True)
argparser.add_argument("--max-mismatch", type=int, default=1)
argparser.add_argument(
    "--barcode-reads",
    choices=["R1", "R2", "either"],
    default="R2",
    help="Which read(s) to look for the barcode at the start of. The ITR read "
    "is the one that carries it - R2 in the usual TagMap layout.",
)
argparser.add_argument("--forward-primer", default=None)
argparser.add_argument("--reverse-primer", default=None)
argparser.add_argument(
    "--require-primer",
    action="store_true",
    help="Only accept a barcode followed by the ITR primer of its own side.",
)
argparser.add_argument("--primer-max-mismatch", type=int, default=2)
argparser.add_argument(
    "--primer-check-length",
    type=int,
    default=15,
    help="How much of the primer, from its 5' end, has to follow the barcode",
)

STATS_COLUMNS = ["source", "category", "plate", "side", "barcode", "read", "mismatches", "n_pairs"]
ASSIGNED = "assigned"
NO_BARCODE = "no_barcode"
CONFLICT = "conflict"
PRIMER_MISMATCH = "primer_mismatch"


def hamming(a, b):
    """Mismatches between two equal-length strings; N never matches."""
    return sum(x != y or x == "N" for x, y in zip(a, b))


class Matcher:
    """Finds which barcode a sequence starts with, if any."""

    def __init__(self, barcodes, max_mismatch):
        self.barcodes = list(barcodes.itertuples(index=False))
        self.max_mismatch = max_mismatch
        # Exact hits are the overwhelming majority, so look those up directly.
        self.exact = {bc.barcode: bc for bc in self.barcodes}
        self.lengths = sorted({len(bc.barcode) for bc in self.barcodes})

    def match(self, seq):
        """(mismatches, barcode row) of the closest barcode within the
        mismatch limit, or None."""
        for length in self.lengths:
            hit = self.exact.get(seq[:length])
            if hit is not None:
                return 0, hit
        if self.max_mismatch == 0:
            return None
        best = None
        for bc in self.barcodes:
            length = len(bc.barcode)
            if len(seq) < length:
                continue
            distance = hamming(seq[:length], bc.barcode)
            if distance <= self.max_mismatch and (best is None or distance < best[0]):
                best = (distance, bc)
        return best


def read_fastq(path):
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                return
            yield header, f.readline().rstrip("\n"), f.readline(), f.readline().rstrip("\n")


def write_record(handle, header, seq, plus, qual, trim=0):
    handle.write(f"{header}{seq[trim:]}\n{plus}{qual[trim:]}\n")


def primer_follows(seq, offset, primer, args):
    if primer is None:
        return True
    primer = primer.upper()[: args.primer_check_length]
    window = seq[offset : offset + len(primer)]
    return len(window) == len(primer) and hamming(window, primer) <= args.primer_max_mismatch


if __name__ == "__main__":
    args = argparser.parse_args()
    if not (len(args.plates) == len(args.out_r1) == len(args.out_r2)):
        raise SystemExit("--plates, --out-r1 and --out-r2 need one entry per plate")
    if args.require_primer and not (args.forward_primer and args.reverse_primer):
        raise SystemExit("--require-primer needs --forward-primer and --reverse-primer")

    barcodes = tagmaplib.read_barcodes(args.barcodes, args.max_mismatch)
    unknown = sorted(set(barcodes["plate"]) - set(args.plates))
    if unknown:
        raise SystemExit(f"Plates {unknown} are in the barcode sheet but not in --plates")
    matcher = Matcher(barcodes, args.max_mismatch)
    primers = {"forward": args.forward_primer, "reverse": args.reverse_primer}
    check_primer = args.require_primer
    reads_to_check = {"R1": (0,), "R2": (1,), "either": (0, 1)}[args.barcode_reads]

    outputs = {
        plate: (gzip.open(r1, "wt"), gzip.open(r2, "wt"))
        for plate, r1, r2 in zip(args.plates, args.out_r1, args.out_r2)
    }
    counts = Counter()

    try:
        for rec1, rec2 in zip(read_fastq(args.r1), read_fastq(args.r2)):
            records = (rec1, rec2)
            hits = []  # (read index, mismatches, barcode row)
            primer_failed = False
            for index in reads_to_check:
                found = matcher.match(records[index][1])
                if found is None:
                    continue
                mismatches, barcode = found
                if check_primer and not primer_follows(
                    records[index][1],
                    len(barcode.barcode),
                    primers[barcode.side],
                    args,
                ):
                    primer_failed = True
                    continue
                hits.append((index, mismatches, barcode))

            if not hits:
                counts[(PRIMER_MISMATCH if primer_failed else NO_BARCODE, "", "", "", "", "")] += 1
                continue
            labels = {(barcode.plate, barcode.side) for _, _, barcode in hits}
            if len(labels) > 1:
                counts[(CONFLICT, "", "", "", "", "")] += 1
                continue

            index, mismatches, barcode = hits[0]
            counts[
                (ASSIGNED, barcode.plate, barcode.side, barcode.barcode, f"R{index + 1}", mismatches)
            ] += 1
            trims = [0, 0]
            for hit_index, _, hit_barcode in hits:
                trims[hit_index] = len(hit_barcode.barcode)
            out1, out2 = outputs[barcode.plate]
            write_record(out1, *rec1, trim=trims[0])
            write_record(out2, *rec2, trim=trims[1])
    finally:
        for out1, out2 in outputs.values():
            out1.close()
            out2.close()

    # Every plate and barcode gets a row, even with nothing assigned, so a
    # barcode that never turns up reads as a zero rather than going missing.
    for barcode in barcodes.itertuples(index=False):
        if not any(
            key[0] == ASSIGNED and key[3] == barcode.barcode for key in counts
        ):
            counts[(ASSIGNED, barcode.plate, barcode.side, barcode.barcode, "", 0)] += 0

    with open(args.stats, "w") as f:
        f.write("\t".join(STATS_COLUMNS) + "\n")
        for key, n in sorted(counts.items(), key=lambda item: tuple(map(str, item[0]))):
            f.write("\t".join([args.source, *map(str, key), str(n)]) + "\n")

    total = sum(counts.values())
    assigned = sum(n for key, n in counts.items() if key[0] == ASSIGNED)
    print(f"{args.source}: {assigned} of {total} pairs assigned to a plate")
