"""Locate the ITR primers in the cassette.

The NGS branch selects read pairs whose cassette side ends exactly at a primer,
so these positions decide what counts as a real tagmentation product; the
Sanger branch uses them to check that a read starts where it should. Getting
them silently wrong would quietly empty the output, so a primer that cannot be
found is an error rather than a -1.
"""

import argparse
import json

import pyfastx

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--genome", required=True)
argparser.add_argument("--cassette-name", required=True)
argparser.add_argument("--forward-primer", required=True)
argparser.add_argument("--reverse-primer", required=True)
argparser.add_argument("--output", "-o", required=True)

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(seq):
    return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))


def find_unique(seq, primer, label, cassette_name):
    """Position of a primer in the cassette, on whichever strand it sits."""
    forward_hits = [i for i in range(len(seq)) if seq.startswith(primer, i)]
    reverse_hits = [i for i in range(len(seq)) if seq.startswith(revcomp(primer), i)]
    hits = [(i, "+") for i in forward_hits] + [(i, "-") for i in reverse_hits]
    if not hits:
        raise ValueError(
            f"The {label} primer {primer!r} does not occur in cassette "
            f"{cassette_name!r}, on either strand. Only the part of the primer "
            "that binds the ITR should be given, without the overhang."
        )
    if len(hits) > 1:
        positions = ", ".join(f"{i}({strand})" for i, strand in hits)
        raise ValueError(
            f"The {label} primer {primer!r} occurs {len(hits)} times in "
            f"cassette {cassette_name!r}, at {positions}. A primer that binds "
            "in more than one place cannot define a unique junction."
        )
    return hits[0]


if __name__ == "__main__":
    args = argparser.parse_args()

    seq = None
    for name, contig in pyfastx.Fasta(args.genome, build_index=False):
        if name == args.cassette_name:
            seq = contig.upper()
            break
    if seq is None:
        raise ValueError(
            f"Cassette {args.cassette_name!r} is not a contig of {args.genome}."
        )

    forward_start, _ = find_unique(
        seq, args.forward_primer.upper(), "forward", args.cassette_name
    )
    reverse_start, _ = find_unique(
        seq, args.reverse_primer.upper(), "reverse", args.cassette_name
    )

    # The forward primer points along the cassette, so reads start at its 5'
    # end; the reverse primer points the other way, so reads start at its 3'
    # end, which is the far side of the match.
    to_save = {
        "forward_ITR_primer_position": forward_start + 1,
        "reverse_ITR_primer_position": reverse_start + len(args.reverse_primer),
    }
    with open(args.output, "w") as f:
        json.dump(to_save, f)
    print(to_save)
