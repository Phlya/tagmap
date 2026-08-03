"""Resolve the genomic position of the original, pre-mobilization insertion.

A founder/control clone can still carry the cassette at its starting locus
rather than a new one - its Sanger reads pass QC and look like a perfectly
good site, but they are not evidence of mobilization. Comparing against this
position downstream (sanger_stats.py) is what tells the two apart. The locus
can be given directly as chrom:pos, or - when it is only known by its
flanking sequence - located once here by searching the genome for the
sequence immediately upstream of it.
"""

import argparse
import json

import pyfastx

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--genome", required=True)
argparser.add_argument(
    "--construct-contigs",
    nargs="*",
    default=[],
    help="Contigs to skip - the cassette itself, not the genome it inserted into",
)
argparser.add_argument(
    "--site", default=None, help="chrom:pos (0-based), given directly"
)
argparser.add_argument(
    "--upstream-seq",
    default=None,
    help="Genome sequence immediately upstream of the site, on the forward strand",
)
argparser.add_argument("--output", "-o", required=True)


def find_upstream(genome_path, construct_contigs, upstream_seq):
    """The one place upstream_seq occurs, as the (chrom, pos) right after it."""
    upstream_seq = upstream_seq.upper()
    hits = []
    for name, seq in pyfastx.Fasta(genome_path, build_index=False):
        if name in construct_contigs:
            continue
        seq = seq.upper()
        i = seq.find(upstream_seq)
        while i != -1:
            hits.append((name, i + len(upstream_seq)))
            i = seq.find(upstream_seq, i + 1)
    if not hits:
        raise ValueError(
            f"{upstream_seq!r} does not occur in {genome_path} outside "
            f"{sorted(construct_contigs)}. original_insertion_upstream_seq must "
            "be exact genome sequence, on the forward strand, immediately "
            "upstream of the original insertion."
        )
    if len(hits) > 1:
        positions = ", ".join(f"{chrom}:{pos}" for chrom, pos in hits)
        raise ValueError(
            f"{upstream_seq!r} occurs {len(hits)} times: {positions}. It must "
            "be unique to define a single original insertion site."
        )
    return hits[0]


if __name__ == "__main__":
    args = argparser.parse_args()

    if bool(args.site) == bool(args.upstream_seq):
        raise ValueError("Give exactly one of --site and --upstream-seq.")

    if args.site:
        chrom, pos = args.site.rsplit(":", 1)
        pos = int(pos)
    else:
        chrom, pos = find_upstream(
            args.genome, set(args.construct_contigs), args.upstream_seq
        )

    with open(args.output, "w") as f:
        json.dump({"chrom": chrom, "pos": pos}, f)
    print(f"Original insertion site: {chrom}:{pos}")
