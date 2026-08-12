"""Build the miniature dataset the workflow is tested against.

Everything is synthetic and tiny: a short genome, a short cassette, and reads
planted at integration sites we know the answer for. Both branches are built
from the same insertions, so the cross-validation step has something to agree
about.

Geometry, which the reads have to respect for the pipeline to find anything:
the forward ITR primer sits near the end of the cassette contig and reads
outwards past it, the reverse primer sits near the start and reads outwards
past position 0. Which way each end runs in the genome depends on the
orientation the cassette integrated in.

Run from the .test directory:  python3 make_test_data.py
"""

import gzip
import os
import random

random.seed(20240117)

HERE = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(HERE, "resources")
DATA = os.path.join(HERE, "data")

CHROM = "chr1"
CHROM_LENGTH = 30000
CASSETTE = "test_cassette"
CASSETTE_LENGTH = 2000

FORWARD_PRIMER = "GGTGATCCTAACTGACCTAAGAC"
REVERSE_PRIMER = "AACGAGTTTTAATGACTCCAACTT"
INSERTION_SEQ = "TA"
# Long enough that the genomic part of a read that crosses out of the cassette
# still aligns on its own, which is what makes the pair chimeric and gives
# pairtools a junction to report.
READ_LENGTH = 100
# How much cassette a read from an ITR primer crosses before reaching the genome
CASSETTE_TAIL = 40

# Where the cassette sits in the genome, and which way round.
INSERTIONS = [
    {"pos": 8000, "strand": "+"},
    {"pos": 17000, "strand": "-"},
]

# A short standalone contig - like SB_launchpad in the real mobilization_FR0_pools
# project - too small to hold a long read. Deliberately independent random
# sequence rather than an excerpt of `genome`, so there is no incidental
# exact-match tie between the two contigs to muddy the test. Exercises
# sanger_sites.py's runs_off_contig_end: a read that crosses the cassette,
# the whole of this contig, and keeps going should still pass QC, landing on
# a second, unrelated, real locus rather than failing as "multiple genomic
# alignments".
LAUNCHPAD_CONTIG = "test_launchpad"
LAUNCHPAD_FLANK = 40
# Where the escaping read's tail lands - far from both INSERTIONS and clear
# of CHROM_LENGTH, so it can't be chained into one contiguous alignment with
# the launchpad contig.
ESCAPE_POS = 26000
ESCAPE_LENGTH = 300

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(seq):
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def random_seq(length):
    return "".join(random.choice("ACGT") for _ in range(length))


def build_genome():
    """A genome unique enough for bwa to place short reads unambiguously."""
    chrom = list(random_seq(CHROM_LENGTH))
    # Guarantee the integration motif right at each planted site, so that
    # pinpointing has something to find.
    for insertion in INSERTIONS:
        chrom[insertion["pos"] : insertion["pos"] + len(INSERTION_SEQ)] = list(
            INSERTION_SEQ
        )
    return "".join(chrom)


def build_cassette():
    """Cassette with an ITR primer close to each end.

    The primers sit in the ITRs, so a read from either of them crosses out of
    the cassette well within its length - which is what the NGS filter's
    min(cassette_length, primer + read_len) is capping at.
    """
    cassette = list(random_seq(CASSETTE_LENGTH))
    forward_at = CASSETTE_LENGTH - CASSETTE_TAIL
    cassette[forward_at : forward_at + len(FORWARD_PRIMER)] = list(FORWARD_PRIMER)
    reverse_revcomp = revcomp(REVERSE_PRIMER)
    reverse_at = CASSETTE_TAIL - len(reverse_revcomp)
    cassette[reverse_at : reverse_at + len(reverse_revcomp)] = list(reverse_revcomp)
    return "".join(cassette)


def build_launchpad():
    """A short contig holding just the TA junction plus a short flank."""
    contig = list(random_seq(len(INSERTION_SEQ) + LAUNCHPAD_FLANK))
    contig[: len(INSERTION_SEQ)] = list(INSERTION_SEQ)
    return "".join(contig)


def make_escape_read(cassette, launchpad, genome):
    """A read that crosses the cassette, the whole of a short contig, and
    keeps going into unrelated real genomic sequence - reproducing a read
    that outruns a short reference contig."""
    forward_start, _ = primer_bounds(cassette)
    cassette_part = cassette[forward_start:]
    return cassette_part + launchpad + genome[ESCAPE_POS : ESCAPE_POS + ESCAPE_LENGTH]


def primer_bounds(cassette):
    """Where a read from each primer starts, in cassette coordinates."""
    forward_start = cassette.index(FORWARD_PRIMER)
    reverse_end = cassette.index(revcomp(REVERSE_PRIMER)) + len(REVERSE_PRIMER)
    return forward_start, reverse_end


def runs_rightwards(side, strand):
    """Whether this end of the cassette reads into the genome to the right."""
    return (side == "forward") == (strand == "+")


def write_fasta(path, contigs):
    with open(path, "w") as f:
        for name, seq in contigs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def write_fastq(path, reads):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt") as f:
        for name, seq in reads:
            f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def genomic_outward(genome, insertion, side, length):
    """Genomic sequence read into on leaving the cassette at one end.

    The insertion motif is duplicated on both sides of an integration, so a
    read from either primer runs into it before continuing into unique
    sequence.
    """
    pos = insertion["pos"]
    if runs_rightwards(side, insertion["strand"]):
        return genome[pos : pos + length]
    hi = pos + len(INSERTION_SEQ)
    return revcomp(genome[hi - length : hi])


def read_from_primer(genome, cassette, insertion, side, length):
    """A read starting at an ITR primer, out through the cassette and beyond."""
    forward_start, reverse_end = primer_bounds(cassette)
    if side == "forward":
        cassette_part = cassette[forward_start:]
    else:
        cassette_part = revcomp(cassette[:reverse_end])
    genomic_part = genomic_outward(
        genome, insertion, side, max(0, length - len(cassette_part))
    )
    return (cassette_part + genomic_part)[:length]


def genomic_mate(genome, insertion, side, offset):
    """The genomic mate, reading back towards the junction from `offset` away.

    The Tn5 cut lands at varying distances from the insertion, which is what
    spreads the reads into the peaks the NGS branch calls.
    """
    pos = insertion["pos"]
    if runs_rightwards(side, insertion["strand"]):
        lo = pos + len(INSERTION_SEQ) + offset
        return revcomp(genome[lo : lo + READ_LENGTH])
    hi = pos - offset
    return genome[hi - READ_LENGTH : hi]


def make_ngs_reads(genome, cassette, sample, n_per_site=40):
    """Pairs of a genomic mate and a cassette mate reading out of the ITR."""
    r1, r2 = [], []
    for i, insertion in enumerate(INSERTIONS):
        for side in ("forward", "reverse"):
            cassette_mate = read_from_primer(
                genome, cassette, insertion, side, READ_LENGTH
            )
            for offset in range(n_per_site):
                name = f"{sample}_{i}_{side}_{offset}"
                r1.append((name, genomic_mate(genome, insertion, side, offset)))
                r2.append((name, cassette_mate))
    return r1, r2


def make_sanger_reads(genome, cassette, plate, direction, wells, length=540):
    letter = "F" if direction == "forward" else "R"
    return [
        (
            f"{well}_{plate}_{letter}",
            read_from_primer(genome, cassette, insertion, direction, length),
        )
        for well, insertion in zip(wells, INSERTIONS)
    ]


if __name__ == "__main__":
    os.makedirs(RESOURCES, exist_ok=True)
    os.makedirs(DATA, exist_ok=True)

    genome = build_genome()
    cassette = build_cassette()
    launchpad = build_launchpad()

    write_fasta(
        os.path.join(RESOURCES, "genome.fa"),
        {CHROM: genome, CASSETTE: cassette, LAUNCHPAD_CONTIG: launchpad},
    )
    with open(os.path.join(RESOURCES, "chromsizes.txt"), "w") as f:
        f.write(
            f"{CHROM}\t{CHROM_LENGTH}\n{CASSETTE}\t{CASSETTE_LENGTH}\n"
            f"{LAUNCHPAD_CONTIG}\t{len(launchpad)}\n"
        )
    with open(os.path.join(RESOURCES, "chromsizes_no_cassette.txt"), "w") as f:
        f.write(f"{CHROM}\t{CHROM_LENGTH}\n")

    r1, r2 = make_ngs_reads(genome, cassette, "clone1")
    write_fastq(os.path.join(DATA, "clone1.R1.fastq.gz"), r1)
    write_fastq(os.path.join(DATA, "clone1.R2.fastq.gz"), r2)

    # One plate sequenced from the forward primer for both clones, and a second
    # that only worked from the reverse primer, and only for one of them - the
    # case where the NGS library has to vouch for a single Sanger read.
    write_fastq(
        os.path.join(DATA, "plate1.fastq.gz"),
        make_sanger_reads(genome, cassette, "001", "forward", ["A01", "B01"])
        + [("C01_001_F", make_escape_read(cassette, launchpad, genome))],
    )
    write_fastq(
        os.path.join(DATA, "plate2.fastq.gz"),
        make_sanger_reads(genome, cassette, "002", "reverse", ["A01", "B01"])[:1],
    )

    # Single-base intervals at the T/A boundary - the same coordinate
    # sanger_sites.py/find_insertion_sites.py report (tagmaplib's
    # insertion_site_from_motif), one base into each planted "TA", rather
    # than the motif's own start.
    expected = "\n".join(
        f"{CHROM}\t{ins['pos'] + len(INSERTION_SEQ) // 2}\t"
        f"{ins['pos'] + len(INSERTION_SEQ) // 2 + 1}\t{ins['strand']}"
        for ins in INSERTIONS
    )
    with open(os.path.join(RESOURCES, "expected_sites.bed"), "w") as f:
        f.write(expected + "\n")
    forward_start, reverse_end = primer_bounds(cassette)
    print(f"forward primer starts at {forward_start}, reverse ends at {reverse_end}")
    print("Planted insertions:")
    print(expected)
