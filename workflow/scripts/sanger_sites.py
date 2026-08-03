"""Find integration sites in individual long Sanger reads.

Each read is sequenced outwards from an ITR primer, so it starts inside the
cassette and crosses into the genome at the integration site. bwa reports that
as a primary plus one or more supplementary alignments; the integration site is
the genomic coordinate at the 5' end - in read orientation - of the genomic
alignment. Unlike the NGS branch there is no peak calling: one read is one
observation of one site, which is what makes this usable on clonal lines.

Writes a table of every read with the QC columns behind the pass/fail decision,
so that reads which were dropped can be inspected rather than just disappearing.
"""

import argparse
import json
import re

import numpy as np
import pandas as pd
import pysam

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--bam", required=True)
argparser.add_argument("--genome", "-g", required=True)
argparser.add_argument(
    "--genome-index",
    default=None,
    help="pyfastx index of --genome. Built next to the FASTA if not given.",
)
argparser.add_argument("--primer-positions", required=True)
argparser.add_argument("--sample-name", required=True)
argparser.add_argument(
    "--construct-contigs",
    nargs="+",
    default=[],
    help="Contigs that are part of the construct rather than the genome",
)
argparser.add_argument("--insertion-seq", default="TA")
argparser.add_argument(
    "--insertion-seq-window",
    type=int,
    default=5,
    help="How far around the junction to look for the insertion motif",
)
argparser.add_argument("--min-mapq", type=int, default=30)
argparser.add_argument(
    "--min-aligned", type=int, default=30, help="Minimum aligned length in the genome"
)
argparser.add_argument(
    "--max-unexplained",
    type=int,
    default=50,
    help="Maximum number of bases before the genomic alignment that match "
    "neither the genome nor the construct",
)
argparser.add_argument(
    "--max-dist",
    type=int,
    default=100,
    help="Tolerance when checking that a read starts at the ITR primer",
)
argparser.add_argument(
    "--direction",
    choices=["forward", "reverse"],
    default="forward",
    help="Which ITR primer these reads were sequenced from, used unless the "
    "read name says otherwise",
)
argparser.add_argument(
    "--direction-map",
    nargs="*",
    default=["F=forward", "R=reverse"],
    help="How to read the direction group of --name-regex, e.g. F=forward",
)
argparser.add_argument(
    "--name-regex",
    default=None,
    help="Regex with named groups matched against read names, e.g. "
    r"'(?P<well>[A-H]\d{2})_(?P<run>\d+)'. Recognised groups are well, clone "
    "and direction; others are kept as extra columns.",
)
argparser.add_argument("--clone", default=None, help="Clone id for all reads")
argparser.add_argument("--output-reads", required=True)
argparser.add_argument("--output-sites", required=True)

CLIP_OPS = (4, 5)  # soft clip, hard clip


def read_coordinates(aln):
    """Where an alignment sits in the read, in original read orientation.

    pysam reports query coordinates against the stored sequence, which is
    reverse complemented for reverse alignments and hard clipped for
    supplementary ones, so work it out from the cigar instead.
    """
    cigar = aln.cigartuples
    lead = 0
    for op, length in cigar:
        if op not in CLIP_OPS:
            break
        lead += length
    trail = 0
    for op, length in reversed(cigar):
        if op not in CLIP_OPS:
            break
        trail += length
    total = aln.infer_read_length()
    aligned = total - lead - trail
    start = trail if aln.is_reverse else lead
    return start, start + aligned, total


def alignment_table(bam_path, construct_contigs):
    """One row per alignment, with read-orientation query coordinates."""
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary:
                continue
            if aln.is_unmapped:
                rows.append(
                    {
                        "readname": aln.query_name,
                        "chrom": None,
                        "start": np.nan,
                        "end": np.nan,
                        "strand": ".",
                        "mapq": 0,
                        "read_start": np.nan,
                        "read_end": np.nan,
                        "read_length": aln.query_length or np.nan,
                        "is_construct": False,
                        "mapped": False,
                    }
                )
                continue
            read_start, read_end, read_length = read_coordinates(aln)
            rows.append(
                {
                    "readname": aln.query_name,
                    "chrom": aln.reference_name,
                    "start": aln.reference_start,
                    "end": aln.reference_end,
                    "strand": "-" if aln.is_reverse else "+",
                    "mapq": aln.mapping_quality,
                    "read_start": read_start,
                    "read_end": read_end,
                    "read_length": read_length,
                    "is_construct": aln.reference_name in construct_contigs,
                    "mapped": True,
                }
            )
    return pd.DataFrame(rows)


def parse_read_name(name, regex, direction_map):
    """Pull well/clone/direction and any extra groups out of a read name."""
    if regex is None:
        return {}
    match = regex.search(name)
    if match is None:
        return {}
    fields = {k: v for k, v in match.groupdict().items() if v is not None}
    if "direction" in fields:
        raw = fields["direction"]
        fields["direction"] = direction_map.get(raw, direction_map.get(raw.upper(), raw))
    return fields


def uncovered_before(segments, position):
    """Bases before `position` in the read covered by no alignment at all."""
    covered = np.zeros(int(position), dtype=bool)
    for segment in segments:
        if not segment["mapped"]:
            continue
        lo = int(min(segment["read_start"], position))
        hi = int(min(segment["read_end"], position))
        covered[lo:hi] = True
    return int((~covered).sum())


def analyse_read(segments, primer_positions, args, default_direction):
    """Turn one read's alignments into a site plus the QC behind it."""
    segments = sorted(segments, key=lambda s: (not s["mapped"], s["read_start"]))
    mapped = [s for s in segments if s["mapped"]]
    construct = [s for s in mapped if s["is_construct"]]
    genomic = [
        s
        for s in mapped
        if not s["is_construct"]
        and s["mapq"] >= args.min_mapq
        and (s["read_end"] - s["read_start"]) >= args.min_aligned
    ]

    result = {
        "n_alignments": len(mapped),
        "n_genome_segments": len(genomic),
        "has_cassette_segment": len(construct) > 0,
        "read_length": segments[0]["read_length"] if segments else np.nan,
    }

    # The read runs outwards from the primer, so the genomic segment nearest the
    # start of the read is the one holding the junction.
    site = genomic[0] if genomic else None

    if construct:
        # Likewise the construct segment nearest the start of the read is the
        # one that should begin at the primer.
        first = construct[0]
        if default_direction == "forward":
            expected = primer_positions["forward_ITR_primer_position"]
            offset = first["start"] - expected
            expected_strand = "+"
        else:
            expected = primer_positions["reverse_ITR_primer_position"]
            offset = first["end"] - expected
            expected_strand = "-"
        result["cassette_start_offset"] = offset
        result["correct_cassette_start"] = bool(
            abs(offset) <= args.max_dist and first["strand"] == expected_strand
        )
    else:
        result["cassette_start_offset"] = np.nan
        result["correct_cassette_start"] = False

    if site is None:
        result.update(
            {
                "chrom": None,
                "start": np.nan,
                "end": np.nan,
                "strand": ".",
                "mapq": max([s["mapq"] for s in mapped], default=0),
                "unexplained_5p_bases": np.nan,
                "pass": False,
                "fail_reason": "no genomic alignment" if mapped else "unmapped",
            }
        )
        return result

    junction = tagmaplib.junction_position(site["start"], site["end"], site["strand"])
    unexplained = uncovered_before(mapped, site["read_start"])

    reasons = []
    if len(genomic) > 1:
        reasons.append("multiple genomic alignments")
    if unexplained > args.max_unexplained:
        reasons.append(f"{unexplained} unexplained bases before the junction")

    result.update(
        {
            "chrom": site["chrom"],
            "start": int(junction),
            "end": int(junction) + 1,
            "strand": site["strand"],
            "mapq": site["mapq"],
            "aligned_length": int(site["read_end"] - site["read_start"]),
            "unexplained_5p_bases": unexplained,
            "pass": not reasons,
            "fail_reason": "; ".join(reasons),
        }
    )
    return result


if __name__ == "__main__":
    args = argparser.parse_args()

    with open(args.primer_positions) as f:
        primer_positions = json.load(f)
    direction_map = dict(item.split("=", 1) for item in args.direction_map)
    name_regex = re.compile(args.name_regex) if args.name_regex else None

    alignments = alignment_table(args.bam, set(args.construct_contigs))

    columns = [
        "sample_name",
        "readname",
        "clone",
        "well",
        "direction",
        "chrom",
        "start",
        "end",
        "strand",
        "mapq",
        "read_length",
        "aligned_length",
        "n_alignments",
        "n_genome_segments",
        "has_cassette_segment",
        "cassette_start_offset",
        "correct_cassette_start",
        "unexplained_5p_bases",
        f"{args.insertion_seq}_found",
        "pass",
        "fail_reason",
    ]

    if alignments.shape[0] == 0:
        print(f"No reads in {args.bam}")
        empty = pd.DataFrame({c: pd.Series(dtype=object) for c in columns})
        empty.to_csv(args.output_reads, sep="\t", index=False)
        empty.head(0).to_csv(args.output_sites, sep="\t", index=False, header=False)
        raise SystemExit(0)

    reads = []
    for readname, group in alignments.groupby("readname", sort=True):
        fields = parse_read_name(readname, name_regex, direction_map)
        direction = fields.get("direction", args.direction)
        if direction not in ("forward", "reverse"):
            raise ValueError(
                f"Read {readname!r} has direction {direction!r}, which is "
                "neither forward nor reverse. Check sanger_direction_map."
            )
        record = {
            "sample_name": args.sample_name,
            "readname": readname,
            "direction": direction,
        }
        record.update(fields)
        record.update(
            analyse_read(
                group.to_dict("records"), primer_positions, args, direction
            )
        )
        record["clone"] = (
            fields.get("clone") or args.clone or fields.get("well") or readname
        )
        reads.append(record)

    reads = pd.DataFrame(reads)

    # Snap onto the motif the transposon integrates into. Reads of one clone
    # then agree exactly rather than to within a base or two.
    motif_start = tagmaplib.find_insertion_seq(
        reads,
        args.genome,
        args.insertion_seq,
        window=args.insertion_seq_window,
        mode="nearest",
        index_file=args.genome_index,
    )
    found = motif_start >= 0
    reads[f"{args.insertion_seq}_found"] = found
    reads.loc[found, "start"] = motif_start[found]
    reads.loc[found, "end"] = motif_start[found] + len(args.insertion_seq)

    # Reads from the two primers run opposite ways along the read, so one of
    # them has to be flipped for a site to have one orientation rather than two.
    # Which one is fixed by the NGS branch: determine_direction in
    # find_insertion_sites.py calls orientation from the order of the forward
    # and reverse peaks, and flipping the forward reads here is what makes the
    # two branches agree about the same insertion.
    reads["strand"] = np.where(
        reads["direction"] == "forward",
        tagmaplib.flip_strand(reads["strand"]),
        reads["strand"],
    )

    for column in columns:
        if column not in reads.columns:
            reads[column] = pd.NA
    extra = [c for c in reads.columns if c not in columns]
    reads = reads[columns + extra].sort_values(["clone", "direction", "readname"])
    reads.to_csv(args.output_reads, sep="\t", index=False)

    sites = reads[reads["pass"].astype(bool)].copy()
    sites[["start", "end"]] = sites[["start", "end"]].astype(int)
    sites[
        ["chrom", "start", "end", "sample_name", "mapq", "strand", "clone", "readname",
         "direction"]
    ].sort_values(["chrom", "start", "end"]).to_csv(
        args.output_sites, sep="\t", index=False, header=False
    )

    print(
        f"{reads.shape[0]} reads, {sites.shape[0]} passing, "
        f"{reads['clone'].nunique()} clones in {args.sample_name}"
    )
