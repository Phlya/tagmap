"""Filter all_sanger_sites.bed down to fully validated sites.

A site counts as fully validated exactly when its clone renders green in the
report's plate grid (see tagmaplib.clone_color): both ITR primers
independently confirm the same site, or one primer passed and NGS
independently confirms the side the other primer left unconfirmed. A clone
whose forward and reverse reads land in different clusters entirely (see
combine_sanger_sites.py) is never green, regardless of how good any single
cluster's own evidence looks on its own - so this can never pick up more than
one site per clone.

Also writes two further copies: one with the construct/landing-pad contigs
dropped (same content as --output-for-ucsc, but a plain BED file with no
track header), and one restricted to one genomic region, for pulling out
validated clones near a specific locus of interest. And, for both the full
set and the region-restricted one, a deduplicated copy collapsing every
clone at one site down to a single row, so a hotspot several clones landed
on independently shows up once rather than as several stacked,
identical-looking features.

Before any of that, a confirmed clone's own Sanger coordinate is replaced
with the matching NGS site's coordinate where --sanger-vs-ngs says one
exists (see pinpoint_with_ngs) - NGS backs a site with many more reads than
the one or few Sanger reads behind any single clone, and snaps it onto its
own motif (find_insertion_sites.py's snap_window), so it's the more
trustworthy of the two, and two clones truly at the same integration will
converge on the identical NGS coordinate regardless of a few bp of Sanger
noise. Deliberately not also merging by raw distance on top of that: unlike
an NGS-backed match, a same-ish Sanger coordinate alone isn't independent
evidence two clones are the same integration - a local re-mobilisation
assay can pack real, distinct integrations within a few bp of each other
(see junction_jitter in find_peaks_junctions.py), so blindly collapsing
anything nearby risks conflating them instead.
"""

import argparse

import bioframe
import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--sites", required=True, help="all_sanger_sites.bed")
argparser.add_argument(
    "--clone-summary", required=True, help="sanger_clone_summary.tsv"
)
argparser.add_argument(
    "--chromsizes",
    default=None,
    help="Chrom sizes of the genome alone (chrom_sizes_path_no_cassette). "
    "Sites outside it - i.e. on the cassette/landing-pad contigs, which "
    "genome browsers don't know about - are dropped from --output-for-ucsc. "
    "Left unfiltered if not given.",
)
argparser.add_argument(
    "--region",
    default=None,
    help="chrom:start-end (0-based, half-open) - only validated sites "
    "overlapping this region go into --output-region. Unset keeps every "
    "validated site there too.",
)
argparser.add_argument(
    "--sanger-vs-ngs",
    default=None,
    help="sanger_vs_ngs.tsv, if NGS data is available - re-pinpoints a "
    "confirmed clone's coordinate onto the matching NGS site's own "
    "coordinate where compare_sanger_ngs.py found one. Unset leaves every "
    "clone at its own Sanger-only coordinate.",
)
argparser.add_argument("--output", "-o", required=True)
argparser.add_argument("--output-for-ucsc", required=True)
argparser.add_argument(
    "--output-no-cassette",
    required=True,
    help="Same sites as --output-for-ucsc (construct/landing-pad contigs "
    "dropped via --chromsizes), but a plain BED file with no track header.",
)
argparser.add_argument("--output-region", required=True)
argparser.add_argument(
    "--output-deduplicated",
    required=True,
    help="--output, collapsed to one row per distinct site - see deduplicate_by_position.",
)
argparser.add_argument(
    "--output-region-deduplicated",
    required=True,
    help="--output-region, collapsed the same way as --output-deduplicated.",
)
args = argparser.parse_args()

BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]
JOIN_COLUMNS = ["sample_name", "clone"]
POSITION_COLUMNS = ["chrom", "start", "end", "strand"]


def write_empty(path):
    empty = pd.DataFrame({c: pd.Series(dtype=object) for c in BED_COLUMNS})
    tagmaplib.write_bed(empty, path, columns=BED_COLUMNS)


def pinpoint_with_ngs(confirmed, sanger_vs_ngs_path):
    """Replace a confirmed clone's own Sanger coordinate with the matching
    NGS site's coordinate, wherever compare_sanger_ngs.py found one
    (confirmed_by_ngs, with ngs_site_chrom/start/end carrying that site's
    own position) - see the module docstring for why that one's preferred.
    A clone with no NGS match (no NGS data configured, or none nearby) is
    left exactly as Sanger reported it.
    """
    if sanger_vs_ngs_path is None:
        return confirmed
    ngs_match = pd.read_csv(
        sanger_vs_ngs_path, sep="\t", dtype={"chrom": str, "ngs_site_chrom": str}
    )
    ngs_match["confirmed_by_ngs"] = tagmaplib.to_bool(ngs_match["confirmed_by_ngs"])
    match_columns = JOIN_COLUMNS + POSITION_COLUMNS[:3]  # chrom, start, end
    ngs_columns = ["confirmed_by_ngs", "ngs_site_chrom", "ngs_site_start", "ngs_site_end"]
    confirmed = confirmed.merge(
        ngs_match[match_columns + ngs_columns], on=match_columns, how="left"
    )
    pinpointed = confirmed["confirmed_by_ngs"].fillna(False) & confirmed[
        "ngs_site_start"
    ].notna()
    for column, ngs_column in zip(
        POSITION_COLUMNS[:3], ["ngs_site_chrom", "ngs_site_start", "ngs_site_end"]
    ):
        confirmed.loc[pinpointed, column] = confirmed.loc[pinpointed, ngs_column]
    confirmed[["start", "end"]] = confirmed[["start", "end"]].astype(int)
    return confirmed.drop(columns=ngs_columns)


def deduplicate_by_position(df):
    """One row per distinct (chrom, start, end, strand) site, keeping the
    alphabetically-first sample_name+clone as the row's own name and
    appending how many other validated clones share that exact site - so a
    hotspot several clones landed on independently collapses to one row
    instead of stacking several identical-looking features on top of each
    other in a genome browser. An exact match only - see pinpoint_with_ngs
    for why two clones at the same true integration should already agree
    down to the base once NGS data can correct Sanger's own noise, and the
    module docstring for why this doesn't also collapse merely nearby sites
    on top of that.
    """
    if df.shape[0] == 0:
        return df
    ordered = df.sort_values(POSITION_COLUMNS + ["sample_name", "clone"])
    counts = ordered.groupby(POSITION_COLUMNS).size().rename("n_clones").reset_index()
    deduped = ordered.drop_duplicates(POSITION_COLUMNS, keep="first").merge(
        counts, on=POSITION_COLUMNS
    )
    n_others = deduped["n_clones"] - 1
    deduped["name"] = deduped["name"] + n_others.map(
        lambda n: f" (+{n} other{'s' if n != 1 else ''})" if n > 0 else ""
    )
    return deduped.sort_values(["chrom", "start", "end"])


sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str})
clones = pd.read_csv(args.clone_summary, sep="\t")

if sites.shape[0] == 0 or clones.shape[0] == 0:
    print("No Sanger sites to filter")
    write_empty(args.output)
    write_empty(args.output_for_ucsc)
    write_empty(args.output_no_cassette)
    write_empty(args.output_region)
    write_empty(args.output_deduplicated)
    write_empty(args.output_region_deduplicated)
    raise SystemExit(0)

has_ngs = "ngs_verification" in clones.columns
clones["color"] = [
    tagmaplib.clone_color(
        row.forward_status,
        row.reverse_status,
        row.both_sides_confirmed,
        row.ngs_verification if has_ngs else None,
    )
    for row in clones.itertuples()
]
validated_clones = clones.loc[clones["color"] == tagmaplib.CLONE_GREEN, JOIN_COLUMNS]

confirmed = sites.merge(validated_clones, on=JOIN_COLUMNS, how="inner")
confirmed = pinpoint_with_ngs(confirmed, args.sanger_vs_ngs)

# sample_name alone repeats across every clone on a plate; plate+clone
# actually identifies which well a browser hit came from.
confirmed["name"] = confirmed["sample_name"] + "_" + confirmed["clone"]
confirmed = confirmed.sort_values(["sample_name", "clone", "chrom", "start"])
tagmaplib.write_bed(confirmed, args.output, columns=BED_COLUMNS)

no_cassette = confirmed[BED_COLUMNS].sort_values(["chrom", "start", "end"])
if args.chromsizes is not None:
    chromsizes = bioframe.read_chromsizes(args.chromsizes)
    no_cassette = bioframe.trim(no_cassette, chromsizes).dropna()
    no_cassette[["start", "end"]] = no_cassette[["start", "end"]].astype(int)
tagmaplib.write_bed(no_cassette, args.output_no_cassette, columns=BED_COLUMNS)

track_name = (
    "confirmed_" + "_".join(sorted(confirmed["sample_name"].unique()))
    if confirmed.shape[0]
    else None
)
tagmaplib.write_bed(
    no_cassette, args.output_for_ucsc, track_name=track_name, color_by_strand=True
)

in_region = bioframe.select(confirmed, args.region) if args.region else confirmed
tagmaplib.write_bed(in_region, args.output_region, columns=BED_COLUMNS)

tagmaplib.write_bed(
    deduplicate_by_position(confirmed), args.output_deduplicated, columns=BED_COLUMNS
)
tagmaplib.write_bed(
    deduplicate_by_position(in_region),
    args.output_region_deduplicated,
    columns=BED_COLUMNS,
)

print(
    f"{confirmed.shape[0]} of {sites.shape[0]} sites fully validated"
    + (f", {in_region.shape[0]} in {args.region}" if args.region else "")
)
