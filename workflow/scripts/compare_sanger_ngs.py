"""Check Sanger integration sites against NGS data from the same material.

The two assays fail in different ways. A Sanger read gives one clean junction
but only from whichever ITR primer happened to work, while an NGS library
covers both sides but has to be pulled out of a pile of reads. So a Sanger site
that was only read from one end can still be trusted if an NGS site sits at the
same place, especially a two-sided one - which is what this reports. The
matching NGS site's own coordinate is carried along too (ngs_site_chrom/
start/end), for pinpointing a confirmed Sanger site onto it later.
"""

import argparse
import json

import bioframe
import numpy as np
import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--sanger", required=True, help="all_sanger_sites.bed")
argparser.add_argument("--ngs-sites", required=True, help="all_sites.bed")
argparser.add_argument("--ngs-peaks", required=True, help="all_peaks.bed")
argparser.add_argument("--max-dist", type=int, default=500)
argparser.add_argument(
    "--ambiguity-dist",
    type=int,
    default=50,
    help="When several NGS sites lie within this distance of a Sanger site and "
    "they disagree on orientation or on which sides they show, the Sanger site "
    "cannot be tied to one of them, and its ngs_verification is 'ambiguous'.",
)
argparser.add_argument(
    "--original-site",
    default=None,
    help="original_site.json, if configured. Sanger sites there are never "
    "called ambiguous: the original locus sits in a repeat that NGS resolves "
    "into several mutually inconsistent sites, and is not a new integration.",
)
argparser.add_argument("--original-max-dist", type=int, default=0)
argparser.add_argument(
    "--sample-pairs",
    nargs="*",
    default=[],
    help="sanger_sample=ngs_sample assignments. A Sanger sample without one is "
    "compared against every NGS sample.",
)
argparser.add_argument("--output", "-o", required=True)
argparser.add_argument("--output-confirmed", required=True)


def nearest(sanger, ngs, prefix, extra_columns, max_dist):
    """Nearest NGS feature for each Sanger site, as a frame of new columns."""
    columns = [f"{prefix}_dist", f"{prefix}_sample"] + [
        f"{prefix}_{c}" for c in extra_columns
    ]
    if sanger.shape[0] == 0 or ngs.shape[0] == 0:
        return pd.DataFrame(
            {c: pd.Series([pd.NA] * sanger.shape[0], index=sanger.index) for c in columns}
        )

    # bioframe.closest returns its hits grouped by chromosome rather than in
    # the order it was handed the query, and drops queries on a chromosome the
    # other frame does not have. Taking the rows as they come and stamping the
    # Sanger index onto them therefore pairs most sites with some other site's
    # nearest neighbour, so carry the row each hit belongs to through the call.
    query = sanger[["chrom", "start", "end"]].copy()
    query["sanger_row"] = np.arange(sanger.shape[0])
    # A query on a chrom missing from `ngs` entirely (e.g. a Sanger site on a
    # short standalone contig NGS never reaches) makes bioframe fill that
    # row's columns with NA - which crashes on a plain bool column (e.g.
    # all_sites.bed's TA_found) rather than silently upcasting, so cast
    # defensively first.
    ngs = ngs.copy()
    bool_columns = ngs.columns[ngs.dtypes == bool]
    ngs[bool_columns] = ngs[bool_columns].astype(object)
    closest = (
        bioframe.closest(query, ngs, suffixes=("", "_ngs"), k=1)
        .set_index("sanger_row")
        .reindex(np.arange(sanger.shape[0]))
    )
    closest.index = sanger.index

    out = pd.DataFrame(index=sanger.index)
    out[f"{prefix}_dist"] = closest["distance"]
    out[f"{prefix}_sample"] = closest["sample_name_ngs"]
    for column in extra_columns:
        out[f"{prefix}_{column}"] = closest[f"{column}_ngs"]
    # Anything further away than the tolerance is not the same site.
    too_far = out[f"{prefix}_dist"] > max_dist
    out.loc[too_far, :] = pd.NA
    return out


def neighbourhood(sanger, ngs, ambiguity_dist):
    """How many NGS sites lie near each Sanger site, and whether they disagree.

    The nearest NGS site is not necessarily the one a Sanger read belongs to:
    a few bp of mapping noise can make it land closer to a neighbouring site
    with a different orientation or side pattern. So every NGS site within
    ambiguity_dist is compared, and the Sanger site is ambiguous if they do
    not all share one (strand, sides) signature.
    """
    n_nearby = np.zeros(sanger.shape[0], dtype=int)
    ambiguous = np.zeros(sanger.shape[0], dtype=bool)
    if sanger.shape[0] == 0 or ngs.shape[0] == 0:
        return pd.DataFrame(
            {"ngs_n_sites_nearby": n_nearby, "ngs_ambiguous": ambiguous},
            index=sanger.index,
        )
    ngs_by_chrom = {chrom: group for chrom, group in ngs.groupby("chrom")}
    for i, row in enumerate(sanger.itertuples()):
        candidates = ngs_by_chrom.get(row.chrom)
        if candidates is None:
            continue
        gap = np.maximum(
            np.maximum(candidates["start"].to_numpy() - row.end, row.start - candidates["end"].to_numpy()),
            0,
        )
        nearby = candidates[gap <= ambiguity_dist]
        n_nearby[i] = nearby.shape[0]
        # all_sites.bed from before site_sides existed has no sides column.
        sides = nearby["sides"] if "sides" in nearby.columns else [None] * len(nearby)
        signatures = set(zip(nearby["strand"], sides))
        ambiguous[i] = len(signatures) > 1
    return pd.DataFrame(
        {"ngs_n_sites_nearby": n_nearby, "ngs_ambiguous": ambiguous},
        index=sanger.index,
    )


if __name__ == "__main__":
    args = argparser.parse_args()

    pairs = dict(item.split("=", 1) for item in args.sample_pairs)

    sanger = pd.read_csv(args.sanger, sep="\t", dtype={"chrom": str})
    ngs_sites = pd.read_csv(args.ngs_sites, sep="\t", dtype={"chrom": str})
    if "site_sides" in ngs_sites.columns:
        # Renamed so the nearest() lookup below produces "ngs_site_sides"
        # rather than the doubled-up "ngs_site_site_sides".
        ngs_sites = ngs_sites.rename(columns={"site_sides": "sides"})
    ngs_peaks = tagmaplib.read_peaks(args.ngs_peaks)

    # chrom/start/end travel through too (as ngs_site_chrom/start/end) so a
    # confirmed Sanger site can be re-pinpointed onto the matching NGS site's
    # own coordinate later (see filter_confirmed_sites.py) - NGS backs a site
    # with many more reads than the one or few Sanger reads behind any single
    # clone, and snaps it onto its motif (see find_insertion_sites.py's
    # snap_window), so its coordinate is the more trustworthy of the two.
    site_columns = ["chrom", "start", "end", "strand", "score", "sides"]
    peak_columns = ["side", "count"]

    if sanger.shape[0] == 0:
        print("No Sanger sites to validate")
        sanger.to_csv(args.output, sep="\t", index=False)
        sanger.head(0).to_csv(
            args.output_confirmed, sep="\t", index=False, header=False
        )
        raise SystemExit(0)

    ngs_sites[["start", "end"]] = ngs_sites[["start", "end"]].astype(int)
    if ngs_peaks.shape[0]:
        ngs_peaks[["start", "end"]] = ngs_peaks[["start", "end"]].astype(int)

    annotated = []
    for sample, group in sanger.groupby("sample_name", sort=True):
        # Only compare against the matching library when one was named, so that
        # an unrelated clone's insertion cannot vouch for this one.
        ngs_sample = pairs.get(sample)
        if ngs_sample is None:
            sites_subset, peaks_subset = ngs_sites, ngs_peaks
        else:
            sites_subset = ngs_sites[ngs_sites["sample_name"] == ngs_sample]
            peaks_subset = ngs_peaks[ngs_peaks["sample_name"] == ngs_sample]
        group = group.copy()
        group["ngs_sample_used"] = ngs_sample if ngs_sample else "all"
        group = pd.concat(
            [
                group,
                nearest(group, sites_subset, "ngs_site", site_columns, args.max_dist),
                nearest(group, peaks_subset, "ngs_peak", peak_columns, args.max_dist),
                neighbourhood(group, sites_subset, args.ambiguity_dist),
            ],
            axis=1,
        )
        annotated.append(group)

    sanger = pd.concat(annotated).sort_values(["sample_name", "clone", "chrom", "start"])

    matched = sanger["ngs_site_dist"].notna()
    # A '.' strand means the NGS side could not be decided, so it neither
    # confirms nor contradicts the Sanger strand.
    strand_known = matched & sanger["ngs_site_strand"].isin(["+", "-"])
    sanger["ngs_strand_agrees"] = np.where(
        strand_known, sanger["ngs_site_strand"] == sanger["strand"], pd.NA
    )
    # find_insertion_sites only assigns a strand when it saw both sides of the
    # cassette, so a stranded NGS site is a two-sided one.
    sanger["ngs_two_sided"] = np.where(matched, strand_known, pd.NA)
    sanger["confirmed_by_ngs"] = matched & (sanger["ngs_strand_agrees"] != False)
    # Only a matched site can be ambiguous: with nothing within max_dist there
    # is no verification to cast doubt on (that is "no NGS site nearby").
    sanger["ngs_ambiguous"] = sanger["ngs_ambiguous"].astype(bool) & matched
    if args.original_site is not None:
        with open(args.original_site) as f:
            original_site = json.load(f)
        at_original_site = (sanger["chrom"] == original_site["chrom"]) & (
            (sanger["start"] - original_site["pos"]).abs() <= args.original_max_dist
        )
        sanger["ngs_ambiguous"] &= ~at_original_site

    def row_ngs_verification(row):
        ngs_sides = set(
            tagmaplib.NGS_SIDE_SETS.get(row.ngs_site_sides, [])
            if row.confirmed_by_ngs
            else []
        )
        # Lacking a called site, a single NGS peak sitting at the exact same
        # base (not just "nearby") is still real, independent evidence - see
        # ngs_verification_label - just from a molecule that never cleared
        # (or was never orientation-resolved into) a full site call.
        exact_peak = pd.notna(row.ngs_peak_dist) and row.ngs_peak_dist == 0
        if exact_peak:
            direction = tagmaplib.PEAK_SIDE_DIRECTIONS.get(row.ngs_peak_side)
            if direction:
                ngs_sides.add(direction)
        return tagmaplib.ngs_verification_label(
            ngs_sides,
            ngs_checked=pd.notna(row.ngs_site_dist) or exact_peak,
            ambiguous=bool(row.ngs_ambiguous),
        )

    # What NGS alone says here (see tagmaplib.ngs_verification_label); which
    # primer(s) Sanger confirmed is in n_forward/n_reverse. "ambiguous" when
    # the NGS sites around the position disagree with each other.
    sanger["ngs_verification"] = [
        row_ngs_verification(row) for row in sanger.itertuples()
    ]

    sanger.to_csv(args.output, sep="\t", index=False)

    confirmed = sanger[sanger["confirmed_by_ngs"]]
    confirmed[tagmaplib.SITE_COLUMNS].sort_values(["chrom", "start", "end"]).to_csv(
        args.output_confirmed, sep="\t", index=False, header=False
    )

    print(
        f"{sanger.shape[0]} Sanger sites, {int(matched.sum())} with a matching "
        f"NGS site, {int(sanger['confirmed_by_ngs'].sum())} confirmed, "
        f"{int((sanger['ngs_two_sided'] == True).sum())} of those two-sided in NGS"
    )
