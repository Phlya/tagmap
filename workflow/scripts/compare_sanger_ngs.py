"""Check Sanger integration sites against NGS data from the same material.

The two assays fail in different ways. A Sanger read gives one clean junction
but only from whichever ITR primer happened to work, while an NGS library
covers both sides but has to be pulled out of a pile of reads. So a Sanger site
that was only read from one end can still be trusted if an NGS site sits at the
same place, especially a two-sided one - which is what this reports.
"""

import argparse

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


if __name__ == "__main__":
    args = argparser.parse_args()

    pairs = dict(item.split("=", 1) for item in args.sample_pairs)

    sanger = pd.read_csv(args.sanger, sep="\t", dtype={"chrom": str})
    ngs_sites = pd.read_csv(args.ngs_sites, sep="\t", dtype={"chrom": str})
    ngs_peaks = tagmaplib.read_peaks(args.ngs_peaks)

    site_columns = ["strand", "score"]
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
