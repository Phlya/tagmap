"""Assemble the pipeline's summary tables into one markdown report.

SECTIONS is also imported by report_pdf.py, so the two reports share the
same section titles and descriptions.
"""

import argparse
import html
import os

import pandas as pd

import tagmaplib

SECTIONS = [
    (
        "ngs_qc",
        "NGS QC",
        "Per library: read pairs anchored at an ITR primer (evidence of "
        "mobilization), and how many of the resulting insertion sites were "
        "seen from both sides of the cassette versus only one - a one-sided "
        "site is weaker evidence, since it has not been confirmed by an "
        "independent primer.",
    ),
    (
        "sanger_qc",
        "Sanger QC",
        "Sanger reads per run and ITR primer direction: how many passed QC, "
        "and why the rest didn't. n_clones/n_clones_pass/frac_clones_pass "
        "give the same breakdown by clone instead of by read - a clone with "
        "several reads on one side (a rerun, or a differently-primed retry "
        "on the same well) counts once, as soon as any of them passed - so "
        "this is the number to read for \"how many clones did this side "
        "work for\", as distinct from n_pass's raw read count.",
    ),
    (
        "sanger_clones",
        "Sanger clones",
        "For each clone, whether the integration site was confirmed from the "
        "forward primer, the reverse primer, both, or neither - and why an "
        "unconfirmed side isn't: never sequenced, or how it failed QC. "
        "position lists the site(s) clustered from that clone's passing "
        "reads (all_sanger_sites.bed), each with its own forward/reverse "
        "read breakdown - usually one site backed by both primers, but a "
        "clone whose forward and reverse reads don't cluster together "
        "shows both, so you can see where they disagree and which primer "
        "supports which locus. orientation pulls the strand(s) already "
        "shown in position out into their own column - which way the "
        "insertion faces, one per site listed there, in the same order; "
        "\".\" marks the rare site where forward and reverse reads land on "
        "the same locus but disagree on which way it faces. Passing QC on "
        "both sides isn't enough on its own: positions_agree checks that "
        "the forward and reverse reads actually cluster onto the same site "
        "rather than "
        "pointing at different loci, and both_sides_confirmed requires "
        "both. Where NGS data for the same material is available, "
        "ngs_verification gives the combined Sanger+NGS call for that "
        "clone's site - \"both\", \"forward only\", \"reverse only\", or "
        "\"not verified\" - using NGS only to fill in a side Sanger's own "
        "primers left unconfirmed, never to override a side Sanger already "
        "confirmed (so a clone Sanger already confirmed from both primers "
        "reads \"both\" regardless of what NGS shows); \"disagreeing "
        "(NGS consistent/inconsistent)\" if the clone's Sanger reads "
        "themselves disagree on locus - disagreement alone already fails "
        "the clone regardless of NGS, so this just reports, for "
        "information only, whether NGS's own side-calling matches which "
        "single primer found each disagreeing site (the original, "
        "pre-mobilization locus is exempt from that check - NGS showing "
        "both sides there is the ordinary, expected result). Where the "
        "original, pre-mobilization insertion "
        "is configured, a clone whose both-sides-agreeing site sits there "
        "instead of a new locus reads \"unmobilized\" in summary rather "
        "than \"both sides confirmed\".",
    ),
    (
        "sanger_positions",
        "Sanger positions",
        "Distinct genomic positions found across the clustered Sanger sites "
        "(all_sanger_sites.bed), per plate and in total - two clones that "
        "land on the very same locus, such as an unmobilized founder control "
        "sequenced more than once, count once rather than twice. n_one_sided "
        "counts how many of those positions have no clone confirming them "
        "from both ITR primers at all - weaker evidence than a clone (or two "
        "agreeing clones) confirmed from both sides would give. Where NGS "
        "data for the same material is available, n_ngs_confirmed separately "
        "counts how many positions were independently backed by NGS, whether "
        "or not they're also one-sided in Sanger.",
    ),
    (
        "sanger_position_counts",
        "Sanger position counts",
        "Every distinct new-integration position from the table above, with "
        "how many clones landed on it - most-found first - and how many of "
        "those were confirmed by only one primer (n_one_sided): weaker "
        "evidence, since the clones just below show that a primer stuck on "
        "the founder locus, when its partner primer is usable at all, often "
        "turns out to disagree and land elsewhere - so a one-sided clone "
        "could equally be a mobilized clone whose other primer simply "
        "failed. A mobilization hotspot stands out here rather than being "
        "folded into a single count. Clones at the configured original "
        "insertion site - one-sided or both-sides-confirmed alike, so this "
        "row's n_times_found can run well ahead of the stricter \"both "
        "primers confirm it\" count noted under the Sanger clones table "
        "above - are pulled out of the per-locus breakdown and tallied in "
        "their own combined row instead (otherwise they'd fragment across "
        "several near-identical founder-locus "
        "coordinates), and likewise clones whose forward and reverse reads "
        "disagree on the locus (e.g. one stuck on a residual copy of the "
        "donor construct) get their own combined row - n_one_sided is N/A "
        "for that row, since both primers there were confirmed, just on "
        "different sites, the opposite situation. Where NGS data is "
        "available, ngs_verification gives the same combined Sanger+NGS "
        "call as the clone summary's own column, pooled across every clone "
        "sharing that row instead of one clone at a time; N/A for the "
        "disagreeing row and the two totals below. The last two rows are "
        "clones with at least one good read (every clone with at least one "
        "QC-passing, motif-snapped site - i.e. represented somewhere above) "
        "and total mobilized positions (the count of real position rows only, "
        "excluding both label rows).",
    ),
    (
        "validation",
        "Sanger vs NGS validation",
        "Sanger sites cross-checked against the NGS data from the same "
        "material, run by run and in total.",
    ),
]

WELL_BORDER = "0.2mm solid #999"
WELL_HIGHLIGHT_BORDER = "0.6mm solid #000"
# Two diagonal bars layered as background images, clipped to the well's own
# border-radius:50% - a CSS-only cross, since inline HTML here has no
# stylesheet to hang a ::before/::after pseudo-element off of.
CROSS_BACKGROUND = (
    "linear-gradient(45deg, transparent 46%, #000 46%, #000 54%, transparent 54%),"
    "linear-gradient(-45deg, transparent 46%, #000 46%, #000 54%, transparent 54%)"
)


def retained_duplicate_names(path):
    """{sample_name}_{clone} names from a deduplicated sites BED (see
    filter_confirmed_sites.py's deduplicate_by_position) - the clone kept to
    represent its coordinate, one entry per row, whether that coordinate had
    one clone or several landing on it. A name with a "(+N other(s))" suffix
    is stripped back down to the bare name first, so it can be matched
    against a clone_summary row's own sample_name/clone. Empty if path is
    None, or the file has no rows (Sanger configured but nothing validated
    yet).
    """
    if path is None or os.path.getsize(path) == 0:
        return set()
    names = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )["name"]
    return {name.split(" (+")[0] for name in names}


def well_grid_html(wells, n_rows, n_cols):
    """The <table> for one color-coded well grid, as raw HTML (markdown
    files pass inline HTML through untouched). `wells` maps (row, col) ->
    (title, color, marker, highlight, duplicate) - shared by plate_html
    (clone status) and read_qc_plate_html (per-well read QC status), which
    only differ in what they put there. marker is a short string (e.g. "u")
    drawn centered in the well, or None for no marker; highlight draws a
    bold ring around the well instead of the usual thin one; duplicate
    draws a black cross over it.

    Cells are WELL_PITCH_MM square - the real plate's well-to-well spacing -
    and each holds a smaller WELL_DIAMETER_MM circle, so the grid is to
    scale with a real plate rather than just a same-shaped table.
    """
    pitch = tagmaplib.WELL_PITCH_MM
    cell_style = f"width:{pitch}mm;height:{pitch}mm;text-align:center;vertical-align:middle;padding:0"
    well_style = (
        f"display:inline-flex;align-items:center;justify-content:center;"
        f"width:{tagmaplib.WELL_DIAMETER_MM}mm;height:{tagmaplib.WELL_DIAMETER_MM}mm;"
        f"border-radius:50%;font-size:3mm;color:#000"
    )

    lines = [
        '<table style="border-collapse:collapse">',
        "<tr><th></th>"
        + "".join(f'<th style="{cell_style}">{c + 1}</th>' for c in range(n_cols))
        + "</tr>",
    ]
    plate_border = "0.4mm solid #333"
    for r in range(n_rows):
        cells = [f'<th style="{cell_style}">{chr(ord("A") + r)}</th>']
        for c in range(n_cols):
            # The plate's own outline, around the wells only - not the
            # row/column labels alongside them - traced by bordering just
            # the outer edge of well cells.
            edges = ";".join(
                filter(
                    None,
                    [
                        f"border-top:{plate_border}" if r == 0 else "",
                        f"border-bottom:{plate_border}" if r == n_rows - 1 else "",
                        f"border-left:{plate_border}" if c == 0 else "",
                        f"border-right:{plate_border}" if c == n_cols - 1 else "",
                    ],
                )
            )
            entry = wells.get((r, c))
            if entry is None:
                cells.append(
                    f'<td style="{cell_style};{edges}">'
                    f'<span style="{well_style};border:{WELL_BORDER};background:#f5f5f5"></span></td>'
                )
            else:
                title, color, marker, highlight, duplicate = entry
                hexcolor = "#%02x%02x%02x" % color
                border = WELL_HIGHLIGHT_BORDER if highlight else WELL_BORDER
                background = f"background-color:{hexcolor}" + (
                    f";background-image:{CROSS_BACKGROUND}" if duplicate else ""
                )
                cells.append(
                    f'<td style="{cell_style};{edges}" title="{html.escape(str(title))}">'
                    f'<span style="{well_style};border:{border};{background}">'
                    f'{html.escape(marker) if marker else ""}</span></td>'
                )
        lines.append("<tr>" + "".join(cells) + "</tr>")
    lines.append("</table>")
    return "\n".join(lines)


def plate_html(sample_name, group, dedup_names=frozenset()):
    """A color-coded well grid for one plate, as raw HTML. Same colours as
    the PDF's clone table rows - see tagmaplib.clone_color(). None if this
    plate's clone ids don't look like well ids at all, since not every
    project's do. A clone whose confirmed site is the configured original
    insertion site (summary == "unmobilized") gets a small "u" marker. A
    validated clone (green) landing inside validated_clones_region gets a
    bold ring, same criterion as report_pdf.py's table row border. A clone
    in dedup_names (see retained_duplicate_names) gets a black cross: it's
    the clone kept to represent its integration coordinate in
    confirmed_sanger_sites_region_deduplicated.bed, whether that coordinate
    had just this one clone or several landing on it independently - see
    that file's own name column for which, and how many others if any.
    """
    wells = {}
    for row in group.itertuples():
        position = tagmaplib.parse_well(row.clone)
        if position is None:
            continue
        color = tagmaplib.clone_color(
            row.forward_status,
            row.reverse_status,
            row.both_sides_confirmed,
            getattr(row, "ngs_verification", None),
        )
        wells[position] = (
            row.clone,
            color,
            "u" if getattr(row, "summary", None) == "unmobilized" else None,
            color == tagmaplib.CLONE_GREEN and bool(getattr(row, "in_region", False)),
            f"{row.sample_name}_{row.clone}" in dedup_names,
        )
    if not wells:
        return None
    n_rows, n_cols = tagmaplib.plate_layout(wells.keys())
    return "\n".join(
        [
            f"**{html.escape(tagmaplib.plate_label(sample_name, group))}**",
            "",
            well_grid_html(wells, n_rows, n_cols),
        ]
    )


def read_qc_plate_html(sample_name, group):
    """Forward and reverse color-coded well grids for one plate's read QC
    (tagmaplib.read_qc_status), side by side on the same page so a well's
    two sides can be compared directly. None if this plate's well ids don't
    look like well ids at all - see plate_html. A well whose passing read(s)
    on that side land at the configured original insertion site
    (at_original_site) get a small "u" marker.
    """
    grids = []
    for direction in ("forward", "reverse"):
        wells = {}
        for row in group[group["direction"] == direction].itertuples():
            position = tagmaplib.parse_well(row.well)
            if position is None:
                continue
            wells[position] = (
                row.well,
                tagmaplib.READ_QC_COLORS[row.status],
                "u" if getattr(row, "at_original_site", False) else None,
                False,
                False,
            )
        if not wells:
            continue
        n_rows, n_cols = tagmaplib.plate_layout(wells.keys())
        grids.append(
            f"<div><strong>{html.escape(direction.capitalize())}</strong><br>"
            + well_grid_html(wells, n_rows, n_cols)
            + "</div>"
        )
    if not grids:
        return None
    return "\n".join(
        [
            f"**{html.escape(tagmaplib.plate_label(sample_name, group))}**",
            "",
            '<div style="display:flex;gap:24px;flex-wrap:wrap;align-items:flex-start">'
            + "".join(grids)
            + "</div>",
        ]
    )


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument("--ngs-qc", default=None)
    argparser.add_argument("--sanger-qc", default=None)
    argparser.add_argument("--sanger-clones", default=None)
    argparser.add_argument("--sanger-positions", default=None)
    argparser.add_argument("--sanger-position-counts", default=None)
    argparser.add_argument("--sanger-read-qc", default=None)
    argparser.add_argument(
        "--deduplicated",
        default=None,
        help="confirmed_sanger_sites_region_deduplicated.bed, to mark the "
        "clone kept to represent each integration coordinate in there with "
        "a black cross in the plate grid below.",
    )
    argparser.add_argument("--validation", default=None)
    argparser.add_argument("--output", "-o", required=True)
    args = argparser.parse_args()

    dedup_names = retained_duplicate_names(args.deduplicated)

    lines = ["# TagMap summary", ""]

    for attr, title, intro in SECTIONS:
        path = getattr(args, attr)
        if path is None:
            continue
        df = tagmaplib.tidy_numeric_dtypes(pd.read_csv(path, sep="\t"))
        lines.append(f"## {title}")
        lines.append("")
        lines.append(intro)
        lines.append("")
        if df.shape[0] == 0:
            lines.append("_No data._")
        else:
            lines.append(tagmaplib.to_markdown(df))
        lines.append("")

        if attr == "sanger_clones" and df.shape[0]:
            n_unmobilized = int((df["summary"] == "unmobilized").sum())
            lines.append(
                f"**{n_unmobilized} of {df.shape[0]} clones unmobilized** "
                "(both primers independently confirm the original insertion "
                "site). A clone confirmed from only one primer there doesn't "
                "count here even though that locus is still its best "
                "evidence - see the \"unmobilized: at the original insertion "
                "site\" row of the Sanger position counts table below for "
                "those too. A small \"u\" in the plate grid below marks a "
                "clone that landed there."
                + (
                    " A well drawn with a bold ring is validated (green) and "
                    "landed inside validated_clones_region."
                    if "in_region" in df.columns
                    else ""
                )
                + (
                    " A well marked with a black cross is the clone kept to "
                    "represent an integration coordinate within "
                    "validated_clones_region - see "
                    "confirmed_sanger_sites_region_deduplicated.bed for "
                    "whether any other validated clone independently landed "
                    "on the same coordinate."
                    if dedup_names
                    else ""
                )
            )
            lines.append("")
            for sample_name, group in df.groupby("sample_name"):
                plate = plate_html(sample_name, group, dedup_names)
                if plate is not None:
                    lines.append(plate)
                    lines.append("")

        if attr == "sanger_qc" and args.sanger_read_qc:
            read_qc = pd.read_csv(args.sanger_read_qc, sep="\t")
            if read_qc.shape[0]:
                lines.append("### Sequencing success by well")
                lines.append("")
                lines.append(
                    "Forward and reverse ITR primer reads, by well: green - "
                    f"worked ({tagmaplib.READ_QC_SHORT_MAX_LENGTH}bp+ of clean, "
                    "trimmed sequence); yellow - very short "
                    f"({tagmaplib.READ_QC_MIN_CLEAN_BASES}-"
                    f"{tagmaplib.READ_QC_SHORT_MAX_LENGTH - 1}bp of clean "
                    "sequence); red - failed (a trace was run, but came back "
                    f"under {tagmaplib.READ_QC_MIN_CLEAN_BASES} clean bases); "
                    "grey - not sequenced (no trace file for this well/side "
                    "at all). A small \"u\" marks a well whose passing "
                    "read(s) on that side land at the configured original "
                    "insertion site."
                )
                lines.append("")
                for sample_name, group in read_qc.groupby("sample_name"):
                    plate = read_qc_plate_html(sample_name, group)
                    if plate is not None:
                        lines.append(plate)
                        lines.append("")

    with open(args.output, "w") as f:
        f.write("\n".join(lines))
