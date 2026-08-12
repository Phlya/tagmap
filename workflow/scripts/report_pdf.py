"""Render the pipeline's summary tables as a single PDF.

Shares its section titles and descriptions with report.py's markdown report.
The Sanger clones table is additionally colour-coded per clone, and followed
by a per-plate well grid in the same colours - see tagmaplib.clone_color()
for the green/yellow/red rules.
"""

import argparse
import math

import pandas as pd
from fpdf import FPDF
from fpdf.enums import TableBorderStyle, TableBordersLayout, TableCellStyle
from fpdf.fonts import FontFace

import tagmaplib
from report import SECTIONS, retained_duplicate_names

WHITE = (255, 255, 255)
HEADING_GREY = (221, 221, 221)
EMPTY_WELL = (245, 245, 245)

MAX_COL_WEIGHT = 45
MIN_COL_WEIGHT = 6
WELL_LABEL = 6

REGION_BORDER_WIDTH = 0.6  # mm - fpdf2's own default hairline is ~0.2mm


class HighlightedRowsBorders(TableBordersLayout):
    """The same grid as TableBordersLayout.ALL, but rows in `bold_rows`
    (0-indexed including the heading row, so a data row's own index needs +1)
    get a thicker box drawn around the row as a whole, rather than around
    each of its cells - the row's own left/right edge is only bold at its
    first/last column, so the columns in between don't come out looking like
    a bold sub-grid of their own.

    A row's top/bottom edge is shared with its neighbour's bottom/top, drawn
    as two separate strokes at the same coordinates - whichever is drawn
    last is the one that actually ends up visible. Thickening that shared
    edge from *either* side, not just the highlighted row's own, means the
    box comes out bold regardless of draw order.
    """

    def __init__(self, bold_rows):
        self.bold_rows = bold_rows

    def cell_style_getter(
        self,
        row_idx,
        col_idx,
        col_pos,
        num_heading_rows,
        num_rows,
        num_col_idx,
        num_col_pos,
    ):
        thick = TableBorderStyle(thickness=REGION_BORDER_WIDTH)
        is_bold = row_idx in self.bold_rows
        return TableCellStyle(
            left=thick if is_bold and col_idx == 0 else True,
            right=thick if is_bold and col_idx == num_col_idx - 1 else True,
            top=thick if is_bold or (row_idx - 1) in self.bold_rows else True,
            bottom=thick if is_bold or (row_idx + 1) in self.bold_rows else True,
        )


def column_weights(df):
    """Relative column widths for pdf.table(), sized off each column's own
    content so long columns (e.g. summary) get proportionally more room."""
    weights = []
    for col in df.columns:
        values = df[col].astype(str)
        typical_len = values.str.len().quantile(0.9) if values.shape[0] else 0
        weights.append(min(max(len(col), typical_len, MIN_COL_WEIGHT), MAX_COL_WEIGHT))
    return weights


def add_table(pdf, df, row_colors=None, bold_rows=None):
    """bold_rows: 0-indexed data row positions (i.e. matching row_colors'
    own indexing, not counting the heading row) to draw a bold box around -
    see HighlightedRowsBorders.
    """
    # astype(object) first: fillna("") on a nullable Int64 column (produced
    # by tagmaplib.tidy_numeric_dtypes for a column with blank cells) raises,
    # since "" isn't a valid Int64 value - object dtype accepts either.
    df = df.astype(object).fillna("")
    # Format before sizing columns, so a long-precision float (e.g.
    # 0.8645833333333334) doesn't widen its column past what the rounded
    # text (0.86) actually needs.
    formatted = df.apply(lambda col: col.map(tagmaplib.format_cell))
    pdf.set_font("helvetica", size=6.5)
    borders_layout = (
        HighlightedRowsBorders({i + 1 for i in bold_rows})
        if bold_rows
        else TableBordersLayout.ALL
    )
    with pdf.table(
        col_widths=column_weights(formatted),
        text_align="LEFT",
        line_height=3.2,
        padding=1,
        headings_style=FontFace(emphasis="BOLD", fill_color=HEADING_GREY),
        borders_layout=borders_layout,
    ) as table:
        header = table.row()
        for col in df.columns:
            header.cell(col)
        for i, (_, data) in enumerate(formatted.iterrows()):
            row = table.row()
            style = FontFace(fill_color=row_colors[i] if row_colors is not None else WHITE)
            for col in formatted.columns:
                row.cell(data[col], style=style)


def add_legend(pdf, has_region=False, has_duplicates=False):
    pdf.set_font("helvetica", size=8)
    for color, label in (
        (tagmaplib.CLONE_GREEN, "both sides confirmed and agree, or only one side passed and NGS confirms both sides there"),
        (tagmaplib.CLONE_YELLOW, "only one side passed, and NGS did not confirm both sides there"),
        (tagmaplib.CLONE_RED, "neither side passed, or both passed but disagree on position"),
    ):
        pdf.set_fill_color(*color)
        pdf.cell(4, 4, "", border=1, fill=True)
        pdf.set_x(pdf.get_x() + 2)
        pdf.cell(0, 4, label, new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("helvetica", style="B", size=8)
    pdf.cell(4, 4, "u", align="C", border=1)
    pdf.set_font("helvetica", size=8)
    pdf.set_x(pdf.get_x() + 2)
    pdf.cell(
        0,
        4,
        "clone's confirmed site is the configured original insertion site",
        new_x="LMARGIN",
        new_y="NEXT",
    )
    if has_region:
        _remember_linewidth = pdf.line_width
        pdf.set_line_width(REGION_BORDER_WIDTH)
        pdf.cell(4, 4, "", border=1)
        pdf.set_line_width(_remember_linewidth)
        pdf.set_x(pdf.get_x() + 2)
        pdf.cell(
            0,
            4,
            "bold border/ring (table row, or well below) - validated, and "
            "inside validated_clones_region",
            new_x="LMARGIN",
            new_y="NEXT",
        )
    if has_duplicates:
        x, y = pdf.get_x(), pdf.get_y()
        pdf.cell(4, 4, "", border=1)
        _remember_linewidth = pdf.line_width
        pdf.set_line_width(REGION_BORDER_WIDTH)
        pdf.set_draw_color(0, 0, 0)
        pdf.line(x + 0.5, y + 0.5, x + 3.5, y + 3.5)
        pdf.line(x + 3.5, y + 0.5, x + 0.5, y + 3.5)
        pdf.set_line_width(_remember_linewidth)
        pdf.set_x(x + 4 + 2)
        pdf.cell(
            0,
            4,
            "black cross in well below - the clone kept to represent an "
            "integration coordinate in validated_clones_region (see "
            "confirmed_sanger_sites_region_deduplicated.bed for whether any "
            "other validated clone independently landed on the same "
            "coordinate)",
            new_x="LMARGIN",
            new_y="NEXT",
        )
    pdf.ln(2)


READ_QC_LEGEND_LABELS = {
    tagmaplib.READ_QC_WORKED: f"worked ({tagmaplib.READ_QC_SHORT_MAX_LENGTH}bp+ of clean, trimmed sequence)",
    tagmaplib.READ_QC_SHORT: (
        f"very short ({tagmaplib.READ_QC_MIN_CLEAN_BASES}-"
        f"{tagmaplib.READ_QC_SHORT_MAX_LENGTH - 1}bp of clean sequence)"
    ),
    tagmaplib.READ_QC_FAILED: (
        f"failed (a trace was run, but under {tagmaplib.READ_QC_MIN_CLEAN_BASES} clean bases)"
    ),
    tagmaplib.READ_QC_NOT_SEQUENCED: "not sequenced (no trace file for this well/side)",
}

READ_QC_LEGEND_ORDER = (
    tagmaplib.READ_QC_WORKED,
    tagmaplib.READ_QC_SHORT,
    tagmaplib.READ_QC_FAILED,
    tagmaplib.READ_QC_NOT_SEQUENCED,
)


def add_read_qc_legend(pdf):
    pdf.set_font("helvetica", size=8)
    for status in READ_QC_LEGEND_ORDER:
        pdf.set_fill_color(*tagmaplib.READ_QC_COLORS[status])
        pdf.cell(4, 4, "", border=1, fill=True)
        pdf.set_x(pdf.get_x() + 2)
        pdf.cell(0, 4, READ_QC_LEGEND_LABELS[status], new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("helvetica", style="B", size=8)
    pdf.cell(4, 4, "u", align="C", border=1)
    pdf.set_font("helvetica", size=8)
    pdf.set_x(pdf.get_x() + 2)
    pdf.cell(
        0,
        4,
        "passing read(s) on this side land at the configured original insertion site",
        new_x="LMARGIN",
        new_y="NEXT",
    )
    pdf.ln(2)


def draw_well_grid(pdf, x0, y0, wells, n_rows, n_cols):
    """The wells, row/column labels and outline of one plate grid at a given
    page position, drawn at WELL_PITCH/WELL_DIAMETER so it's to scale with a
    real plate rather than just a same-shaped table. Shared by draw_plate
    (clone status) and draw_read_qc_plate (per-well read QC status), which
    only differ in what `wells` maps each position to: (color, marker,
    highlight, duplicate), where marker is a short string (e.g. "u") drawn
    centered in the well, or None for no marker; highlight draws a bold ring
    around the well - see HighlightedRowsBorders for the same idea applied
    to a table row; duplicate draws a black cross over it.
    """
    pitch = tagmaplib.WELL_PITCH_MM
    diameter = tagmaplib.WELL_DIAMETER_MM
    inset = (pitch - diameter) / 2
    pdf.set_font("helvetica", size=6)
    for c in range(n_cols):
        pdf.set_xy(x0 + WELL_LABEL + c * pitch, y0)
        pdf.cell(pitch, WELL_LABEL, str(c + 1), align="C")
    for r in range(n_rows):
        y = y0 + WELL_LABEL + r * pitch
        pdf.set_xy(x0, y)
        pdf.cell(WELL_LABEL, pitch, chr(ord("A") + r), align="C")
        for c in range(n_cols):
            color, marker, highlight, duplicate = wells.get(
                (r, c), (EMPTY_WELL, None, False, False)
            )
            well_x = x0 + WELL_LABEL + c * pitch + inset
            well_y = y + inset
            pdf.set_fill_color(*color)
            pdf.ellipse(well_x, well_y, diameter, diameter, style="DF")
            if highlight:
                _remember_linewidth = pdf.line_width
                pdf.set_line_width(REGION_BORDER_WIDTH)
                pdf.ellipse(well_x, well_y, diameter, diameter, style="D")
                pdf.set_line_width(_remember_linewidth)
            if duplicate:
                # Inset to the largest square that fits inside the well's
                # own circle, so the cross's arms don't poke out past it.
                _remember_linewidth = pdf.line_width
                pdf.set_line_width(REGION_BORDER_WIDTH)
                pdf.set_draw_color(0, 0, 0)
                cx, cy = well_x + diameter / 2, well_y + diameter / 2
                half = diameter / 2 / math.sqrt(2)
                pdf.line(cx - half, cy - half, cx + half, cy + half)
                pdf.line(cx + half, cy - half, cx - half, cy + half)
                pdf.set_line_width(_remember_linewidth)
            if marker:
                pdf.set_font("helvetica", style="B", size=5)
                pdf.set_text_color(0, 0, 0)
                pdf.set_xy(x0 + WELL_LABEL + c * pitch, y)
                pdf.cell(pitch, pitch, marker, align="C")
                pdf.set_font("helvetica", size=6)
    # The plate's own outline, around the wells only - not the row/column
    # labels floating outside it.
    pdf.rect(x0 + WELL_LABEL, y0 + WELL_LABEL, n_cols * pitch, n_rows * pitch, style="D")


def draw_plate(pdf, sample_name, group, dedup_names=frozenset(), new_page=True):
    """One color-coded 96-well (or larger) grid for a single plate, in the
    same colours as the clone table rows. Wells whose clone id doesn't parse
    (see tagmaplib.parse_well) are simply left off the grid; if none parse,
    nothing is drawn at all - not every project's clone ids are well ids. A
    clone whose confirmed site is the configured original insertion site
    (summary == "unmobilized") gets a small "u" marker. A validated clone
    (green) landing inside validated_clones_region gets a bold ring, same
    criterion as the table's own bold row border. A clone in dedup_names
    (see report.retained_duplicate_names) gets a black cross: it's the
    clone kept to represent its integration coordinate in
    confirmed_sanger_sites_region_deduplicated.bed, whether that coordinate
    had just this one clone or several landing on it independently - see
    that file's own name column for which, and how many others if any.

    new_page=False packs this plate below whatever's already on the current
    page (the caller's way of fitting two plates per page) - but only if it
    actually fits there; draw_well_grid has no pagination logic of its own,
    so a plate too big to fit isn't given the choice and starts a fresh page
    regardless, rather than being clipped by the page bottom.

    Returns True if a plate was actually drawn, False if this plate's clone
    ids don't look like well ids at all - so callers pairing up plates two
    to a page can skip a group that used up no space.
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
            color,
            "u" if getattr(row, "summary", None) == "unmobilized" else None,
            color == tagmaplib.CLONE_GREEN and bool(getattr(row, "in_region", False)),
            f"{row.sample_name}_{row.clone}" in dedup_names,
        )
    if not wells:
        return False
    n_rows, n_cols = tagmaplib.plate_layout(wells.keys())
    block_height = 8 + 1 + WELL_LABEL + n_rows * tagmaplib.WELL_PITCH_MM

    if new_page or pdf.get_y() + block_height > pdf.page_break_trigger:
        pdf.add_page()
    pdf.set_font("helvetica", style="B", size=12)
    pdf.cell(0, 8, tagmaplib.plate_label(sample_name, group), new_x="LMARGIN", new_y="NEXT")
    pdf.ln(1)

    x0, y0 = pdf.get_x(), pdf.get_y()
    draw_well_grid(pdf, x0, y0, wells, n_rows, n_cols)
    pdf.set_xy(x0, y0 + WELL_LABEL + n_rows * tagmaplib.WELL_PITCH_MM + 4)
    return True


def draw_read_qc_plate(pdf, sample_name, group):
    """Forward and reverse color-coded well grids for one plate's read QC
    (tagmaplib.read_qc_status), side by side on one landscape page so a
    well's two sides can be compared directly. Wells whose id doesn't parse
    (see tagmaplib.parse_well) are left off the grid; if neither direction
    has any, nothing is drawn. A well whose passing read(s) on that side
    land at the configured original insertion site (at_original_site) get a
    small "u" marker.
    """
    direction_wells = {}
    for direction in ("forward", "reverse"):
        wells = {}
        for row in group[group["direction"] == direction].itertuples():
            position = tagmaplib.parse_well(row.well)
            if position is None:
                continue
            wells[position] = (
                tagmaplib.READ_QC_COLORS[row.status],
                "u" if getattr(row, "at_original_site", False) else None,
                False,
                False,
            )
        if wells:
            direction_wells[direction] = wells
    if not direction_wells:
        return

    pdf.add_page(orientation="L")
    pdf.set_font("helvetica", style="B", size=12)
    pdf.cell(0, 8, tagmaplib.plate_label(sample_name, group), new_x="LMARGIN", new_y="NEXT")
    pdf.ln(1)
    add_read_qc_legend(pdf)

    pitch = tagmaplib.WELL_PITCH_MM
    x0, y0 = pdf.get_x(), pdf.get_y()
    x = x0
    for direction, wells in direction_wells.items():
        n_rows, n_cols = tagmaplib.plate_layout(wells.keys())
        pdf.set_xy(x, y0)
        pdf.set_font("helvetica", style="B", size=10)
        pdf.cell(0, 6, direction.capitalize())
        draw_well_grid(pdf, x, y0 + 7, wells, n_rows, n_cols)
        x += WELL_LABEL + n_cols * pitch + 10
    pdf.set_xy(x0, y0)


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

    pdf = FPDF(orientation="P", unit="mm", format="A4")
    pdf.set_auto_page_break(True, margin=10)

    pdf.add_page()
    pdf.set_font("helvetica", style="B", size=20)
    pdf.cell(0, 14, "TagMap summary", new_x="LMARGIN", new_y="NEXT")

    first_section = True
    for attr, title, intro in SECTIONS:
        path = getattr(args, attr)
        if path is None:
            continue
        df = tagmaplib.tidy_numeric_dtypes(pd.read_csv(path, sep="\t"))

        if not first_section:
            pdf.add_page()
        first_section = False

        pdf.set_font("helvetica", style="B", size=14)
        pdf.cell(0, 9, title, new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("helvetica", size=9)
        pdf.multi_cell(0, 4.5, intro, new_x="LMARGIN", new_y="NEXT")
        pdf.ln(2)

        if df.shape[0] == 0:
            pdf.set_font("helvetica", style="I", size=9)
            pdf.cell(0, 6, "No data.")
            continue

        row_colors = None
        bold_rows = None
        if attr == "sanger_clones":
            add_legend(
                pdf,
                has_region="in_region" in df.columns,
                has_duplicates=bool(dedup_names),
            )
            row_colors = [
                tagmaplib.clone_color(
                    row.forward_status,
                    row.reverse_status,
                    row.both_sides_confirmed,
                    getattr(row, "ngs_verification", None),
                )
                for row in df.itertuples()
            ]
            if "in_region" in df.columns:
                bold_rows = {
                    i
                    for i, (color, row) in enumerate(zip(row_colors, df.itertuples()))
                    if color == tagmaplib.CLONE_GREEN and row.in_region
                }

        add_table(pdf, df, row_colors, bold_rows)

        if attr == "sanger_qc" and args.sanger_read_qc:
            read_qc = pd.read_csv(args.sanger_read_qc, sep="\t")
            if read_qc.shape[0]:
                for sample_name, group in read_qc.groupby("sample_name"):
                    draw_read_qc_plate(pdf, sample_name, group)

        if attr == "sanger_clones":
            n_unmobilized = int((df["summary"] == "unmobilized").sum())
            pdf.ln(2)
            pdf.set_font("helvetica", style="B", size=9)
            pdf.multi_cell(
                0,
                5,
                f"{n_unmobilized} of {df.shape[0]} clones unmobilized (both "
                "primers independently confirm the original insertion "
                "site) - see the Sanger position counts table's own "
                "\"unmobilized\" row for clones confirmed there from only "
                "one primer",
                new_x="LMARGIN",
                new_y="NEXT",
            )
            # Two plates to a page: force a fresh page for every other
            # plate actually drawn, and let draw_plate pack the one in
            # between below it if it fits.
            plates_drawn = 0
            for sample_name, group in df.groupby("sample_name"):
                if draw_plate(
                    pdf, sample_name, group, dedup_names, new_page=plates_drawn % 2 == 0
                ):
                    plates_drawn += 1

    pdf.output(args.output)
