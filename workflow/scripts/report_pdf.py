"""Render the pipeline's summary tables as a single PDF.

Shares its section titles and descriptions with report.py's markdown report.
The Sanger clones table is additionally colour-coded per clone - see
clone_color() for the green/yellow/red rules.
"""

import argparse

import pandas as pd
from fpdf import FPDF
from fpdf.fonts import FontFace

from report import SECTIONS

GREEN = (198, 239, 206)
YELLOW = (255, 235, 156)
RED = (255, 199, 206)
WHITE = (255, 255, 255)
HEADING_GREY = (221, 221, 221)

MAX_COL_WEIGHT = 45
MIN_COL_WEIGHT = 6


def clone_color(forward_status, reverse_status, both_sides_confirmed):
    """Which of GREEN/YELLOW/RED a clone's row should be filled with.

    - green: both sides confirmed and agree on position; or only one side
      was even sequenced, and it passed - a single clean read with nothing
      to contradict it.
    - yellow: one side passed, but the other was sequenced and its read
      just didn't come out clean (failed QC) - a real attempt that didn't
      pan out, short of positive evidence of a problem.
    - red: neither side passed - not sequenced, failed QC, or both - or
      both sides passed individually but disagree on position, which is
      positive evidence of a conflict.
    """
    if both_sides_confirmed:
        return GREEN
    forward_pass = forward_status == "PASSED"
    reverse_pass = reverse_status == "PASSED"
    if forward_pass and reverse_pass:
        return RED  # both passed, but both_sides_confirmed is False: positions disagree
    if forward_pass or reverse_pass:
        other_status = reverse_status if forward_pass else forward_status
        return GREEN if other_status == "not sequenced" else YELLOW
    return RED


def column_weights(df):
    """Relative column widths for pdf.table(), sized off each column's own
    content so long columns (e.g. summary) get proportionally more room."""
    weights = []
    for col in df.columns:
        values = df[col].astype(str)
        typical_len = values.str.len().quantile(0.9) if values.shape[0] else 0
        weights.append(min(max(len(col), typical_len, MIN_COL_WEIGHT), MAX_COL_WEIGHT))
    return weights


def add_table(pdf, df, row_colors=None):
    df = df.fillna("")
    pdf.set_font("helvetica", size=6.5)
    with pdf.table(
        col_widths=column_weights(df),
        text_align="LEFT",
        line_height=3.2,
        padding=1,
        headings_style=FontFace(emphasis="BOLD", fill_color=HEADING_GREY),
    ) as table:
        header = table.row()
        for col in df.columns:
            header.cell(col)
        for i, (_, data) in enumerate(df.iterrows()):
            row = table.row()
            style = FontFace(fill_color=row_colors[i] if row_colors is not None else WHITE)
            for col in df.columns:
                row.cell(str(data[col]), style=style)


def add_legend(pdf):
    pdf.set_font("helvetica", size=8)
    for color, label in (
        (GREEN, "both sides confirmed and agree, or only one side was sequenced and it passed"),
        (YELLOW, "one side passed, but the other was sequenced and failed QC"),
        (RED, "neither side passed, or both passed but disagree on position"),
    ):
        pdf.set_fill_color(*color)
        pdf.cell(4, 4, "", border=1, fill=True)
        pdf.set_x(pdf.get_x() + 2)
        pdf.cell(0, 4, label, new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument("--ngs-qc", default=None)
    argparser.add_argument("--sanger-qc", default=None)
    argparser.add_argument("--sanger-clones", default=None)
    argparser.add_argument("--validation", default=None)
    argparser.add_argument("--output", "-o", required=True)
    args = argparser.parse_args()

    pdf = FPDF(orientation="L", unit="mm", format="A4")
    pdf.set_auto_page_break(True, margin=10)

    pdf.add_page()
    pdf.set_font("helvetica", style="B", size=20)
    pdf.cell(0, 14, "TagMap summary", new_x="LMARGIN", new_y="NEXT")

    first_section = True
    for attr, title, intro in SECTIONS:
        path = getattr(args, attr)
        if path is None:
            continue
        df = pd.read_csv(path, sep="\t")

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
        if attr == "sanger_clones":
            add_legend(pdf)
            row_colors = [
                clone_color(row.forward_status, row.reverse_status, row.both_sides_confirmed)
                for row in df.itertuples()
            ]

        add_table(pdf, df, row_colors)

    pdf.output(args.output)
