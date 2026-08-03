"""Assemble the pipeline's summary tables into one markdown report.

SECTIONS is also imported by report_pdf.py, so the two reports share the
same section titles and descriptions.
"""

import argparse

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
        "and why the rest didn't.",
    ),
    (
        "sanger_clones",
        "Sanger clones",
        "For each clone, whether the integration site was confirmed from the "
        "forward primer, the reverse primer, both, or neither - and why an "
        "unconfirmed side isn't: never sequenced, or how it failed QC. "
        "position lists the site(s) clustered from that clone's passing "
        "reads (all_sanger_sites.bed) - usually one, but a clone whose "
        "forward and reverse reads don't cluster together shows both, so "
        "you can see where they disagree. Passing QC on both sides isn't "
        "enough on its own: positions_agree checks that the forward and "
        "reverse reads actually cluster onto the same site rather than "
        "pointing at different loci, and both_sides_confirmed requires "
        "both. Where NGS data for the same material is available, also "
        "whether that clone's site was independently confirmed there, and "
        "from which side(s).",
    ),
    (
        "validation",
        "Sanger vs NGS validation",
        "Sanger sites cross-checked against the NGS data from the same "
        "material, run by run and in total.",
    ),
]

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument("--ngs-qc", default=None)
    argparser.add_argument("--sanger-qc", default=None)
    argparser.add_argument("--sanger-clones", default=None)
    argparser.add_argument("--validation", default=None)
    argparser.add_argument("--output", "-o", required=True)
    args = argparser.parse_args()

    lines = ["# TagMap summary", ""]

    for attr, title, intro in SECTIONS:
        path = getattr(args, attr)
        if path is None:
            continue
        df = pd.read_csv(path, sep="\t")
        lines.append(f"## {title}")
        lines.append("")
        lines.append(intro)
        lines.append("")
        if df.shape[0] == 0:
            lines.append("_No data._")
        else:
            lines.append(tagmaplib.to_markdown(df))
        lines.append("")

    with open(args.output, "w") as f:
        f.write("\n".join(lines))
