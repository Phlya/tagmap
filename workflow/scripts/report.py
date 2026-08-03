"""Assemble the pipeline's summary tables into one markdown report."""

import argparse

import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--ngs-qc", default=None)
argparser.add_argument("--sanger-qc", default=None)
argparser.add_argument("--sanger-clones", default=None)
argparser.add_argument("--validation", default=None)
argparser.add_argument("--output", "-o", required=True)
args = argparser.parse_args()

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
        "forward primer, the reverse primer, both, or neither - and whether an "
        "unconfirmed side failed QC or was simply never sequenced. Where NGS "
        "data for the same material is available, also whether that clone's "
        "site was independently confirmed there, and from which side(s).",
    ),
    (
        "validation",
        "Sanger vs NGS validation",
        "Sanger sites cross-checked against the NGS data from the same "
        "material, run by run and in total.",
    ),
]

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
