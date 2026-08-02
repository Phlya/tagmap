"""Assemble the pipeline's summary tables into one markdown report."""

import argparse

import pandas as pd

import tagmaplib

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--ngs-mapping", default=None)
argparser.add_argument("--ngs-sidedness", default=None)
argparser.add_argument("--sanger-reads", default=None)
argparser.add_argument("--sanger-fail-reasons", default=None)
argparser.add_argument("--sanger-clones", default=None)
argparser.add_argument("--validation", default=None)
argparser.add_argument("--output", "-o", required=True)
args = argparser.parse_args()

SECTIONS = [
    (
        "ngs_mapping",
        "NGS mobilization",
        "Read pairs per library, and how many were anchored at an ITR primer - "
        "the evidence that a pair actually captured a transposition event.",
    ),
    (
        "ngs_sidedness",
        "NGS insertion site sidedness",
        "Insertion sites found from both sides of the cassette versus only one. "
        "A one-sided site is weaker evidence, since it has not been confirmed by "
        "an independent primer.",
    ),
    (
        "sanger_reads",
        "Sanger reads",
        "Sanger reads per run and ITR primer direction, and how many passed QC.",
    ),
    (
        "sanger_fail_reasons",
        "Sanger QC failure reasons",
        "Why the reads above that did not pass QC failed.",
    ),
    (
        "sanger_clones",
        "Sanger clones",
        "For each clone, whether the integration site was confirmed from the "
        "forward primer, the reverse primer, both, or neither - and whether an "
        "unconfirmed side failed QC or was simply never sequenced.",
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
