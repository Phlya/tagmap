# Snakemake workflow: `TagMap`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for mapping transposon integration sites from tagmentation
mapping data (https://www.biorxiv.org/content/10.1101/037762v5.full), either
from pooled NGS libraries or from individual Sanger reads of clonal lines.

Two complementary ways of finding an integration site are supported, and can
be run together or on their own:

- **NGS**: paired-end tagmentation libraries are mapped, filtered down to
  read pairs anchored at an ITR primer, and turned into peaks and then
  single-base insertion sites. Suited to pooled populations, and reports both
  sides of an insertion when both were captured.
- **Sanger**: individual long reads from a clonal line - typically `.ab1`
  traces from one well per clone - are mapped directly, and the exact
  cassette/genome junction is read off each read. Suited to identifying the
  integration site in a specific clone.

When both are configured for the same material, Sanger sites can be
cross-validated against the NGS data (`validate_sanger_with_ngs`): a Sanger
read that only worked from one ITR primer can still be trusted if the NGS
library shows an insertion at the same position from both sides.

Alongside the site calls, the workflow writes a set of summary tables to
`stats_folder` - NGS mobilization rate and insertion-site sidedness, Sanger
read/QC counts and per-clone forward/reverse coverage, and the
Sanger-vs-NGS validation outcome - and assembles them into one `report.md`.
A founder/control clone that still carries the cassette at its original,
pre-mobilization locus can be flagged as such (`original_insertion_site`),
rather than counted as a new integration.

## Usage

See `config/README.md` for how to configure the workflow, and
`config/example_config.yaml` for a fully commented example. A minimal test
dataset lives under `.test/`.

If you use this workflow in a paper, don't forget to give credits to the
authors by citing the URL of this (original) repository and its DOI (see
above).
