# Configuration

Copy `example_config.yaml` to `config.yaml` and edit it for your experiment.
Every key is documented inline; only `refgen_path`, `chrom_sizes_path`,
`cassette_name`, `forward_primer_sequence` and `reverse_primer_sequence` are
required, everything else has a default (see
`workflow/schemas/config.schema.yaml` for the authoritative list, which
`snakemake` validates the config against on every run).

## Sample sheets

Give either or both, depending on what you're analysing:

- **`samples.tsv`** (`samples_path`) - paired-end NGS tagmentation libraries.
  Columns: `name`, `fastq1`, `fastq2`. Rows sharing a `name` are merged, so a
  library split across lanes can be listed as several rows.
- **`sanger_samples.tsv`** (`sanger_samples_path`) - Sanger runs of individual
  long reads, for pinpointing the integration site in a clonal line. Columns:
  `name`, and either `ab1` (a folder or glob of `.ab1` traces) or `fastq`
  (already-basecalled reads); optionally `direction` (`forward`/`reverse`),
  `clone`, and `ngs_sample` (an NGS sample from `samples.tsv` covering the same
  material, used to cross-validate the Sanger site - see
  `validate_sanger_with_ngs`).

Read-level metadata such as well, clone and primer direction is normally
pulled from the trace file names via `sanger_name_regex`; sample-sheet columns
fill in whatever the regex doesn't capture, or override it for a whole row.

Both sheets take `#`-prefixed comment lines and blank rows.

## Output paths

Every stage's output folder (`fastq_folder`, `bams_folder`, `pairs_folder`, ...)
and the `primer_position_file`/`fasta_index_file` can be set individually, but
default to living under `results_folder` (itself defaulting to `results`), so
usually just setting `results_folder` is enough to get a self-contained,
namespaced results tree - only override an individual path if you want that
one file somewhere else.

## Summary tables

Besides the site calls, the workflow writes small TSV summaries to
`stats_folder`, and assembles them into a single `report.md`:

- `ngs_qc_stats.tsv` - per NGS library: read pairs anchored at an ITR primer
  (evidence of mobilization), and how many of the resulting insertion sites
  were seen from both sides of the cassette versus only one.
- `sanger_qc_stats.tsv` - Sanger reads per run and primer direction, how many
  passed QC, and a summary of why the rest didn't.
- `sanger_clone_summary.tsv` - per clone, whether the forward and reverse
  primers each gave a confirmed site, failed QC, or were never sequenced, and
  - when NGS data for the same material is available - whether that clone's
  site was independently confirmed there, and from which side(s). When
  `original_insertion_site`/`original_insertion_upstream_seq` is configured, a
  side whose site lands there instead of a new locus is classified
  `unmobilized` rather than `confirmed`, and the clone gets an `unmobilized`
  column.
- `validation_summary.tsv` - Sanger sites confirmed by the NGS data, per run
  and in total (only written when `validate_sanger_with_ngs` applies).
