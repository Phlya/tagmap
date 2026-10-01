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
  An optional `plates` column (comma-separated) restricts which plates to
  expect in a library when demultiplexing by barcode (see below).
- **`barcodes.tsv`** (`barcodes_path`) - only for barcoded libraries. Columns:
  `plate`, `side` (`forward` or `reverse`) and `barcode`. See
  [Barcoded libraries](#barcoded-libraries).
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

## Barcoded libraries

If the ITR primers carry a barcode in front of the ITR sequence - to label
plates, or just to add sequence diversity at the start of the library - set
`barcodes_path` and every library in `samples.tsv` is split by it first, into
one sample per plate named `{plate}_{name}`. A typical layout is one FASTQ pair
per well (from the Illumina indexes) and the in-read barcode saying which plate
the clone came from.

```
plate   side     barcode
plate1  forward  TTATTCGAGG
plate1  forward  AAGATTGGATA
plate1  reverse  GTAACATGCG
plate1  reverse  ACAGTAACTAT
plate2  forward  TGTAGTCATTGG
...
```

The barcode starts the read that starts at the ITR primer (R2 in the usual
TagMap layout; `barcode_reads`). Forward-side barcodes sit on forward-ITR reads
and reverse-side ones on reverse-ITR reads. A plate can list any number of
barcodes per side, of different lengths if you like; a read carries any one of
them and they all label that plate. The matched barcode is trimmed, so the rest
of the workflow sees a read that starts at the ITR primer, exactly as without
barcodes. `barcode_max_mismatch` mismatches are tolerated, and the barcodes are
checked when the workflow starts to be distinguishable at that tolerance (a
barcode that is a prefix of another, or only a base or two away from it, is an
error). With `barcode_require_primer`, a barcode must also be followed by the
ITR primer of its own side, which also checks the read layout.

Reads with no barcode, with conflicting barcodes, or (with
`barcode_require_primer`) with the wrong primer are not assigned to any plate;
they are counted in `stats/demux_stats.tsv` and the report.

## Clonal libraries

`ngs_mode: clone` is for libraries of single clonal lines, each expected to
carry one insertion at a few hundred to a few thousand reads. Sites are called
exactly as for a pool; each library is then judged in
`stats/ngs_clone_summary.tsv` on whether one position dominates and nothing
else is present above background. Abundance is counted in independent
molecules (distinct tagmentation positions), not reads, so PCR copies do not
count. Statuses: `clean`, `contaminated` (a second site at or above
`clone_max_contamination_frac`; `shared_with` names the clones whose dominant
site it is, which points at cross-talk rather than a second insertion),
`multiple` (a second site at or above `clone_multi_site_frac`), `unmobilized`
(clean, but at the original insertion site), `weak`, `no_insertion` and
`too_few_reads`. The thresholds depend on how much cross-talk a run has, so
tune them on real data. The report adds a colour-coded plate grid per plate.

A library with no reads at all (e.g. a plate absent from a well) is fine - it
comes out as `too_few_reads`.

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
- `demux_stats.tsv` - with `barcodes_path`: read pairs per library, plate and
  barcode, and how many were left unassigned and why.
- `ngs_clone_summary.tsv` - with `ngs_mode: clone`: per clonal library, the
  dominant insertion site, its share of all molecules, the runner-up and a
  clean/contaminated/multiple/... status.
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
- `sanger_positions.tsv` - distinct genomic positions found across the
  clustered Sanger sites, per plate and in total (a locus seen on more than
  one plate, e.g. an unmobilized founder control resequenced, still counts
  once).
- `validation_summary.tsv` - Sanger sites confirmed by the NGS data, per run
  and in total (only written when `validate_sanger_with_ngs` applies).
