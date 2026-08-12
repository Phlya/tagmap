import glob
import json
import os
from os.path import normpath

import numpy as np
import pandas as pd
from snakemake.utils import validate

if "use_only_read_junctions" in config:
    raise ValueError(
        "use_only_read_junctions has been folded into peak_caller, so that one "
        "setting decides which pairs are used and how, instead of two that had "
        "to agree. Replace it with peak_caller: 'coverage' (was False), "
        "'coverage_junctions' (was True), or 'junction_tiered' - the default, "
        "which needs no such choice because it tiers the pairs itself. See "
        "workflow/schemas/config.schema.yaml for what each one does."
    )

validate(config, schema="../schemas/config.schema.yaml")

# Only the pairs that read through the cassette/genome junction are of any use
# to the coverage caller's junction mode; the tiered caller sorts them out for
# itself and so needs to be handed both kinds.
only_read_junctions = config["peak_caller"] == "coverage_junctions"

# Every output path below can be set individually, but defaults to living
# under results_folder, so that in the common case one path is enough for a
# whole self-contained, namespaced results tree.
config.setdefault("results_folder", "results")
for _key, _subpath in {
    "fastq_folder": "fastq",
    "bams_folder": "bams",
    "pairs_folder": "pairs",
    "coverage_folder": "coverage",
    "peaks_folder": "peaks",
    "insertion_sites_folder": "insertion_sites",
    "sanger_folder": "sanger",
    "validation_folder": "validation",
    "stats_folder": "stats",
    "primer_position_file": "primer_positions.json",
    "sanger_primer_position_file": "sanger_primer_positions.json",
    "original_site_file": "original_site.json",
    "fasta_index_file": "refgen.fai",
}.items():
    config.setdefault(_key, os.path.join(config["results_folder"], _subpath))

# The Sanger branch defaults to the same primers as the NGS branch, but can be
# pointed at a different pair - see sanger_forward_primer_sequence in the
# schema for why the two sometimes differ. Neither forward_primer_sequence
# nor reverse_primer_sequence is required at this point - a Sanger-only
# project may skip them entirely and give sanger_forward_primer_sequence/
# sanger_reverse_primer_sequence instead (checked once sample_list/
# sanger_sample_list are known, further down).
config.setdefault("sanger_forward_primer_sequence", config.get("forward_primer_sequence"))
config.setdefault("sanger_reverse_primer_sequence", config.get("reverse_primer_sequence"))


# Anchored to the workflow rather than the working directory, so that the
# pipeline can be run from anywhere - including from .test, as CI does.
scripts_dir = os.path.join(workflow.basedir, "scripts")

import sys

sys.path.insert(0, scripts_dir)
import tagmaplib

refgen_path = normpath(config["refgen_path"])

fastq_folder = normpath(config["fastq_folder"])
bams_folder = normpath(config["bams_folder"])
pairs_folder = normpath(config["pairs_folder"])
coverage_folder = normpath(config["coverage_folder"])
peaks_folder = normpath(config["peaks_folder"])
insertion_sites_folder = normpath(config["insertion_sites_folder"])
sanger_folder = normpath(config["sanger_folder"])
validation_folder = normpath(config["validation_folder"])
stats_folder = normpath(config["stats_folder"])
original_site_file = normpath(config["original_site_file"])

# Two ways of pointing at the same pre-mobilization locus - see
# original_insertion_site/original_insertion_upstream_seq in the schema.
if config.get("original_insertion_site") and config.get(
    "original_insertion_upstream_seq"
):
    raise ValueError(
        "Only one of original_insertion_site and original_insertion_upstream_seq "
        "should be given - they are two ways of specifying the same locus."
    )
has_original_site = bool(
    config.get("original_insertion_site")
    or config.get("original_insertion_upstream_seq")
)


def cassette_length():
    """Length of the cassette contig, in chrom_sizes_path.

    Deferred rather than read as soon as the config is parsed, because
    chrom_sizes_path need not exist yet at that point - the make_chromsizes
    rule can still be the one to build it, as an ordinary dependency of
    whichever rule first calls this, rather than something the user has to
    generate by hand before the DAG can even be built. Sanger-only runs never
    call this at all, so for them chrom_sizes_path need not exist on disk,
    only be set. Deliberately not cached: chrom_sizes_path is itself often a
    make_chromsizes output, so a value read from it during this same run's
    DAG-building pass (before that rule has actually executed) can be stale -
    every caller needs its own fresh read, not one frozen for the whole run.
    """
    chromsizes = pd.read_table(
        config["chrom_sizes_path"],
        header=None,
        sep="\t",
        index_col=0,
        names=["chrom", "size"],
    )
    if config["cassette_name"] not in chromsizes.index:
        raise ValueError(
            f"Cassette {config['cassette_name']!r} is not in "
            f"{config['chrom_sizes_path']}. cassette_name must name a contig of "
            "the reference genome."
        )
    return chromsizes.loc[config["cassette_name"]]["size"]


# Depending on the mapper, the index files will be different
if config["mapper"] == "bwa-mem":
    idx = multiext(refgen_path, ".amb", ".ann", ".bwt", ".pac", ".sa")
elif config["mapper"] == "bwa-mem2":
    idx = multiext(refgen_path, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
elif config["mapper"] == "bwa-meme":
    idx = multiext(
        refgen_path,
        ".0123",
        ".amb",
        ".ann",
        ".pac",
        ".pos_packed",
        ".suffixarray_uint64",
        ".suffixarray_uint64_L0_PARAMETERS",
        ".suffixarray_uint64_L1_PARAMETERS",
        ".suffixarray_uint64_L2_PARAMETERS",
    )
elif config["mapper"] == "minibwa":
    idx = multiext(refgen_path, ".mbw", ".l2b")


def read_sample_sheet(path, schema):
    """Read a sample sheet, or return an empty frame if it is not configured."""
    if not path:
        return pd.DataFrame()
    samples = pd.read_table(path, comment="#", dtype=str)
    samples = samples.dropna(how="all")
    if samples.shape[0] == 0:
        return pd.DataFrame()
    validate(samples, schema=schema)
    return samples


samples = read_sample_sheet(
    config.get("samples_path"), "../schemas/samples.schema.yaml"
)
sanger_samples = read_sample_sheet(
    config.get("sanger_samples_path"), "../schemas/sanger_samples.schema.yaml"
)

sample_list = sorted(samples["name"].unique()) if samples.shape[0] else []
sanger_sample_list = (
    sorted(sanger_samples["name"].unique()) if sanger_samples.shape[0] else []
)

if not sample_list and not sanger_sample_list:
    raise ValueError(
        "No samples to process. Point samples_path at a TSV of NGS libraries, "
        "sanger_samples_path at a TSV of Sanger runs, or both."
    )

if sample_list and not config.get("chrom_sizes_path_no_cassette"):
    raise ValueError(
        "chrom_sizes_path_no_cassette is needed to write genome browser tracks "
        "for NGS samples. It should be the chrom sizes of the genome alone, "
        "without the cassette contig."
    )

if sample_list and not (
    config.get("forward_primer_sequence") and config.get("reverse_primer_sequence")
):
    raise ValueError(
        "forward_primer_sequence and reverse_primer_sequence are required "
        "when samples_path is set - the NGS branch uses them to tell real "
        "tagmentation pairs from the rest."
    )

if sanger_sample_list and not (
    config.get("sanger_forward_primer_sequence")
    and config.get("sanger_reverse_primer_sequence")
):
    raise ValueError(
        "forward_primer_sequence and reverse_primer_sequence (or, if the "
        "Sanger primer differs from the NGS one, sanger_forward_primer_sequence "
        "and sanger_reverse_primer_sequence) are required when "
        "sanger_samples_path is set."
    )

if sanger_samples.shape[0]:
    has_input = pd.Series(False, index=sanger_samples.index)
    for column in ("ab1", "fastq"):
        if column in sanger_samples.columns:
            has_input |= sanger_samples[column].notna()
    if not has_input.all():
        missing = sorted(sanger_samples.loc[~has_input, "name"].unique())
        raise ValueError(
            f"Sanger samples {missing} have neither an ab1 nor a fastq value. "
            "One of the two is needed to know what to read."
        )

# Cross-validation needs both kinds of data.
do_validation = bool(
    config["validate_sanger_with_ngs"] and sample_list and sanger_sample_list
)


def sanger_rows(sample):
    return sanger_samples[sanger_samples["name"] == sample]


def sanger_row_value(sample, column):
    """A per-sample column of the Sanger sheet, or None if unset/inconsistent."""
    if column not in sanger_samples.columns:
        return None
    values = sanger_rows(sample)[column].dropna().unique()
    return values[0] if len(values) == 1 else None


def ab1_files(sample):
    """Expand the ab1 column of a Sanger sample into actual trace files.

    The traces exist before the run, so a plain glob at DAG build time is
    enough and the files still get tracked as inputs.
    """
    rows = sanger_rows(sample)
    if "ab1" not in rows.columns:
        return []
    files = []
    for pattern in rows["ab1"].dropna():
        pattern = str(pattern)
        if os.path.isdir(pattern):
            pattern = os.path.join(pattern, "*.ab1")
        files.extend(sorted(glob.glob(pattern)))
    return files


def sanger_input_files(wildcards):
    """Whatever a Sanger sample is built from: its traces, or ready fastqs."""
    files = ab1_files(wildcards.sample)
    if files:
        return files
    rows = sanger_rows(wildcards.sample)
    return sorted(str(f) for f in rows.get("fastq", pd.Series(dtype=str)).dropna())


def get_filter(side, primer_positions_file):
    """pairtools expression selecting pairs anchored at one ITR primer."""
    with open(primer_positions_file) as f:
        primer_positions = json.load(f)
    walk_pair_type = (
        '(walk_pair_type in ["R1", "R2", "R1&2"]) and' if only_read_junctions else ""
    )
    forward_ITR_primer_position = primer_positions["forward_ITR_primer_position"]
    reverse_ITR_primer_position = primer_positions["reverse_ITR_primer_position"]
    selection_3prime_forward_readthrough = f"(abs(pos2-{cassette_length()})<=2)"
    selection_3prime_reverse_readthrough = f"(pos2<=2)"
    # "All" has to be a superset of the readthrough case. Writing it as the
    # single position min(cassette end, primer + read length) is only that when
    # the reads are long enough to reach the end of the cassette; with shorter
    # ones the min picks where a read that stopped inside the cassette ends,
    # and silently throws away every pair that did read through - the pairs the
    # junction_tiered caller relies on most.
    selection_3prime_forward_all = (
        f"((abs(pos2-min({cassette_length()}, {forward_ITR_primer_position}+read_len2))<=2)"
        f" or {selection_3prime_forward_readthrough})"
    )
    selection_3prime_reverse_all = (
        f"((abs(pos2-max(0, {reverse_ITR_primer_position}-read_len2))<=2)"
        f" or {selection_3prime_reverse_readthrough})"
    )
    selection_3prime_forward = (
        selection_3prime_forward_readthrough
        if only_read_junctions
        else selection_3prime_forward_all
    )
    selection_3prime_reverse = (
        selection_3prime_reverse_readthrough
        if only_read_junctions
        else selection_3prime_reverse_all
    )
    selection_3prime = (
        selection_3prime_forward if side == "forward" else selection_3prime_reverse
    )
    return (
        f"""((pair_type in ["UR", "UU", "RU"]) and (chrom1!=chrom2) and """
        f"""{walk_pair_type} (chrom2=="{config['cassette_name']}") and """
        f"""{selection_3prime})"""
    )


def workflow_targets():
    """Final outputs, depending on which kinds of samples were configured."""
    targets = []
    if sample_list:
        if config["do_fastqc"]:
            targets += expand(
                f"{fastq_folder}/{{sample}}.{{read}}_fastqc.html",
                sample=sample_list,
                read=["R1", "R2"],
            )
        targets += expand(f"{pairs_folder}/{{sample}}_stats.yml", sample=sample_list)
        targets += expand(
            f"{coverage_folder}/{{sample}}_{{side}}_coverage_for_ucsc.bedgraph",
            sample=sample_list,
            side=["forward", "reverse"],
        )
        targets += [
            f"{peaks_folder}/all_peaks.bed",
            f"{insertion_sites_folder}/all_sites.bed",
            f"{insertion_sites_folder}/confirmed_ngs_sites.bed",
            f"{insertion_sites_folder}/confirmed_ngs_sites_no_cassette.bed",
            f"{insertion_sites_folder}/sample_summary.tsv",
            f"{stats_folder}/ngs_qc_stats.tsv",
        ]
    if sanger_sample_list:
        targets += expand(
            f"{sanger_folder}/{{sample}}_reads.tsv", sample=sanger_sample_list
        )
        targets += [
            f"{sanger_folder}/all_sanger_sites.bed",
            f"{sanger_folder}/confirmed_sanger_sites.bed",
            f"{sanger_folder}/confirmed_sanger_sites_no_cassette.bed",
            f"{sanger_folder}/confirmed_sanger_sites_region.bed",
            f"{sanger_folder}/confirmed_sanger_sites_deduplicated.bed",
            f"{sanger_folder}/confirmed_sanger_sites_region_deduplicated.bed",
            f"{stats_folder}/sanger_qc_stats.tsv",
            f"{stats_folder}/sanger_clone_summary.tsv",
            f"{stats_folder}/sanger_positions.tsv",
            f"{stats_folder}/sanger_position_counts.tsv",
        ]
    if do_validation:
        targets += [
            f"{validation_folder}/sanger_vs_ngs.tsv",
            f"{validation_folder}/confirmed_sites.bed",
            f"{stats_folder}/validation_summary.tsv",
        ]
    targets += [f"{stats_folder}/report.md", f"{stats_folder}/report.pdf"]
    return targets


def sanger_ngs_pairs():
    """sanger_sample=ngs_sample assignments taken from the Sanger sheet.

    The left-hand side has to match the "sample_name" values compare_sanger_ngs.py
    actually sees in all_sanger_sites.bed - sanger_sites.py renames those from
    the raw sheet name via tagmaplib.rename_clone_id, so the same rename is
    applied here rather than passing the pre-rename sheet name through.
    """
    pairs = []
    for sample in sanger_sample_list:
        ngs_sample = sanger_row_value(sample, "ngs_sample")
        if ngs_sample:
            if ngs_sample not in sample_list:
                raise ValueError(
                    f"Sanger sample {sample!r} names NGS sample {ngs_sample!r}, "
                    f"which is not in {config['samples_path']}."
                )
            renamed_sample, _ = tagmaplib.rename_clone_id(sample, "")
            pairs.append(f"{renamed_sample}={ngs_sample}")
    return pairs
