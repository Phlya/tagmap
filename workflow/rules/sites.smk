localrules:
    ngs_stats,


rule for_ucsc:
    input:
        bg=f"{coverage_folder}/{{sample}}_{{side}}_coverage.bedgraph",
        script=f"{scripts_dir}/for_ucsc.py",
    output:
        bg=f"{coverage_folder}/{{sample}}_{{side}}_coverage_for_ucsc.bedgraph",
        bw=f"{coverage_folder}/{{sample}}_{{side}}_coverage_for_ucsc.bw",
    log:
        "logs/for_ucsc/{sample}_{side}.log",
    conda:
        "../envs/all.yaml"
    params:
        chromsizes=config.get("chrom_sizes_path_no_cassette", ""),
        name=lambda wildcards: f"{wildcards.sample}_{wildcards.side}",
    shell:
        """
        python3 {input.script} -i {input.bg} --chromsizes {params.chromsizes} \
            --name {params.name} -o {output.bg} --output-bigwig {output.bw} \
            >{log[0]} 2>&1
        """


# Only the junction_tiered caller records a per-peak tier, so only it has
# evidence files - for sample_summary.py's own junction/anchor split, and for
# find_insertion_sites.py's --min-orientation-support-confirmed.
_junction_evidence = (
    expand(
        f"{peaks_folder}/{{sample}}_{{side}}_evidence.tsv",
        sample=sample_list,
        side=["forward", "reverse"],
    )
    if config["peak_caller"] == "junction_tiered"
    else []
)


rule find_insertion_sites:
    input:
        peaks=f"{peaks_folder}/all_peaks.bed",
        evidence=_junction_evidence,
        genome=refgen_path,
        genome_index=config["fasta_index_file"],
        chromsizes=config.get("chrom_sizes_path_no_cassette", []),
        script=f"{scripts_dir}/find_insertion_sites.py",
    output:
        output=f"{insertion_sites_folder}/all_sites.bed",
        for_ucsc=f"{insertion_sites_folder}/all_sites_for_ucsc.bed",
        confirmed=f"{insertion_sites_folder}/confirmed_ngs_sites.bed",
        confirmed_no_cassette=f"{insertion_sites_folder}/confirmed_ngs_sites_no_cassette.bed",
    log:
        "logs/find_insertion_sites/all.log",
    benchmark:
        "benchmarks/find_insertion_sites/all.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        insertion_seq=config["insertion_seq"],
        max_dist_between_sides=config["max_dist_between_sides"],
        snap_window=config["snap_window"],
        min_orientation_support=config["min_orientation_support"],
        min_orientation_support_confirmed=config["min_orientation_support_confirmed"],
        chromsizes_arg=lambda wildcards, input: (
            f"--chromsizes {input.chromsizes}" if input.chromsizes else ""
        ),
        evidence_arg=lambda wildcards, input: (
            f"--evidence {' '.join(input.evidence)}" if input.evidence else ""
        ),
    shell:
        """
        python3 {input.script} --peaks {input.peaks} {params.evidence_arg} \
            --max-dist {params.max_dist_between_sides} \
            --genome {input.genome} --genome-index {input.genome_index} \
            --insertion-seq {params.insertion_seq} \
            --snap-window {params.snap_window} \
            --min-orientation-support {params.min_orientation_support} \
            --min-orientation-support-confirmed {params.min_orientation_support_confirmed} \
            {params.chromsizes_arg} \
            -o {output.output} --output-for-ucsc {output.for_ucsc} \
            --output-confirmed {output.confirmed} \
            --output-confirmed-no-cassette {output.confirmed_no_cassette} \
            >{log[0]} 2>&1
        """


rule sample_summary:
    input:
        sites=f"{insertion_sites_folder}/all_sites.bed",
        peaks=f"{peaks_folder}/all_peaks.bed",
        evidence=_junction_evidence,
        coverage=expand(
            f"{coverage_folder}/{{sample}}_{{side}}_coverage.bedgraph",
            sample=sample_list,
            side=["forward", "reverse"],
        ),
        # Not read directly - construct_contigs() reads these paths itself
        # (see common.smk) - declared as inputs only so Snakemake builds
        # chrom_sizes_path (a make_chromsizes output) first, and reruns this
        # rule if either chromsizes file's contents change.
        chrom_sizes_path=config["chrom_sizes_path"],
        chrom_sizes_path_no_cassette=config.get("chrom_sizes_path_no_cassette", []),
        script=f"{scripts_dir}/sample_summary.py",
    output:
        f"{insertion_sites_folder}/sample_summary.tsv",
    log:
        "logs/sample_summary/log.log",
    benchmark:
        "benchmarks/sample_summary/benchmark.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        construct_contigs=lambda wildcards, input: " ".join(construct_contigs()),
        max_dist_between_sides=config["max_dist_between_sides"],
        evidence_arg=lambda wildcards, input: (
            f"--evidence {' '.join(input.evidence)}" if input.evidence else ""
        ),
    shell:
        """
        python3 {input.script} --sites {input.sites} --peaks {input.peaks} \
            {params.evidence_arg} --coverage {input.coverage} \
            --construct-contigs {params.construct_contigs} \
            --max-dist {params.max_dist_between_sides} \
            -o {output} >{log[0]} 2>&1
        """


rule ngs_stats:
    input:
        stats=expand(f"{pairs_folder}/{{sample}}_stats.yml", sample=sample_list),
        raw_pairs=expand(f"{pairs_folder}/{{sample}}_sorted.pairs", sample=sample_list),
        dedup_pairs=expand(
            (
                f"{pairs_folder}/{{sample}}_dupmarked.pairs"
                if config["dedup"]
                else f"{pairs_folder}/{{sample}}_sorted.pairs"
            ),
            sample=sample_list,
        ),
        pairs=expand(
            f"{pairs_folder}/{{sample}}_{{side}}.pairs",
            sample=sample_list,
            side=["forward", "reverse"],
        ),
        peaks=f"{peaks_folder}/all_peaks.bed",
        sites=f"{insertion_sites_folder}/all_sites.bed",
        primer_positions=config["primer_position_file"],
        script=f"{scripts_dir}/ngs_stats.py",
    output:
        f"{stats_folder}/ngs_qc_stats.tsv",
    log:
        "logs/ngs_stats/log.log",
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        only_junctions_arg="--only-junctions" if only_read_junctions else "",
        cassette_name=config["cassette_name"],
        cassette_length=cassette_length(),
    shell:
        """
        python3 {input.script} --stats-yml {input.stats} --raw-pairs {input.raw_pairs} \
            --dedup-pairs {input.dedup_pairs} \
            --pairs {input.pairs} --peaks {input.peaks} --sites {input.sites} \
            --primer-positions {input.primer_positions} \
            --cassette-name {params.cassette_name} --cassette-length {params.cassette_length} \
            {params.only_junctions_arg} \
            -o {output} \
            >{log[0]} 2>&1
        """
