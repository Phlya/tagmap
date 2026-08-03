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


rule find_insertion_sites:
    input:
        peaks=f"{peaks_folder}/all_peaks.bed",
        genome=refgen_path,
        genome_index=config["fasta_index_file"],
        script=f"{scripts_dir}/find_insertion_sites.py",
    output:
        output=f"{insertion_sites_folder}/all_sites.bed",
        for_ucsc=f"{insertion_sites_folder}/all_sites_for_ucsc.bed",
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
    shell:
        """
        python3 {input.script} --peaks {input.peaks} \
            --max-dist {params.max_dist_between_sides} \
            --genome {input.genome} --genome-index {input.genome_index} \
            --insertion-seq {params.insertion_seq} \
            -o {output.output} --output-for-ucsc {output.for_ucsc} \
            >{log[0]} 2>&1
        """


# Only the junction_tiered caller records a per-peak tier, so only it has
# evidence files for sample_summary.py to read a junction/anchor split from.
_summary_evidence = (
    expand(
        f"{peaks_folder}/{{sample}}_{{side}}_evidence.tsv",
        sample=sample_list,
        side=["forward", "reverse"],
    )
    if config["peak_caller"] == "junction_tiered"
    else []
)


rule sample_summary:
    input:
        sites=f"{insertion_sites_folder}/all_sites.bed",
        peaks=f"{peaks_folder}/all_peaks.bed",
        evidence=_summary_evidence,
        coverage=expand(
            f"{coverage_folder}/{{sample}}_{{side}}_coverage.bedgraph",
            sample=sample_list,
            side=["forward", "reverse"],
        ),
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
        construct_contigs=" ".join(config["construct_contigs"]),
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
        sites=f"{insertion_sites_folder}/all_sites.bed",
        script=f"{scripts_dir}/ngs_stats.py",
    output:
        f"{stats_folder}/ngs_qc_stats.tsv",
    log:
        "logs/ngs_stats/log.log",
    conda:
        "../envs/all.yaml"
    threads: 1
    shell:
        """
        python3 {input.script} --stats-yml {input.stats} --sites {input.sites} \
            -o {output} \
            >{log[0]} 2>&1
        """
