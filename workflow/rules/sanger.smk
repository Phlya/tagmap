localrules:
    combine_sanger_sites,
    filter_confirmed_sites,
    sanger_stats,


rule ab1_to_fastq:
    input:
        traces=sanger_input_files,
        script=f"{scripts_dir}/ab1_to_fastq.py",
    output:
        f"{sanger_folder}/{{sample}}.fastq.gz",
    log:
        "logs/ab1_to_fastq/{sample}.log",
    benchmark:
        "benchmarks/ab1_to_fastq/{sample}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        basecaller=config["sanger_basecaller"],
    shell:
        """
        python3 {input.script} --basecaller {params.basecaller} \
            -o {output} {input.traces} \
            >{log[0]} 2>&1
        """


if config["mapper"] == "minibwa":

    rule minibwa_sanger_map:
        input:
            reads=f"{sanger_folder}/{{sample}}.fastq.gz",
            reference=refgen_path,
            idx=idx,
        output:
            f"{sanger_folder}/{{sample}}.bam",
        log:
            "logs/sanger_map/{sample}.log",
        benchmark:
            "benchmarks/sanger_map/{sample}.tsv"
        conda:
            "../envs/all.yaml"
        threads: 4
        params:
            extra=config["sanger_map_args"],
        shell:
            # Unsorted: sanger_sites reads the file straight through, and
            # sorting would only add a dependency.
            """
            minibwa mem -t {threads} {params.extra} {input.reference} {input.reads} 2>{log[0]} \
                | samtools view -b -o {output} - 2>>{log[0]}
            """

else:

    rule sanger_map:
        input:
            reads=[f"{sanger_folder}/{{sample}}.fastq.gz"],
            reference=refgen_path,
            idx=idx,
        output:
            f"{sanger_folder}/{{sample}}.bam",
        log:
            "logs/sanger_map/{sample}.log",
        benchmark:
            "benchmarks/sanger_map/{sample}.tsv"
        threads: 4
        params:
            bwa=config["mapper"],
            # Unsorted: sanger_sites reads the file straight through, and sorting
            # would only add a dependency.
            sort="none",
            dedup="none",
            extra=config["sanger_map_args"],
        wrapper:
            "v3.3.3/bio/bwa-memx/mem"


rule sanger_sites:
    input:
        bam=f"{sanger_folder}/{{sample}}.bam",
        genome=refgen_path,
        genome_index=config["fasta_index_file"],
        primer_positions=config["sanger_primer_position_file"],
        script=f"{scripts_dir}/sanger_sites.py",
    output:
        reads=f"{sanger_folder}/{{sample}}_reads.tsv",
        sites=f"{sanger_folder}/{{sample}}_sites.bed",
    log:
        "logs/sanger_sites/{sample}.log",
    benchmark:
        "benchmarks/sanger_sites/{sample}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        cassette_name=config["cassette_name"],
        insertion_seq=config["insertion_seq"],
        min_mapq=config["sanger_min_mapq"],
        min_aligned=config["sanger_min_aligned"],
        max_unexplained=config["sanger_max_unexplained"],
        max_dist=config["sanger_max_dist"],
        name_regex_arg=lambda wildcards: (
            f"--name-regex '{config['sanger_name_regex']}'"
            if config["sanger_name_regex"]
            else ""
        ),
        direction_map=" ".join(
            f"{key}={value}" for key, value in config["sanger_direction_map"].items()
        ),
        direction=lambda wildcards: (
            sanger_row_value(wildcards.sample, "direction")
            or config["sanger_default_direction"]
        ),
        clone_arg=lambda wildcards: (
            f"--clone {sanger_row_value(wildcards.sample, 'clone')}"
            if sanger_row_value(wildcards.sample, "clone")
            else ""
        ),
    shell:
        """
        python3 {input.script} --bam {input.bam} --genome {input.genome} \
            --genome-index {input.genome_index} \
            --primer-positions {input.primer_positions} \
            --sample-name {wildcards.sample} \
            --construct-contigs {params.cassette_name} \
            --insertion-seq {params.insertion_seq} \
            --min-mapq {params.min_mapq} --min-aligned {params.min_aligned} \
            --max-unexplained {params.max_unexplained} \
            --max-dist {params.max_dist} \
            --direction {params.direction} \
            --direction-map {params.direction_map} \
            {params.name_regex_arg} {params.clone_arg} \
            --output-reads {output.reads} --output-sites {output.sites} \
            >{log[0]} 2>&1
        """


rule combine_sanger_sites:
    input:
        sites=expand(f"{sanger_folder}/{{sample}}_sites.bed", sample=sanger_sample_list),
        blacklist=config.get("blacklist", []),
        script=f"{scripts_dir}/combine_sanger_sites.py",
    output:
        sites=f"{sanger_folder}/all_sanger_sites.bed",
        for_ucsc=f"{sanger_folder}/all_sanger_sites_for_ucsc.bed",
    log:
        "logs/combine_sanger_sites/log.log",
    benchmark:
        "benchmarks/combine_sanger_sites/benchmark.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        max_dist=config["sanger_max_dist"],
        blacklist_arg=lambda wildcards, input: (
            f"--blacklist {input.blacklist}" if input.blacklist else ""
        ),
        chromsizes_arg=(
            f"--chromsizes {config['chrom_sizes_path_no_cassette']}"
            if config.get("chrom_sizes_path_no_cassette")
            else ""
        ),
    shell:
        """
        python3 {input.script} --sites {input.sites} \
            --max-dist {params.max_dist} {params.blacklist_arg} \
            {params.chromsizes_arg} \
            -o {output.sites} --output-for-ucsc {output.for_ucsc} \
            >{log[0]} 2>&1
        """


rule filter_confirmed_sites:
    input:
        sites=f"{sanger_folder}/all_sanger_sites.bed",
        clones=f"{stats_folder}/sanger_clone_summary.tsv",
        chromsizes=config.get("chrom_sizes_path_no_cassette", []),
        sanger_vs_ngs=f"{validation_folder}/sanger_vs_ngs.tsv" if do_validation else [],
        script=f"{scripts_dir}/filter_confirmed_sites.py",
    output:
        sites=f"{sanger_folder}/confirmed_sanger_sites.bed",
        for_ucsc=f"{sanger_folder}/confirmed_sanger_sites_for_ucsc.bed",
        no_cassette=f"{sanger_folder}/confirmed_sanger_sites_no_cassette.bed",
        region=f"{sanger_folder}/confirmed_sanger_sites_region.bed",
        deduplicated=f"{sanger_folder}/confirmed_sanger_sites_deduplicated.bed",
        region_deduplicated=f"{sanger_folder}/confirmed_sanger_sites_region_deduplicated.bed",
    log:
        "logs/filter_confirmed_sites/log.log",
    benchmark:
        "benchmarks/filter_confirmed_sites/benchmark.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        chromsizes_arg=lambda wildcards, input: (
            f"--chromsizes {input.chromsizes}" if input.chromsizes else ""
        ),
        region_arg=lambda wildcards: (
            f"--region {config['validated_clones_region']}"
            if config.get("validated_clones_region")
            else ""
        ),
        sanger_vs_ngs_arg=lambda wildcards, input: (
            f"--sanger-vs-ngs {input.sanger_vs_ngs}" if input.sanger_vs_ngs else ""
        ),
    shell:
        """
        python3 {input.script} --sites {input.sites} --clone-summary {input.clones} \
            {params.chromsizes_arg} {params.region_arg} {params.sanger_vs_ngs_arg} \
            -o {output.sites} --output-for-ucsc {output.for_ucsc} \
            --output-no-cassette {output.no_cassette} \
            --output-region {output.region} \
            --output-deduplicated {output.deduplicated} \
            --output-region-deduplicated {output.region_deduplicated} \
            >{log[0]} 2>&1
        """


rule sanger_stats:
    input:
        reads=expand(f"{sanger_folder}/{{sample}}_reads.tsv", sample=sanger_sample_list),
        sites=f"{sanger_folder}/all_sanger_sites.bed",
        validation=f"{validation_folder}/sanger_vs_ngs.tsv" if do_validation else [],
        original_site=original_site_file if has_original_site else [],
        script=f"{scripts_dir}/sanger_stats.py",
    output:
        qc=f"{stats_folder}/sanger_qc_stats.tsv",
        clones=f"{stats_folder}/sanger_clone_summary.tsv",
        positions=f"{stats_folder}/sanger_positions.tsv",
        position_counts=f"{stats_folder}/sanger_position_counts.tsv",
        read_qc=f"{stats_folder}/sanger_read_qc.tsv",
    log:
        "logs/sanger_stats/log.log",
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        validation_arg=lambda wildcards, input: (
            f"--validation {input.validation}" if input.validation else ""
        ),
        original_site_arg=lambda wildcards, input: (
            f"--original-site {input.original_site}" if input.original_site else ""
        ),
        max_dist=config["original_insertion_max_dist"],
        region_arg=lambda wildcards: (
            f"--region {config['validated_clones_region']}"
            if config.get("validated_clones_region")
            else ""
        ),
    shell:
        """
        python3 {input.script} --reads {input.reads} --sites {input.sites} \
            {params.validation_arg} \
            {params.original_site_arg} --original-max-dist {params.max_dist} \
            {params.region_arg} \
            --output-qc {output.qc} \
            --output-clone-summary {output.clones} \
            --output-positions {output.positions} \
            --output-position-counts {output.position_counts} \
            --output-read-qc {output.read_qc} \
            >{log[0]} 2>&1
        """
