localrules:
    combine_sanger_sites,
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
            "../envs/minibwa.yaml"
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
        primer_positions=config["primer_position_file"],
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
    shell:
        """
        python3 {input.script} --sites {input.sites} \
            --max-dist {params.max_dist} {params.blacklist_arg} \
            -o {output.sites} --output-for-ucsc {output.for_ucsc} \
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
    shell:
        """
        python3 {input.script} --reads {input.reads} --sites {input.sites} \
            {params.validation_arg} \
            {params.original_site_arg} --original-max-dist {params.max_dist} \
            --output-qc {output.qc} \
            --output-clone-summary {output.clones} \
            --output-positions {output.positions} \
            >{log[0]} 2>&1
        """
