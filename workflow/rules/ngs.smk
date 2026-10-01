localrules:
    combine_peaks,
    combine_all_peaks,


def library_fastqs(wildcards):
    """The fastq files of one library as listed in samples_path - R1 or R2
    according to the read wildcard (the sheet's columns are fastq1/fastq2)."""
    return samples.loc[
        samples["name"] == wildcards.source, f"fastq{wildcards.read[-1]}"
    ]


if demultiplex:

    # Each library is merged first and then split by in-read barcode into one
    # sample per plate, named "{plate}_{name}" - which is what everything below
    # sees, so from here on a plate of a library is an ordinary NGS sample.

    rule merge_library_fastq:
        input:
            library_fastqs,
        output:
            temp(f"{fastq_folder}/raw/{{source}}.{{read}}.fastq.gz"),
        log:
            "logs/merge_fastq/{source}_{read}.log",
        conda:
            "../envs/all.yaml"
        shell:
            """
            cat {input} >{output} 2>{log[0]}
            """

    rule demux:
        input:
            r1=f"{fastq_folder}/raw/{{source}}.R1.fastq.gz",
            r2=f"{fastq_folder}/raw/{{source}}.R2.fastq.gz",
            barcodes=config["barcodes_path"],
            script=f"{scripts_dir}/demux_barcodes.py",
            # tagmaplib validates the barcode sheet, so a change to it reruns
            # the split too.
            lib=f"{scripts_dir}/tagmaplib.py",
        output:
            r1=temp([f"{fastq_folder}/{plate}_{{source}}.R1.fastq.gz" for plate in all_plates]),
            r2=temp([f"{fastq_folder}/{plate}_{{source}}.R2.fastq.gz" for plate in all_plates]),
            stats=f"{stats_folder}/demux/{{source}}_demux.tsv",
        log:
            "logs/demux/{source}.log",
        benchmark:
            "benchmarks/demux/{source}.tsv"
        conda:
            "../envs/all.yaml"
        threads: 1
        params:
            plates=" ".join(all_plates),
            max_mismatch=config["barcode_max_mismatch"],
            barcode_reads=config["barcode_reads"],
            primer_args=(
                f"--require-primer --forward-primer {config.get('forward_primer_sequence')} "
                f"--reverse-primer {config.get('reverse_primer_sequence')} "
                f"--primer-max-mismatch {config['barcode_primer_max_mismatch']}"
                if config["barcode_require_primer"]
                else ""
            ),
        shell:
            """
            python3 {input.script} --r1 {input.r1} --r2 {input.r2} \
                --barcodes {input.barcodes} --source {wildcards.source} \
                --plates {params.plates} \
                --out-r1 {output.r1} --out-r2 {output.r2} \
                --stats {output.stats} \
                --max-mismatch {params.max_mismatch} \
                --barcode-reads {params.barcode_reads} \
                {params.primer_args} \
                >{log[0]} 2>&1
            """

    rule combine_demux_stats:
        input:
            expand(f"{stats_folder}/demux/{{source}}_demux.tsv", source=source_list),
        output:
            f"{stats_folder}/demux_stats.tsv",
        log:
            "logs/combine_demux_stats/log.log",
        conda:
            "../envs/all.yaml"
        shell:
            """
            awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} >{output} 2>{log[0]}
            """

    localrules:
        combine_demux_stats,

else:

    rule merge_fastq:
        input:
            lambda wildcards: samples.loc[
                samples["name"] == wildcards.sample, f"fastq{wildcards.read[-1]}"
            ],
        output:
            temp(f"{fastq_folder}/{{sample}}.{{read}}.fastq.gz"),
        log:
            "logs/merge_fastq/{sample}_{read}.log",
        conda:
            "../envs/all.yaml"
        shell:
            """
            cat {input} >{output} 2>{log[0]}
            """


rule fastqc:
    input:
        f"{fastq_folder}/{{sample}}.{{read}}.fastq.gz",
    output:
        html=f"{fastq_folder}/{{sample}}.{{read}}_fastqc.html",
        zip=f"{fastq_folder}/{{sample}}.{{read}}_fastqc.zip",
    log:
        "logs/fastqc/{sample}_{read}.log",
    benchmark:
        "benchmarks/fastqc/{sample}_{read}.tsv"
    threads: 1
    params:
        extra="--quiet",
    wrapper:
        "v3.9.0/bio/fastqc"


rule trim:
    input:
        sample=[
            f"{fastq_folder}/{{sample}}.R1.fastq.gz",
            f"{fastq_folder}/{{sample}}.R2.fastq.gz",
        ],
    output:
        trimmed=[
            f"{fastq_folder}/{{sample}}_trimmed.R1.fastq.gz",
            f"{fastq_folder}/{{sample}}_trimmed.R2.fastq.gz",
        ],
        json=f"{fastq_folder}/{{sample}}.fastp.json",
        html=f"{fastq_folder}/{{sample}}.fastp.html",
    log:
        "logs/fastp/{sample}.log",
    threads: 2
    params:
        extra=config["trim_args"],
    wrapper:
        "v3.9.0/bio/fastp"


if config["mapper"] == "minibwa":

    rule minibwa_map:
        input:
            reads=(
                [
                    f"{fastq_folder}/{{sample}}_trimmed.R1.fastq.gz",
                    f"{fastq_folder}/{{sample}}_trimmed.R2.fastq.gz",
                ]
                if config["trim"]
                else [
                    f"{fastq_folder}/{{sample}}.R1.fastq.gz",
                    f"{fastq_folder}/{{sample}}.R2.fastq.gz",
                ]
            ),
            reference=refgen_path,
            idx=idx,
        output:
            f"{bams_folder}/{{sample}}.bam",
        log:
            "logs/minibwa/{sample}.log",
        benchmark:
            "benchmarks/minibwa/{sample}.tsv"
        conda:
            "../envs/all.yaml"
        threads: 12
        params:
            # --hic: independent mate mapping, no proper-pair assumption
            # (equivalent to -5P) - matches the intent of the other mappers'
            # -SP, since these are chimeric tagmentation fragments rather
            # than conventional inserts.
            extra="--hic -s 30",
        shell:
            """
            minibwa map {params.extra} -t {threads} {input.reference} {input.reads} 2>{log[0]} \
                | samtools view -b -o {output} - 2>>{log[0]}
            """

else:

    rule bwamap:
        input:
            reads=(
                [
                    f"{fastq_folder}/{{sample}}_trimmed.R1.fastq.gz",
                    f"{fastq_folder}/{{sample}}_trimmed.R2.fastq.gz",
                ]
                if config["trim"]
                else [
                    f"{fastq_folder}/{{sample}}.R1.fastq.gz",
                    f"{fastq_folder}/{{sample}}.R2.fastq.gz",
                ]
            ),
            reference=refgen_path,
            idx=idx,
        output:
            f"{bams_folder}/{{sample}}.bam",
        log:
            "logs/bwa_memx/{sample}.log",
        benchmark:
            "benchmarks/bwa_memx/{sample}.tsv"
        threads: 12
        params:
            bwa=config["mapper"],
            sort="none",
            dedup="none",
            # Lower minimal alignment score for bwa-mem to increase sensitivity for
            # short reads and short alignments
            extra="-SP -T 30",
        wrapper:
            "v3.3.3/bio/bwa-memx/mem"


rule parse2:
    input:
        bam=f"{bams_folder}/{{sample}}.bam",
        chromsizes=config["chrom_sizes_path"],
    output:
        pairs=f"{pairs_folder}/{{sample}}_sorted.pairs",
    log:
        "logs/parse2/{sample}.log",
    benchmark:
        "benchmarks/parse2/{sample}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 8
    params:
        # Every flag but the input, which differs for an empty BAM (below).
        parse2="--drop-sam --flip --min-mapq 1 --max-insert-size 5000 "
        "--add-columns pos5,pos3,read_len,mapq --add-pair-index "
        "--report-position junction --report-orientation pair",
    shell:
        """
        # pairtools parse2 crashes on a BAM with no alignments at all (it never
        # works out its own column list), which is what a sample with no reads
        # is - e.g. a plate with nothing in one well, once demultiplexed.
        # Giving it one unmapped placeholder pair makes it write the same
        # header as for any other sample, and the placeholder is then dropped,
        # leaving a valid, empty pairs file for the rest of the workflow.
        if [ "$(samtools view -c {input.bam})" = "0" ]; then
            ( samtools view -h {input.bam}; \
              printf 'placeholder\t77\t*\t0\t0\t*\t*\t0\t0\tN\tI\n' ; \
              printf 'placeholder\t141\t*\t0\t0\t*\t*\t0\t0\tN\tI\n' ) \
                | pairtools parse2 -c {input.chromsizes} {params.parse2} - 2>{log[0]} \
                | pairtools select 'False' 2>>{log[0]} \
                | pairtools sort --nproc {threads} -o {output.pairs} >>{log[0]} 2>&1
        else
            pairtools parse2 -c {input.chromsizes} {params.parse2} {input.bam} \
                | pairtools sort --nproc {threads} -o {output.pairs} \
                    >{log[0]} 2>&1
        fi
        """


if config["dedup"]:

    ruleorder: dedup > stats

else:

    ruleorder: stats > dedup


rule stats:
    input:
        pairs=f"{pairs_folder}/{{sample}}_sorted.pairs",
        primer_positions=config["primer_position_file"],
        chromsizes=config["chrom_sizes_path"],
    output:
        stats=f"{pairs_folder}/{{sample}}_stats.yml",
    log:
        "logs/stats/{sample}.log",
    benchmark:
        "benchmarks/stats/{sample}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        filter_forward=lambda wildcards, input: get_filter(
            "forward", input.primer_positions
        ),
        filter_reverse=lambda wildcards, input: get_filter(
            "reverse", input.primer_positions
        ),
    shell:
        """
        pairtools stats --engine python --yaml \
            --filter 'forward:{params.filter_forward}' \
            --filter 'reverse:{params.filter_reverse}' \
            -o {output.stats} {input.pairs} \
            >{log[0]} 2>&1
        """


rule dedup:
    input:
        pairs=f"{pairs_folder}/{{sample}}_sorted.pairs",
        primer_positions=config["primer_position_file"],
        chromsizes=config["chrom_sizes_path"],
    output:
        pairs=f"{pairs_folder}/{{sample}}_dupmarked.pairs",
        stats=f"{pairs_folder}/{{sample}}_stats.yml",
    log:
        "logs/dedup/{sample}.log",
    benchmark:
        "benchmarks/dedup/{sample}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        filter_forward=lambda wildcards, input: get_filter(
            "forward", input.primer_positions
        ),
        filter_reverse=lambda wildcards, input: get_filter(
            "reverse", input.primer_positions
        ),
    shell:
        """
        pairtools dedup --backend cython --mark-dups --max-mismatch 2 \
            --p1 pos51 --p2 pos31 \
            --engine python --yaml \
            --filter 'forward:{params.filter_forward}' \
            --filter 'reverse:{params.filter_reverse}' \
            --output-stats {output.stats} \
            --output-dups {output.pairs} -o {output.pairs} {input.pairs} \
            >{log[0]} 2>&1
        """


rule get_trans_side_pairs:
    input:
        pairs=(
            f"{pairs_folder}/{{sample}}_dupmarked.pairs"
            if config["dedup"]
            else f"{pairs_folder}/{{sample}}_sorted.pairs"
        ),
        primer_positions=config["primer_position_file"],
        chromsizes=config["chrom_sizes_path"],
    output:
        f"{pairs_folder}/{{sample}}_{{side}}.pairs",
    log:
        "logs/get_trans_pairs/{sample}_{side}.log",
    benchmark:
        "benchmarks/get_trans_pairs/{sample}_{side}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        filter=lambda wildcards, input: get_filter(
            wildcards.side, input.primer_positions
        ),
    shell:
        """
        pairtools select -t pos51 int -t pos31 int -t pos52 int -t pos32 int \
            -t read_len1 int -t read_len2 int \
            '{params.filter}' \
            -o {output} {input.pairs} \
            >{log[0]} 2>&1
        """


rule coverage:
    input:
        pairs=f"{pairs_folder}/{{sample}}_{{side}}.pairs",
        script=f"{scripts_dir}/coverage.py",
    output:
        bg=f"{coverage_folder}/{{sample}}_{{side}}_coverage.bedgraph",
        bw=f"{coverage_folder}/{{sample}}_{{side}}_coverage.bw",
    log:
        "logs/coverage/{sample}_{side}.log",
    benchmark:
        "benchmarks/coverage/{sample}_{side}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        side=1,
    shell:
        """
        python3 {input.script} -i {input.pairs} --side {params.side} \
            --itr-side {wildcards.side} -o {output.bg} \
            -t {threads} --output-bigwig {output.bw} \
            >{log[0]} 2>&1
        """


if config["peak_caller"] == "junction_tiered":

    rule find_peaks:
        input:
            pairs=f"{pairs_folder}/{{sample}}_{{side}}.pairs",
            genome=refgen_path,
            genome_index=config["fasta_index_file"],
            script=f"{scripts_dir}/find_peaks_junctions.py",
        output:
            peaks=f"{peaks_folder}/{{sample}}_{{side}}.bed",
            evidence=f"{peaks_folder}/{{sample}}_{{side}}_evidence.tsv",
        log:
            "logs/find_peaks/{sample}_{side}.log",
        benchmark:
            "benchmarks/find_peaks/{sample}_{side}.tsv"
        conda:
            "../envs/all.yaml"
        threads: 1
        params:
            junction_jitter=config["junction_jitter"],
            max_fragment=config["max_fragment"],
            min_junction_frags=config["min_junction_frags"],
            min_junction_mapq=config["min_junction_mapq"],
            motif_max_support=config["motif_max_support"],
            insertion_seq=config["insertion_seq"],
            min_peak_width=config["min_peak_width"],
            min_peak_reads=config["min_peak_reads"],
            min_peak_frac=config["min_peak_frac"],
            min_peak_dist=config["min_peak_dist"],
            min_peak_positions=config["min_peak_positions"],
            cassette_name=config["cassette_name"],
        shell:
            """
            python3 {input.script} -i {input.pairs} --side {wildcards.side} \
                -o {output.peaks} --output-detail {output.evidence} \
                --junction-jitter {params.junction_jitter} \
                --max-fragment {params.max_fragment} \
                --min-junction-frags {params.min_junction_frags} \
                --min-junction-mapq {params.min_junction_mapq} \
                --motif-max-support {params.motif_max_support} \
                --insertion-seq {params.insertion_seq} \
                --genome {input.genome} --genome-index {input.genome_index} \
                --min-peak-width {params.min_peak_width} --min-peak-frac {params.min_peak_frac} \
                --min-peak-reads {params.min_peak_reads} --min-peak-dist {params.min_peak_dist} \
                --min-peak-positions {params.min_peak_positions} \
                --ignore-chrom {params.cassette_name} \
                >{log[0]} 2>&1
            """

else:

    rule find_peaks:
        input:
            coverage=f"{coverage_folder}/{{sample}}_{{side}}_coverage.bedgraph",
            script=f"{scripts_dir}/find_peaks.py",
        output:
            f"{peaks_folder}/{{sample}}_{{side}}.bed",
        log:
            "logs/find_peaks/{sample}_{side}.log",
        benchmark:
            "benchmarks/find_peaks/{sample}_{side}.tsv"
        conda:
            "../envs/all.yaml"
        threads: 1
        params:
            cluster_arg=("--no-cluster" if only_read_junctions else "--cluster"),
            auto_li_arg=("--auto-li" if only_read_junctions else "--no-auto-li"),
            min_peak_width=config["min_peak_width"],
            min_peak_reads=config["min_peak_reads"],
            min_peak_frac=config["min_peak_frac"],
            min_peak_dist=config["min_peak_dist"],
            min_peak_positions=config["min_peak_positions"],
            cassette_name=config["cassette_name"],
        shell:
            """
            python3 {input.script} -i {input.coverage} -o {output} \
                {params.cluster_arg} {params.auto_li_arg} \
                --min-peak-width {params.min_peak_width} --min-peak-frac {params.min_peak_frac} \
                --min-peak-reads {params.min_peak_reads} --min-peak-dist {params.min_peak_dist} \
                --min-peak-positions {params.min_peak_positions} \
                --ignore-chrom {params.cassette_name} \
                >{log[0]} 2>&1
            """


rule combine_peaks:
    input:
        peaks_fwd=f"{peaks_folder}/{{sample}}_forward.bed",
        peaks_rev=f"{peaks_folder}/{{sample}}_reverse.bed",
        blacklist=config.get("blacklist", []),
        script=f"{scripts_dir}/combine_peaks.py",
    output:
        peaks=f"{peaks_folder}/{{sample}}_peaks.bed",
        genome_browser=f"{peaks_folder}/{{sample}}_peaks_for_genome_browser.bed",
    log:
        "logs/combine_peaks/{sample}.log",
    benchmark:
        "benchmarks/combine_peaks/{sample}.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        blacklist_arg=lambda wildcards, input: (
            f"--blacklist {input.blacklist}" if input.blacklist else ""
        ),
    shell:
        """
        python3 {input.script} --fwd {input.peaks_fwd} --rev {input.peaks_rev} \
            {params.blacklist_arg} \
            --sample-name {wildcards.sample} -o {output.peaks} \
            --output-genome-browser {output.genome_browser} >{log[0]} 2>&1
        """


rule combine_all_peaks:
    input:
        peaks=expand(f"{peaks_folder}/{{sample}}_peaks.bed", sample=sample_list),
    output:
        f"{peaks_folder}/all_peaks.bed",
    log:
        "logs/combine_all_peaks/log.log",
    benchmark:
        "benchmarks/combine_all_peaks/benchmark.tsv"
    conda:
        "../envs/all.yaml"
    threads: 1
    shell:
        """
        cat {input.peaks} | sort -k4,4V -k1,1V -k2,2n -k3,3n -k6,6 >{output}
        """
