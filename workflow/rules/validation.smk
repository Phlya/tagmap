localrules:
    compare_sanger_ngs,


rule compare_sanger_ngs:
    input:
        sanger=f"{sanger_folder}/all_sanger_sites.bed",
        ngs_sites=f"{insertion_sites_folder}/all_sites.bed",
        ngs_peaks=f"{peaks_folder}/all_peaks.bed",
        script=f"{scripts_dir}/compare_sanger_ngs.py",
    output:
        table=f"{validation_folder}/sanger_vs_ngs.tsv",
        confirmed=f"{validation_folder}/confirmed_sites.bed",
    log:
        "logs/compare_sanger_ngs/log.log",
    benchmark:
        "benchmarks/compare_sanger_ngs/benchmark.tsv"
    conda:
        "../envs/sanger.yaml"
    threads: 1
    params:
        max_dist=config["validation_max_dist"],
        pairs_arg=lambda wildcards: (
            f"--sample-pairs {' '.join(sanger_ngs_pairs())}"
            if sanger_ngs_pairs()
            else ""
        ),
    shell:
        """
        python3 {input.script} --sanger {input.sanger} \
            --ngs-sites {input.ngs_sites} --ngs-peaks {input.ngs_peaks} \
            --max-dist {params.max_dist} {params.pairs_arg} \
            -o {output.table} --output-confirmed {output.confirmed} \
            >{log[0]} 2>&1
        """
