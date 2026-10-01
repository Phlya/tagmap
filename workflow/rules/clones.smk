localrules:
    ngs_clone_stats,


rule ngs_clone_stats:
    input:
        support=f"{insertion_sites_folder}/all_sites_support.tsv",
        ngs_qc=f"{stats_folder}/ngs_qc_stats.tsv",
        original_site=original_site_file if has_original_site else [],
        script=f"{scripts_dir}/ngs_clone_stats.py",
    output:
        f"{stats_folder}/ngs_clone_summary.tsv",
    log:
        "logs/ngs_clone_stats/log.log",
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        samples=" ".join(ngs_sample_map()),
        min_reads=config["clone_min_reads"],
        min_molecules=config["clone_min_molecules"],
        max_contamination_frac=config["clone_max_contamination_frac"],
        multi_site_frac=config["clone_multi_site_frac"],
        merge_dist=config["clone_merge_dist"],
        shared_dist=config["max_dist_between_sides"],
        require_both_sides_arg=(
            "--require-both-sides" if config["clone_require_both_sides"] else ""
        ),
        original_site_arg=lambda wildcards, input: (
            f"--original-site {input.original_site} "
            f"--original-max-dist {config['original_insertion_max_dist']}"
            if input.original_site
            else ""
        ),
    shell:
        """
        python3 {input.script} --support {input.support} --ngs-qc {input.ngs_qc} \
            --samples {params.samples} \
            --min-reads {params.min_reads} --min-molecules {params.min_molecules} \
            --max-contamination-frac {params.max_contamination_frac} \
            --multi-site-frac {params.multi_site_frac} \
            --merge-dist {params.merge_dist} --shared-dist {params.shared_dist} \
            {params.require_both_sides_arg} {params.original_site_arg} \
            -o {output} >{log[0]} 2>&1
        """
