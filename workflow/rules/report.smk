localrules:
    report,


rule report:
    input:
        script=f"{scripts_dir}/report.py",
        ngs_mapping=f"{stats_folder}/ngs_mapping_stats.tsv" if sample_list else [],
        ngs_sidedness=f"{stats_folder}/ngs_site_sidedness.tsv" if sample_list else [],
        sanger_reads=f"{stats_folder}/sanger_read_stats.tsv" if sanger_sample_list else [],
        sanger_fail_reasons=(
            f"{stats_folder}/sanger_fail_reasons.tsv" if sanger_sample_list else []
        ),
        sanger_clones=f"{stats_folder}/sanger_clone_summary.tsv" if sanger_sample_list else [],
        validation=f"{stats_folder}/validation_summary.tsv" if do_validation else [],
    output:
        f"{stats_folder}/report.md",
    log:
        "logs/report/log.log",
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        ngs_mapping_arg=lambda wildcards, input: (
            f"--ngs-mapping {input.ngs_mapping}" if input.ngs_mapping else ""
        ),
        ngs_sidedness_arg=lambda wildcards, input: (
            f"--ngs-sidedness {input.ngs_sidedness}" if input.ngs_sidedness else ""
        ),
        sanger_reads_arg=lambda wildcards, input: (
            f"--sanger-reads {input.sanger_reads}" if input.sanger_reads else ""
        ),
        sanger_fail_reasons_arg=lambda wildcards, input: (
            f"--sanger-fail-reasons {input.sanger_fail_reasons}"
            if input.sanger_fail_reasons
            else ""
        ),
        sanger_clones_arg=lambda wildcards, input: (
            f"--sanger-clones {input.sanger_clones}" if input.sanger_clones else ""
        ),
        validation_arg=lambda wildcards, input: (
            f"--validation {input.validation}" if input.validation else ""
        ),
    shell:
        """
        python3 {input.script} {params.ngs_mapping_arg} {params.ngs_sidedness_arg} \
            {params.sanger_reads_arg} {params.sanger_fail_reasons_arg} \
            {params.sanger_clones_arg} {params.validation_arg} \
            -o {output} \
            >{log[0]} 2>&1
        """
