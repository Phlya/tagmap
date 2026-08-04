localrules:
    report,
    report_pdf,


rule report:
    input:
        script=f"{scripts_dir}/report.py",
        ngs_qc=f"{stats_folder}/ngs_qc_stats.tsv" if sample_list else [],
        sanger_qc=f"{stats_folder}/sanger_qc_stats.tsv" if sanger_sample_list else [],
        sanger_clones=f"{stats_folder}/sanger_clone_summary.tsv"
        if sanger_sample_list
        else [],
        sanger_positions=f"{stats_folder}/sanger_positions.tsv"
        if sanger_sample_list
        else [],
        validation=f"{stats_folder}/validation_summary.tsv" if do_validation else [],
    output:
        f"{stats_folder}/report.md",
    log:
        "logs/report/log.log",
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        ngs_qc_arg=lambda wildcards, input: (
            f"--ngs-qc {input.ngs_qc}" if input.ngs_qc else ""
        ),
        sanger_qc_arg=lambda wildcards, input: (
            f"--sanger-qc {input.sanger_qc}" if input.sanger_qc else ""
        ),
        sanger_clones_arg=lambda wildcards, input: (
            f"--sanger-clones {input.sanger_clones}" if input.sanger_clones else ""
        ),
        sanger_positions_arg=lambda wildcards, input: (
            f"--sanger-positions {input.sanger_positions}"
            if input.sanger_positions
            else ""
        ),
        validation_arg=lambda wildcards, input: (
            f"--validation {input.validation}" if input.validation else ""
        ),
    shell:
        """
        python3 {input.script} {params.ngs_qc_arg} {params.sanger_qc_arg} \
            {params.sanger_clones_arg} {params.sanger_positions_arg} \
            {params.validation_arg} \
            -o {output} \
            >{log[0]} 2>&1
        """


rule report_pdf:
    input:
        script=f"{scripts_dir}/report_pdf.py",
        # report_pdf.py imports SECTIONS from report.py, so it depends on it too.
        report_script=f"{scripts_dir}/report.py",
        ngs_qc=f"{stats_folder}/ngs_qc_stats.tsv" if sample_list else [],
        sanger_qc=f"{stats_folder}/sanger_qc_stats.tsv" if sanger_sample_list else [],
        sanger_clones=f"{stats_folder}/sanger_clone_summary.tsv"
        if sanger_sample_list
        else [],
        sanger_positions=f"{stats_folder}/sanger_positions.tsv"
        if sanger_sample_list
        else [],
        validation=f"{stats_folder}/validation_summary.tsv" if do_validation else [],
    output:
        f"{stats_folder}/report.pdf",
    log:
        "logs/report_pdf/log.log",
    conda:
        "../envs/all.yaml"
    threads: 1
    params:
        ngs_qc_arg=lambda wildcards, input: (
            f"--ngs-qc {input.ngs_qc}" if input.ngs_qc else ""
        ),
        sanger_qc_arg=lambda wildcards, input: (
            f"--sanger-qc {input.sanger_qc}" if input.sanger_qc else ""
        ),
        sanger_clones_arg=lambda wildcards, input: (
            f"--sanger-clones {input.sanger_clones}" if input.sanger_clones else ""
        ),
        sanger_positions_arg=lambda wildcards, input: (
            f"--sanger-positions {input.sanger_positions}"
            if input.sanger_positions
            else ""
        ),
        validation_arg=lambda wildcards, input: (
            f"--validation {input.validation}" if input.validation else ""
        ),
    shell:
        """
        python3 {input.script} {params.ngs_qc_arg} {params.sanger_qc_arg} \
            {params.sanger_clones_arg} {params.sanger_positions_arg} \
            {params.validation_arg} \
            -o {output} \
            >{log[0]} 2>&1
        """
