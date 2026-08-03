if config["mapper"] == "minibwa":

    rule minibwa_index:
        input:
            refgen_path,
        output:
            idx,
        log:
            "logs/minibwa_index/{}.log".format(
                os.path.splitext(os.path.basename(refgen_path))[0]
            ),
        priority: 100
        conda:
            "../envs/minibwa.yaml"
        threads: 8
        shell:
            """
            minibwa index -t {threads} {input} >{log[0]} 2>&1
            """

else:

    rule bwaindex:
        input:
            refgen_path,
        output:
            idx,
        log:
            "logs/bwa-memx_index/{}.log".format(
                os.path.splitext(os.path.basename(refgen_path))[0]
            ),
        priority: 100
        threads: 8  # Only affects bwa-meme
        params:
            bwa=config["mapper"],
        wrapper:
            "v3.8.0/bio/bwa-memx/index"


rule fasta_index:
    input:
        refgen_path,
    output:
        config["fasta_index_file"],
    log:
        "logs/fasta_index/log.log",
    conda:
        "../envs/all.yaml"
    shell:
        """
        python3 -c "import pysam, sys; pysam.faidx(sys.argv[1], '-o', sys.argv[2])" \
            {input} {output} >{log[0]} 2>&1
        """


rule make_chromsizes:
    input:
        refgen_path=refgen_path,
    output:
        chrom_sizes_path=config["chrom_sizes_path"],
    log:
        "logs/make_chromsizes/log.log",
    conda:
        "../envs/all.yaml"
    threads: 8
    shell:
        """
        chromsize --sequence {input.refgen_path} \
            -o $(dirname {output.chrom_sizes_path}) \
            -p $(basename {output.chrom_sizes_path}) \
            -t {threads} \
            >{log[0]} 2>&1
        """


rule find_original_site:
    input:
        script=f"{scripts_dir}/find_original_site.py",
        refgen_path=refgen_path,
    output:
        original_site_file,
    log:
        "logs/find_original_site/log.log",
    conda:
        "../envs/all.yaml"
    params:
        cassette_name=config["cassette_name"],
        site_arg=(
            f"--site {config['original_insertion_site']}"
            if config.get("original_insertion_site")
            else ""
        ),
        upstream_seq_arg=(
            f"--upstream-seq {config['original_insertion_upstream_seq']}"
            if config.get("original_insertion_upstream_seq")
            else ""
        ),
    shell:
        """
        python3 {input.script} --genome {input.refgen_path} \
            --construct-contigs {params.cassette_name} \
            {params.site_arg} {params.upstream_seq_arg} \
            -o {output} \
            >{log[0]} 2>&1
        """


rule get_primer_positions:
    input:
        script=f"{scripts_dir}/get_primer_positions.py",
        refgen_path=refgen_path,
    output:
        config["primer_position_file"],
    log:
        "logs/get_primer_positions/log.log",
    conda:
        "../envs/all.yaml"
    params:
        cassette_name=config["cassette_name"],
        forward_primer=config["forward_primer_sequence"],
        reverse_primer=config["reverse_primer_sequence"],
    shell:
        """
        python3 {input.script} --genome {input.refgen_path} \
            --cassette-name {params.cassette_name} \
            --forward-primer {params.forward_primer} \
            --reverse-primer {params.reverse_primer} \
            -o {output} \
            >{log[0]} 2>&1
        """
