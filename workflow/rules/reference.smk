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
        "../envs/pyfastx.yaml"
    shell:
        """
        python3 -c "import pyfastx, sys; pyfastx.Fasta(sys.argv[1], index_file=sys.argv[2])" \
            {input} {output} >{log[0]} 2>&1
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
        "../envs/pyfastx.yaml"
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
