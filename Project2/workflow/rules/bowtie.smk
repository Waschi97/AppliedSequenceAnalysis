rule mapping:
    input:
        fqs = fq_assembly_input,
        idx = multiext(
            str(result_dir / "ref_idx" / "reference"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        result_dir / "sam" / "{sample}.sam"
    log:
        result_dir / "log" / "bowtie2_mapping" / "{sample}.log"
    threads:
        10
    params:
        min_frag_len = config["bowtie2_mapping"]['I'],
        max_frag_len = config["bowtie2_mapping"]['X'],
        idx_base = result_dir / "ref_idx" / "reference"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2 -x {params.idx_base} -1 {input.fqs[0]} -2 {input.fqs[1]} -p {threads} -X {params.max_frag_len} -I {params.min_frag_len} -S {output} > {log} 2>&1"

rule genome_index:
    input:
        result_dir / "best_reference" / "reference.fasta"
    output:
        idx=multiext(
            str(result_dir / "ref_idx" / "reference"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        result_dir / "log" / "reference_indexing" / "bowtie2_indexing.log"
    threads:
        10
    params:
        idx_base = result_dir / "ref_idx" / "reference"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2-build {reference} {params.idx_base} -p {threads} > {log} 2>&1"