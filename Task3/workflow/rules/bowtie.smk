rule mapping:
    input:
        fq1 = lambda wildcards: samples.at[wildcards.sample, "fq1"],
        fq2 = lambda wildcards: samples.at[wildcards.sample, "fq2"],
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
        4
    params:
        min_frag_len = config["bowtie2_mapping"]['I'],
        max_frag_len = config["bowtie2_mapping"]['X'],
        idx_base = result_dir / "ref_idx" / "reference"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2 -x {params.idx_base} -1 {input.fq1} -2 {input.fq2} -p {threads} -X {params.max_frag_len} -I {params.min_frag_len} -S {output} 2> {log}"

rule genome_index:
    input:
        reference
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
        result_dir / "log" / "bowtie2_indexing.log"
    threads:
        4
    params:
        idx_base = result_dir / "ref_idx" / "reference"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2-build {reference} {params.idx_base} -p {threads} 2> {log}"