rule mapping:
    input:
        fq1 = lambda wildcards: samples.at[wildcards.sample, "fq1"],
        fq2 = lambda wildcards: samples.at[wildcards.sample, "fq2"],
        idx = multiext(
            ref_base,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        result_dir / "sam" / "{sample}.sam"
    threads:
        4
    conda:
        Path.cwd() / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2 -x {ref_base} -1 {input.fq1} -2 {input.fq2} -p {threads} -S {output}"

rule genome_index:
    input:
        reference
    output:
        idx=multiext(
            ref_base,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads:
        4
    conda:
        Path.cwd() / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2-build {reference} {ref_base} -p {threads}"