rule mapping:
    input:
        fqs = fq_assembly_input,
        idx = expand(result_dir / "final_references" / "{{sample}}_idx" / "genome{end}", end=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"])
    output:
        result_dir / "sam" / "{sample}.sam"
    log:
        result_dir / "log" / "bowtie2_mapping" / "{sample}.log"
    threads:
        10
    params:
        min_frag_len = config["bowtie2_mapping"]['I'],
        max_frag_len = config["bowtie2_mapping"]['X'],
        idx_base = lambda wildcards: result_dir / "final_references" / f"{wildcards.sample}_idx" / "genome"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2 -x {params.idx_base} -1 {input.fqs[0]} -2 {input.fqs[1]} -p {threads} -X {params.max_frag_len} -I {params.min_frag_len} -S {output} > {log} 2>&1"

rule genome_index:
    input:
        result_dir / "final_references" / "{sample}_ref.fasta"
    output:
        expand(result_dir / "final_references" / "{{sample}}_idx" / "genome{end}", end=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"])
    log:
        result_dir / "log" / "reference_indexing" / "{sample}_ref_idx.log"
    threads:
        10
    params:
        idx_base = lambda wildcards: result_dir / "final_references" / f"{wildcards.sample}_idx" / "genome"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2-build {reference} {params.idx_base} -p {threads} > {log} 2>&1"