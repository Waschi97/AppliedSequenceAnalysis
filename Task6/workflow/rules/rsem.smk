rule calculate_expression:
    input:
        reads = fq_assembly_input,
        ref = result_dir / "denovo" / "denovo_index" / "trinity.seq"
    output:
        result_dir / "denovo" / "gene_expression" / "{sample}.genes.results",
        result_dir / "denovo" / "gene_expression" / "{sample}.isoforms.results"
    log:
        result_dir / "logs" / "denovo" / "rsem" / "{sample}_expression.log"
    conda:
        Path("..") / "envs" / "rsem_env.yaml"
    threads:
        10
    params:
        base_ref = result_dir / "denovo" / "denovo_index" / "trinity",
        base_out = lambda wildcards: result_dir / "denovo" / "gene_expression" / f"{wildcards.sample}"
    shell:
        "rsem-calculate-expression --hisat2-hca --num-threads {threads} --paired-end {input.reads} {params.base_ref} {params.base_out} > {log} 2>&1"

rule prepare_ref:
    input:
        result_dir / "denovo" / "trinity_reference" / "Trinity.fasta"
    output:
        result_dir / "denovo" / "denovo_index" / "trinity.seq"
    log:
        result_dir / "logs" / "denovo" / "rsem" / "rsem_prepare.log"
    conda:
        Path("..") / "envs" / "rsem_env.yaml"
    threads:
        10
    params:
        base_name = result_dir / "denovo" / "denovo_index" / "trinity"
    shell:
        "rsem-prepare-reference --hisat2-hca --num-threads {threads} {input} {params.base_name} > {log} 2>&1"