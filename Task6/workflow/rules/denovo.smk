rule kallisto_quant:
    input:
        fastq = fq_assembly_input,
        index = result_dir / "trinity_idx" / "Trinity.idx"
    output:
        directory(result_dir / "quant_results" / "{sample}")
    params:
        extra = ""
    log:
        result_dir / "logs" / "kallisto_quant" / "{sample}.log"
    threads: 10
    wrapper:
        "v1.3.2/bio/kallisto/quant"

rule kallisto_index:
    input:
        fasta = result_dir / "trinity_reference" / "Trinity.fasta"
    output:
        index = result_dir / "trinity_idx" / "Trinity.idx"
    params:
        extra = "--kmer-size=5"
    log:
        result_dir / "logs" / "kallisto_index.log"
    threads: 10
    wrapper:
        "v1.3.2/bio/kallisto/index"


rule trinity:
    input:
        left = fq_denovo_input(True),
        right = fq_denovo_input(False)
    output:
        result_dir / "trinity_reference" / "Trinity.fasta"
    log:
        result_dir / "logs" / "trinity" / "trinity.log"
    params:
        extra=""
    threads: 10
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_gb=20
    wrapper:
        "v1.7.0/bio/trinity"