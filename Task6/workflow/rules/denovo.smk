rule trinity:
    input:
        left = fq_denovo_input(True),
        right = fq_denovo_input(False)
    output:
        result_dir / "denovo" / "trinity_reference" / "Trinity.fasta"
    log:
        result_dir / "logs" / "denovo" / "trinity" / "trinity.log"
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