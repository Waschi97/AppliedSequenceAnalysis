rule core_genome:
    input:
        expand(result_dir / "annotation" / "{sample}.gff", sample = list(samples.index))
    output:
        result_dir / "core_genome" / "core_gene_alignment.aln"
    log:
        log = result_dir / "log" / "roary" / "core_genome.log",
        err = result_dir / "log" / "roary" / "core_genome.err"
    params:
        out_dir = result_dir / "core_genome",
        threshold = config["roary"]["threshold"]
    threads:
        config['max_threads']
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    shell:
        "rm -rf {params.out_dir}; roary -e -n -p {threads} -cd {params.threshold} -f {params.out_dir} -v {input} 2> {log.err} > {log.log}"


