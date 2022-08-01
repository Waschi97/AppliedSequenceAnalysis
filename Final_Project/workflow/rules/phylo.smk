rule core_genome:
    input:
        expand(result_dir / "annotation" / "{sample}.gff", sample = list(samples.index))
    output:
        result_dir / "core_genome" / "core_gene_alignment.aln"
    log:
        log = result_dir / "log" / "roary" / "core_genome.log",
        err = result_dir / "log" / "roary" / "core_genome.err"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    threads:
        config['max_threads']
    params:
        out_dir = result_dir / "core_genome",
        threshold = config["roary"]["threshold"]
    shell:
        "rm -rf {params.out_dir}; roary -e -n -p {threads} -cd {params.threshold} -f {params.out_dir} -v {input} 2> {log.err} > {log.log}"

rule tree_ML:
    input:
        result_dir / "core_genome" / "core_gene_alignment.aln"
    output:
        result_dir / "RAxML_files" / "RAxML_bestTree.tree"
    log:
        result_dir / "log" / "raxml.log"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    threads:
        config["max_threads"]
    params:
        out_dir = result_dir / "RAxML_files",
        model = raxml_model
    shell:
        "raxmlHPC -s {input} -n tree -m {params.model} -w {params.out_dir} > {log} 2>&1"