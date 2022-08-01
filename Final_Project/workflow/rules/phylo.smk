# Prepare blacklist before core genome
sample_names = list(samples.index)

if config['blacklist'] != "" and config['blacklist_core_genome']:
    blacklist = list(pd.read_csv(Path(config["blacklist"]), header=None)[0])
    
    for name in blacklist:
        if name in sample_names:
            sample_names.remove(name)

# function for blacklisting only the tree (post alignment calculation)
def alignment_input():
    if config['blacklist_core_genome'] or config['blacklist'] == "":
        return result_dir / "core_genome" / "core_gene_alignment.aln"
    
    return result_dir / "core_genome_alignment_blacklisted" / "alignment.aln"

rule core_genome:
    input:
        expand(result_dir / "annotation" / "{sample}.gff", sample = sample_names)
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

rule blacklist_alignment:
    input:
        algn = result_dir / "core_genome" / "core_gene_alignment.aln",
        blk_lst = config['blacklist']
    output:
        result_dir / "core_genome_alignment_blacklisted" / "alignment.aln"
    conda:
        Path("..") / "envs" / "python_env.yaml"
    script:
        str(Path("..") / "scripts" / "blacklist_alignment.py")

rule tree_ML:
    input:
        alignment_input()
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

rule newick_convert:
    input:
        result_dir / "RAxML_files" / "RAxML_bestTree.tree"
    output:
        result_dir / "core_genome_phylogeny" / "RAxML_bestTree.png"
    conda:
        Path("..") / "envs" / "python_env.yaml"
    script:
        str(Path("..") / "scripts" / "draw_tree.py")