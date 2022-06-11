rule tree_convert:
    input:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.txt"
    output:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.pdf"
    log:
        result_dir / "log" / "phylogeny" / "text2pdf.log"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    params:
        script = Path("workflow") / "scripts" / "txt2pdf.py"
    shell:
        "python3 {params.script} -o {output} {input} 2> {log}"

rule tree_topology:
    input:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.tree"
    output:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.txt"
    log:
        result_dir / "log" / "phylogeny" / "nw_utils.log"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    shell:
        "nw_topology {input} | nw_display - > {output} 2> {log}"

rule tree:
    input:
        result_dir / "mult_align" / "alignment.fa"
    output:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.tree"
    log:
        result_dir / "log" / "phylogeny" / "fasttree.log"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    shell:
        "fasttree -gtr -nt {input} > {output} 2> {log}"

rule mult_align:
    input:
        result_dir / "assembled_genomes" / "all_genomes.fa"
    output:
        result_dir / "mult_align" / "alignment.fa"
    log:
        result_dir / "log" / "phylogeny" / "mafft.log"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    threads:
        30
    shell:
        "mafft --thread {threads} {input} > {output} 2> {log}"

rule merge_fasta:
    input:
        expand(result_dir / "assembled_genomes" / "{sample}.fa", sample=list(samples.index))
    output:
        result_dir / "assembled_genomes" / "all_genomes.fa"
    shell:
        "cat {input} > {output}"