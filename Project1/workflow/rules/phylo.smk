rule tree_convert:
    input:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.txt"
    output:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.pdf"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    params:
        script = Path("workflow") / "scripts" / "txt2pdf.py"
    shell:
        "python3 {params.script} -o {output} {input}"

rule tree_topology:
    input:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.tree"
    output:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.txt"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    shell:
        "nw_topology {input} | nw_display - > {output}"

rule tree:
    input:
        result_dir / "mult_align" / "alignment.fa"
    output:
        result_dir / "phylogenetic_tree" / "cov_phylogeny.tree"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    shell:
        "fasttree -gtr -nt {input} > {output}"

rule mult_align:
    input:
        result_dir / "assembled_genomes" / "all_genomes.fa"
    output:
        result_dir / "mult_align" / "alignment.fa"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    threads:
        30
    shell:
        "mafft --thread {threads} {input} > {output}"

rule merge_fasta:
    input:
        expand(result_dir / "assembled_genomes" / "{sample}.fa", sample=list(samples.index))
    output:
        result_dir / "assembled_genomes" / "all_genomes.fa"
    shell:
        "cat {input} > {output}"