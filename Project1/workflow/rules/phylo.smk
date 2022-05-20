rule mult_align:
    input:
        result_dir / "assembled_genomes" / "all_genomes.fa"
    output:
        result_dir / "mult_align" / "alignment.fa"
    conda:
        Path("..") / "envs" / "phylo_env.yaml"
    params:
    shell:
        "clustalo -i {input} -o {output}"

rule merge_fasta:
    input:
        expand(result_dir / "assembled_genomes" / "{sample}.fa", sample=list(samples.index))
    output:
        result_dir / "assembled_genomes" / "all_genomes.fa"
    shell:
        "cat {input} > {output}"