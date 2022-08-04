# DeNovo genome assembly

rule assembly:
    input:
        Illumina = fq_assembly_input,
        ONT = ONT_assembly_input
    output:
        result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta"
    log:
        result_dir / "log" / "spades_denovo_assembly" / "{sample}.log"
    conda:
        Path("..") / "envs" / "assembly_env.yaml"
    threads:
        config['max_threads']
    params:
        out_dir = lambda wildcards: result_dir / "denovo_assembly" / f"{wildcards.sample}",
        ONT_flag = "--nanopore " if ONT_correction else ""
    shell:
        "spades.py {params.ONT_flag}{input.ONT} -t {threads} -1 {input.Illumina[0]} -2 {input.Illumina[1]} -o {params.out_dir} > {log} 2>&1"

# Quality Control of the assembly

rule assembly_qc:
    input:
        expand(result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta", sample=list(samples.index))
    output:
        result_dir / "qc_files" / "assembly" / "report.html"
    log:
        result_dir / "log" / "assembly_qc" / "quast.log"
    conda:
        Path("..") / "envs" / "assembly_env.yaml"
    threads:
        min(10, config['max_threads'])
    params:
        out_dir = result_dir / "qc_files" / "assembly"
    shell:
        "quast -o {params.out_dir} --threads {threads} --conserved-genes-finding --k-mer-stats {input} > {log} 2>&1"

# annotation of the assembled genome

rule annotation:
    input:
        result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta"
    output:
        result_dir / "annotation" / "{sample}.gff"
    log:
        result_dir / "log" / "prokka_annotation" / "{sample}.log"
    conda:
        Path("..") / "envs" / "assembly_env.yaml"
    threads:
        config['max_threads']
    params:
        out_dir = result_dir / "annotation",
        prefix = lambda wildcards: {wildcards.sample}
    shell:
        "prokka --force --outdir {params.out_dir} --prefix {params.prefix} {input} > {log} 2>&1"