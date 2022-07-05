rule kallisto_index:
    input:
        fasta = result_dir / "denovo" / "trinity_reference" / "Trinity.fasta"
    output:
        index = result_dir / "denovo" / "trinity_idx" / "Trinity.idx"
    params:
        extra = "--kmer-size=5"
    log:
        result_dir / "logs" / "denovo" / "kallisto_index.log"
    threads: 10
    wrapper:
        "v1.3.2/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq = fq_assembly_input,
        index = result_dir / "denovo" / "trinity_idx" / "Trinity.idx"
    output:
        directory(result_dir / "denovo" / "quant_results" / "{sample}")
    params:
        extra = "-b 3"
    log:
        result_dir / "logs" / "denovo" / "kallisto_quant" / "{sample}.log"
    threads: 10
    wrapper:
        "v1.3.2/bio/kallisto/quant"

rule sleuth:
    input:
        sample_files,
        expand(rules.kallisto_quant.output, sample=list(samples.index))
    output:
        result_dir / "sleuth" / "results.csv",
        result_dir / "sleuth" / "pca.pdf"
    log:
        result_dir / "logs" / "sleuth.log"
    params:
        str(result_dir / "denovo" / "quant_results")
    conda:
        Path("..") / "envs" / "dge.yaml"
    script:
        str(Path("..") / "scripts" / "sleuth.R")