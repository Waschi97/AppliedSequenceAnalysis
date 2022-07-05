rule kallisto_quant:
    input:
        fastq = fq_assembly_input,
        index = result_dir / "denovo" / "trinity_idx" / "Trinity.idx"
    output:
        directory(result_dir / "denovo" / "quant_results" / "{sample}")
    params:
        extra = ""
    log:
        result_dir / "logs" / "denovo" / "kallisto_quant" / "{sample}.log"
    threads: 10
    wrapper:
        "v1.3.2/bio/kallisto/quant"

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