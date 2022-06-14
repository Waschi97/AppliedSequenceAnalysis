rule best_reference:
    input:
        scores = expand(result_dir / "blast" / "{sample}.tsv", sample=list(samples.index)),
        refs = ref_files
    output:
        result_dir / "best_reference" / "reference.fasta"
    log:
        result_dir / "log" / "best_reference.log"
    conda:
        Path("..") / "envs" / "python_env.yaml"
    script:
        str(Path("..") / "scripts" / "best_reference.py")

rule blastn:
    input:
        seqs = result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta",
        dbs = expand(result_dir / "potential_refs" / "blast_db" / "{ref_id}.ndb", ref_id=ref_ids),
    output:
        result_dir / "blast" / "{sample}.tsv"
    log:
        result_dir / "log" / "blast" / "{sample}.log"
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    threads:
        10
    params:
        base_names = expand(result_dir / "potential_refs" / "blast_db" / "{ref_id}", ref_id=ref_ids)
    shell:
        "blastn -num_threads {threads} -db {params.base_names} -query {input.seqs} -outfmt '6' -out {output}"

rule blast_db:
    input:
        result_dir / "potential_refs" / "{ref_id}.fasta"
    output:
        result_dir / "potential_refs" / "blast_db" / "{ref_id}.ndb"
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    threads:
        10
    params:
        base_name = lambda wildcards: result_dir / "potential_refs" / "blast_db" / f"{wildcards.ref_id}"
    shell:
        "makeblastdb -in {input} -dbtype nucl -parse_seqids -out {params.base_name}"

rule denovo_assembly:
    input:
        nfq1 = result_dir / "normalization" / "{sample}_norm_1.fastq.gz",
        nfq2 = result_dir / "normalization" / "{sample}_norm_2.fastq.gz"
    output:
        result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta",
        result_dir / "denovo_assembly" / "{sample}" / "scaffolds.fasta"
    log:
        result_dir / "log" / "denovo_assembly" / "{sample}.log"
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    threads:
        10
    params:
        out_dir = lambda wildcards: result_dir / "denovo_assembly" / f"{wildcards.sample}"
    shell:
        "spades.py --careful -t {threads} -1 {input.nfq1} -2 {input.nfq2} -o {params.out_dir} > {log} 2>&1"    

rule cov_normalization:
    input:
        fq_assembly_input
    output:
        nfq1 = result_dir / "normalization" / "{sample}_norm_1.fastq.gz",
        nfq2 = result_dir / "normalization" / "{sample}_norm_2.fastq.gz"
    log:
        result_dir / "log" / "normalization" / "{sample}.log"
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    threads:
        10
    shell:
        "bbnorm.sh in={input[0]} in2={input[1]} out={output.nfq1} out2={output.nfq2} threads={threads}> {log} 2>&1"