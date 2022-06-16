rule final_reference:
    input:
        calls = result_dir / "bctftools" / "mpileup" / "{sample}" / "calls.vcf.gz",
        ref = result_dir / "best_references" / "{sample}_ref.fasta"
    output:
        result_dir / "final_references" / "{sample}_ref.fasta"
    log:
        idx = result_dir / "log" / "bcftools" / "index" / "{sample}.log",
        cns = result_dir / "log" / "bcftools" / "consensus" / "{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bcftools index {input.calls} --threads {threads} 2> {log.idx} && cat {input.ref} | bcftools consensus {input.calls} > {output} 2> {log.cns}"

rule contig_mpileup:
    input:
        bam = result_dir / "contig_map_sorted" / "{sample}_sorted.bam",
        idx = result_dir / "contig_map_sorted" / "{sample}_sorted.bam.bai",
        ref = result_dir / "best_references" / "{sample}_ref.fasta"
    output:
        result_dir / "bctftools" / "mpileup" / "{sample}" / "calls.vcf.gz"
    log:
        mpile = result_dir / "log" / "contig_mpileup" / "{sample}_mpileup.log",
        call = result_dir / "log" / "contig_mpileup" / "{sample}_call.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bcftools mpileup -Ou -f {input.ref} {input.bam} --threads {threads} 2> {log.mpile} | bcftools call --threads {threads} -mv -Oz -o {output} 2> {log.call}"


rule contig_idx:
    input:
        result_dir / "contig_map_sorted" / "{sample}_sorted.bam"
    output:
        result_dir / "contig_map_sorted" / "{sample}_sorted.bam.bai"
    log:
        result_dir / "log" / "contig_idx" / "{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools index {input} -@ {threads} 2> {log}"


rule contig_sort:
    input:
        result_dir / "contig_mapping" / "{sample}.bam"
    output:
        result_dir / "contig_map_sorted" / "{sample}_sorted.bam"
    log:
        result_dir / "log" / "contig_sort" / "{sample}.log"
    threads:
        10
    conda:
       Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools sort -o {output} {input} -@ {threads} 2> {log}"


rule contig_mapping:
    input:
        ref = result_dir / "best_references" / "{sample}_ref.fasta",
        contigs = result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta"
    output:
        result_dir / "contig_mapping" / "{sample}.bam"
    log:
        samtools = result_dir / "log" / "contig_mapping" / "{sample}_sam.log",
        minimap = result_dir / "log" / "contig_mapping" / "{sample}_map.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    shell:
        "minimap2 -t {threads} -a {input.ref} {input.contigs} 2> {log.minimap} | samtools view -h  --output-fmt BAM - > {output} 2> {log.samtools}"
    
rule best_reference:
    input:
        scores = result_dir / "blast" / "{sample}.csv",
        refs = config["reference"]
    output:
        result_dir / "best_references" / "{sample}_ref.fasta"
    conda:
        Path("..") / "envs" / "python_env.yaml"
    script:
        str(Path("..") / "scripts" / "best_reference.py")

rule blastn:
    input:
        query = result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta",
        blastdb = multiext(str(result_dir / "potential_refs" / "blast_db" / "reference"),
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    output:
        result_dir / "blast" / "{sample}.csv"
    log:
        result_dir / "log" / "blast" / "{sample}.log"
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    threads:
        10
    params:
        format = "10 qseqid sseqid evalue",
        extra=""
    wrapper:
        "v1.5.0/bio/blast/blastn"

rule blast_db:
    input:
        config["reference"]
    output:
        multiext(str(result_dir / "potential_refs" / "blast_db" / "reference"),
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    log:
        result_dir / "log" / "blast_db" / "reference.log"
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    threads:
        10
    params:
        base_name = lambda wildcards: result_dir / "potential_refs" / "blast_db" / "reference"
    shell:
        "makeblastdb -in {input} -dbtype nucl -parse_seqids -out {params.base_name} > {log} 2>&1"

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
        "spades.py -t {threads} -1 {input.nfq1} -2 {input.nfq2} -o {params.out_dir} > {log} 2>&1"    

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