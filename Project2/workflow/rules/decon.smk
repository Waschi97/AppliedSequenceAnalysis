rule screen:
    input:
        expand(result_dir / "kraken_screening" / "{sample}_kraken.report", sample=list(samples.index))
    output:
        result_dir / "kraken_screening" / "final_kraken_report.html"
    log:
        result_dir / "log" / "screening" / "multiqc_kraken.log"
    params:
        dir = result_dir / "kraken_screening"
    conda:
        Path("..") / "envs" / "decon_env.yaml"
    shell:
        "multiqc {params.dir} -n {output} 2> {log}"

rule kraken:
    input:
        fq_trimmed_input
    output:
        report = result_dir / "kraken_screening" / "{sample}_kraken.report",
        tsv = result_dir / "kraken_screening" / "{sample}_kraken.tsv"
    log:
        result_dir / "log" / "screening" / "{sample}_kraken.log"
    params:
        db = config["contamination_db"]
    conda:
        Path("..") / "envs" / "decon_env.yaml"
    shell:
        "kraken2 --db {params.db} --report {output.report} --paired {input[0]} {input[1]} > {output.tsv} 2> {log}"

rule filtered_fastq:
    input:
        result_dir / "contaminant_files" / "bam_sorted" / "{sample}.bam"
    output:
        fq1 = result_dir / "filtered_fastq" / "{sample}_1.fastq.gz",
        fq2 = result_dir / "filtered_fastq" / "{sample}_2.fastq.gz"
    log:
        result_dir / "log" / "filtered_fastq" / "{sample}.log"
    threads:
        10
    shell:
        "samtools fastq -@ {threads} {input} -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -n 2> {log}"

rule contaminant_sort:
    input:
        result_dir / "contaminant_files" / "bam" / "{sample}.bam"
    output:
        result_dir / "contaminant_files" / "bam_sorted" / "{sample}.bam"
    log:
        result_dir / "log" / "contaminant_sort" / "bam_sort_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools sort -n -m 5G -@ {threads} -o {output} {input} 2> {log}"

rule filter:
    input:
        result_dir / "contaminant_files" / "bam" / "{sample}.bam"
    output:
        result_dir / "contaminant_files" / "bam_filtered" / "{sample}.bam"
    log:
        result_dir / "log" / "filter" / "filter_{sample}.log"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools view -b -f 12 -F 256 {input} > {output} 2> {log}"

rule contaminant_convert:
    input:
        result_dir / "contaminant_files" / "sam" / "{sample}.sam"
    output:
        result_dir / "contaminant_files" / "bam" / "{sample}.bam"
    log:
        result_dir / "log" / "contaminant_convert" / "sam2bam_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools view -@ {threads} -b -o {output} {input} 2> {log}"

rule contaminant_mapping:
    input:
        fqs = fq_trimmed_input,
        idx = expand(result_dir / "contaminant_files" / "genome_idx" / "{{sample}}" / "genome{end}", end=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"])
    output:
        result_dir / "contaminant_files" / "sam" / "{sample}.sam"
    log:
        result_dir / "log" / "bowtie2_mapping" / "{sample}.log"
    threads:
        10
    params:
        min_frag_len = config["bowtie2_mapping"]['I'],
        max_frag_len = config["bowtie2_mapping"]['X'],
        idx_base = lambda wildcards: result_dir / "contaminant_files" / "genome_idx" / f"{wildcards.sample}" / "genome"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2 -x {params.idx_base} -1 {input.fqs[0]} -2 {input.fqs[1]} -p {threads} -X {params.max_frag_len} -I {params.min_frag_len} -S {output} 2> {log}"

rule contaminant_indexing:
    input:
        get_contaminant_db
    output:
        idx = expand(result_dir / "contaminant_files" / "genome_idx" / "{{sample}}" / "genome{end}", end=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"])
    log:
        result_dir / "log" / "contaminant_indexing" / "{sample}_bowtie2_indexing.log"
    threads:
        10
    params:
        idx_base = lambda wildcards: result_dir / "contaminant_files" / "genome_idx" / f"{wildcards.sample}" / "genome"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "bowtie2-build {input} {params.idx_base} -p {threads} 2> {log}"
