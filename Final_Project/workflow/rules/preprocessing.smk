# Trimming Illumina

rule trimm_short:
    input:
        fq1 = lambda wildcards: samples.at[wildcards.sample, "fq1"],
        fq2 = lambda wildcards: samples.at[wildcards.sample, "fq2"]
    output:
        p1 = result_dir / "trimmed" / "{sample}_1P.fq.gz",
        u1 = result_dir / "trimmed" / "{sample}_1U.fq.gz",
        p2 = result_dir / "trimmed" / "{sample}_2P.fq.gz",
        u2 = result_dir / "trimmed" / "{sample}_2U.fq.gz"
    log:
        result_dir / "log" / "trimmomatic" / "{sample}.log"
    params:
        output_dir = result_dir / "trimmed",
        base_name = lambda wildcards: result_dir / "trimmed" / f"{wildcards.sample}.fq.gz",
        mismatches = config["trimmomatic"]["mismatches"],
        score_threshold = config["trimmomatic"]["score_threshold"],
    threads:
        min(10, config['max_threads'])
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "mkdir -p {params.output_dir} && trimmomatic PE -threads {threads} {input.fq1} {input.fq2} -baseout {params.base_name} ILLUMINACLIP:{adapter_seqs}:{params.mismatches}:{params.score_threshold}:10 > {log} 2>&1"

# Trimming Oxford Nanopore

rule trimm_long:
    input:
        lambda wildcards: samples.at[wildcards.sample, "ONT"]
    output:
        result_dir / "trimmed_ONT" / "{sample}_ONT_trimmed.fastq.gz"
    log:
        result_dir / "log" / "porechop" / "{sample}.log"
    threads:
        min(10, config['max_threads'])
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "porechop -i {input} -o {output} --threads {threads} > {log} 2>&1"

# Contaminant Filtering

rule contaminant_indexing:
    input:
        config['contamination_db']
    output:
        idx = expand(result_dir / "contaminant_files" / "genome_idx" / "{{sample}}" / "genome{end}", end=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"])
    log:
        result_dir / "log" / "contaminant_indexing" / "{sample}_bowtie2_indexing.log"
    threads:
        config['max_threads']
    params:
        idx_base = lambda wildcards: result_dir / "contaminant_files" / "genome_idx" / f"{wildcards.sample}" / "genome"
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "bowtie2-build {input} {params.idx_base} -p {threads} > {log} 2>&1"

rule contaminant_mapping:
    input:
        fqs = fq_trimmed_input,
        idx = expand(result_dir / "contaminant_files" / "genome_idx" / "{{sample}}" / "genome{end}", end=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"])
    output:
        result_dir / "contaminant_files" / "sam" / "{sample}.sam"
    log:
        result_dir / "log" / "bowtie2_mapping" / "{sample}.log"
    threads:
        config['max_threads']
    params:
        idx_base = lambda wildcards: result_dir / "contaminant_files" / "genome_idx" / f"{wildcards.sample}" / "genome"
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "bowtie2 -x {params.idx_base} -1 {input.fqs[0]} -2 {input.fqs[1]} -p {threads} -S {output} > {log} 2>&1"

rule contaminant_convert:
    input:
        result_dir / "contaminant_files" / "sam" / "{sample}.sam"
    output:
        result_dir / "contaminant_files" / "bam" / "{sample}.bam"
    log:
        result_dir / "log" / "contaminant_convert" / "sam2bam_{sample}.log"
    threads:
        min(10, config['max_threads'])
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "samtools view -@ {threads} -b -o {output} {input} > {log} 2>&1"

rule filter:
    input:
        result_dir / "contaminant_files" / "bam" / "{sample}.bam"
    output:
        result_dir / "contaminant_files" / "bam_filtered" / "{sample}.bam"
    log:
        result_dir / "log" / "filter" / "filter_{sample}.log"
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "samtools view -b -f 12 -F 256 {input} > {output} > {log} 2>&1"

rule contaminant_sort:
    input:
        result_dir / "contaminant_files" / "bam" / "{sample}.bam"
    output:
        result_dir / "contaminant_files" / "bam_sorted" / "{sample}.bam"
    log:
        result_dir / "log" / "contaminant_sort" / "bam_sort_{sample}.log"
    threads:
        min(10, config['max_threads'])
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "samtools sort -n -m 5G -@ {threads} -o {output} {input} > {log} 2>&1"

rule filtered_fastq:
    input:
        result_dir / "contaminant_files" / "bam_sorted" / "{sample}.bam"
    output:
        fq1 = result_dir / "filtered_fastq" / "{sample}_1.fastq.gz",
        fq2 = result_dir / "filtered_fastq" / "{sample}_2.fastq.gz"
    log:
        result_dir / "log" / "filtered_fastq" / "{sample}.log"
    threads:
        min(10, config['max_threads'])
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "samtools fastq -@ {threads} {input} -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -n > {log} 2>&1"

# Coverage Depth Normalization

rule cov_normalization:
    input:
        fq_assembly_input
    output:
        nfq1 = result_dir / "normalization" / "{sample}_norm_1.fastq.gz",
        nfq2 = result_dir / "normalization" / "{sample}_norm_2.fastq.gz"
    log:
        result_dir / "log" / "normalization" / "{sample}.log"
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    threads:
        min(10, config['max_threads'])
    shell:
        "bbnorm.sh in={input[0]} in2={input[1]} out={output.nfq1} out2={output.nfq2} threads={threads}> {log} 2>&1"

# Quality Control Illumina

rule fastqc_raw:
    input:
        lambda wildcards: samples.at[wildcards.sample, "fq1"] if wildcards.i == '1' else samples.at[wildcards.sample, "fq2"]
    output:
        html = result_dir / "qc_files" / "preprocessing" / "qc_raw" / "{sample}_{i}.html",
        zip = result_dir / "qc_files" / "preprocessing" / "qc_raw" / "{sample}_{i}_fastqc.zip"
        # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        result_dir / "log" / "qc_raw" / "{sample}_{i}_fqc_raw.log"
    params: "--quiet"
    threads: min(5, config['max_threads'])
    wrapper:
        "v1.4.0/bio/fastqc"

rule fastqc_prep:
    input:
        fq_preproccessed
    output:
        html = result_dir / "qc_files" / "preprocessing" / "qc_prep" / "{sample}_{i}.html",
        zip = result_dir / "qc_files" / "preprocessing" / "qc_prep" / "{sample}_{i}_fastqc.zip"
        # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        result_dir / "log" / "qc_prep" / "{sample}_{i}_fqc_prep.log"
    params: "--quiet"
    threads: min(5, config['max_threads'])
    wrapper:
        "v1.4.0/bio/fastqc"

# Quality Control Oxford Nanopore

rule nanoplot_raw:
    input:
        lambda wildcards: samples.at[wildcards.sample, "ONT"]
    output:
        result_dir / "qc_files" / "preprocessing" / "nanoplot" / "{sample}_raw_NanoStats.txt"
    log:
        result_dir / "log" / "nanoplot"/ "{sample}_raw.log"
    params:
        out_dir = lambda wildcards: result_dir / "qc_files" / "preprocessing" / "nanoplot",
        prefix = lambda wildcards: f"{wildcards.sample}_raw_"
    threads:
        min(5, config['max_threads'])
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "NanoPlot -t {threads} -o {params.out_dir} -p {params.prefix} --fastq {input} 2> {log}"

rule nanoplot_prep:
    input:
        rules.trimm_long.output
    output:
        result_dir / "qc_files" / "preprocessing" / "nanoplot" / "{sample}_prep_NanoStats.txt"
    log:
        result_dir / "log" / "nanoplot"/ "{sample}_prep.log"
    params:
        out_dir = lambda wildcards: result_dir / "qc_files" / "preprocessing" / "nanoplot",
        prefix = lambda wildcards: f"{wildcards.sample}_prep_"
    threads:
        min(5, config['max_threads'])
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "NanoPlot -t {threads} -o {params.out_dir} -p {params.prefix} --fastq {input} 2> {log}"

# MultiQC

rule multiqc:
    input:
        fastqc_files = determine_fastqc_exe(),
        nanoplot_files = determine_nanoplot_exe()
    output:
        result_dir / "qc_files" / "preprocessing" / "final_qc_report.html"
    log:
        result_dir / "log" / "multiqc.log"
    params:
        dir = result_dir / "qc_files" / "preprocessing"
    conda:
        Path("..") / "envs" / "preprocessing_env.yaml"
    shell:
        "multiqc {params.dir} -n {output} > {log} 2>&1"