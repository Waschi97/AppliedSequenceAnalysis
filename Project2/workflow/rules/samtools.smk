rule stats:
    input:
        data = result_dir / "bam_sorted" / "{sample}.bam",
        idx = result_dir / "bam_sorted" / "{sample}.bam.bai"
    output:
        result_dir / "stats" / "{sample}.stats"
    log:
        result_dir / "log" / "stats" / "statistics_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools idxstats -@ {threads} {input.data} > {output} 2> {log}"

rule index:
    input:
        result_dir / "bam_sorted" / "{sample}.bam"
    output:
        result_dir / "bam_sorted" / "{sample}.bam.bai"
    log:
        result_dir / "log" / "index" / "indexing_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"

rule assembly:
    input:
        reads = result_dir / "bam_sorted" / "{sample}.bam"
    output:
        result_dir / "assembled_genomes" / "{sample}.fa"
    log:
        mpileup = result_dir / "log" / "ref_assembly" / "mpileup_{sample}.log",
        ivar = result_dir / "log" / "ref_assembly" / "ivar_{sample}.log"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    params:
        out_prefix = lambda wildcards: result_dir / "assembled_genomes" / f"{wildcards.sample}",
        min_depth = config["ivar_assembly"]["min_depth"]
    shell:
        "samtools mpileup -aa -A -d 0 -Q 0 {input.reads} 2> {log.mpileup} | ivar consensus -m {params.min_depth} -p {params.out_prefix} 2> {log.ivar}"

rule sort:
    input:
        result_dir / "bam" / "{sample}.bam"
    output:
        result_dir / "bam_sorted" / "{sample}.bam"
    log:
        result_dir / "log" / "sort" / "bam_sort_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"

rule convert:
    input:
        result_dir / "sam" / "{sample}.sam"
    output:
        result_dir / "bam" / "{sample}.bam"
    log:
        result_dir / "log" / "convert" / "sam2bam_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools view -@ {threads} -b -o {output} {input} 2> {log}"