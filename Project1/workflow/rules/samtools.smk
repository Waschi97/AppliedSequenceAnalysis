rule stats:
    input:
        data = result_dir / "bam_sorted" / "{sample}.bam",
        idx = result_dir / "bam_sorted" / "{sample}.bam.bai"
    output:
        result_dir / "stats" / "{sample}.stats"
    threads:
        10
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools idxstats -@ {threads} {input.data} > {output}"

rule index:
    input:
        result_dir / "bam_sorted" / "{sample}.bam"
    output:
        result_dir / "bam_sorted" / "{sample}.bam.bai"
    threads:
        30
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools index -@ {threads} {input}"

rule assembly:
    input:
        reads = result_dir / "bam_sorted" / "{sample}.bam"
    output:
        result_dir / "assembled_genomes" / "{sample}.fa"
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    params:
        out_prefix = lambda wildcards: result_dir / "assembled_genomes" / f"{wildcards.sample}",
        min_depth = config["ivar_assembly"]["min_depth"]
    shell:
        "samtools mpileup -aa -A -d 0 -Q 0 {input.reads} | ivar consensus -m {params.min_depth} -p {params.out_prefix}"

rule sort:
    input:
        result_dir / "bam" / "{sample}.bam"
    output:
        result_dir / "bam_sorted" / "{sample}.bam"
    threads:
        30
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule convert:
    input:
        result_dir / "sam" / "{sample}.sam"
    output:
        result_dir / "bam" / "{sample}.bam"
    threads:
        30
    conda:
        Path("..") / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools view -@ {threads} -b -o {output} {input}"