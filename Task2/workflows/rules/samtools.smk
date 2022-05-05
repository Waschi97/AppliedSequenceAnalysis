rule stats:
    input:
        data = result_dir / "bam_sorted" / "{sample}.bam",
        idx = result_dir / "bam_sorted" / "{sample}.bam.bai"
    output:
        result_dir / "stats" / "{sample}.stats"
    threads:
        4
    conda:
        Path.cwd() / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools idxstats -@ {threads} {input.data} > {output}"

rule index:
    input:
        result_dir / "bam_sorted" / "{sample}.bam"
    output:
        result_dir / "bam_sorted" / "{sample}.bam.bai"
    threads:
        4
    conda:
        Path.cwd() / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools index -@ {threads} {input}"

rule sort:
    input:
        result_dir / "bam" / "{sample}.bam"
    output:
        result_dir / "bam_sorted" / "{sample}.bam"
    threads:
        4
    conda:
        Path.cwd() / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule convert:
    input:
        result_dir / "sam" / "{sample}.sam"
    output:
        result_dir / "bam" / "{sample}.bam"
    threads:
        4
    conda:
        Path.cwd() / "envs" / "samtools_bowtie_env.yaml"
    shell:
        "samtools view -@ {threads} -b -o {output} {input}"