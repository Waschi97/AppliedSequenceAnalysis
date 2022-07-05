rule transcripts:
    input:
        reads = result_dir / "ref_mapping" / "bam_sorted" / "{sample}.bam",
        gff = ref_annotation
    output:
        result_dir / "ref_mapping" / "transcripts" / "{sample}.gft"
    log:
        result_dir / "logs" / "ref_mapping" / "transcripts" / "stringtie_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "ref_assembly.yaml"
    shell:
        "stringtie -o {output} -G {input.gff} {input.reads}"

rule ref_sort:
    input:
        result_dir / "ref_mapping" / "mapped" / "{sample}.bam"
    output:
        result_dir / "ref_mapping" / "bam_sorted" / "{sample}.bam"
    log:
        result_dir / "logs" / "ref_mapping" / "sort" / "bam_sort_{sample}.log"
    threads:
        10
    conda:
        Path("..") / "envs" / "ref_assembly.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} > {log} 2>&1"

rule ref_align:
    input:
      reads = fq_assembly_input,
      idx_base = result_dir / "ref_mapping" / "ref_index"
    output:
      result_dir / "ref_mapping" / "mapped" / "{sample}.bam"
    log:
        result_dir / "logs" / "ref_mapping" / "hisat2_align_{sample}.log"
    params:
      extra = "",
      idx = result_dir / "ref_mapping" / "ref_index" / "genome",
    threads: 10
    wrapper:
      "v1.7.0/bio/hisat2/align"

rule ref_index:
    input:
        fasta = reference
    output:
        directory(result_dir / "ref_mapping" / "ref_index")
    params:
        prefix = result_dir / "ref_mapping" / "ref_index" / "genome"
    log:
        result_dir / "logs" / "ref_mapping" / "hisat2_index.log"
    threads: 10
    wrapper:
        "v1.7.0/bio/hisat2/index"