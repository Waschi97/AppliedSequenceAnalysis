rule mlst:
    input:
        assembly = result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta"
    output:
        mlst = result_dir / "db_screening" / "mlst" / "{sample}.tsv",
    log:
        result_dir / "log" / "mlst" / "{sample}.log",
    params:
        extra = lambda wildcards: f"--nopath --quiet --label {wildcards.sample}",
    threads: 
        config["max_threads"]
    wrapper:
        "v1.7.0/bio/mlst"

rule mlst_summary:
    input:
        expand(result_dir / "db_screening" / "mlst" / "{sample}.tsv", sample=list(samples.index))
    output:
        result_dir / "db_screening" / "mlst_summary.tsv"
    shell:
        "cat {input} > {output}"

rule abricate_card:
    input:
        expand(result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta", sample=list(samples.index))
    output:
        result_dir / "db_screening" / "abricate" / "card.tsv"
    log:
        result_dir / "log" / "abricate_card.log"
    params:
        db = "card"
    conda:
        Path("..") / "envs" / "abricate_env.yaml"
    shell:
        "abricate --db {params.db} {input} > {output} 2> {log}"


rule abricate_vfdb:
    input:
        expand(result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta", sample=list(samples.index))
    output:
        result_dir / "db_screening" / "abricate" / "vfdb.tsv"
    log:
        result_dir / "log" / "abricate_vfdb.log"
    params:
        db = "vfdb"
    conda:
        Path("..") / "envs" / "abricate_env.yaml"
    shell:
        "abricate --db {params.db} {input} > {output} 2> {log}"

rule abricate_plasmid:
    input:
        expand(result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta", sample=list(samples.index))
    output:
        result_dir / "db_screening" / "abricate" / "plasmidfinder.tsv"
    log:
        result_dir / "log" / "abricate_plasmidfinder.log"
    params:
        db = "plasmidfinder"
    conda:
        Path("..") / "envs" / "abricate_env.yaml"
    shell:
        "abricate --db {params.db} {input} > {output} 2> {log}"

rule abricate_summary:
    input:
        determine_db_screening_exe()
    output:
        result_dir / "db_screening" / "abricate_summary.tsv"
    log:
        result_dir / "log" / "abricate_summary.log"
    conda:
        Path("..") / "envs" / "abricate_env.yaml"
    shell:
        "abricate --summary {input} > {output} 2> {log}"

