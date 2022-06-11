rule multiqc:
    input:
        fastqc_files = determine_fastqc_exe(),
        mapping_html = determine_qualimap_exe()
    output:
        result_dir / "qc_files" / "final_qc_report.html"
    log:
        result_dir / "log" / "multiqc.log"
    params:
        dir = result_dir / "qc_files"
    conda:
        Path("..") / "envs" / "qc_tools_env.yaml"
    shell:
        "multiqc {params.dir} -n {output} > {log} 2>&1"

rule fastqc_trimm:
    input:
        lambda wildcards: result_dir / "trimmed" / f"{wildcards.sample}_{wildcards.j}P.fq.gz"
    output:
        html = result_dir / "qc_files" / "qc_trimmed" / "{sample}_{j}P.html",
        zip = result_dir / "qc_files" / "qc_trimmed" / "{sample}_{j}P_fastqc.zip"
    log:
        result_dir / "log" / "qc_trimmed" / "{sample}_{j}_fqc_trimm.log"
    params: "--quiet"
    threads: 1
    wrapper:
        "v1.4.0/bio/fastqc"

rule mapping_stats:
    input:
        result_dir / "bam_sorted" / "{sample}.bam"
    output:
        result_dir / "qc_files" / "qualimap" / "{sample}" / "qualimapReport.html"
    log:
        result_dir / "log" / "qualimap" / "{sample}.log"
    params:
        dir = lambda wildcards: result_dir / "qc_files" / "qualimap" / f"{wildcards.sample}"
    conda:
        Path("..") / "envs" / "qc_tools_env.yaml"
    shell:
        "qualimap bamqc -bam {input} -outdir {params.dir} > {log} 2>&1"

rule trimm:
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
        10
    conda:
        Path("..") / "envs" / "qc_tools_env.yaml"
    shell:
        "mkdir -p {params.output_dir} && trimmomatic PE -threads {threads} {input.fq1} {input.fq2} -baseout {params.base_name} ILLUMINACLIP:{adapter_seqs}:{params.mismatches}:{params.score_threshold}:10 > {log} 2>&1"

rule fastqc_raw:
    input:
        lambda wildcards: samples.at[wildcards.sample, "fq1"] if wildcards.i == '1' else samples.at[wildcards.sample, "fq2"]
    output:
        html = result_dir / "qc_files" / "qc_raw" / "{sample}_{i}.html",
        zip = result_dir / "qc_files" / "qc_raw" / "{sample}_{i}_fastqc.zip"
        # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        result_dir / "log" / "qc_raw" / "{sample}_{i}_fqc_raw.log"
    params: "--quiet"
    threads: 1
    wrapper:
        "v1.4.0/bio/fastqc"