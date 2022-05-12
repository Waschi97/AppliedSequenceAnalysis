rule multiqc:
    input:
        raw_html = expand(result_dir / "qc_files" / "qc_raw" / "{sample}_{i}.html", i=[1,2], sample=list(samples.index)),
        raw_zip = expand(result_dir / "qc_files" / "qc_raw" / "{sample}_{i}_fastqc.zip", i=[1,2], sample=list(samples.index)),
        trimm_html = expand(result_dir / "qc_files" / "qc_trimmed" / "{sample}_{j}{k}.html", k=['P','U'], j=[1,2], sample=list(samples.index)),
        trimm_zip = expand(result_dir / "qc_files" / "qc_trimmed" / "{sample}_{j}{k}_fastqc.zip", k=['P','U'], j=[1,2], sample=list(samples.index))
    output:
        result_dir / "qc_files" / "final_qc_report.html"
    params:
        dir = result_dir / "qc_files"
    conda:
        Path.cwd() / "envs" / "qc_tools_env.yaml"
    shell:
        "multiqc {params.dir} -n {output}"

rule fastqc_trimm:
    input:
        lambda wildcards: result_dir / "trimmed" / f"{wildcards.sample}_{wildcards.j}{wildcards.k}.fq.gz"
    output:
        html = result_dir / "qc_files" / "qc_trimmed" / "{sample}_{j}{k}.html",
        zip = result_dir / "qc_files" / "qc_trimmed" / "{sample}_{j}{k}_fastqc.zip"
    params: "--quiet"
    log:
        Path.cwd() / "logs" / "fastqc_raw" / "{sample}_{j}{k}.log"
    threads: 1
    wrapper:
        "v1.4.0/bio/fastqc"

rule trimm:
    input:
        fq1 = lambda wildcards: samples.at[wildcards.sample, "fq1"],
        fq2 = lambda wildcards: samples.at[wildcards.sample, "fq2"]
    output:
        p1 = result_dir / "trimmed" / "{sample}_1P.fq.gz",
        u1 = result_dir / "trimmed" / "{sample}_1U.fq.gz",
        p2 = result_dir / "trimmed" / "{sample}_2P.fq.gz",
        u2 = result_dir / "trimmed" / "{sample}_2U.fq.gz"
    params:
        output_dir = result_dir / "trimmed",
        base_name = lambda wildcards: result_dir / "trimmed" / f"{wildcards.sample}.fq.gz"
        # TODO: set actual parameters
    threads:
        4
    conda:
        Path.cwd() / "envs" / "qc_tools_env.yaml"
    shell:
        "mkdir -p {params.output_dir} && trimmomatic PE -threads {threads} {input.fq1} {input.fq2} -baseout {params.base_name} ILLUMINACLIP:{adapter_seqs}:2:30:10"

rule fastqc_raw:
    input:
        lambda wildcards: samples.at[wildcards.sample, "fq1"] if wildcards.i == '1' else samples.at[wildcards.sample, "fq2"]
    output:
        html = result_dir / "qc_files" / "qc_raw" / "{sample}_{i}.html",
        zip = result_dir / "qc_files" / "qc_raw" / "{sample}_{i}_fastqc.zip"
        # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        result_dir / "logs" / "fastqc_raw" / "{sample}_{i}.log"
    threads: 1
    wrapper:
        "v1.4.0/bio/fastqc"