rule kraken:
    input:
        fq_trimmed_input
    output:
        report = result_dir / "kraken_screening" / "{sample}_kraken.report",
        tsv = result_dir / "kraken_screening" / "{sample}_kraken.tsv"
    log:
        result_dir / "log" / "screening" / "{sample}_kraken.log"
    params:
        db = config["kraken_db"]
    conda:
        Path("..") / "envs" / "tax_screen_env.yaml"
    shell:
        "kraken2 --db {params.db} --report {output.report} --paired {input[0]} {input[1]} > {output.tsv} 2> {log}"

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
        Path("..") / "envs" / "tax_screen_env.yaml"
    shell:
        "multiqc {params.dir} -n {output} > {log} 2>&1"