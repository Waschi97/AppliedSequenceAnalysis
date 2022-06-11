rule cov_normalization:
    input:
        fq_assembly_input
    output:
        nfq1 = result_dir / "normalization" / "{sample}_norm_1.fastq.gz",
        nfq2 = result_dir / "normalization" / "{sample}_norm_2.fastq.gz"
    log:
        result_dir / "log" / "normalization" / "{sample}_norm.log"
    conda:
        Path("..") / "envs" / "denovo_env.yaml"
    threads:
        10
    shell:
        "bbnorm.sh in={input[0]} in2={input[1]} out={output.nfq1} out2={output.nfq2} threads={threads}> {log} 2>&1"