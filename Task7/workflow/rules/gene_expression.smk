rule star_index:
    input:
        fasta = reference
    output:
        directory(result_dir / "star_index"),
    threads: 
        10
    params:
        extra = f"--sjdbGTFfile {ref_annotation}",
    log:
        result_dir / "logs" / "star_index.log"
    wrapper:
        "v1.7.0/bio/star/index"

rule star_pe_multi:
    input:
        rules.star_index.output,
        fq1 = fq_denovo_input(True),
        fq2 = fq_denovo_input(False)
    output:
        sam = result_dir / "star" / "{sample}_aligned.out.sam",
        log = result_dir / "star" / "{sample}_log.out"
    log:
        result_dir/ "logs" / "star" / "{sample}.log"
    params:
        idx = str(result_dir / "star_index"),
        extra = "",
    threads: 
        10
    wrapper:
        "v1.7.0/bio/star/align"

rule feature_counts:
    input:
        sam = result_dir / "star" / "{sample}_aligned.out.sam",
        annotation = ref_annotation
    output:
        expand(result_dir / "featurecounts" / "{{sample}}.{end}", end=[".featureCounts", ".featureCounts.summary", ".featureCounts.jcounts"])
    threads:
        10
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -p"
    log:
        result_dir / "logs" / "featurecounts" / "{sample}.log"
    wrapper:
        "0.72.0/bio/subread/featurecounts"
        