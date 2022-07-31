rule genbank_conversion:
    input:
        seqs = result_dir / "denovo_assembly" / "{sample}" / "contigs.fasta",
        gff = result_dir / "annotation" / "{sample}.gff",
    output:
        result_dir / "genbank_genomes" / "{sample}.gbk"
    log:
        result_dir / "log" / "genbank_conversion" / "{sample}.log"
    conda:
        Path("..") / "envs" / "pan_genome_env.yaml"
    params:
        base_name = lambda wildcards: result_dir / "genbank_genomes" / f"{wildcards.sample}"
    shell:
        "seqret -sequence {input.seqs} -feature -fformat gff -fopenfile {input.gff} -osformat genbank -osname_outseq {params.base_name} -ofdirectory_outseq gbk_file -auto > {log} 2>&1"