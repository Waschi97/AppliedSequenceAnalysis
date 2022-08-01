from Bio import SeqIO
import pandas as pd

blacklist = list(pd.read_csv(str(snakemake.input.blk_lst), header=None)[0])

with open(str(snakemake.output), "w") as outfile:
    for record in SeqIO.parse(str(snakemake.input.algn), "fasta"):
        if str(record.id) in blacklist:
            continue
        SeqIO.write(record, outfile, "fasta")
