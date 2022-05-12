import pandas as pd 
import sys
from pathlib import Path

if __name__ == "__main__":
    stat_file = snakemake.input[0]
    out_file = snakemake.output[0]
    df = pd.read_csv(stat_file, header=None, sep='\t')
    rpk = df.iloc[:,2].div(df.iloc[:,1]) * 1000
    df.insert(4, None ,rpk)
    df.to_csv(out_file, header=None, index=False, sep='\t')