import pandas as pd
import os

quant_files = snakemake.input
output_file = snakemake.output[0]

dfs = []
for qf in quant_files:
    # Get patient name from the directory structure
    sample_name = qf.split("/")[-2] 
    
    # Read Salmon output
    df = pd.read_csv(qf, sep='\t')
    
    # Keep only Gene ID and Transcripts Per Million (TPM)
    df = df[['Name', 'TPM']].rename(columns={'TPM': sample_name})
    df.set_index('Name', inplace=True)
    dfs.append(df)

# Concatenate all patients into a single matrix and transpose
expr_matrix = pd.concat(dfs, axis=1).T
expr_matrix.index.name = "sample"

expr_matrix.to_csv(output_file)