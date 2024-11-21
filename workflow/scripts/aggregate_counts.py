import os
import pandas as pd

''' 
Join count files to create count table

sgRNA names are expected to be in the format:
GENE_sgGENE_1, etc
'''

# Prepare data frame with gene and sgRNA names
fasta = pd.read_csv(snakemake.input["fasta"], header=None)
df = fasta[fasta[0].str.contains(">")].copy()
df.loc[:,0] = fasta.loc[:,0].str.replace(">", "", regex=False)
df.columns = ["sgRNA"]
df.loc[:,"gene"] = df.loc[:,"sgRNA"].str.split(pat="_", n=1, expand=True)[0]

# Merge all count files to df to create count table
for file in snakemake.input["files"]:
    # Read count file
    tmp = pd.read_csv(file, sep=" ", header=None)
    
    # Rename headers to wildcard value (for counts) and sgRNA 
    wildcard_value = os.path.basename(file).replace(".guidecounts.txt", "")
    tmp.columns = [wildcard_value, "sgRNA"]
    
    # Merge to df
    df = pd.merge(df, tmp, on="sgRNA", how="left")

# Replace missing values with zero and round to whole numbers
df = df.fillna(0)
numeric_cols = df.select_dtypes(include="number").columns
df[numeric_cols] = df[numeric_cols].round(decimals=0).astype(int)

# Remove duplicate rows (observed in some cases)
df = df.drop_duplicates()

# Save data frame to file
df.to_csv(snakemake.output[0], sep='\t', index=False)
