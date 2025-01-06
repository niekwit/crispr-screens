import os
import pandas as pd

''' 
Join count files to create count table
'''

# Prepare data frame with gene and sgRNA names
gene_column = snakemake.config["csv"]["gene_column"]
name_column = snakemake.config["csv"]["name_column"]
csv = pd.read_csv(snakemake.input["csv"])
df = csv.iloc[:, [name_column, gene_column]]
df.columns = ["sgRNA", "gene"]

# Merge all count files to df to create count table
for file in snakemake.input["files"]:
    # Read count file
    try:
        tmp = pd.read_csv(file, sep=" ", header=None)
    except pd.errors.EmptyDataError:
        # Create empty data frame if file is empty
        tmp = pd.DataFrame(columns=[0, 1])
    
    # Rename headers to wildcard value (for counts) and sgRNA 
    wildcard_value = os.path.basename(file).replace(".guidecounts.txt", "")
    tmp.columns = [wildcard_value, "sgRNA"]
    
    # Merge to df
    df = pd.merge(df, tmp, on="sgRNA", how="left")

# Check if columns other than 'gene' and 'sgRNA' contain only NAs
count_columns = df.columns.difference(["gene", "sgRNA"])
na_columns = df[count_columns].isna().all()

if na_columns.all():
    raise ValueError("All count columns are empty. Check trim/alignment settings.")

# Replace missing values with zero and round to whole numbers
df = df.fillna(0)
numeric_cols = df.select_dtypes(include="number").columns
df[numeric_cols] = df[numeric_cols].round(decimals=0).astype(int)

# Remove duplicate rows (observed for some libraries)
df = df.drop_duplicates()

# Save data frame to file
df.to_csv(snakemake.output[0], sep='\t', index=False)
