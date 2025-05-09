import pandas as pd

# Snakemake variables
count_file = snakemake.input["counts"]
matrix_file = snakemake.input["matrix"]
output_file = snakemake.output[0]

# Read matrix and retrieve sample names (first column)
matrix = pd.read_csv(matrix_file, sep="\t", low_memory=False)
sample_names = matrix.iloc[:, 0].tolist()
to_keep = ["sgRNA", "gene"] + sample_names

# Read count file
counts = pd.read_csv(count_file, sep="\t", low_memory=False)

# Make sure that the order of the sample columns is the same as
# in the design matrix
counts = counts[to_keep]

# Save to file
counts.to_csv(output_file, sep="\t", index=False)
