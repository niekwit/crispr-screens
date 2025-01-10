import pandas as pd

"""
Convert csv file to fasta file. Barcode name and sequence 
column numbers must be specified in config.yml.
"""

sgrna_column = snakemake.config["csv"]["name_column"]
sequence_column = snakemake.config["csv"]["sequence_column"]

# Load csv file
df = pd.read_csv(snakemake.input["csv"])

# Only keep sgRNA ID and sequence columns
df = df.iloc[:, [sgrna_column, sequence_column]]

# Check if sequence column is a DNA sequence
if not df.iloc[:, 1].str.match(r"^[ACGTacgt]+$").any():
    raise ValueError("Sequence column must contain valid DNA sequences.")

# Prepend ">" to the sgRNA ID
df.iloc[:, 0] = ">" + df.iloc[:, 0]

# Write fasta file
df.to_csv(snakemake.output["fasta"],
          sep="\n",
          index=False,
          header=False)