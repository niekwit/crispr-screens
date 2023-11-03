import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import general_functions as utils

mpl.use('agg')

# write stdout/stderr to log file
sys.stderr = open(snakemake.log[0], "w")
    
# get number of sgRNAs in CRISPR library
fasta = pd.read_table(snakemake.params["fasta"], header = None)
lib_size = len(fasta) / 2

# extract number of single mapped aligned reads from counts-aggregated.tsv
df = pd.read_table(snakemake.input[0], sep = "\t")
column_names = list(df.columns)
del column_names[0:2] # delete sgRNA and gene columns

counts = {}
for i in column_names:
    count_sum = []
    count_sum.append(df[i].sum())
    counts[i] = count_sum

# convert counts to sequence coverage
df = pd.DataFrame(counts)
df = df / lib_size

# order columns alphabetically
df = df.reindex(sorted(df.columns), axis=1)

# transpose data frame
df = df.transpose()
df["samples"] = df.index
names = ["coverage", "samples"]
df.columns = names
df = df[["samples","coverage"]]

# plot coverage per sample
utils.plot(df,"Fold sequence coverage per sample",snakemake.output[0])


