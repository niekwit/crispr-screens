import os
import pandas as pd

''' 
Join count files to create count table
'''

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

counts = {}
    
for file in snakemake.input:
    # Check if file size is larger than 0
    assert os.path.getsize(file) > 0, f"{file} is empty (alignment failed?)"
    
    # Add counts to counts dict
    df = pd.read_csv(file,sep=" ",header = None)
    key = os.path.basename(file).replace(".guidecounts.txt", "")
    df.columns = [key, "sgRNA"]
    df[key].astype(int)
    counts.update({key:df})

df = pd.read_csv(snakemake.params[0], header = None)
df.columns = ["sgRNA"]

df = df[df["sgRNA"].str.contains(">")]
df = df.reset_index(drop = True)
df["sgRNA"] = df["sgRNA"].str.replace(">", "")
df["gene"] = df["sgRNA"].str.split(pat = "_",n = 1,expand = True)[0]

# Perform left join on count files
for key, value in counts.items():
    df = pd.merge(df, value, on='sgRNA', how='left')
df["sgRNA"] = df["sgRNA"].str.split(pat = "_",n = 1,expand = True)[1]

# Replace nan with zero
df = df.fillna(0)
df = df.sort_values(by = ["sgRNA"])

# Convert floats back to int after pandas merge (bug in pandas)
index_range = range(2, len(df.columns))
index_list = []
for i in index_range:
    index_list.append(i)
df[df.columns[index_list]] = df[df.columns[index_list]].astype(int)

# Save data frame to file
df.to_csv(snakemake.output[0], sep = '\t', index = False)
