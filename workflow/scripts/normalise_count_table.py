import pandas as pd

df = pd.read_table(snakemake.input[0])
column_range = range(2,len(df.columns))
for i in column_range:
    column_sum = df.iloc[:,i].sum()
    df.iloc[:,i] = df.iloc[:,i] / column_sum * 1E8
    df.iloc[:,i] = df.iloc[:,i].astype(int)
df.to_csv(snakemake.output[0],
        index = False,
        header = True)