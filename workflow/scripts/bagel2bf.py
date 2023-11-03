import pandas as pd
from snakemake.shell import shell

b2dir = snakemake.input["b2dir"]
species = snakemake.params["species"]
fc = snakemake.input["fc"]
bf = snakemake.output["bf"]
comparison = snakemake.wildcards[0]
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

#load gene sets
if species == "human":
    eg = f"{b2dir}/CEGv2.txt" # essential genes
    neg = f"{b2dir}/NEGv1.txt" # non-essential genes

elif species == "mouse":
    eg = f"{b2dir}/CEG_mouse.txt" # essential genes
    neg = f"{b2dir}/NEG_mouse.txt" # non-essential genes
    
    
# create dictionary to store the sample column numbers in foldchange file
# (samples cannot be refered to by column name)
fc_table = pd.read_csv(fc, sep="\t")
column_names = list(fc_table.columns)
column_dict = {key: i for i, key in enumerate(column_names)}
column_dict = {key: column_dict[key] - 1 for key in column_dict} # first sample column should have value 1

sample = comparison.split("_vs_")[0]

if not "-" in sample:
    sample_column = column_dict[sample]

else: # multiple samples will be pooled
    # get individual sample names
    samples = sample.split("-")
    
    # get column number for each sample name
    sample_columns = []
    for i in samples:
        number = str(column_dict[i])
        sample_columns.append(number)
    
    # convert to comma-seperated string
    sample_column = ",".join(sample_columns)

shell(
    "python {b2dir}/BAGEL.py bf -i {fc} "
    " -o {bf} " 
    "-e {eg} " 
    "-n {neg} "
    "-c {sample_column} {log}"
)   






