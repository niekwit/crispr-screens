import sys
import os
import pandas as pd
import numpy as np
import general_functions as utils

# write stdout/stderr to log file
sys.stderr = open(snakemake.log[0], "w")

# create df to store alignment rates
df = pd.DataFrame(columns=["sample","alignment_rate"], index=np.arange(len(snakemake.input)))
samples = []
rates = []

for i in sorted(snakemake.input):
    # get sample name from file name
    sample = os.path.basename(i).replace(".log","")
    samples.append(sample)
    
    # extract alignment rate from file
    with open(i) as f:
        lines = f.readlines()
    rate = float([x for x in lines if "overall alignment rate" in x][0].replace("% overall alignment rate\n",""))
    rates.append(rate)

# add values to empty df for plotting
df["sample"] = samples
df["alignment_rate"] = rates

# plot alignment rate
utils.plot(df,"Overall alignment rate (%)", snakemake.output[0])



