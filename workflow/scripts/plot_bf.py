#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

comparison = snakemake.wildcards[0]

# load Bayes factors
df = pd.read_table(snakemake.input[0])

# remove genes with extreme negative BF values (controls) for better plotting
df = df[df["BF"] > -50 ]

#create BF plot
sns.set_style("white")
sns.set_style("ticks")
sns.despine()
df.hist("BF", 
        bins=50,
        color="seagreen")
plt.xlabel('Bayes Factor')
plt.ylabel('Number of Genes')
plt.savefig(snakemake.output[0])
plt.clf()





