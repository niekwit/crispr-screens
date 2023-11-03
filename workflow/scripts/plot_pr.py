#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

if os.path.exists(f"bagel2/{comparison}/{comparison}.no_bf"):
    
    #create empty file
    open(pr, 'w').close()

comparison = snakemake.wildcards[0]

print(f"BAGEL2: plotting precision-recall for {comparison}")

#load PR
df = pd.read_table(snakemake.input[0])

#create plot
sns.set_style("white")
sns.set_style("ticks")
sns.despine()
plt.plot(df.Recall, df.Precision, 
         linewidth=2,
         color="seagreen")
plt.xlim(0,1.01)
plt.ylim(0,1.02)
plt.xlabel('Recall')
plt.ylabel('Precision (1-FDR)')
plt.savefig(snakemake.output[0])
plt.clf()




