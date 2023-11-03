import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

import general_functions as utils

mpl.use('agg')


if snakemake.params[0] == "plot_alignment_rate":

    #create df to store alignment rates
    df = pd.DataFrame(columns=["sample","alignment_rate"],index=np.arange(len(snakemake.input)))
    samples = []
    rates = []

    for i in sorted(snakemake.input):
        #get sample name from file name
        sample = os.path.basename(i).replace(".log","")
        samples.append(sample)
        
        #extract alignment rate from file
        with open(i) as f:
            lines = f.readlines()
        rate = float([x for x in lines if "overall alignment rate" in x][0].replace("% overall alignment rate\n",""))
        rates.append(rate)

    #add values to empty df for plotting
    df["sample"] = samples
    df["alignment_rate"] = rates

    #plot alignment rate
    plot(df,"Overall alignment rate (%)", snakemake.output[0])

elif snakemake.params[0] == "plot_coverage":

    #get number of sgRNAs in CRISPR library
    fasta = pd.read_table(snakemake.params[1], header = None)
    lib_size = len(fasta) / 2

    #extract number of single mapped aligned reads from counts-aggregated.tsv
    df = pd.read_table("count/counts-aggregated.tsv", sep = "\t")
    column_names = list(df.columns)
    del column_names[0:2] #delete sgRNA and gene columns
    
    counts = {}
    for i in column_names:
        count_sum = []
        count_sum.append(df[i].sum())
        counts[i] = count_sum
    
    #convert counts to sequence coverage
    df = pd.DataFrame(counts)
    df = df / lib_size
    
    #order columns alphabetically
    df = df.reindex(sorted(df.columns), axis=1)
    
    #transpose data frame
    df = df.transpose()
    df["samples"] = df.index
    names = ["coverage", "samples"]
    df.columns = names
    df = df[["samples","coverage"]]
    
    #plot coverage per sample
    plot(df,"Fold sequence coverage per sample",snakemake.output[0])






