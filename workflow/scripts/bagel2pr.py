#!/usr/bin/env python

import os
from snakemake.shell import shell

b2dir = snakemake.input["b2dir"]
species = snakemake.params["species"]
bf = snakemake.input["bf"]
pr = snakemake.output[0]
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

#get comparison
comparison = os.path.basename(bf.replace(".bf",""))

#load gene sets
if species == "human":
    eg = f"{b2dir}/CEGv2.txt" #essential genes
    neg = f"{b2dir}/NEGv1.txt" #non-essential genes

elif species == "mouse":
    eg = f"{b2dir}/CEG_mouse.txt" #essential genes
    neg = f"{b2dir}/NEG_mouse.txt" #non-essential genes

# run BAGEL2 pr
shell(
    "python {b2dir}/BAGEL.py pr -i {bf} " 
    "-o {pr} " 
    "-e {eg} " 
    "-n {neg} {log}"
)





