#!/usr/bin/env python

import os
import logging
from snakemake.shell import shell

b2dir = snakemake.input["b2dir"]
species = snakemake.params["species"]
bf = snakemake.input["bf"]
pr = snakemake.output[0]
ceg = snakemake.params["ceg"]
cneg = snakemake.params["cneg"]
extra = snakemake.params["extra"]

# Set up logging
log_stdout = snakemake.log["stdout"]
log_stderr = snakemake.log["stderr"]
logging.basicConfig(format='%(levelname)s:%(message)s', 
                    level=logging.DEBUG,
                    handlers=[logging.FileHandler(snakemake.log["command"])])

# Get comparison
comparison = os.path.basename(bf.replace(".bf",""))

# Load gene sets
if ceg == "none" and cneg == "none":
    if species == "human":
        eg = f"{b2dir}/CEGv2.txt" # essential genes
        neg = f"{b2dir}/NEGv1.txt" # non-essential genes

    elif species == "mouse":
        eg = f"{b2dir}/CEG_mouse.txt" # essential genes
        neg = f"{b2dir}/NEG_mouse.txt" # non-essential genes
else:
    eg = ceg
    neg = cneg

# Prepare command and log, and run
command = f"python {b2dir}/BAGEL.py pr -i {bf} -o {pr} -e {eg} -n {neg} {extra} > {log_stdout} 2> {log_stderr}"
logging.debug(command)
shell(command)
