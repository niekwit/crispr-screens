import logging
from snakemake.shell import shell

b2dir = snakemake.input["b2dir"]
species = snakemake.params["species"]
fc = snakemake.input["fc"]
bf = snakemake.output["bf"]
comparison = snakemake.wildcards[0]
ceg = snakemake.params["ceg"]
cneg = snakemake.params["cneg"]
extra = snakemake.params["extra"]

# Set up logging
log_stdout = snakemake.log["stdout"]
log_stderr = snakemake.log["stderr"]
logging.basicConfig(
    format="%(levelname)s:%(message)s",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(snakemake.log["command"])],
)

# Load gene sets
if ceg == "none" and cneg == "none":
    if species == "human":
        eg = f"{b2dir}/CEGv2.txt"  # essential genes
        neg = f"{b2dir}/NEGv1.txt"  # non-essential genes

    elif species == "mouse":
        eg = f"{b2dir}/CEG_mouse.txt"  # essential genes
        neg = f"{b2dir}/NEG_mouse.txt"  # non-essential genes
else:
    eg = ceg
    neg = cneg

# Prepare command and log, and run
# Column number for -c is always 1, as CRISPRcleanR
# has averaged the samples, if multiple replicates were provided
command = f"python {b2dir}/BAGEL.py bf -i {fc} -o {bf} -e {eg} -n {neg} -c 1 {extra} > {log_stdout} 2> {log_stderr}"
logging.debug(command)
shell(command)
