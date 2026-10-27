import logging
from snakemake.shell import shell

b2dir = snakemake.input["b2dir"]
species = snakemake.params["species"]
fc = snakemake.input["fc"]
bf = snakemake.output["bf"]
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

# Load gene sets: BAGEL2 default list per species, unless a custom list is set
default_lists = {
    "human": ("CEGv2.txt", "NEGv1.txt"),  # essential, non-essential genes
    "mouse": ("CEG_mouse.txt", "NEG_mouse.txt"),
}
if "none" in (ceg, cneg) and species not in default_lists:
    raise ValueError(
        f"No BAGEL2 default gene lists for species {species!r}: "
        "provide custom_gene_lists in config.yml"
    )
eg = f"{b2dir}/{default_lists[species][0]}" if ceg == "none" else ceg
neg = f"{b2dir}/{default_lists[species][1]}" if cneg == "none" else cneg

# Prepare command and log, and run
# Column number for -c is always 1, as CRISPRcleanR
# has averaged the samples, if multiple replicates were provided
command = f"python {b2dir}/BAGEL.py bf -i {fc} -o {bf} -e {eg} -n {neg} -c 1 {extra} > {log_stdout} 2> {log_stderr}"
logging.debug(command)
shell(command)
