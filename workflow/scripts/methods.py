"""
Methods section of the run (snakemake-methods): the prose is in
workflow/methods/methods.md.j2. This script adds the values of the library
file (number of sgRNAs and genes) and the NCBI taxon of the species.
"""

import csv
import os

from snakemake_methods.snakemake import run

config = snakemake.config
library_file = snakemake.params["library_file"]

# sgRNAs and genes of the library (columns as in csv_to_fasta.py)
name_column = config["csv"]["name_column"]
gene_column = config["csv"]["gene_column"]
sgrnas, genes = set(), set()
with open(library_file, newline="") as fh:
    rows = csv.reader(fh)
    next(rows)
    for row in rows:
        sgrnas.add(row[name_column])
        if row[gene_column]:
            genes.add(row[gene_column])

# NCBI taxon IDs as in string_db.py
TAXA = {"human": 9606, "mouse": 10090}

PREAMBLE = (
    "*Generated on {date} from the configuration of the run (appended below). "
    "Software versions are those of the run: installed in its conda "
    "environments (including those of the Snakemake wrappers), git commits of "
    "tools cloned at run time, and versions reported by online services.*"
)

run(
    snakemake,
    preamble=PREAMBLE,
    context={
        "library_name": os.path.basename(library_file),
        "n_sgrnas": len(sgrnas),
        "n_genes": len(genes),
        "taxon": TAXA.get(config["lib_info"]["species"]),
    },
)
