# Snakemake workflow: `crispr-screens`

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10286661.svg)](https://doi.org/10.5281/zenodo.10286661)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/crispr-screens/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/crispr-screens/actions/workflows/main.yml)


<p align="center">
  <img src="docs/_static/logo2_small.png" width="200" alt="GPSW Logo" />
</p>


A Snakemake workflow for the analysis of CRISPR screens.

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

Instructions of how to use `crispr-screens` can be found here:

https://crispr-screens.readthedocs.io/en/latest/

## Annotating sgRNAs with genomic coordinates

CRISPRcleanR, which is always run before BAGEL2 (and optionally before MAGeCK and DrugZ), needs the genomic coordinates of every sgRNA. Some libraries do not provide them. The standalone script `annotate_sgrna_coordinates.py` (not part of the Snakemake workflow) finds them by searching the sgRNA sequences in a genome FASTA file and writes a CRISPRcleanR library file.

See the [documentation](https://crispr-screens.readthedocs.io/en/latest/annotate_sgrna_coordinates/annotate_sgrna_coordinates.html) for usage, options, output format, and test results.
