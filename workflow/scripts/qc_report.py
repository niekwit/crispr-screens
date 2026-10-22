"""Build a single interactive HTML QC report (overview, alignment rate,
sequence coverage, Gini index, missed sgRNAs, read count distribution,
sgRNA abundance) and the alignment-rates.csv file consumed by ngs-tracker.

Each metric lives in its own module under qc/ (one dict with id, title,
description, fig and df) so new QC metrics can be added without touching
this file beyond adding one entry to `sections`.
"""

import logging
import os

import pandas as pd

from qc import (
    abundance,
    alignment,
    coverage,
    count_distribution,
    gini,
    missed_sgrnas,
    report_html,
    summary,
)
from qc.common import clean_sample_name

log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

log_files = snakemake.input["logs"]
counts_file = snakemake.input["counts"]
fasta_file = snakemake.params["fasta"]
output_html = snakemake.output["html"]
output_alignment_csv = snakemake.output["csv"]
output_count_dist_csv = snakemake.output["count_csv"]

# ── Alignment rate (also written to alignment-rates.csv for ngs-tracker) ─────
logging.info("Parsing Bowtie alignment logs")
alignment_section = alignment.build(log_files, output_alignment_csv)
logging.info(f"Written alignment rates to {output_alignment_csv}")

# ── Load count table (shared by coverage, Gini, missed sgRNAs, distribution) ─
logging.info(f"Loading {counts_file}")
counts = pd.read_csv(counts_file, sep="\t")
sample_cols = list(counts.columns[2:])
samples = [clean_sample_name(c) for c in sample_cols]

# ── Assemble report ───────────────────────────────────────────────────────────
coverage_section = coverage.build(counts, sample_cols, samples, fasta_file)
gini_section = gini.build(counts, sample_cols, samples)
missed_section = missed_sgrnas.build(counts, sample_cols, samples)
count_dist_section = count_distribution.build(
    counts, sample_cols, samples, output_count_dist_csv
)
logging.info(f"Written count distribution stats to {output_count_dist_csv}")
abundance_section = abundance.build(counts, sample_cols, samples)

summary_section = summary.build(
    alignment_section,
    coverage_section,
    gini_section,
    missed_section,
    count_dist_section,
)

sections = [
    summary_section,
    alignment_section,
    coverage_section,
    gini_section,
    missed_section,
    count_dist_section,
    abundance_section,
]

html = report_html.render(sections)

os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as fh:
    fh.write(html)
logging.info(f"Written QC report to {output_html}")
