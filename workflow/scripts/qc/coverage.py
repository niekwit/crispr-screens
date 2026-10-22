"""Fold sequence coverage per sample (total aligned reads / library size)."""

import pandas as pd

from .common import bar_fig


def build(counts, sample_cols, samples, fasta_file):
    with open(fasta_file) as fh:
        n_lines = sum(1 for _ in fh)
    sgrnas = n_lines / 2

    df = pd.DataFrame(
        {
            "sample": samples,
            "coverage": [counts[c].sum() / sgrnas for c in sample_cols],
        }
    )
    fig = bar_fig(
        df, "sample", "coverage", "Fold sequence coverage", "Sequence coverage"
    )
    return dict(
        id="sequence-coverage",
        title="Sequence coverage",
        description=(
            "Fold sequence coverage of the sgRNA library per sample (total aligned "
            "reads / library size). 200-500x is typically recommended."
        ),
        fig=fig,
        df=df,
        filename="sequence-coverage",
    )
