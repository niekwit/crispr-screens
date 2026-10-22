"""Gini index per sample: evenness of sgRNA representation."""

import numpy as np
import pandas as pd

from .common import bar_fig


def gini_index(values):
    """Uncorrected Gini coefficient, equivalent to R's ineq::Gini(corr = FALSE),
    which was used by the previous version of this report."""
    x = np.sort(np.asarray(values, dtype=float))
    n = len(x)
    total = x.sum()
    if total == 0:
        return 0.0
    ranks = np.arange(1, n + 1)
    return (2 * np.sum(ranks * x)) / (n * total) - (n + 1) / n


def build(counts, sample_cols, samples):
    df = pd.DataFrame(
        {
            "sample": samples,
            "gini_index": [gini_index(counts[c].values) for c in sample_cols],
        }
    )
    fig = bar_fig(df, "sample", "gini_index", "Gini index", "Evenness of sgRNAs")
    return dict(
        id="gini-index",
        title="Gini index",
        description=(
            "Evenness of sgRNA representation per sample. Values close to 0 "
            "indicate an even distribution; values close to 1 indicate strong skew."
        ),
        fig=fig,
        df=df,
        filename="gini-index",
    )
