"""Number of sgRNAs with zero reads per sample."""

import pandas as pd

from .common import bar_fig


def build(counts, sample_cols, samples):
    df = pd.DataFrame(
        {
            "sample": samples,
            "missed_sgrnas": [int((counts[c] == 0).sum()) for c in sample_cols],
        }
    )
    fig = bar_fig(df, "sample", "missed_sgrnas", "Missed sgRNAs", "Missed sgRNAs")
    return dict(
        id="missed-sgrnas",
        title="Missed sgRNAs",
        description="Number of sgRNAs with zero reads per sample.",
        fig=fig,
        df=df,
        filename="missed-sgrnas",
    )
