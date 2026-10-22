"""Replicate correlation heatmap: Pearson correlation of log-transformed
sgRNA counts between every pair of samples."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .common import FLAG_COLOR, style_fig

# Diverging scale: rose for negative correlation, near-white at 0, indigo for
# positive — matches the report's palette instead of a default D3/ColorBrewer
# scale.
CORRELATION_COLORSCALE = [
    [0.0, FLAG_COLOR],
    [0.5, "#f8fafc"],
    [1.0, "#6366f1"],
]


def build(counts, sample_cols, samples):
    log_counts = np.log2(counts[sample_cols].astype(float) + 1)
    log_counts.columns = samples
    corr = log_counts.corr(method="pearson")

    fig = go.Figure(
        go.Heatmap(
            z=corr.values,
            x=corr.columns.tolist(),
            y=corr.columns.tolist(),
            zmin=-1,
            zmax=1,
            colorscale=CORRELATION_COLORSCALE,
            colorbar=dict(title="Pearson r"),
            text=corr.round(3).values,
            texttemplate="%{text}",
            textfont=dict(size=13),
            hovertemplate="<b>%{x} vs %{y}</b><br>r = %{z:.3f}<extra></extra>",
        )
    )
    fig.update_xaxes(tickangle=45, showgrid=False, constrain="domain")
    fig.update_yaxes(
        autorange="reversed",
        showgrid=False,
        scaleanchor="x",
        scaleratio=1,
        constrain="domain",
    )
    style_fig(
        fig,
        height=max(450, 80 * len(samples) + 150),
        margin=dict(t=60, b=80, l=110, r=30),
    )

    corr_df = corr.reset_index().rename(columns={"index": "sample"})

    return dict(
        id="sample-correlation",
        title="Sample correlation",
        description=(
            "Pearson correlation of log2(sgRNA count + 1) between every pair of "
            "samples. High correlation (close to 1) between replicates indicates "
            "good reproducibility; a replicate that correlates poorly with its "
            "counterparts may point to a technical problem with that sample."
        ),
        fig=fig,
        df=corr_df,
        filename="sample-correlation",
    )
