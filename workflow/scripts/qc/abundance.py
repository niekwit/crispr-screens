"""Smoothed sgRNA read-count density per sample, overlaid, with the P90/P10
ratio (evenness of representation) shown in the legend."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import gaussian_kde

from .common import hex_to_rgba, sample_color, style_fig


def p90_p10(values):
    """Ratio of the 90th to the 10th percentile read count: a standard
    CRISPR-screen evenness metric (e.g. used by Broad GPP's PoolQ). Close to
    1 means even sgRNA representation; higher values mean a wider spread
    between the best- and worst-represented sgRNAs."""
    values = np.asarray(values, dtype=float)
    p10 = np.percentile(values, 10)
    p90 = np.percentile(values, 90)
    return p90 / p10 if p10 else float("inf")


def build(counts, sample_cols, samples):
    log_counts = {
        col: np.log10(counts[col].values.astype(float) + 1) for col in sample_cols
    }
    x_max = max(v.max() for v in log_counts.values()) * 1.05
    grid = np.linspace(0, x_max, 400)

    stats_rows = []
    fig = go.Figure()
    for i, (col, sample) in enumerate(zip(sample_cols, samples)):
        ratio = p90_p10(counts[col].values)
        stats_rows.append({"sample": sample, "p90_p10": ratio})

        density = gaussian_kde(log_counts[col])(grid)
        color = sample_color(i)
        fig.add_trace(
            go.Scatter(
                x=grid,
                y=density,
                mode="lines",
                line=dict(color=color, width=2),
                fill="tozeroy",
                fillcolor=hex_to_rgba(color, 0.12),
                name=f"{sample} | P90/P10 = {ratio:.2f}",
                hovertemplate=(
                    f"<b>{sample}</b><br>"
                    "log10(count+1): %{x:.2f}<br>"
                    "Density: %{y:.3f}<extra></extra>"
                ),
            )
        )
    stats_df = pd.DataFrame(stats_rows)

    fig.update_xaxes(title="log10(read count + 1)")
    fig.update_yaxes(title="Density", rangemode="tozero")
    style_fig(
        fig,
        height=600,
        showlegend=True,
        margin=dict(t=60, b=80, l=70, r=30),
        legend=dict(title="Sample (P90/P10)"),
    )

    return dict(
        id="sgrna-abundance",
        title="sgRNA abundance",
        description=(
            "Smoothed distribution (kernel density estimate) of sgRNA read counts "
            "per sample, overlaid for comparison. P90/P10 is the ratio of the 90th "
            "to the 10th percentile read count: a standard evenness metric for "
            "CRISPR screens, shown per sample in the legend. A value close to 1 "
            "means almost all sgRNAs have similar counts (even representation); a "
            "high value means a wide gap between the best- and worst-represented "
            "sgRNAs (uneven representation). Click a legend entry to hide or show "
            "that sample's curve."
        ),
        fig=fig,
        df=stats_df,
        filename="sgrna-abundance",
    )
