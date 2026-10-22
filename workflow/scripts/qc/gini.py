"""Lorenz curve per sample: evenness of sgRNA representation, with the Gini
index for each sample shown in the legend."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .common import CUTOFF_COLOR, sample_color, style_fig


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


def lorenz_curve(values):
    """Cumulative share of sgRNAs (x) vs. cumulative share of reads (y),
    ranked from the lowest to the highest count."""
    x = np.sort(np.asarray(values, dtype=float))
    n = len(x)
    total = x.sum()
    cum_reads = np.cumsum(x) / total if total else np.zeros(n)
    cum_sgrnas = np.arange(1, n + 1) / n
    return np.concatenate(([0], cum_sgrnas)), np.concatenate(([0], cum_reads))


def build(counts, sample_cols, samples):
    gini_df = pd.DataFrame(
        {
            "sample": samples,
            "gini_index": [gini_index(counts[c].values) for c in sample_cols],
        }
    )

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            line=dict(color=CUTOFF_COLOR, dash="dash", width=1.5),
            name="Perfect equality",
            hoverinfo="skip",
        )
    )
    for i, (col, sample) in enumerate(zip(sample_cols, samples)):
        cum_sgrnas, cum_reads = lorenz_curve(counts[col].values)
        gini = gini_df.loc[gini_df["sample"] == sample, "gini_index"].iloc[0]
        fig.add_trace(
            go.Scatter(
                x=cum_sgrnas,
                y=cum_reads,
                mode="lines",
                line=dict(color=sample_color(i), width=2.5),
                name=f"{sample} (Gini={gini:.3f})",
                hovertemplate=(
                    f"<b>{sample}</b><br>"
                    "Cumulative sgRNAs: %{x:.1%}<br>"
                    "Cumulative reads: %{y:.1%}<extra></extra>"
                ),
            )
        )

    fig.update_xaxes(title="Cumulative share of sgRNAs", range=[0, 1], tickformat=".0%")
    fig.update_yaxes(title="Cumulative share of reads", range=[0, 1], tickformat=".0%")
    style_fig(
        fig,
        height=600,
        showlegend=True,
        margin=dict(t=60, b=80, l=70, r=30),
        legend=dict(title="Sample (Gini index)"),
    )

    return dict(
        id="gini-index",
        title="Gini index",
        description=(
            "Lorenz curve of sgRNA read counts per sample: the cumulative share of "
            "reads plotted against the cumulative share of sgRNAs, ranked from "
            "lowest to highest count. The dashed diagonal is perfect equality — "
            "the further a curve bows below it, the more uneven the sgRNA "
            "representation. The Gini index, shown per sample in the legend, "
            "summarises this as a single value between 0 (perfectly even) and 1 "
            "(maximally uneven)."
        ),
        fig=fig,
        df=gini_df,
        filename="gini-index",
    )
