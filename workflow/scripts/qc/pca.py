"""PCA of per-sample sgRNA count profiles: samples that cluster together
have similar sgRNA representation; an outlier may indicate a batch effect,
sample swap, or other technical problem."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sklearn.decomposition import PCA

from .common import sample_color, style_fig

MIN_SAMPLES = 3  # fewer than this and a 2D PCA isn't meaningful


def build(counts, sample_cols, samples):
    if len(samples) < MIN_SAMPLES:
        return _too_few_samples(len(samples))

    log_counts = np.log2(counts[sample_cols].astype(float) + 1).T.values

    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(log_counts)
    explained = 100 * pca.explained_variance_ratio_

    df = pd.DataFrame({"sample": samples, "PC1": coords[:, 0], "PC2": coords[:, 1]})

    fig = go.Figure(
        go.Scatter(
            x=df["PC1"],
            y=df["PC2"],
            mode="markers+text",
            text=df["sample"],
            textposition="top center",
            marker=dict(
                size=14,
                color=[sample_color(i) for i in range(len(samples))],
                line=dict(width=1.5, color="#1e293b"),
            ),
            hovertemplate="<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>",
        )
    )
    fig.update_xaxes(title=f"PC1 ({explained[0]:.1f}% variance)", zeroline=True)
    fig.update_yaxes(title=f"PC2 ({explained[1]:.1f}% variance)", zeroline=True)
    style_fig(fig, height=550, margin=dict(t=60, b=80, l=70, r=30))

    return dict(
        id="pca",
        title="PCA",
        description=(
            "Principal component analysis (scikit-learn) of log2(sgRNA count + 1) "
            "profiles per sample. Samples that cluster together have similar sgRNA "
            "representation; a sample far from its replicates may indicate a "
            "batch effect, sample swap, or other technical problem."
        ),
        fig=fig,
        df=df,
        filename="pca",
    )


def _too_few_samples(n_samples):
    fig = go.Figure()
    fig.add_annotation(
        text=f"PCA requires at least {MIN_SAMPLES} samples (found {n_samples}).",
        showarrow=False,
        font=dict(size=15, color="#64748b"),
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    style_fig(fig, height=300)

    return dict(
        id="pca",
        title="PCA",
        description=(
            "Principal component analysis of log2(sgRNA count + 1) profiles per "
            f"sample. Not enough samples ({n_samples}) to compute a meaningful 2D "
            f"projection — at least {MIN_SAMPLES} are needed."
        ),
        fig=fig,
        df=pd.DataFrame(columns=["sample", "PC1", "PC2"]),
        filename="pca",
    )
