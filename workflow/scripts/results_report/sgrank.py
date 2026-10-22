"""sgRNA-rank view per MAGeCK comparison: for the top N enriched/depleted
genes (by MAGeCK gene-level rank), shows each gene's individual sgRNA log
fold changes against the overall sgRNA LFC distribution.

Ported from workflow/scripts/plot_sgrank.R (itself inspired by
MAGeCKFlute::sgRankView and https://www.jci.org/articles/view/201639
Figure 1D) as an interactive Plotly version. Unlike the enriched/depleted
scatter plots, this one has no STRING-db panel or gene search — it's a
fixed top-N-genes view, not a searchable scatter of every gene — so it's
built as a "simple" plot (see common.build_simple_plot).
"""

import logging

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde

from . import common

TOP_N = 5  # number of top genes to show for enriched/depleted
ENRICHED_COLOR = "#e41a1c"
DEPLETED_COLOR = "#377eb8"
ROW_HEIGHT = 0.8
ROW_GAP = 1.0  # centre-to-centre spacing between gene rows


def _rank_plot(df_all, df_top, genes, tick_color, title):
    n_genes = len(genes)
    x_min, x_max = df_all["LFC"].min(), df_all["LFC"].max()
    pad = (x_max - x_min) * 0.02
    x_range = [x_min - pad, x_max + pad]
    grid = np.linspace(x_range[0], x_range[1], 300)
    density = gaussian_kde(df_all["LFC"].dropna().values)(grid)

    row_of = {gene: i + 1 for i, gene in enumerate(genes)}  # 1-indexed, bottom-up

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.25, 0.75],
        vertical_spacing=0.03,
    )

    # Reference density of the full sgRNA LFC distribution.
    fig.add_trace(
        go.Scatter(
            x=grid,
            y=density,
            mode="lines",
            fill="tozeroy",
            line=dict(color="#888888", width=1),
            fillcolor="rgba(224, 224, 224, 0.6)",
            hoverinfo="skip",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Same reference density repeated as a background shade behind every
    # gene row (white = low density, black = high density).
    fig.add_trace(
        go.Heatmap(
            x=grid,
            y=[row_of[g] for g in genes],
            z=[density for _ in genes],
            colorscale=[[0, "white"], [1, "black"]],
            showscale=False,
            hoverinfo="skip",
        ),
        row=2,
        col=1,
    )

    # One sgRNA tick (vertical segment) per row, drawn as a single traced
    # polyline with None-separated breaks between segments.
    xs, ys, hover = [], [], []
    for _, r in df_top.iterrows():
        y0 = row_of[r["Gene"]] - ROW_HEIGHT / 2
        y1 = row_of[r["Gene"]] + ROW_HEIGHT / 2
        xs += [r["LFC"], r["LFC"], None]
        ys += [y0, y1, None]
        hover.append(f"{r['sgrna']} ({r['Gene']})<br>LFC: {r['LFC']:.3f}")
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="lines",
            line=dict(color=tick_color, width=1.5),
            hoverinfo="skip",
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    # Invisible markers at each tick's midpoint purely to carry hover text.
    fig.add_trace(
        go.Scatter(
            x=df_top["LFC"],
            y=[row_of[g] for g in df_top["Gene"]],
            mode="markers",
            marker=dict(size=10, opacity=0),
            hovertext=hover,
            hoverinfo="text",
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    for gene, row in row_of.items():
        fig.add_shape(
            type="rect",
            xref="x2",
            yref="y2",
            x0=x_range[0],
            x1=x_range[1],
            y0=row - ROW_HEIGHT / 2,
            y1=row + ROW_HEIGHT / 2,
            line=dict(color="#444444", width=1),
            fillcolor="rgba(0,0,0,0)",
        )

    fig.update_xaxes(range=x_range, row=1, col=1, showticklabels=False)
    fig.update_xaxes(title="Log2(Fold change)", range=x_range, row=2, col=1)
    fig.update_yaxes(title="Density", row=1, col=1)
    fig.update_yaxes(
        tickvals=[row_of[g] for g in genes],
        ticktext=genes,
        range=[0.5, n_genes + 0.5],
        row=2,
        col=1,
    )
    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=15)),
        template="plotly_white",
        height=180 + 45 * n_genes,
        margin=dict(t=60, b=60, l=100, r=20),
        showlegend=False,
    )
    return fig


def build_plots(gene_txt, sgrna_txt, key, comparison, cnv):
    """Returns a list of 0-2 "simple" plot dicts (enriched/depleted) to
    append to the MAGeCK entry for this (comparison, cnv)."""
    logging.info(f"Loading {sgrna_txt}")
    gene_df = pd.read_csv(gene_txt, sep="\t")[["id", "neg|rank", "pos|rank"]]
    gene_df = gene_df.rename(
        columns={"id": "Gene", "neg|rank": "neg_rank", "pos|rank": "pos_rank"}
    )
    sg_df = pd.read_csv(sgrna_txt, sep="\t")
    data = sg_df.merge(gene_df, on="Gene", how="left")

    df_enriched = data[data["pos_rank"] <= TOP_N].sort_values("pos_rank")
    df_depleted = data[data["neg_rank"] <= TOP_N].sort_values("neg_rank")

    # R's rev(unique(...)): rank order (best first), then reversed so the
    # best-ranked gene ends up in the last (topmost, highest y) row.
    genes_enriched = list(dict.fromkeys(df_enriched["Gene"]))[::-1]
    genes_depleted = list(dict.fromkeys(df_depleted["Gene"]))[::-1]

    plots = []
    if genes_enriched:
        fig = _rank_plot(
            data,
            df_enriched,
            genes_enriched,
            ENRICHED_COLOR,
            f"sgRNA rank — enriched genes — {comparison} ({cnv})",
        )
        plots.append(
            common.build_simple_plot(
                div_id=f"{key}-sgrank-enriched",
                title="sgRNA rank (enriched)",
                fig=fig,
                direction="enriched",
            )
        )
    if genes_depleted:
        fig = _rank_plot(
            data,
            df_depleted,
            genes_depleted,
            DEPLETED_COLOR,
            f"sgRNA rank — depleted genes — {comparison} ({cnv})",
        )
        plots.append(
            common.build_simple_plot(
                div_id=f"{key}-sgrank-depleted",
                title="sgRNA rank (depleted)",
                fig=fig,
                direction="depleted",
            )
        )
    return plots
