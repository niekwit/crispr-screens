"""sgRNA-rank view per MAGeCK comparison: for a chosen set of genes (the
top N enriched/depleted by MAGeCK gene-level rank, initially), shows each
gene's individual sgRNA log fold changes against the overall sgRNA LFC
distribution.

Ported from workflow/scripts/plot_sgrank.R (itself inspired by
MAGeCKFlute::sgRankView and https://www.jci.org/articles/view/201639
Figure 1D) as an interactive Plotly version. Unlike the enriched/depleted
scatter plots, this one has no STRING-db panel or gene search in the usual
sense — it's a fixed-genes view, not a searchable scatter of every gene —
so it's built as a "simple" plot (see common.build_simple_plot), with the
gene list itself editable client-side (see rank_data/default_genes below
and the rebuildSgRank() JS in report_html.py): the figure only renders the
genes typed into the "Genes to show" box above it, not a fixed top N, so
the data for every gene has to be available for the client to redraw from,
not just the N shown by default.
"""

import logging

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde

from . import common

DEFAULT_TOP_N = 5  # number of top genes shown initially, for enriched/depleted
# A gene's rank can disagree with some of its own individual sgRNAs — ticks
# are coloured by each sgRNA's own LFC sign, not by which direction
# (enriched/depleted) the surrounding plot is for.
POS_COLOR = "#e41a1c"  # LFC > 0
NEG_COLOR = "#377eb8"  # LFC < 0
ROW_HEIGHT = 0.8
ROW_GAP = 1.0  # centre-to-centre spacing between gene rows


def _tick_trace(df_subset, row_of, color, name):
    xs, ys, hover = [], [], []
    for _, r in df_subset.iterrows():
        y0 = row_of[r["Gene"]] - ROW_HEIGHT / 2
        y1 = row_of[r["Gene"]] + ROW_HEIGHT / 2
        xs += [r["LFC"], r["LFC"], None]
        ys += [y0, y1, None]
        hover.append(f"{r['sgrna']} ({r['Gene']})<br>LFC: {r['LFC']:.3f}")
    return go.Scatter(
        x=xs,
        y=ys,
        mode="lines",
        line=dict(color=color, width=1.5),
        name=name,
        hoverinfo="skip",
        showlegend=True,
    )


def _rank_plot(grid, density, x_range, df_top, genes, title):
    n_genes = len(genes)
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
    # gene row (white = low density, black = high density). Without ygap,
    # Plotly sizes each cell to the full row-to-row spacing (ROW_GAP),
    # which paints straight through the gap between adjacent rows' outline
    # rectangles instead of stopping at each one's own border.
    fig.add_trace(
        go.Heatmap(
            x=grid,
            y=[row_of[g] for g in genes],
            z=[density for _ in genes],
            colorscale=[[0, "white"], [1, "black"]],
            showscale=False,
            hoverinfo="skip",
            ygap=4,
        ),
        row=2,
        col=1,
    )

    # One sgRNA tick (vertical segment) per row, as two traced polylines
    # (None-separated breaks between segments) split by each sgRNA's own
    # LFC sign, not by the plot's overall enriched/depleted direction.
    pos = df_top[df_top["LFC"] >= 0]
    neg = df_top[df_top["LFC"] < 0]
    if len(pos):
        fig.add_trace(
            _tick_trace(pos, row_of, POS_COLOR, "sgRNA LFC > 0"), row=2, col=1
        )
    if len(neg):
        fig.add_trace(
            _tick_trace(neg, row_of, NEG_COLOR, "sgRNA LFC < 0"), row=2, col=1
        )
    hover = [
        f"{r['sgrna']} ({r['Gene']})<br>LFC: {r['LFC']:.3f}"
        for _, r in df_top.iterrows()
    ]
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
        showlegend=True,
        legend=dict(orientation="h", x=0.5, xanchor="center", y=1.0, yanchor="bottom"),
    )
    return fig


def _build_rank_data(data):
    """The payload rebuildSgRank() (report_html.py) needs to redraw this
    plot client-side for an arbitrary gene list: the reference density
    curve (identical for every gene, computed once here rather than
    per-direction) and every gene's individual sgRNA LFCs, as a flat row
    list for compact JSON — not just the default top N, since the whole
    point is letting the gene selection change after the report is built."""
    x_min, x_max = data["LFC"].min(), data["LFC"].max()
    pad = (x_max - x_min) * 0.02
    x_range = [float(x_min - pad), float(x_max + pad)]
    grid = np.linspace(x_range[0], x_range[1], 300)
    density = gaussian_kde(data["LFC"].dropna().values)(grid)
    rows = data.dropna(subset=["LFC"])
    return {
        "grid": [round(float(v), 4) for v in grid],
        "density": [round(float(v), 6) for v in density],
        "xRange": x_range,
        # Parallel arrays rather than a list of {gene, sgrna, lfc} dicts —
        # one gene/sgrna/lfc triple per sgRNA across the WHOLE library (not
        # just the default top N, see module docstring), so repeating the
        # three key names per row would otherwise double report size on a
        # genome-wide screen.
        "rowGenes": [str(g) for g in rows["Gene"]],
        "rowSgrnas": [str(s) for s in rows["sgrna"]],
        "rowLfcs": [round(float(v), 4) for v in rows["LFC"]],
        "genes": sorted(data["Gene"].dropna().unique().tolist()),
    }


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

    rank_data = _build_rank_data(data)
    grid = np.array(rank_data["grid"])
    density = np.array(rank_data["density"])
    x_range = rank_data["xRange"]

    df_enriched = data[data["pos_rank"] <= DEFAULT_TOP_N].sort_values("pos_rank")
    df_depleted = data[data["neg_rank"] <= DEFAULT_TOP_N].sort_values("neg_rank")

    # Natural rank order (best first) — shown as-typed in the "Genes to
    # show" box, and in that order reversed internally so the best-ranked
    # gene ends up in the last (topmost, highest y) row, matching R's
    # rev(unique(...)).
    default_enriched = list(dict.fromkeys(df_enriched["Gene"]))
    default_depleted = list(dict.fromkeys(df_depleted["Gene"]))

    plots = []
    if default_enriched:
        fig = _rank_plot(
            grid,
            density,
            x_range,
            df_enriched,
            default_enriched[::-1],
            f"sgRNA rank — enriched genes — {comparison} ({cnv})",
        )
        plots.append(
            common.build_simple_plot(
                div_id=f"{key}-sgrank-enriched",
                title="sgRNA rank (enriched)",
                fig=fig,
                direction="enriched",
                rank_data=rank_data,
                default_genes=default_enriched,
            )
        )
    if default_depleted:
        fig = _rank_plot(
            grid,
            density,
            x_range,
            df_depleted,
            default_depleted[::-1],
            f"sgRNA rank — depleted genes — {comparison} ({cnv})",
        )
        plots.append(
            common.build_simple_plot(
                div_id=f"{key}-sgrank-depleted",
                title="sgRNA rank (depleted)",
                fig=fig,
                direction="depleted",
                rank_data=rank_data,
                default_genes=default_depleted,
            )
        )
    return plots
