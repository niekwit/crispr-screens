"""Shared helpers for building per-comparison interactive scatter plots
(MAGeCK, DrugZ, BAGEL2) that get assembled into one combined results report.
"""

import logging
import os
import re

import pandas as pd
import plotly.graph_objects as go

# Palette used for STRING-db term highlighting, shared across all three tools
# (identical to what each standalone report used before).
COLORS = [
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
    "#f781bf",
    "#17becf",
    "#bcbd22",
    "#8c564b",
]


def slugify(text):
    return re.sub(r"[^A-Za-z0-9_-]", "-", str(text))


def load_string_terms(path):
    if not path or not os.path.exists(path):
        return []
    try:
        df = pd.read_csv(path)
        if df.empty or "preferredNames" not in df.columns:
            return []
        terms = []
        for _, row in df.iterrows():
            genes = [
                g.strip() for g in str(row["preferredNames"]).split(",") if g.strip()
            ]
            if genes:
                terms.append(
                    {
                        "label": f"{row['description']} (FDR={float(row['fdr']):.3f})",
                        "genes": genes,
                    }
                )
        logging.info(f"Loaded {len(df)} STRING-db terms from {path}")
        return terms
    except Exception as exc:  # pragma: no cover - defensive, mirrors prior scripts
        logging.warning(f"Could not load STRING-db terms from {path}: {exc}")
        return []


def build_plot(
    div_id,
    title,
    term_panel_title,
    term_list_id,
    df,
    y_col,
    y_title,
    color_col,
    color_title,
    id_col,
    customdata_cols,
    hovertemplate,
    terms,
    direction,
):
    """One Scattergl figure + the JS payload needed to drive it (gene
    search, click-to-label, STRING-db term highlighting)."""
    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=df["x"],
            y=df[y_col],
            mode="markers",
            marker=dict(
                color=df[color_col],
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title=color_title),
                size=8,
                line=dict(width=0.5, color="black"),
            ),
            customdata=list(zip(*(df[c] for c in customdata_cols))),
            hovertemplate=hovertemplate + "<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=15)),
        height=500,
        template="plotly_white",
        margin=dict(t=70, b=60, l=70, r=20),
        xaxis_title="Random Index",
        yaxis_title=y_title,
        showlegend=True,  # needed so term legend traces show when added
    )

    return {
        "div_id": div_id,
        "panel_title": title,
        "term_panel_title": term_panel_title,
        "term_list_id": term_list_id,
        "fig": fig,
        "genes": df[id_col].astype(str).tolist(),
        "color_values": [round(float(v), 4) for v in df[color_col].tolist()],
        "gene_data": [
            {"id": str(r[id_col]), "x": r["x"], "y": r[y_col]} for _, r in df.iterrows()
        ],
        "terms": terms,
        "direction": direction,
    }


def build_simple_plot(div_id, title, fig, direction):
    """A plain zoomable/pannable Plotly chart with no STRING-db panel and
    no gene-search wiring (e.g. a fixed top-N-genes view rather than a
    searchable scatter of every gene)."""
    return {
        "div_id": div_id,
        "panel_title": title,
        "simple": True,
        "fig": fig,
        "direction": direction,
    }


def build_svg_panel(div_id, title, svg_path, direction):
    """An inline-embedded STRING-db network SVG — not a Plotly figure, so
    it's rendered as raw markup rather than via Plotly.newPlot."""
    if not svg_path or not os.path.exists(svg_path):
        return None
    try:
        with open(svg_path) as fh:
            svg = fh.read()
    except OSError as exc:  # pragma: no cover - defensive
        logging.warning(f"Could not read STRING-db SVG {svg_path}: {exc}")
        return None
    return {
        "div_id": div_id,
        "panel_title": title,
        "kind": "svg",
        "svg": svg,
        "direction": direction,
    }
