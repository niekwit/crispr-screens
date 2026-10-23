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


# NCBI taxon IDs for the species this workflow supports — same mapping used
# by workflow/scripts/string_db.py when it originally fetched the network.
TAXON_IDS = {"human": 9606, "mouse": 10090}


def taxon_id_for(species):
    return TAXON_IDS.get(species)


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


def load_id_mapping(path):
    """query_gene -> preferred_name, as resolved by STRING (see
    string_db.py) — the preferred name is what labels a node in the
    network SVG, which can differ from our submitted identifier (e.g. a
    deprecated alias resolving to its current official symbol)."""
    if not path or not os.path.exists(path):
        return {}
    try:
        df = pd.read_csv(path)
        return dict(zip(df["query_gene"].astype(str), df["preferred_name"].astype(str)))
    except Exception as exc:  # pragma: no cover - defensive
        logging.warning(f"Could not load STRING-db id mapping from {path}: {exc}")
        return {}


def _render_hover_html(hovertemplate, customdata_cols, row):
    """Renders the same per-point hover text used by the Plotly scatter
    (minus the %{customdata[i]} placeholders) as plain HTML, so the STRING
    network's nodes can show identical stats in a custom tooltip."""
    html = hovertemplate
    for i, col in enumerate(customdata_cols):
        html = html.replace(f"%{{customdata[{i}]}}", str(row[col]))
    return html


def build_plot(
    div_id,
    title,
    term_panel_title,
    term_list_id,
    df,
    y_col,
    y_title,
    color_options,
    id_col,
    customdata_cols,
    hovertemplate,
    terms,
    direction,
    id_mapping=None,
):
    """One Scattergl figure + the JS payload needed to drive it (gene
    search, click-to-label, STRING-db term highlighting, and switching
    which column governs marker fill colour).

    color_options: ordered list of (key, label, series) tuples — series is
    a column of df, same length/order. The first option is the initial
    marker colour; all are embedded so the "Color by" dropdown can switch
    between them client-side with no server round-trip.
    """
    default_key, default_label, default_series = color_options[0]
    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=df["x"],
            y=df[y_col],
            mode="markers",
            marker=dict(
                color=default_series,
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title=default_label),
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

    hover_by_gene = {
        str(r[id_col]): _render_hover_html(hovertemplate, customdata_cols, r)
        for _, r in df.iterrows()
    }
    # Also key by the name STRING actually displays on the node, when it
    # differs from our submitted identifier (see load_id_mapping).
    for query_gene, preferred_name in (id_mapping or {}).items():
        if query_gene in hover_by_gene and preferred_name not in hover_by_gene:
            hover_by_gene[preferred_name] = hover_by_gene[query_gene]

    color_values_by_option = {
        key: [round(float(v), 4) for v in series.tolist()]
        for key, _, series in color_options
    }

    return {
        "div_id": div_id,
        "panel_title": title,
        "term_panel_title": term_panel_title,
        "term_list_id": term_list_id,
        "fig": fig,
        "genes": df[id_col].astype(str).tolist(),
        "color_options": [(key, label) for key, label, _ in color_options],
        "color_values_by_option": color_values_by_option,
        "default_color_key": default_key,
        "gene_data": [
            {"id": str(r[id_col]), "x": r["x"], "y": r[y_col]} for _, r in df.iterrows()
        ],
        "terms": terms,
        "direction": direction,
        "hover_by_gene": hover_by_gene,
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


def build_svg_panel(
    div_id, title, svg_path, direction, stats_div_id=None, taxon_id=None
):
    """An inline-embedded STRING-db network SVG — not a Plotly figure, so
    it's rendered as raw markup rather than via Plotly.newPlot.

    stats_div_id: the div_id of the sibling scatter plot (same direction)
    whose per-gene hover stats should be shown when hovering a network
    node — STRING tags each node with the gene's name
    (data-safe_div_label), matched against that plot's hover_by_gene map.

    taxon_id: NCBI taxon ID (see taxon_id_for), used to query STRING-db's
    own API live, client-side, for a gene info popup when a node is
    clicked — not pre-fetched/embedded, so it needs internet access when
    the report is viewed (unlike everything else in it), but keeps the
    report from having to carry per-gene description text for every gene
    in every network up front.
    """
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
        "stats_div_id": stats_div_id,
        "taxon_id": taxon_id,
    }
