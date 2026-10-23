"""Shared helpers for building per-comparison interactive scatter plots
(MAGeCK, DrugZ, BAGEL2) that get assembled into one combined results report.
"""

import logging
import os
import re
import xml.etree.ElementTree as ET

import pandas as pd
import plotly.graph_objects as go

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

# STRING-db's node radius is always 20 in the SVGs it serves (confirmed
# across every network this workflow has generated) — the baseline for
# _resize_nodes_by_rank's scale factor.
DEFAULT_NODE_RADIUS = 20
MIN_NODE_RADIUS = 12
MAX_NODE_RADIUS = 32

# STRING-db's edge evidence-channel colours, verified against the live
# legend on a real string-db.org network page (/network/<id>, "Edges"
# section) rather than guessed — the exact same colours the /api/svg/network
# endpoint this workflow calls renders with.
EDGE_COLOR_LEGEND = [
    ("rgb(31,120,180)", "Databases"),
    ("rgb(166,206,227)", "Textmining"),
    ("rgb(227,26,28)", "Experiments"),
    ("rgb(254,178,76)", "Co-expression"),
    ("rgb(28,140,30)", "Gene neighborhood"),
    ("rgb(118,201,98)", "Gene co-occurrence"),
    ("rgb(178,223,138)", "Gene fusions"),
]

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
    rank_col=None,
):
    """One Scattergl figure + the JS payload needed to drive it (gene
    search, click-to-label, STRING-db term highlighting, and switching
    which column governs marker fill colour).

    rank_col: optional column holding this gene's rank (lower = more
    significant) — if given, a gene-name-keyed rank_by_gene dict is
    returned alongside hover_by_gene, for the paired STRING network panel
    to size nodes by (see build_svg_panel).

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

    rank_by_gene = {}
    if rank_col:
        rank_by_gene = {str(r[id_col]): int(r[rank_col]) for _, r in df.iterrows()}
        for query_gene, preferred_name in (id_mapping or {}).items():
            if query_gene in rank_by_gene and preferred_name not in rank_by_gene:
                rank_by_gene[preferred_name] = rank_by_gene[query_gene]

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
        "rank_by_gene": rank_by_gene,
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


def _resize_nodes_by_rank(svg_text, rank_by_gene):
    """Scale each node's circle, and its label's offset (so the text stays
    flush against the new edge), by that gene's rank — lower rank (more
    significant) draws bigger. Scaled relative to the best/worst rank
    actually present among the nodes drawn in THIS network, not rank's
    absolute scale (which can run into the thousands for genes that aren't
    shown here), so the size spread is always visually meaningful."""
    if not rank_by_gene:
        return svg_text
    try:
        root = ET.fromstring(svg_text)
    except ET.ParseError as exc:  # pragma: no cover - defensive
        logging.warning(f"Could not parse STRING-db SVG for node resizing: {exc}")
        return svg_text
    ns = {"svg": SVG_NS}
    nodes = root.findall(".//svg:g[@class='nwnodecontainer']", ns)
    ranks_present = [
        rank_by_gene[g.get("data-safe_div_label")]
        for g in nodes
        if g.get("data-safe_div_label") in rank_by_gene
    ]
    if not ranks_present:
        return svg_text
    min_rank, max_rank = min(ranks_present), max(ranks_present)

    for g in nodes:
        rank = rank_by_gene.get(g.get("data-safe_div_label"))
        if rank is None:
            continue
        norm = (
            1.0 if max_rank == min_rank else (max_rank - rank) / (max_rank - min_rank)
        )
        new_r = MIN_NODE_RADIUS + norm * (MAX_NODE_RADIUS - MIN_NODE_RADIUS)
        scale = new_r / DEFAULT_NODE_RADIUS
        cx = cy = None
        for circle in g.findall("svg:circle", ns):
            cx, cy = float(circle.get("cx")), float(circle.get("cy"))
            circle.set("r", f"{new_r:.2f}")
        if cx is None:
            continue
        for text in g.findall("svg:text", ns):
            x, y = float(text.get("x")), float(text.get("y"))
            text.set("x", f"{cx + (x - cx) * scale:.2f}")
            text.set("y", f"{cy + (y - cy) * scale:.2f}")

    return ET.tostring(root, encoding="unicode")


def _edge_legend(svg_text):
    """Which of STRING-db's edge evidence-channel colours actually appear
    in this network, in EDGE_COLOR_LEGEND's order — only those get shown,
    since most networks don't use all seven channels."""
    return [(color, label) for color, label in EDGE_COLOR_LEGEND if color in svg_text]


def build_svg_panel(
    div_id,
    title,
    svg_path,
    direction,
    stats_div_id=None,
    taxon_id=None,
    rank_by_gene=None,
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

    rank_by_gene: the paired scatter plot's rank_by_gene (see build_plot),
    used to size each node by its rank among the genes shown here.
    """
    if not svg_path or not os.path.exists(svg_path):
        return None
    try:
        with open(svg_path) as fh:
            svg = fh.read()
    except OSError as exc:  # pragma: no cover - defensive
        logging.warning(f"Could not read STRING-db SVG {svg_path}: {exc}")
        return None
    svg = _resize_nodes_by_rank(svg, rank_by_gene or {})
    return {
        "div_id": div_id,
        "panel_title": title,
        "kind": "svg",
        "svg": svg,
        "direction": direction,
        "stats_div_id": stats_div_id,
        "taxon_id": taxon_id,
        "edge_legend": _edge_legend(svg),
    }
