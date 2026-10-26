"""CRISPRcleanR QC diagnostics (ROC, precision-recall, depletion-rank), all
three grouped under one "CRISPRcleanR" tab, each switchable by comparison
via its own in-figure dropdown — built from the curve/rank data
crisprcleaner.R exports per comparison, not pre-rendered PDFs (see
workflow/scripts/crisprcleaner.R; the statistics themselves — ROC/PR
curves, the 5% FDR rank threshold — are computed there by the
CRISPRcleanR/pROC packages, not re-derived here).

build() returns one section, only when the crisprcleanr rule actually ran
for at least one comparison (see general_functions.smk's
crisprcleanr_comparisons())."""

import logging

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .common import style_fig

GENE_COLOR = "#6366f1"  # indigo-500
SGRNA_COLOR = "#f59e0b"  # amber-500
DIAGONAL_COLOR = "#94a3b8"  # slate-400
MAIN_COLOR = "#334155"  # slate-700
HIT_COLOR = "#2563eb"  # blue-600
MISS_COLOR = "#cbd5e1"  # slate-300
THRESHOLD_COLOR = "#dc2626"  # red-600

SIGNATURE_ORDER = [
    "Ribosomal_Proteins",
    "DNA_Replication",
    "RNA_polymerase",
    "Proteasome",
    "Spliceosome",
    "Common_essential",
    "Non_essential",
]

CITATION = "Iorio et al. 2018, Genome Biology, doi:10.1186/s13059-018-1459-6"


def _comparison_dropdown(comparisons, n_traces_per_comparison):
    """A dropdown (Plotly updatemenu) toggling trace visibility between
    comparisons — every comparison's traces are pre-built and all but the
    first are hidden; switching just restyles "visible", no server
    round-trip."""
    total = len(comparisons) * n_traces_per_comparison
    buttons = []
    for i, comparison in enumerate(comparisons):
        visible = [False] * total
        start = i * n_traces_per_comparison
        for j in range(start, start + n_traces_per_comparison):
            visible[j] = True
        buttons.append(
            dict(label=comparison, method="restyle", args=[{"visible": visible}])
        )
    return dict(
        type="dropdown",
        direction="down",
        buttons=buttons,
        x=0,
        xanchor="left",
        y=1.15,
        yanchor="top",
        showactive=True,
        bgcolor="#ffffff",
        bordercolor="#cbd5e1",
        font=dict(size=13),
    )


def _hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))


def _gradient(base_hex, n):
    """n shades of base_hex from light to full colour, one per comparison
    (comparisons are usually successive timepoints) — lets every curve be
    plotted on one axis at once instead of behind a dropdown, so an
    improving/worsening AUC across timepoints reads directly off the
    colour progression."""
    if n == 1:
        return [base_hex]
    r, g, b = _hex_to_rgb(base_hex)
    colors = []
    for i in range(n):
        t = i / (n - 1)  # 0 = lightest (first comparison), 1 = full colour (last)
        frac_white = 0.75 * (1 - t)
        rr = round(r + (255 - r) * frac_white)
        gg = round(g + (255 - g) * frac_white)
        bb = round(b + (255 - b) * frac_white)
        colors.append(f"#{rr:02x}{gg:02x}{bb:02x}")
    return colors


def _build_curve_fig(comparisons, dfs, x_col, y_col, x_title, y_title, show_diagonal):
    """Shared by ROC (x=1-specificity, y=sensitivity) and PR (x=recall,
    y=precision): two lines per comparison (gene- and sgRNA-level), all
    plotted together — gene lines shaded light-to-dark indigo and sgRNA
    lines light-to-dark amber, in comparison order, so multiple comparisons
    (e.g. timepoints) are directly comparable on one axis rather than
    switched between."""
    gene_colors = _gradient(GENE_COLOR, len(comparisons))
    sgrna_colors = _gradient(SGRNA_COLOR, len(comparisons))

    traces = []
    for i, (comparison, df) in enumerate(zip(comparisons, dfs)):
        for cls, color in (("gene", gene_colors[i]), ("sgRNA", sgrna_colors[i])):
            sub = df[df["class"] == cls].sort_values(x_col)
            auc = sub["auc"].iloc[0] if len(sub) else float("nan")
            traces.append(
                go.Scatter(
                    x=sub[x_col],
                    y=sub[y_col],
                    mode="lines",
                    name=f"{comparison} — {cls} (AUC={auc:.3f})",
                    line=dict(color=color, width=2),
                    hovertemplate=(
                        f"<b>{comparison} — {cls}</b><br>{x_title}: "
                        "%{x:.3f}<br>"
                        f"{y_title}: "
                        "%{y:.3f}<extra></extra>"
                    ),
                )
            )

    fig = go.Figure(data=traces)
    if show_diagonal:
        fig.add_trace(
            go.Scatter(
                x=[0, 1],
                y=[0, 1],
                mode="lines",
                line=dict(color=DIAGONAL_COLOR, width=1, dash="dash"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    fig.update_xaxes(title=x_title, range=[0, 1])
    fig.update_yaxes(title=y_title, range=[0, 1])
    style_fig(
        fig,
        height=550,
        showlegend=True,
        margin=dict(t=60, b=80, l=70, r=30),
    )
    return fig, pd.concat(dfs, ignore_index=True)


def build_roc(comparisons, csv_paths):
    logging.info(f"Building CRISPRcleanR ROC plot for {len(comparisons)} comparison(s)")
    dfs = []
    for comparison, path in zip(comparisons, csv_paths):
        df = pd.read_csv(path)
        # The exported "specificity" column is the true specificity; the
        # plot wants false-positive rate (1 - specificity) on x.
        df["fpr"] = 1 - df["specificity"]
        df["comparison"] = comparison
        dfs.append(df)
    fig, df = _build_curve_fig(
        comparisons,
        dfs,
        x_col="fpr",
        y_col="sensitivity",
        x_title="1 - Specificity",
        y_title="Sensitivity (Recall)",
        show_diagonal=True,
    )
    df = df.drop(columns="fpr")
    return dict(
        id="crisprcleanr-roc",
        title="CRISPRcleanR ROC",
        description=(
            "ROC (Receiver Operating Characteristic) curve: at every "
            "possible corrected log fold-change cutoff, plots sensitivity "
            "(the fraction of essential genes correctly called depleted) "
            "against 1 - specificity (the fraction of non-essential genes "
            "incorrectly called depleted). The closer the curve bows "
            "toward the top-left (AUC near 1), the better essential and "
            "non-essential genes are separated. Unlike the PR curve below, "
            "ROC doesn't depend on how rare essential genes are among all "
            f"genes screened ({CITATION})."
        ),
        fig=fig,
        df=df,
        filename="crisprcleanr-roc",
    )


def build_pr(comparisons, csv_paths):
    logging.info(f"Building CRISPRcleanR PR plot for {len(comparisons)} comparison(s)")
    dfs = []
    for comparison, path in zip(comparisons, csv_paths):
        df = pd.read_csv(path)
        df["comparison"] = comparison
        dfs.append(df)
    fig, df = _build_curve_fig(
        comparisons,
        dfs,
        x_col="recall",
        y_col="precision",
        x_title="Recall",
        y_title="Precision",
        show_diagonal=False,
    )
    return dict(
        id="crisprcleanr-pr",
        title="CRISPRcleanR PR",
        description=(
            "Precision-recall (PR) curve: at the same cutoffs as the ROC "
            "curve above, plots precision (of the genes called depleted, "
            "the fraction that are truly essential) against recall (the "
            "fraction of all essential genes successfully called "
            "depleted). Essential genes are a small minority of all genes "
            "screened, and PR — unlike ROC — is sensitive to that "
            "imbalance, so it penalises false positives more visibly and "
            "tends to look less optimistic than ROC for the same "
            f"classifier ({CITATION})."
        ),
        fig=fig,
        df=df,
        filename="crisprcleanr-pr",
    )


def build_depletion_rank(comparisons, gene_csv_paths, signature_csv_paths):
    logging.info(
        f"Building CRISPRcleanR depletion-rank plot for {len(comparisons)} comparison(s)"
    )
    n_cols = 1 + len(SIGNATURE_ORDER)
    fig = make_subplots(
        rows=1,
        cols=n_cols,
        shared_yaxes=True,
        column_widths=[0.5] + [0.5 / len(SIGNATURE_ORDER)] * len(SIGNATURE_ORDER),
        horizontal_spacing=0.008,
        subplot_titles=["Depletion profile"]
        + [s.replace("_", " ") for s in SIGNATURE_ORDER],
    )
    base_annotations = [a.to_plotly_json() for a in fig.layout.annotations]
    # The narrow signature columns can't fit a horizontal title without
    # overlapping their neighbours — angle them steeply (near-vertical, like
    # the original R plot's srt=80 rotated column labels in
    # ccr.VisDepAndSig) so each one's horizontal footprint stays within its
    # own column instead of spilling into the next.
    for annotation in base_annotations[1:]:
        annotation.update(
            textangle=-80, xanchor="left", yanchor="bottom", font=dict(size=13)
        )
    # Centre x (paper fraction) of each signature column, for the hit-%
    # labels placed below the plot — make_subplots already fixed every
    # column's domain via column_widths, so these don't change per
    # comparison.
    col_centers = {
        col: sum(fig.layout[f"xaxis{col}"].domain) / 2 for col in range(2, n_cols + 1)
    }

    traces_per_comparison = None
    comparison_annotations = []
    comparison_shapes = []
    all_gene_dfs = []
    for comparison, gpath, spath in zip(
        comparisons, gene_csv_paths, signature_csv_paths
    ):
        genes = pd.read_csv(gpath).sort_values("fc").reset_index(drop=True)
        genes["rank"] = genes.index + 1
        genes["comparison"] = comparison
        all_gene_dfs.append(genes)
        fdr_rank = genes["fdr5perc_rank"].iloc[0]
        rank_by_gene = dict(zip(genes["gene"], genes["rank"]))
        sigs = pd.read_csv(spath)

        is_first = comparison == comparisons[0]
        start = len(fig.data)

        fig.add_trace(
            go.Scattergl(
                x=genes["fc"],
                y=genes["rank"],
                mode="markers",
                marker=dict(size=4, color=MAIN_COLOR),
                showlegend=False,
                visible=is_first,
                customdata=genes["gene"],
                hovertemplate=(
                    "<b>%{customdata}</b><br>LFC: %{x:.3f}<br>Rank: %{y}<extra></extra>"
                ),
            ),
            row=1,
            col=1,
        )
        # A layout shape (not a trace) spanning the full figure width via
        # xref="paper", so the threshold line crosses every signature
        # column too, not just the main depletion-profile panel — yref="y"
        # still lines it up correctly everywhere since shared_yaxes=True
        # (see make_subplots call above) matches every column's y-axis
        # range/type to this one.
        comparison_shapes.append(
            [
                dict(
                    type="line",
                    xref="paper",
                    yref="y",
                    x0=0,
                    x1=1,
                    y0=fdr_rank,
                    y1=fdr_rank,
                    line=dict(color=THRESHOLD_COLOR, width=2, dash="dash"),
                )
            ]
        )

        comparison_annotations.append([dict(a) for a in base_annotations])
        for j, sig in enumerate(SIGNATURE_ORDER):
            col = j + 2
            sig_genes = sigs.loc[sigs["signature"] == sig, "gene"]
            ranks = sorted(rank_by_gene[g] for g in sig_genes if g in rank_by_gene)
            hits = [r for r in ranks if r <= fdr_rank]
            misses = [r for r in ranks if r > fdr_rank]
            for group, color in ((hits, HIT_COLOR), (misses, MISS_COLOR)):
                xs, ys = [], []
                for r in group:
                    xs += [0, 1, None]
                    ys += [r, r, None]
                fig.add_trace(
                    go.Scatter(
                        x=xs,
                        y=ys,
                        mode="lines",
                        line=dict(color=color, width=1.5),
                        showlegend=False,
                        visible=is_first,
                        hoverinfo="skip",
                    ),
                    row=1,
                    col=col,
                )
            pct = round(100 * len(hits) / len(ranks)) if ranks else 0
            comparison_annotations[-1].append(
                dict(
                    text=f"{pct}%",
                    x=col_centers[col],
                    xref="paper",
                    xanchor="center",
                    y=-0.04,
                    yref="paper",
                    yanchor="top",
                    showarrow=False,
                    font=dict(size=13, color=HIT_COLOR),
                )
            )

        end = len(fig.data)
        if traces_per_comparison is None:
            traces_per_comparison = end - start

    fig.update_yaxes(type="log", autorange="reversed", title="Depletion rank", col=1)
    # shared_yaxes=True only syncs the numeric range across columns (via
    # matches="y"), not type/autorange — every signature column needs its
    # own type="log"/reversed, or its ticks render on a mismatched linear
    # scale, visually misplacing which genes look "above" vs "below" the
    # threshold line relative to the main panel.
    for col in range(2, n_cols + 1):
        fig.update_yaxes(
            type="log", autorange="reversed", showticklabels=False, col=col
        )
        fig.update_xaxes(visible=False, range=[0, 1], col=col)
    fig.update_xaxes(title="Log2(Fold change)", col=1)

    updatemenus = []
    if len(comparisons) > 1:
        menu = _comparison_dropdown(comparisons, traces_per_comparison)
        for button, annotations, shapes in zip(
            menu["buttons"], comparison_annotations, comparison_shapes
        ):
            button["method"] = "update"
            button["args"] = [
                button["args"][0],
                {"annotations": annotations, "shapes": shapes},
            ]
        updatemenus = [menu]
    # The dropdown buttons' "update" args only apply on click — the figure's
    # own initial state (matching the first comparison's traces, the ones
    # left visible above) needs the same annotations/shapes set directly,
    # or the default view renders with make_subplots' unrotated titles and
    # no threshold line until a comparison is explicitly picked.
    fig.layout.annotations = comparison_annotations[0]
    fig.layout.shapes = comparison_shapes[0]

    style_fig(
        fig,
        height=850,
        showlegend=False,
        updatemenus=updatemenus,
        margin=dict(t=220 if updatemenus else 190, b=90, l=70, r=20),
    )

    return dict(
        id="crisprcleanr-depletion-rank",
        title="CRISPRcleanR depletion rank",
        description=(
            "Gene-level depletion rank (by corrected log fold-change) against "
            "known essential gene-set signatures. Blue ticks mark genes within "
            "the 5% FDR threshold (dashed red line); the percentage per "
            "signature is the share of that signature's genes meeting it "
            f"({CITATION})."
        ),
        fig=fig,
        df=pd.concat(all_gene_dfs, ignore_index=True),
        filename="crisprcleanr-depletion-rank",
        # No legend, but several extra gene-signature columns each with a
        # rotated header — those overlap/overflow at kaleido's narrower
        # default PDF export width (see report_html.py's _encode_pdf).
        wide_pdf=True,
    )


def build(
    comparisons,
    roc_csv_paths,
    pr_csv_paths,
    drnk_gene_csv_paths,
    drnk_signature_csv_paths,
):
    """All three plots under one "CRISPRcleanR" tab (see qc/report_html.py's
    "multi"-kind section handling) rather than three separate top-level
    tabs."""
    return dict(
        id="crisprcleanr",
        title="CRISPRcleanR",
        kind="multi",
        subsections=[
            build_roc(comparisons, roc_csv_paths),
            build_pr(comparisons, pr_csv_paths),
            build_depletion_rank(
                comparisons, drnk_gene_csv_paths, drnk_signature_csv_paths
            ),
        ],
    )
