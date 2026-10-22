"""Shared styling helpers for the QC report figures."""

import re

import plotly.graph_objects as go

# ── Modern colour palette (replaces the old "snakemake green" #419179) ───────
BAR_COLOR = "#6366f1"  # indigo-500
BAR_LINE_COLOR = "#3730a3"  # indigo-800
MEDIAN_COLOR = "#0ea5e9"  # sky-500
CUTOFF_COLOR = "#64748b"  # slate-500
FLAG_COLOR = "#f43f5e"  # rose-500
TEMPLATE = "plotly_white"
FONT = dict(family="Inter, -apple-system, 'Segoe UI', Roboto, sans-serif", size=15)

# Qualitative palette for multi-series charts (one colour per sample).
SAMPLE_PALETTE = [
    "#6366f1",  # indigo-500
    "#0ea5e9",  # sky-500
    "#f43f5e",  # rose-500
    "#f59e0b",  # amber-500
    "#10b981",  # emerald-500
    "#8b5cf6",  # violet-500
    "#ec4899",  # pink-500
    "#14b8a6",  # teal-500
]


def clean_sample_name(name):
    # remove prepended X from sample names (only happens when they start
    # with a number)
    return re.sub(r"^X", "", name)


def sample_color(i):
    return SAMPLE_PALETTE[i % len(SAMPLE_PALETTE)]


def style_fig(fig, height=450, margin=None, showlegend=False, **layout_kwargs):
    fig.update_layout(
        template=TEMPLATE,
        font=FONT,
        height=height,
        margin=margin or dict(t=60, b=80, l=60, r=30),
        showlegend=showlegend,
        **layout_kwargs,
    )
    return fig


def _sort_menu(categories, values, bar_color, bar_line_color):
    """A dropdown (Plotly updatemenu) that restyles the bar order in place,
    between the original sample order and ascending/descending by value."""
    is_list_color = isinstance(bar_color, list)
    orig_idx = list(range(len(categories)))
    orders = {
        "Sample order": orig_idx,
        "Value: high → low": sorted(orig_idx, key=lambda i: values[i], reverse=True),
        "Value: low → high": sorted(orig_idx, key=lambda i: values[i]),
    }

    def restyle_args(idx):
        update = {
            "x": [[categories[i] for i in idx]],
            "y": [[values[i] for i in idx]],
        }
        if is_list_color:
            update["marker.color"] = [[bar_color[i] for i in idx]]
            update["marker.line.color"] = [[bar_line_color[i] for i in idx]]
        return update

    return dict(
        type="dropdown",
        direction="down",
        buttons=[
            dict(label=label, method="restyle", args=[restyle_args(idx)])
            for label, idx in orders.items()
        ],
        x=0,
        xanchor="left",
        y=1.2,
        yanchor="top",
        showactive=True,
        bgcolor="#ffffff",
        bordercolor="#cbd5e1",
        font=dict(size=13),
    )


def bar_fig(df, x, y, y_title, title, y_range=None, flagged=None, sortable=True):
    categories = list(df[x])
    values = list(df[y])
    if flagged is None:
        bar_color = BAR_COLOR
        bar_line_color = BAR_LINE_COLOR
    else:
        flagged = list(flagged)
        bar_color = [FLAG_COLOR if f else BAR_COLOR for f in flagged]
        bar_line_color = [FLAG_COLOR if f else BAR_LINE_COLOR for f in flagged]
    fig = go.Figure(
        go.Bar(
            x=categories,
            y=values,
            marker=dict(color=bar_color, line=dict(color=bar_line_color, width=1)),
            hovertemplate=f"<b>%{{x}}</b><br>{y_title}: %{{y:.3f}}<extra></extra>",
        )
    )
    fig.update_xaxes(tickangle=45, title=None)
    fig.update_yaxes(title=y_title, range=y_range, rangemode="tozero")

    updatemenus = []
    if sortable and len(categories) > 1:
        updatemenus = [_sort_menu(categories, values, bar_color, bar_line_color)]

    return style_fig(
        fig,
        title=dict(text=title, x=0.5),
        margin=dict(t=110, b=80, l=60, r=30),
        updatemenus=updatemenus,
    )
