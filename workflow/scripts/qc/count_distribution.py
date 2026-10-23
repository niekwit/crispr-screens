"""Read count distribution of sgRNAs per sample, faceted as a subplot grid."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .common import (
    BAR_COLOR,
    BAR_LINE_COLOR,
    CUTOFF_COLOR,
    FLAG_COLOR,
    MEDIAN_COLOR,
    style_fig,
)

LOW_FRACTION = 0.1  # sgRNAs below this fraction of the sample mean are "low"
FLAG_THRESHOLD = 0.1  # samples with more low sgRNAs than this are flagged


def build(counts, sample_cols, samples, output_csv, sorted_samples=frozenset()):
    dist_rows = []
    log_counts = {}
    for col, sample in zip(sample_cols, samples):
        values = counts[col].values.astype(float)
        mean = values.mean()
        median = float(np.median(values))
        dist_rows.append(
            {
                "sample": sample,
                "mean": mean,
                "median": median,
                "median_over_mean": median / mean if mean else float("nan"),
                "fraction_low": float((values < LOW_FRACTION * mean).mean()),
                "fraction_zero": float((values == 0).mean()),
                "sorted_sample": sample in sorted_samples,
            }
        )
        log_counts[sample] = np.log10(values + 1)
    df = pd.DataFrame(dist_rows)
    # A FACS-sorted sample is expected to have a bottlenecked sgRNA
    # distribution from the sort itself, not from a technical problem, so
    # it's exempted from the dropout flag (see sorted_samples() in
    # general_functions.smk) even if its fraction_low would otherwise
    # cross FLAG_THRESHOLD.
    df["flagged"] = (df["fraction_low"] > FLAG_THRESHOLD) & ~df["sorted_sample"]
    df.to_csv(output_csv, index=False)

    ncol_facet = min(4, len(samples))
    nrow_facet = -(-len(samples) // ncol_facet)  # ceiling division
    fig = make_subplots(
        rows=nrow_facet,
        cols=ncol_facet,
        subplot_titles=[
            f"<b>{row['sample']}</b>"
            + (" (sorted)" if row["sorted_sample"] else "")
            + f"<br>{100 * row['fraction_low']:.0f}% low<br>"
            f"median/mean {row['median_over_mean']:.2f}"
            for _, row in df.iterrows()
        ],
        vertical_spacing=max(0.5 / nrow_facet, 0.12),
    )
    for i, (_, row) in enumerate(df.iterrows()):
        r, c = i // ncol_facet + 1, i % ncol_facet + 1
        sample = row["sample"]
        low_cutoff_log = np.log10(LOW_FRACTION * row["mean"] + 1)
        median_log = np.log10(row["median"] + 1)

        if row["flagged"]:
            axis_n = "" if i == 0 else str(i + 1)
            fig.add_shape(
                type="rect",
                xref=f"x{axis_n} domain",
                yref=f"y{axis_n} domain",
                x0=0,
                x1=1,
                y0=0,
                y1=1,
                fillcolor=FLAG_COLOR,
                opacity=0.12,
                line_width=0,
                layer="below",
                row=r,
                col=c,
            )

        fig.add_trace(
            go.Histogram(
                x=log_counts[sample],
                nbinsx=60,
                marker=dict(
                    color=BAR_COLOR, line=dict(color=BAR_LINE_COLOR, width=0.5)
                ),
                showlegend=False,
                hovertemplate="log10(count+1): %{x:.2f}<br>sgRNAs: %{y}<extra></extra>",
            ),
            row=r,
            col=c,
        )
        fig.add_vline(
            x=low_cutoff_log, line_dash="dash", line_color=CUTOFF_COLOR, row=r, col=c
        )
        fig.add_vline(x=median_log, line_color=MEDIAN_COLOR, line_width=2, row=r, col=c)

    fig.update_annotations(font_size=12)
    fig.update_xaxes(title_text="log10(count + 1)", title_font_size=13)
    fig.update_yaxes(title_text="sgRNAs", title_font_size=13)
    style_fig(fig, height=320 * nrow_facet + 80, margin=dict(t=60, b=80, l=60, r=30))

    return dict(
        id="count-distribution",
        title="Read count distribution",
        description=(
            "Read count distribution of sgRNAs per sample. A single peak indicates "
            "a healthy sample. A second peak at low counts (or a large fraction of "
            "sgRNAs below 10% of the mean count) indicates sgRNA dropout, typically "
            "caused by a bottleneck when too few cells were sorted or collected. In "
            "such samples MAGeCK gene-level LFCs can show discrete bands, because "
            "the median of the sgRNA LFCs depends mainly on how many sgRNAs per "
            "gene were lost by chance. Blue line: median. Dashed line: low-count "
            "cut-off. Shaded panel: sample flagged for sgRNA dropout. Samples "
            "marked (sorted) in config/samples.csv are exempt from the flag, "
            "since a bottlenecked distribution is expected from the sort itself."
        ),
        fig=fig,
        df=df,
        filename="count-distribution",
    )
