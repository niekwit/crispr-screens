"""Per-sample overview table combining the flags raised by the other QC
metrics, so problem samples are visible without checking every tab."""

import math

import plotly.graph_objects as go

from .alignment import MIN_ALIGNMENT_RATE
from .common import style_fig
from .count_distribution import FLAG_THRESHOLD, LOW_FRACTION

FLAG_ROW_COLOR = "#fff1f2"  # rose-50
OK_ROW_COLOR = "#ffffff"
HEADER_COLOR = "#0f172a"

CELL_FONT_SIZE = 14
# go.Table doesn't expose the wrapped line count it actually renders, and a
# fixed per-row cell height (the old 32px) only fits a single line — any row
# whose text wraps to a second line (e.g. a long "Flags" value) still only
# gets 32px, so go.Table silently clips it instead of wrapping the figure's
# height to fit: rows run off the bottom of the fixed-height SVG canvas with
# no scrollbar to reveal them. Estimating the wrapped line count per cell
# (from the column's pixel width and an average character width) and sizing
# the header/cell/figure heights to the worst case avoids that clipping.
AVG_CHAR_WIDTH_FACTOR = 0.56  # of CELL_FONT_SIZE; Arial-ish average glyph width
CELL_LINE_HEIGHT = 18  # px per wrapped line at CELL_FONT_SIZE
CELL_PADDING = 16  # go.Table's default horizontal cell padding (both sides)


def _wrapped_lines(text, column_width):
    """How many lines `text` wraps to inside a go.Table column of this
    pixel width, at CELL_FONT_SIZE."""
    chars_per_line = max(
        1, int((column_width - CELL_PADDING) / (CELL_FONT_SIZE * AVG_CHAR_WIDTH_FACTOR))
    )
    return max(1, math.ceil(len(str(text)) / chars_per_line))


def build(
    alignment_section,
    coverage_section,
    gini_section,
    missed_section,
    count_dist_section,
):
    align_df = alignment_section["df"].rename(columns={"flagged": "low_alignment"})
    cov_df = coverage_section["df"]
    gini_df = gini_section["df"]
    missed_df = missed_section["df"]
    dist_df = count_dist_section["df"][
        ["sample", "fraction_low", "flagged", "sorted_sample"]
    ].rename(
        columns={"fraction_low": "fraction_low_sgrnas", "flagged": "sgrna_dropout"}
    )

    summary = (
        align_df[["sample", "mapping_rate", "low_alignment"]]
        .merge(cov_df[["sample", "coverage"]], on="sample")
        .merge(gini_df[["sample", "gini_index"]], on="sample")
        .merge(missed_df[["sample", "missed_sgrnas"]], on="sample")
        .merge(dist_df, on="sample")
    )
    summary["flagged"] = summary["low_alignment"] | summary["sgrna_dropout"]

    def flag_text(row):
        flags = []
        if row["low_alignment"]:
            flags.append("Low alignment rate")
        if row["sgrna_dropout"]:
            flags.append("sgRNA dropout")
        return ", ".join(flags) if flags else "—"

    summary["flags"] = summary.apply(flag_text, axis=1)

    header = [
        "Sample",
        "Sorted",
        "Mapping rate (%)",
        "Coverage (x)",
        "Gini index",
        "Missed sgRNAs",
        "Low-count sgRNAs (%)",
        "Flags",
    ]
    cell_values = [
        summary["sample"],
        summary["sorted_sample"].map({True: "Yes", False: ""}),
        summary["mapping_rate"].round(1),
        summary["coverage"].round(1),
        summary["gini_index"].round(3),
        summary["missed_sgrnas"],
        (100 * summary["fraction_low_sgrnas"]).round(1),
        summary["flags"],
    ]
    row_colors = [
        FLAG_ROW_COLOR if flagged else OK_ROW_COLOR for flagged in summary["flagged"]
    ]

    # Sample column sized to the longest actual sample name (a fixed 90
    # truncated/clipped longer ones); the rest widened so each header title
    # fits on one line instead of wrapping (go.Table wraps header text that
    # doesn't fit the column's width, there's no "nowrap" option).
    #
    # go.Table's columnwidth values are *relative* weights, rescaled to
    # fill whatever pixel width the table actually renders at — not
    # absolute pixels. Since the live report's .plot div is responsive
    # (width: 100% of its ~1100px card), just raising some columns'
    # weights while leaving others unchanged redistributes the *same*
    # total width, shrinking the untouched columns (e.g. "Sorted") instead
    # of actually growing the ones that needed it. So this table also gets
    # an explicit figure width (plot_width, below) equal to the sum of
    # these, so the values below really do land as pixels.
    sample_col_width = max(100, 10 * summary["sample"].astype(str).str.len().max() + 40)
    columnwidth = [sample_col_width, 90, 210, 160, 140, 170, 250, 260]
    table_width = sum(columnwidth)

    # "Flags" is the only column whose values can run long enough to wrap
    # (e.g. "Low alignment rate, sgRNA dropout") — go.Table has a single
    # fixed cell height for the whole table, so any row that needs a second
    # line must get it, or go.Table clips that row's text (and any rows
    # pushed past the figure's fixed total height) instead of wrapping to
    # fit. Size cell_height/table height for the worst case across rows.
    max_flag_lines = max(
        (_wrapped_lines(flags, columnwidth[-1]) for flags in summary["flags"]),
        default=1,
    )
    cell_height = (
        32 if max_flag_lines <= 1 else 32 + (max_flag_lines - 1) * CELL_LINE_HEIGHT
    )

    fig = go.Figure(
        data=[
            go.Table(
                columnwidth=columnwidth,
                header=dict(
                    values=header,
                    fill_color=HEADER_COLOR,
                    font=dict(color="white", size=CELL_FONT_SIZE),
                    align="left",
                    height=36,
                ),
                cells=dict(
                    values=cell_values,
                    fill_color=[row_colors] * len(header),
                    align="left",
                    height=cell_height,
                    font=dict(size=CELL_FONT_SIZE),
                ),
            )
        ]
    )
    style_fig(
        fig,
        height=110 + (cell_height + 6) * len(summary),
        margin=dict(t=20, b=20, l=20, r=20),
        width=table_width,
    )

    return dict(
        id="summary",
        title="Overview",
        # Fixes the .plot div's on-screen width to table_width (scrolling
        # horizontally if that's wider than its card) instead of 100% of
        # the card — see the columnwidth comment above for why.
        plot_width=table_width,
        description=(
            "Per-sample summary across all QC metrics. Rows highlighted in rose "
            f"are flagged for a low alignment rate (<{MIN_ALIGNMENT_RATE}%) or "
            f"sgRNA dropout (>{100 * FLAG_THRESHOLD:.0f}% of sgRNAs below "
            f"{100 * LOW_FRACTION:.0f}% of the sample mean count) — samples "
            "marked sorted in config/samples.csv are exempt from the dropout "
            "flag."
        ),
        fig=fig,
        df=summary.drop(columns=["flagged"]),
        filename="qc-summary",
        # Wider columns (above) need the wider PDF canvas too, or the same
        # header-wrapping this fixes on-screen just reappears in export.
        wide_pdf=True,
    )
