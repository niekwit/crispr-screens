"""Per-sample overview table combining the flags raised by the other QC
metrics, so problem samples are visible without checking every tab."""

import plotly.graph_objects as go

from .alignment import MIN_ALIGNMENT_RATE
from .common import style_fig
from .count_distribution import FLAG_THRESHOLD, LOW_FRACTION

FLAG_ROW_COLOR = "#fff1f2"  # rose-50
OK_ROW_COLOR = "#ffffff"
HEADER_COLOR = "#0f172a"


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

    fig = go.Figure(
        data=[
            go.Table(
                columnwidth=[90, 70, 110, 90, 80, 100, 130, 220],
                header=dict(
                    values=header,
                    fill_color=HEADER_COLOR,
                    font=dict(color="white", size=14),
                    align="left",
                    height=36,
                ),
                cells=dict(
                    values=cell_values,
                    fill_color=[row_colors] * len(header),
                    align="left",
                    height=32,
                    font=dict(size=14),
                ),
            )
        ]
    )
    style_fig(fig, height=110 + 38 * len(summary), margin=dict(t=20, b=20, l=20, r=20))

    return dict(
        id="summary",
        title="Overview",
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
    )
