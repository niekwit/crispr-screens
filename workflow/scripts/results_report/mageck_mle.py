"""Per-(matrix, cnv) MAGeCK MLE scatter entries.

Unlike MAGeCK test (one enriched/depleted pair per comparison), MLE reports
one beta/z/p-value/fdr/wald-p-value/wald-fdr block per non-baseline column
of the design matrix — the "direction" tabs here are per-condition rather
than enriched/depleted. There's no STRING-db pathway analysis for MLE (no
rule computes it), so these are plain scatter plots: no STRING term panel,
no network SVG, but gene search/click-to-label and the STRING-db gene-info
popup on click still work (see common.build_plot/report_html.py).
"""

import logging

import numpy as np
import pandas as pd

from . import common


def _conditions(columns):
    """Condition names are the non-baseline design-matrix columns, detected
    from the "|beta" columns MAGeCK MLE reports one of per condition."""
    return [c[: -len("|beta")] for c in columns if c.endswith("|beta")]


def _color_options(df, cond):
    log_pval = -np.log10(df[f"{cond}|p-value"].clip(lower=1e-300))
    log_fdr = -np.log10(df[f"{cond}|fdr"].clip(lower=1e-300))
    log_wald_fdr = -np.log10(df[f"{cond}|wald-fdr"].clip(lower=1e-300))
    return [
        ("neg_log10_fdr", "-log10(FDR)", log_fdr),
        ("neg_log10_pvalue", "-log10(p-value)", log_pval),
        ("neg_log10_wald_fdr", "-log10(Wald FDR)", log_wald_fdr),
        ("beta", "Beta score", df[f"{cond}|beta"]),
    ]


def build_entries(txt_files, pairs, taxon_id=None):
    """pairs: list of (matrix, cnv) tuples, positionally matching
    txt_files."""
    entries = []
    for i, (matrix, cnv) in enumerate(pairs):
        txt = txt_files[i]
        logging.info(f"Loading {txt}")
        data = pd.read_csv(txt, sep="\t")
        np.random.seed(42)

        key = common.slugify(f"mageck_mle__{cnv}__{matrix}")
        conditions = _conditions(data.columns)
        logging.info(f"{matrix} ({cnv}): conditions {conditions}")

        plots = []
        for cond in conditions:
            df = data.copy()
            df["x"] = np.random.permutation(len(df))
            df["_beta_str"] = df[f"{cond}|beta"].round(4).astype(str)
            df["_z_str"] = df[f"{cond}|z"].round(3).astype(str)
            df["_fdr_str"] = df[f"{cond}|fdr"].round(6).astype(str)
            df["_wald_fdr_str"] = df[f"{cond}|wald-fdr"].round(6).astype(str)

            hovertemplate = (
                "<b>%{customdata[0]}</b><br>"
                "Beta: %{customdata[1]}<br>"
                "Z: %{customdata[2]}<br>"
                "FDR: %{customdata[3]}<br>"
                "Wald FDR: %{customdata[4]}"
            )
            customdata_cols = [
                "Gene",
                "_beta_str",
                "_z_str",
                "_fdr_str",
                "_wald_fdr_str",
            ]

            plots.append(
                common.build_plot(
                    div_id=common.slugify(f"{key}-{cond}"),
                    title=f"MAGeCK MLE — {matrix} ({cnv}) — {cond}",
                    term_panel_title="STRING-db terms",
                    term_list_id=common.slugify(f"{key}-{cond}-terms"),
                    df=df,
                    y_col=f"{cond}|beta",
                    y_title="Beta score",
                    color_options=_color_options(df, cond),
                    id_col="Gene",
                    customdata_cols=customdata_cols,
                    hovertemplate=hovertemplate,
                    terms=[],
                    direction=cond,
                    taxon_id=taxon_id,
                )
            )

        entries.append(
            dict(
                key=key,
                label=f"{matrix} ({cnv})",
                plots=plots,
            )
        )
    return entries
