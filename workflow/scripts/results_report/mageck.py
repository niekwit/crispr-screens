"""Per-(comparison, cnv) MAGeCK enriched/depleted scatter entries."""

import logging

import numpy as np
import pandas as pd

from . import common


def build_entries(
    txt_files,
    string_enriched_files,
    string_depleted_files,
    string_enriched_svg_files,
    string_depleted_svg_files,
    pairs,
):
    """pairs: list of (comparison, cnv) tuples, positionally matching
    txt_files / string_enriched_files / string_depleted_files (the latter
    two, and their _svg counterparts, are empty lists when STRING-db is
    disabled for that direction)."""
    entries = []
    for i, (comparison, cnv) in enumerate(pairs):
        txt = txt_files[i]
        logging.info(f"Loading {txt}")
        data = pd.read_csv(txt, sep="\t")
        np.random.seed(42)

        enriched = data[data["pos|lfc"] > 0].copy()
        enriched["x"] = np.random.permutation(len(enriched))
        enriched["lfc"] = enriched["pos|lfc"]
        enriched["log_pval"] = -np.log10(enriched["pos|p-value"].clip(lower=1e-300))
        enriched["_lfc_str"] = enriched["lfc"].round(4).astype(str)
        enriched["_logpval_str"] = enriched["log_pval"].round(3).astype(str)
        enriched["_fdr_str"] = enriched["pos|fdr"].round(6).astype(str)

        depleted = data[data["neg|lfc"] < 0].copy()
        depleted["x"] = np.random.permutation(len(depleted))
        depleted["lfc"] = depleted["neg|lfc"]
        depleted["log_pval"] = -np.log10(depleted["neg|p-value"].clip(lower=1e-300))
        depleted["_lfc_str"] = depleted["lfc"].round(4).astype(str)
        depleted["_logpval_str"] = depleted["log_pval"].round(3).astype(str)
        depleted["_fdr_str"] = depleted["neg|fdr"].round(6).astype(str)

        logging.info(
            f"Enriched: {len(enriched)} genes, Depleted: {len(depleted)} genes"
        )

        enriched_terms = (
            common.load_string_terms(string_enriched_files[i])
            if string_enriched_files
            else []
        )
        depleted_terms = (
            common.load_string_terms(string_depleted_files[i])
            if string_depleted_files
            else []
        )

        key = common.slugify(f"mageck__{cnv}__{comparison}")
        hovertemplate = (
            "<b>%{customdata[0]}</b><br>"
            "LFC: %{customdata[1]}<br>"
            "-log10(p-value): %{customdata[2]}<br>"
            "FDR: %{customdata[3]}"
        )
        customdata_cols = ["id", "_lfc_str", "_logpval_str", "_fdr_str"]

        plots = [
            common.build_plot(
                div_id=f"{key}-enriched",
                title=f"Enriched genes — {comparison} ({cnv})",
                term_panel_title="STRING-db enriched terms",
                term_list_id=f"{key}-enriched-terms",
                df=enriched,
                y_col="lfc",
                y_title="Log₂ Fold Change",
                color_col="log_pval",
                color_title="-log10(p-value)",
                id_col="id",
                customdata_cols=customdata_cols,
                hovertemplate=hovertemplate,
                terms=enriched_terms,
                direction="enriched",
            ),
            common.build_svg_panel(
                div_id=f"{key}-enriched-svg",
                title="STRING-db enriched network",
                svg_path=(
                    string_enriched_svg_files[i] if string_enriched_svg_files else None
                ),
                direction="enriched",
            ),
            common.build_plot(
                div_id=f"{key}-depleted",
                title=f"Depleted genes — {comparison} ({cnv})",
                term_panel_title="STRING-db depleted terms",
                term_list_id=f"{key}-depleted-terms",
                df=depleted,
                y_col="lfc",
                y_title="Log₂ Fold Change",
                color_col="log_pval",
                color_title="-log10(p-value)",
                id_col="id",
                customdata_cols=customdata_cols,
                hovertemplate=hovertemplate,
                terms=depleted_terms,
                direction="depleted",
            ),
            common.build_svg_panel(
                div_id=f"{key}-depleted-svg",
                title="STRING-db depleted network",
                svg_path=(
                    string_depleted_svg_files[i] if string_depleted_svg_files else None
                ),
                direction="depleted",
            ),
        ]
        plots = [p for p in plots if p is not None]

        entries.append(
            dict(
                key=key,
                label=f"{comparison} ({cnv})",
                plots=plots,
            )
        )
    return entries
