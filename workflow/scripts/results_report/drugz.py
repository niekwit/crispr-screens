"""Per-comparison DrugZ depleted/enriched scatter entries."""

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
    string_enriched_id_mapping_files,
    string_depleted_id_mapping_files,
    comparisons,
):
    entries = []
    for i, comparison in enumerate(comparisons):
        txt = txt_files[i]
        logging.info(f"Loading {txt}")
        data = pd.read_csv(txt, sep="\t")

        # DrugZ reports a single signed normZ per gene: negative normZ are
        # synthetic-lethal/depleted interactions (significance in
        # fdr_synth), positive normZ are suppressor/enriched interactions
        # (significance in fdr_supp) — same split plot_drugz_results.R uses.
        depleted = data[data["normZ"] < 0].copy()
        depleted["x"] = np.random.RandomState(42).permutation(len(depleted))
        depleted["log_fdr"] = -np.log10(depleted["fdr_synth"].clip(lower=1e-300))
        depleted["_normz_str"] = depleted["normZ"].round(4).astype(str)
        depleted["_logfdr_str"] = depleted["log_fdr"].round(3).astype(str)
        depleted["_numobs_str"] = depleted["numObs"].astype(str)

        enriched = data[data["normZ"] > 0].copy()
        enriched["x"] = np.random.RandomState(42).permutation(len(enriched))
        enriched["log_fdr"] = -np.log10(enriched["fdr_supp"].clip(lower=1e-300))
        enriched["_normz_str"] = enriched["normZ"].round(4).astype(str)
        enriched["_logfdr_str"] = enriched["log_fdr"].round(3).astype(str)
        enriched["_numobs_str"] = enriched["numObs"].astype(str)

        logging.info(
            f"Depleted: {len(depleted)} genes, Enriched: {len(enriched)} genes"
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
        enriched_id_mapping = (
            common.load_id_mapping(string_enriched_id_mapping_files[i])
            if string_enriched_id_mapping_files
            else {}
        )
        depleted_id_mapping = (
            common.load_id_mapping(string_depleted_id_mapping_files[i])
            if string_depleted_id_mapping_files
            else {}
        )

        key = common.slugify(f"drugz__{comparison}")
        hovertemplate = (
            "<b>%{customdata[0]}</b><br>"
            "normZ: %{customdata[1]}<br>"
            "-log10(FDR): %{customdata[2]}<br>"
            "numObs: %{customdata[3]}"
        )
        customdata_cols = ["GENE", "_normz_str", "_logfdr_str", "_numobs_str"]

        plots = [
            common.build_plot(
                div_id=f"{key}-depleted",
                title=f"Depleted (synthetic-lethal) genes — {comparison}",
                term_panel_title="STRING-db depleted terms",
                term_list_id=f"{key}-depleted-terms",
                df=depleted,
                y_col="normZ",
                y_title="normZ",
                color_col="log_fdr",
                color_title="-log10(FDR)",
                id_col="GENE",
                customdata_cols=customdata_cols,
                hovertemplate=hovertemplate,
                terms=depleted_terms,
                direction="depleted",
                id_mapping=depleted_id_mapping,
            ),
            common.build_svg_panel(
                div_id=f"{key}-depleted-svg",
                title="STRING-db depleted network",
                svg_path=(
                    string_depleted_svg_files[i] if string_depleted_svg_files else None
                ),
                direction="depleted",
                stats_div_id=f"{key}-depleted",
            ),
            common.build_plot(
                div_id=f"{key}-enriched",
                title=f"Enriched (suppressor) genes — {comparison}",
                term_panel_title="STRING-db enriched terms",
                term_list_id=f"{key}-enriched-terms",
                df=enriched,
                y_col="normZ",
                y_title="normZ",
                color_col="log_fdr",
                color_title="-log10(FDR)",
                id_col="GENE",
                customdata_cols=customdata_cols,
                hovertemplate=hovertemplate,
                terms=enriched_terms,
                direction="enriched",
                id_mapping=enriched_id_mapping,
            ),
            common.build_svg_panel(
                div_id=f"{key}-enriched-svg",
                title="STRING-db enriched network",
                svg_path=(
                    string_enriched_svg_files[i] if string_enriched_svg_files else None
                ),
                direction="enriched",
                stats_div_id=f"{key}-enriched",
            ),
        ]
        plots = [p for p in plots if p is not None]

        entries.append(
            dict(
                key=key,
                label=comparison,
                plots=plots,
            )
        )
    return entries
