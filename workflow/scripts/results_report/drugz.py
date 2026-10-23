"""Per-comparison DrugZ depleted/enriched scatter entries."""

import logging

import numpy as np
import pandas as pd

from . import common


def _color_options(df, suffix):
    """DrugZ's output columns for one direction (_synth/_supp), offered as
    "Color by" choices for the dot plot's marker fill."""
    log_fdr = -np.log10(df[f"fdr_{suffix}"].clip(lower=1e-300))
    log_pval = -np.log10(df[f"pval_{suffix}"].clip(lower=1e-300))
    return [
        ("neg_log10_fdr", "-log10(FDR)", log_fdr),
        ("neg_log10_pvalue", "-log10(p-value)", log_pval),
        ("rank", "Rank", df[f"rank_{suffix}"]),
        ("num_obs", "numObs", df["numObs"]),
    ]


def build_entries(
    txt_files,
    string_enriched_files,
    string_depleted_files,
    string_enriched_svg_files,
    string_depleted_svg_files,
    string_enriched_id_mapping_files,
    string_depleted_id_mapping_files,
    comparisons,
    taxon_id=None,
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
        depleted["_normz_str"] = depleted["normZ"].round(4).astype(str)
        depleted["_logfdr_str"] = (
            (-np.log10(depleted["fdr_synth"].clip(lower=1e-300))).round(3).astype(str)
        )
        depleted["_numobs_str"] = depleted["numObs"].astype(str)
        depleted["_rank_str"] = depleted["rank_synth"].astype(str)

        enriched = data[data["normZ"] > 0].copy()
        enriched["x"] = np.random.RandomState(42).permutation(len(enriched))
        enriched["_normz_str"] = enriched["normZ"].round(4).astype(str)
        enriched["_logfdr_str"] = (
            (-np.log10(enriched["fdr_supp"].clip(lower=1e-300))).round(3).astype(str)
        )
        enriched["_numobs_str"] = enriched["numObs"].astype(str)
        enriched["_rank_str"] = enriched["rank_supp"].astype(str)

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
            "numObs: %{customdata[3]}<br>"
            "Rank: %{customdata[4]}"
        )
        customdata_cols = [
            "GENE",
            "_normz_str",
            "_logfdr_str",
            "_numobs_str",
            "_rank_str",
        ]

        depleted_plot = common.build_plot(
            div_id=f"{key}-depleted",
            title=f"Depleted (synthetic-lethal) genes — {comparison}",
            term_panel_title="STRING-db depleted terms",
            term_list_id=f"{key}-depleted-terms",
            df=depleted,
            y_col="normZ",
            y_title="normZ",
            color_options=_color_options(depleted, "synth"),
            id_col="GENE",
            customdata_cols=customdata_cols,
            hovertemplate=hovertemplate,
            terms=depleted_terms,
            direction="depleted",
            id_mapping=depleted_id_mapping,
            rank_col="rank_synth",
        )
        enriched_plot = common.build_plot(
            div_id=f"{key}-enriched",
            title=f"Enriched (suppressor) genes — {comparison}",
            term_panel_title="STRING-db enriched terms",
            term_list_id=f"{key}-enriched-terms",
            df=enriched,
            y_col="normZ",
            y_title="normZ",
            color_options=_color_options(enriched, "supp"),
            id_col="GENE",
            customdata_cols=customdata_cols,
            hovertemplate=hovertemplate,
            terms=enriched_terms,
            direction="enriched",
            id_mapping=enriched_id_mapping,
            rank_col="rank_supp",
        )
        plots = [
            depleted_plot,
            common.build_svg_panel(
                div_id=f"{key}-depleted-svg",
                title="STRING-db depleted network",
                svg_path=(
                    string_depleted_svg_files[i] if string_depleted_svg_files else None
                ),
                direction="depleted",
                stats_div_id=f"{key}-depleted",
                taxon_id=taxon_id,
                rank_by_gene=depleted_plot["rank_by_gene"],
            ),
            enriched_plot,
            common.build_svg_panel(
                div_id=f"{key}-enriched-svg",
                title="STRING-db enriched network",
                svg_path=(
                    string_enriched_svg_files[i] if string_enriched_svg_files else None
                ),
                direction="enriched",
                stats_div_id=f"{key}-enriched",
                taxon_id=taxon_id,
                rank_by_gene=enriched_plot["rank_by_gene"],
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
