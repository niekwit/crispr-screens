"""Per-(matrix, cnv) MAGeCK MLE scatter entries.

Unlike MAGeCK test (one enriched/depleted pair per comparison), MLE reports
one beta/z/p-value/fdr/wald-p-value/wald-fdr block per non-baseline column
of the design matrix — the "direction" tabs here are per-condition rather
than enriched/depleted. Each condition gets a single scatter of every gene
(not split into separate enriched/depleted scatters, since beta can be
either sign within one condition), paired with two STRING-db network
panels — enriched (beta > 0) and depleted (beta < 0) — fed by gene subsets
chosen by sign + FDR rather than MAGeCK's own pos|rank/neg|rank columns
(MLE has no equivalent column), computed in string_db.py's "mageck_mle"
branch.
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


def _rank_by_gene(df, id_mapping):
    """gene -> rank (1 = best), by ascending FDR within the sign-filtered
    subset already selected by the caller — independent of
    common.build_plot's own rank_by_gene (which assumes one rank column
    for the whole df, not a sign-split one), also alias-resolved to
    STRING's preferred names like common.build_plot does."""
    rank_by_gene = {gene: i + 1 for i, gene in enumerate(df["Gene"])}
    for query_gene, preferred_name in id_mapping.items():
        if query_gene in rank_by_gene and preferred_name not in rank_by_gene:
            rank_by_gene[preferred_name] = rank_by_gene[query_gene]
    return rank_by_gene


def build_entries(
    txt_files,
    pairs,
    triples,
    string_enriched_files,
    string_depleted_files,
    string_enriched_svg_files,
    string_depleted_svg_files,
    string_enriched_id_mapping_files,
    string_depleted_id_mapping_files,
    taxon_id=None,
):
    """pairs: list of (matrix, cnv) tuples, positionally matching txt_files
    (one gene_summary.txt per design matrix x cnv). triples: list of
    (matrix, cnv, condition), positionally matching the string_* file
    lists (one condition's STRING-db files per triple; the string_*
    lists are empty altogether when STRING-db/that direction is disabled)."""
    data_by_pair = {}
    for i, (matrix, cnv) in enumerate(pairs):
        logging.info(f"Loading {txt_files[i]}")
        data_by_pair[(matrix, cnv)] = pd.read_csv(txt_files[i], sep="\t")

    plots_by_pair = {}
    for j, (matrix, cnv, cond) in enumerate(triples):
        data = data_by_pair[(matrix, cnv)]
        np.random.seed(42)
        key = common.slugify(f"mageck_mle__{cnv}__{matrix}")
        div_key = common.slugify(f"{key}-{cond}")

        df = data.copy()
        df["x"] = np.random.permutation(len(df))
        df["_beta_str"] = df[f"{cond}|beta"].round(4).astype(str)
        df["_z_str"] = df[f"{cond}|z"].round(3).astype(str)
        df["_fdr_str"] = df[f"{cond}|fdr"].round(6).astype(str)
        df["_wald_fdr_str"] = df[f"{cond}|wald-fdr"].round(6).astype(str)

        enriched_terms = (
            common.load_string_terms(string_enriched_files[j])
            if string_enriched_files
            else []
        )
        depleted_terms = (
            common.load_string_terms(string_depleted_files[j])
            if string_depleted_files
            else []
        )
        terms = [{**t, "label": f"[Enriched] {t['label']}"} for t in enriched_terms] + [
            {**t, "label": f"[Depleted] {t['label']}"} for t in depleted_terms
        ]

        enriched_id_mapping = (
            common.load_id_mapping(string_enriched_id_mapping_files[j])
            if string_enriched_id_mapping_files
            else {}
        )
        depleted_id_mapping = (
            common.load_id_mapping(string_depleted_id_mapping_files[j])
            if string_depleted_id_mapping_files
            else {}
        )
        id_mapping = {**enriched_id_mapping, **depleted_id_mapping}

        hovertemplate = (
            "<b>%{customdata[0]}</b><br>"
            "Beta: %{customdata[1]}<br>"
            "Z: %{customdata[2]}<br>"
            "FDR: %{customdata[3]}<br>"
            "Wald FDR: %{customdata[4]}"
        )
        customdata_cols = ["Gene", "_beta_str", "_z_str", "_fdr_str", "_wald_fdr_str"]

        plot = common.build_plot(
            div_id=div_key,
            title=f"MAGeCK MLE — {matrix} ({cnv}) — {cond}",
            term_panel_title="STRING-db terms",
            term_list_id=common.slugify(f"{div_key}-terms"),
            df=df,
            y_col=f"{cond}|beta",
            y_title="Beta score",
            color_options=_color_options(df, cond),
            id_col="Gene",
            customdata_cols=customdata_cols,
            hovertemplate=hovertemplate,
            terms=terms,
            direction=cond,
            id_mapping=id_mapping or None,
            taxon_id=taxon_id,
        )

        entry_plots = [plot]
        if string_enriched_svg_files:
            enriched_df = df[df[f"{cond}|beta"] > 0].sort_values(f"{cond}|fdr")
            entry_plots.append(
                common.build_svg_panel(
                    div_id=common.slugify(f"{div_key}-enriched-svg"),
                    title=f"STRING-db enriched network — {cond}",
                    svg_path=string_enriched_svg_files[j],
                    direction=cond,
                    stats_div_id=div_key,
                    taxon_id=taxon_id,
                    rank_by_gene=_rank_by_gene(enriched_df, enriched_id_mapping),
                )
            )
        if string_depleted_svg_files:
            depleted_df = df[df[f"{cond}|beta"] < 0].sort_values(f"{cond}|fdr")
            entry_plots.append(
                common.build_svg_panel(
                    div_id=common.slugify(f"{div_key}-depleted-svg"),
                    title=f"STRING-db depleted network — {cond}",
                    svg_path=string_depleted_svg_files[j],
                    direction=cond,
                    stats_div_id=div_key,
                    taxon_id=taxon_id,
                    rank_by_gene=_rank_by_gene(depleted_df, depleted_id_mapping),
                )
            )
        entry_plots = [p for p in entry_plots if p is not None]
        plots_by_pair.setdefault((matrix, cnv), []).extend(entry_plots)

    entries = []
    for matrix, cnv in pairs:
        key = common.slugify(f"mageck_mle__{cnv}__{matrix}")
        entries.append(
            dict(
                key=key,
                label=f"{matrix} ({cnv})",
                plots=plots_by_pair.get((matrix, cnv), []),
            )
        )
    return entries
