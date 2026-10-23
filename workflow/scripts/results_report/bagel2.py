"""Per-comparison BAGEL2 depleted (essential) scatter entries.

BAGEL2 only ever reports depleted/essential genes — there is no enriched
set to split out, unlike MAGeCK and DrugZ.
"""

import logging

import numpy as np
import pandas as pd

from . import common


def _color_options(df):
    """BAGEL2's .pr columns, offered as "Color by" choices for the dot
    plot's marker fill."""
    log_fdr = -np.log10(df["FDR"].clip(lower=1e-300))
    return [
        ("neg_log10_fdr", "-log10(FDR)", log_fdr),
        ("recall", "Recall", df["Recall"]),
        ("precision", "Precision", df["Precision"]),
    ]


def build_entries(
    txt_files,
    string_depleted_files,
    string_depleted_svg_files,
    string_depleted_id_mapping_files,
    comparisons,
    taxon_id=None,
):
    entries = []
    for i, comparison in enumerate(comparisons):
        txt = txt_files[i]
        logging.info(f"Loading {txt}")
        data = pd.read_csv(txt, sep="\t")
        data["x"] = np.random.RandomState(42).permutation(len(data))
        data["_bf_str"] = data["BF"].round(3).astype(str)
        data["_fdr_str"] = data["FDR"].round(6).astype(str)
        # BAGEL2's .pr has no rank column — higher BF means more confidently
        # essential/depleted, so rank 1 is the highest BF.
        data["rank"] = data["BF"].rank(ascending=False, method="min").astype(int)
        data["_rank_str"] = data["rank"].astype(str)

        logging.info(f"Loaded {len(data)} genes")

        terms = (
            common.load_string_terms(string_depleted_files[i])
            if string_depleted_files
            else []
        )
        id_mapping = (
            common.load_id_mapping(string_depleted_id_mapping_files[i])
            if string_depleted_id_mapping_files
            else {}
        )

        key = common.slugify(f"bagel2__{comparison}")
        plot = common.build_plot(
            div_id=f"{key}-depleted",
            title=f"Depleted (essential) genes — {comparison}",
            term_panel_title="STRING-db depleted terms",
            term_list_id=f"{key}-depleted-terms",
            df=data,
            y_col="BF",
            y_title="Bayes Factor (BF)",
            color_options=_color_options(data),
            id_col="Gene",
            customdata_cols=["Gene", "_bf_str", "_fdr_str", "_rank_str"],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>BF: %{customdata[1]}<br>"
                "FDR: %{customdata[2]}<br>Rank: %{customdata[3]}"
            ),
            terms=terms,
            direction="depleted",
            id_mapping=id_mapping,
            rank_col="rank",
        )
        svg_panel = common.build_svg_panel(
            div_id=f"{key}-depleted-svg",
            title="STRING-db depleted network",
            svg_path=(
                string_depleted_svg_files[i] if string_depleted_svg_files else None
            ),
            direction="depleted",
            stats_div_id=f"{key}-depleted",
            taxon_id=taxon_id,
            rank_by_gene=plot["rank_by_gene"],
        )
        plots = [p for p in [plot, svg_panel] if p is not None]

        entries.append(
            dict(
                key=key,
                label=comparison,
                plots=plots,
            )
        )
    return entries
