"""Per-comparison BAGEL2 depleted (essential) scatter entries.

BAGEL2 only ever reports depleted/essential genes — there is no enriched
set to split out, unlike MAGeCK and DrugZ.
"""

import logging

import numpy as np
import pandas as pd

from . import common


def build_entries(txt_files, string_depleted_files, comparisons):
    entries = []
    for i, comparison in enumerate(comparisons):
        txt = txt_files[i]
        logging.info(f"Loading {txt}")
        data = pd.read_csv(txt, sep="\t")
        data["x"] = np.random.RandomState(42).permutation(len(data))
        data["log_fdr"] = -np.log10(data["FDR"].clip(lower=1e-300))
        data["_bf_str"] = data["BF"].round(3).astype(str)
        data["_fdr_str"] = data["FDR"].round(6).astype(str)

        logging.info(f"Loaded {len(data)} genes")

        terms = (
            common.load_string_terms(string_depleted_files[i])
            if string_depleted_files
            else []
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
            color_col="log_fdr",
            color_title="-log10(FDR)",
            id_col="Gene",
            customdata_cols=["Gene", "_bf_str", "_fdr_str"],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>BF: %{customdata[1]}<br>FDR: %{customdata[2]}"
            ),
            terms=terms,
        )

        entries.append(
            dict(
                key=key,
                label=comparison,
                plots=[plot],
            )
        )
    return entries
