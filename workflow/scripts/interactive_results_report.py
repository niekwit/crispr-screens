"""Build one combined interactive HTML report for every MAGeCK, DrugZ and
BAGEL2 comparison (MAGeCK also gets its sgRNA-rank view), replacing the old
per-comparison interactive report files. Navigate via tool tabs and a
per-tool comparison dropdown.
"""

import logging
import os

from results_report import (
    bagel2,
    common,
    drugz,
    mageck,
    mageck_mle,
    report_html,
    sgrank,
)

log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

output_html = snakemake.output["html"]

mageck_pairs = snakemake.params["mageck_pairs"]
mageck_mle_pairs = snakemake.params["mageck_mle_pairs"]
drugz_comparisons = snakemake.params["drugz_comparisons"]
bagel2_comparisons = snakemake.params["bagel2_comparisons"]
taxon_id = common.taxon_id_for(snakemake.params["species"])

mageck_entries = (
    mageck.build_entries(
        snakemake.input["mageck_txt"],
        snakemake.input["mageck_string_enriched"],
        snakemake.input["mageck_string_depleted"],
        snakemake.input["mageck_string_enriched_svg"],
        snakemake.input["mageck_string_depleted_svg"],
        snakemake.input["mageck_string_enriched_id_mapping"],
        snakemake.input["mageck_string_depleted_id_mapping"],
        mageck_pairs,
        taxon_id,
    )
    if mageck_pairs
    else []
)
logging.info(f"Built {len(mageck_entries)} MAGeCK comparison entries")

mageck_sgrna_txt = snakemake.input["mageck_sgrna_txt"]
for i, entry in enumerate(mageck_entries):
    comparison, cnv = mageck_pairs[i]
    entry["plots"] += sgrank.build_plots(
        snakemake.input["mageck_txt"][i],
        mageck_sgrna_txt[i],
        entry["key"],
        comparison,
        cnv,
    )

mageck_mle_entries = (
    mageck_mle.build_entries(
        snakemake.input["mageck_mle_txt"],
        mageck_mle_pairs,
        taxon_id,
    )
    if mageck_mle_pairs
    else []
)
logging.info(f"Built {len(mageck_mle_entries)} MAGeCK MLE entries")

drugz_entries = (
    drugz.build_entries(
        snakemake.input["drugz_txt"],
        snakemake.input["drugz_string_enriched"],
        snakemake.input["drugz_string_depleted"],
        snakemake.input["drugz_string_enriched_svg"],
        snakemake.input["drugz_string_depleted_svg"],
        snakemake.input["drugz_string_enriched_id_mapping"],
        snakemake.input["drugz_string_depleted_id_mapping"],
        drugz_comparisons,
        taxon_id,
    )
    if drugz_comparisons
    else []
)
logging.info(f"Built {len(drugz_entries)} DrugZ comparison entries")

bagel2_entries = (
    bagel2.build_entries(
        snakemake.input["bagel2_txt"],
        snakemake.input["bagel2_string_depleted"],
        snakemake.input["bagel2_string_depleted_svg"],
        snakemake.input["bagel2_string_depleted_id_mapping"],
        bagel2_comparisons,
        taxon_id,
    )
    if bagel2_comparisons
    else []
)
logging.info(f"Built {len(bagel2_entries)} BAGEL2 comparison entries")

html = report_html.render(
    [
        ("mageck", "MAGeCK", mageck_entries),
        ("mageck_mle", "MAGeCK MLE", mageck_mle_entries),
        ("drugz", "DrugZ", drugz_entries),
        ("bagel2", "BAGEL2", bagel2_entries),
    ]
)

os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as fh:
    fh.write(html)
logging.info(f"Written combined results report to {output_html}")
