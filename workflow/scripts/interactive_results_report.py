"""Build one combined interactive HTML report for every MAGeCK, DrugZ and
BAGEL2 comparison, replacing the old per-comparison interactive report
files. Navigate via tool tabs and a per-tool comparison dropdown.
"""

import logging
import os

from results_report import bagel2, drugz, mageck, report_html

log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

output_html = snakemake.output["html"]

mageck_pairs = snakemake.params["mageck_pairs"]
drugz_comparisons = snakemake.params["drugz_comparisons"]
bagel2_comparisons = snakemake.params["bagel2_comparisons"]

mageck_entries = (
    mageck.build_entries(
        snakemake.input["mageck_txt"],
        snakemake.input["mageck_string_enriched"],
        snakemake.input["mageck_string_depleted"],
        mageck_pairs,
    )
    if mageck_pairs
    else []
)
logging.info(f"Built {len(mageck_entries)} MAGeCK comparison entries")

drugz_entries = (
    drugz.build_entries(
        snakemake.input["drugz_txt"],
        snakemake.input["drugz_string_enriched"],
        snakemake.input["drugz_string_depleted"],
        drugz_comparisons,
    )
    if drugz_comparisons
    else []
)
logging.info(f"Built {len(drugz_entries)} DrugZ comparison entries")

bagel2_entries = (
    bagel2.build_entries(
        snakemake.input["bagel2_txt"],
        snakemake.input["bagel2_string_depleted"],
        bagel2_comparisons,
    )
    if bagel2_comparisons
    else []
)
logging.info(f"Built {len(bagel2_entries)} BAGEL2 comparison entries")

html = report_html.render(
    [
        ("mageck", "MAGeCK", mageck_entries),
        ("drugz", "DrugZ", drugz_entries),
        ("bagel2", "BAGEL2", bagel2_entries),
    ]
)

os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as fh:
    fh.write(html)
logging.info(f"Written combined results report to {output_html}")
