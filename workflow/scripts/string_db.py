import io
import logging
import re

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Set up logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

# Load Snakemake variables
pathway_data = snakemake.wildcards["pathway_data"]
results_txt = snakemake.input["txt"]
svg_output = snakemake.output["svg"]
csv_output = snakemake.output["csv"]
id_mapping_output = snakemake.output["id_mapping"]
fdr_cutoff = snakemake.config["stats"]["string_db"]["fdr"]
top_genes = snakemake.config["stats"]["string_db"]["top_genes"]
organism = snakemake.config["lib_info"]["species"]
data_source = snakemake.params["data"]

# Column names differ per upstream tool
if data_source == "mageck":
    id_column = "id"
    rank_column = "pos|rank" if pathway_data == "enriched" else "neg|rank"
    fdr_column = "pos|fdr" if pathway_data == "enriched" else "neg|fdr"
elif data_source == "drugz":
    id_column = "GENE"
    rank_column = "rank_supp" if pathway_data == "enriched" else "rank_synth"
    fdr_column = "fdr_supp" if pathway_data == "enriched" else "fdr_synth"
elif data_source == "bagel2":
    # BAGEL2's .bf file only ever represents depleted (essential) genes,
    # as a raw BF value rather than a rank/FDR column: a higher BF means
    # higher confidence the gene is essential/depleted.
    if pathway_data != "depleted":
        message = (
            f"BAGEL2 only supports 'depleted' STRING-db analysis, got {pathway_data!r}."
        )
        logging.error(message)
        raise ValueError(message)
    id_column = "GENE"
    bf_cutoff = snakemake.config["stats"]["string_db"]["bf_cutoff"]
elif data_source == "mageck_mle":
    # MLE reports one beta/fdr block per design-matrix condition (no
    # separate pos|/neg| columns like MAGeCK test), and beta can be
    # either sign within that one column — enriched/depleted here means
    # beta > 0 / beta < 0, not a pre-split column.
    condition = snakemake.wildcards["condition"]
    id_column = "Gene"
    beta_column = f"{condition}|beta"
    fdr_column = f"{condition}|fdr"
else:
    message = (
        f"Unsupported data source: {data_source}. "
        "Please use 'mageck', 'mageck_mle', 'drugz', or 'bagel2'."
    )
    logging.error(message)
    raise ValueError(message)

# Set NCBI taxon ID based on species
if organism == "human":
    taxon_id = 9606
elif organism == "mouse":
    taxon_id = 10090
else:
    message = f"Unsupported species: {organism}. Please use 'human' or 'mouse'."
    logging.error(message)
    raise ValueError(message)

# Load upstream results
data = pd.read_csv(results_txt, sep="\t")

# Subset to top genes based on FDR/BF cutoff or top N genes
if data_source == "bagel2":
    if top_genes > 0:
        logging.info(f"Selecting top {top_genes} genes based on BF ranking.")
        data = data.sort_values(by="BF", ascending=False).head(top_genes)
    else:
        logging.info(f"Filtering genes based on BF cutoff of {bf_cutoff}.")
        data = data[data["BF"] > bf_cutoff]
elif data_source == "mageck_mle":
    data = (
        data[data[beta_column] > 0]
        if pathway_data == "enriched"
        else data[data[beta_column] < 0]
    )
    if top_genes > 0:
        logging.info(f"Selecting top {top_genes} genes based on {pathway_data} FDR.")
        data = data.sort_values(by=fdr_column, ascending=True).head(top_genes)
    else:
        logging.info(f"Filtering genes based on FDR threshold of {fdr_cutoff}.")
        data = data[data[fdr_column] < fdr_cutoff]
elif top_genes > 0:
    logging.info(f"Selecting top {top_genes} genes based on {pathway_data} ranking.")
    data = data.sort_values(by=rank_column, ascending=True).head(top_genes)
else:
    logging.info(f"Filtering genes based on FDR threshold of {fdr_cutoff}.")
    data = data[data[fdr_column] < fdr_cutoff]

# Get genes for STRING-db analysis
genes = data[id_column].tolist()

# Raise error if no genes are found after filtering
if not genes:
    message = (
        f"No genes found for {pathway_data} analysis after filtering. "
        f"Please check your {data_source} results and filtering criteria."
    )
    logging.error(message)
    raise ValueError(message)

# Format input genes as a carriage-return (\r) separated payload
genes_payload = "%0d".join(genes)

# Session with retry: backoff on 400 (STRING returns this, not 429, when
# concurrent CI jobs exceed its request rate), 429/500/502/503/504, up to 5
# attempts. backoff_factor=2 gives waits of 2, 4, 8, 16 s between retries.
# allowed_methods must include POST explicitly — urllib3's Retry excludes
# non-idempotent methods by default, so without this every retry setting
# here is a silent no-op for our session.post() calls.
retry = Retry(
    total=5,
    backoff_factor=2,
    status_forcelist=[400, 429, 500, 502, 503, 504],
    allowed_methods=frozenset(["POST"]),
    raise_on_status=False,
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retry))

params = {
    "identifiers": genes_payload,
    "species": taxon_id,
}

# STRING resolves each submitted identifier to its current canonical gene
# symbol (e.g. the deprecated alias FAM192A resolves to PSME3IP1), and that
# resolved name — not our submitted identifier — is what labels the node in
# the network SVG. Fetch that mapping so results_report can look up a
# node's stats by the name STRING actually displays.
id_mapping_url = "https://string-db.org/api/tsv/get_string_ids"
id_mapping_response = session.post(id_mapping_url, data=params, timeout=60)

if id_mapping_response.status_code == 200:
    mapping_df = pd.read_csv(io.StringIO(id_mapping_response.text), sep="\t")
    mapping_df["query_gene"] = mapping_df["queryIndex"].apply(lambda i: genes[i])
    mapping_df = mapping_df[["query_gene", "preferredName"]].rename(
        columns={"preferredName": "preferred_name"}
    )
    mapping_df.to_csv(id_mapping_output, index=False)
    logging.info(
        f"Saved {len(mapping_df)} resolved gene identifiers to {id_mapping_output}."
    )
else:
    raise RuntimeError(
        f"STRING-db identifier resolution request failed with status "
        f"{id_mapping_response.status_code}: {id_mapping_response.text}"
    )

# Query the STRING Enrichment endpoint
enrichment_url = "https://string-db.org/api/tsv/enrichment"

response = session.post(enrichment_url, data=params, timeout=60)

if response.status_code == 200:
    df = pd.read_csv(io.StringIO(response.text), sep="\t")
    df_significant = df[df["fdr"] < 0.05]
    df_significant.to_csv(csv_output, index=False)
    logging.info(
        f"Saved {len(df_significant)} significant functional pathways to {csv_output}."
    )
else:
    raise RuntimeError(
        f"STRING-db enrichment request failed with status {response.status_code}: {response.text}"
    )

# Download network SVG. Nodes get a structure-thumbnail image whenever
# STRING happens to have 3D structure data for that gene, which looks like
# a bordered node versus a flat circle for genes without one — purely a
# per-gene coincidence, not something meaningful to this analysis. Disable
# it so node appearance is consistent across all networks.
network_url = "https://string-db.org/api/svg/network"
svg_params = {**params, "block_structure_pics_in_bubbles": 1}
svg_response = session.post(network_url, data=svg_params, timeout=60)

if svg_response.status_code == 200:
    # STRING always draws a radial-gradient "3D bubble" highlight circle on
    # top of each flat node circle (no API parameter disables it, unlike
    # the structure thumbnails above), which reads as a border/ring around
    # otherwise-flat nodes. Strip those overlay circles for a truly flat
    # look; the unused gradient <defs> they reference are harmless to leave.
    svg_text = re.sub(
        r"<circle[^>]*fill='url\(#bubble_gradient[^']*'[^>]*/>",
        "",
        svg_response.text,
    )
    with open(svg_output, "w") as f:
        f.write(svg_text)
    logging.info(f"Saved network SVG to {svg_output}.")
else:
    raise RuntimeError(
        f"STRING-db network SVG request failed with status {svg_response.status_code}: {svg_response.text}"
    )
