import io
import logging
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
else:
    message = f"Unsupported data source: {data_source}. Please use 'mageck' or 'drugz'."
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

# Subset to top genes based on FDR or top N genes
if top_genes > 0:
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

# Session with retry: backoff on 429/500/502/503/504, up to 5 attempts.
# backoff_factor=2 gives waits of 2, 4, 8, 16 s between retries.
retry = Retry(
    total=5,
    backoff_factor=2,
    status_forcelist=[429, 500, 502, 503, 504],
    raise_on_status=False,
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retry))

# Query the STRING Enrichment endpoint
enrichment_url = "https://string-db.org/api/tsv/enrichment"
params = {
    "identifiers": genes_payload,
    "species": taxon_id,
}

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

# Download network SVG
network_url = "https://string-db.org/api/svg/network"
svg_response = session.post(network_url, data=params, timeout=60)

if svg_response.status_code == 200:
    with open(svg_output, "w") as f:
        f.write(svg_response.text)
    logging.info(f"Saved network SVG to {svg_output}.")
else:
    raise RuntimeError(
        f"STRING-db network SVG request failed with status {svg_response.status_code}: {svg_response.text}"
    )
