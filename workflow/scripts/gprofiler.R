# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(tidyverse)
library(gprofiler2)
library(cowplot)

# Load Snakemake variables
txt <- snakemake@input[["txt"]]
csv <- snakemake@output[["csv"]]
pdf <- snakemake@output[["pdf"]]
top_genes <- snakemake@params[["top_genes"]]
fdr <- snakemake@params[["fdr"]]
comparison <- snakemake@wildcards[["comparison"]]
dt <- snakemake@wildcards[["pathway_data"]]
data <- snakemake@params[["data"]]

### Load data and extract genes
df <- read.delim(txt)

if (data == "mageck") {
  if (dt == "enriched") {
    rank <- "pos.rank"
    dt_fdr <- "pos.fdr"
  } else {
    rank <- "neg.rank"
    dt_fdr <- "neg.fdr"
  } 
  id_column <- "id"
} else if (data == "bagel2") {
  # BAGEL2 will only have data for depleted genes
  rank <- "BF"
  dt_fdr <- "FDR"
  id_column <- "Gene"
} else if (data == "drugz") {
  if (dt == "enriched") {
    rank <- "rank_supp"
    dt_fdr <- "fdr_supp"
  } else {
    rank <- "rank_synth"
    dt_fdr <- "fdr_synth"
  } 
  id_column <- "GENE"
}

if (top_genes > 0) {
  genes <- df %>%
    arrange(get(rank)) %>%
    head(top_genes) %>%
    pull(get(id_column))
} else {
  genes <- df %>%
    filter(get(dt_fdr) < fdr) %>%
    pull(get(id_column))
}

# Print genes to stdout on a separate line
cat("Genes for", dt, "analysis:\n", paste(genes, collapse = "\n"), "\n")

### Run gprofiler2
res <- gost(genes,
            sources = c("GO:BP", "GO:MF", "KEGG", "REAC"),
            correction_method = "fdr",
            significant = TRUE,
            evcodes = TRUE)

# Flatten and clean up results before saving to CSV
if (!is.null(res)) {
  result_df <- as.data.frame(res$result)
  
  # Ensure all list columns are converted to strings
  clean_result_df <- result_df %>%
    mutate(across(where(is.list), ~ map_chr(.x, ~ paste(.x, collapse = ";"))))
  
  # Remove evidence_codes column
  clean_result_df <- clean_result_df %>%
    dplyr::select(-evidence_codes)
  
  # Save results
  write.csv(clean_result_df, csv, row.names = FALSE)
} else {
  warning("No significant results returned by gprofiler2.")
  write.csv(data.frame(), csv, row.names = FALSE)  # Write an empty CSV if no results
}

### Plot results
if (!is.null(res)) {
  # Get the top 3 terms of each category
  top_terms <- res$result %>%
    group_by(source) %>%
    top_n(3, -log10(p_value)) %>%
    ungroup() %>%
    arrange(-log10(p_value)) %>%
    pull(term_id)
  
  # Create plot and save
  p <- gostplot(res, 
                capped = FALSE, 
                interactive = FALSE) +
    theme_cowplot(18) +
    theme(legend.position = "none") +
    ggtitle(paste(dt, "in", comparison))
  
  publish_gostplot(p,
                   highlight_terms = top_terms,
                   width = 12,
                   height = 16,
                   filename = pdf)
} else {
  warning("No significant results returned by gprofiler2.")
  # Write empty PDF
  pdf(file = pdf, width = 12, height = 16)
  dev.off()
}
file.remove("Rplots.pdf", showWarnings = FALSE)  # Remove the default Rplots.pdf

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
