# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(tidyverse)
library(enrichR)
library(cowplot)

# Load Snakemake variables
txt <- snakemake@input[["txt"]]
csvs <- snakemake@output[["csv"]]
plots <- snakemake@output[["plots"]]
terms <- snakemake@params[["terms"]]
dbs <- snakemake@params[["dbs"]]
fdr <- snakemake@params[["fdr"]]
top_genes <- snakemake@params[["top_genes"]]
data_type <- snakemake@params[["data_type"]]

# Create output directories
plots_out_dir <- dirname(plots[1])
dir.create(plots_out_dir, showWarnings = FALSE)
csv_out_dir <- dirname(csvs[1])
dir.create(csv_out_dir, showWarnings = FALSE)

# Load available databases
available_dbs <- listEnrichrDbs()

# Check if the requested databases are available
dbs_check <- dbs %in% available_dbs$libraryName
if (any(!dbs_check)) {
  print("Following databases are not available:")
  print(dbs[!dbs_check])
  print(paste("Available databases are:", paste(available_dbs$libraryName, collapse = ", ")))
  stop()
} else {
  print("All requested databases are available...")
}

# Load MAGeCK results
df <- read.delim(txt) 

pathway_analysis <- function(dt) {
  if (dt == "enriched") {
      rank <- "pos.rank"
      dt_fdr <- "pos.fdr"
  } else {
      rank <- "neg.rank"
      dt_fdr <- "neg.fdr"
  } 
  
  if (top_genes > 0) {
      genes <- df %>%
          arrange(desc(rank)) %>%
          head(top_genes) %>%
          pull(id)
  } else {
      genes <- df %>%
          filter(dt_fdr < fdr) %>%
          pull(id)
    if (length(genes) == 0) {
      stop("No genes found with FDR < ", fdr)
    }
  }
  # Run EnrichR
  enrichr_results <- enrichr(genes, dbs)
  
  # Plot results for each db
  for (db in dbs) {
    tmp <- enrichr_results[[db]] %>%
      mutate(neg.log.P.value = -log10(Adjusted.P.value)) %>%
      arrange(desc(neg.log.P.value))
    
    # Save results to csv file
    write.csv(enrichr_results[[db]], 
              file = paste0(csv_out_dir, "/", db, "_", dt,".csv"),
              quote = FALSE,
              row.names = FALSE)
    
    # Calculate ratio of genes found and total genes in the term
    tmp$Ratio <- str_split(tmp$Overlap, "/") %>%
      map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2]))
    
    # Wrap terms longer than 60 characters over multiple lines
    tmp$Term <- str_wrap(tmp$Term, 60)
    
    # Relevel terms to avoid alphabetical sorting
    tmp$Term <- fct_rev(factor(tmp$Term, levels = tmp$Term))
    
    # Check whether the number of terms is less than the amout of results
    if (nrow(tmp) < terms) {
      terms <- nrow(tmp)
    }
    
    # Plot
    p <- ggplot(tmp[1:terms,], aes(x = neg.log.P.value,
                                  y = Term)) +
      geom_point(aes(fill = Ratio), 
                 alpha = 1,
                 size = 12,
                 shape = 21,
                 colour = "black") +
      scale_fill_gradient(low = "white", 
                          high = "forestgreen",
                          guide = guide_colorbar(frame.colour = "black", 
                                                 ticks.colour = "black")) +
      theme_cowplot(18) +
      theme(axis.text.x = element_text(angle = 45, 
                                       hjust = 1)) +
      labs(x = "-log10(Adjusted P value)", 
           y = NULL)
    
    # Save plot
    ggsave(paste0(plots_out_dir, "/", db, "_", dt, ".pdf"),
           p, 
           width = 12, 
           height = terms * 0.9)
  }
}

lapply(data_type, pathway_analysis)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")