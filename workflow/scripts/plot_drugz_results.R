# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(tidyverse)
library(cowplot)
library(ggrepel)

# Load Snakemake variables
txt <- snakemake@input[["txt"]]
pdf <- snakemake@output[["pdf"]]
fdr <- snakemake@params[["fdr"]]

# Load data
df <- read.delim(txt) %>%
  # if normZ is greater than 0, fdr is fdr_synthetic, else fdr_supp
  mutate(fdr = ifelse(normZ < 0, fdr_synth, fdr_supp))

# Data for labels (top 10 genes of lowest and highest normZ)
lowest <- df %>% top_n(10, -normZ) %>% arrange(normZ)
highest <- df %>% top_n(10, normZ) %>% arrange(normZ)

# Plot normZ against fdr
p <- ggplot(df, aes(x = normZ, y = -log10(fdr))) +
  geom_point(size = 8,
             alpha = 0.4) +
  theme_cowplot(18) +
  geom_hline(yintercept = -log10(fdr), linetype = "dashed") +
  geom_label_repel(data = rbind(lowest, highest), 
                  aes(label = GENE))

ggsave(pdf, p)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
