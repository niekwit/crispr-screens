# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# ADAPTED (ANNOTATED AND SIMPLIFIED) FROM https://github.com/WubingZhang/MAGeCKFlute/blob/master/R/sgRankView.R

# load required libraries
library(tidyverse)

# load data
data <- read.delim(snakemake@input[["sg"]])

# load gene rank data
rank <- read.delim(snakemake@input[["rank"]]) %>%
  dplyr::select(id, neg.rank, pos.rank) %>%
  rename(Gene = id)

# Add Gene rank to data
data <- data %>% left_join(rank, by = "Gene") 

# set parameters
select <- 5 # number of genes to plot for enriched and depleted
binwidth <- 0.3 # 
interval <- 0.1

# select top x genes for enrichment
df.enriched <- data %>% filter(pos.rank <= select) 

# select top x genes with lowest median lfc
df.depleted <- data %>% filter(neg.rank <= select) 

# add top and bottom genes data together and get genes to plot
df <- rbind(df.enriched, df.depleted)
genes <- unique(df$Gene)

# add index, y values for plotting and colour param
df <- df %>%
  mutate(Gene = factor(Gene, levels = genes),
         index = rep(1:length(genes), as.numeric(table(Gene)[genes])),
         y = (binwidth+interval)*(index-1),
         yend = (binwidth+interval)*index-interval,
         colour = ifelse(LFC > 0, "pos", "neg")) %>%
  dplyr::select(c("sgrna", "Gene", "LFC", "y", "yend", "colour", "index")) %>%
  as.data.frame()
  
# set the scale of x-axis
a <- -Inf
b <- Inf

# set values for rectangle dimensions (x/y values)/colour fill to plot sgRNAs lfc inside
bgcol <- as.vector(sapply(seq(1, max(df$index), 1), function(x){rep(x, 4)})) %>%
  as.data.frame() %>%
  rename("id" = 1) %>%
  mutate(x = rep(c(a,b,b,a),max(df$index)),
         y = as.vector(sapply(seq(1,max(df$index),1), function(x){
           c((interval + binwidth)*(x-1), (interval + binwidth)*(x-1),
             (interval + binwidth)*x-interval,(interval + binwidth)*x-interval)
         })))

# create plot
p <- ggplot() +
  geom_polygon(
    # Use aes() with unquoted column names
    aes(x = x, y = y, group = id),
    fill = "#dedede",
    color = "gray20",
    data = bgcol
  ) +
  geom_segment(
    # Use aes() with unquoted column names
    aes(x = LFC, y = y, xend = LFC, yend = yend, color = colour),
    data = df
  ) +
  scale_color_manual(values = c("pos" = "#e41a1c", "neg" = "#377eb8")) +
  scale_y_continuous(
    breaks = bgcol$y[seq(1, nrow(bgcol), 4)] + binwidth / 2,
    labels = genes, expand = c(0, 0)
  ) +
  labs(
    x = "Log2(Fold change)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), panel.background = element_blank()
  )

# save plot to file
ggsave(plot = p, 
       filename = snakemake@output[[1]], 
       units = "in", 
       width = 10, 
       height = 7)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
