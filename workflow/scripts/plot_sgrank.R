# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# ADAPTED (ANNOTATED AND SIMPLIFIED) FROM https://github.com/WubingZhang/MAGeCKFlute/blob/master/R/sgRankView.R

# load required libraries
library(tidyverse)
library(cowplot)

# load data
data <- read.delim(snakemake@input[["sg"]])

# load gene rank data
rank <- read.delim(snakemake@input[["gene"]]) %>%
  dplyr::select(id, neg.rank, pos.rank) %>%
  rename(Gene = id)

# Add Gene rank to data
data <- data %>% left_join(rank, by = "Gene")

# set parameters
select <- 5 # number of genes to plot for enriched and depleted
binwidth <- 0.3
interval <- 0.1

# shared x-axis range across all panels, based on the full sgRNA LFC distribution
x.range <- range(data$LFC, na.rm = TRUE)
x.pad <- diff(x.range) * 0.02
x.limits <- c(x.range[1] - x.pad, x.range[2] + x.pad)

# density of the full sgRNA LFC distribution (background reference)
dens <- density(
  data$LFC,
  na.rm = TRUE,
  n = 512,
  from = x.limits[1],
  to = x.limits[2]
)
dens.df <- tibble(x = dens$x, dens = dens$y)
max.dens <- max(dens.df$dens)
# slightly overlapping tile width avoids hairline rendering seams between
# adjacent density tiles in the gene-row shading
dens.tile.width <- diff(range(dens.df$x)) / (nrow(dens.df) - 1) * 1.02

# density panel: just the curve, no separate gradient strip (the density
# shading now lives inside the gene rows themselves, see build_sgrank_plot)
p.density <- ggplot(dens.df, aes(x = x, y = dens)) +
  geom_area(fill = "grey88", colour = "grey20", linewidth = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, max.dens * 1.05),
    breaks = scales::breaks_pretty(n = 3)(c(0, max.dens)),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = x.limits) +
  labs(y = "Density") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# builds the gene-rank row panel for a set of top genes and stacks it under
# the shared density panel (colour: red for enriched genes, blue for depleted)
build_sgrank_plot <- function(df, genes, colour) {
  n.genes <- length(genes)
  row.centre <- function(i) (interval + binwidth) * (i - 1) + binwidth / 2

  df <- df %>%
    mutate(
      Gene = factor(Gene, levels = genes),
      # as.integer(Gene) (not row position) so each sgRNA lines up with its
      # gene's row regardless of the row order genes/df happen to be in
      index = as.integer(Gene),
      y = (binwidth + interval) * (index - 1),
      yend = (binwidth + interval) * index - interval
    ) %>%
    dplyr::select(c("sgrna", "Gene", "LFC", "y", "yend", "index")) %>%
    as.data.frame()

  # density-shaded background per gene row (white = low density, black =
  # high density), echoing the overall LFC distribution behind each gene
  shade.df <- do.call(
    rbind,
    lapply(seq_len(n.genes), function(i) {
      dens.df %>% mutate(y = row.centre(i))
    })
  )

  # row outline, drawn on top of the shading with no fill of its own
  box.df <- tibble(
    index = seq_len(n.genes),
    ymin = (interval + binwidth) * (index - 1),
    ymax = (interval + binwidth) * index - interval
  )

  p.ranks <- ggplot() +
    geom_tile(
      aes(x = x, y = y, height = binwidth, fill = dens),
      data = shade.df,
      width = dens.tile.width
    ) +
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    geom_rect(
      aes(xmin = x.limits[1], xmax = x.limits[2], ymin = ymin, ymax = ymax),
      data = box.df,
      fill = NA,
      colour = "gray20"
    ) +
    geom_segment(
      aes(x = LFC, y = y, xend = LFC, yend = yend),
      colour = colour,
      data = df
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
      breaks = row.centre(seq_len(n.genes)),
      labels = genes,
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = x.limits) +
    labs(x = "Log2(Fold change)", y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )

  plot_grid(
    p.density,
    p.ranks,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(3, length(genes))
  )
}

# select top x genes for enrichment/depletion, ordered by rank
df.enriched <- data %>% filter(pos.rank <= select) %>% arrange(pos.rank)
df.depleted <- data %>% filter(neg.rank <= select) %>% arrange(neg.rank)

# reversed so the top-ranked gene (rank 1) is drawn at the top of the plot
genes.enriched <- rev(unique(df.enriched$Gene))
genes.depleted <- rev(unique(df.depleted$Gene))

p.pos <- build_sgrank_plot(df.enriched, genes.enriched, "#e41a1c")
p.neg <- build_sgrank_plot(df.depleted, genes.depleted, "#377eb8")

# save plots to file
ggsave(
  plot = p.pos,
  filename = snakemake@output[["pos"]],
  units = "in",
  width = 10,
  height = 2 + length(genes.enriched) * 0.6
)
ggsave(
  plot = p.neg,
  filename = snakemake@output[["neg"]],
  units = "in",
  width = 10,
  height = 2 + length(genes.depleted) * 0.6
)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
