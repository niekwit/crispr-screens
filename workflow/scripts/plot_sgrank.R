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
dens.df <- tibble(x = dens$x, y = dens$y)
max.y <- max(dens.df$y)

# gradient strip sits just below the density curve's baseline
strip.height <- max.y * 0.18
strip.y <- -strip.height / 2

# density panel: curve on top, density-shaded gradient strip below the baseline
p.density <- ggplot(dens.df, aes(x = x)) +
  geom_area(aes(y = y), fill = "grey88", colour = "grey20", linewidth = 0.3) +
  geom_tile(aes(y = strip.y, height = strip.height, fill = y)) +
  scale_fill_gradient(low = "white", high = "black", guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(-strip.height, max.y * 1.05),
    breaks = scales::breaks_pretty(n = 3)(c(0, max.y)),
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

  # rectangle background per gene row
  bgcol <- tibble(
    id = rep(seq(1, max(df$index)), each = 4),
    x = rep(
      c(x.limits[1], x.limits[2], x.limits[2], x.limits[1]),
      max(df$index)
    ),
    y = unlist(lapply(seq(1, max(df$index)), function(i) {
      c(
        (interval + binwidth) * (i - 1),
        (interval + binwidth) * (i - 1),
        (interval + binwidth) * i - interval,
        (interval + binwidth) * i - interval
      )
    }))
  )

  p.ranks <- ggplot() +
    geom_polygon(
      aes(x = x, y = y, group = id),
      fill = "#dedede",
      colour = "gray20",
      data = bgcol
    ) +
    geom_segment(
      aes(x = LFC, y = y, xend = LFC, yend = yend),
      colour = colour,
      data = df
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
      breaks = bgcol$y[seq(1, nrow(bgcol), 4)] + binwidth / 2,
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
