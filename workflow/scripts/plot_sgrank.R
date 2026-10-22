# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Inspired by
# https://github.com/WubingZhang/MAGeCKFlute/blob/master/R/sgRankView.R
# https://www.jci.org/articles/view/201639 Figure 1D

# load required libraries
library(tidyverse)
library(cowplot)

# load data
data <- read.delim(snakemake@input[["sg"]])

# load gene rank data
rank <- read.delim(snakemake@input[["gene"]]) |>
  dplyr::select(id, neg.rank, pos.rank) |>
  rename(Gene = id)

# Add Gene rank to data
data <- data |> left_join(rank, by = "Gene")

# set parameters
select <- 5 # number of genes to plot for enriched and depleted
binwidth <- 0.3
interval <- 0.1

# shared x-axis range across all panels, based on the full sgRNA LFC
# distribution
x_range <- range(data$LFC, na.rm = TRUE)
x_pad <- diff(x_range) * 0.02
x_limits <- c(x_range[1] - x_pad, x_range[2] + x_pad)

# density of the full sgRNA LFC distribution (background reference)
dens <- density(
  data$LFC,
  na.rm = TRUE,
  n = 512,
  from = x_limits[1],
  to = x_limits[2]
)
dens_df <- tibble(x = dens$x, dens = dens$y)
max.dens <- max(dens_df$dens)
# slightly overlapping tile width avoids hairline rendering seams between
# adjacent density tiles in the gene-row shading
dens_tile_width <- diff(range(dens_df$x)) / (nrow(dens_df) - 1) * 1.02

# density panel: just the curve, no separate gradient strip (the density
# shading now lives inside the gene rows themselves, see build_sgrank_plot)
p_density <- ggplot(dens_df, aes(x = x, y = dens)) +
  geom_area(fill = "grey88", colour = "grey20", linewidth = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, max.dens * 1.05),
    breaks = scales::breaks_pretty(n = 3)(c(0, max.dens)),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = x_limits) +
  labs(y = "Density") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black", linewidth = 0.3),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# builds the gene-rank row panel for a set of top genes and stacks it under
# the shared density panel (colour: red for enriched genes, blue for depleted)
build_sgrank_plot <- function(df, genes, colour) {
  .data <- rlang::.data
  n_genes <- length(genes)
  row_centre <- function(i) (interval + binwidth) * (i - 1) + binwidth / 2

  # as.integer(Gene) (not row position) so each sgRNA lines up with its
  # gene's row regardless of the row order genes/df happen to be in
  df$Gene <- factor(df$Gene, levels = genes)
  df$index <- as.integer(df$Gene)
  df$y <- (binwidth + interval) * (df$index - 1)
  df$yend <- (binwidth + interval) * df$index - interval
  df <- as.data.frame(df[, c("sgrna", "Gene", "LFC", "y", "yend", "index")])

  # density-shaded background per gene row (white = low density, black =
  # high density), echoing the overall LFC distribution behind each gene
  shade_df <- do.call(
    rbind,
    lapply(seq_len(n_genes), function(i) {
      dplyr::mutate(dens_df, y = row_centre(i))
    })
  )

  # row outline, drawn on top of the shading with no fill of its own
  row_index <- seq_len(n_genes)
  box_df <- tibble::tibble(
    index = row_index,
    ymin = (interval + binwidth) * (row_index - 1),
    ymax = (interval + binwidth) * row_index - interval
  )

  p_ranks <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        height = binwidth,
        fill = .data$dens
      ),
      data = shade_df,
      width = dens_tile_width
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "black",
      guide = "none"
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = x_limits[1],
        xmax = x_limits[2],
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      data = box_df,
      fill = NA,
      colour = "gray20"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$LFC,
        y = .data$y,
        xend = .data$LFC,
        yend = .data$yend
      ),
      colour = colour,
      data = df
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = row_centre(seq_len(n_genes)),
      labels = genes,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = x_limits) +
    ggplot2::labs(x = "Log2(Fold change)", y = NULL) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )

  plot_grid(
    p_density,
    p_ranks,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(3, length(genes))
  )
}

# select top x genes for enrichment/depletion, ordered by rank
df_enriched <- data |> filter(pos.rank <= select) |> arrange(pos.rank)
df_depleted <- data |> filter(neg.rank <= select) |> arrange(neg.rank)

# reversed so the top-ranked gene (rank 1) is drawn at the top of the plot
genes_enriched <- rev(unique(df_enriched$Gene))
genes_depleted <- rev(unique(df_depleted$Gene))

p_pos <- build_sgrank_plot(df_enriched, genes_enriched, "#e41a1c")
p_neg <- build_sgrank_plot(df_depleted, genes_depleted, "#377eb8")

# save plots to file
ggsave(
  plot = p_pos,
  filename = snakemake@output[["pos"]],
  units = "in",
  width = 10,
  height = 2 + length(genes_enriched) * 0.6
)
ggsave(
  plot = p_neg,
  filename = snakemake@output[["neg"]],
  units = "in",
  width = 10,
  height = 2 + length(genes_depleted) * 0.6
)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
