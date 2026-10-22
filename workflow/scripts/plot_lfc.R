# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# load required libraries
library(tidyverse)
library(viridis)
library(ggrepel)
library(cowplot)

# load data
data <- read.delim(snakemake@input[[1]])

# dot plot function
dot_plot <- function(df, df_label, filename) {
  # create plot
  .data <- rlang::.data
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$x, y = .data$log2fc, fill = .data$log.p.value)
  ) +
    ggplot2::geom_point(size = 4, shape = 21) +
    ggplot2::labs(
      x = "Random Index",
      y = "Log2(Fold Change)",
      fill = "-log10(p value)"
    ) +
    cowplot::theme_cowplot(16) +
    ggrepel::geom_text_repel(
      data = df_label,
      ggplot2::aes(x = .data$x, y = .data$pos.lfc, label = .data$id)
    ) +
    viridis::scale_fill_viridis(
      guide = ggplot2::guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black"
      )
    )
  # save to file
  ggplot2::ggsave(filename, p, width = 12, height = 10)
}

# enriched genes
enriched <- data |>
  filter(pos.lfc > 0) |>
  mutate(
    log2fc = pos.lfc,
    log.p.value = -log10(pos.p.value),
    x = sample.int(n(), n())
  ) |>
  arrange(x)

df_label <- enriched |>
  arrange(desc(log2fc)) |>
  slice(1:10)

dot_plot(enriched, df_label, snakemake@output[["pos"]])

# depleted genes
depleted <- data |>
  filter(neg.lfc < 0) |>
  mutate(
    log2fc = neg.lfc,
    log.p.value = -log10(neg.p.value),
    x = sample.int(n(), n())
  ) |>
  arrange(x)

df_label <- depleted |>
  arrange(log2fc) |>
  slice(1:10)

dot_plot(depleted, df_label, snakemake@output[["neg"]])


# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
