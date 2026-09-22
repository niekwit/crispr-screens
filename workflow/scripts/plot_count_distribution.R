# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# fraction of sgRNAs below this fraction of the sample mean is called "depleted"
low.fraction <- 0.1
# samples with more than this fraction of depleted sgRNAs are flagged
flag.threshold <- 0.1

# load count table
counts <- read.delim(snakemake@input[[1]], check.names = FALSE) %>%
  dplyr::select(-c(1, 2))

# per sample summary
summary.df <- counts %>%
  pivot_longer(everything(), names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  summarise(
    mean = mean(count),
    median = median(count),
    median_over_mean = median / mean,
    fraction_low = mean(count < low.fraction * mean),
    fraction_zero = mean(count == 0),
    .groups = "drop"
  ) %>%
  mutate(
    flagged = fraction_low > flag.threshold,
    sample = factor(sample, levels = names(counts)),
    label = sprintf(
      "%s\n%.0f%% sgRNAs < %.0f%% of mean | median/mean %.2f",
      sample,
      100 * fraction_low,
      100 * low.fraction,
      median_over_mean
    ),
    label = factor(label, levels = unique(label[order(sample)]))
  )

write.csv(
  summary.df %>% dplyr::select(-label),
  snakemake@output[["csv"]],
  row.names = FALSE
)

# log10 counts (+1 to keep zeros)
df <- counts %>%
  pivot_longer(everything(), names_to = "sample", values_to = "count") %>%
  mutate(
    sample = factor(sample, levels = names(counts)),
    log.count = log10(count + 1)
  ) %>%
  left_join(summary.df %>% dplyr::select(sample, label), by = "sample")

# vertical lines: median (solid), and low count cut-off (dashed)
lines.df <- summary.df %>%
  mutate(
    low.cutoff = log10(low.fraction * mean + 1),
    median = log10(median + 1)
  )

ncol.facet <- min(4, nlevels(df$sample))
flagged.df <- dplyr::filter(summary.df, flagged)

p <- ggplot(df, aes(x = log.count))

# A zero-row layer under facet_wrap crashes gtable rendering in ggplot2
# 3.5.2, so only add the flagged-sample highlight when something is flagged
if (nrow(flagged.df) > 0) {
  p <- p +
    geom_rect(
      data = flagged.df,
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "#D55E00",
      alpha = 0.12,
      inherit.aes = FALSE
    )
}

p <- p +
  geom_histogram(
    bins = 60,
    fill = "#419179",
    colour = "black",
    linewidth = 0.1
  ) +
  geom_vline(
    data = lines.df,
    aes(xintercept = low.cutoff),
    linetype = "dashed",
    colour = "grey30"
  ) +
  geom_vline(
    data = lines.df,
    aes(xintercept = median),
    colour = "#0072B2",
    linewidth = 0.8
  ) +
  facet_wrap(~label, ncol = ncol.facet) +
  theme_cowplot(12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9)
  ) +
  xlab(expression(log[10](read ~ count + 1))) +
  ylab("Number of sgRNAs") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ggtitle(
    "Read count distribution per sample",
    subtitle = sprintf(
      "Blue: median. Dashed: %.0f%% of mean count. Shaded: more than %.0f%% of sgRNAs below dashed line",
      100 * low.fraction,
      100 * flag.threshold
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

# save plot
n.rows <- ceiling(nlevels(df$sample) / ncol.facet)
ggsave(
  snakemake@output[[1]],
  p,
  width = 3.2 * ncol.facet,
  height = 2.6 * n.rows + 0.8,
  limitsize = FALSE
)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
