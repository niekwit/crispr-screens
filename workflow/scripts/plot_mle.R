# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# load required libraries
library(tidyverse)
library(viridis)
library(ggrepel)
library(cowplot)

# load data (check.names = FALSE keeps the "|" in e.g. "condition|beta")
data <- read.delim(snakemake@input[[1]], check.names = FALSE)

# MAGeCK mle reports one <condition>|beta, <condition>|z, <condition>|p-value,
# <condition>|fdr, <condition>|wald-p-value, <condition>|wald-fdr block per
# non-baseline column of the design matrix (the baseline/intercept column
# itself is never reported), so detect condition names from the "|beta"
# columns instead of hard-coding them
conditions <- sub(
  "\\|beta$",
  "",
  grep("\\|beta$", colnames(data), value = TRUE)
)

# reshape into one row per gene x condition, with a random index per
# condition to spread out points (same approach as plot_lfc.R)
long <- map_dfr(conditions, function(cond) {
  tibble(
    condition = cond,
    Gene = data[["Gene"]],
    beta = data[[paste0(cond, "|beta")]],
    log.fdr = -log10(data[[paste0(cond, "|fdr")]]),
    x = sample.int(nrow(data), nrow(data))
  )
}) %>%
  arrange(condition, x)

# a gene can get fdr = 0 (common with few permutation rounds), which turns
# -log10(fdr) into Inf and breaks the colour scale; clip it to the largest
# finite value instead of dropping/greying out the point
finite.max <- max(long$log.fdr[is.finite(long$log.fdr)], na.rm = TRUE)
long <- long %>%
  mutate(log.fdr = ifelse(is.finite(log.fdr), log.fdr, finite.max))

# label the top 5 most enriched and top 5 most depleted genes per condition
df.label <- long %>%
  group_by(condition) %>%
  group_modify(
    ~ bind_rows(
      slice_max(.x, beta, n = 5),
      slice_min(.x, beta, n = 5)
    )
  ) %>%
  ungroup()

# create plot: one facet per design matrix condition
p <- ggplot(long, aes(x = x, y = beta, fill = log.fdr)) +
  geom_point(size = 4, shape = 21) +
  geom_text_repel(data = df.label, aes(x = x, y = beta, label = Gene)) +
  labs(
    x = "Random Index",
    y = "Beta score",
    fill = "-log10(FDR)"
  ) +
  facet_wrap(~condition, scales = "free") +
  theme_cowplot(16) +
  scale_fill_viridis(
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )

# scale canvas size to the number of conditions (facets)
n.cond <- length(conditions)
n.col <- min(3, n.cond)
n.row <- ceiling(n.cond / 3)

# save to file
ggsave(snakemake@output[[1]], p, width = 6 * n.col, height = 5 * n.row)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
