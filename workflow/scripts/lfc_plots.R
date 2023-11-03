# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
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
dot_plot <- function(df, df.label, filename) {
  # create plot
  p <- ggplot(df, aes(x = x, y = log2fc, fill=log.p.value),
              ) +
    geom_point(size=4,
               shape=21) +
    labs(x="Random Index", 
         y="Log2(Fold Change)",
         fill="-log10(p value)") +
    theme_cowplot(16) +
    geom_text_repel(data=df.label,
                    aes(x=`x`,
                        y=`pos.lfc`,
                        label=`id`)) +
    scale_fill_viridis(guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black"))
  # save to file
   ggsave(filename, p, width = 12, height = 10)
}

# enriched genes
enriched <- data %>% 
  filter(pos.lfc > 0) %>%
  mutate(log2fc = pos.lfc,
         log.p.value = -log10(pos.p.value),
         x = sample.int(nrow(.), nrow(.))) %>%
  arrange(x) 

df.label <- enriched %>% 
  arrange(desc(log2fc)) %>% 
  slice(1:10)

dot_plot(enriched, df.label, snakemake@output[["pos"]])

# depleted genes
depleted <- data %>% 
  filter(neg.lfc < 0) %>%
  mutate(log2fc = neg.lfc,
         log.p.value = -log10(neg.p.value),
         x = sample.int(nrow(.), nrow(.))) %>%
  arrange(x)

df.label <- depleted %>% 
  arrange(log2fc) %>% 
  slice(1:10)

dot_plot(depleted, df.label, snakemake@output[["neg"]])


# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

