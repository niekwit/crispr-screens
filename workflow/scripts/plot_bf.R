# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# load libraries
library(ggplot2)
library(cowplot)

# load data
data <- read.delim(snakemake@input[[1]])

# plot histogram for BF values in data
p <- ggplot(data, aes(x = BF)) +
  geom_histogram(bins = 100,
                 fill = "aquamarine4") +
  labs(x = "Bayes Factor", 
       y = "Number of genes") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot(18)

# save plot
ggsave(snakemake@output[[1]], p)


# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

