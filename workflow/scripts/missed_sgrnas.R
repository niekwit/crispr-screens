# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(reshape2)
library(cowplot)

# load count table
counts <- read.delim(snakemake@input[[1]]) %>%
  dplyr::select(-c("sgRNA", "gene"))

# get number of sgRNAs with zero counts
df <- colSums(counts == 0) %>%
  melt() %>%
  rownames_to_column(var = "sample") 

# create bar graph
p <- ggplot(data = df, aes(x = sample, y = value)) + 
  geom_bar(stat = "identity",
           fill = "aquamarine4",
           colour = "black") +
  theme_cowplot(16) +
  xlab(NULL) +
  ylab("Missed sgRNAs") +
  scale_x_discrete(guide = guide_axis(angle = 45))

# save plot
ggsave(snakemake@output[[1]], p)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")




