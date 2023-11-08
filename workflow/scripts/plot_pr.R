# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# load libraries
library(ggplot2)
library(cowplot)

# load data
data <- read.delim(snakemake@input[[1]])

# plot Precision and Recall from data
p <- ggplot(data, aes(x = Recall, y = Precision)) +
  geom_line() +
  labs(x = "Recall", 
       y = "Precision") +
  theme_cowplot(16) +
  scale_y_continuous(limits = c(0, 1))

# save plot
ggsave(snakemake@output[[1]], p)


# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

