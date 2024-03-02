# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# get data
fasta <- snakemake@params[["fasta"]]
data <- read.delim(snakemake@input[[1]]) %>% # count table
  dplyr::select(-c(1,2))

# create df to store coverage
# remove prepended X from samples names (only happens when they start with a number)
df <- as.data.frame(matrix(data = NA,
                           ncol = 2,
                           nrow = ncol(data))) %>%
  rename(sample = 1,
         coverage = 2) %>%
  mutate(sample = str_remove(names(data), "^X")) 

# get number of sgRNAs from fasta file (library size)
sgrnas <- length(readLines(fasta)) / 2

# calculate sequence coverage for each sample
for (n in seq(names(data))){
  tmp <- data %>%
    dplyr::select(n)
  
  sum.reads <- sum(tmp)
  coverage <- sum.reads / sgrnas
  
  df[n,"coverage"] <- coverage
}

# plot sequence coverage
p <- ggplot(df, aes(x = sample, y = coverage)) +
  geom_bar(stat = "identity", 
           fill = "aquamarine4",
           colour = "black") +
  theme_cowplot(16) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ylab("Fold sequence coverage") +
  xlab(NULL)

# save plot
ggsave(snakemake@output[[1]], p)


# close log file
sink(log, type = "output")
sink(log, type = "message")

