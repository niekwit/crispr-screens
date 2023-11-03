# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# get HISAT2 alignment rates (in log files)
files <- snakemake@input

# create df for storing alignment rates
df <- as.data.frame(matrix(ncol = 2, nrow = 0))
names(df) <- c("sample","alignment.rate")

# get sample name and mapping rates from log files
for (i in seq(files)){
  sample <- system(paste0("echo ", files[i], "| sed 's/.log//'"), intern=TRUE)
  sample <- basename(sample)
  
  rate <- system(paste0('grep "aligned exactly 1 time" ', files[i], " | awk '{print $2}' | sed 's/(*[%)]*//g'"), intern=TRUE)
  rate <- as.numeric(rate)

  # add to df
  df[i,"sample"] <- sample
  df[i,"mapping.rate"] <- rate
}

# remove prepended X from samples names (only happens when they start with a number)
df$sample <- str_remove(df$sample, "^X")

# create plot
p <- ggplot(df, aes(x = sample, y = mapping.rate)) +
  geom_bar(stat = "identity", 
           fill = "aquamarine4",
           colour = "black") +
  theme_cowplot(16) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("Mapping rate (%)") +
  xlab(NULL)

# save plot
ggsave(snakemake@output[[1]], p)

# close log file
sink(log, type = "output")
sink(log, type = "message")
