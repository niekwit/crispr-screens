#redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)

# create correlation matrix
# https://stackoverflow.com/a/22282852/11329736
norm.counts <- read.csv(snakemake@input[[1]]) %>%
  select(-c("sgRNA","gene"))%>% 
  as.matrix %>%
  cor(method = "spearman") %>%
  as.data.frame %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1)

# plot heatmap
p <- ggplot(norm.counts, aes(var1, var2, fill = value)) + 
  geom_tile() + 
  scale_fill_viridis_c(guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = NULL) +
  ggtitle("Spearman correlation\nof normalised counts") +
  theme_bw(14) + 
  theme(aspect.ratio=1,
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.title = element_text(hjust = 0.5)) 

#save plot
ggsave(snakemake@output[[1]], p)


# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")


