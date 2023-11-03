#redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

#load required libraries
library(tidyverse)
library(viridis)
library(MAGeCKFlute)
library(clusterProfiler)
library(yaml)

# path to gene summary file 
file1 <- snakemake@input[[1]]

# path to sgRNA summary file 
file2 <- snakemake@input[[2]]

# get fdr cutoff
fdr <- snakemake@params["fdr"]

# load data
gdata <- ReadRRA(file1)
gdata$LogFDR <- -log10(gdata$FDR)

# load plot settings
settings <- read_yaml("workflow/envs/plot_settings.yaml",readLines.warn=FALSE)
theme <- paste0(settings["ggplot2_theme"][[1]],"(base_size = ",settings["font_size"][[1]],")")

# volcano plot
pdf(snakemake@output[[1]])
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id")
dev.off()

# dot plots
gdata$RandomIndex <- sample(1:nrow(gdata), nrow(gdata))
gdata <- gdata[order(-gdata$Score), ]

gg <- gdata[gdata$Score > 0, ]
pdf(snakemake@output[[2]])
ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,2), 
                 groups = "top", top = 10, ylab = "Log2FC")
dev.off()

gg <- gdata[gdata$Score < 0, ]
pdf(snakemake@output[[3]])
ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,2), 
                 groups = "bottom", top = 10, ylab = "Log2FC")
dev.off()

# sgRNA rank plot
sdata <- ReadsgRRA(file2)

pdf(snakemake@output[[4]])
sgRankView(sdata, top = 5, bottom = 5)
dev.off()

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")



