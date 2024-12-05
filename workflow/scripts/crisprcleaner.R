# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(CRISPRcleanR)

# Load Snakemake variable
counts <- read.table(snakemake@input[["counts"]], header = TRUE)
fasta <- snakemake@params[["fasta"]]
control.name <- snakemake@params[["control"]]
test.name <- snakemake@params[["test"]]
library.annotation <- read.table(snakemake@input[["annotation"]], header = TRUE)
comparison <- snakemake@wildcards[["mcomparison"]]
cell.line <- snakemake@params[["cell_line"]]
corrected.fc.file <- snakemake@output[["corr_fc"]]
corrected.counts <- snakemake@output[["corr_counts"]]
out.dir <- dirname(corrected.fc.file)
dir.create(out.dir, showWarnings = FALSE)

# Remove duplicated rows in counts (this can cause BAGEL2 to fail)
counts <- counts %>% distinct()

# Only keep test and control columns in count table
# REASON: the average LFCs will be calculated later of all count columns that are not control samples
all.samples <- c(str_split(test.name, ",")[[1]], str_split(control.name, ",")[[1]])
counts <- counts %>%
  select(1, 2, all_of(all.samples))

# Get column number(s) of control(s)
if (str_detect(control.name, ",")) {
  control.columns <- which(colnames(counts) %in% strsplit(control.name, ",")[[1]])
  } else {
  control.columns <- which(colnames(counts) == control.name)
}

# Get all columns numbers except control(s) starting from 3rd column
non.control.columns <- setdiff(3:ncol(counts), control.columns)

# Move control column(s) after 2nd column
counts <- counts %>%
  select(1, 2, control.columns, all_of(non.control.columns))

# Load fasta file and create library annotation file
sequences <- read.table(fasta, 
                        header = FALSE,
                        comment.char = ">")
names <- read.table(fasta, header = FALSE) %>%
  # Remove lines not starting with ">"
  filter(grepl("^>", V1)) %>%
  # Remove ">" from the beginning of the line
  mutate(V1 = gsub("^>", "", V1))
  
full.annotations <- data.frame(CODE = names$V1, seq = sequences$V1) %>%
  # Create GENE column by removing everything after the first underscore
  mutate(GENES = gsub("_.*", "", CODE)) %>%
  left_join(library.annotation, by = c("seq", "GENES")) %>%
  # Row names as CODE
  column_to_rownames(var = "CODE") %>%
  mutate(CODE = rownames(.)) %>%
  # Remove chr from chromosome name
  mutate(CHRM = gsub("chr", "", CHRM))

# Determine how many control samples are in the count table
ncontrols <- length(control.columns)

# Perform normalisation of raw counts and compute sgRNA log fold-changes
normANDfcs <- ccr.NormfoldChanges(Dframe=counts,
                                  min_reads = 30,
                                  EXPname = comparison,
                                  libraryAnnotation = full.annotations,
                                  display = FALSE,
                                  ncontrols = ncontrols,
                                  outdir = out.dir)

# Remove sgRNAs from full.annotations that are not in normANDfcs$logFCs$sgRNA
full.annotations <- full.annotations[rownames(full.annotations) %in% normANDfcs$logFCs$sgRNA,] 

# Correct fold changes according to position
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs, full.annotations)

# Correct biases due to gene independent responses to CRISPR-Cas9 targeting
correctedFCs<-ccr.GWclean(gwSortedFCs, 
                          display = FALSE, 
                          label = comparison)

# Correct counts based on corrected fold changes (new input file for MAGeCK)
correctedCounts <- ccr.correctCounts(comparison,
                                     normANDfcs$norm_counts,
                                     correctedFCs,
                                     full.annotations,
                                     minTargetedGenes = 3,
                                     OutDir = out.dir)

# Count table already in MAGeCK format
write.table(correctedFCs$corrected_logFCs, 
            corrected.counts, 
            sep = "\t", 
            quote = FALSE)

# Save corrected fold changes in BAGEL2 format
lfc.bagel <- correctedFCs$corrected_logFCs %>%
  # column CODE with rownames
  mutate(CODE = rownames(.)) %>%
  # keep CODE, corrrectedFC, genes columns
  select(CODE, correctedFC, genes) %>%
  # Add seq column from full.annotations with join
  left_join(full.annotations %>% select(seq, CODE), by = "CODE") %>%
  # Rename columns: seq to REAGENT_ID, correctedFC to test.name variable
  rename(REAGENT_ID = seq,
         GENE = genes,
         !!test.name := correctedFC) %>%
  # Remove CODE column and reorder other columns: REAGENT_ID, GENE, test.name
  select(-CODE) %>%
  select(REAGENT_ID, GENE, !!test.name)


# Prepare data for QC plots
FCs <- correctedFCs$corrected_logFCs$avgFC
names(FCs) <- rownames(correctedFCs$corrected_logFCs)
BAGEL_essential_sgRNAs <- ccr.genes2sgRNAs(full.annotations, BAGEL_essential)
BAGEL_nonEssential_sgRNAs <- ccr.genes2sgRNAs(full.annotations, BAGEL_nonEssential)

### ROC curves
# Calculate ROC curve data for gene and sgRNA level
sgRNA_level_ROC <- ccr.ROC_Curve(FCs,
                                 BAGEL_essential_sgRNAs,
                                 BAGEL_nonEssential_sgRNAs,
                                 display = FALSE) 
df.sgRNA <- as.data.frame(sgRNA_level_ROC$curve) %>%
  mutate(class = "sgRNA")

geneFCs <- ccr.geneMeanFCs(FCs, full.annotations)
gene_level_ROC <- ccr.ROC_Curve(geneFCs,
                                BAGEL_essential,
                                BAGEL_nonEssential,
                                FDRth = 0.05,
                                display = FALSE)
df.gene <- as.data.frame(gene_level_ROC$curve) %>%
  mutate(class = "gene")

# Combine data for gene and sgRNA level ROC curves to plot
df <- rbind(df.sgRNA, df.gene)

roc.plot <- ggplot(df, aes(x = 1 - specificity, 
                       y = sensitivity, 
                       color = class)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = paste0("ROC curve: gene and sgRNA level\n(", comparison, ")"),
       x = "TNR",
       y = "Recall") +
  theme_cowplot(18) +
  annotate("text", 
           x = 0.5, 
           y = 0.05,
           size = 6,
           label = paste("AUC = ", round(gene_level_ROC$AUC, 3), 
                         " (gene), ", 
                         round(sgRNA_level_ROC$AUC, 3), 
                         " (sgRNA)")) +
  theme(legend.position = c(0.85, 0.1),
        legend.title = element_blank())

ggsave(filename = snakemake@output[["roc"]],
       roc.plot, 
       width = 8, 
       height = 6)

### PR curves
# Calculate PR curve data for gene and sgRNA level
sgRNA_level_PR <- ccr.PrRc_Curve(FCs,
                               BAGEL_essential_sgRNAs,
                               BAGEL_nonEssential_sgRNAs,
                               display = FALSE)
df.sgRNA <- as.data.frame(sgRNA_level_PR$curve) %>%
  mutate(class = "sgRNA")

gene_level_PR <- ccr.PrRc_Curve(geneFCs,
                             BAGEL_essential,
                             BAGEL_nonEssential,
                             display = FALSE,
                             FDRth = 0.05)
df.gene <- as.data.frame(gene_level_PR$curve) %>%
  mutate(class = "gene")

# Combine data for gene and sgRNA level PR curves to plot
df <- rbind(df.sgRNA, df.gene)

pr.plot <- ggplot(df, aes(x = recall, 
                      y = precision, 
                      color = class)) +
  geom_line(linewidth = 1) +
  labs(title = paste0("PR curve: gene and sgRNA level\n(", comparison, ")"),
       x = "Recall",
       y = "Precision") +
  theme_cowplot(18) +
  annotate("text", 
           x = 0.35, 
           y = 0.25,
           size = 6,
           label = paste("AUC = ", round(gene_level_PR$AUC, 3), 
                         " (gene), ", 
                         round(sgRNA_level_PR$AUC, 3), 
                         " (sgRNA)")) +
  theme(legend.position = c(0.025, 0.15),
        legend.title = element_blank())

ggsave(snakemake@output[["pr"]], 
       pr.plot,
       width = 8, 
       height = 6)

### Depletion profile of gene signatures
# Load gene sets
data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)

SIGNATURES <- list (Ribosomal_Proteins=EssGenes.ribosomalProteins,
                   DNA_Replication = EssGenes.DNA_REPLICATION_cons,
                   RNA_polymerase = EssGenes.KEGG_rna_polymerase,
                   Proteasome = EssGenes.PROTEASOME_cons,
                   Spliceosome = EssGenes.SPLICEOSOME_cons,
                   Common_essential = BAGEL_essential,
                   Non_essential = BAGEL_nonEssential)

pdf(snakemake@output[["drnk"]])
Recall_scores <- ccr.VisDepAndSig(FCsprofile = geneFCs,
                                  SIGNATURES = SIGNATURES,
                                  TITLE = comparison,
                                  pIs = 6,
                                  nIs = 7)
dev.off()

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")