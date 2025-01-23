# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(CRISPRcleanR)

# Load Snakemake variable
counts <- read.table(snakemake@input[["counts"]], header = TRUE)
control.name <- snakemake@params[["control"]]
test.name <- snakemake@params[["test"]]
library.name <- snakemake@params[["lib_name"]]
comparison <- snakemake@wildcards[["comparison"]]
corrected.fc.file <- snakemake@output[["corr_lfc"]]
corrected.counts.file <- snakemake@output[["corr_counts"]]
ceg <- snakemake@params[["ceg"]]
cneg <- snakemake@params[["cneg"]]

print(paste("Running CRISPRcleanR for", comparison))

# Create output directory
print("Creating output directory...")
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
  select(1, 2, all_of(control.columns), all_of(non.control.columns))

# Check if library name is in available libraries
available.libs <- c("AVANA_Library", "Brunello_Library", "GeCKO_Library_v2", "KY_Library_v1.0",
                    "KY_Library_v1.1", "MiniLibCas9_Library", "Whitehead_Library")
if (!library.name %in% available.libs) {
  print(paste("Loading sgRNA information from", snakemake@params[["lib"]]))
  # Library data not available so load from provided file
  full.annotations <- read.csv(snakemake@params[["lib"]])
  
  # Rown names should be same as CODE column
  rownames(full.annotations) <- full.annotations$CODE
  
  # Check if required columns are present in library file
  required.columns <- c("seq", "GENES", "CODE", "CHRM", "STARTpos", "ENDpos", "STRAND")
  if (!all(required.columns %in% colnames(full.annotations))) {
    stop("Library file does not contain all required columns: ", 
         paste(setdiff(required.columns, colnames(full.annotations)), collapse = ", "))
  }
} else {
  print(paste("Loading sgRNA information from CRISPRcleanR:", library.name))
  # Load library data from CRISPRcleanR package
  data(list = library.name)
  full.annotations <- get(library.name)
  rm(list = library.name)
}

# Remove chr from chromosome name
full.annotations$CHRM <- gsub("chr", "", full.annotations$CHRM)

# Check if any non-existing chromosome exist in annotations
# This would happen with sgRNAs targeting non-genomic sequences, eg. EGFP
non.real.chr <- setdiff(unique(full.annotations$CHRM), c(1:24,"X","Y"))
non.real.chr.count <- length(non.real.chr)
if (non.real.chr.count > 0) {
  # Convert any non-real chromosome to an integer
  # Assign these chromosomes a number higher than 24
  # This is to avoid confusing CRISPRcleanR
  print(paste("Non-existing chromosomes found in annotations:", paste(non.real.chr, collapse = ", ")))
  print("Assigning them a chromosome number higher than 24")
  new.chr.names <- 25:(24 + non.real.chr.count)
  for (i in 1:non.real.chr.count) {
    print(paste("Converting", non.real.chr[i], "to", new.chr.names[i]))
    full.annotations$CHRM <- gsub(non.real.chr[i], new.chr.names[i], full.annotations$CHRM)
  }
}

# Determine how many control samples are in the count table
ncontrols <- length(control.columns)

# Perform normalisation of raw counts and compute sgRNA log fold-changes
print(paste("Normalising raw counts and computing sgRNA log fold-changes against ", control.name), "...")
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
print("Sorting sgRNAs based on there chromosomal location...")
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs, full.annotations)

print("Correcting fold changes based on chromosomal position...")
# Correct biases due to gene independent responses to CRISPR-Cas9 targeting
correctedFCs<-ccr.GWclean(gwSortedFCs, 
                          display = FALSE, 
                          label = comparison)

# Correct counts based on corrected fold changes (new input file for MAGeCK)
print("Correcting counts based on corrected fold changes...")
correctedCounts <- ccr.correctCounts(comparison,
                                     normANDfcs$norm_counts,
                                     correctedFCs,
                                     full.annotations,
                                     minTargetedGenes = 3,
                                     OutDir = out.dir)

# Count table already in MAGeCK/DrugZ format
print("Saving corrected count table...")
write.table(correctedCounts, 
            corrected.counts.file, 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE)

# Save corrected fold changes in BAGEL2 format
print("Creating BAGEL2 format LFC table...")
lfc.bagel <- correctedFCs$corrected_logFCs %>%
  mutate(CODE = rownames(.)) %>%
  select(CODE, correctedFC, -genes) %>%
  left_join(full.annotations %>% select(seq, CODE, GENES), by = "CODE") %>%
  rename_with(~ "GENE", GENES, .cols = GENES) %>%
  rename_with(~ "REAGENT_ID", seq, .cols = seq) %>%
  rename_with(~ test.name, everything(), .cols = correctedFC) %>%
  select(-CODE) %>%
  select(REAGENT_ID, GENE, !!test.name)

print("Saving BAGEL2 format LFC table...")
write.table(lfc.bagel, 
            corrected.fc.file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# Prepare data for QC plots
if (ceg == "none") {
  data(BAGEL_essential)
} else {
  BAGEL_essential <- scan(ceg, what = "character")
}
if (cneg == "none") {
  data(BAGEL_nonEssential)
} else {
  BAGEL_nonEssential <- scan(cneg, what = "character")
}

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