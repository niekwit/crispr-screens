# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(CRISPRcleanR)

# Some CRISPRcleanR functions open a graphics device even with
# display = FALSE, which in a non-interactive session creates an empty
# Rplots.pdf in the working directory: send that to a null device instead
pdf(NULL)

# Load Snakemake variable
counts <- read.delim(snakemake@input[["counts"]], header = TRUE)
control_name <- snakemake@params[["control"]]
test_name <- snakemake@params[["test"]]
library_name <- snakemake@params[["lib_name"]]
comparison <- snakemake@wildcards[["comparison"]]
corrected_fc_file <- snakemake@output[["corr_lfc"]]
corrected_counts_file <- snakemake@output[["corr_counts"]]
ceg <- snakemake@params[["ceg"]]
cneg <- snakemake@params[["cneg"]]
min_reads <- snakemake@params[["min_reads"]]

print(paste("Running CRISPRcleanR for", comparison))

# Create output directory
print("Creating output directory...")
out_dir <- dirname(corrected_fc_file)
dir.create(out_dir, showWarnings = FALSE)

# ccr.NormfoldChanges()/ccr.correctCounts() build their .RData paths via
# paste0(outdir, EXPname, "_suffix.RData") — no path separator inserted —
# so outdir must end in "/" or the filename gets mashed into out_dir's own
# name (e.g. "results/crisprcleanrL15_vs_plasmid_..." instead of landing
# inside out_dir). Point them at their own subdirectory instead of
# out_dir itself, both to satisfy that trailing-slash requirement and to
# keep these package side-effect files out of out_dir's "real" outputs.
rdata_dir <- file.path(out_dir, "rdata")
dir.create(rdata_dir, showWarnings = FALSE)
rdata_outdir <- paste0(rdata_dir, "/")

# Remove duplicated rows in counts (this can cause BAGEL2 to fail)
counts <- counts |> distinct()

# Only keep test and control columns in count table
# REASON: the average LFCs will be calculated later of all count
# columns that are not control samples
all.samples <- c(
  str_split(test_name, ",")[[1]],
  str_split(control_name, ",")[[1]]
)
counts <- counts |>
  select(1, 2, all_of(all.samples))

# Get column number(s) of control(s)
if (str_detect(control_name, ",")) {
  # multiple control columns
  control_columns <- which(
    colnames(counts) %in% strsplit(control_name, ",")[[1]]
  )
} else {
  control_columns <- which(colnames(counts) == control_name)
}

# Get all columns numbers except control(s) starting from 3rd column
non_control_columns <- setdiff(3:ncol(counts), control_columns)

# Move control column(s) after 2nd column
counts <- counts |>
  select(1, 2, all_of(control_columns), all_of(non_control_columns))

# Check if library name is in available libraries
available_libs <- c(
  "AVANA_Library",
  "Brunello_Library",
  "GeCKO_Library_v2",
  "KY_Library_v1.0",
  "KY_Library_v1.1",
  "MiniLibCas9_Library",
  "Whitehead_Library"
)
if (!library_name %in% available_libs) {
  print(paste("Loading sgRNA information from", snakemake@params[["lib"]]))
  # Library data not available so load from provided file
  full_annotations <- read.csv(snakemake@params[["lib"]])

  # Rown names should be same as CODE column
  rownames(full_annotations) <- full_annotations$CODE

  # Check if required columns are present in library file
  required_columns <- c(
    "seq",
    "GENES",
    "CODE",
    "CHRM",
    "STARTpos",
    "ENDpos",
    "STRAND"
  )
  if (!all(required_columns %in% colnames(full_annotations))) {
    stop(
      "Library file does not contain all required columns: ",
      paste(
        setdiff(required_columns, colnames(full_annotations)),
        collapse = ", "
      )
    )
  }
} else {
  print(paste("Loading sgRNA information from CRISPRcleanR:", library_name))
  # Load library data from CRISPRcleanR package
  data(list = library_name)
  full_annotations <- get(library_name)
  rm(list = library_name)
}

# Remove chr from chromosome name
full_annotations$CHRM <- gsub("chr", "", full_annotations$CHRM)

# Check if any non-existing chromosome exist in annotations
# This would happen with sgRNAs targeting non-genomic sequences (eg. EGFP),
# unplaced/alt genome assembly contigs (eg. KI270*, GL000*), or sgRNAs with
# no chromosomal location at all (CHRM is NA)
non_real_chr <- setdiff(unique(full_annotations$CHRM), c(1:24, "X", "Y"))
if (length(non_real_chr) > 0) {
  # Pool ALL non-real chromosomes/contigs into a single fake chromosome,
  # rather than giving each one its own number. Several of these contigs
  # (particularly unplaced/alt scaffolds) often have only a handful of
  # sgRNAs mapped to them, sometimes just one, and CRISPRcleanR's CBS-based
  # smoothing (ccr.GWclean -> smooth.CNA) crashes with "NA/NaN/Inf in
  # foreign function call" on a chromosome with a single data point, so
  # they need to be pooled together to have enough sgRNAs to smooth over.
  new_chr_name <- 25
  print(paste(
    "Non-existing chromosomes/contigs found in annotations:",
    paste(non_real_chr, collapse = ", ")
  ))
  print(paste("Assigning them all to chromosome", new_chr_name))
  full_annotations$CHRM[full_annotations$CHRM %in% non_real_chr] <-
    new_chr_name
}

# Check if gene and sgRNA names fully match between counts and annotations
# This is important for CRISPRcleanR to work properly
if (!all(full_annotations$GENES == counts$gene)) {
  stop("Gene names in counts and annotations do not match.")
}
if (!all(full_annotations$CODE == counts$sgRNA)) {
  stop("sgRNA names in counts and annotations do not match.")
}

# Check if there are any rows with NA in STARTpos/ENDpos in full_annotations
# (only these columns: NAs elsewhere, e.g. the optional EXONE column, are fine)
no_position <- is.na(full_annotations$STARTpos) |
  is.na(full_annotations$ENDpos)
if (any(no_position)) {
  print("NA values detected in annotation file for STARTpos or ENDpos.")
  print(
    "These might be non-targeting sgRNAs with no chromosomal location assigned."
  )
  print("Assigning random locations...")

  # Get rows with NA values and add some values
  na_rows <- full_annotations[no_position, ]
  na_rows$STARTpos <- seq(from = 10, to = nrow(na_rows) * 10, by = 10)
  na_rows$ENDpos <- na_rows$STARTpos + 10

  # Add gene name
  na_rows$GENES <- "CONTROL_GENE"

  # Remove rows with NA values from full_annotations
  full_annotations <- full_annotations[!no_position, ]

  # Combine the two data frames
  full_annotations <- rbind(full_annotations, na_rows)
}

# Check if there are any rows with NA in gene column in counts
# replace with CONTROL_GENE
if (any(is.na(counts$gene))) {
  print("NA values detected in counts file for gene column.")
  print("Replacing with CONTROL_GENE...")

  # Get rows with NA values and add some values
  na_rows <- counts[is.na(counts$gene), ]
  na_rows$gene <- "CONTROL_GENE"

  # Remove rows with NA values from counts
  counts <- counts[!is.na(counts$gene), ]

  # Combine the two data frames
  counts <- rbind(counts, na_rows)
}

# Determine how many control samples are in the count table
ncontrols <- length(control_columns)

# Perform normalisation of raw counts and compute sgRNA log fold-changes
print(paste(
  "Normalising raw counts and computing sgRNA log fold-changes against",
  control_name,
  "..."
))
print(paste(
  "Filtering sgRNAs for minimum number of reads in control:",
  min_reads
))
norm_and_fcs <- ccr.NormfoldChanges(
  Dframe = counts,
  min_reads = min_reads,
  EXPname = comparison,
  libraryAnnotation = full_annotations,
  display = FALSE,
  ncontrols = ncontrols,
  outdir = rdata_outdir
)

# Check if norm_and_fcs is empty
if (nrow(norm_and_fcs$logFCs) == 0) {
  stop(
    "ccr.NormfoldChanges command returned an empty data frame.\n
       Please check the library annotation file or min_reads value."
  )
}

# Remove sgRNAs from full_annotations that are not in norm_and_fcs$logFCs$sgRNA
full_annotations <- full_annotations[
  rownames(full_annotations) %in% norm_and_fcs$logFCs$sgRNA,
]

# Correct fold changes according to position
print("Sorting sgRNAs based on there chromosomal location...")
gw_sorted_fcs <- ccr.logFCs2chromPos(norm_and_fcs$logFCs, full_annotations)

print("Correcting fold changes based on chromosomal position...")
# Correct biases due to gene independent responses to CRISPR-Cas9 targeting
corrected_fcs <- ccr.GWclean(
  gw_sorted_fcs,
  display = FALSE,
  label = comparison
)

# Correct counts based on corrected fold changes (new input file for MAGeCK)
print("Correcting counts based on corrected fold changes...")
corrected_counts <- ccr.correctCounts(
  comparison,
  norm_and_fcs$norm_counts,
  corrected_fcs,
  full_annotations,
  minTargetedGenes = 3,
  OutDir = rdata_outdir
)

# Count table already in MAGeCK/DrugZ format
print("Saving corrected count table...")
write.table(
  corrected_counts,
  corrected_counts_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Save corrected fold changes in BAGEL2 format
print("Creating BAGEL2 format LFC table...")
corrected_log_fcs <- corrected_fcs$corrected_logFCs
annotations_subset <- full_annotations |> select(seq, CODE, GENES)
lfc_bagel <- corrected_log_fcs |>
  mutate(CODE = rownames(corrected_log_fcs)) |>
  select(CODE, correctedFC, -genes) |>
  left_join(annotations_subset, by = "CODE") |>
  rename_with(~"GENE", GENES, .cols = GENES) |>
  rename_with(~"REAGENT_ID", seq, .cols = seq) |>
  rename_with(~test_name, everything(), .cols = correctedFC) |>
  select(-CODE) |>
  select(REAGENT_ID, GENE, !!test_name)

print("Saving BAGEL2 format LFC table...")
write.table(
  lfc_bagel,
  corrected_fc_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Prepare data for QC plots
if (ceg == "none") {
  data(BAGEL_essential)
  bagel_essential <- BAGEL_essential
} else {
  bagel_essential <- scan(ceg, what = "character")
}
if (cneg == "none") {
  data(BAGEL_nonEssential)
  bagel_non_essential <- BAGEL_nonEssential
} else {
  bagel_non_essential <- scan(cneg, what = "character")
}

fcs <- corrected_fcs$corrected_logFCs$avgFC
names(fcs) <- rownames(corrected_fcs$corrected_logFCs)
bagel_essential_sgrnas <- ccr.genes2sgRNAs(full_annotations, bagel_essential)
bagel_non_essential_sgrnas <- ccr.genes2sgRNAs(
  full_annotations,
  bagel_non_essential
)

# These feed the QC report's interactive plots (workflow/scripts/qc/crisprcleanr.py)
# as plain data, not pre-rendered PDFs — Plotly draws them client-side the
# same way every other QC section does, so only the underlying curve/rank
# data needs exporting here, not a finished plot.

### ROC curves
# Calculate ROC curve data for gene and sgRNA level
sgrna_level_roc <- ccr.ROC_Curve(
  fcs,
  bagel_essential_sgrnas,
  bagel_non_essential_sgrnas,
  display = FALSE
)
df_sgrna <- as.data.frame(sgrna_level_roc$curve) |>
  mutate(class = "sgRNA", auc = sgrna_level_roc$AUC)

gene_fcs <- ccr.geneMeanFCs(fcs, full_annotations)
gene_level_roc <- ccr.ROC_Curve(
  gene_fcs,
  bagel_essential,
  bagel_non_essential,
  FDRth = 0.05,
  display = FALSE
)
df_gene <- as.data.frame(gene_level_roc$curve) |>
  mutate(class = "gene", auc = gene_level_roc$AUC)

print("Saving ROC curve data...")
write.csv(
  rbind(df_sgrna, df_gene),
  snakemake@output[["roc"]],
  row.names = FALSE
)

### PR curves
# Calculate PR curve data for gene and sgRNA level
sgrna_level_pr <- ccr.PrRc_Curve(
  fcs,
  bagel_essential_sgrnas,
  bagel_non_essential_sgrnas,
  display = FALSE
)
df_sgrna <- as.data.frame(sgrna_level_pr$curve) |>
  mutate(class = "sgRNA", auc = sgrna_level_pr$AUC)

gene_level_pr <- ccr.PrRc_Curve(
  gene_fcs,
  bagel_essential,
  bagel_non_essential,
  display = FALSE,
  FDRth = 0.05
)
df_gene <- as.data.frame(gene_level_pr$curve) |>
  mutate(class = "gene", auc = gene_level_pr$AUC)

print("Saving PR curve data...")
write.csv(rbind(df_sgrna, df_gene), snakemake@output[["pr"]], row.names = FALSE)

### Depletion profile of gene signatures
# Load gene sets
data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)

signatures <- list(
  Ribosomal_Proteins = EssGenes.ribosomalProteins,
  DNA_Replication = EssGenes.DNA_REPLICATION_cons,
  RNA_polymerase = EssGenes.KEGG_rna_polymerase,
  Proteasome = EssGenes.PROTEASOME_cons,
  Spliceosome = EssGenes.SPLICEOSOME_cons,
  Common_essential = bagel_essential,
  Non_essential = bagel_non_essential
)

# The rank at which gene-level FDR crosses 5% — ccr.VisDepAndSig computes
# this internally (via the same unexported helper) purely to draw its
# threshold line; call it directly rather than re-deriving it in Python; the
# gene/signature data is otherwise plain enough to re-plot client-side.
fdr5perc <- CRISPRcleanR:::ccr.fixedFDRthreshold(
  FCsprofile = gene_fcs,
  TruePositives = bagel_essential,
  TrueNegatives = bagel_non_essential,
  th = 0.05
)

print("Saving depletion rank gene data...")
write.csv(
  data.frame(
    gene = names(gene_fcs),
    fc = as.numeric(gene_fcs),
    fdr5perc_rank = fdr5perc$RANK
  ),
  snakemake@output[["drnk_genes"]],
  row.names = FALSE
)

print("Saving depletion rank signature membership...")
write.csv(
  bind_rows(lapply(names(signatures), function(name) {
    data.frame(signature = name, gene = signatures[[name]])
  })),
  snakemake@output[["drnk_signatures"]],
  row.names = FALSE
)

# Close redirection of output/messages
sink(type = "message")
sink(type = "output")
close(log)
