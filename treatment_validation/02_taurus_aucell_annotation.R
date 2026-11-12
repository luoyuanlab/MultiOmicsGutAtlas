library(anndata)
library(Seurat)
library(AUCell)
library(GSEABase)
library(data.table)
library(ggplot2)
library(dplyr)
set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRNA_DIR  <- "/path/to/scrna/cd"
TAURUS_DIR <- "/path/to/taurus"
ANNO_DIR   <- file.path(TAURUS_DIR, "anno")
PLOT_DIR   <- file.path(TAURUS_DIR, "plots")

dir.create(ANNO_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Shared gene set collection (50 cell types)
load(file.path(SCRNA_DIR, "cell_type_markers_gene_sets.RData"))

# ── PART 1: AUCell cell type annotation ───────────────────────────────────────
# Input: normalized CD counts exported from 01_taurus_preprocessing.ipynb
cd <- read_h5ad(file.path(TAURUS_DIR, "xnorm_counts_for_aucell.h5ad.gz"))
counts_mat <- t(cd$X)

gut_sub      <- subsetGeneSets(gsc, rownames(counts_mat))
all_markers  <- unique(unlist(lapply(gut_sub, geneIds)))
counts_subset <- counts_mat[rownames(counts_mat) %in% all_markers, ]

cells_AUC <- AUCell_run(counts_subset, gut_sub)
save(cells_AUC, file = file.path(ANNO_DIR, "norm_aucell.Rdata"))

cells_AUC_mat   <- t(getAUC(cells_AUC))
max_col_indices <- max.col(cells_AUC_mat, "first")
cell_labels     <- colnames(cells_AUC_mat)[max_col_indices]

result <- data.frame(cells = rownames(cells_AUC_mat), label = cell_labels)
write.csv(result, file.path(ANNO_DIR, "norm_aucell.csv"))
message("Done: AUCell annotation")
