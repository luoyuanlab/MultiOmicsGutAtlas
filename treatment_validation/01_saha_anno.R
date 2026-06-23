library(Seurat)
library(AUCell)
library(GSEABase)
library(data.table)
library(ggplot2)
library(dplyr)
set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRNA_DIR  <- "/path/to/scrna/cd"
SAHA_DIR   <- "/path/to/saha"
CSV_DIR    <- file.path(SAHA_DIR, "csv")
ANNO_DIR   <- file.path(SAHA_DIR, "anno")
PLOT_DIR   <- file.path(SAHA_DIR, "plots")

dir.create(ANNO_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,  recursive = TRUE, showWarnings = FALSE)

# Shared gene set collection (50 cell types)
load(file.path(SCRNA_DIR, "cell_type_markers_gene_sets.RData"))

# ── PART 1: AUCell cell type annotation ───────────────────────────────────────
norm_6k <- fread(file.path(CSV_DIR, "6k_trt_norm_counts.csv"))
norm_6k <- as.data.frame(norm_6k)
rownames(norm_6k) <- norm_6k[[1]]
norm_6k <- norm_6k[, -1]

norm_so <- CreateSeuratObject(counts = t(norm_6k))
counts  <- GetAssayData(object = norm_so, layer = "counts")
gut_sub <- subsetGeneSets(gsc, rownames(counts))

cells_AUC <- AUCell_run(counts, gut_sub)
save(cells_AUC, file = file.path(ANNO_DIR, "6k_trt_norm_aucell.Rdata"))

cells_AUC_mat   <- t(getAUC(cells_AUC))
max_col_indices <- max.col(cells_AUC_mat, "first")
cell_labels     <- colnames(cells_AUC_mat)[max_col_indices]

result <- data.frame(cells = rownames(cells_AUC_mat), label = cell_labels)
write.csv(result, file.path(ANNO_DIR, "6k_trt_norm_aucell.csv"))
message("Done: AUCell annotation")

# ── PART 2: Treg pseudotime bin label transfer ─────────────────────────────────
# Reference: scRNA Treg pseudotime bins (built in cd_scrna/12_cd_treg_pseudotime.ipynb)
scrna <- readRDS(file.path(SCRNA_DIR, "treg_pt_bin.RDS"))

cosmx_counts <- read.csv(file.path(CSV_DIR, "6k_trt_treg_raw_counts.csv"), row.names = 1)
cosmx_counts <- t(cosmx_counts)

cosmx <- CreateSeuratObject(counts = as.matrix(cosmx_counts), assay = "RNA")
cosmx <- NormalizeData(cosmx)
cosmx <- FindVariableFeatures(cosmx, selection.method = "vst", nfeatures = 2000)
cosmx <- ScaleData(cosmx)
cosmx <- RunPCA(cosmx, npcs = 30, verbose = FALSE)
cosmx <- RunUMAP(cosmx, reduction = "pca", dims = 1:30)

transfer_anchors <- FindTransferAnchors(
  reference            = scrna,
  query                = cosmx,
  dims                 = 1:30,
  reference.reduction  = "pca",
  normalization.method = "LogNormalize"
)

scrna$pt_bin <- as.character(scrna$pt_bin)
predictions  <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = scrna$pt_bin,
  dims             = 1:30,
  weight.reduction = cosmx[["pca"]]
)

cosmx <- AddMetaData(cosmx, metadata = predictions)
saveRDS(cosmx, file.path(CSV_DIR, "6k_trt_treg_ptbin_anno.RDS"))

# QC plots
ggplot(cosmx@meta.data, aes(x = prediction.score.max)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "SAHA: Label Transfer Prediction Scores",
       x = "Max Prediction Score", y = "Number of Cells") +
  theme_classic()
ggsave(file.path(PLOT_DIR, "6k_trt_tregs_prediction_scores_histogram.pdf"), width = 8, height = 6)

DimPlot(cosmx, group.by = "predicted.id", label = TRUE) +
  ggtitle("SAHA CosMx: Transferred Pseudotime Bins")
ggsave(file.path(PLOT_DIR, "6k_trt_tregs_cosmx_transferred_labels.pdf"), width = 8, height = 6)

ggplot(cosmx@meta.data, aes(x = predicted.id, y = prediction.score.max)) +
  geom_boxplot(fill = "steelblue") +
  labs(title = "SAHA: Prediction Confidence by Pseudotime Bin",
       x = "Predicted Bin", y = "Prediction Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, "6k_trt_tregs_prediction_scores_by_bin.pdf"), width = 8, height = 6)

write.csv(cosmx@meta.data, file.path(CSV_DIR, "6k_trt_reg_pt_transfer_labels.csv"))
message("Done: Treg label transfer")
