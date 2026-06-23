library(Seurat)
library(dplyr)
library(ggplot2)

# ── Set cohort here ────────────────────────────────────────────────────────────
# cohort <- "biopsy"    # CD WTX biopsies
cohort <- "resection"   # CD surgical resection

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRNA_DIR <- "/path/to/scrna/cd"

if (cohort == "biopsy") {
  COSMX_DIR  <- "/path/to/cosmx_data/cd_6k_wtx/whole_trans/Processed_merged"
  cosmx_csv  <- file.path(COSMX_DIR, "count_csv/all_treg_raw_counts.csv")
  out_rds    <- file.path(COSMX_DIR, "rds/treg_cosmx_with_labels.rds")
  out_csv    <- file.path(COSMX_DIR, "anno/alltreg_pt_transfer_labels.csv")
  plot_dir   <- file.path(COSMX_DIR, "plots")
  tag        <- "alltregs"
} else if (cohort == "resection") {
  COSMX_DIR  <- "/path/to/cosmx_data/cd_resection/combined"
  cosmx_csv  <- file.path(COSMX_DIR, "raw_counts_strict_tregs.csv")
  out_rds    <- file.path(COSMX_DIR, "strict_treg_ptbin_anno.RDS")
  out_csv    <- file.path(COSMX_DIR, "strict_treg_pt_transfer_labels.csv")
  plot_dir   <- file.path(COSMX_DIR, "plots")
  tag        <- "strict_tregs"
} else {
  stop("cohort must be 'biopsy' or 'resection'")
}

# ── PART 1: Load and prepare scRNA-seq reference ───────────────────────────────
scrna_counts   <- read.csv(file.path(SCRNA_DIR, "treg_raw_counts.csv"), row.names = 1)
scrna_counts   <- t(scrna_counts)
scrna_metadata <- read.csv(file.path(SCRNA_DIR, "treg_pt_bin.csv"), row.names = 1)

print(paste("scRNA genes:", nrow(scrna_counts)))
print(paste("scRNA cells:", ncol(scrna_counts)))
print(paste("Metadata rows:", nrow(scrna_metadata)))

scrna <- CreateSeuratObject(
  counts    = as.matrix(scrna_counts),
  meta.data = scrna_metadata,
  assay     = "RNA"
)

print("Pseudotime bin distribution:")
print(table(scrna$pt_bin))

scrna <- NormalizeData(scrna)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000)
scrna <- ScaleData(scrna)
scrna <- RunPCA(scrna, npcs = 30, verbose = FALSE)
scrna <- RunUMAP(scrna, reduction = "pca", dims = 1:30)
saveRDS(scrna, file.path(SCRNA_DIR, "treg_pt_bin.RDS"))

# ── PART 2: Load and prepare CosMx query ──────────────────────────────────────
cosmx_counts <- read.csv(cosmx_csv, row.names = 1)
cosmx_counts <- t(cosmx_counts)

print(paste("CosMx genes:", nrow(cosmx_counts)))
print(paste("CosMx cells:", ncol(cosmx_counts)))

cosmx <- CreateSeuratObject(counts = as.matrix(cosmx_counts), assay = "RNA")
cosmx <- NormalizeData(cosmx)
cosmx <- FindVariableFeatures(cosmx, selection.method = "vst", nfeatures = 2000)
cosmx <- ScaleData(cosmx)
cosmx <- RunPCA(cosmx, npcs = 30, verbose = FALSE)
cosmx <- RunUMAP(cosmx, reduction = "pca", dims = 1:30)

# ── PART 3: Gene overlap ───────────────────────────────────────────────────────
common_genes <- intersect(rownames(scrna), rownames(cosmx))
print(paste("Common genes:", length(common_genes)))
print(paste("scRNA-only genes:", length(setdiff(rownames(scrna), rownames(cosmx)))))
print(paste("CosMx-only genes:", length(setdiff(rownames(cosmx), rownames(scrna)))))

# ── PART 4: Label transfer ─────────────────────────────────────────────────────
transfer_anchors <- FindTransferAnchors(
  reference           = scrna,
  query               = cosmx,
  dims                = 1:30,
  reference.reduction = "pca",
  normalization.method = "LogNormalize"
)

print(paste("Number of anchors found:", nrow(transfer_anchors@anchors)))

scrna$pt_bin <- as.character(scrna$pt_bin)
predictions <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = scrna$pt_bin,
  dims             = 1:30,
  weight.reduction = cosmx[["pca"]]
)

cosmx <- AddMetaData(cosmx, metadata = predictions)
saveRDS(cosmx, out_rds)

print("Prediction score summary:")
print(summary(cosmx$prediction.score.max))

# ── PART 5: QC plots ───────────────────────────────────────────────────────────
ggplot(cosmx@meta.data, aes(x = prediction.score.max)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "Label Transfer Prediction Scores",
       x = "Max Prediction Score", y = "Number of Cells") +
  theme_classic()
ggsave(file.path(plot_dir, paste0(tag, "_prediction_scores_histogram.pdf")), width = 8, height = 6)

DimPlot(cosmx, group.by = "predicted.id", label = TRUE) +
  ggtitle("CosMx: Transferred Pseudotime Bins")
ggsave(file.path(plot_dir, paste0(tag, "_cosmx_transferred_labels.pdf")), width = 8, height = 6)

ggplot(cosmx@meta.data, aes(x = predicted.id, y = prediction.score.max)) +
  geom_boxplot(fill = "steelblue") +
  labs(title = "Prediction Confidence by Pseudotime Bin",
       x = "Predicted Bin", y = "Prediction Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, paste0(tag, "_prediction_scores_by_bin.pdf")), width = 8, height = 6)

print(table(cosmx$predicted.id))

# Marker validation
early_markers <- c("BATF", "TNFRSF4", "TNFRSF18", "CD27", "ARID5B", "CTLA4")
late_markers  <- c("KLRD1", "NKG7", "CCL5", "CCL4", "CCL3", "CD69")

for (markers in list(early_markers, late_markers)) {
  tag_m <- if (identical(markers, early_markers)) "early" else "late"
  avail <- markers[markers %in% rownames(cosmx)]
  if (length(avail) > 0) {
    VlnPlot(cosmx, features = avail, group.by = "predicted.id", ncol = 3)
    ggsave(file.path(plot_dir, paste0(tag, "_", tag_m, "_marker_expression_by_bin.pdf")),
           width = 12, height = 8)
    FeaturePlot(cosmx, features = avail, reduction = "umap", ncol = 3)
    ggsave(file.path(plot_dir, paste0(tag, "_", tag_m, "_marker_expression_umap.pdf")),
           width = 12, height = 8)
  }
}

# ── PART 6: Export ─────────────────────────────────────────────────────────────
write.csv(cosmx@meta.data, out_csv)

summary_stats <- cosmx@meta.data %>%
  group_by(predicted.id) %>%
  summarise(
    n_cells       = n(),
    mean_score    = mean(prediction.score.max),
    median_score  = median(prediction.score.max),
    pct_high_conf = sum(prediction.score.max > 0.5) / n() * 100
  ) %>%
  arrange(predicted.id)

print(summary_stats)

message("Label transfer complete!")
message("Total CosMx cells: ", ncol(cosmx))
message("High confidence (>0.5): ", sum(cosmx$prediction.score.max > 0.5))
message("Low confidence (<=0.5): ", sum(cosmx$prediction.score.max <= 0.5))
