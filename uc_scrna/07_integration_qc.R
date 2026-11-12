library(Seurat)
library(SingleCellExperiment)
library(bluster)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/uc/scrna/output"
OUTPUT_DIR <- "/path/to/uc/scrna/output/qc"

so <- readRDS(file.path(DATA_DIR, "uc_atlas_harmony_integrated.RDS"))

# ── 1. Export embeddings for scIB silhouette scoring (Supplementary Table 3) ──
# scIB lineage-stratified silhouette analysis is run externally in Python;
# embeddings exported here serve as input.
pca_df <- as.data.frame(Embeddings(so, reduction = "pca"))
pca_df$cell_id <- rownames(pca_df)
write.csv(pca_df, file.path(OUTPUT_DIR, "seurat_pca_embeddings.csv"), row.names = FALSE)

harmony_df <- as.data.frame(Embeddings(so, reduction = "harmony"))
harmony_df$cell_id <- rownames(harmony_df)
write.csv(harmony_df, file.path(OUTPUT_DIR, "seurat_harmony_embeddings.csv"), row.names = FALSE)

# ── 2. Approximate silhouette scores on Harmony embeddings ────────────────────
# Measures within-cluster cohesion vs. between-cluster separation after integration
sce <- as.SingleCellExperiment(so)
colLabels(sce) <- so$RNA_snn_res.0.1

sil <- approxSilhouette(reducedDim(sce, "HARMONY"), clusters = colLabels(sce))
sil_df <- as.data.frame(sil)
sil_df$cluster <- colLabels(sce)
sil_df$closest <- ifelse(sil_df$width > 0, colLabels(sce), sil_df$other)
cat("Overall mean silhouette width:", mean(sil_df$width), "\n")
write.csv(sil_df, file.path(OUTPUT_DIR, "silhouette_scores.csv"), row.names = FALSE)

# ── 3. Neighborhood purity ────────────────────────────────────────────────────
pure <- neighborPurity(reducedDim(sce, "HARMONY"), colLabels(sce))
pure_df <- as.data.frame(pure)
pure_df$cluster <- colLabels(sce)
cat("Overall mean purity:", mean(pure_df$purity), "\n")
write.csv(pure_df, file.path(OUTPUT_DIR, "neighborhood_purity.csv"), row.names = FALSE)
