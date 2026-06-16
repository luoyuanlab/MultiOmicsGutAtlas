# CD Atlas: Myeloid Cell Subclustering (Fig. X)
#
# Iterative subclustering of CD68+ myeloid cells from the CD atlas.
# Four rounds of HVG selection → ScaleData → PCA → Harmony → UMAP →
# clustering, with manual removal of doublet/contaminating clusters
# at each iteration.
#
# Inputs:
#   - 9p_harmony_integ.RDS            : Harmony-integrated atlas
# Outputs (per iteration):
#   - myeloid/mye_iter<N>.RDS         : subclustered Seurat object
#   - myeloid/mye_iter4_umap.csv      : final UMAP + annotation for Python
#   - myeloid/mye_harmony_iter4.csv   : Harmony embedding for pseudotime

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/cd/scrna/output"
MYE_DIR    <- file.path(DATA_DIR, "myeloid")
dir.create(MYE_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Helper: remove CC, RP, MT, HSP genes from HVG list ────────────────────────
filter_hvg <- function(seurat_obj) {
  hvgs <- VariableFeatures(seurat_obj)
  hvgs[!grepl("^RP[SL]|^MT-|^HSP|^DNAJ|^MKI67|^TOP2A|^UBB|^UBC", hvgs)]
}

# ── Load full atlas and subset myeloid cells ───────────────────────────────────
scrna   <- readRDS(file.path(DATA_DIR, "9p_harmony_integ.RDS"))
myeloid <- scrna[, scrna$RNA_snn_res.0.1 %in% c("5")]  # myeloid cluster(s)

# Add binary condition and tissue labels
myeloid$condition_binary <- case_when(
  myeloid$condition %in% c("Inactive (Trt Naive)", "Inactive") ~ "Inactive",
  myeloid$condition == "Control"                               ~ "Healthy",
  myeloid$condition == "CD"                                    ~ "CD",
  TRUE                                                         ~ "Active"
)
myeloid$tissue_region_binary <- case_when(
  myeloid$tissue_region %in% c("Duodenum", "Ileum", "Ileum (Terminal)",
                                "Jejunum", "Small Bowel")     ~ "Small Intestine",
  TRUE                                                         ~ "Colon"
)

# ── Iteration 0 ───────────────────────────────────────────────────────────────
myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid, nfeatures = 2000)
VariableFeatures(myeloid) <- filter_hvg(myeloid)
myeloid <- ScaleData(myeloid)
myeloid <- RunPCA(myeloid, npcs = 25)
myeloid <- RunHarmony(myeloid,
                      group.by.vars = c("orig.ident", "patient"),
                      dims.use = 1:25,
                      theta = c(1, 1), lambda = c(2, 2))
myeloid <- RunUMAP(myeloid,   reduction = "harmony", dims = 1:25)
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:25)
myeloid <- FindClusters(myeloid, resolution = 0.3)
saveRDS(myeloid, file.path(MYE_DIR, "mye_iter0.RDS"))

# Remove doublet / contaminating clusters (non-myeloid)
remove_iter0 <- c("7", "9", "11", "13", "14", "15", "16")
mye1 <- myeloid[, !myeloid$seurat_clusters %in% remove_iter0]

# ── Iteration 1 ───────────────────────────────────────────────────────────────
mye1 <- NormalizeData(mye1)
mye1 <- FindVariableFeatures(mye1, nfeatures = 2000)
VariableFeatures(mye1) <- filter_hvg(mye1)
mye1 <- ScaleData(mye1)
mye1 <- RunPCA(mye1, npcs = 25)
mye1 <- RunHarmony(mye1,
                   group.by.vars = c("orig.ident", "patient"),
                   dims.use = 1:25,
                   theta = c(1, 1), lambda = c(2, 2))
mye1 <- RunUMAP(mye1, reduction = "harmony", dims = 1:25)
mye1 <- FindNeighbors(mye1, reduction = "harmony", dims = 1:25)
mye1 <- FindClusters(mye1, resolution = 0.3)
saveRDS(mye1, file.path(MYE_DIR, "mye_iter1.RDS"))

remove_iter1 <- c("12", "14")
mye2 <- mye1[, !mye1$seurat_clusters %in% remove_iter1]

# ── Iteration 2 ───────────────────────────────────────────────────────────────
mye2 <- NormalizeData(mye2)
mye2 <- FindVariableFeatures(mye2, nfeatures = 2000)
VariableFeatures(mye2) <- filter_hvg(mye2)
mye2 <- ScaleData(mye2)
mye2 <- RunPCA(mye2, npcs = 25)
mye2 <- RunHarmony(mye2,
                   group.by.vars = c("orig.ident", "patient"),
                   dims.use = 1:25,
                   theta = c(1, 1), lambda = c(2, 2))
mye2 <- RunUMAP(mye2, reduction = "harmony", dims = 1:25)
mye2 <- FindNeighbors(mye2, reduction = "harmony", dims = 1:25)
mye2 <- FindClusters(mye2, resolution = 0.3)
saveRDS(mye2, file.path(MYE_DIR, "mye_iter2.RDS"))

remove_iter2 <- c("12", "13")
mye3 <- mye2[, !mye2$seurat_clusters %in% remove_iter2]

# ── Iteration 3 ───────────────────────────────────────────────────────────────
mye3 <- NormalizeData(mye3)
mye3 <- FindVariableFeatures(mye3, nfeatures = 2000)
VariableFeatures(mye3) <- filter_hvg(mye3)
mye3 <- ScaleData(mye3)
mye3 <- RunPCA(mye3, npcs = 22)
mye3 <- RunHarmony(mye3,
                   group.by.vars = c("orig.ident", "patient"),
                   dims.use = 1:22,
                   theta = c(1, 1), lambda = c(2, 2))
mye3 <- RunUMAP(mye3, reduction = "harmony", dims = 1:22)
mye3 <- FindNeighbors(mye3, reduction = "harmony", dims = 1:22)
mye3 <- FindClusters(mye3, resolution = 0.3)

remove_iter3 <- c("12")
mye4 <- mye3[, !mye3$seurat_clusters %in% remove_iter3]

# ── Iteration 4: final clustering ─────────────────────────────────────────────
mye4 <- NormalizeData(mye4)
mye4 <- FindVariableFeatures(mye4, nfeatures = 2000)
VariableFeatures(mye4) <- filter_hvg(mye4)
mye4 <- ScaleData(mye4)
mye4 <- RunPCA(mye4, npcs = 22)
mye4 <- RunHarmony(mye4,
                   group.by.vars = c("orig.ident", "patient"),
                   dims.use = 1:22,
                   theta = c(1, 1), lambda = c(2, 2))
mye4 <- RunUMAP(mye4, reduction = "harmony", dims = 1:22)
mye4 <- FindNeighbors(mye4, reduction = "harmony", dims = 1:22)
mye4 <- FindClusters(mye4, resolution = 0.3)

# ── Manual annotation: iteration 4 ───────────────────────────────────────────
# Based on marker gene expression (SELENOP, SIGLEC1, XIST, MMP9, PLA2G2D,
# CD9, FCGR3A/B, CCR7, cDC markers, cycling)
mye4$iter4_anno <- case_when(
  mye4$seurat_clusters == "0"  ~ "SELENOP+SIGLEC1+ Mac",
  mye4$seurat_clusters == "1"  ~ "SELENOP+XIST- Mac",
  mye4$seurat_clusters == "2"  ~ "Mo-Mac",
  mye4$seurat_clusters == "3"  ~ "Inflammatory Mo-Mac",
  mye4$seurat_clusters == "4"  ~ "cDC2",
  mye4$seurat_clusters == "5"  ~ "Transitory cDC2/Mac",
  mye4$seurat_clusters == "6"  ~ "SELENOP+MMP9+PLA2G2D+ Mac",
  mye4$seurat_clusters == "7"  ~ "cDC1",
  mye4$seurat_clusters == "8"  ~ "CCR7+ Mature DC",
  mye4$seurat_clusters == "9"  ~ "Cycling",
  mye4$seurat_clusters == "10" ~ "Neutrophil",
  mye4$seurat_clusters == "11" ~ "SELENOP-CD9+ M2 Mac",
  TRUE                         ~ "Unknown"
)

saveRDS(mye4, file.path(MYE_DIR, "mye_iter4.RDS"))

# ── SELENOP+/- DEG analysis ───────────────────────────────────────────────────
mye4[["RNA"]] <- JoinLayers(mye4[["RNA"]])
mye4$selenop_group <- case_when(
  grepl("^SELENOP\\+", mye4$iter4_anno) ~ "SELENOP+",
  mye4$iter4_anno %in% c("Mo-Mac", "Inflammatory Mo-Mac",
                          "SELENOP-CD9+ M2 Mac")  ~ "SELENOP-",
  TRUE                                            ~ NA_character_
)
mye_selenop <- mye4[, !is.na(mye4$selenop_group)]
Idents(mye_selenop) <- "selenop_group"
selenop_markers <- FindMarkers(mye_selenop,
                               ident.1 = "SELENOP+",
                               ident.2 = "SELENOP-",
                               min.pct = 0.25,
                               logfc.threshold = 0.5)
write.csv(selenop_markers,
          file.path(MYE_DIR, "selenop_pos_vs_neg_markers.csv"))

# ── Export for Python notebooks ───────────────────────────────────────────────
umap_mye              <- as.data.frame(mye4@reductions$umap@cell.embeddings)
colnames(umap_mye)    <- c("umapharmonymyeiter4_1", "umapharmonymyeiter4_2")
umap_mye$cell         <- rownames(umap_mye)
umap_mye$iter4_anno   <- mye4$iter4_anno
write.csv(umap_mye, file.path(MYE_DIR, "mye_iter4_umap.csv"))

harm_emb <- as.data.frame(mye4@reductions$harmony@cell.embeddings)
write.csv(harm_emb, file.path(MYE_DIR, "mye_harmony_iter4.csv"))
