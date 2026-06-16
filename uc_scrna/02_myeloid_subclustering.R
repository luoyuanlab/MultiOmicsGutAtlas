library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(RColorBrewer)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/uc/scrna/output"
OUTPUT_DIR <- "/path/to/uc/scrna/output/myeloid"

# ── Helper: fix patient IDs across datasets ────────────────────────────────────
fix_patient_ids <- function(seur) {
  seur$Patient <- seur$sample
  # P1: paired inflamed / non-inflamed — strip condition suffix
  seur$Patient[seur$orig.ident == "P1"] <- gsub("R|S|UN|INF|MAR", "",
    seur@meta.data[seur$orig.ident == "P1", "sample"])
  # P4 (Smillie): strip sample suffix
  seur$Patient[seur$orig.ident == "P4"] <- gsub("\\..*", "",
    gsub(".*#_", "", seur@meta.data[seur$orig.ident == "P4", "sample"]))
  # P6: strip tissue label
  seur$Patient[seur$orig.ident == "P6"] <- gsub("_M.*|_F.*", "",
    gsub(".*#_", "", seur@meta.data[seur$orig.ident == "P6", "sample"]))
  seur
}

# ── Helper: iterative Harmony sub-clustering ───────────────────────────────────
recluster <- function(seur, pca_name, harmony_name, umap_name,
                      npcs = 30, resolution = 0.5, cluster_name = "cluster") {
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur)
  # Remove cell cycle, ribosomal, mitochondrial, heat-shock genes from HVGs
  rg <- which(VariableFeatures(seur) %in% c(cc.genes$s.genes, cc.genes$g2m.genes))
  if (length(rg)) VariableFeatures(seur) <- VariableFeatures(seur)[-rg]
  rg <- grep("^MT-|^RPL|^RPS|^MT1|^MT2|^MRPL|^MRPS|^HSP", VariableFeatures(seur))
  if (length(rg)) VariableFeatures(seur) <- VariableFeatures(seur)[-rg]
  seur <- ScaleData(seur)
  seur <- RunPCA(seur, assay = "RNA", npcs = npcs, reduction.name = pca_name)
  seur <- RunHarmony(seur,
                     reduction.use  = pca_name,
                     reduction.save = harmony_name,
                     dims.use       = 1:npcs,
                     group.by.vars  = c("orig.ident", "Patient"),
                     theta = c(1, 1), lambda = c(2, 2))
  seur <- RunUMAP(seur, assay = "RNA", reduction = harmony_name,
                  reduction.name = umap_name, dims = 1:npcs,
                  min.dist = 0.5, n.neighbors = 50)
  seur <- FindNeighbors(seur, assay = "RNA", reduction = harmony_name, dims = 1:npcs)
  seur <- FindClusters(seur, resolution = resolution, cluster.name = cluster_name)
  seur
}

# ── 1. Load myeloid subset from primary atlas ──────────────────────────────────
`%ni%` <- Negate(`%in%`)

seur <- readRDS(file.path(DATA_DIR, "subcluster_myeloid_raw.RDS"))
seur[["SCT"]] <- NULL
seur <- fix_patient_ids(seur)

# ── 2. Iteration 1: initial sub-clustering ─────────────────────────────────────
seur <- recluster(seur,
                  pca_name     = "pca_mye_iter1",
                  harmony_name = "harmony_mye_iter1",
                  umap_name    = "umap_harmony_mye_iter1",
                  npcs         = 30, resolution = 0.5,
                  cluster_name = "mye_iter1_cluster")

seur[["mye_iter1_cluster_num"]] <- as.numeric(as.character(seur[["mye_iter1_cluster"]][, 1]))
markers_iter1 <- FindAllMarkers(seur, assay = "RNA", only.pos = TRUE,
                                min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_iter1 %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "mye_iter1_markers.csv"))
saveRDS(seur, file.path(OUTPUT_DIR, "mye_iter1.RDS"))

# ── 3. Iteration 2: remove doublet / contaminating clusters ───────────────────
# Clusters 7 and 15 identified as non-myeloid doublets by marker inspection
seur_iter2 <- subset(seur, subset = mye_iter1_cluster_num %ni% c(7, 15))
seur_iter2 <- recluster(seur_iter2,
                        pca_name     = "pca_mye_iter2",
                        harmony_name = "harmony_mye_iter2",
                        umap_name    = "umap_harmony_mye_iter2",
                        npcs         = 30, resolution = 0.5,
                        cluster_name = "mye_iter2_cluster")

seur_iter2[["mye_iter2_cluster_num"]] <- as.numeric(as.character(seur_iter2[["mye_iter2_cluster"]][, 1]))
markers_iter2 <- FindAllMarkers(seur_iter2, assay = "RNA", only.pos = TRUE,
                                min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_iter2 %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "mye_iter2_markers.csv"))
saveRDS(seur_iter2, file.path(OUTPUT_DIR, "mye_iter2.RDS"))

# ── 4. Iteration 3: final clean sub-clustering ────────────────────────────────
# Remove donor-specific cluster (patient N661) and transitional cluster 14
uc_mye_iter3 <- subset(seur_iter2,
                        subset = iter2_anno != "Transitory" & Patient != "N661")

uc_mye_iter3 <- recluster(uc_mye_iter3,
                           pca_name     = "pca_mye_iter3",
                           harmony_name = "harmony_mye_iter3",
                           umap_name    = "umap_harmony_mye_iter3",
                           npcs         = 30, resolution = 0.6,
                           cluster_name = "mye_iter3_cluster")

uc_mye_iter3[["mye_iter3_cluster_num"]] <- as.numeric(as.character(uc_mye_iter3[["mye_iter3_cluster"]][, 1]))
markers_iter3 <- FindAllMarkers(uc_mye_iter3, assay = "RNA", only.pos = TRUE,
                                min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_iter3 %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "mye_iter3_markers.csv"))

# ── 5. Final annotation (Fig. 1B, D; Ext. Fig. 2) ────────────────────────────
uc_mye_meta <- uc_mye_iter3@meta.data %>%
  mutate(iter3_anno = case_when(
    mye_iter3_cluster == 0  ~ "Neutrophil",
    mye_iter3_cluster == 1  ~ "SELENOP+SIGLEC1+ Mac",
    mye_iter3_cluster == 2  ~ "SELENOP+XIST- Mac",
    mye_iter3_cluster == 3  ~ "Mo-Mac",
    mye_iter3_cluster == 4  ~ "cDC2",
    mye_iter3_cluster == 5  ~ "Inflammatory Mo-Mac",
    mye_iter3_cluster == 6  ~ "Transitory cDC2/Mac",
    mye_iter3_cluster == 7  ~ "cDC1",
    mye_iter3_cluster == 8  ~ "SELENOP+MMP9+SPP1+ Mac",
    mye_iter3_cluster == 9  ~ "Cycling",
    mye_iter3_cluster == 10 ~ "CCR7+ Mature DC",
    mye_iter3_cluster == 11 ~ "SELENOP+MMP9+PLA2G2D+ Mac"
  ))
uc_mye_iter3[["iter3_anno"]] <- uc_mye_meta$iter3_anno
uc_mye_iter3 <- SetIdent(uc_mye_iter3, value = "iter3_anno")

saveRDS(uc_mye_iter3, file.path(OUTPUT_DIR, "mye_iter3.RDS"))

# ── 6. SELENOP+ vs SELENOP- DEG in CD68+ cells (input to GSEA; Fig. 1F) ──────
CD68_cells  <- subset(uc_mye_iter3, subset = CD68 > 0)
CD68_cells[["CellName"]] <- colnames(CD68_cells)
selenop_pos <- Cells(subset(CD68_cells, subset = SELENOP > 0))
CD68_cells@meta.data <- CD68_cells@meta.data %>%
  mutate(SELENOP_exp = ifelse(CellName %in% selenop_pos, "Pos", "Neg"))
CD68_cells <- SetIdent(CD68_cells, value = "SELENOP_exp")

deg_selenop <- FindMarkers(CD68_cells,
                           ident.1           = "Pos",
                           ident.2           = "Neg",
                           min.pct           = -Inf,
                           logfc.threshold   = -Inf,
                           min.cells.feature = 1,
                           min.cells.group   = 1)
write.csv(deg_selenop, file.path(OUTPUT_DIR, "cd68_selenop_pos_neg_degs.csv"))
