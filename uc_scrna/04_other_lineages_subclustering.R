library(dplyr)
library(Seurat)
library(harmony)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/uc/scrna/output"
OUTPUT_DIR <- "/path/to/uc/scrna/output"

`%ni%` <- Negate(`%in%`)

# ── Helper: fix patient IDs across datasets ────────────────────────────────────
fix_patient_ids <- function(seur) {
  seur$Patient <- seur$sample
  seur$Patient[seur$orig.ident == "P1"] <- gsub("R|S|UN|INF|MAR", "",
    seur@meta.data[seur$orig.ident == "P1", "sample"])
  seur$Patient[seur$orig.ident == "P4"] <- gsub("\\..*", "",
    gsub(".*#_", "", seur@meta.data[seur$orig.ident == "P4", "sample"]))
  seur$Patient[seur$orig.ident == "P6"] <- gsub("_M.*|_F.*", "",
    gsub(".*#_", "", seur@meta.data[seur$orig.ident == "P6", "sample"]))
  seur
}

# ── Helper: standard iterative Harmony sub-clustering ─────────────────────────
recluster <- function(seur, pca_name, harmony_name, umap_name,
                      npcs = 25, resolution = 0.5, cluster_name = "cluster") {
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur)
  rg <- which(VariableFeatures(seur) %in% c(cc.genes$s.genes, cc.genes$g2m.genes))
  if (length(rg)) VariableFeatures(seur) <- VariableFeatures(seur)[-rg]
  rg <- grep("^MT-|^RPL|^RPS|^MT1|^MT2|^MRPL|^MRPS|^HSP", VariableFeatures(seur))
  if (length(rg)) VariableFeatures(seur) <- VariableFeatures(seur)[-rg]
  seur[["percent_rp"]] <- PercentageFeatureSet(seur, pattern = "^RPL|^RPS")
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

# ─────────────────────────────────────────────────────────────────────────────
# T CELLS (Ext. Fig. 1C)
# ─────────────────────────────────────────────────────────────────────────────
uc_t <- readRDS(file.path(DATA_DIR, "subcluster_t_cd3.RDS"))
uc_t <- fix_patient_ids(uc_t)

# Remove doublet clusters (9, 10, 11) identified by QC and marker inspection
t_iter1 <- subset(uc_t, subset = RNA_snn_res.0.2 %ni% c(9, 10, 11))
t_iter1[["RNA_snn_res.0.1"]] <- NULL; t_iter1[["RNA_snn_res.0.3"]] <- NULL
t_iter1[["SCT"]] <- NULL
t_iter1@reductions$umap <- NULL; t_iter1@reductions$pca <- NULL
t_iter1@reductions$harmony <- NULL

t_iter1 <- recluster(t_iter1,
                     pca_name     = "pca_cd3t_iter1",
                     harmony_name = "harmony_cd3t_iter1",
                     umap_name    = "umap_harmony_cd3t_iter1",
                     npcs = 25, resolution = 0.4,
                     cluster_name = "cd3t_iter1_cluster")

markers_t <- FindAllMarkers(t_iter1, assay = "RNA", only.pos = TRUE,
                            min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_t %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "t_iter1_markers.csv"))
saveRDS(t_iter1, file.path(OUTPUT_DIR, "t_iter1.RDS"))

# ─────────────────────────────────────────────────────────────────────────────
# B CELLS (Ext. Fig. 1C)
# ─────────────────────────────────────────────────────────────────────────────
uc_b  <- readRDS(file.path(DATA_DIR, "subcluster_b_ms4a1.RDS"))
uc_b2 <- readRDS(file.path(DATA_DIR, "subcluster_b_ms4a1_neg.RDS"))
uc_bm <- merge(uc_b, uc_b2)
for (col in c("RNA_snn_res.0.2", "RNA_snn_res.0.3", "RNA_snn_res.0.05",
              "RNA_snn_res.0.1", "SCT")) uc_bm[[col]] <- NULL

uc_bm <- fix_patient_ids(uc_bm)
# Remove PBMC-derived samples from P7 that lack mucosal representation
pbmc_p7 <- paste0("P7_#_", c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30))
uc_bm   <- subset(uc_bm, subset = sample %ni% pbmc_p7)

b_iter1 <- recluster(uc_bm,
                     pca_name     = "pca_b_iter1",
                     harmony_name = "harmony_b_iter1",
                     umap_name    = "umap_harmony_b_iter1",
                     npcs = 20, resolution = 0.4,
                     cluster_name = "b_iter1_cluster")

markers_b <- FindAllMarkers(b_iter1, assay = "RNA", only.pos = TRUE,
                            min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_b %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "b_iter1_markers.csv"))
saveRDS(b_iter1, file.path(OUTPUT_DIR, "b_iter1.RDS"))

# ─────────────────────────────────────────────────────────────────────────────
# PLASMA CELLS (Ext. Fig. 1C)
# ─────────────────────────────────────────────────────────────────────────────
uc_plasma <- readRDS(file.path(DATA_DIR, "subcluster_plasma_raw.RDS"))
uc_plasma <- fix_patient_ids(uc_plasma)

plasma_iter1 <- recluster(uc_plasma,
                          pca_name     = "pca_plasma_iter1",
                          harmony_name = "harmony_plasma_iter1",
                          umap_name    = "umap_harmony_plasma_iter1",
                          npcs = 20, resolution = 0.4,
                          cluster_name = "plasma_iter1_cluster")

markers_plasma <- FindAllMarkers(plasma_iter1, assay = "RNA", only.pos = TRUE,
                                 min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_plasma %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "plasma_iter1_markers.csv"))
saveRDS(plasma_iter1, file.path(OUTPUT_DIR, "plasma_iter1.RDS"))

# ─────────────────────────────────────────────────────────────────────────────
# NK CELLS (Ext. Fig. 1C)
# ─────────────────────────────────────────────────────────────────────────────
uc_nk <- readRDS(file.path(DATA_DIR, "subcluster_nk_raw.RDS"))
uc_nk <- fix_patient_ids(uc_nk)

nk_iter1 <- recluster(uc_nk,
                      pca_name     = "pca_nk_iter1",
                      harmony_name = "harmony_nk_iter1",
                      umap_name    = "umap_harmony_nk_iter1",
                      npcs = 20, resolution = 0.4,
                      cluster_name = "nk_iter1_cluster")

markers_nk <- FindAllMarkers(nk_iter1, assay = "RNA", only.pos = TRUE,
                             min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_nk %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "nk_iter1_markers.csv"))
saveRDS(nk_iter1, file.path(OUTPUT_DIR, "nk_iter1.RDS"))

# ─────────────────────────────────────────────────────────────────────────────
# EPITHELIAL CELLS (Ext. Fig. 1D)
# ─────────────────────────────────────────────────────────────────────────────
uc_epi <- readRDS(file.path(DATA_DIR, "subcluster_epi_raw.RDS"))
uc_epi <- fix_patient_ids(uc_epi)

epi_iter1 <- recluster(uc_epi,
                       pca_name     = "pca_epi_iter1",
                       harmony_name = "harmony_epi_iter1",
                       umap_name    = "umap_harmony_epi_iter1",
                       npcs = 25, resolution = 0.4,
                       cluster_name = "epi_iter1_cluster")

markers_epi <- FindAllMarkers(epi_iter1, assay = "RNA", only.pos = TRUE,
                              min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_epi %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "epi_iter1_markers.csv"))
saveRDS(epi_iter1, file.path(OUTPUT_DIR, "epi_iter1.RDS"))

# ─────────────────────────────────────────────────────────────────────────────
# ENDOTHELIAL CELLS (Ext. Fig. 1F)
# ─────────────────────────────────────────────────────────────────────────────
uc_endo <- readRDS(file.path(DATA_DIR, "subcluster_endo_raw.RDS"))
uc_endo <- fix_patient_ids(uc_endo)

endo_iter1 <- recluster(uc_endo,
                        pca_name     = "pca_endo_iter1",
                        harmony_name = "harmony_endo_iter1",
                        umap_name    = "umap_harmony_endo_iter1",
                        npcs = 20, resolution = 0.3,
                        cluster_name = "endo_iter1_cluster")

markers_endo <- FindAllMarkers(endo_iter1, assay = "RNA", only.pos = TRUE,
                               min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_endo %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "endo_iter1_markers.csv"))
saveRDS(endo_iter1, file.path(OUTPUT_DIR, "endo_iter1.RDS"))
