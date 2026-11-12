library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/uc/scrna/output"
OUTPUT_DIR <- "/path/to/uc/scrna/output/fibroblast"

`%ni%` <- Negate(`%in%`)

# ── Helper: derive patient IDs from sample names (shared across lineage scripts)
harmonise_patient_ids <- function(seur) {
  seur$Patient <- seur$sample
  seur$Patient[seur$orig.ident == "P1"] <- gsub("R|S|UN|INF|MAR", "",
    seur@meta.data[seur$orig.ident == "P1", "sample"])
  seur$Patient[seur$orig.ident == "P4"] <- gsub("\\..*", "",
    gsub(".*#_", "", seur@meta.data[seur$orig.ident == "P4", "sample"]))
  seur$Patient[seur$orig.ident == "P6"] <- gsub("_M.*|_F.*", "",
    gsub(".*#_", "", seur@meta.data[seur$orig.ident == "P6", "sample"]))
  seur
}

# ── Helper: iterative Harmony sub-clustering ───────────────────────────────────
recluster <- function(seur, pca_name, harmony_name, umap_name,
                      npcs = 20, resolution = 0.5, cluster_name = "cluster") {
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur)
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

# ── 1. Load fibroblast / connective tissue subset from primary atlas ───────────
seur <- readRDS(file.path(DATA_DIR, "subcluster_fibroblast_raw.RDS"))
seur[["SCT"]]          <- NULL
seur[["RNA_snn_res.0.05"]] <- NULL
seur[["RNA_snn_res.0.2"]]  <- NULL
seur[["RNA_snn_res.0.3"]]  <- NULL
seur[["pca"]] <- NULL; seur[["harmony"]] <- NULL; seur[["umap"]] <- NULL
seur <- harmonise_patient_ids(seur)

# ── 2. Iteration 1: initial sub-clustering ─────────────────────────────────────
# Doublet clusters 6, 7, 9 identified by marker inspection and QC metrics
seur <- subset(seur, subset = RNA_snn_res.0.1 %ni% c(6, 7, 9))
seur <- CellCycleScoring(seur, s.features = cc.genes$s.genes,
                         g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

seur <- recluster(seur,
                  pca_name     = "pca_con_iter1",
                  harmony_name = "harmony_con_iter1",
                  umap_name    = "umap_harmony_con_iter1",
                  npcs         = 20, resolution = 0.5,
                  cluster_name = "con_iter1_cluster")

seur[["con_iter1_cluster_num"]] <- as.numeric(as.character(seur[["con_iter1_cluster"]][, 1]))
markers_iter1 <- FindAllMarkers(seur, assay = "RNA", only.pos = TRUE,
                                min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers_iter1 %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC),
          file.path(OUTPUT_DIR, "con_iter1_markers.csv"))
saveRDS(seur, file.path(OUTPUT_DIR, "con_iter1.RDS"))

# ── 3. Subsequent iterations: repeat recluster → inspect → annotate ───────────
# (Repeated until the fine-grained fibroblast subtypes are stable; the annotated
# object is saved as con_iter4.RDS)
# Annotations include: ADAMDEC1+ Fib, OGN+RSPO3+ Fib, CXCL5+ Activated Fib,
# Activated Fib, GREM2+ Myofibroblast, HHIP+ Myofibroblast,
# RERGL+ Contractile Pericyte, CD36+ Pericyte, Cycling Conn
# (see Supp Fig 3A and Fig 2H for the fibroblast UMAP and pseudotime)

# ── 4. Supp Fig 3A: fibroblast marker feature plots on the iteration-4 object ──────
# Run after con_iter4.RDS exists; iteration-4 fibroblast UMAP (as in Fig 2C and 2H). Output: ed3a_fib_markers_iter4.png
O <- '/share/fsmresfiles/UC/revision_jenny/SuppFig3A_fib_marker_umaps'
uc_fib_iter4 <- readRDS('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/subcluster4/con_iter4.RDS')
ecm <- c('OGN','SFRP2','IGF1','GPNMB','RSPO3','GREM1','C3','THBS2')   # 8 marker genes, in Supp Fig 3A panel order
p <- FeaturePlot(uc_fib_iter4, features = ecm, raster = FALSE, reduction = 'umap_harmony_con_iter4', ncol = 4) &
     NoLegend() & labs(x = 'UMAP1', y = 'UMAP2')
ggsave(file.path(O, 'ed3a_fib_markers_iter4.png'), p, width = 12, height = 6, dpi = 300)
