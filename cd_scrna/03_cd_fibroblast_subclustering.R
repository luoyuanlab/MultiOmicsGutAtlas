# CD Atlas: Fibroblast / Connective Tissue Subclustering
#
# Iterative subclustering of connective tissue cells from the CD atlas.
# Two rounds of HVG selection → ScaleData → PCA → Harmony → UMAP →
# clustering, with manual removal of contaminating clusters at each
# iteration.
#
# Inputs:
#   - 9p_harmony_integ.RDS              : Harmony-integrated atlas
# Outputs:
#   - connective/con_iter2.RDS          : final subclustered object
#   - connective/con_iter2_umap.csv     : UMAP + annotation for Python
#   - connective/con_harmony_iter2.csv  : Harmony embedding for pseudotime

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/cd/scrna/output"
CON_DIR    <- file.path(DATA_DIR, "connective")
dir.create(CON_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Helper: remove CC, RP, MT, HSP genes from HVG list ────────────────────────
filter_hvg <- function(seurat_obj) {
  hvgs <- VariableFeatures(seurat_obj)
  hvgs[!grepl("^RP[SL]|^MT-|^HSP|^DNAJ|^MKI67|^TOP2A|^UBB|^UBC", hvgs)]
}

# ── Load atlas and subset connective tissue cells ──────────────────────────────
scrna <- readRDS(file.path(DATA_DIR, "9p_harmony_integ.RDS"))
con   <- scrna[, scrna$RNA_snn_res.0.1 %in% c("6")]  # connective cluster(s)

# Add binary condition and tissue labels
con$condition_binary <- case_when(
  con$condition %in% c("Inactive (Trt Naive)", "Inactive") ~ "Inactive",
  con$condition == "Control"                               ~ "Healthy",
  con$condition == "CD"                                    ~ "CD",
  TRUE                                                     ~ "Active"
)
con$tissue_region_binary <- case_when(
  con$tissue_region %in% c("Duodenum", "Ileum", "Ileum (Terminal)",
                            "Jejunum", "Small Bowel")      ~ "Small Intestine",
  TRUE                                                     ~ "Colon"
)

# ── Iteration 0 ───────────────────────────────────────────────────────────────
con <- NormalizeData(con)
con <- FindVariableFeatures(con, nfeatures = 2000)
VariableFeatures(con) <- filter_hvg(con)
con <- ScaleData(con)
con <- RunPCA(con, npcs = 20)
con <- RunHarmony(con,
                  group.by.vars = c("orig.ident", "patient"),
                  dims.use = 1:20,
                  theta = c(1, 1), lambda = c(2, 2))
con <- RunUMAP(con,   reduction = "harmony", dims = 1:20)
con <- FindNeighbors(con, reduction = "harmony", dims = 1:20)
con <- FindClusters(con, resolution = 0.3)
saveRDS(con, file.path(CON_DIR, "con_iter0.RDS"))

# Remove contaminating clusters
remove_iter0 <- c("9", "10", "11")
con1 <- con[, !con$seurat_clusters %in% remove_iter0]

# ── Iteration 1 ───────────────────────────────────────────────────────────────
con1 <- NormalizeData(con1)
con1 <- FindVariableFeatures(con1, nfeatures = 2000)
VariableFeatures(con1) <- filter_hvg(con1)
con1 <- ScaleData(con1)
con1 <- RunPCA(con1, npcs = 20)
con1 <- RunHarmony(con1,
                   group.by.vars = c("orig.ident", "patient"),
                   dims.use = 1:20,
                   theta = c(1, 1), lambda = c(2, 2))
con1 <- RunUMAP(con1, reduction = "harmony", dims = 1:20)
con1 <- FindNeighbors(con1, reduction = "harmony", dims = 1:20)
con1 <- FindClusters(con1, resolution = 0.3)
saveRDS(con1, file.path(CON_DIR, "con_iter1.RDS"))

remove_iter1 <- c("9", "10")
con2 <- con1[, !con1$seurat_clusters %in% remove_iter1]

# ── Iteration 2: final clustering ─────────────────────────────────────────────
con2 <- NormalizeData(con2)
con2 <- FindVariableFeatures(con2, nfeatures = 2000)
VariableFeatures(con2) <- filter_hvg(con2)
con2 <- ScaleData(con2)
con2 <- RunPCA(con2, npcs = 20)
con2 <- RunHarmony(con2,
                   group.by.vars = c("orig.ident", "patient"),
                   dims.use = 1:20,
                   theta = c(1, 1), lambda = c(2, 2))
con2 <- RunUMAP(con2, reduction = "harmony", dims = 1:20)
con2 <- FindNeighbors(con2, reduction = "harmony", dims = 1:20)
con2 <- FindClusters(con2, resolution = 0.3)

# ── Manual annotation: iteration 2 ───────────────────────────────────────────
# Based on marker genes: ADAMDEC1, OGN, RSPO3, PLA2G2A, VSTM2A,
# HHIP, GREM2, RERGL, CD36, SELENOP, activated fib markers (ACTA2, etc.)
con2$iter2_anno <- case_when(
  con2$seurat_clusters == "0"  ~ "ADAMDEC1+ Fib",
  con2$seurat_clusters == "1"  ~ "VSTM2A+ Crypt Top Fib",
  con2$seurat_clusters == "2"  ~ "OGN+RSPO3+ Fib",
  con2$seurat_clusters == "3"  ~ "ADAMDEC1+ Activated Fib",
  con2$seurat_clusters == "4"  ~ "SELENOP+ Fib",
  con2$seurat_clusters == "5"  ~ "HHIP+ Myofibroblast",
  con2$seurat_clusters == "6"  ~ "GREM2+ Myofibroblast",
  con2$seurat_clusters == "7"  ~ "PLA2G2A+ ECM Fib",
  con2$seurat_clusters == "8"  ~ "RERGL+ Contractile Pericyte",
  con2$seurat_clusters == "9"  ~ "CD36+ Pericyte",
  con2$seurat_clusters == "10" ~ "Activated Fib",
  con2$seurat_clusters == "11" ~ "T cell-Interacting Fib",
  TRUE                         ~ "Unknown"
)

saveRDS(con2, file.path(CON_DIR, "con_iter2.RDS"))

# ── Export for Python notebooks ───────────────────────────────────────────────
umap_con           <- as.data.frame(con2@reductions$umap@cell.embeddings)
colnames(umap_con) <- c("umapharmonyconiter2_1", "umapharmonyconiter2_2")
umap_con$cell      <- rownames(umap_con)
umap_con$iter2_anno <- con2$iter2_anno
write.csv(umap_con, file.path(CON_DIR, "con_iter2_umap.csv"))

harm_emb <- as.data.frame(con2@reductions$harmony@cell.embeddings)
write.csv(harm_emb, file.path(CON_DIR, "con_harmony_iter2.csv"))
