# CD Atlas: T Cell and NK Cell Subclustering
#
# Separates CD3E+/CD3D+ T cells from NK cells, then iteratively
# subcluster T cells through three rounds of refinement.
# NK cells are saved after a single-iteration subclustering in
# 05_cd_other_lineages_subclustering.R.
#
# Inputs:
#   - 9p_harmony_integ.RDS               : Harmony-integrated atlas
# Outputs:
#   - t_cells/tcell_iter3.RDS            : final T cell object
#   - t_cells/nk_iter0.RDS              : NK cells (for further subclustering)
#   - t_cells/tcell_iter3_umap.csv       : UMAP + annotation for Python
#   - t_cells/t_harmony_iter3.csv        : Harmony embedding for pseudotime

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR <- "/path/to/cd/scrna/output"
T_DIR    <- file.path(DATA_DIR, "t_cells")
dir.create(T_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Helper: remove CC, RP, MT, HSP genes from HVG list ────────────────────────
filter_hvg <- function(seurat_obj) {
  hvgs <- VariableFeatures(seurat_obj)
  hvgs[!grepl("^RP[SL]|^MT-|^HSP|^DNAJ|^MKI67|^TOP2A|^UBB|^UBC", hvgs)]
}

# ── Load atlas and subset T/NK cells ──────────────────────────────────────────
scrna <- readRDS(file.path(DATA_DIR, "9p_harmony_integ.RDS"))
tcell <- scrna[, scrna$RNA_snn_res.0.1 %in% c("0")]  # T/NK cluster(s)

# Add binary condition and tissue labels
tcell$condition_binary <- case_when(
  tcell$condition %in% c("Inactive (Trt Naive)", "Inactive") ~ "Inactive",
  tcell$condition == "Control"                               ~ "Healthy",
  tcell$condition == "CD"                                    ~ "CD",
  TRUE                                                       ~ "Active"
)
tcell$tissue_region_binary <- case_when(
  tcell$tissue_region %in% c("Duodenum", "Ileum", "Ileum (Terminal)",
                              "Jejunum", "Small Bowel")       ~ "Small Intestine",
  TRUE                                                        ~ "Colon"
)

# ── Iteration 0: separate T cells from NK cells ────────────────────────────────
tcell <- NormalizeData(tcell)
tcell <- FindVariableFeatures(tcell, nfeatures = 2000)
VariableFeatures(tcell) <- filter_hvg(tcell)
tcell <- ScaleData(tcell)
tcell <- RunPCA(tcell, npcs = 25)
tcell <- RunHarmony(tcell,
                    group.by.vars = c("source", "patient"),
                    dims.use = 1:25,
                    theta = c(1, 1), lambda = c(2, 2))
tcell <- RunUMAP(tcell,   reduction = "harmony", dims = 1:25)
tcell <- FindNeighbors(tcell, reduction = "harmony", dims = 1:25)
tcell <- FindClusters(tcell, resolution = 0.3)

# Identify NK cluster by NCAM1/KLRD1/FCGR3A expression, no CD3E/CD3D
# Save NK cells separately for subclustering in 05_cd_other_lineages_subclustering.R
nk_clusters <- c("11")   # NK-defining cluster(s)
nk   <- tcell[, tcell$seurat_clusters %in% nk_clusters]
saveRDS(nk, file.path(T_DIR, "nk_iter0.RDS"))

# Retain only CD3E+/CD3D+ T cells
tcell1 <- tcell[, !tcell$seurat_clusters %in% nk_clusters]

# ── Iteration 1 ───────────────────────────────────────────────────────────────
tcell1 <- NormalizeData(tcell1)
tcell1 <- FindVariableFeatures(tcell1, nfeatures = 2000)
VariableFeatures(tcell1) <- filter_hvg(tcell1)
tcell1 <- ScaleData(tcell1)
tcell1 <- RunPCA(tcell1, npcs = 25)
tcell1 <- RunHarmony(tcell1,
                     group.by.vars = c("source", "patient"),
                     dims.use = 1:25,
                     theta = c(1, 1), lambda = c(2, 2))
tcell1 <- RunUMAP(tcell1, reduction = "harmony", dims = 1:25)
tcell1 <- FindNeighbors(tcell1, reduction = "harmony", dims = 1:25)
tcell1 <- FindClusters(tcell1, resolution = 0.3)
saveRDS(tcell1, file.path(T_DIR, "tcell_iter1.RDS"))

remove_iter1 <- c("13", "14")
tcell2 <- tcell1[, !tcell1$seurat_clusters %in% remove_iter1]

# ── Iteration 2 ───────────────────────────────────────────────────────────────
tcell2 <- NormalizeData(tcell2)
tcell2 <- FindVariableFeatures(tcell2, nfeatures = 2000)
VariableFeatures(tcell2) <- filter_hvg(tcell2)
tcell2 <- ScaleData(tcell2)
tcell2 <- RunPCA(tcell2, npcs = 25)
tcell2 <- RunHarmony(tcell2,
                     group.by.vars = c("source", "patient"),
                     dims.use = 1:25,
                     theta = c(1, 1), lambda = c(2, 2))
tcell2 <- RunUMAP(tcell2, reduction = "harmony", dims = 1:25)
tcell2 <- FindNeighbors(tcell2, reduction = "harmony", dims = 1:25)
tcell2 <- FindClusters(tcell2, resolution = 0.3)
saveRDS(tcell2, file.path(T_DIR, "tcell_iter2.RDS"))

remove_iter2 <- c("11")
tcell3_input <- tcell2[, !tcell2$seurat_clusters %in% remove_iter2]

# ── Iteration 3: source-corrected final clustering ────────────────────────────
# P19 contributes non-overlapping gene sets in some layers;
# correct by using only the source-level Harmony grouping
tcell3_input[["RNA"]] <- JoinLayers(tcell3_input[["RNA"]])
tcell3_input[["RNA"]] <- split(tcell3_input[["RNA"]], f = tcell3_input$source)

tcell3_input <- NormalizeData(tcell3_input)
tcell3_input <- FindVariableFeatures(tcell3_input, nfeatures = 2000)
VariableFeatures(tcell3_input) <- filter_hvg(tcell3_input)
tcell3_input <- ScaleData(tcell3_input)
tcell3_input <- RunPCA(tcell3_input, npcs = 22)
tcell3_input <- RunHarmony(tcell3_input,
                           group.by.vars = c("source", "patient"),
                           dims.use = 1:22,
                           theta = c(1, 1), lambda = c(2, 2))
tcell3_input <- RunUMAP(tcell3_input, reduction = "harmony", dims = 1:22)
tcell3_input <- FindNeighbors(tcell3_input, reduction = "harmony", dims = 1:22)
tcell3_input <- FindClusters(tcell3_input, resolution = 0.3)

# ── Manual annotation: iteration 3 ───────────────────────────────────────────
# Based on marker genes: GZMK, GZMB, FOXP3, CXCR5, TRGV, CCR7, etc.
tcell3_input$iter3_anno <- case_when(
  tcell3_input$seurat_clusters == "0"  ~ "CD8+GZMK- Trm(Prolif)",
  tcell3_input$seurat_clusters == "1"  ~ "CD4+ Trm",
  tcell3_input$seurat_clusters == "2"  ~ "Resting T",
  tcell3_input$seurat_clusters == "3"  ~ "CD4+ Effector",
  tcell3_input$seurat_clusters == "4"  ~ "CD8+GZMK+ Trm",
  tcell3_input$seurat_clusters == "5"  ~ "Treg",
  tcell3_input$seurat_clusters == "6"  ~ "Gd T",
  tcell3_input$seurat_clusters == "7"  ~ "Tfh",
  tcell3_input$seurat_clusters == "8"  ~ "CD8+GZMK- Trm",
  tcell3_input$seurat_clusters == "9"  ~ "CD4+ Trm(Prolif)",
  tcell3_input$seurat_clusters == "10" ~ "Exhausted CD4+ Tfh",
  tcell3_input$seurat_clusters == "11" ~ "Gd T(Vd29g+)",
  TRUE                                 ~ "Unknown"
)

saveRDS(tcell3_input, file.path(T_DIR, "tcell_iter3.RDS"))

# ── Export for Python notebooks ───────────────────────────────────────────────
umap_t             <- as.data.frame(tcell3_input@reductions$umap@cell.embeddings)
colnames(umap_t)   <- c("umapharmonytcelliter3_1", "umapharmonytcelliter3_2")
umap_t$cell        <- rownames(umap_t)
umap_t$iter3_anno  <- tcell3_input$iter3_anno
write.csv(umap_t, file.path(T_DIR, "tcell_iter3_umap.csv"))

harm_emb <- as.data.frame(tcell3_input@reductions$harmony@cell.embeddings)
write.csv(harm_emb, file.path(T_DIR, "t_harmony_iter3.csv"))
