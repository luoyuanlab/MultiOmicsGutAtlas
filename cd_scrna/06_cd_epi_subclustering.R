# CD Atlas: Epithelial Cell Subclustering
#
# Iterative subclustering of epithelial cells from the CD atlas.
# Two rounds of HVG selection → ScaleData → PCA → Harmony → UMAP →
# clustering, with removal of contaminating clusters between iterations.
# A second pass (epi2) removes low-quality / doublet clusters for
# the final clean epithelial annotation.
#
# Inputs:
#   - 9p_harmony_integ.RDS            : Harmony-integrated atlas
# Outputs:
#   - epithelial/epi_iter1.RDS        : iter1 annotated object
#   - epithelial/epi_iter2.RDS        : final clean object (iter1_anno column)
#   - epithelial/epi_iter2_umap.csv   : UMAP + annotation for Python

library(Seurat)
library(harmony)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR <- "/path/to/cd/scrna/output"
EPI_DIR  <- file.path(DATA_DIR, "epithelial")
dir.create(EPI_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Helper: remove CC, RP, MT, HSP genes from HVG list ────────────────────────
filter_hvg <- function(seurat_obj) {
  hvgs <- VariableFeatures(seurat_obj)
  hvgs[!grepl("^RP[SL]|^MT-|^HSP|^DNAJ|^MKI67|^TOP2A|^UBB|^UBC", hvgs)]
}

# ── Load atlas and subset epithelial cells ─────────────────────────────────────
scrna <- readRDS(file.path(DATA_DIR, "9p_harmony_integ.RDS"))
epi   <- scrna[, scrna$RNA_snn_res.0.1 %in% c("1", "4")]  # epithelial cluster(s)

# Add binary condition and tissue labels
epi$condition_binary <- case_when(
  epi$condition %in% c("Inactive (Trt Naive)", "Inactive") ~ "Inactive",
  epi$condition == "Control"                               ~ "Healthy",
  epi$condition == "CD"                                    ~ "CD",
  TRUE                                                     ~ "Active"
)
epi$tissue_region_binary <- case_when(
  epi$tissue_region %in% c("Duodenum", "Ileum", "Ileum (Terminal)",
                            "Jejunum", "Small Bowel")       ~ "Small Intestine",
  TRUE                                                      ~ "Colon"
)

# ── Iteration 0 ───────────────────────────────────────────────────────────────
epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi, nfeatures = 2000)
VariableFeatures(epi) <- filter_hvg(epi)
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 25)
epi <- RunHarmony(epi,
                  group.by.vars = c("source", "patient"),
                  dims.use = 1:25,
                  theta = c(1, 1), lambda = c(2, 2))
epi <- RunUMAP(epi,   reduction = "harmony", dims = 1:25)
epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:25)
epi <- FindClusters(epi, resolution = 0.3)
saveRDS(epi, file.path(EPI_DIR, "epi_iter0.RDS"))

remove_epi0 <- c("12", "13", "14")
epi1_input  <- epi[, !epi$seurat_clusters %in% remove_epi0]

# ── Iteration 1: source-corrected clustering ───────────────────────────────────
# P19 source correction: re-split layers by source before re-normalizing
epi1_input[["RNA"]] <- JoinLayers(epi1_input[["RNA"]])
epi1_input[["RNA"]] <- split(epi1_input[["RNA"]], f = epi1_input$source)

epi1_input <- NormalizeData(epi1_input)
epi1_input <- FindVariableFeatures(epi1_input, nfeatures = 2000)
VariableFeatures(epi1_input) <- filter_hvg(epi1_input)
epi1_input <- ScaleData(epi1_input)
epi1_input <- RunPCA(epi1_input, npcs = 22)
epi1_input <- RunHarmony(epi1_input,
                         group.by.vars = c("source", "patient"),
                         dims.use = 1:22,
                         theta = c(1, 1), lambda = c(2, 2))
epi1_input <- RunUMAP(epi1_input, reduction = "harmony", dims = 1:22)
epi1_input <- FindNeighbors(epi1_input, reduction = "harmony", dims = 1:22)
epi1_input <- FindClusters(epi1_input, resolution = 0.3)

# ── Manual annotation: iteration 1 ───────────────────────────────────────────
# Based on canonical epithelial markers: EPCAM, OLFM4 (Stem), MUC2 (Goblet),
# CHGA (EEC), BEST4/OTOP2 (BEST4+), DCLK1 (Tuft), LYZ (Paneth),
# SI/FABP1 (Enterocyte), CA1/CA4/ITGA6 (Progenitors)
epi1_input$iter1_anno <- case_when(
  epi1_input$seurat_clusters == "0"  ~ "Enterocyte",
  epi1_input$seurat_clusters == "1"  ~ "Stem",
  epi1_input$seurat_clusters == "2"  ~ "Goblet",
  epi1_input$seurat_clusters == "3"  ~ "Immature Enterocyte",
  epi1_input$seurat_clusters == "4"  ~ "Absorptive Progenitor (SI)",
  epi1_input$seurat_clusters == "5"  ~ "Absorptive Progenitor (Col)",
  epi1_input$seurat_clusters == "6"  ~ "BEST4+ Cells",
  epi1_input$seurat_clusters == "7"  ~ "Tuft",
  epi1_input$seurat_clusters == "8"  ~ "Enteroendocrine",
  epi1_input$seurat_clusters == "9"  ~ "Paneth",
  epi1_input$seurat_clusters == "10" ~ "Microfold",
  TRUE                               ~ "Unknown"
)

saveRDS(epi1_input, file.path(EPI_DIR, "epi_iter1.RDS"))

# ── Iteration 2: remove low-quality / doublet clusters ────────────────────────
remove_epi1 <- c("10", "12")
epi2 <- epi1_input[, !epi1_input$seurat_clusters %in% remove_epi1]
saveRDS(epi2, file.path(EPI_DIR, "epi_iter2.RDS"))

# ── Export for Python notebooks ───────────────────────────────────────────────
# NOTE: atlas assembly uses iter1_anno column from epi_iter2 object
umap_epi            <- as.data.frame(epi2@reductions$umap@cell.embeddings)
colnames(umap_epi)  <- c("umapharmonyepi1iter1_1", "umapharmonyepi1iter1_2")
umap_epi$cell       <- rownames(umap_epi)
umap_epi$iter1_anno <- epi2$iter1_anno
write.csv(umap_epi, file.path(EPI_DIR, "epi_iter2_umap.csv"))
