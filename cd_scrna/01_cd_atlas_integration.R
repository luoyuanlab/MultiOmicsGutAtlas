# CD Atlas: Harmony Integration
#
# Loads the pre-merged, QC'd Seurat object produced by 00_cd_dataset_preparation.Rmd,
# runs Harmony integration at the source x patient level, and exports the integrated
# object for downstream lineage subclustering.
#
# Input:
#   - merged_cd/9p_merge_qc.RDS         : post-QC merged Seurat object
#
# Outputs:
#   - merged_cd/9p_harmony_integ.RDS    : Harmony-integrated atlas
#   - merged_cd/plots/                  : integration UMAP PDFs

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(harmony)

# -- Paths ---------------------------------------------------------------------
DATA_DIR   <- "/path/to/cd/scrna/data"
OUTPUT_DIR <- "/path/to/cd/scrna/output"

dir.create(file.path(OUTPUT_DIR, "merged_cd/plots"), showWarnings = FALSE, recursive = TRUE)

# -- Load pre-merged QC'd object -----------------------------------------------
so_qc <- readRDS(file.path(OUTPUT_DIR, "merged_cd/9p_merge_qc.RDS"))

# -- Gene QC: retain genes expressed in >=0.01% of cells ----------------------
counts       <- GetAssayData(so_qc, slot = "counts", assay = "RNA")
genes.pct    <- rowMeans(counts > 0) * 100
genes.filter <- names(genes.pct[genes.pct > 0.01])
so_qc        <- so_qc[genes.filter, ]

# -- QC violin plots -----------------------------------------------------------
so_qc <- SetIdent(so_qc, value = "source")
pdf(file.path(OUTPUT_DIR, "merged_cd/plots/post_qc_vlnplot.pdf"))
print(VlnPlot(so_qc, features = "nCount_RNA",   pt.size = 0) + NoLegend())
print(VlnPlot(so_qc, features = "nFeature_RNA", pt.size = 0) + NoLegend())
print(VlnPlot(so_qc, features = "percent.mt",   pt.size = 0) + NoLegend())
dev.off()

# -- Normalize and variable features -------------------------------------------
so_qc[["RNA"]] <- split(so_qc[["RNA"]], f = so_qc$source)
so_norm        <- NormalizeData(so_qc)
so_norm        <- FindVariableFeatures(so_norm)

# Augment HVGs with curated gut marker genes
gut_markers       <- read.delim(file.path(DATA_DIR, "cell_markers/gut_markers.txt"), header = FALSE)
var_gut_features  <- unique(c(VariableFeatures(so_norm), gut_markers$V1))
VariableFeatures(so_norm) <- var_gut_features

so_norm <- ScaleData(so_norm)
so_norm <- RunPCA(object = so_norm, assay = "RNA", npcs = 30, features = var_gut_features)

# -- Pre-integration UMAP (baseline) ------------------------------------------
so_norm_uninteg <- RunUMAP(so_norm, assay = "RNA", reduction = "pca", dims = 1:30)
coul            <- colorRampPalette(brewer.pal(8, "Set1"))(9)
pdf(file.path(OUTPUT_DIR, "merged_cd/plots/uninteg_umap_9p.pdf"))
print(DimPlot(so_norm_uninteg, cols = alpha(coul, 0.55),
              group.by = "source", raster = FALSE, repel = TRUE))
dev.off()

# -- Harmony integration: source x patient ------------------------------------
so_harm <- IntegrateLayers(
  object         = so_norm,
  method         = HarmonyIntegration,
  assay          = "RNA",
  features       = var_gut_features,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  group.by.vars  = c("source", "patient"),
  theta          = c(1, 1),
  lambda         = c(2, 2),
  verbose        = TRUE
)
so_harm <- FindNeighbors(so_harm, assay = "RNA", reduction = "harmony", dims = 1:30)
so_harm <- FindClusters(so_harm,  resolution = 0.1)
so_harm <- RunUMAP(so_harm, assay = "RNA", reduction = "harmony", dims = 1:30)

# -- Integration UMAP ---------------------------------------------------------
pdf(file.path(OUTPUT_DIR, "merged_cd/plots/integ_umap_9p.pdf"))
print(DimPlot(so_harm, cols = alpha(coul, 0.55),
              group.by = "source", raster = FALSE, repel = TRUE))
dev.off()

# Export UMAP coordinates for Python visualization
umap_df       <- as.data.frame(so_harm@reductions$umap@cell.embeddings)
umap_df$cell  <- rownames(umap_df)
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
meta_cols     <- c("source", "sample", "patient", "condition", "tissue_region",
                   "disease_location", "RNA_snn_res.0.1")
umap_df       <- cbind(umap_df, so_harm@meta.data[rownames(umap_df), meta_cols])
write.csv(umap_df, file.path(OUTPUT_DIR, "merged_cd/9p_integ_p19_upd_anno_umap.csv"))

# Save integrated object
so_harm[["RNA"]] <- JoinLayers(so_harm[["RNA"]])
saveRDS(so_harm, file.path(OUTPUT_DIR, "merged_cd/9p_harmony_integ.RDS"))
