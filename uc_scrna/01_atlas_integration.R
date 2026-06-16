library(Matrix)
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(RColorBrewer)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/uc/scrna"          # directory with per-dataset SCT objects
OUTPUT_DIR <- "/path/to/uc/scrna/output"   # where integrated objects are saved
MARKER_DIR <- "/path/to/cell_markers"      # PanglaoDB and CellMarker files

# ── 1. Build gut-enriched feature set ─────────────────────────────────────────
# Combine HVGs from an initial integration pass with curated gut/immune markers
# from CellMarker and PanglaoDB to improve cell type resolution

cell_marker <- read.csv(file.path(MARKER_DIR, "Cell_marker_Human.csv"))
panglao     <- read.table(file.path(MARKER_DIR, "PanglaoDB_markers.tsv"),
                          sep = "\t", header = TRUE)

cell_marker_gut <- cell_marker %>%
  filter(tissue_class %in% c("Gastrointestinal tract", "Intestine")) %>%
  pull(Symbol)
cell_marker_gut <- cell_marker_gut[nzchar(cell_marker_gut)]

panglao_gut <- panglao %>%
  filter(organ %in% c("GI tract", "Immune system")) %>%
  pull(official.gene.symbol)

merged_gut_marker <- unique(c(panglao_gut, cell_marker_gut))

# Load previously merged Seurat list (per-dataset SCT objects after QC)
# Each dataset was independently: normalized, log-transformed (nFeature >= 500,
# nCount >= 500, percent.mt < 25%), variable features selected
plist <- readRDS(file.path(DATA_DIR, "Seurat_paper_sct.RDS"))

# Select integration features: HVGs shared across datasets
features <- SelectIntegrationFeatures(plist, nfeatures = 3000, method = "glmGamPoi")

# Combine HVGs with gut marker genes (4601 total features)
merged_var_gene_new_marker <- unique(c(features, merged_gut_marker))

# Keep only genes present in the merged object
so_merged <- merge(plist[[1]], y = plist[2:length(plist)],
                   project = "uc", merge.data = TRUE)
feature_upd <- intersect(merged_var_gene_new_marker, rownames(so_merged))

# ── 2. Harmony integration ─────────────────────────────────────────────────────
# Batch correction across datasets (orig.ident) and samples simultaneously
DefaultAssay(so_merged) <- "RNA"
so_merged <- NormalizeData(so_merged)
VariableFeatures(so_merged) <- feature_upd
so_merged <- ScaleData(so_merged)
so_merged <- RunPCA(so_merged, assay = "RNA", npcs = 30, features = feature_upd)
so_merged <- RunHarmony(so_merged,
                        assay.use    = "RNA",
                        reduction    = "pca",
                        dims.use     = 1:30,
                        group.by.vars = c("orig.ident", "sample"),
                        theta        = c(1, 1),
                        lambda       = c(2, 2))

# ── 3. Dimensionality reduction and primary clustering ─────────────────────────
so_merged <- RunUMAP(so_merged, assay = "RNA", reduction = "harmony", dims = 1:30)
so_merged <- FindNeighbors(so_merged, assay = "RNA", reduction = "harmony", dims = 1:30)
so_merged <- FindClusters(so_merged, resolution = c(0.1, 0.2, 0.3, 0.4))

saveRDS(so_merged, file.path(OUTPUT_DIR, "uc_atlas_harmony_integrated.RDS"))
