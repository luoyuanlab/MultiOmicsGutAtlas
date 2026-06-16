# CD Atlas: Final Assembly of All Lineages
#
# Reads per-lineage subclustered Seurat objects, selects a common
# set of metadata columns, assigns final cell_type labels, and
# merges everything into the annotated CD atlas.
#
# Inputs (from subclustering scripts 02–06):
#   - myeloid/mye_iter4.RDS
#   - connective/con_iter2.RDS
#   - t_cells/tcell_iter3.RDS
#   - nk/nk_iter1.RDS
#   - b_cells/b_iter1.RDS
#   - plasma/plasma_iter1.RDS
#   - endothelial/endo_iter1.RDS
#   - epithelial/epi_iter2.RDS    (iter1_anno column used)
# Output:
#   - merged_cd/cleaned_annoed_all_cell_types.RDS  : annotated atlas
#   - merged_cd/meta_cleaned_annoed_all_cell_types.csv

library(Seurat)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR <- "/path/to/cd/scrna/output"

# ── Load all lineage objects ───────────────────────────────────────────────────
myeloid <- readRDS(file.path(DATA_DIR, "myeloid/mye_iter4.RDS"))
con     <- readRDS(file.path(DATA_DIR, "connective/con_iter2.RDS"))
tcell   <- readRDS(file.path(DATA_DIR, "t_cells/tcell_iter3.RDS"))
nk      <- readRDS(file.path(DATA_DIR, "nk/nk_iter1.RDS"))
b       <- readRDS(file.path(DATA_DIR, "b_cells/b_iter1.RDS"))
plasma  <- readRDS(file.path(DATA_DIR, "plasma/plasma_iter1.RDS"))
endo    <- readRDS(file.path(DATA_DIR, "endothelial/endo_iter1.RDS"))
epi     <- readRDS(file.path(DATA_DIR, "epithelial/epi_iter2.RDS"))

# ── Common metadata columns to retain ─────────────────────────────────────────
KEEP_COLS <- c("orig.ident", "nCount_RNA", "nFeature_RNA",
               "source", "sample", "patient",
               "condition", "tissue_region", "disease_location",
               "treatment", "percent.mt", "label")

# ── Helper: slim down a Seurat object for merging ────────────────────────────
slim <- function(so, anno_col) {
  so  <- DietSeurat(so, assays = "RNA")
  so$label <- so@meta.data[[anno_col]]
  # Keep only standard metadata
  keep <- intersect(KEEP_COLS, colnames(so@meta.data))
  so@meta.data <- so@meta.data[, keep, drop = FALSE]
  so
}

myeloid <- slim(myeloid, "iter4_anno")
con     <- slim(con,     "iter2_anno")
tcell   <- slim(tcell,   "iter3_anno")
nk      <- slim(nk,      "cell_type")
b       <- slim(b,       "cell_type")
plasma  <- slim(plasma,  "cell_type")
endo    <- slim(endo,    "cell_type")
epi     <- slim(epi,     "iter1_anno")

# ── Merge all lineages ────────────────────────────────────────────────────────
atlas <- merge(myeloid, con)
atlas <- merge(atlas, tcell)
atlas <- merge(atlas, nk)
atlas <- merge(atlas, b)
atlas <- merge(atlas, plasma)
atlas <- merge(atlas, endo)
atlas <- merge(atlas, epi)
atlas[["RNA"]] <- JoinLayers(atlas[["RNA"]])

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(atlas, file.path(DATA_DIR, "merged_cd/cleaned_annoed_all_cell_types.RDS"))
write.csv(atlas@meta.data,
          file.path(DATA_DIR, "merged_cd/meta_cleaned_annoed_all_cell_types.csv"))

cat("Atlas assembled:", ncol(atlas), "cells,", nrow(atlas), "genes\n")
print(table(atlas$label))
