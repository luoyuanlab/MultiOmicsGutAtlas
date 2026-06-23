library(Seurat)
library(AUCell)
library(GSEABase)
library(data.table)
set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRNA_DIR  <- "/share/fsmresfiles/UC/scRNA-seq/merged_cd"
COSMX_DIR  <- "/share/fsmresfiles/UC/scOmics_model/data/Spatial_Protein_CosMx"

# Shared gene set collection (50 cell types)
load(file.path(SCRNA_DIR, "cell_type_markers_gene_sets.RData"))

# UC multi-omics samples to annotate
samples <- list(
  list(
    csv    = file.path(COSMX_DIR, "count_csv/Cos4_UC_M_batch1_S2_norm_nolog.csv"),
    outdir = file.path(COSMX_DIR, "aucell"),
    tag    = "Cos4_UC_M_batch1_S2"
  ),
  list(
    csv    = file.path(COSMX_DIR, "count_csv/Cos4_UC_M_batch2_S2_norm_nolog.csv"),
    outdir = file.path(COSMX_DIR, "aucell"),
    tag    = "Cos4_UC_M_batch2_S2"
  )
)

for (s in samples) {
  message("Starting: ", s$tag)

  rna <- fread(s$csv)
  rna <- as.data.frame(rna)
  rownames(rna) <- rna[[1]]
  rna <- rna[, -1]

  rna_so <- CreateSeuratObject(counts = t(rna))

  counts  <- GetAssayData(object = rna_so, slot = "counts")
  gut_sub <- subsetGeneSets(gsc, rownames(counts))

  set.seed(123)
  cells_AUC <- AUCell_run(counts, gut_sub)

  cells_AUC_mat   <- t(getAUC(cells_AUC))
  max_col_indices <- max.col(cells_AUC_mat, "first")
  cell_labels     <- colnames(cells_AUC_mat)[max_col_indices]

  result <- data.frame(cells = rownames(cells_AUC_mat), label = cell_labels)
  print(table(result$label))

  write.csv(result, file.path(s$outdir, paste0(s$tag, "_norm_nolog_anno.csv")))

  message("Done: ", s$tag)
}
