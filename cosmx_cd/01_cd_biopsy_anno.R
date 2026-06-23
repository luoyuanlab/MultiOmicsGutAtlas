library(Seurat)
library(AUCell)
library(GSEABase)
library(data.table)
set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
# All CD 1K biopsy samples available: CD_B1, CD_B3, CD_B4 (primary for figures),
# CD_B5 (validation), plus non-inflamed/HC controls.
# Primary sample for figures: CD_B4; Validation: CD_B5
SCRNA_DIR  <- "/path/to/scrna/cd"
COSMX_DIR  <- "/path/to/cosmx_data/uc_serial_cuts/1000plex"

load(file.path(SCRNA_DIR, "cell_type_markers_gene_sets.RData"))

# Input files: edit this vector to add/remove samples
files <- c(
  file.path(COSMX_DIR, "CD_batch1_dir/h5ad_obj/CD_B1_norm_counts.csv"),
  file.path(COSMX_DIR, "CD_batch2_dir/h5ad_obj/CD_B3_norm_counts.csv")
)

annotate_one_file <- function(file, gsc) {
  message("Starting: ", file)

  # Read and preprocess
  norm_1k <- fread(file)
  norm_1k <- as.data.frame(norm_1k)
  rownames(norm_1k) <- norm_1k[[1]]
  norm_1k <- norm_1k[, -1, drop = FALSE]

  # Create Seurat object (cells are columns in Seurat -> transpose)
  norm_1k_so <- CreateSeuratObject(counts = t(norm_1k))

  # AUCell
  counts <- GetAssayData(object = norm_1k_so, slot = "counts")
  gut_sub <- subsetGeneSets(gsc, rownames(counts))
  cells_AUC <- AUCell_run(counts, gut_sub)

  # Output directory = same folder as the input CSV
  anno_dir <- dirname(file)

  # Save AUCell object
  save(cells_AUC, file = file.path(anno_dir, "aucell_anno.Rdata"))

  # Extract AUC matrix and assign top label per cell
  cells_AUC_mat <- t(getAUC(cells_AUC))  # cells x geneSets
  max_col_indices <- max.col(cells_AUC_mat, ties.method = "first")
  cell_labels <- colnames(cells_AUC_mat)[max_col_indices]

  result <- data.frame(
    cells = rownames(cells_AUC_mat),
    label = cell_labels,
    stringsAsFactors = FALSE
  )

  write.csv(result, file.path(anno_dir, "aucell_anno.csv"), row.names = FALSE)

  message("Done: ", file)
  invisible(result)
}

# Run for all files
for (f in files) {
  annotate_one_file(f, gsc)
}
