library(Seurat)
library(AUCell)
library(GSEABase)
library(data.table)
set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRNA_DIR   <- "/share/fsmresfiles/UC/scRNA-seq/merged_cd"
RESECT_DIR  <- "/share/fsmresfiles/UC/AtoMx/CD_surgical_resection"

load(file.path(SCRNA_DIR, "cell_type_markers_gene_sets.RData"))

# All 12 resection samples
samples <- c(
  'Cos1_6KDC_strict18888_A5_0715_11_08_2025_14_21_12_760',
  'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
  'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437',
  'Cos1_6KDC_strictCD_R1_B4_0715_11_08_2025_14_20_23_133',
  'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455',
  'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
  'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615',
  'Cos3_6KDC_strictCD_R1_B5_0716_11_08_2025_14_17_56_277',
  'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727',
  'Cos4_6KDC_strictCD_R1_B7_0716_11_08_2025_14_18_41_738',
  'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775',
  'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952'
)

for (sample in samples) {
  message("Starting: ", sample)

  file     <- file.path(RESECT_DIR, sample, "Processed", "norm_counts.csv")
  anno_dir <- file.path(RESECT_DIR, sample)

  # Read and preprocess data
  norm_1k <- fread(file)
  norm_1k <- as.data.frame(norm_1k)
  rownames(norm_1k) <- norm_1k[[1]]
  norm_1k <- norm_1k[, -1]

  # Create Seurat object
  norm_1k_so <- CreateSeuratObject(counts = t(norm_1k))

  # AUCell
  counts <- GetAssayData(object = norm_1k_so, slot = "counts")
  gut_sub <- subsetGeneSets(gsc, rownames(counts))
  cells_AUC <- AUCell_run(counts, gut_sub)

  # Save AUCell results
  save(cells_AUC, file = file.path(anno_dir, "aucell_anno.Rdata"))

  # Extract AUC matrix and get top-scoring labels
  cells_AUC_mat <- t(getAUC(cells_AUC))
  max_col_indices <- max.col(cells_AUC_mat, "first")
  cell_labels <- colnames(cells_AUC_mat)[max_col_indices]

  result <- data.frame(cells = rownames(cells_AUC_mat), label = cell_labels)
  write.csv(result, file.path(anno_dir, "aucell_anno.csv"))

  message("Processed: ", sample)
}
