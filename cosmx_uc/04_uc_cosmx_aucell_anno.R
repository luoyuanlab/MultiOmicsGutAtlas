library(Seurat)
library(AUCell)
library(GSEABase)
library(data.table)
set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
# Annotates multiple UC CosMx panels for imputation method comparison (Fig 2)
# Final annotation used for biology: YWL-B v7 imputed 1K (auc_cell_type_ywl-b_anno.csv)
SCRNA_DIR  <- "/path/to/scrna/cd"
WTX_DIR    <- "/path/to/cosmx_data/uc_1k_6k_wtx/whole_trans/Processed_merged"
K1_DIR     <- "/path/to/cosmx_data/uc_1k_6k_wtx/1k/Processed_merged"
K6_DIR     <- "/path/to/cosmx_data/uc_1k_6k_wtx/6k/Processed_merged_upd"
IMP_DIR    <- "/path/to/imputation_model/results/Jenny_cosmx_geneimputation/results/1K_imputation"
BENCH_DIR  <- file.path(WTX_DIR, "v7")

load(file.path(SCRNA_DIR, "cell_type_markers_gene_sets.RData"))

run_aucell <- function(count_csv, save_rdata, save_csv, msg = NULL) {
  if (!is.null(msg)) message(msg)
  dat <- fread(count_csv)
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[[1]]
  dat <- dat[, -1]
  so <- CreateSeuratObject(counts = t(dat))
  cts <- GetAssayData(so, slot = "counts")
  gut_sub <- subsetGeneSets(gsc, rownames(cts))
  set.seed(123)
  cells_AUC <- AUCell_run(cts, gut_sub)
  save(cells_AUC, file = save_rdata)
  mat <- t(getAUC(cells_AUC))
  cell_labels <- colnames(mat)[max.col(mat, "first")]
  result <- data.frame(cells = rownames(mat), label = cell_labels)
  print(table(result$label))
  write.csv(result, save_csv)
  invisible(result)
}

# ── 1. WTX norm counts (no log1p) ─────────────────────────────────────────────
run_aucell(
  count_csv   = file.path(WTX_DIR, "v7/h5ad_obj/norm_counts.csv"),
  save_rdata  = file.path(WTX_DIR, "v7/aucell/norm_counts_aucell.RData"),
  save_csv    = file.path(WTX_DIR, "v7/aucell/norm_counts_aucell.csv"),
  msg = "Annotating: WTX norm counts"
)

# ── 2. YWL-B v7 imputed 1K (final annotation for biology) ────────────────────
log1p_1k <- fread(file.path(K1_DIR, "count_csv/log1p_normalized_count.csv"))
ywl_f1   <- fread(file.path(IMP_DIR, "1K_pred_UC_batch1_f_log1p.csv"))
ywl_f2   <- fread(file.path(IMP_DIR, "1K_pred_UC_batch2_f_log1p.csv"))
ywl_f    <- rbind(as.data.frame(ywl_f1), as.data.frame(ywl_f2))
ywl_f    <- ywl_f[, -1]
ywl_f_full <- cbind(log1p_1k, ywl_f)

ywl_f_full_so <- CreateSeuratObject(counts = t(ywl_f_full))
cts <- GetAssayData(ywl_f_full_so, slot = "counts")
gut_sub <- subsetGeneSets(gsc, rownames(cts))
cells_AUC <- AUCell_run(cts, gut_sub)
save(cells_AUC, file = file.path(K1_DIR, "anno/ywl_imp_f_log_aucell.RData"))
mat <- t(getAUC(cells_AUC))
cell_labels <- colnames(mat)[max.col(mat, "first")]
result <- data.frame(cells = rownames(mat), label = cell_labels)
print(table(result$label))
write.csv(result, file.path(K1_DIR, "anno/auc_cell_type_ywl-b_anno.csv"))
message("Done: YWL-B imputed 1K annotation")

# ── 3. Updated 6K annotation ──────────────────────────────────────────────────
run_aucell(
  count_csv  = file.path(K6_DIR, "count_csv/norm_counts_6175genes.csv"),
  save_rdata = file.path(K6_DIR, "aucell/norm_counts_6175genes_aucell.RData"),
  save_csv   = file.path(K6_DIR, "aucell/norm_counts_6175genes_aucell.csv"),
  msg = "Annotating: updated 6K counts"
)

# ── 4. Imputation benchmarking (negative control validation, Fig 2 / Ext Fig 6G) ──
# Annotates pseudo-6K benchmarks: GimVI, Tangram, STPlus, SPRITE, CellPLM,
# YWL-B, YWL-F, and original 1K/6K/12K panels for comparison.
bench_files <- list(
  list(csv = file.path(BENCH_DIR, "count_csv/half_ori_half_gimvi_imp_log1p.csv"),        tag = "gimvi"),
  list(csv = file.path(BENCH_DIR, "count_csv/half_ori_half_tangram_imp_log1p.csv"),      tag = "tangram"),
  list(csv = file.path(BENCH_DIR, "count_csv/half_ori_half_stplus_imp_log1p.csv"),       tag = "stplus"),
  list(csv = file.path(BENCH_DIR, "count_csv/half_ori_half_sprite_spage_imp_log1p.csv"), tag = "sprite"),
  list(csv = file.path(BENCH_DIR, "count_csv/half_ori_half_cellplm_imp_log1p.csv"),      tag = "cellplm"),
  list(csv = file.path(BENCH_DIR, "count_csv/half_ori_half_ywlb_imp_log1p.csv"),         tag = "ywlb"),
  list(csv = file.path(BENCH_DIR, "count_csv/half_ori_half_ywlf_imp_log1p.csv"),         tag = "ywlf"),
  list(csv = file.path(BENCH_DIR, "count_csv/psuedo_1k.csv"),                            tag = "pseudo_1k"),
  list(csv = file.path(BENCH_DIR, "count_csv/psuedo_6k.csv"),                            tag = "pseudo_6k"),
  list(csv = file.path(BENCH_DIR, "count_csv/psuedo_12k.csv"),                           tag = "pseudo_12k")
)

for (b in bench_files) {
  run_aucell(
    count_csv  = b$csv,
    save_rdata = file.path(BENCH_DIR, "aucell", paste0(b$tag, "_aucell.Rdata")),
    save_csv   = file.path(BENCH_DIR, "aucell", paste0(b$tag, "_aucell.csv")),
    msg = paste("Annotating benchmark:", b$tag)
  )
}
