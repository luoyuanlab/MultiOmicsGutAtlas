# CD CosMx and CITE-seq: Seurat label transfer of scRNA-seq Treg pseudotime bins.
#
# Panels made
#   none (Treg pseudotime-bin labels for the WTX biopsies, the 6K resections and the CITE-seq Tregs; used by
#   10_treg_suppression_score.py and cite-tcr/05_treg_tcr_repertoire.py)
#
# Inputs
#   /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer/ (reference CSVs from 08_treg_reference_export.py)
#   /share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged/count_csv/all_treg_raw_counts.csv (WTX biopsy Treg counts,
#     02_wtx_biopsy_treg_labels_liana.ipynb)
#   /share/fsmresfiles/UC/AtoMx/CD_surgical_resection/combined/raw_counts_strict_tregs.csv (resection strict-Treg counts, 05_resection_treg_counts.ipynb)
#   /share/fsmresfiles/UC/Schipma/cite-seq_combined/annotation_for_cosmx/third_round/tcell/notstrict_treg_raw_counts.csv
#     (CITE-seq Treg counts, cite-tcr/04_treg_subset_export.ipynb)
#
# Outputs
#   /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer/ (reference RDS and <cohort>_treg_pt_transfer_labels_rowfix.csv per query)
#
# Run order of sections
#   08_treg_reference_export.py -> treg_reference_build -> treg_pt_bin_transfer -> 10_treg_suppression_score.py
#
# Notes
#   treg_pt_bin_transfer takes optional cohort names as arguments (default: all cohorts).
#
# How to run
#   Rscript <this file>                          run every section in the order below
#   RUN_SECTIONS=<id>[,<id>...] Rscript <file>   run only the listed sections (still in file order)
#   Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
#   upstream outputs exist.
#   Paths are the lab's absolute paths; edit them to match the local data location.
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {   # all sections run in file order unless RUN_SECTIONS=id1,id2,... is set
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}


# ======================================================================================================================
# Section treg_reference_build: Seurat reference of the scRNA-seq Tregs with pseudotime bins
# Panels: label-transfer reference
# ======================================================================================================================
if (.run("treg_reference_build")) {
# Builds the Seurat reference from the scRNA-seq Treg raw counts and pseudotime bins (0-4):
#   inputs  objects/label_transfer/treg_raw_counts_rowfix.csv + treg_pt_bin_rowfix.csv (08_treg_reference_export.py)
#   steps   NormalizeData -> vst 2000 HVG -> ScaleData -> PCA 30 -> UMAP
#   output  objects/label_transfer/treg_pt_bin_rowfix.RDS
library(Seurat)
library(dplyr)
library(ggplot2)

SCRNA_DIR <- "/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer"

# ── Load and prepare the scRNA-seq reference ───────────────────────────────────
scrna_counts   <- read.csv(file.path(SCRNA_DIR, "treg_raw_counts_rowfix.csv"), row.names = 1)
scrna_counts   <- t(scrna_counts)
scrna_metadata <- read.csv(file.path(SCRNA_DIR, "treg_pt_bin_rowfix.csv"), row.names = 1)

print(paste("scRNA genes:", nrow(scrna_counts)))
print(paste("scRNA cells:", ncol(scrna_counts)))
print(paste("Metadata rows:", nrow(scrna_metadata)))
stopifnot(identical(colnames(scrna_counts), make.names(rownames(scrna_metadata))) ||
          identical(colnames(scrna_counts), rownames(scrna_metadata)))
colnames(scrna_counts) <- rownames(scrna_metadata)   # read.csv keeps row names verbatim; guard only

scrna <- CreateSeuratObject(
  counts    = as.matrix(scrna_counts),
  meta.data = scrna_metadata,
  assay     = "RNA"
)

print("Pseudotime bin distribution:")
print(table(scrna$pt_bin))

scrna <- NormalizeData(scrna)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000)
scrna <- ScaleData(scrna)
scrna <- RunPCA(scrna, npcs = 30, verbose = FALSE)
scrna <- RunUMAP(scrna, reduction = "pca", dims = 1:30)
saveRDS(scrna, file.path(SCRNA_DIR, "treg_pt_bin_rowfix.RDS"))
message("saved ", file.path(SCRNA_DIR, "treg_pt_bin_rowfix.RDS"))
}  # end section treg_reference_build


# ======================================================================================================================
# Section treg_pt_bin_transfer: transfer of Treg pseudotime bins (0-4) onto the CosMx and CITE-seq Treg queries
# Panels: Treg pseudotime-bin labels for WTX biopsies, resections and CITE-seq
# ======================================================================================================================
if (.run("treg_pt_bin_transfer")) {
# Transfers the Treg pseudotime bins (0-4) of the scRNA-seq reference onto every query:
#   query    NormalizeData -> vst 2000 HVG -> ScaleData -> PCA 30;
#   anchors  FindTransferAnchors(dims 1:30, reference.reduction "pca", LogNormalize);
#   labels   TransferData(refdata = pt_bin, dims 1:30, weight.reduction = query pca).
# Queries (raw count CSVs):
#   biopsy_alltreg   CD WTX biopsies, all 1,204 AUCell Tregs
#   resection_strict CD resections, 6,554 strict Tregs
#   cite_notstrict   CITE-seq, 745 AUCell Tregs (TCR analysis)
# Outputs (objects/label_transfer/): <cohort>_treg_pt_transfer_labels_rowfix.csv (+ <cohort>_treg_ptbin_anno_rowfix.RDS)
# Arguments: optional cohort names, e.g. Rscript <this file> resection_strict
library(Seurat)
library(dplyr)

OUT_DIR <- "/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer"
WTX     <- "/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged"
RES     <- "/share/fsmresfiles/UC/AtoMx/CD_surgical_resection/combined"
CITE    <- "/share/fsmresfiles/UC/Schipma/cite-seq_combined/annotation_for_cosmx/third_round/tcell"

queries <- list(
  biopsy_alltreg   = list(csv = file.path(WTX, "count_csv/all_treg_raw_counts.csv"), subset = NULL),
  resection_strict = list(csv = file.path(RES, "raw_counts_strict_tregs.csv"), subset = NULL),
  cite_notstrict   = list(csv = file.path(CITE, "notstrict_treg_raw_counts.csv"), subset = NULL)
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) queries <- queries[args]

scrna <- readRDS(file.path(OUT_DIR, "treg_pt_bin_rowfix.RDS"))
scrna$pt_bin <- as.character(scrna$pt_bin)

for (cohort in names(queries)) {
  q <- queries[[cohort]]
  message("\n==== ", cohort, " ====")
  # ── Load and prepare the query ──
  cosmx_counts <- read.csv(q$csv, row.names = 1)
  if (!is.null(q$subset)) {
    keep <- rownames(read.csv(q$subset, row.names = 1))
    stopifnot(all(keep %in% rownames(cosmx_counts)))
    cosmx_counts <- cosmx_counts[keep, ]
  }
  cosmx_counts <- t(cosmx_counts)
  print(paste("Query genes:", nrow(cosmx_counts)))
  print(paste("Query cells:", ncol(cosmx_counts)))

  cosmx <- CreateSeuratObject(counts = as.matrix(cosmx_counts), assay = "RNA")
  cosmx <- NormalizeData(cosmx)
  cosmx <- FindVariableFeatures(cosmx, selection.method = "vst", nfeatures = 2000)
  cosmx <- ScaleData(cosmx)
  cosmx <- RunPCA(cosmx, npcs = 30, verbose = FALSE)
  cosmx <- RunUMAP(cosmx, reduction = "pca", dims = 1:30, verbose = FALSE)

  # ── Gene overlap ──
  common_genes <- intersect(rownames(scrna), rownames(cosmx))
  print(paste("Common genes:", length(common_genes)))
  print(paste("Reference HVGs present in query:", sum(VariableFeatures(scrna) %in% rownames(cosmx)), "/ 2000"))

  # ── Label transfer ──
  transfer_anchors <- FindTransferAnchors(
    reference            = scrna,
    query                = cosmx,
    dims                 = 1:30,
    reference.reduction  = "pca",
    normalization.method = "LogNormalize"
  )
  print(paste("Number of anchors found:", nrow(transfer_anchors@anchors)))

  predictions <- TransferData(
    anchorset        = transfer_anchors,
    refdata          = scrna$pt_bin,
    dims             = 1:30,
    weight.reduction = cosmx[["pca"]]
  )
  cosmx <- AddMetaData(cosmx, metadata = predictions)
  saveRDS(cosmx, file.path(OUT_DIR, paste0(cohort, "_treg_ptbin_anno_rowfix.RDS")))

  print("Prediction score summary:")
  print(summary(cosmx$prediction.score.max))
  print(table(cosmx$predicted.id))

  # ── Export ──
  write.csv(cosmx@meta.data, file.path(OUT_DIR, paste0(cohort, "_treg_pt_transfer_labels_rowfix.csv")))
  summary_stats <- cosmx@meta.data %>%
    group_by(predicted.id) %>%
    summarise(n_cells = n(), mean_score = mean(prediction.score.max),
              median_score = median(prediction.score.max),
              pct_high_conf = sum(prediction.score.max > 0.5) / n() * 100) %>%
    arrange(predicted.id)
  print(as.data.frame(summary_stats))
  message("High confidence (>0.5): ", sum(cosmx$prediction.score.max > 0.5), " / ", ncol(cosmx))
}
}  # end section treg_pt_bin_transfer
