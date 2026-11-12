# CD CosMx and CITE-seq: Seurat label transfer of scRNA-seq Treg pseudotime bins (0-4).
# Inputs:  objects/label_transfer/ reference CSVs (08); Treg raw counts of the WTX biopsies (02), resections (05) and CITE-seq (cite-tcr/04)
# Outputs: objects/label_transfer/treg_pt_bin.RDS (reference) and <cohort>_treg_pt_transfer_labels.csv per query
#          (read by 10_treg_suppression_score.py and cite-tcr/05_treg_tcr_repertoire.py)
# Sections: treg_reference_build -> treg_pt_bin_transfer (optional cohort names as arguments; default all)
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {   # all sections run in file order unless RUN_SECTIONS=id1,id2,... is set
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}
library(Seurat)
OUT_DIR <- "/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer"


# ---- Section treg_reference_build: Seurat reference of the scRNA-seq Tregs (NormalizeData, 2000 vst HVGs, PCA 30, UMAP)
if (.run("treg_reference_build")) {
scrna_counts   <- t(read.csv(file.path(OUT_DIR, "treg_raw_counts.csv"), row.names = 1))
scrna_metadata <- read.csv(file.path(OUT_DIR, "treg_pt_bin.csv"), row.names = 1)

scrna <- CreateSeuratObject(
  counts    = as.matrix(scrna_counts),
  meta.data = scrna_metadata,
  assay     = "RNA"
)
scrna <- NormalizeData(scrna)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000)
scrna <- ScaleData(scrna)
scrna <- RunPCA(scrna, npcs = 30, verbose = FALSE)
scrna <- RunUMAP(scrna, reduction = "pca", dims = 1:30)
saveRDS(scrna, file.path(OUT_DIR, "treg_pt_bin.RDS"))
}  # end section treg_reference_build


# ---- Section treg_pt_bin_transfer: pseudotime bins onto each query (query PCA 30, UMAP; FindTransferAnchors dims 1:30, reference PCA,
#      LogNormalize; TransferData weighted by the query PCA). Queries: 1,204 biopsy AUCell Tregs, 6,554 resection strict Tregs,
#      745 CITE-seq AUCell Tregs.
if (.run("treg_pt_bin_transfer")) {
WTX  <- "/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged"
RES  <- "/share/fsmresfiles/UC/AtoMx/CD_surgical_resection/combined"
CITE <- "/share/fsmresfiles/UC/Schipma/cite-seq_combined/annotation_for_cosmx/third_round/tcell"

queries <- list(
  biopsy_alltreg   = file.path(WTX, "count_csv/all_treg_raw_counts.csv"),
  resection_strict = file.path(RES, "raw_counts_strict_tregs.csv"),
  cite_notstrict   = file.path(CITE, "notstrict_treg_raw_counts.csv")
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) queries <- queries[args]

scrna <- readRDS(file.path(OUT_DIR, "treg_pt_bin.RDS"))
scrna$pt_bin <- as.character(scrna$pt_bin)

for (cohort in names(queries)) {
  message("==== ", cohort, " ====")
  cosmx_counts <- t(read.csv(queries[[cohort]], row.names = 1))
  cosmx <- CreateSeuratObject(counts = as.matrix(cosmx_counts), assay = "RNA")
  cosmx <- NormalizeData(cosmx)
  cosmx <- FindVariableFeatures(cosmx, selection.method = "vst", nfeatures = 2000)
  cosmx <- ScaleData(cosmx)
  cosmx <- RunPCA(cosmx, npcs = 30, verbose = FALSE)
  cosmx <- RunUMAP(cosmx, reduction = "pca", dims = 1:30, verbose = FALSE)   # kept for identical random state

  transfer_anchors <- FindTransferAnchors(
    reference            = scrna,
    query                = cosmx,
    dims                 = 1:30,
    reference.reduction  = "pca",
    normalization.method = "LogNormalize"
  )
  predictions <- TransferData(
    anchorset        = transfer_anchors,
    refdata          = scrna$pt_bin,
    dims             = 1:30,
    weight.reduction = cosmx[["pca"]]
  )
  cosmx <- AddMetaData(cosmx, metadata = predictions)
  print(table(cosmx$predicted.id))
  write.csv(cosmx@meta.data, file.path(OUT_DIR, paste0(cohort, "_treg_pt_transfer_labels.csv")))
}
}  # end section treg_pt_bin_transfer
