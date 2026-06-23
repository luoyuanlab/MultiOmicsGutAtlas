library(Seurat)
library(dplyr)
library(plyr)
library(harmony)
library(SeuratDisk)
library(viridis)
library(CellChat)

# ── Paths ──────────────────────────────────────────────────────────────────────
K1_DIR   <- "/path/to/cosmx_data/uc_1k_6k_wtx/1k/Processed_merged"
OUT_DIR  <- file.path(K1_DIR, "cci_plots")
RDS_DIR  <- file.path(K1_DIR, "seurat_obj")

# ── Load imputed 1K data ───────────────────────────────────────────────────────
cosmx_1k <- readRDS(file.path(RDS_DIR, "1k_4pat_norm_dimreduc_noneg.RDS"))
anno      <- read.csv(file.path(K1_DIR, "anno/auc_cell_type_ywl-b_anno.csv"))
mapper    <- read.csv("/path/to/scrna/uc/cell_type_name_mapper.csv")

anno <- left_join(anno, mapper, by = c("label" = "cell_type_short"))
cosmx_1k$cell_index     <- anno$cells
cosmx_1k$cell_type_1k_fine   <- anno$label
cosmx_1k$cell_type_1k_coarse <- anno$cell_category

# Collapse broad categories; keep Myeloid and Fibroblast subtypes
cosmx_fibmye <- cosmx_1k
cosmx_fibmye@meta.data <- cosmx_fibmye@meta.data %>%
  mutate(label_cci = case_when(
    cell_type_1k_coarse == "Myeloid"     ~ cell_type_1k_fine,
    cell_type_1k_coarse == "Connective"  ~ cell_type_1k_fine,
    cell_type_1k_coarse == "Epithelial"  ~ "Epi",
    cell_type_1k_coarse == "Endothelial" ~ "Endo",
    cell_type_1k_coarse == "T"           ~ "T cell",
    cell_type_1k_coarse == "B"           ~ "B cell",
    TRUE ~ cell_type_1k_coarse
  ))

data.input <- Seurat::GetAssayData(cosmx_fibmye, slot = "data")
meta       <- data.frame(
  cosmx_fibmye@meta.data %>% dplyr::select(fov, label_cci, slide)
)

# Fix missing misc slot in FOV images (Seurat compatibility)
for (nm in names(cosmx_fibmye@images)) {
  img <- cosmx_fibmye@images[[nm]]
  if (inherits(img, "FOV")) {
    tryCatch(slot(img, "misc"),
             error = function(e) { slot(img, "misc") <<- list() })
    cosmx_fibmye@images[[nm]] <- img
  }
}
validObject(cosmx_fibmye)

cosmx_tmp          <- cosmx_fibmye
cosmx_tmp@images   <- list()  # remove FOV payload for spatial locs

# ── Build per-slide spatial coordinates ───────────────────────────────────────
sample_names <- c(
  "UC_batch1_slide1", "UC_batch1_slide2",
  "UC_batch1_slide3", "UC_batch1_slide4",
  "UC_batch2_slide1", "UC_batch2_slide2",
  "UC_batch2_slide3", "UC_batch2_slide4"
)

for (sn in sample_names) assign(sn, subset(cosmx_tmp, subset = slide == sn))

conversion.factor      <- 0.12
combined_spatial_locs  <- NULL
combined_spatial_factors <- NULL

for (sample_name in sample_names) {
  seurat_obj   <- get(sample_name)
  spatial_locs <- data.frame(
    x      = seurat_obj@meta.data$CenterX_global_px,
    y      = seurat_obj@meta.data$CenterY_global_px,
    sample = sample_name,
    row.names = seurat_obj@meta.data$colname
  )
  d         <- computeCellDistance(spatial_locs[, c("x", "y")])
  spot.size <- min(d) * conversion.factor
  spatial_factors <- data.frame(
    sample = sample_name,
    ratio  = conversion.factor,
    tol    = spot.size / 2
  )
  combined_spatial_locs    <- rbind(combined_spatial_locs, spatial_locs)
  combined_spatial_factors <- rbind(combined_spatial_factors, spatial_factors)
}
combined_spatial_locs$sample <- NULL

# ── Create CellChat object ────────────────────────────────────────────────────
cellchat <- createCellChat(
  object          = data.input,
  meta            = meta,
  group.by        = "label_cci",
  datatype        = "spatial",
  coordinates     = as.matrix(combined_spatial_locs),
  spatial.factors = combined_spatial_factors
)

CellChatDB     <- CellChatDB.human
cellchat@DB    <- subsetDB(CellChatDB)
cellchat       <- subsetData(cellchat)

future::plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Main CellChat run (interaction range = 200px, contact range = 10px)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = FALSE, interaction.range = 200,
                              scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 10)

# Subset to Fib/Mye cell types and rename for display
idents_keep <- c(
  "cDC1", "Neutrophil", "HHIP+ Myofib", "GREM2+ Myofib",
  "Mac S+XS-", "cDC2", "Mac S+M+S+", "Activ Fib",
  "Mac S+M+P+", "OGN+RSPO3+ Fib", "RERGL+ Contr Peri", "SELENOP+ Fib",
  "Mac S+SG+", "ADAMDEC1+ Fib", "Mo-Mac", "Inf Mo-Mac", "Cycl Fib",
  "CCR7+ DC", "CD36+ Peri", "Cycl Myeloid", "Trans cDC2/Mac",
  "CXCL5+ Activ Fib", "NRG1+ Crypt Top Fib", "T-interact Fib",
  "VSTM2A+ Crypt Top Fib"
)
cellchat <- subsetCellChat(cellchat, idents.use = idents_keep)
cellchat@meta$label <- gsub("OGN\\+RSPO3\\+ Fib",          "RSPO3+ Fib",       cellchat@meta$label)
cellchat@meta$label <- gsub("VSTM2A\\+ Crypt Top Fib",     "VSTM2A+ CT Fib",  cellchat@meta$label)
cellchat@meta$label <- gsub("NRG1\\+ Crypt Top Fib",       "NRG1+ CT Fib",    cellchat@meta$label)
cellchat@meta$label <- gsub("CXCL5\\+ Activ Fib",          "CXCL5+ Act Fib",  cellchat@meta$label)
cellchat@meta$label <- gsub("RERGL\\+ Contr Peri",          "RERGL+ Peri",     cellchat@meta$label)
cellchat <- setIdent(cellchat, ident.use = "label")

cellchat <- filterCommunication(cellchat, min.cells = 20)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, file.path(RDS_DIR, "cellchat_fib_mye_ywl_v7_imp.RDS"))

# ── IL1 pathway plots (main figure) ───────────────────────────────────────────
pathways.show <- "IL1"

weight.by <- rowSums(cellchat@net$count) + colSums(cellchat@net$count)
pdf(file.path(OUT_DIR, paste0("cci_mye_fib_", pathways.show, "_imp_weighted_chord.pdf")))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord",
                    vertex.size = weight.by, top = 0.2)
dev.off()

pdf(file.path(OUT_DIR, "cci_mye_fib_imp_heatmap.pdf"))
netAnalysis_signalingRole_network(cellchat, width = 12, height = 3, font.size = 9)
dev.off()

# ── Distance sensitivity analysis ─────────────────────────────────────────────
# Tests: 250px range, 15px contact range, original 200px/10px
cellchat_imp250 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                     distance.use = FALSE, interaction.range = 250,
                                     scale.distance = NULL,
                                     contact.dependent = TRUE, contact.range = 10)
cellchat_imp15  <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                     distance.use = FALSE, interaction.range = 200,
                                     scale.distance = NULL,
                                     contact.dependent = TRUE, contact.range = 15)
cellchat_ori    <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                     distance.use = FALSE, interaction.range = 200,
                                     scale.distance = NULL,
                                     contact.dependent = TRUE, contact.range = 10)

relabel <- function(cc) {
  cc <- subsetCellChat(cc, idents.use = idents_keep)
  cc@meta$label <- gsub("OGN\\+RSPO3\\+ Fib",      "RSPO3+ Fib",    cc@meta$label)
  cc@meta$label <- gsub("VSTM2A\\+ Crypt Top Fib", "VSTM2A+ CT Fib", cc@meta$label)
  cc@meta$label <- gsub("NRG1\\+ Crypt Top Fib",   "NRG1+ CT Fib",  cc@meta$label)
  cc@meta$label <- gsub("CXCL5\\+ Activ Fib",       "CXCL5+ Act Fib", cc@meta$label)
  cc@meta$label <- gsub("RERGL\\+ Contr Peri",      "RERGL+ Peri",   cc@meta$label)
  setIdent(cc, ident.use = "label")
}

finalize <- function(cc) {
  cc <- filterCommunication(cc, min.cells = 20)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  netAnalysis_computeCentrality(cc, slot.name = "netP")
}

cellchat_imp250 <- finalize(relabel(cellchat_imp250))
cellchat_imp15  <- finalize(relabel(cellchat_imp15))
cellchat_ori    <- finalize(relabel(cellchat_ori))

# Sensitivity chord plots
for (cfg in list(
  list(obj = cellchat_imp250, tag = "250um_int"),
  list(obj = cellchat_imp15,  tag = "15um_cont"),
  list(obj = cellchat_ori,    tag = "ori")
)) {
  weight.by <- rowSums(cfg$obj@net$count) + colSums(cfg$obj@net$count)
  pdf(file.path(OUT_DIR, paste0("cci_mye_fib_", pathways.show, "_imp_weighted_chord_", cfg$tag, ".pdf")))
  netVisual_aggregate(cfg$obj, signaling = pathways.show, layout = "chord",
                      vertex.size = weight.by, top = 0.2)
  dev.off()
}
