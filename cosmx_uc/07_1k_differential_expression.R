# UC CosMx 1K: differential expression between inflamed rectal UC and other tissue.
#
# Panels made
#   Fig 3A  DE between inflamed rectal UC (4 samples) and other tissue (12 samples), CosMx 1K
#
# Inputs
#   DATA_DIR_1K = /share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/
#     sample_list_16.txt                        16 CosMx 1K samples (sample_id, block, tissue_type, disease, disease_state)
#     comparison/<sample_id>_smi_so_qcnorm.RDS  QC'd, normalised Seurat objects, one per sample (made by a separate QC
#                                               pipeline that is not part of this repository)
#
# Outputs
#   DATA_DIR_1K/comparison/cosmx_merged_category_DGE.txt   Seurat FindAllMarkers (SCT) table, one block per category
#   The volcano plot is drawn from the "Rectum, Inflamed UC" rows of that table; no image file is written.
#
# Run order of sections
#   F3A_1k_de
#
# How to run
#   Rscript <this file>                          run every section in the order below
#   RUN_SECTIONS=<id>[,<id>...] Rscript <file>   run only the listed sections (still in file order)
#   Each section runs one analysis step inside `if (.run("<id>")) { ... }` and reads its inputs from disk.
#   Paths are the lab's absolute paths (to be replaced by placeholders before public release).
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {   # all sections run in file order unless RUN_SECTIONS=id1,id2,... is set
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}


# ======================================================================================================================
# Section F3A_1k_de: SCTransform, merge and FindAllMarkers by tissue/disease category across the 16 CosMx 1K samples
# Panels: Fig 3A
# Note: the section reads and writes paths relative to DATA_DIR_1K (set below); setwd() makes it independent of the launch folder
# ======================================================================================================================
if (.run("F3A_1k_de")) {
DATA_DIR_1K <- "/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex"   # folder holding sample_list_16.txt and comparison/
setwd(DATA_DIR_1K)
# cd /share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex

library(data.table)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

sample_list <- fread("sample_list_16.txt")
sample_list$category <- paste0(sample_list$tissue_type, ", ", sample_list$disease_state, " ", sample_list$disease)
table(sample_list$category)
#      Ileum, Inflamed CD  Ileum, Non-inflamed CD  Ileum, Non-inflamed UC
#                       2                       2                       2
#     Rectum, Inflamed CD     Rectum, Inflamed UC Rectum, Non-inflamed CD
#                       2                       4                       2
# Rectum, Non-inflamed UC
#                       2

subsetList <- function(myList, elementNames) {
	sapply(elementNames, FUN=function(x) myList[[x]])
}

so_list <- list()
for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	so <- readRDS(paste0("comparison/", as.character(sample_id), "_smi_so_qcnorm.RDS"))
	so@meta.data$sample_id <- sample_id
	so@meta.data$block <- sample_list$block[i]
	so@meta.data$tissue_type <- sample_list$tissue_type[i]
	so@meta.data$disease <- sample_list$disease[i]
	so@meta.data$disease_state <- sample_list$disease_state[i]
	so@meta.data$category <- sample_list$category[i]
	
	so <- SCTransform(so, assay = "Nanostring", return.only.var.genes = FALSE, clip.range = c(-10, 10), verbose = FALSE)
	
	so_list[[sample_id]] <- so
}
so <- merge(x = so_list[[1]], y = c(subsetList(so_list, seq(2, length(so_list)))), add.cell.ids = as.character(1:16), project = "CosMx")

so <- PrepSCTFindMarkers(object = so)
VariableFeatures(so) <- rownames(so) # set all the genes to be variable because we want to use all the genes
so <- RunPCA(so, assay = "SCT", verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = 1:40, verbose = FALSE)
so <- FindNeighbors(so, reduction = "pca", dims = 1:40, verbose = FALSE)
so <- FindClusters(so, dims = 1:40, resolution = 0.3, verbose = FALSE)

### Category (tissue type, disease state, disease type) DGE
annotations <- so[["category"]]$category
names(annotations) <- rownames(so[["category"]])
Idents(so) <- annotations
category_DGE <- FindAllMarkers(so, assay = "SCT")
fwrite(category_DGE, "comparison/cosmx_merged_category_DGE.txt", sep = "\t")



### Volcano plot (Fig 3A) uses this subset
category_DGE <- category_DGE[category_DGE$cluster == "Rectum, Inflamed UC",]
}  # end section F3A_1k_de
