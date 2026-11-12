# UC CosMx 1K: differential expression between inflamed rectal UC (4 samples) and other tissue (12 samples) (Fig 3A).
#
# Inputs:  /share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/ sample_list_16.txt and comparison/<sample_id>_smi_so_qcnorm.RDS
#          (QC'd Seurat objects from a QC pipeline that is not in this repository)
# Output:  comparison/cosmx_merged_category_DGE.txt (FindAllMarkers, SCT); the Fig 3A volcano uses its "Rectum, Inflamed UC" rows
# Run:     Rscript <file>, or RUN_SECTIONS=<id>[,<id>] Rscript <file>
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}


# Section F3A_1k_de: SCTransform per sample, merge, FindAllMarkers by tissue/disease category
if (.run("F3A_1k_de")) {
DATA_DIR_1K <- "/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex"   # folder holding sample_list_16.txt and comparison/
setwd(DATA_DIR_1K)

library(data.table)
library(Seurat)

sample_list <- fread("sample_list_16.txt")
sample_list$category <- paste0(sample_list$tissue_type, ", ", sample_list$disease_state, " ", sample_list$disease)

subsetList <- function(myList, elementNames) {
	sapply(elementNames, FUN=function(x) myList[[x]])
}

so_list <- list()
for (i in 1:nrow(sample_list)) {
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

### Category (tissue type, disease state, disease type) DGE
annotations <- so[["category"]]$category
names(annotations) <- rownames(so[["category"]])
Idents(so) <- annotations
category_DGE <- FindAllMarkers(so, assay = "SCT")
fwrite(category_DGE, "comparison/cosmx_merged_category_DGE.txt", sep = "\t")
}  # end section F3A_1k_de
