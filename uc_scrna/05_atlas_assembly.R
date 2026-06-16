library(dplyr)
library(Seurat)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/uc/scrna/output"
OUTPUT_DIR <- "/path/to/uc/scrna/output"

# ── 1. Load final sub-clustered and annotated lineage objects ─────────────────
# Each object has been DietSeurat'd to RNA assay only and carries a 'label' column

uc_mye <- readRDS(file.path(DATA_DIR, "myeloid/mye_iter3.RDS"))
uc_mye <- DietSeurat(uc_mye, assays = "RNA")
uc_mye@meta.data <- uc_mye@meta.data %>% select(orig.ident, sample, Patient, condition)
uc_mye[["label"]] <- uc_mye[["iter3_anno"]]
# Shorten labels for CCI visualizations
uc_mye$label <- dplyr::recode(uc_mye$label,
  "Inflammatory Mo-Mac"        = "Inf Mo-Mac",
  "SELENOP+MMP9+SPP1+ Mac"     = "Mac S+M+S+",
  "SELENOP+XIST- Mac"          = "Mac S+XS-",
  "Transitory cDC2/Mac"        = "Trans cDC2/Mac",
  "SELENOP+MMP9+PLA2G2D+ Mac"  = "Mac S+M+P+",
  "SELENOP+SIGLEC1+ Mac"       = "Mac S+SG+",
  "CCR7+ Mature DC"            = "CCR7+ Mat DC",
  "Cycling"                    = "Cycl Myeloid"
)

uc_fib <- readRDS(file.path(DATA_DIR, "fibroblast/con_iter4.RDS"))
uc_fib <- DietSeurat(uc_fib, assays = "RNA")
uc_fib@meta.data <- uc_fib@meta.data %>% select(orig.ident, sample, Patient, condition, label)

uc_t <- readRDS(file.path(DATA_DIR, "t_iter2.RDS"))
uc_t <- DietSeurat(uc_t, assays = "RNA")
uc_t@meta.data <- uc_t@meta.data %>% select(orig.ident, sample, Patient, condition)
uc_t[["label"]] <- "T cell"

uc_nk <- readRDS(file.path(DATA_DIR, "nk_iter1.RDS"))
uc_nk <- DietSeurat(uc_nk, assays = "RNA")
uc_nk@meta.data <- uc_nk@meta.data %>% select(orig.ident, sample, Patient, condition)
uc_nk[["label"]] <- "NK"

uc_b <- readRDS(file.path(DATA_DIR, "b_iter1.RDS"))
uc_b <- DietSeurat(uc_b, assays = "RNA")
uc_b@meta.data <- uc_b@meta.data %>% select(orig.ident, sample, Patient, condition)
uc_b[["label"]] <- "B cell"

uc_plasma <- readRDS(file.path(DATA_DIR, "plasma_iter1.RDS"))
uc_plasma <- DietSeurat(uc_plasma, assays = "RNA")
uc_plasma@meta.data <- uc_plasma@meta.data %>% select(orig.ident, sample, Patient, condition)
uc_plasma[["label"]] <- "Plasma"

uc_epi <- readRDS(file.path(DATA_DIR, "epi_iter3.RDS"))
uc_epi <- DietSeurat(uc_epi, assays = "RNA")
uc_epi@meta.data <- uc_epi@meta.data %>% select(orig.ident, sample, Patient, condition)
uc_epi[["label"]] <- "Epi"

uc_endo <- readRDS(file.path(DATA_DIR, "endo_iter1.RDS"))
uc_endo <- DietSeurat(uc_endo, assays = "RNA")
uc_endo@meta.data <- uc_endo@meta.data %>% select(orig.ident, sample, Patient, condition)
uc_endo[["label"]] <- "Endo"

# Granulocytes and neurons pulled from the original atlas (not sub-clustered)
uc_orig <- readRDS(file.path(DATA_DIR, "uc_atlas_harmony_integrated.RDS"))

gran <- subset(uc_orig, subset = primary_annotation == "Granulocyte")
gran <- DietSeurat(gran, assays = "RNA")
gran@meta.data <- gran@meta.data %>% select(orig.ident, sample, Patient, condition)
gran[["label"]] <- "Granulocyte"

neu <- subset(uc_orig, subset = primary_annotation == "Neuron")
neu <- DietSeurat(neu, assays = "RNA")
neu@meta.data <- neu@meta.data %>% select(orig.ident, sample, Patient, condition)
neu[["label"]] <- "Neuron"

rm(uc_orig); gc()

# ── 2. Merge all lineages into final annotated atlas ──────────────────────────
uc_atlas <- merge(uc_mye,
                  y = list(uc_fib, uc_t, uc_nk, uc_b, uc_plasma, uc_epi, uc_endo,
                           gran, neu))

saveRDS(uc_atlas, file.path(OUTPUT_DIR, "uc_atlas_annotated.RDS"))
