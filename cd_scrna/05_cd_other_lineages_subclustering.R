# CD Atlas: NK, B Cell, Plasma, and Endothelial Subclustering
#
# Iterative subclustering for four minor lineages:
#   - NK cells      (1 iteration from nk_iter0.RDS)
#   - B cells       (2 iterations)
#   - Plasma cells  (2 iterations)
#   - Endothelial   (2 iterations)
#
# Inputs:
#   - 9p_harmony_integ.RDS       : Harmony-integrated atlas
#   - t_cells/nk_iter0.RDS       : NK cells from T cell subclustering
# Outputs:
#   - nk/nk_iter1.RDS
#   - b_cells/b_iter1.RDS
#   - plasma/plasma_iter1.RDS
#   - endothelial/endo_iter1.RDS

library(Seurat)
library(harmony)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR <- "/path/to/cd/scrna/output"
NK_DIR   <- file.path(DATA_DIR, "nk")
B_DIR    <- file.path(DATA_DIR, "b_cells")
PL_DIR   <- file.path(DATA_DIR, "plasma")
EN_DIR   <- file.path(DATA_DIR, "endothelial")
T_DIR    <- file.path(DATA_DIR, "t_cells")
for (d in c(NK_DIR, B_DIR, PL_DIR, EN_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ── Helper: remove CC, RP, MT, HSP genes from HVG list ────────────────────────
filter_hvg <- function(seurat_obj) {
  hvgs <- VariableFeatures(seurat_obj)
  hvgs[!grepl("^RP[SL]|^MT-|^HSP|^DNAJ|^MKI67|^TOP2A|^UBB|^UBC", hvgs)]
}

# ── Load full atlas ────────────────────────────────────────────────────────────
scrna <- readRDS(file.path(DATA_DIR, "9p_harmony_integ.RDS"))

# ============================================================
# 1. NK CELLS
# ============================================================
nk <- readRDS(file.path(T_DIR, "nk_iter0.RDS"))

nk <- NormalizeData(nk)
nk <- FindVariableFeatures(nk, nfeatures = 2000)
VariableFeatures(nk) <- filter_hvg(nk)
nk <- ScaleData(nk)
nk <- RunPCA(nk, npcs = 20)
nk <- RunHarmony(nk,
                 group.by.vars = c("source", "patient"),
                 dims.use = 1:20,
                 theta = c(1, 1), lambda = c(2, 2))
nk <- RunUMAP(nk,   reduction = "harmony", dims = 1:20)
nk <- FindNeighbors(nk, reduction = "harmony", dims = 1:20)
nk <- FindClusters(nk, resolution = 0.2)

remove_nk0 <- c("6")
nk1 <- nk[, !nk$seurat_clusters %in% remove_nk0]

# Final NK subclustering
nk1 <- NormalizeData(nk1)
nk1 <- FindVariableFeatures(nk1, nfeatures = 2000)
VariableFeatures(nk1) <- filter_hvg(nk1)
nk1 <- ScaleData(nk1)
nk1 <- RunPCA(nk1, npcs = 20)
nk1 <- RunHarmony(nk1,
                  group.by.vars = c("source", "patient"),
                  dims.use = 1:20,
                  theta = c(1, 1), lambda = c(2, 2))
nk1 <- RunUMAP(nk1, reduction = "harmony", dims = 1:20)
nk1 <- FindNeighbors(nk1, reduction = "harmony", dims = 1:20)
nk1 <- FindClusters(nk1, resolution = 0.2)
nk1$cell_type <- "NK"
saveRDS(nk1, file.path(NK_DIR, "nk_iter1.RDS"))

# ============================================================
# 2. B CELLS
# ============================================================
b <- scrna[, scrna$RNA_snn_res.0.1 %in% c("2")]  # B cell cluster

b <- NormalizeData(b)
b <- FindVariableFeatures(b, nfeatures = 2000)
VariableFeatures(b) <- filter_hvg(b)
b <- ScaleData(b)
b <- RunPCA(b, npcs = 20)
b <- RunHarmony(b,
                group.by.vars = c("source", "patient"),
                dims.use = 1:20,
                theta = c(1, 1), lambda = c(2, 2))
b <- RunUMAP(b,   reduction = "harmony", dims = 1:20)
b <- FindNeighbors(b, reduction = "harmony", dims = 1:20)
b <- FindClusters(b, resolution = 0.2)
saveRDS(b, file.path(B_DIR, "b_iter0.RDS"))

remove_b0 <- c("7", "8")
b1 <- b[, !b$seurat_clusters %in% remove_b0]

b1 <- NormalizeData(b1)
b1 <- FindVariableFeatures(b1, nfeatures = 2000)
VariableFeatures(b1) <- filter_hvg(b1)
b1 <- ScaleData(b1)
b1 <- RunPCA(b1, npcs = 20)
b1 <- RunHarmony(b1,
                 group.by.vars = c("source", "patient"),
                 dims.use = 1:20,
                 theta = c(1, 1), lambda = c(2, 2))
b1 <- RunUMAP(b1, reduction = "harmony", dims = 1:20)
b1 <- FindNeighbors(b1, reduction = "harmony", dims = 1:20)
b1 <- FindClusters(b1, resolution = 0.2)
b1$cell_type <- "B cell"
saveRDS(b1, file.path(B_DIR, "b_iter1.RDS"))

# ============================================================
# 3. PLASMA CELLS
# ============================================================
plasma <- scrna[, scrna$RNA_snn_res.0.1 %in% c("3")]  # plasma cluster

plasma <- NormalizeData(plasma)
plasma <- FindVariableFeatures(plasma, nfeatures = 2000)
VariableFeatures(plasma) <- filter_hvg(plasma)
plasma <- ScaleData(plasma)
plasma <- RunPCA(plasma, npcs = 20)
plasma <- RunHarmony(plasma,
                     group.by.vars = c("source", "patient"),
                     dims.use = 1:20,
                     theta = c(1, 1), lambda = c(2, 2))
plasma <- RunUMAP(plasma,   reduction = "harmony", dims = 1:20)
plasma <- FindNeighbors(plasma, reduction = "harmony", dims = 1:20)
plasma <- FindClusters(plasma, resolution = 0.2)
saveRDS(plasma, file.path(PL_DIR, "plasma_iter0.RDS"))

remove_pl0 <- c("7")
plasma1 <- plasma[, !plasma$seurat_clusters %in% remove_pl0]

plasma1 <- NormalizeData(plasma1)
plasma1 <- FindVariableFeatures(plasma1, nfeatures = 2000)
VariableFeatures(plasma1) <- filter_hvg(plasma1)
plasma1 <- ScaleData(plasma1)
plasma1 <- RunPCA(plasma1, npcs = 20)
plasma1 <- RunHarmony(plasma1,
                      group.by.vars = c("source", "patient"),
                      dims.use = 1:20,
                      theta = c(1, 1), lambda = c(2, 2))
plasma1 <- RunUMAP(plasma1, reduction = "harmony", dims = 1:20)
plasma1 <- FindNeighbors(plasma1, reduction = "harmony", dims = 1:20)
plasma1 <- FindClusters(plasma1, resolution = 0.2)
plasma1$cell_type <- "Plasma"
saveRDS(plasma1, file.path(PL_DIR, "plasma_iter1.RDS"))

# ============================================================
# 4. ENDOTHELIAL CELLS
# ============================================================
endo <- scrna[, scrna$RNA_snn_res.0.1 %in% c("9")]  # endothelial cluster

endo <- NormalizeData(endo)
endo <- FindVariableFeatures(endo, nfeatures = 2000)
VariableFeatures(endo) <- filter_hvg(endo)
endo <- ScaleData(endo)
endo <- RunPCA(endo, npcs = 20)
endo <- RunHarmony(endo,
                   group.by.vars = c("source", "patient"),
                   dims.use = 1:20,
                   theta = c(1, 1), lambda = c(2, 2))
endo <- RunUMAP(endo,   reduction = "harmony", dims = 1:20)
endo <- FindNeighbors(endo, reduction = "harmony", dims = 1:20)
endo <- FindClusters(endo, resolution = 0.2)
saveRDS(endo, file.path(EN_DIR, "endo_iter0.RDS"))

remove_en0 <- c("6")
endo1 <- endo[, !endo$seurat_clusters %in% remove_en0]

endo1 <- NormalizeData(endo1)
endo1 <- FindVariableFeatures(endo1, nfeatures = 2000)
VariableFeatures(endo1) <- filter_hvg(endo1)
endo1 <- ScaleData(endo1)
endo1 <- RunPCA(endo1, npcs = 20)
endo1 <- RunHarmony(endo1,
                    group.by.vars = c("source", "patient"),
                    dims.use = 1:20,
                    theta = c(1, 1), lambda = c(2, 2))
endo1 <- RunUMAP(endo1, reduction = "harmony", dims = 1:20)
endo1 <- FindNeighbors(endo1, reduction = "harmony", dims = 1:20)
endo1 <- FindClusters(endo1, resolution = 0.2)
endo1$cell_type <- "Endothelial"
saveRDS(endo1, file.path(EN_DIR, "endo_iter1.RDS"))
