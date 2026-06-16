library(dplyr)
library(plyr)
library(Seurat)
library(CellChat)
library(future)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/uc/scrna/output"
OUTPUT_DIR <- "/path/to/uc/scrna/output/cci"

plan("multisession", workers = 8)
options(future.globals.maxSize = 50 * 1024^3)

# ── 1. Load final annotated atlas and add cell category metadata ───────────────
uc_scrna <- readRDS(file.path(DATA_DIR, "uc_atlas_annotated.RDS"))
mapper   <- read.csv(file.path(DATA_DIR, "cell_type_name_mapper.csv"))

uc_scrna@meta.data <- uc_scrna@meta.data %>%
  mutate(
    cell_category   = mapvalues(label, mapper$cell_type, mapper$cell_category,    warn_missing = FALSE),
    cell_type_short = mapvalues(label, mapper$cell_type, mapper$cell_type_short,  warn_missing = FALSE)
  )

# ── 2. Myeloid–fibroblast focused CCI (Ext. Fig. 4C-D, 5A-B) ─────────────────
uc_sub <- subset(uc_scrna, subset = cell_category %in% c("Myeloid", "Connective"))

# Shorten labels for readability in chord/circle plots
uc_sub$label_cci <- dplyr::recode(uc_sub$cell_type_short,
  "OGN+RSPO3+ Fib"       = "RSPO3+ Fib",
  "VSTM2A+ Crypt Top Fib" = "VSTM2A+ CT Fib",
  "NRG1+ Crypt Top Fib"  = "NRG1+ CT Fib",
  "CXCL5+ Activ Fib"     = "CXCL5+ Act Fib",
  "RERGL+ Contr Peri"    = "RERGL+ Peri"
)

cellchat_mye_fib <- createCellChat(object = uc_sub,
                                   meta   = data.frame(labels = uc_sub$label_cci,
                                                       row.names = colnames(uc_sub)),
                                   group.by = "label_cci", assay = "RNA")
cellchat_mye_fib@DB <- CellChatDB.human

cellchat_mye_fib <- subsetData(cellchat_mye_fib)
cellchat_mye_fib <- identifyOverExpressedGenes(cellchat_mye_fib)
cellchat_mye_fib <- identifyOverExpressedInteractions(cellchat_mye_fib)
cellchat_mye_fib <- computeCommunProb(cellchat_mye_fib, type = "truncatedMean", trim = 0.1)
cellchat_mye_fib <- filterCommunication(cellchat_mye_fib, min.cells = 50)
cellchat_mye_fib <- computeCommunProbPathway(cellchat_mye_fib)
cellchat_mye_fib <- aggregateNet(cellchat_mye_fib)
cellchat_mye_fib <- netAnalysis_computeCentrality(cellchat_mye_fib, slot.name = "netP")

saveRDS(cellchat_mye_fib, file.path(OUTPUT_DIR, "mye_fib_cellchat.RDS"))

# Export pairwise pathway interactions
net_df <- subsetCommunication(cellchat_mye_fib, slot.name = "netP") %>%
  arrange(pathway, -prob)
write.csv(net_df, file.path(OUTPUT_DIR, "mye_fib_cci_all_pathways.csv"))

# ── 3. All-cell atlas CCI (Ext. Fig. 4C-D incoming/outgoing patterns) ─────────
# Includes all major cell types; used for NMF-based communication pattern discovery
cellchat_all <- createCellChat(object = uc_scrna, group.by = "label", assay = "RNA")
cellchat_all@DB <- CellChatDB.human

cellchat_all <- subsetData(cellchat_all)
cellchat_all <- identifyOverExpressedGenes(cellchat_all)
cellchat_all <- identifyOverExpressedInteractions(cellchat_all)
cellchat_all <- computeCommunProb(cellchat_all, type = "triMean")
cellchat_all <- filterCommunication(cellchat_all, min.cells = 50)
cellchat_all <- computeCommunProbPathway(cellchat_all)
cellchat_all <- aggregateNet(cellchat_all)
cellchat_all <- netAnalysis_computeCentrality(cellchat_all, slot.name = "netP")

saveRDS(cellchat_all, file.path(OUTPUT_DIR, "all_cell_cellchat.RDS"))
