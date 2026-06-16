# CD Atlas: NicheNet Immune → Treg Ligand–Receptor Prioritization
#
# Prioritizes ligand–receptor pairs from immune sender cells
# (B, DC, Granulocyte, Macrophage, NK, Other T, Plasma) to
# Tregs binned along a pseudotime trajectory (Early vs. Late).
#
# Priority score: geometric mean of five percentile-ranked components
#   1. Ligand activity (AUPR against Late-Treg target genes, α = 2)
#   2. Destabilization activity (AUPR against Early-Treg targets, β = 1.5)
#   3. Ligand availability in senders (γ = 1)
#   4. Receptor early prevalence in Early-Tregs (δ = 1)
#   5. Receptor late up-regulation in Late- vs Early-Tregs (η = 1.5)
#
# Outputs:
#   - nichenet_destab_late_treg_ligand_strength_imm_logpr.csv
#   - nichenet_best_sender_per_ligand.csv
#   - plots/nichenet_chord_diagram.pdf
#   - plots/nichenet_heatmap.pdf
#   - plots/nichenet_weight_sensitivity.pdf

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(purrr)
library(PRROC)
library(scales)
library(nichenetr)
library(tidytext)
library(ggforce)
library(stringr)
library(forcats)
library(ggplot2)
library(arrow)
library(circlize)
library(grid)

eps <- 1e-6

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR    <- "/path/to/cd/scrna/output"
NICHENET_DB <- "/path/to/nichenet"
OUTPUT_DIR  <- file.path(DATA_DIR, "nichenet")
dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)

# ── Load Treg + immune sender expression data ──────────────────────────────────
# Cells stratified: Tregs binned into pseudotime bins B0–B4;
# immune senders: B, DC, Granulocyte, Macrophage, NK, Other T, Plasma
meta   <- as.data.frame(read_parquet(file.path(DATA_DIR, "strict_treg_imm_meta.parquet")))
counts <- as.data.frame(read_parquet(file.path(DATA_DIR, "strict_treg_imm_counts.parquet")))

if ("__index_level_0__" %in% colnames(counts)) counts <- counts[, colnames(counts) != "__index_level_0__"]
if ("__index_level_0__" %in% colnames(meta))   meta   <- meta[, colnames(meta)   != "__index_level_0__"]

rownames(counts) <- meta$index

# ── Build Seurat object and assign Treg Early / Late labels ───────────────────
obj <- CreateSeuratObject(counts = t(counts))
obj$ct_lr <- meta$ct_lr

treg_mask <- grepl("^Treg_B[0-4]$", obj$ct_lr)
pt_bin    <- sub("^Treg_(B[0-4])$", "\\1", obj$ct_lr)
pt_bin[!treg_mask] <- NA_character_

EARLY_BINS <- c("B0", "B1")
LATE_BINS  <- c("B2", "B3", "B4")

nn_id <- obj$ct_lr
nn_id[treg_mask & pt_bin %in% EARLY_BINS] <- "Treg_Early"
nn_id[treg_mask & pt_bin %in% LATE_BINS]  <- "Treg_Late"

obj$nn_id <- factor(nn_id)
Idents(obj) <- obj$nn_id

# Normalize
obj <- NormalizeData(obj, layer = "counts", new.layer = "data")
obj <- ScaleData(obj, layer = "data", features = rownames(obj))

mat          <- GetAssayData(obj, slot = "data")
genes_avail  <- rownames(mat)
early_cells  <- colnames(obj)[obj$nn_id == "Treg_Early"]
late_cells   <- colnames(obj)[obj$nn_id == "Treg_Late"]

# ── Load NicheNet databases ────────────────────────────────────────────────────
lr_network          <- readRDS(file.path(NICHENET_DB, "lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file.path(NICHENET_DB, "ligand_target_matrix_nsga2r_final.rds"))
weighted_networks   <- readRDS(file.path(NICHENET_DB, "weighted_networks_nsga2r_final.rds"))

ligands_use   <- intersect(unique(lr_network$from), genes_avail)
receptors_use <- intersect(unique(lr_network$to),   genes_avail)

SENDER_TYPES <- c("B", "DC", "Granulocyte", "Macrophage", "NK", "Other T", "Plasma")

rm(counts); gc()

# ── Helper functions ───────────────────────────────────────────────────────────
beta_shrink <- function(k, n, alpha = 1, beta = 1) (k + alpha) / (n + alpha + beta)

get_pct <- function(cells, genes) {
  if (length(cells) == 0) return(setNames(rep(0, length(genes)), genes))
  m    <- mat[genes, cells, drop = FALSE]
  rowMeans(m > 0)
}

get_mean <- function(cells, genes) {
  if (length(cells) == 0) return(setNames(rep(0, length(genes)), genes))
  rowMeans(mat[genes, cells, drop = FALSE])
}

# ── DEGs: Late-Treg vs Early-Treg (targets of ligand activity) ────────────────
treg_de <- FindMarkers(obj,
                       ident.1 = "Treg_Late",
                       ident.2 = "Treg_Early",
                       min.pct = 0.1,
                       logfc.threshold = 0.1,
                       test.use = "wilcox")
treg_de$gene <- rownames(treg_de)

late_targets  <- treg_de %>% filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>% pull(gene)
early_targets <- treg_de %>% filter(avg_log2FC < -0.25, p_val_adj < 0.05) %>% pull(gene)

# ── NicheNet ligand activity (AUPR): late targets ─────────────────────────────
lt_mat <- ligand_target_matrix[intersect(rownames(ligand_target_matrix), late_targets),
                                intersect(colnames(ligand_target_matrix), ligands_use),
                                drop = FALSE]
background_genes <- intersect(genes_avail, rownames(ligand_target_matrix))

aupr_late <- sapply(colnames(lt_mat), function(lig) {
  scores <- lt_mat[, lig]
  fg     <- scores[intersect(names(scores), late_targets)]
  bg     <- scores[setdiff(names(scores), late_targets)]
  if (length(fg) < 3 || length(bg) < 3) return(NA_real_)
  tryCatch(pr.curve(scores.class1 = fg, scores.class0 = bg)$auc.integral,
           error = function(e) NA_real_)
})

# ── NicheNet ligand destabilization (AUPR): early targets ─────────────────────
et_mat <- ligand_target_matrix[intersect(rownames(ligand_target_matrix), early_targets),
                                intersect(colnames(ligand_target_matrix), ligands_use),
                                drop = FALSE]

aupr_destab <- sapply(colnames(et_mat), function(lig) {
  scores <- et_mat[, lig]
  fg     <- scores[intersect(names(scores), early_targets)]
  bg     <- scores[setdiff(names(scores), early_targets)]
  if (length(fg) < 3 || length(bg) < 3) return(NA_real_)
  tryCatch(pr.curve(scores.class1 = fg, scores.class0 = bg)$auc.integral,
           error = function(e) NA_real_)
})

# ── Per-sender ligand availability (mean expression) ──────────────────────────
sender_lig_avail <- lapply(SENDER_TYPES, function(stype) {
  cells <- colnames(obj)[obj$ct_lr == stype]
  get_mean(cells, intersect(ligands_use, genes_avail))
})
names(sender_lig_avail) <- SENDER_TYPES

# ── Receptor prevalence: early Tregs ──────────────────────────────────────────
rec_early_pct <- get_pct(early_cells, intersect(receptors_use, genes_avail))

# ── Receptor late up-regulation ───────────────────────────────────────────────
rec_late_mean  <- get_mean(late_cells,  intersect(receptors_use, genes_avail))
rec_early_mean <- get_mean(early_cells, intersect(receptors_use, genes_avail))
rec_late_fc    <- log2((rec_late_mean + eps) / (rec_early_mean + eps))

# ── LR pair table ─────────────────────────────────────────────────────────────
lr_use <- lr_network %>%
  filter(from %in% names(aupr_late), to %in% names(rec_early_pct)) %>%
  rename(ligand = from, receptor = to) %>%
  select(ligand, receptor) %>%
  distinct()

# ── Priority score ────────────────────────────────────────────────────────────
alpha_w <- 2; beta_w <- 1.5; gamma_w <- 1; delta_w <- 1; eta_w <- 1.5

score_df <- lr_use %>%
  rowwise() %>%
  mutate(
    aupr_act    = aupr_late[ligand],
    aupr_dstb   = aupr_destab[ligand],
    lig_avail   = max(sapply(SENDER_TYPES, function(s) sender_lig_avail[[s]][ligand]), na.rm = TRUE),
    rec_prev    = rec_early_pct[receptor],
    rec_fc_late = rec_late_fc[receptor]
  ) %>%
  ungroup() %>%
  filter(!is.na(aupr_act)) %>%
  mutate(
    pct_act    = percent_rank(aupr_act),
    pct_dstb   = percent_rank(aupr_dstb),
    pct_avail  = percent_rank(lig_avail),
    pct_prev   = percent_rank(rec_prev),
    pct_fc     = percent_rank(rec_fc_late),
    priority   = (pct_act^alpha_w * pct_dstb^beta_w * pct_avail^gamma_w *
                    pct_prev^delta_w * pct_fc^eta_w)^
                   (1 / (alpha_w + beta_w + gamma_w + delta_w + eta_w))
  ) %>%
  arrange(desc(priority))

write.csv(score_df,
          file.path(OUTPUT_DIR, "nichenet_destab_late_treg_ligand_strength_imm_logpr.csv"),
          row.names = FALSE)

# ── Best sender per ligand ─────────────────────────────────────────────────────
best_sender_df <- lapply(unique(score_df$ligand), function(lig) {
  avail <- sapply(SENDER_TYPES, function(s) {
    v <- sender_lig_avail[[s]][lig]
    if (is.null(v) || is.na(v)) 0 else v
  })
  data.frame(ligand = lig, best_sender = names(which.max(avail)),
             sender_mean_expr = max(avail, na.rm = TRUE))
}) %>% bind_rows()

write.csv(best_sender_df,
          file.path(OUTPUT_DIR, "nichenet_best_sender_per_ligand.csv"),
          row.names = FALSE)

# ── Chord diagram: top LR pairs ───────────────────────────────────────────────
top_pairs <- score_df %>%
  left_join(best_sender_df, by = "ligand") %>%
  slice_max(priority, n = 30)

chord_df <- top_pairs %>%
  transmute(
    from  = best_sender,
    to    = "Treg",
    value = priority,
    label = paste0(ligand, "–", receptor)
  )

sender_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(SENDER_TYPES)),
  SENDER_TYPES
)

pdf(file.path(OUTPUT_DIR, "plots/nichenet_chord_diagram.pdf"), width = 8, height = 8)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(chord_df[, c("from", "to", "value")],
             grid.col  = c(sender_colors, Treg = "#d62728"),
             transparency = 0.4,
             annotationTrack = "grid",
             preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)
dev.off()

# ── Heatmap: top ligands × senders ────────────────────────────────────────────
top_ligs    <- top_pairs %>% pull(ligand) %>% unique() %>% head(20)
heat_mat    <- sapply(SENDER_TYPES, function(s) {
  sapply(top_ligs, function(lig) {
    v <- sender_lig_avail[[s]][lig]
    if (is.null(v) || is.na(v)) 0 else v
  })
})

pdf(file.path(OUTPUT_DIR, "plots/nichenet_heatmap.pdf"), width = 6, height = 6)
pheatmap::pheatmap(heat_mat,
                   color = colorRampPalette(c("white", "#1f77b4"))(50),
                   cluster_cols = FALSE,
                   fontsize_row = 9, fontsize_col = 9,
                   main = "Ligand availability by sender")
dev.off()

# ── Parameter sensitivity: Spearman rank correlation across weight schemes ─────
weight_schemes <- list(
  base      = c(alpha_w = 2,   beta_w = 1.5, gamma_w = 1, delta_w = 1, eta_w = 1.5),
  act_only  = c(alpha_w = 4,   beta_w = 0,   gamma_w = 1, delta_w = 1, eta_w = 0),
  balanced  = c(alpha_w = 1,   beta_w = 1,   gamma_w = 1, delta_w = 1, eta_w = 1),
  rec_focus = c(alpha_w = 1,   beta_w = 1,   gamma_w = 0.5, delta_w = 2, eta_w = 2),
  dstb_focus= c(alpha_w = 1,   beta_w = 3,   gamma_w = 1, delta_w = 1, eta_w = 1.5)
)

compute_priority <- function(df, weights) {
  total <- sum(weights)
  df %>%
    mutate(priority_new =
             (pct_act^weights["alpha_w"] * pct_dstb^weights["beta_w"] *
              pct_avail^weights["gamma_w"] * pct_prev^weights["delta_w"] *
              pct_fc^weights["eta_w"])^(1 / total)) %>%
    pull(priority_new)
}

rank_mat <- sapply(weight_schemes, function(w) {
  rank(-compute_priority(score_df, w))
})

cor_mat <- cor(rank_mat, method = "spearman")

pdf(file.path(OUTPUT_DIR, "plots/nichenet_weight_sensitivity.pdf"), width = 5, height = 5)
pheatmap::pheatmap(cor_mat,
                   color       = colorRampPalette(c("#f0f0f0", "#2166ac"))(50),
                   breaks      = seq(0.8, 1, length.out = 51),
                   display_numbers = TRUE,
                   number_format = "%.3f",
                   fontsize    = 10,
                   main        = "Spearman correlation of priority rankings\nacross weight schemes")
dev.off()

cat("NicheNet analysis complete.\n")
cat("Top 10 LR pairs by priority:\n")
print(score_df %>% left_join(best_sender_df, by = "ligand") %>%
        select(ligand, receptor, best_sender, priority) %>% head(10))
