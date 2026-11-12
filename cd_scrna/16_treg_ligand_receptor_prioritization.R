# CD scRNA-seq: NicheNet ligand-receptor prioritization of signals linked to loss of Treg suppression.
#
# Panels made
#   Fig 4H     ligand availability by function from immune, stromal and epithelial senders to Tregs
#   Fig 4I     top 25 sender-ligand-receptor interactions (composite priority; chord diagram)
#   Supp 8G-H  stability of the LR ranking under five priority weighting schemes
#   Section immune_sender_ranking ranks immune-sender ligands; cosmx_cd/10_treg_suppression_score.py takes the
#   antigen-presentation LR module from this ranking and leaves its genes out of the Treg Suppression Score.
#
# Inputs
#   /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/nichenet_inputs/ (15_treg_ligand_receptor_inputs.py)
#   /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/nichenet_outputs_lossact/ (two precomputed gene-set
#     tables read by F4HI_nichenet_lr, see the section note)
#   /share/fsmresfiles/ydn4687/nichenet/ (NicheNet v2 lr_network, ligand_target_matrix, weighted_networks)
#
# Outputs
#   /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/nichenet_outputs/ (immune_sender_ranking) and
#     /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/nichenet_outputs_general_destab/ (F4HI_nichenet_lr)
#   /share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/ (Fig 4H-I panels and CSVs)
#   /share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness/ (Supp 8G-H panels and CSVs)
#
# Run order of sections
#   15_treg_ligand_receptor_inputs.py -> immune_sender_ranking -> F4HI_nichenet_lr -> S8GH_lr_weighting_schemes
#   immune_sender_ranking must run before cosmx_cd/10_treg_suppression_score.py.
#
# How to run
#   Rscript <this file>                          run every section in the order below
#   RUN_SECTIONS=<id>[,<id>...] Rscript <file>   run only the listed sections (still in file order)
#   Each section runs one analysis step inside `if (.run("<id>")) { ... }` and reads its inputs from disk, so a single
#   section can be run once its upstream outputs exist.
#   Paths are the lab's absolute paths (to be replaced by placeholders before public release).
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {   # all sections run in file order unless RUN_SECTIONS=id1,id2,... is set
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}

# ---- helper (always loaded): LR chord diagram (Fig 4I) with the title centred over the whole figure and a legend panel.
draw_lr_chord <- function(chord_data, grid_cols, file, main) {
  ragg::agg_png(file, width = 11, height = 8.5, units = "in", res = 300)
  par(oma = c(0, 0, 3, 0))
  layout(matrix(1:2, nrow = 1), widths = c(0.68, 0.32))
  par(mar = c(1, 1, 1, 1))
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.01, 0.01), cell.padding = c(0.02, 0, 0.02, 0))
  chordDiagram(chord_data[, c("from", "to", "priority_rank_geom")], grid.col = grid_cols,
               transparency = chord_data$transparency, directional = 1, direction.type = "arrows",
               link.arr.type = "big.arrow", link.arr.length = 0.1, annotationTrack = "grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(chord_data)))), track.margin = c(0.01, 0)))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), cex = 0.95, font = 2)
  }, bg.border = NA)
  circos.clear()
  mtext(main, outer = TRUE, side = 3, line = 0.5, cex = 1.6, font = 2, adj = 0.5)
  par(mar = c(1, 0, 1, 0.5)); plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  y <- 0.90; text(0.02, y, "Node colour: total priority", font = 2, cex = 1.0, adj = 0)
  for (k in 1:2) {
    y <- y - 0.08
    pal <- if (k == 1) c("#fee5d9", "#fc9272", "#a50f15") else c("#deebf7", "#9ecae1", "#3182bd")
    for (j in 1:3) rect(0.30 + 0.12 * (j - 1), y - 0.02, 0.42 + 0.12 * (j - 1), y + 0.02, col = pal[j], border = "grey30")
    text(0.02, y, if (k == 1) "Ligands" else "Receptors", cex = 0.9, adj = 0)
    text(0.30, y - 0.04, "Low", cex = 0.75, adj = 0); text(0.66, y - 0.04, "High", cex = 0.75, adj = 1)
    text(0.68, y, if (k == 1) "(immune sender)" else "(Tregs)", cex = 0.8, adj = 0, font = 3)
  }
  y <- y - 0.13; text(0.02, y, "Link width: LR priority", font = 2, cex = 1.0, adj = 0)
  for (lw in c(1, 4, 8)) { y <- y - 0.045; segments(0.02, y, 0.36, y, lwd = lw, col = "grey50")
    text(0.40, y, c("1" = "Low", "4" = "Medium", "8" = "High")[as.character(lw)], cex = 0.85, adj = 0) }
  y <- y - 0.11; text(0.02, y, "Link opacity: receptor log2FC\n(late vs early Treg)", font = 2, cex = 1.0, adj = 0); y <- y - 0.03
  for (a in c(0.3, 0.6, 0.9)) { y <- y - 0.045; rect(0.02, y - 0.018, 0.36, y + 0.018, col = rgb(0.5, 0.5, 0.5, a), border = "grey30")
    text(0.40, y, c("Low", "Medium", "High")[match(a, c(0.3, 0.6, 0.9))], cex = 0.85, adj = 0) }
  y <- y - 0.11; text(0.02, y, "Direction", font = 2, cex = 1.0, adj = 0)
  y <- y - 0.05; arrows(0.02, y, 0.30, y, length = 0.12, lwd = 2, col = "grey30"); text(0.34, y, "Ligand → Receptor", cex = 0.85, adj = 0)
  dev.off()
}


# ======================================================================================================================
# Section immune_sender_ranking: NicheNet immune -> Treg LR prioritization (immune senders).
# Panels: Treg Suppression Score prerequisite (Fig 5, Supp 9, Fig 6A): its ranking (objects/nichenet_outputs/) defines the
#   antigen-presentation Treg LR module of cosmx_cd/10_treg_suppression_score.py (section antigen_presentation_module); the
#   genes of that module are excluded from the Treg Suppression Score signature.
# Note: draw_lr_chord() is defined in the helper block at the top of this file
# ======================================================================================================================
if (.run("immune_sender_ranking")) {
# NicheNet immune -> Treg LR prioritization. Senders: B, DC, Granulocyte, Macrophage, NK, Other T, Plasma
#   (objects/nichenet_inputs/strict_treg_imm_*_rowfix.parquet). The state-transition gene set is the hand-picked NF-kB / inflammatory
#   destabilisation gene set defined in step 4.
# Receivers: pseudotime Tregs, Early = bins 0-1, Late = bins 2-4. Steps: (1) receiver DE (late vs early) and receptor
#   detection; (2) sender ligand availability; (3) NicheNet ligand activity on late-Treg genes per sender; (4) ligand
#   activity on the Treg state-transition gene set; (5) LR pairs per sender; (6) composite priority = weighted geometric
#   mean of the percentile ranks r_destab, r_late, r_avail, r_recEarly, r_recLateUp (alpha = 2, beta = 1.5, gamma = 1,
#   delta = 1, eta = 1.5).
# The dense parquet is converted to a sparse matrix in column chunks (avoids a ~100 GB dense copy); ScaleData() is not
#   run because only the 'data' layer is used (FindMarkers / get_expressed_genes / predict_ligand_activities).
# Outputs (objects/nichenet_outputs/): priority table nichenet_destab_late_treg_ligand_strength_imm_logpr.csv, its best
#   sender per LR pair (_best_sender.csv, read by cosmx_cd/10_treg_suppression_score.py), component tables, RDS.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble)
  library(Matrix); library(purrr); library(PRROC); library(scales)
  library(nichenetr); library(stringr); library(forcats); library(ggplot2)
  library(arrow); library(circlize); library(grid)
})
set.seed(1)
eps <- 1e-6
BASE <- "/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix"
OUT  <- file.path(BASE, "objects", "nichenet_outputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
meta_f   <- file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_meta_rowfix.parquet")
counts_f <- file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_counts_rowfix.parquet")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# read in data
meta <- as.data.frame(read_parquet(meta_f))
tbl  <- read_parquet(counts_f, as_data_frame = FALSE)
bad_gene <- "__index_level_0__"
genes <- setdiff(names(tbl), bad_gene)
if (bad_gene %in% colnames(meta)) meta <- meta[, colnames(meta) != bad_gene, drop = FALSE]
msg("parquet", nrow(tbl), "cells x", length(genes), "genes")
chunks <- split(genes, ceiling(seq_along(genes) / 1000))
counts_sp <- do.call(cbind, lapply(chunks, function(g) {
  m <- as.matrix(as.data.frame(tbl[, g]))
  Matrix(m, sparse = TRUE)
}))
rownames(counts_sp) <- meta$index
rm(tbl); invisible(gc())
msg("sparse counts built", dim(counts_sp))

obj <- CreateSeuratObject(counts = t(counts_sp))
obj$ct_lr <- meta$ct_lr
rm(counts_sp); invisible(gc())

# Map Treg_B* -> Early/Late; keep non-Treg types unchanged
treg_mask <- grepl("^Treg_B[0-4]$", obj$ct_lr)
pt_bin <- rep(NA_character_, ncol(obj))
pt_bin[treg_mask] <- sub("^Treg_(B[0-4])$", "\\1", obj$ct_lr[treg_mask])
EARLY_BINS <- c("B0","B1"); LATE_BINS <- c("B2","B3","B4")
nn_id <- obj$ct_lr
nn_id[treg_mask & pt_bin %in% EARLY_BINS] <- "Treg_Early"
nn_id[treg_mask & pt_bin %in% LATE_BINS]  <- "Treg_Late"
obj$nn_id <- factor(nn_id)
Idents(obj) <- obj$nn_id
idt <- as.data.frame(table(Idents(obj))); colnames(idt) <- c("identity", "n_cells")
write.csv(idt, file.path(OUT, "identity_counts.csv"), row.names = FALSE)
print(idt)

obj <- NormalizeData(obj, layer = "counts", new.layer = "data")
mat <- GetAssayData(obj, layer = "data")  # log-normalised expression
genes_avail <- rownames(mat)
early_cells <- colnames(obj)[obj$nn_id == "Treg_Early"]
late_cells  <- colnames(obj)[obj$nn_id == "Treg_Late"]

NN <- "/share/fsmresfiles/ydn4687/nichenet"
lr_network <- readRDS(file.path(NN, "lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file.path(NN, "ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(file.path(NN, "weighted_networks_nsga2r_final.rds"))
ligands_use   <- intersect(unique(lr_network$from), genes_avail)
receptors_use <- intersect(unique(lr_network$to),   genes_avail)
TREG_LABELS   <- c("Treg_Early","Treg_Late")
OTHER_SENDERS <- c('B', 'DC','Granulocyte','Macrophage', 'Neuron','NK','Other T' ,'Plasma')

beta_shrink <- function(k, n, alpha = 1, beta = 1) (k + alpha) / (n + alpha + beta)
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

## 1) Receiver (Tregs): DE & detection
seurat_obj_receiver <- subset(obj, idents = TREG_LABELS)
Idents(seurat_obj_receiver) <- factor(Idents(seurat_obj_receiver), levels = c("Treg_Early","Treg_Late"))
de_rec <- FindMarkers(seurat_obj_receiver, ident.1 = "Treg_Late", ident.2 = "Treg_Early",
                      logfc.threshold = 0, min.pct = 0.05, only.pos = FALSE) %>% tibble::rownames_to_column("gene")
write.csv(de_rec, file.path(OUT, "treg_late_vs_early_DE.csv"), row.names = FALSE)
rec_bin_early <- mat[receptors_use, early_cells, drop = FALSE] > 0
rec_bin_late  <- mat[receptors_use, late_cells,  drop = FALSE] > 0
k_early <- Matrix::rowSums(rec_bin_early); k_late <- Matrix::rowSums(rec_bin_late)
n_early <- length(early_cells); n_late <- length(late_cells)
rec_frac_early_shrunk <- beta_shrink(k_early, n_early, alpha = 1, beta = 1)
names(rec_frac_early_shrunk) <- rownames(rec_bin_early)
rec_tbl <- de_rec %>%
  dplyr::filter(gene %in% receptors_use) %>%
  dplyr::mutate(log2fc_late = avg_log2FC, sig_late = p_val_adj < 0.1 & log2fc_late >= 0.25) %>%
  dplyr::left_join(tibble::tibble(gene = names(rec_frac_early_shrunk),
                                  rec_frac_early = as.numeric(rec_frac_early_shrunk)), by = "gene")
# raw % of early / late Tregs expressing each receptor (reported in the table, not used in scoring)
rec_pct <- tibble::tibble(gene = names(k_early), pct_early = 100 * k_early / n_early, pct_late = 100 * k_late / n_late)
write.csv(rec_tbl %>% left_join(rec_pct, by = "gene"), file.path(OUT, "receiver_receptor_table.csv"), row.names = FALSE)

## 2) Sender ligand availability
sender_expr <- dplyr::bind_rows(lapply(OTHER_SENDERS, function(s) {
  cells <- colnames(obj)[Idents(obj) == s]
  m <- mat[ligands_use, cells, drop = FALSE]
  k <- Matrix::rowSums(m > 0); n <- length(cells)
  pi_hat <- beta_shrink(k, n, alpha = 1, beta = 1)
  pos_means <- vapply(seq_len(nrow(m)), function(i) { v <- m[i, ]; vp <- v[v > 0]; if (length(vp) == 0) 0 else mean(vp) }, numeric(1))
  names(pos_means) <- rownames(m)
  mu_scaled <- scale01(log1p(pos_means))
  tibble::tibble(sender = s, ligand = rownames(m), L_k = as.numeric(k), L_n = n, L_pi_hat = as.numeric(pi_hat),
                 L_mu_pos = as.numeric(pos_means), L_mu_pos01 = as.numeric(mu_scaled), Avail_L = as.numeric(pi_hat * mu_scaled))
}))
write.csv(sender_expr, file.path(OUT, "sender_ligand_availability.csv"), row.names = FALSE)
msg("sender availability done")

## 3) NicheNet ligand activity (Late-focused)
run_nn_sender <- function(sender_label) {
  nichenet_seuratobj_cluster_de(
    seurat_obj = obj, receiver_reference = "Treg_Early", receiver_affected = "Treg_Late",
    sender = sender_label, ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
    weighted_networks = weighted_networks, expression_pct = 0.05, lfc_cutoff = 0.25,
    geneset = "up", cutoff_visualization = 0.33)
}
nn_results <- setNames(lapply(OTHER_SENDERS, function(s) { msg("nichenet sender", s); run_nn_sender(s) }), OTHER_SENDERS)
saveRDS(lapply(nn_results, function(x) list(ligand_activities = x$ligand_activities, top_ligands = x$top_ligands,
                                            geneset_oi = x$geneset_oi, background_expressed_genes = x$background_expressed_genes)),
        file.path(OUT, "nn_results_ligand_activities.rds"))
lig_act_by_sender <- lapply(names(nn_results), function(s) {
  nn_results[[s]]$ligand_activities %>%
    dplyr::transmute(from = .data$test_ligand, A_late = pmin(pmax(.data$aupr_corrected, 0), 1)) %>%
    dplyr::mutate(sender = s)
}) |> dplyr::bind_rows()
write.csv(lig_act_by_sender, file.path(OUT, "ligand_late_activity_by_sender.csv"), row.names = FALSE)
cat("late-Treg up geneset size per sender:", sapply(nn_results, function(x) length(x$geneset_oi)), "\n")
writeLines(nn_results[[1]]$geneset_oi, file.path(OUT, "late_treg_up_geneset.txt"))

## 4) Destabilization activity
receiver_cells <- c(early_cells, late_cells)
expr_frac_receiver <- Matrix::rowMeans(mat[, receiver_cells, drop = FALSE] > 0)
bg_expressed_genes <- names(expr_frac_receiver)[expr_frac_receiver >= 0.01]
destab_prog_use <- c("NFKBIA","NFKBIZ","TNFAIP3","RELB","TNF","IFNG","OSM","CXCL8",
                     "CCL3","CCL4","CCL5","CCL20","IL23A","TNFSF9","TNFSF14","TNFSF12","TNFSF13B")
geneset_destab <- intersect(destab_prog_use, bg_expressed_genes)
cat("destab geneset in background:", geneset_destab, "\n")
all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, bg_expressed_genes)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
lig_destab_activity_raw <- nichenetr::predict_ligand_activities(
  geneset = geneset_destab, background_expressed_genes = bg_expressed_genes,
  ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
lig_destab_activity <- tibble::tibble(ligand = potential_ligands) %>%
  dplyr::left_join(lig_destab_activity_raw %>% dplyr::select(ligand = test_ligand, A_destab = aupr_corrected), by = "ligand")
write.csv(lig_destab_activity, file.path(OUT, "ligand_destab_activity.csv"), row.names = FALSE)

## 5) LR pairs per sender
lr_by_sender <- tidyr::crossing(sender = OTHER_SENDERS, from = ligands_use) %>%
  dplyr::left_join(lr_network %>% dplyr::distinct(from, to) %>% dplyr::filter(from %in% ligands_use, to %in% receptors_use), by = "from") %>%
  dplyr::left_join(lig_act_by_sender, by = c("sender","from")) %>%
  dplyr::mutate(A_late = dplyr::coalesce(A_late, eps)) %>%
  dplyr::left_join(lig_destab_activity, by = c("from" = "ligand")) %>%
  dplyr::mutate(A_destab = dplyr::coalesce(A_destab, eps))

## 6) Priority
pctrank <- function(x) {
  n <- sum(!is.na(x)); if (n <= 1) return(rep(0.5, length(x)))
  rx <- rank(x, na.last = "keep", ties.method = "average"); as.numeric((rx - 1) / pmax(n - 1, 1))
}
alpha <- 2; beta <- 1.5; gamma <- 1; delta <- 1; eta <- 1.5
OTHER_SENDERS <- setdiff(OTHER_SENDERS, "Neuron")
lr_by_sender <- lr_by_sender %>% filter(sender != 'Neuron')
df <- lr_by_sender %>%
  dplyr::left_join(sender_expr %>% dplyr::select(sender, ligand, Avail_L, L_pi_hat, L_mu_pos01), by = c("sender","from" = "ligand")) %>%
  dplyr::left_join(rec_tbl %>% dplyr::select(to = gene, log2fc_late, rec_frac_early), by = "to") %>%
  filter(!is.na(rec_frac_early)) %>%
  mutate(eps = 1e-6, lfc_cap = 2, log2fc_cap = pmin(pmax(log2fc_late, -lfc_cap), lfc_cap),
         RecLateUp_R = stats::plogis(1.0 * log2fc_cap)) %>%
  filter(log2fc_late > 0) %>%
  mutate(A_late = dplyr::coalesce(A_late, eps), A_destab = dplyr::coalesce(A_destab, eps),
         A_late = pmin(pmax(A_late, eps), 1), A_destab = pmin(pmax(A_destab, eps), 1),
         Avail_L = pmin(pmax(Avail_L, eps), 1), RecEarly_R = pmin(pmax(rec_frac_early, eps), 1),
         RecLateUp_R = pmin(pmax(RecLateUp_R, eps), 1),
         r_destab = pctrank(A_destab), r_late = pctrank(A_late), r_avail = pctrank(Avail_L),
         r_recEarly = pctrank(RecEarly_R), r_recLateUp = pctrank(RecLateUp_R),
         rank_geom_num = (r_destab ^ alpha) * (r_late ^ beta) * (r_avail ^ gamma) * (r_recEarly ^ delta) * (r_recLateUp ^ eta),
         rank_geom_den = alpha + beta + gamma + delta + eta,
         priority_rank_geom = rank_geom_num ^ (1 / rank_geom_den)) %>%
  arrange(desc(priority_rank_geom)) %>%
  select(sender, from, to, A_late, A_destab, L_pi_hat, L_mu_pos01, Avail_L, rec_frac_early, log2fc_late, RecLateUp_R,
         r_destab, r_late, r_avail, r_recEarly, r_recLateUp, priority_rank_geom)
write.csv(df, file.path(OUT, "nichenet_destab_late_treg_ligand_strength_imm_logpr.csv"))
best_sender_per_LR <- df %>% group_by(from, to) %>% slice_max(order_by = priority_rank_geom, n = 1, with_ties = FALSE) %>% ungroup()
write.csv(best_sender_per_LR, file.path(OUT, "nichenet_destab_late_treg_ligand_strength_imm_logpr_best_sender.csv"))
msg("priority table:", nrow(df), "rows;", nrow(best_sender_per_LR), "LR pairs")

saveRDS(list(df = df, best_sender_per_LR = best_sender_per_LR, sender_expr = sender_expr, rec_tbl = rec_tbl,
             lig_act_by_sender = lig_act_by_sender, lig_destab_activity = lig_destab_activity,
             geneset_destab = geneset_destab, identity_counts = idt), file.path(OUT, "nichenet_outputs.rds"))
msg("done")
}  # end section immune_sender_ranking


# ======================================================================================================================
# Section F4HI_nichenet_lr: NicheNet LR prioritization with immune, stromal and epithelial senders (Fig 4H-I).
# Panels: Fig 4H, 4I
# Note: draw_lr_chord() is defined in the helper block at the top of this file
# Note: reads two precomputed tables in objects/nichenet_outputs_lossact/ (loss_geneset_lineage_check.csv: per-gene ratio of
#   Treg expression to the highest other lineage, used to drop lineage contaminants; destab_markers_bin4_vs_bin0.csv:
#   bin 4 vs bin 0 destabilisation markers). The script that writes them is not part of this repository.
# ======================================================================================================================
if (.run("F4HI_nichenet_lr")) {
# NicheNet LR prioritization for Fig 4H-I. Senders = general cell types: the immune groups of immune_sender_ranking plus
#   Fibroblast, Endothelial and Epithelial (objects/nichenet_inputs/strict_treg_imm_general_*_rowfix.parquet,
#   15_treg_ligand_receptor_inputs.py section F4HI_nichenet_inputs).
# State-transition gene set (step 4) = data-derived loss-of-activation program (genes down in late vs early Tregs:
#   FindMarkers p_val_adj < 0.05, avg_log2FC <= -0.25, expressed in >= 1% of Tregs; up to 250 by most negative log2FC)
#   minus lineage contaminants (Treg / max other lineage expression ratio < 0.2) plus destabilisation markers up in
#   bin 4 vs bin 0 (Wilcoxon p < 0.05).
# Receivers: pseudotime Tregs, Early = bins 0-1, Late = bins 2-4. Steps: (1) receiver DE (late vs early) and receptor
#   detection; (2) sender ligand availability; (3) NicheNet ligand activity on late-Treg genes per sender; (4) ligand
#   activity on the Treg state-transition gene set; (5) LR pairs per sender; (6) composite priority = weighted geometric
#   mean of the percentile ranks r_destab, r_late, r_avail, r_recEarly, r_recLateUp (alpha = 2, beta = 1.5, gamma = 1,
#   delta = 1, eta = 1.5).
# The dense parquet is converted to a sparse matrix in column chunks (avoids a ~100 GB dense copy); ScaleData() is not
#   run because only the 'data' layer is used (FindMarkers / get_expressed_genes / predict_ligand_activities).
# Outputs: objects/nichenet_outputs_general_destab/ (priority tables, RDS); Fig3H_sender_ligand_families.png/.csv (Fig 4H)
#   and Fig3I_LR_chord.png, Fig3I_top25_LR.csv (Fig 4I) in the Fig3D-I folder.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble)
  library(Matrix); library(purrr); library(PRROC); library(scales)
  library(nichenetr); library(stringr); library(forcats); library(ggplot2)
  library(arrow); library(circlize); library(grid)
})
set.seed(1)
eps <- 1e-6
BASE <- "/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix"
OUT  <- file.path(BASE, "objects", "nichenet_outputs_general_destab")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG  <- "/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals"   # Fig 4H-I panels and CSVs
meta_f   <- file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_general_meta_rowfix.parquet")
counts_f <- file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_general_counts_rowfix.parquet")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# read in data
meta <- as.data.frame(read_parquet(meta_f))
tbl  <- read_parquet(counts_f, as_data_frame = FALSE)
bad_gene <- "__index_level_0__"
genes <- setdiff(names(tbl), bad_gene)
if (bad_gene %in% colnames(meta)) meta <- meta[, colnames(meta) != bad_gene, drop = FALSE]
msg("parquet", nrow(tbl), "cells x", length(genes), "genes")
chunks <- split(genes, ceiling(seq_along(genes) / 1000))
counts_sp <- do.call(cbind, lapply(chunks, function(g) {
  m <- as.matrix(as.data.frame(tbl[, g]))
  Matrix(m, sparse = TRUE)
}))
rownames(counts_sp) <- meta$index
rm(tbl); invisible(gc())
msg("sparse counts built", dim(counts_sp))

obj <- CreateSeuratObject(counts = t(counts_sp))
obj$ct_lr <- meta$ct_lr
print(table(obj$ct_lr))
rm(counts_sp); invisible(gc())

# Map Treg_B* -> Early/Late; keep non-Treg types unchanged
treg_mask <- grepl("^Treg_B[0-4]$", obj$ct_lr)
pt_bin <- rep(NA_character_, ncol(obj))
pt_bin[treg_mask] <- sub("^Treg_(B[0-4])$", "\\1", obj$ct_lr[treg_mask])
EARLY_BINS <- c("B0","B1"); LATE_BINS <- c("B2","B3","B4")
nn_id <- obj$ct_lr
nn_id[treg_mask & pt_bin %in% EARLY_BINS] <- "Treg_Early"
nn_id[treg_mask & pt_bin %in% LATE_BINS]  <- "Treg_Late"
obj$nn_id <- factor(nn_id)
Idents(obj) <- obj$nn_id
idt <- as.data.frame(table(Idents(obj))); colnames(idt) <- c("identity", "n_cells")
write.csv(idt, file.path(OUT, "identity_counts.csv"), row.names = FALSE)
print(idt)

obj <- NormalizeData(obj, layer = "counts", new.layer = "data")
mat <- GetAssayData(obj, layer = "data")  # log-normalised expression
genes_avail <- rownames(mat)
early_cells <- colnames(obj)[obj$nn_id == "Treg_Early"]
late_cells  <- colnames(obj)[obj$nn_id == "Treg_Late"]

NN <- "/share/fsmresfiles/ydn4687/nichenet"
lr_network <- readRDS(file.path(NN, "lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file.path(NN, "ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(file.path(NN, "weighted_networks_nsga2r_final.rds"))
ligands_use   <- intersect(unique(lr_network$from), genes_avail)
receptors_use <- intersect(unique(lr_network$to),   genes_avail)
TREG_LABELS   <- c("Treg_Early","Treg_Late")
STROMAL_SENDERS <- c('Fibroblast', 'Endothelial', 'Epithelial')   # general non-immune senders
OTHER_SENDERS <- c('B', 'DC','Granulocyte','Macrophage', 'Neuron','NK','Other T' ,'Plasma', STROMAL_SENDERS)

beta_shrink <- function(k, n, alpha = 1, beta = 1) (k + alpha) / (n + alpha + beta)
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

## 1) Receiver (Tregs): DE & detection
seurat_obj_receiver <- subset(obj, idents = TREG_LABELS)
Idents(seurat_obj_receiver) <- factor(Idents(seurat_obj_receiver), levels = c("Treg_Early","Treg_Late"))
de_rec <- FindMarkers(seurat_obj_receiver, ident.1 = "Treg_Late", ident.2 = "Treg_Early",
                      logfc.threshold = 0, min.pct = 0.05, only.pos = FALSE) %>% tibble::rownames_to_column("gene")
write.csv(de_rec, file.path(OUT, "treg_late_vs_early_DE.csv"), row.names = FALSE)
rec_bin_early <- mat[receptors_use, early_cells, drop = FALSE] > 0
rec_bin_late  <- mat[receptors_use, late_cells,  drop = FALSE] > 0
k_early <- Matrix::rowSums(rec_bin_early); k_late <- Matrix::rowSums(rec_bin_late)
n_early <- length(early_cells); n_late <- length(late_cells)
rec_frac_early_shrunk <- beta_shrink(k_early, n_early, alpha = 1, beta = 1)
names(rec_frac_early_shrunk) <- rownames(rec_bin_early)
rec_tbl <- de_rec %>%
  dplyr::filter(gene %in% receptors_use) %>%
  dplyr::mutate(log2fc_late = avg_log2FC, sig_late = p_val_adj < 0.1 & log2fc_late >= 0.25) %>%
  dplyr::left_join(tibble::tibble(gene = names(rec_frac_early_shrunk),
                                  rec_frac_early = as.numeric(rec_frac_early_shrunk)), by = "gene")
# raw % of early / late Tregs expressing each receptor (reported in the table, not used in scoring)
rec_pct <- tibble::tibble(gene = names(k_early), pct_early = 100 * k_early / n_early, pct_late = 100 * k_late / n_late)
write.csv(rec_tbl %>% left_join(rec_pct, by = "gene"), file.path(OUT, "receiver_receptor_table.csv"), row.names = FALSE)

## 2) Sender ligand availability
sender_expr <- dplyr::bind_rows(lapply(OTHER_SENDERS, function(s) {
  cells <- colnames(obj)[Idents(obj) == s]
  m <- mat[ligands_use, cells, drop = FALSE]
  k <- Matrix::rowSums(m > 0); n <- length(cells)
  pi_hat <- beta_shrink(k, n, alpha = 1, beta = 1)
  pos_means <- vapply(seq_len(nrow(m)), function(i) { v <- m[i, ]; vp <- v[v > 0]; if (length(vp) == 0) 0 else mean(vp) }, numeric(1))
  names(pos_means) <- rownames(m)
  mu_scaled <- scale01(log1p(pos_means))
  tibble::tibble(sender = s, ligand = rownames(m), L_k = as.numeric(k), L_n = n, L_pi_hat = as.numeric(pi_hat),
                 L_mu_pos = as.numeric(pos_means), L_mu_pos01 = as.numeric(mu_scaled), Avail_L = as.numeric(pi_hat * mu_scaled))
}))
write.csv(sender_expr, file.path(OUT, "sender_ligand_availability.csv"), row.names = FALSE)
msg("sender availability done")

## 3) NicheNet ligand activity (Late-focused)
run_nn_sender <- function(sender_label) {
  nichenet_seuratobj_cluster_de(
    seurat_obj = obj, receiver_reference = "Treg_Early", receiver_affected = "Treg_Late",
    sender = sender_label, ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
    weighted_networks = weighted_networks, expression_pct = 0.05, lfc_cutoff = 0.25,
    geneset = "up", cutoff_visualization = 0.33)
}
nn_results <- setNames(lapply(OTHER_SENDERS, function(s) { msg("nichenet sender", s); run_nn_sender(s) }), OTHER_SENDERS)
saveRDS(lapply(nn_results, function(x) list(ligand_activities = x$ligand_activities, top_ligands = x$top_ligands,
                                            geneset_oi = x$geneset_oi, background_expressed_genes = x$background_expressed_genes)),
        file.path(OUT, "nn_results_ligand_activities.rds"))
lig_act_by_sender <- lapply(names(nn_results), function(s) {
  nn_results[[s]]$ligand_activities %>%
    dplyr::transmute(from = .data$test_ligand, A_late = pmin(pmax(.data$aupr_corrected, 0), 1)) %>%
    dplyr::mutate(sender = s)
}) |> dplyr::bind_rows()
write.csv(lig_act_by_sender, file.path(OUT, "ligand_late_activity_by_sender.csv"), row.names = FALSE)
cat("late-Treg up geneset size per sender:", sapply(nn_results, function(x) length(x$geneset_oi)), "\n")
writeLines(nn_results[[1]]$geneset_oi, file.path(OUT, "late_treg_up_geneset.txt"))

## 4) Destabilization activity
receiver_cells <- c(early_cells, late_cells)
expr_frac_receiver <- Matrix::rowMeans(mat[, receiver_cells, drop = FALSE] > 0)
bg_expressed_genes <- names(expr_frac_receiver)[expr_frac_receiver >= 0.01]
geneset_loss <- de_rec %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC <= -0.25, gene %in% bg_expressed_genes) %>%
  dplyr::arrange(avg_log2FC) %>% dplyr::slice_head(n = 250) %>% dplyr::pull(gene)   # loss-of-activation genes
lin <- read.csv(file.path(BASE, "objects/nichenet_outputs_lossact/loss_geneset_lineage_check.csv"), row.names = 1)
contam <- rownames(lin)[lin$ratio_treg_over_max_other < 0.2]
dm <- read.csv(file.path(BASE, "objects/nichenet_outputs_lossact/destab_markers_bin4_vs_bin0.csv"))
destab_up <- intersect(dm$gene[dm$wilcox_p < 0.05 & dm$median_delta > 0], bg_expressed_genes)
geneset_destab <- union(setdiff(geneset_loss, contam), destab_up)
writeLines(geneset_destab, file.path(OUT, "state_transition_geneset.txt"))
cat("loss genes", length(geneset_loss), "| contaminants removed", length(intersect(geneset_loss, contam)), "| destab markers added", paste(destab_up, collapse = ","), "| total", length(geneset_destab), "\n")
cat("destab geneset in background:", geneset_destab, "\n")
all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, bg_expressed_genes)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
lig_destab_activity_raw <- nichenetr::predict_ligand_activities(
  geneset = geneset_destab, background_expressed_genes = bg_expressed_genes,
  ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
lig_destab_activity <- tibble::tibble(ligand = potential_ligands) %>%
  dplyr::left_join(lig_destab_activity_raw %>% dplyr::select(ligand = test_ligand, A_destab = aupr_corrected), by = "ligand")
write.csv(lig_destab_activity, file.path(OUT, "ligand_destab_activity.csv"), row.names = FALSE)

## 5) LR pairs per sender
lr_by_sender <- tidyr::crossing(sender = OTHER_SENDERS, from = ligands_use) %>%
  dplyr::left_join(lr_network %>% dplyr::distinct(from, to) %>% dplyr::filter(from %in% ligands_use, to %in% receptors_use), by = "from") %>%
  dplyr::left_join(lig_act_by_sender, by = c("sender","from")) %>%
  dplyr::mutate(A_late = dplyr::coalesce(A_late, eps)) %>%
  dplyr::left_join(lig_destab_activity, by = c("from" = "ligand")) %>%
  dplyr::mutate(A_destab = dplyr::coalesce(A_destab, eps))

## 6) Priority
pctrank <- function(x) {
  n <- sum(!is.na(x)); if (n <= 1) return(rep(0.5, length(x)))
  rx <- rank(x, na.last = "keep", ties.method = "average"); as.numeric((rx - 1) / pmax(n - 1, 1))
}
alpha <- 2; beta <- 1.5; gamma <- 1; delta <- 1; eta <- 1.5
OTHER_SENDERS <- setdiff(OTHER_SENDERS, "Neuron")
lr_by_sender <- lr_by_sender %>% filter(sender != 'Neuron')
df <- lr_by_sender %>%
  dplyr::left_join(sender_expr %>% dplyr::select(sender, ligand, Avail_L, L_pi_hat, L_mu_pos01), by = c("sender","from" = "ligand")) %>%
  dplyr::left_join(rec_tbl %>% dplyr::select(to = gene, log2fc_late, rec_frac_early), by = "to") %>%
  filter(!is.na(rec_frac_early)) %>%
  mutate(eps = 1e-6, lfc_cap = 2, log2fc_cap = pmin(pmax(log2fc_late, -lfc_cap), lfc_cap),
         RecLateUp_R = stats::plogis(1.0 * log2fc_cap)) %>%
  filter(log2fc_late > 0) %>%
  mutate(A_late = dplyr::coalesce(A_late, eps), A_destab = dplyr::coalesce(A_destab, eps),
         A_late = pmin(pmax(A_late, eps), 1), A_destab = pmin(pmax(A_destab, eps), 1),
         Avail_L = pmin(pmax(Avail_L, eps), 1), RecEarly_R = pmin(pmax(rec_frac_early, eps), 1),
         RecLateUp_R = pmin(pmax(RecLateUp_R, eps), 1),
         r_destab = pctrank(A_destab), r_late = pctrank(A_late), r_avail = pctrank(Avail_L),
         r_recEarly = pctrank(RecEarly_R), r_recLateUp = pctrank(RecLateUp_R),
         rank_geom_num = (r_destab ^ alpha) * (r_late ^ beta) * (r_avail ^ gamma) * (r_recEarly ^ delta) * (r_recLateUp ^ eta),
         rank_geom_den = alpha + beta + gamma + delta + eta,
         priority_rank_geom = rank_geom_num ^ (1 / rank_geom_den)) %>%
  arrange(desc(priority_rank_geom)) %>%
  select(sender, from, to, A_late, A_destab, L_pi_hat, L_mu_pos01, Avail_L, rec_frac_early, log2fc_late, RecLateUp_R,
         r_destab, r_late, r_avail, r_recEarly, r_recLateUp, priority_rank_geom)
write.csv(df, file.path(OUT, "nichenet_destab_late_treg_ligand_strength_imm_logpr.csv"))
best_sender_per_LR <- df %>% group_by(from, to) %>% slice_max(order_by = priority_rank_geom, n = 1, with_ties = FALSE) %>% ungroup()
write.csv(best_sender_per_LR, file.path(OUT, "nichenet_destab_late_treg_ligand_strength_imm_logpr_best_sender.csv"))
msg("priority table:", nrow(df), "rows;", nrow(best_sender_per_LR), "LR pairs")

## ---------------- Fig 4H: ligand availability by function per sender ----------------
family_of <- function(g) {
  dplyr::case_when(
    grepl("^HLA-D", g) ~ "Antigen presentation",
    grepl("^(CCL|CXCL|CX3CL)", g) ~ "Chemokines",
    grepl("^(TNF|TNFSF|LTB)$", g) ~ "TNF/TNFSF",
    grepl("^COL[0-9A-Z]", g) | g %in% c("FN1","VTN","VCAN","COMP","TNC","TNXB","LAMA1","LAMA2","LAMA3","LAMA4","LAMA5",
      "LAMB1","LAMB2","LAMB3","LAMC1","LAMC2","LAMC3","FBLN1","FBN1","RELN","OGN","MATN1","CHAD","NID1","MIA","TNN","CCN1") ~ "ECM / matrix",
    grepl("^TGFB", g) | grepl("^ANGPT", g) | grepl("^ANGPTL", g) | grepl("^BMP", g) | grepl("^WNT", g) |
      g %in% c("MDK","TSLP","VSTM1") ~ "Growth / angiogenic",
    g %in% c("HMGB1","S100A8","MIF","IFNG","IFNE","GZMB","CRISP3") ~ "Stress / cytotoxic",
    grepl("^LGALS", g) | grepl("^CLEC", g) | grepl("^SIGLEC", g) | g %in% c("ZG16B") ~ "Lectin",
    grepl("^ADAM", g) | grepl("^MMP", g) | g %in% c("PLA2G2A","FURIN","FAP","PCSK9","REN","TGM2") ~ "Protease / sheddase",
    grepl("^IL[0-9A-Z]", g) | grepl("^IFN", g) ~ "Cytokines",
    TRUE ~ "Others")
}
cols <- c("Antigen presentation" = "#377EB8", "Chemokines" = "#556B2F", "Cytokines" = "#984EA3", "TNF/TNFSF" = "#F781BF",
          "ECM / matrix" = "#E41A1C", "Growth / angiogenic" = "#FF7F00", "Stress / cytotoxic" = "#A65628",
          "Lectin" = "#A6CEE3", "Protease / sheddase" = "#FDBF6F", "Others" = "grey70")
lr_pairs_all <- lr_network %>% distinct(from, to) %>% filter(from %in% ligands_use, to %in% receptors_use)
edges_all <- tidyr::crossing(sender = OTHER_SENDERS, from = ligands_use) %>%
  left_join(lr_pairs_all, by = "from") %>%
  left_join(sender_expr %>% select(sender, ligand, Avail_L), by = c("sender","from" = "ligand")) %>%
  filter(!is.na(to))
stack_all <- edges_all %>%
  mutate(sender = stringr::str_squish(as.character(sender)), family = family_of(from)) %>%
  group_by(sender, family) %>% summarise(Avail_L = sum(Avail_L, na.rm = TRUE), .groups = "drop") %>%
  mutate(mode = "All LR pairs") %>% filter(family != 'Others')
write.csv(stack_all %>% mutate(Avail_L = signif(Avail_L, 3)), file.path(FIG, "Fig3H_sender_ligand_families.csv"), row.names = FALSE)
sender_order <- stack_all %>% group_by(sender) %>% summarise(tot = sum(Avail_L, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(tot)) %>% pull(sender)
plot_df <- stack_all %>% mutate(sender = factor(sender, levels = sender_order))
p <- ggplot(plot_df, aes(x = sender, y = Avail_L, fill = family)) +
  geom_col(width = 0.8, na.rm = TRUE) + coord_flip() +
  scale_x_discrete(drop = FALSE) + scale_fill_manual(values = cols, drop = FALSE) +
  labs(x = "Sender Cell", y = "Sum of Ligand Availability in LRs", fill = "Ligand Family",
       title = "Signal to Tregs (Immune + Stromal + Epithelial)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right", text = element_text(family = "DejaVu Sans"), axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
ggsave(file.path(FIG, "Fig3H_sender_ligand_families.png"), p, width = 8, height = 7, dpi = 300, device = ragg::agg_png)

## ---------------- Fig 4I: chord, top 25 LR (best sender per LR) ----------------
df_fac <- df %>% group_by(from, to) %>% slice_max(priority_rank_geom, n = 1, with_ties = FALSE) %>% ungroup()
chord_data <- df_fac %>% select(from, to, sender, priority_rank_geom, Avail_L, log2fc_late) %>%
  filter(from != to) %>% arrange(desc(priority_rank_geom)) %>% slice_head(n = 25)
write.csv(chord_data %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), file.path(FIG, "Fig3I_top25_LR.csv"), row.names = FALSE)
ligand_importance <- chord_data %>% group_by(from) %>% summarise(total_priority = sum(priority_rank_geom)) %>% arrange(total_priority)
receptor_importance <- chord_data %>% group_by(to) %>% summarise(total_priority = sum(priority_rank_geom)) %>% arrange(total_priority)
ligand_cols <- colorRampPalette(c("#fee5d9", "#a50f15"))(nrow(ligand_importance)); names(ligand_cols) <- ligand_importance$from
receptor_cols <- colorRampPalette(c("#deebf7", "#3182bd"))(nrow(receptor_importance)); names(receptor_cols) <- receptor_importance$to
grid_cols <- c(ligand_cols, receptor_cols)
chord_data$transparency <- 1 - scales::rescale(chord_data$log2fc_late, to = c(0.3, 0.9))
# a gene that is both a top ligand and a top receptor would share one sector; flag it
both <- intersect(chord_data$from, chord_data$to); if (length(both)) cat("WARNING genes both ligand and receptor in chord:", both, "\n")

draw_lr_chord(chord_data, grid_cols, file.path(FIG, "Fig3I_LR_chord.png"), "Treg State LRs (All Sender Types)")
saveRDS(list(chord_data = chord_data, grid_cols = grid_cols), file.path(OUT, "fig3H_chord_inputs.rds"))

saveRDS(list(df = df, best_sender_per_LR = best_sender_per_LR, sender_expr = sender_expr, rec_tbl = rec_tbl,
             lig_act_by_sender = lig_act_by_sender, lig_destab_activity = lig_destab_activity,
             geneset_destab = geneset_destab, identity_counts = idt), file.path(OUT, "nichenet_outputs.rds"))
msg("done")
}  # end section F4HI_nichenet_lr


# ======================================================================================================================
# Section S8GH_lr_weighting_schemes: LR priority ranking under five weighting schemes.
# Panels: Supp Fig 8G-H
# ======================================================================================================================
if (.run("S8GH_lr_weighting_schemes")) {
# Stability of the Treg state-transition LR priority ranking under five weighting schemes, on the F4HI_nichenet_lr
# priority table (objects/nichenet_outputs_general_destab/nichenet_destab_late_treg_ligand_strength_imm_logpr.csv;
# component percentile ranks r_destab, r_late, r_avail, r_recEarly, r_recLateUp). NicheNet ligand activities are not
# recomputed: only the composite score is recomputed from the saved component ranks (weight_schemes, compute_priority).
#   Supp 8G  overlap of the top-50 LR rows across the 5 schemes (UpSetR; one row = sender x ligand x receptor), with the
#            weight table as inset (w1..w5 = alpha..eta of compute_priority)
#   Supp 8H  Spearman correlation of ranks across the 5 schemes
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(reshape2); library(gridExtra)
  library(grid); library(UpSetR)})
D <- "/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix"; OUT <- "/share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness"
s3 <- function(df) { df[] <- lapply(df, function(x) if (is.numeric(x)) signif(x, 3) else x); df }

df <- read.csv(file.path(D, "objects/nichenet_outputs_general_destab/nichenet_destab_late_treg_ligand_strength_imm_logpr.csv"))
cat("rows (sender x ligand x receptor):", nrow(df), "| LR pairs:", nrow(distinct(df, from, to)), "\n")
weight_schemes <- list(
  balanced      = c(alpha=1,   beta=1,   gamma=1, delta=1, eta=1),
  activity_heavy= c(alpha=2.5, beta=2,   gamma=1, delta=1, eta=1),
  expr_heavy    = c(alpha=1,   beta=1,   gamma=2, delta=2, eta=1),
  late_bias     = c(alpha=1.5, beta=2.5, gamma=1, delta=1, eta=2),
  conservative  = c(alpha=1.5, beta=1.5, gamma=1, delta=1, eta=1)
)
compute_priority <- function(df, w) {
  num <- (df$r_destab    ^ w["alpha"]) *
         (df$r_late      ^ w["beta"])  *
         (df$r_avail     ^ w["gamma"]) *
         (df$r_recEarly  ^ w["delta"]) *
         (df$r_recLateUp ^ w["eta"])
  den <- sum(w)
  num ^ (1 / den)
}
# the main ranking (priority_rank_geom, Fig 4I) uses alpha=2, beta=1.5, gamma=1, delta=1, eta=1.5 (not one of the 5 schemes)

rankings <- lapply(names(weight_schemes), function(nm) {
  w <- weight_schemes[[nm]]
  df %>% mutate(priority = compute_priority(., w)) %>% arrange(desc(priority)) %>% mutate(rank = row_number()) %>%
    select(sender, from, to, rank) %>% mutate(scheme = nm)
}) %>% bind_rows()
rank_wide <- rankings %>% unite(pair, sender, from, to, remove = FALSE) %>% pivot_wider(names_from = scheme, values_from = rank)
write.csv(s3(rank_wide %>% arrange(balanced)), file.path(OUT, "ext8GH_LR_ranks_all_schemes.csv"), row.names = FALSE)

## ---------------- Supp 8G inset: parameter table
tab <- do.call(rbind, lapply(names(weight_schemes), function(nm) data.frame(scheme = nm, t(weight_schemes[[nm]]))))
colnames(tab) <- c("scheme", "w1", "w2", "w3", "w4", "w5")
tabp <- tab; tabp[, -1] <- lapply(tabp[, -1], function(x) format(x, drop0trailing = TRUE))
tabcsv <- tab; tabcsv$w1_term <- "r_destab (ligand activity on destabilization genes)"; tabcsv$w2_term <- "r_late (ligand activity on late-Treg genes)"
tabcsv$w3_term <- "r_avail (sender ligand availability)"; tabcsv$w4_term <- "r_recEarly (receptor expression in early Tregs)"
tabcsv$w5_term <- "r_recLateUp (receptor induction in late Tregs)"
write.csv(tabcsv, file.path(OUT, "ext8G_LR_parameter_schemes.csv"), row.names = FALSE)
tg <- tableGrob(tabp, rows = NULL, theme = ttheme_minimal(base_size = 13, core = list(fg_params = list(hjust = 0, x = 0.08)),
                                                          colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.08))))
tg <- gtable::gtable_add_grob(tg, grobs = segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(lwd = 1.5)), t = 1, l = 1, r = ncol(tg))
tg <- gtable::gtable_add_grob(tg, grobs = segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(lwd = 1.5)), t = nrow(tg), l = 1, r = ncol(tg))
png(file.path(OUT, "ext8G_LR_parameter_schemes.png"), width = 4.6, height = 2.6, units = "in", res = 300)
grid.arrange(tg, top = textGrob("LR Parameter Schemes", gp = gpar(fontsize = 15)))
dev.off()

## ---------------- Supp 8H: Spearman heatmap (all 5 schemes)
rank_mat5 <- rank_wide %>% dplyr::select(all_of(names(weight_schemes))) %>% as.matrix()
cor5 <- cor(rank_mat5, method = "spearman", use = "pairwise.complete.obs")
write.csv(s3(data.frame(scheme = rownames(cor5), cor5, check.names = FALSE)), file.path(OUT, "ext8H_LR_rank_spearman.csv"), row.names = FALSE)
print(round(cor5, 3))
rank_mat <- rank_wide %>% dplyr::select(all_of(names(weight_schemes))) %>% as.matrix()
cor_mat <- cor(rank_mat, method = "spearman", use = "pairwise.complete.obs")
cor_df <- melt(cor_mat)
lo <- min(0.75, floor(min(cor_mat) * 20) / 20)
p <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 4) +
  scale_fill_gradient(low = "white", high = "#d73027", limits = c(lo, 1)) +   # all rho > 0: sequential scale
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(fill = "Spearman ρ", title = "Spearman Correlation of LR Ranks",
       caption = sprintf("n = %d sender-ligand-receptor rows", nrow(df)))
ggsave(file.path(OUT, "ext8H_LR_rank_spearman.png"), p, width = 5.4, height = 4.6, dpi = 300, bg = "white")

## ---------------- Supp 8G: top-50 overlap (UpSet)
top50 <- lapply(names(weight_schemes), function(nm) rank_wide$pair[rank_wide[[nm]] <= 50]); names(top50) <- names(weight_schemes)
core <- Reduce(intersect, top50)
cat("top-50 rows shared by all 5 schemes:", length(core), "\n")
memb <- data.frame(pair = unique(unlist(top50)))
for (nm in names(top50)) memb[[nm]] <- as.integer(memb$pair %in% top50[[nm]])
memb <- memb %>% left_join(rank_wide %>% select(pair, sender, from, to), by = "pair") %>%
  mutate(n_schemes = rowSums(across(all_of(names(top50))))) %>% arrange(desc(n_schemes), pair)
write.csv(memb, file.path(OUT, "ext8G_top50_LR_overlap.csv"), row.names = FALSE)
# LR-pair level (best sender per ligand-receptor, as Fig 4I)
top50_lr <- lapply(names(weight_schemes), function(nm) {
  rank_wide %>% group_by(from, to) %>% summarise(r = min(.data[[nm]]), .groups = "drop") %>% arrange(r) %>%
    slice_head(n = 50) %>% mutate(lr = paste(from, to, sep = "_")) %>% pull(lr) })
cat("LR-pair level (best sender) top-50 shared by all 5 schemes:", length(Reduce(intersect, top50_lr)), "\n")
jac <- outer(names(top50), names(top50), Vectorize(function(a, b) length(intersect(top50[[a]], top50[[b]]))))
dimnames(jac) <- list(names(top50), names(top50))
write.csv(data.frame(scheme = rownames(jac), jac, check.names = FALSE), file.path(OUT, "ext8G_top50_pairwise_shared.csv"), row.names = FALSE)
print(jac)
# each scheme vs the main ranking (Fig 4H-I weights)
main_rank <- rank(-df$priority_rank_geom, ties.method = "first"); names(main_rank) <- paste(df$sender, df$from, df$to, sep = "_")
# reference top-25 = the Fig 4I chord table
fig3h <- read.csv("/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/Fig3I_top25_LR.csv"); fig3h_lr <- paste(fig3h$from, fig3h$to, sep = "_")
vs_main <- do.call(rbind, lapply(names(weight_schemes), function(nm) {
  r <- rank_wide[[nm]]; m <- main_rank[rank_wide$pair]
  data.frame(scheme = nm, spearman_vs_fig3H_ranking = cor(r, m, method = "spearman"),
             top50_rows_shared_with_fig3H_ranking = length(intersect(rank_wide$pair[r <= 50], names(main_rank)[main_rank <= 50])),
             fig3H_top25_LR_in_scheme_top50_LR = sum(fig3h_lr %in% top50_lr[[which(names(weight_schemes) == nm)]]))
}))
print(vs_main); write.csv(s3(vs_main), file.path(OUT, "ext8GH_schemes_vs_fig3H_ranking.csv"), row.names = FALSE)
png(file.path(OUT, "ext8G_top50_LR_overlap.png"), width = 6.6, height = 5.2, units = "in", res = 300)
print(upset(fromList(top50), sets = rev(names(weight_schemes)), keep.order = TRUE, order.by = "freq", nintersects = NA,
            mainbar.y.label = "Intersection Size", sets.x.label = "Top-50 Rows", text.scale = c(1.4, 1.2, 1.1, 1.1, 1.3, 1.2),
            point.size = 2.8, line.size = 0.9, mainbar.y.max = ceiling(length(core) * 1.35)))
grid.text("Top 50 LR Pairs Overlap", x = 0.62, y = 0.975, gp = gpar(fontsize = 15))
grid.text(sprintf("%d of 50 shared by all 5 schemes", length(core)), x = 0.62, y = 0.93, gp = gpar(fontsize = 10.5))
dev.off()
cat("done\n")
}  # end section S8GH_lr_weighting_schemes
