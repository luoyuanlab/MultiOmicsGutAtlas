# CD scRNA-seq: NicheNet ligand-receptor prioritization of signals linked to loss of Treg suppression.
#
# Panels: Fig 4H (ligand availability by function per sender), Fig 4I (top 25 LR chord), Supp 8G-H (weighting schemes).
#   Section immune_sender_ranking writes the immune-sender ranking from which cosmx_cd/10_treg_suppression_score.py
#   takes the antigen-presentation LR module.
# Inputs:  objects/nichenet_inputs/ (15_treg_ligand_receptor_inputs.py); objects/nichenet_outputs_lossact/ (two gene-set
#   tables for F4HI_nichenet_lr, not in repo); NicheNet v2 model files in /share/fsmresfiles/ydn4687/nichenet/.
# Outputs: objects/nichenet_outputs/ and objects/nichenet_outputs_general_destab/ (priority tables); Fig 4H-I panels and
#   CSVs in Fig3D-I_treg_pseudotime_signals/; Supp 8G-H panels and CSVs in SuppFig8A-I_treg_pseudotime_robustness/.
# Run:     Rscript <file>  or  RUN_SECTIONS=<id>[,<id>] Rscript <file>
#   (immune_sender_ranking -> F4HI_nichenet_lr -> S8GH_lr_weighting_schemes)
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}
BASE <- "/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix"
NN   <- "/share/fsmresfiles/ydn4687/nichenet"
eps  <- 1e-6

load_nn_libs <- function() suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble)
  library(Matrix); library(purrr); library(PRROC); library(scales)
  library(nichenetr); library(stringr); library(forcats); library(ggplot2)
  library(arrow); library(circlize); library(grid)
})

beta_shrink <- function(k, n, alpha = 1, beta = 1) (k + alpha) / (n + alpha + beta)
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}
pctrank <- function(x) {
  n <- sum(!is.na(x)); if (n <= 1) return(rep(0.5, length(x)))
  rx <- rank(x, na.last = "keep", ties.method = "average"); as.numeric((rx - 1) / pmax(n - 1, 1))
}

# Seurat object from the dense parquet (built in column chunks); Tregs bins 0-1 = Treg_Early, bins 2-4 = Treg_Late
read_nn_obj <- function(meta_f, counts_f) {
  meta <- as.data.frame(read_parquet(meta_f))
  tbl  <- read_parquet(counts_f, as_data_frame = FALSE)
  bad_gene <- "__index_level_0__"
  genes <- setdiff(names(tbl), bad_gene)
  if (bad_gene %in% colnames(meta)) meta <- meta[, colnames(meta) != bad_gene, drop = FALSE]
  chunks <- split(genes, ceiling(seq_along(genes) / 1000))
  counts_sp <- do.call(cbind, lapply(chunks, function(g) Matrix(as.matrix(as.data.frame(tbl[, g])), sparse = TRUE)))
  rownames(counts_sp) <- meta$index
  rm(tbl); invisible(gc())
  obj <- CreateSeuratObject(counts = t(counts_sp))
  obj$ct_lr <- meta$ct_lr
  rm(counts_sp); invisible(gc())
  treg_mask <- grepl("^Treg_B[0-4]$", obj$ct_lr)
  pt_bin <- rep(NA_character_, ncol(obj))
  pt_bin[treg_mask] <- sub("^Treg_(B[0-4])$", "\\1", obj$ct_lr[treg_mask])
  nn_id <- obj$ct_lr
  nn_id[treg_mask & pt_bin %in% c("B0", "B1")] <- "Treg_Early"
  nn_id[treg_mask & pt_bin %in% c("B2", "B3", "B4")] <- "Treg_Late"
  obj$nn_id <- factor(nn_id)
  Idents(obj) <- obj$nn_id
  NormalizeData(obj, layer = "counts", new.layer = "data")
}

# Steps 1-3: receiver DE and receptor detection, sender ligand availability, NicheNet ligand activity on late-Treg genes
nn_receiver_sender <- function(obj, senders, lr_network, ligand_target_matrix, weighted_networks) {
  mat <- GetAssayData(obj, layer = "data")
  early_cells <- colnames(obj)[obj$nn_id == "Treg_Early"]
  late_cells  <- colnames(obj)[obj$nn_id == "Treg_Late"]
  ligands_use   <- intersect(unique(lr_network$from), rownames(mat))
  receptors_use <- intersect(unique(lr_network$to),   rownames(mat))

  rec_obj <- subset(obj, idents = c("Treg_Early", "Treg_Late"))
  Idents(rec_obj) <- factor(Idents(rec_obj), levels = c("Treg_Early", "Treg_Late"))
  de_rec <- FindMarkers(rec_obj, ident.1 = "Treg_Late", ident.2 = "Treg_Early",
                        logfc.threshold = 0, min.pct = 0.05, only.pos = FALSE) %>% tibble::rownames_to_column("gene")
  k_early <- Matrix::rowSums(mat[receptors_use, early_cells, drop = FALSE] > 0)
  rec_frac_early_shrunk <- beta_shrink(k_early, length(early_cells), alpha = 1, beta = 1)
  rec_tbl <- de_rec %>%
    dplyr::filter(gene %in% receptors_use) %>%
    dplyr::mutate(log2fc_late = avg_log2FC, sig_late = p_val_adj < 0.1 & log2fc_late >= 0.25) %>%
    dplyr::left_join(tibble::tibble(gene = names(rec_frac_early_shrunk),
                                    rec_frac_early = as.numeric(rec_frac_early_shrunk)), by = "gene")

  sender_expr <- dplyr::bind_rows(lapply(senders, function(s) {
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

  nn_results <- setNames(lapply(senders, function(s) nichenet_seuratobj_cluster_de(
    seurat_obj = obj, receiver_reference = "Treg_Early", receiver_affected = "Treg_Late",
    sender = s, ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
    weighted_networks = weighted_networks, expression_pct = 0.05, lfc_cutoff = 0.25,
    geneset = "up", cutoff_visualization = 0.33)), senders)
  lig_act_by_sender <- lapply(names(nn_results), function(s) {
    nn_results[[s]]$ligand_activities %>%
      dplyr::transmute(from = .data$test_ligand, A_late = pmin(pmax(.data$aupr_corrected, 0), 1)) %>%
      dplyr::mutate(sender = s)
  }) |> dplyr::bind_rows()

  expr_frac_receiver <- Matrix::rowMeans(mat[, c(early_cells, late_cells), drop = FALSE] > 0)
  list(de_rec = de_rec, rec_tbl = rec_tbl, sender_expr = sender_expr, lig_act_by_sender = lig_act_by_sender,
       ligands_use = ligands_use, receptors_use = receptors_use,
       bg_expressed_genes = names(expr_frac_receiver)[expr_frac_receiver >= 0.01])
}

# Step 4: NicheNet ligand activity on the Treg state-transition gene set
nn_destab_activity <- function(geneset_destab, bg_expressed_genes, lr_network, ligand_target_matrix) {
  expressed_receptors <- intersect(unique(lr_network$to), bg_expressed_genes)
  potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
  raw <- nichenetr::predict_ligand_activities(
    geneset = geneset_destab, background_expressed_genes = bg_expressed_genes,
    ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  tibble::tibble(ligand = potential_ligands) %>%
    dplyr::left_join(raw %>% dplyr::select(ligand = test_ligand, A_destab = aupr_corrected), by = "ligand")
}

# Steps 5-6: LR pairs per sender; composite priority = weighted geometric mean of percentile ranks
#   (alpha = 2, beta = 1.5, gamma = 1, delta = 1, eta = 1.5); writes the priority table and its best sender per LR pair
nn_priority <- function(senders, r, lig_destab_activity, lr_network, OUT) {
  lr_by_sender <- tidyr::crossing(sender = senders, from = r$ligands_use) %>%
    dplyr::left_join(lr_network %>% dplyr::distinct(from, to) %>% dplyr::filter(from %in% r$ligands_use, to %in% r$receptors_use), by = "from") %>%
    dplyr::left_join(r$lig_act_by_sender, by = c("sender", "from")) %>%
    dplyr::mutate(A_late = dplyr::coalesce(A_late, eps)) %>%
    dplyr::left_join(lig_destab_activity, by = c("from" = "ligand")) %>%
    dplyr::mutate(A_destab = dplyr::coalesce(A_destab, eps)) %>%
    filter(sender != "Neuron")
  alpha <- 2; beta <- 1.5; gamma <- 1; delta <- 1; eta <- 1.5
  df <- lr_by_sender %>%
    dplyr::left_join(r$sender_expr %>% dplyr::select(sender, ligand, Avail_L, L_pi_hat, L_mu_pos01), by = c("sender", "from" = "ligand")) %>%
    dplyr::left_join(r$rec_tbl %>% dplyr::select(to = gene, log2fc_late, rec_frac_early), by = "to") %>%
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
  df
}

# LR chord diagram (Fig 4I) with the title centred over the whole figure and a key panel
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


# ---- Section immune_sender_ranking: immune senders -> Treg LR priority; state-transition gene set = NF-kB /
#      inflammatory destabilisation genes. Its best-sender table is read by cosmx_cd/10_treg_suppression_score.py.
if (.run("immune_sender_ranking")) {
load_nn_libs()
set.seed(1)
OUT <- file.path(BASE, "objects", "nichenet_outputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
obj <- read_nn_obj(file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_meta.parquet"),
                   file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_counts.parquet"))
lr_network <- readRDS(file.path(NN, "lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file.path(NN, "ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(file.path(NN, "weighted_networks_nsga2r_final.rds"))
SENDERS <- c('B', 'DC', 'Granulocyte', 'Macrophage', 'Neuron', 'NK', 'Other T', 'Plasma')

r <- nn_receiver_sender(obj, SENDERS, lr_network, ligand_target_matrix, weighted_networks)
destab_prog_use <- c("NFKBIA","NFKBIZ","TNFAIP3","RELB","TNF","IFNG","OSM","CXCL8",
                     "CCL3","CCL4","CCL5","CCL20","IL23A","TNFSF9","TNFSF14","TNFSF12","TNFSF13B")
geneset_destab <- intersect(destab_prog_use, r$bg_expressed_genes)
lig_destab_activity <- nn_destab_activity(geneset_destab, r$bg_expressed_genes, lr_network, ligand_target_matrix)
nn_priority(SENDERS, r, lig_destab_activity, lr_network, OUT)
}


# ---- Section F4HI_nichenet_lr: immune, fibroblast, endothelial and epithelial senders (Fig 4H-I). State-transition gene
#      set = loss-of-activation genes (late vs early: p_adj < 0.05, log2FC <= -0.25, >= 1% of Tregs; up to 250) minus lineage
#      contaminants (Treg / max other lineage ratio < 0.2) plus destabilisation markers up in bin 4 vs bin 0 (Wilcoxon p < 0.05).
if (.run("F4HI_nichenet_lr")) {
load_nn_libs()
set.seed(1)
OUT <- file.path(BASE, "objects", "nichenet_outputs_general_destab")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG <- "/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals"   # Fig 4H-I panels and CSVs
obj <- read_nn_obj(file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_general_meta.parquet"),
                   file.path(BASE, "objects/nichenet_inputs/strict_treg_imm_general_counts.parquet"))
lr_network <- readRDS(file.path(NN, "lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file.path(NN, "ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(file.path(NN, "weighted_networks_nsga2r_final.rds"))
SENDERS <- c('B', 'DC', 'Granulocyte', 'Macrophage', 'Neuron', 'NK', 'Other T', 'Plasma', 'Fibroblast', 'Endothelial', 'Epithelial')

r <- nn_receiver_sender(obj, SENDERS, lr_network, ligand_target_matrix, weighted_networks)
bg <- r$bg_expressed_genes
geneset_loss <- r$de_rec %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC <= -0.25, gene %in% bg) %>%
  dplyr::arrange(avg_log2FC) %>% dplyr::slice_head(n = 250) %>% dplyr::pull(gene)
lin <- read.csv(file.path(BASE, "objects/nichenet_outputs_lossact/loss_geneset_lineage_check.csv"), row.names = 1)
contam <- rownames(lin)[lin$ratio_treg_over_max_other < 0.2]
dm <- read.csv(file.path(BASE, "objects/nichenet_outputs_lossact/destab_markers_bin4_vs_bin0.csv"))
destab_up <- intersect(dm$gene[dm$wilcox_p < 0.05 & dm$median_delta > 0], bg)
geneset_destab <- union(setdiff(geneset_loss, contam), destab_up)
writeLines(geneset_destab, file.path(OUT, "state_transition_geneset.txt"))
lig_destab_activity <- nn_destab_activity(geneset_destab, bg, lr_network, ligand_target_matrix)
df <- nn_priority(SENDERS, r, lig_destab_activity, lr_network, OUT)
SENDERS <- setdiff(SENDERS, "Neuron")

## Fig 4H: ligand availability by function per sender
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
lr_pairs_all <- lr_network %>% distinct(from, to) %>% filter(from %in% r$ligands_use, to %in% r$receptors_use)
stack_all <- tidyr::crossing(sender = SENDERS, from = r$ligands_use) %>%
  left_join(lr_pairs_all, by = "from") %>%
  left_join(r$sender_expr %>% select(sender, ligand, Avail_L), by = c("sender", "from" = "ligand")) %>%
  filter(!is.na(to)) %>%
  mutate(sender = stringr::str_squish(as.character(sender)), family = family_of(from)) %>%
  group_by(sender, family) %>% summarise(Avail_L = sum(Avail_L, na.rm = TRUE), .groups = "drop") %>%
  mutate(mode = "All LR pairs") %>% filter(family != 'Others')
write.csv(stack_all %>% mutate(Avail_L = signif(Avail_L, 3)), file.path(FIG, "Fig3H_sender_ligand_families.csv"), row.names = FALSE)
sender_order <- stack_all %>% group_by(sender) %>% summarise(tot = sum(Avail_L, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(tot)) %>% pull(sender)
p <- ggplot(stack_all %>% mutate(sender = factor(sender, levels = sender_order)), aes(x = sender, y = Avail_L, fill = family)) +
  geom_col(width = 0.8, na.rm = TRUE) + coord_flip() +
  scale_x_discrete(drop = FALSE) + scale_fill_manual(values = cols, drop = FALSE) +
  labs(x = "Sender Cell", y = "Sum of Ligand Availability in LRs", fill = "Ligand Family",
       title = "Signal to Tregs (Immune + Stromal + Epithelial)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right", text = element_text(family = "DejaVu Sans"), axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
ggsave(file.path(FIG, "Fig3H_sender_ligand_families.png"), p, width = 8, height = 7, dpi = 300, device = ragg::agg_png)

## Fig 4I: chord, top 25 LR (best sender per LR)
chord_data <- df %>% group_by(from, to) %>% slice_max(priority_rank_geom, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(from, to, sender, priority_rank_geom, Avail_L, log2fc_late) %>%
  filter(from != to) %>% arrange(desc(priority_rank_geom)) %>% slice_head(n = 25)
write.csv(chord_data %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), file.path(FIG, "Fig3I_top25_LR.csv"), row.names = FALSE)
ligand_importance <- chord_data %>% group_by(from) %>% summarise(total_priority = sum(priority_rank_geom)) %>% arrange(total_priority)
receptor_importance <- chord_data %>% group_by(to) %>% summarise(total_priority = sum(priority_rank_geom)) %>% arrange(total_priority)
ligand_cols <- colorRampPalette(c("#fee5d9", "#a50f15"))(nrow(ligand_importance)); names(ligand_cols) <- ligand_importance$from
receptor_cols <- colorRampPalette(c("#deebf7", "#3182bd"))(nrow(receptor_importance)); names(receptor_cols) <- receptor_importance$to
chord_data$transparency <- 1 - scales::rescale(chord_data$log2fc_late, to = c(0.3, 0.9))
draw_lr_chord(chord_data, c(ligand_cols, receptor_cols), file.path(FIG, "Fig3I_LR_chord.png"), "Treg State LRs (All Sender Types)")
}


# ---- Section S8GH_lr_weighting_schemes: composite priority recomputed from the saved component ranks of
#      F4HI_nichenet_lr under five weighting schemes (NicheNet activities are not recomputed).
#      Supp 8G: top-50 overlap (UpSet) + weight table inset; Supp 8H: Spearman correlation of ranks.
if (.run("S8GH_lr_weighting_schemes")) {
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(reshape2); library(gridExtra)
  library(grid); library(UpSetR)})
OUT <- "/share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness"
s3 <- function(df) { df[] <- lapply(df, function(x) if (is.numeric(x)) signif(x, 3) else x); df }

df <- read.csv(file.path(BASE, "objects/nichenet_outputs_general_destab/nichenet_destab_late_treg_ligand_strength_imm_logpr.csv"))
weight_schemes <- list(
  balanced      = c(alpha=1,   beta=1,   gamma=1, delta=1, eta=1),
  activity_heavy= c(alpha=2.5, beta=2,   gamma=1, delta=1, eta=1),
  expr_heavy    = c(alpha=1,   beta=1,   gamma=2, delta=2, eta=1),
  late_bias     = c(alpha=1.5, beta=2.5, gamma=1, delta=1, eta=2),
  conservative  = c(alpha=1.5, beta=1.5, gamma=1, delta=1, eta=1)
)
compute_priority <- function(df, w) {
  num <- (df$r_destab ^ w["alpha"]) * (df$r_late ^ w["beta"]) * (df$r_avail ^ w["gamma"]) *
         (df$r_recEarly ^ w["delta"]) * (df$r_recLateUp ^ w["eta"])
  num ^ (1 / sum(w))
}
rankings <- lapply(names(weight_schemes), function(nm) {
  w <- weight_schemes[[nm]]
  df %>% mutate(priority = compute_priority(., w)) %>% arrange(desc(priority)) %>% mutate(rank = row_number()) %>%
    select(sender, from, to, rank) %>% mutate(scheme = nm)
}) %>% bind_rows()
rank_wide <- rankings %>% unite(pair, sender, from, to, remove = FALSE) %>% pivot_wider(names_from = scheme, values_from = rank)
write.csv(s3(rank_wide %>% arrange(balanced)), file.path(OUT, "ext8GH_LR_ranks_all_schemes.csv"), row.names = FALSE)

## Supp 8G inset: parameter table
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

## Supp 8H: Spearman heatmap (all 5 schemes)
cor_mat <- cor(rank_wide %>% dplyr::select(all_of(names(weight_schemes))) %>% as.matrix(), method = "spearman", use = "pairwise.complete.obs")
write.csv(s3(data.frame(scheme = rownames(cor_mat), cor_mat, check.names = FALSE)), file.path(OUT, "ext8H_LR_rank_spearman.csv"), row.names = FALSE)
lo <- min(0.75, floor(min(cor_mat) * 20) / 20)
p <- ggplot(melt(cor_mat), aes(Var1, Var2, fill = value)) +
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

## Supp 8G: top-50 overlap (UpSet)
top50 <- lapply(names(weight_schemes), function(nm) rank_wide$pair[rank_wide[[nm]] <= 50]); names(top50) <- names(weight_schemes)
core <- Reduce(intersect, top50)
memb <- data.frame(pair = unique(unlist(top50)))
for (nm in names(top50)) memb[[nm]] <- as.integer(memb$pair %in% top50[[nm]])
memb <- memb %>% left_join(rank_wide %>% select(pair, sender, from, to), by = "pair") %>%
  mutate(n_schemes = rowSums(across(all_of(names(top50))))) %>% arrange(desc(n_schemes), pair)
write.csv(memb, file.path(OUT, "ext8G_top50_LR_overlap.csv"), row.names = FALSE)
jac <- outer(names(top50), names(top50), Vectorize(function(a, b) length(intersect(top50[[a]], top50[[b]]))))
dimnames(jac) <- list(names(top50), names(top50))
write.csv(data.frame(scheme = rownames(jac), jac, check.names = FALSE), file.path(OUT, "ext8G_top50_pairwise_shared.csv"), row.names = FALSE)
# each scheme vs the main ranking (Fig 4H-I weights; top-25 = Fig 4I chord table)
top50_lr <- lapply(names(weight_schemes), function(nm) {
  rank_wide %>% group_by(from, to) %>% summarise(r = min(.data[[nm]]), .groups = "drop") %>% arrange(r) %>%
    slice_head(n = 50) %>% mutate(lr = paste(from, to, sep = "_")) %>% pull(lr) })
main_rank <- rank(-df$priority_rank_geom, ties.method = "first"); names(main_rank) <- paste(df$sender, df$from, df$to, sep = "_")
fig3h <- read.csv("/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/Fig3I_top25_LR.csv"); fig3h_lr <- paste(fig3h$from, fig3h$to, sep = "_")
vs_main <- do.call(rbind, lapply(names(weight_schemes), function(nm) {
  r <- rank_wide[[nm]]; m <- main_rank[rank_wide$pair]
  data.frame(scheme = nm, spearman_vs_fig3H_ranking = cor(r, m, method = "spearman"),
             top50_rows_shared_with_fig3H_ranking = length(intersect(rank_wide$pair[r <= 50], names(main_rank)[main_rank <= 50])),
             fig3H_top25_LR_in_scheme_top50_LR = sum(fig3h_lr %in% top50_lr[[which(names(weight_schemes) == nm)]]))
}))
write.csv(s3(vs_main), file.path(OUT, "ext8GH_schemes_vs_fig3H_ranking.csv"), row.names = FALSE)
png(file.path(OUT, "ext8G_top50_LR_overlap.png"), width = 6.6, height = 5.2, units = "in", res = 300)
print(upset(fromList(top50), sets = rev(names(weight_schemes)), keep.order = TRUE, order.by = "freq", nintersects = NA,
            mainbar.y.label = "Intersection Size", sets.x.label = "Top-50 Rows", text.scale = c(1.4, 1.2, 1.1, 1.1, 1.3, 1.2),
            point.size = 2.8, line.size = 0.9, mainbar.y.max = ceiling(length(core) * 1.35)))
grid.text("Top 50 LR Pairs Overlap", x = 0.62, y = 0.975, gp = gpar(fontsize = 15))
grid.text(sprintf("%d of 50 shared by all 5 schemes", length(core)), x = 0.62, y = 0.93, gp = gpar(fontsize = 10.5))
dev.off()
}
