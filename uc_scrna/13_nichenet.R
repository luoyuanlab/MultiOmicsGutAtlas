library(nichenetr)
library(Seurat)
library(tidyverse)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR      <- "/path/to/uc/scrna/output"
NICHENET_DIR  <- "/path/to/nichenet/resources"   # NicheNet model files
OUTPUT_DIR    <- "/path/to/uc/scrna/output/nichenet"

# ── Load data ──────────────────────────────────────────────────────────────────
seurat_obj <- readRDS(file.path(DATA_DIR, "uc_atlas_annotated.RDS"))
Idents(seurat_obj) <- "label"

ligand_target_matrix <- readRDS(file.path(NICHENET_DIR, "ligand_target_matrix_nsga2r_final.rds"))
lr_network            <- readRDS(file.path(NICHENET_DIR, "lr_network_human_21122021.rds"))
ligand_tf_matrix      <- readRDS(file.path(NICHENET_DIR, "ligand_tf_matrix.rds"))
gr_network            <- readRDS(file.path(NICHENET_DIR, "gr_network.rds"))

# ── 1. Define sender and receiver cell populations (Fig. 1I) ──────────────────
sender_cells   <- "SELENOP+MMP9+SPP1+ Mac"
receiver_cells <- "OGN+RSPO3+ Fib"

expressed_genes_sender   <- get_expressed_genes(sender_cells,   seurat_obj, pct = 0.05)
expressed_genes_receiver <- get_expressed_genes(receiver_cells, seurat_obj, pct = 0.05)
background_expressed_genes <- expressed_genes_receiver[
  expressed_genes_receiver %in% rownames(ligand_target_matrix)]

# ── 2. Gene set of interest: genes upregulated in OGN+RSPO3+ vs ADAMDEC1+ Fib ─
deg_results <- FindMarkers(seurat_obj,
                           ident.1         = "OGN+RSPO3+ Fib",
                           ident.2         = "ADAMDEC1+ Fib",
                           logfc.threshold = 0.25,
                           min.pct         = 0.05)
geneset_oi <- rownames(deg_results)[deg_results$p_val_adj < 0.1 & deg_results$avg_log2FC > 0]

# ── 3. Identify potential ligands ─────────────────────────────────────────────
expressed_ligands   <- intersect(lr_network$from, expressed_genes_sender)
expressed_receptors <- intersect(lr_network$to,   expressed_genes_receiver)
potential_ligands   <- lr_network %>%
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>% unique()

# ── 4. Ligand activity analysis ───────────────────────────────────────────────
ligand_activities <- predict_ligand_activities(
  geneset                  = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix     = ligand_target_matrix,
  potential_ligands        = potential_ligands
) %>%
  arrange(-pearson) %>%
  mutate(rank = rank(desc(pearson)))

best_upstream_ligands <- ligand_activities %>%
  top_n(50, pearson) %>%
  pull(test_ligand)

# ── 5. Downstream TF predictions ──────────────────────────────────────────────
active_ligand_tf <- best_upstream_ligands %>%
  lapply(function(ligand) {
    tf_scores <- ligand_tf_matrix[ligand, ] %>% sort(decreasing = TRUE) %>% head(50)
    data.frame(ligand = ligand, tf = names(tf_scores),
               score = as.numeric(tf_scores), stringsAsFactors = FALSE)
  }) %>%
  bind_rows() %>%
  filter(score > quantile(score, 0.5))

predicted_tfs_expressed <- intersect(unique(active_ligand_tf$tf),
                                     expressed_genes_receiver)

# ── 6. TF-target gene network ─────────────────────────────────────────────────
tf_target_network <- gr_network %>%
  filter(from %in% predicted_tfs_expressed) %>%
  group_by(from) %>%
  slice_head(n = 50) %>%
  ungroup()

# ── 7. Save outputs ───────────────────────────────────────────────────────────
write.csv(ligand_activities,   file.path(OUTPUT_DIR, "nichenet_ligand_activities.csv"),  row.names = FALSE)
write.csv(active_ligand_tf,    file.path(OUTPUT_DIR, "nichenet_ligand_tf_predictions.csv"), row.names = FALSE)
write.csv(data.frame(tf = predicted_tfs_expressed),
          file.path(OUTPUT_DIR, "nichenet_predicted_tfs.csv"), row.names = FALSE)
write.csv(tf_target_network,   file.path(OUTPUT_DIR, "nichenet_tf_target_network.csv"),  row.names = FALSE)
