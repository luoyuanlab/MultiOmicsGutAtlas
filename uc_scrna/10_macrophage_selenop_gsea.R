# UC scRNA-seq atlas: GSEA of CD68+ SELENOP+ vs SELENOP- myeloid cells.
#
# Panels: Fig 2F (GSEA, CD68+ SELENOP+ vs SELENOP- myeloid cells); Supp Table 4 (all GSEA results).
# Input: subcluster5/mye_iter3_symfix.RDS (myeloid object, gene symbols harmonised across studies).
# Outputs: Fig1F_gsea_selenop_symfix/ (DEG CSV, Supp Table 4, plotted values, panel).
# Run: Rscript <file>, or RUN_SECTIONS=<id>[,<id>...] Rscript <file>.
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {   # all sections run in file order unless RUN_SECTIONS=id1,id2,... is set
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}


# ======================================================================================================================
# Section F2F_selenop_degs: DEGs of SELENOP+ vs SELENOP- cells among CD68+ myeloid cells
# Panels: Fig 2F (DEGs)
# ======================================================================================================================
if (.run("F2F_selenop_degs")) {
# CD68+ myeloid cells (CD68 > 0) are split into SELENOP+ (SELENOP > 0) and SELENOP- groups; Seurat FindMarkers with no
# expression or fold-change pre-filter. Output: cd68_selenop_pos_neg_degs.csv
suppressMessages({library(Seurat); library(dplyr)})
OUT <- "/share/fsmresfiles/UC/revision_jenny/Fig1F_gsea_selenop_symfix"
uc_mye_iter3 <- readRDS("/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/subcluster5/mye_iter3_symfix.RDS")

CD68_cells  <- subset(uc_mye_iter3, subset = CD68 > 0)
CD68_cells[["CellName"]] <- colnames(CD68_cells)
selenop_pos <- Cells(subset(CD68_cells, subset = SELENOP > 0))
CD68_cells@meta.data <- CD68_cells@meta.data %>%
  mutate(SELENOP_exp = ifelse(CellName %in% selenop_pos, "Pos", "Neg"))
CD68_cells <- SetIdent(CD68_cells, value = "SELENOP_exp")
deg_selenop <- FindMarkers(CD68_cells, ident.1 = "Pos", ident.2 = "Neg",
                           min.pct = -Inf, logfc.threshold = -Inf,
                           min.cells.feature = 1, min.cells.group = 1)
write.csv(deg_selenop, file.path(OUT, "cd68_selenop_pos_neg_degs.csv"))
}  # end section F2F_selenop_degs


# ======================================================================================================================
# Section F2F_gsea_supp_table4: preranked GSEA (GO, KEGG, WikiPathways, Reactome), Supp Table 4 and the Fig 2F pathway values
# Panels: Fig 2F (values), Supp Table 4
# ======================================================================================================================
if (.run("F2F_gsea_supp_table4")) {
# Genes with |log2FC| > 0.5 are ranked by log2FC (Entrez IDs) for gseGO / gseKEGG / gseWP / gsePathway; pvalueCutoff = 1
# so every pathway shown in Fig 2F has an NES (significance is shown by colour). Leading-edge genes are given as symbols.
# Outputs: supp_table4_gsea_all.csv (Supp Table 4), fig1f_gsea_plotted.csv (Fig 2F values)
suppressMessages({library(dplyr); library(DOSE); library(org.Hs.eg.db); library(AnnotationDbi)
  library(clusterProfiler); library(ReactomePA)})
OUT <- "/share/fsmresfiles/UC/revision_jenny/Fig1F_gsea_selenop_symfix"
degs <- read.csv(file.path(OUT, "cd68_selenop_pos_neg_degs.csv"))
degs <- degs %>% filter(avg_log2FC > 0.5 | avg_log2FC < -0.5) %>% select(X, avg_log2FC) %>% rename(gene = X)
degs$entrez <- as.vector(mapIds(org.Hs.eg.db, keys = degs$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
degs <- degs %>% filter(!is.na(entrez)) %>% arrange(-avg_log2FC)
entrez_degs <- setNames(degs$avg_log2FC, degs$entrez)
set.seed(1)
ego   <- gseGO(entrez_degs, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, verbose = FALSE, seed = TRUE)
mkk   <- gseKEGG(entrez_degs, organism = "hsa", pvalueCutoff = 1, verbose = FALSE, seed = TRUE)
wp    <- gseWP(entrez_degs, organism = "Homo sapiens", pvalueCutoff = 1, verbose = FALSE, seed = TRUE)
react <- gsePathway(entrez_degs, pvalueCutoff = 1, verbose = FALSE, seed = TRUE)
sel <- function(x, db) x@result %>% select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, core_enrichment) %>% mutate(Database = db)
all_gsea <- bind_rows(sel(react, "Reactome"), sel(ego, "GO"), sel(mkk, "KEGG"), sel(wp, "WikiPathways")) %>% arrange(enrichmentScore)
ue <- unique(unlist(strsplit(all_gsea$core_enrichment, "/")))
sym <- mapIds(org.Hs.eg.db, keys = ue, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
all_gsea$leading_edge_genes <- vapply(strsplit(all_gsea$core_enrichment, "/"), function(i) paste(na.omit(sym[i]), collapse = "/"), "")
supp4 <- all_gsea %>%
  select(ID, Database, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, leading_edge_genes) %>%
  mutate(across(c(enrichmentScore, NES, pvalue, p.adjust, qvalue), ~ signif(.x, 3))) %>% arrange(NES)
write.csv(supp4, file.path(OUT, "supp_table4_gsea_all.csv"), row.names = FALSE)

fig2f_ids <- c("WP623", "WP28", "WP3931", "R-HSA-917937", "hsa04979", "hsa04668", "WP382", "R-HSA-181438",
  "R-HSA-6783783", "GO:0007249", "WP4496", "GO:0032615", "R-HSA-913531", "R-HSA-2408557", "GO:0034097")  # pathways shown in Fig 2F
f <- all_gsea %>% filter(ID %in% fig2f_ids) %>%
  mutate(Description = recode(Description, "response to cytokine" = "Response to cytokine",
    "interleukin-12 production" = "Interleukin-12 production",
    "canonical NF-kappaB signal transduction" = "Canonical NF-kappaB signalling",
    "Toll Like Receptor 2 (TLR2) Cascade" = "Toll Like Receptor 2 Cascade",
    "Embryonic stem cell pluripotency pathways" = "Stem cell pluripotency pathways"),
    Description = paste0(ID, ": ", Description)) %>%
  arrange(NES)
write.csv(f %>% select(ID, Description, Database, setSize, NES, pvalue, p.adjust) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))),
          file.path(OUT, "fig1f_gsea_plotted.csv"), row.names = FALSE)
}  # end section F2F_gsea_supp_table4


# ======================================================================================================================
# Section F2F_gsea_barplot: NES bar plot of the Fig 2F pathways, coloured by direction and FDR q < 0.1
# Panels: Fig 2F
# ======================================================================================================================
if (.run("F2F_gsea_barplot")) {
# Reads fig1f_gsea_plotted.csv (section F2F_gsea_supp_table4). Output: fig1f_gsea_selenop.png
suppressMessages({library(dplyr); library(ggplot2)})
OUT <- "/share/fsmresfiles/UC/revision_jenny/Fig1F_gsea_selenop_symfix"
f <- read.csv(file.path(OUT, "fig1f_gsea_plotted.csv")) %>% arrange(NES) %>%
  mutate(Pathway = factor(Description, levels = Description[order(NES)]),
         Color_Group = case_when(NES > 0 & p.adjust < 0.1 ~ "Higher in SELENOP⁺ (FDR q < 0.1)", NES > 0 ~ "Higher in SELENOP⁺ (q ≥ 0.1)",
                                 NES < 0 & p.adjust < 0.1 ~ "Higher in SELENOP⁻ (FDR q < 0.1)", TRUE ~ "Higher in SELENOP⁻ (q ≥ 0.1)"))
color_map <- c("Higher in SELENOP⁺ (FDR q < 0.1)" = "darkred", "Higher in SELENOP⁺ (q ≥ 0.1)" = "lightcoral",
               "Higher in SELENOP⁻ (FDR q < 0.1)" = "darkblue", "Higher in SELENOP⁻ (q ≥ 0.1)" = "lightblue")
offset <- 0.02 * max(abs(f$NES))
p <- ggplot(f, aes(x = NES, y = Pathway, fill = Color_Group)) + geom_col() +
  scale_fill_manual(values = color_map, breaks = names(color_map), drop = FALSE) +
  geom_text(data = filter(f, NES > 0), aes(label = Pathway), x = -offset, hjust = 1, size = 3) +
  geom_text(data = filter(f, NES < 0), aes(label = Pathway), x = offset, hjust = 0, size = 3) +
  scale_x_continuous(limits = c(-3, 3)) + coord_cartesian(clip = "off") +
  theme_classic() + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank(), legend.position = "top", panel.grid = element_blank(), legend.title = element_text(hjust = 0.5, size = 10), plot.margin = margin(0.5, 1, 0.5, 1, "cm")) +
  guides(fill = guide_legend(nrow = 2, title.position = "top")) +
  labs(x = "Normalized Enrichment Score (NES)\nSELENOP⁺ vs SELENOP⁻ CD68⁺ myeloid cells", fill = "Pathway enrichment (CD68⁺ myeloid cells, GSEA)")
ggsave(file.path(OUT, "fig1f_gsea_selenop.png"), plot = p, width = 9, height = 5, dpi = 300)
}  # end section F2F_gsea_barplot
