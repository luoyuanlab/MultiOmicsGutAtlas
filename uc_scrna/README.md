# UC single-cell atlas

Integration and analysis of nine public scRNA-seq studies of healthy colon and ulcerative colitis (441,649 cells, 97 patients, 234 samples). The atlas resolves macrophage and fibroblast states, shows how they shift with inflammation, and predicts signals between SELENOP⁺ macrophages and OGN⁺RSPO3⁺ fibroblasts.

Panel letters refer to the figures (Fig 2 = UC scRNA-seq atlas).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_dataset_preparation.Rmd` | Loads each study, harmonises metadata and the source cell labels | per-study count matrices → per-study Seurat objects | Supp 3B-C (label harmonisation) |
| `02_atlas_integration.R` | Merges the studies, Harmony integration, coarse clustering | 01 objects → integrated atlas | Fig 2A, Supp 1 |
| `03_myeloid_subclustering.R` | Iterative myeloid sub-clustering and myeloid annotation | atlas → `mye_iter3.RDS` | Fig 2B, 2D, Supp 2 |
| `04_fibroblast_subclustering.R` | Stromal sub-clustering; OGN⁺RSPO3⁺ marker feature plots on the stromal UMAP (section 4) | atlas → `con_iter4.RDS`, UMAP and Harmony tables | Fig 2C, Supp 3A |
| `05_other_lineages_subclustering.R` | Sub-clustering of T cells, epithelium and the remaining lineages | atlas → lineage objects | Supp 1 |
| `06_atlas_assembly.R` | Puts the fine labels on the full atlas; writes the RDS and h5ad used downstream | 03-05 → annotated atlas | Fig 2A, Supp 3B-C |
| `07_integration_qc.R` | Exports PCA and Harmony embeddings | atlas → embedding CSVs | Supp Table 3 |
| `08_umaps_integration_qc.ipynb` | Atlas and lineage UMAPs; lineage-wise silhouette and graph connectivity before and after Harmony | 02-07 exports → plots, scores | Fig 2A-C, Supp 1, Supp 2, Supp Table 3 |
| `09_cell_type_composition.py` | scCODA differential abundance (est. FDR 0.1) by condition; macrophage subtype composition per patient and condition (Kruskal-Wallis per subtype, BH across 6 subtypes) | atlas h5ad → scCODA effects and table, `macrophage_composition_per_patient.csv`, panels | Supp 3D, Fig 2E |
| `10_macrophage_selenop_gsea.R` | DEGs and preranked GSEA of CD68⁺ SELENOP⁺ vs SELENOP⁻ myeloid cells, on the myeloid object with gene symbols harmonised across studies | `mye_iter3_symfix.RDS` → DEG table, Supp Table 4, panel | Fig 2F, Supp Table 4 |
| `11_macrophage_selenium_regulators.ipynb` | Stacked violin of selenium-influenced regulators across macrophage subtypes | atlas h5ad, myeloid embedding → panel | Fig 2G |
| `12_fibroblast_activation_trajectory.py` | DPT and Palantir from four candidate roots; ADAMDEC1⁺ fraction by condition (mixed model); disease alignment of each root (patient bootstrap); PROGENy TGF-β / NF-κB / TNF-α along the trajectory | atlas h5ad, stromal embedding, `adamdec1_sample_proportions.csv` → pseudotime CSV, panels, supplementary tables | Fig 2H-J |
| `13_cellchat_myeloid_fibroblast.R` | CellChat between myeloid and fibroblast populations; incoming and outgoing communication patterns | atlas RDS → CellChat objects, plots | Supp 4C-D, Supp 5A-B |
| `14_nichenet_macrophage_to_fibroblast.R` | NicheNet: Mac S+M+S+ ligands for genes up in OGN⁺RSPO3⁺ vs ADAMDEC1⁺ fibroblasts; ligand-regulator links | atlas RDS → ligand activity and regulator tables | Fig 3D (scores), Fig 3G (TNC targets) |
| `15_fibroblast_tnc_targets.ipynb` | Dot plot of TNC and its NicheNet targets across fibroblast subtypes | atlas h5ad, stromal embedding → panel | Fig 3G |
| `16_nichenet_fibroblast_to_macrophage.py` | Patient pseudobulk DESeq2 (Mac S+M+S+ vs resident SELENOP⁺ macrophages), NicheNet ranking of fibroblast ligands, gene-set-size sensitivity, patient bootstrap, ligand-response test | atlas h5ad → panel, supplementary tables | Supp 5C |

## Run order

`01 → 02 → 03, 04, 05 → 06 → 07 → 08`, then any of `09`-`16`. `15` uses the NicheNet scores of `14`, and `14` also feeds `cosmx_uc/10`.

`09` section `S3D_sccoda` runs in the separate pertpy 1.0.2 environment (see the main README).

## Not in this repository

- The gene-symbol harmonisation of the myeloid object read by `10` (`mye_iter3_symfix.RDS`; renamed symbols such as SEPP1 → SELENOP brought to current names from the per-study source data, which are no longer available).
- The per-sample table `adamdec1_sample_proportions.csv` read by `12` (ADAMDEC1⁺ fraction of connective cells per sample).
