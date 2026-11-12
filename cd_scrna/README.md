# CD single-cell atlas and Treg states

Integration of public scRNA-seq studies of control, inactive and active Crohn's disease (793,762 cells, 267 patients, 146 with CD), followed by an analysis of regulatory T cells (Tregs). The question is how Treg states change along a pseudotime trajectory in CD, which suppressive programs are lost, and which cells send the signals that late, low-suppressive Tregs receive.

Panel letters refer to the figures (Fig 4 = CD Treg pseudotime).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_dataset_preparation.Rmd` | Loads each study and harmonises metadata | per-study data → Seurat objects | Supp 7 |
| `02_atlas_integration.R` | Merges the studies, Harmony integration, coarse clustering | 01 objects → integrated atlas | Fig 4A, Supp 7A-B |
| `03_myeloid_subclustering.R` | Myeloid sub-clustering and annotation | atlas → myeloid object, UMAP tables | Supp 7C |
| `04_fibroblast_subclustering.R` | Stromal sub-clustering and annotation | atlas → stromal object, UMAP tables | Supp 7E |
| `05_t_cell_subclustering.R` | T/NK sub-clustering and annotation, including Treg subtypes | atlas → T cell object, UMAP tables | Fig 4B |
| `06_other_lineages_subclustering.R` | Remaining lineages | atlas → lineage objects | Fig 4A |
| `07_epithelial_subclustering.R` | Epithelial sub-clustering and annotation | atlas → epithelial object | Supp 7D |
| `08_atlas_assembly.R` | Merges the annotated lineages into the atlas | 03-07 → atlas and count matrix | Fig 4A |
| `09_atlas_h5ad_conversion.ipynb` | Builds the h5ad (X rows and `.obs` in sorted-barcode order) and attaches the labels | 08 export → `cd_atlas_annotated.h5ad` | input of 10-15 |
| `10_umaps.ipynb` | Atlas and lineage UMAPs | UMAP tables, 09 h5ad → plots | Fig 4A-B, Supp 7A-E |
| `11_sccoda_composition.ipynb` | scCODA of cell-type abundance by location and activity (overall, T cells, myeloid cells) | 09 h5ad → effects, forest plots | Fig 4C, Supp 7G |
| `12_treg_pseudotime.py` | Treg input objects; strict Treg selection, Harmony, Leiden, diffusion pseudotime from the core Treg sub-cluster; module UMAPs, markers by bin, Suppressive Core Score vs pseudotime (LMM, sign test) | 09 h5ad → `treg_pseudotime.h5ad`, panels | Fig 4D-F |
| `13_treg_pseudotime_root_comparison.py` | Th17-like Tregs vs conventional Th17 cells; Leiden sub-clusters and root choice; Palantir from the main root; ten alternative roots with DPT and Palantir | `treg_pseudotime.h5ad` → panels, CSVs | Supp 8A-E |
| `14_treg_pathway_tf_activity.py` | Pathway activity per pseudotime quintile (patient-paired); CollecTRI TF activity along pseudotime | `treg_pseudotime.h5ad` → panels, TF table | Fig 4G, Supp 8F, Supp Table 10 (pathway gene sets) |
| `15_treg_ligand_receptor_inputs.py` | Writes sender and Treg expression tables for NicheNet | Treg input objects, `treg_pseudotime.h5ad` → parquet tables | input of 16 |
| `16_treg_ligand_receptor_prioritization.R` | NicheNet ligand activity for late- vs early-pseudotime Treg genes; composite LR priority; weighting schemes; immune-sender ranking used by the Treg Suppression Score | 15 tables, NicheNet model files → panels, ranking tables | Fig 4H-I, Supp 8G-H |
| `utils/http_user_agent.py` | Helper imported by `14`: sends a browser User-Agent with `requests.get` so decoupler can download its resources | | |

## Run order

`01 → 02 → 03-07 → 08 → 09 → 10, 11`, then `12 → 13, 14, 15 → 16`.

`16` (section `immune_sender_ranking`) runs before `cosmx_cd/10_treg_suppression_score.py`, which takes the antigen-presentation LR module from that ranking. `12` and `14` also feed `treatment_validation/03`.

### Why there is a `.py` and an `.R` file for the ligand-receptor step

The atlas and the Treg object are AnnData files, so `15` (Python) writes the sender and Treg expression tables. NicheNet runs in R, so `16` reads those tables and does the prioritization. They are two consecutive steps.

## Not in this repository

- `16` (section `F4HI_nichenet_lr`) reads two tables in `objects/nichenet_outputs_lossact/` (a lineage table of the late-vs-early gene set and the bin 4 vs bin 0 markers); the code that writes them is not included.
