# CD CosMx spatial transcriptomics: Treg suppression niches

CosMx data of Crohn's disease: WTX ileal biopsies and 6K + custom resections (11 sections from 4 patients, 6,554 strict Tregs). Treg pseudotime labels and a Treg Suppression Score are carried over from the scRNA-seq atlas (`cd_scrna/`). The question is where in tissue Tregs are more or less suppressive: which cellular neighbourhoods, which nearby antigen-presenting cells, and how this changes in CXCL12-High Regions next to OGN⁺RSPO3⁺ fibroblasts.

Panel letters refer to the figures (Fig 5 = CD CosMx niches; Supp 9).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_wtx_biopsy_processing.ipynb` | WTX biopsies (two slides, four patients): sections, QC, normalisation, ComBat embedding; exports counts for AUCell and maps the AUCell labels; adds patient and condition | AtoMx flat files, AUCell labels → `norm_anno.h5ad` | input of 02 |
| `02_wtx_biopsy_treg_labels_liana.ipynb` | Exports the raw counts of all biopsy Tregs; flags core Tregs; writes the biopsy object with Treg labels; LIANA bivariate ligand-receptor scores of the inflamed slide (NicheNet network, Moran's I) | `norm_anno.h5ad` → `all_treg_raw_counts.csv`, `norm_anno_tregbin.h5ad`, `lr_liana_i0476_1_nn_db.h5ad` | input of 09-11 |
| `03_resection_processing.ipynb` | Resection QC, normalisation and tissue-section assignment (12 runs) | AtoMx flat files → per-run `Processed/norm_log1p.h5ad`, count CSVs | input of 04-06, 11 |
| `04_resection_aucell_annotation.R` | AUCell annotation of every resection run | 03 counts → `aucell_anno.csv` | resection labels |
| `05_resection_treg_counts.ipynb` | Raw counts of AUCell Tregs, then strict Tregs (core Treg score ≥ 0, cell-cycle score < 10; 6,554 cells) | 03 objects, 04 labels → `raw_counts_tregs.csv`, `raw_counts_strict_tregs.csv` | input of 09-11 |
| `06_resection_liana.ipynb` | Per-FOV LIANA bivariate ligand-receptor scores (Gaussian kernel, 100 permutations) | 03 objects, 04 labels → per-FOV h5ads | input of 07 |
| `07_resection_liana_merge.ipynb` | Merges the per-FOV LIANA output of all sections | 06 output → `combined/spot_lr_df.parquet` | input of 11, 13, 17 |
| `08_treg_reference_export.py` | Writes the scRNA-seq Treg reference (raw counts, pseudotime bins) for label transfer | `treg_pseudotime.h5ad` → reference CSVs | input of 09 |
| `09_treg_label_transfer.R` | Seurat label transfer of Treg pseudotime bins onto WTX biopsy, resection and CITE-seq Tregs | 08 reference, Treg count CSVs of 02, 05 and `cite-tcr/04` → label tables | input of 10, `cite-tcr/05` |
| `10_treg_suppression_score.py` | Antigen-presentation Treg LR module; Treg Suppression Score per CosMx Treg (late- vs early-pseudotime signature, sign reversed, z within section) | Treg object, NicheNet ranking, 09 labels, 02 objects, 05 counts → per-cell score table, signature genes | Fig 5A-F, Supp 9, Fig 6A, Supp Table 10 |
| `11_neighbourhoods_local_tables.py` | Cellular neighbourhoods (20 nearest cells, k = 10); per-Treg local ligand, receptor and cell-type tables; local CXCL12 from expected stromal producers and its models | cell metadata, LIANA scores, 10 → CN and per-Treg tables | input of 12-17 |
| `12_neighbourhood_treg_score.py` | CN composition table; Treg Suppression Score inside vs outside each CN | 10-11 tables → CN tables, forest plot | Fig 5C (left), inputs of Supp 9A, 9C |
| `13_neighbourhood_signalling.py` | Ligand-receptor and cytokine signalling per CN (CellChatDB modules, inflammatory axes); CN characterisation; combined score and signalling figure; per-niche error bars | LIANA scores, 10-12 → panels, module tables | Fig 5C, Supp 9A-B, Supp Table 10 |
| `14_treg_score_maps.py` | Treg Suppression Score maps of the WTX biopsies with zoom windows; resection section coloured by CN with high- and low-score windows | WTX polygons, 10-11 → panels | Fig 5A-B |
| `15_apc_local_models.py` | Score vs local fraction of each cell type; joint APC-type model; DC effect by OGN⁺RSPO3⁺ tertile; summary and interaction tables | 11 per-Treg tables, 12 → panels, CSVs | Fig 5D, Supp 9D, numbers in the Results |
| `16_cxcl12_regions.py` | CXCL12-High Regions and APC-rich areas, section maps; score by APC quintile inside and outside the regions; CXCL12 by cell type; APC abundance per CN | 10-12, 15 → panels, CSVs | Fig 5E-F, Supp 9E-F |
| `17_neighbourhood_apc_mhcii.py` | Antigen-presenting cells and MHC-II per CN, and the APC-matched inside-minus-outside score | 11-13, 16 → panel | Supp 9C |

## Run order

`01` up to its count export → AUCell → rest of `01` → `02`. `03 → 04 → 05 → 06 → 07`. Then `08 → 09 → 10 → 11 → 12 → 13 → 14 → 15 → 16 → 17`.

Before `08-10`: `cd_scrna/12` (Treg object) and `cd_scrna/16` section `immune_sender_ranking` (antigen-presentation LR module).

`17` reads the APC table that `16` writes with `RUN_FIG2=1 RUN_SECTIONS=S9E_cxcl12_by_celltype`. In `16`, section `F5EF_section_maps` runs once per figure: `T94_SEC=43563_A4 T94_MAKE=row` (Fig 5E) and `T94_SEC=43563_A8 T94_MAKE=piece` (Fig 5F). Fig 5C is drawn with `VM=1.2 WIDE=1 RUN_SECTIONS=F5C_score_signalling python 13_neighbourhood_signalling.py`.

### Why the Treg step has a `.py` and an `.R` file

They are three consecutive steps. The Treg pseudotime object is an AnnData file, so `08` (Python) writes the reference as CSV. Seurat label transfer runs in R (`09`). The Treg Suppression Score and everything after it are in Python (`10`).

## Not in this repository

- The AUCell run that labels the WTX biopsies (`anno/norm_counts_aucell.csv`, read by `01`).
- The resection cell metadata table `combined/all_samples.parquet`.
- `combined/strict_treg_pt_transfer_labels.csv`, from which `07` takes the strict-Treg cell IDs (the same 6,554 cells as the label table of `09`).
- The CellChatDB snapshot `objects/cellchatdb_human_interactions_snapshot.csv` (read by `13` and `treatment_validation/05`) and the per-section table `ligand_receptor_state/c_top_apc_quintile_per_section_exp.csv` (read by `16`).
