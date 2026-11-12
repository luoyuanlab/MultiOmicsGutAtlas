# UC CosMx spatial transcriptomics

CosMx data of inflamed UC and control colon: whole-transcriptome (WTX), 6K with custom probes and 1K panels cut from the same blocks, an independent 1K validation set, and additional CosMx 6K sections of UC and control colon. The question is whether the macrophage-fibroblast signals predicted from the scRNA-seq atlas are seen in tissue: do SELENOP⁺MMP9⁺SPP1⁺ macrophages sit next to OGN⁺RSPO3⁺ and activated fibroblasts, and do ligand and target genes co-vary with distance?

Panel letters refer to the figures (Fig 3 = UC spatial; Supp 6 = panels and imputation).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_wtx_processing.ipynb` | WTX QC, normalisation, ComBat; exports counts for AUCell; reads the labels back | AtoMx flat files → `qc_wtx_v7_merged_anno.h5ad` | Fig 3B-C, Supp 6A |
| `02_6k_processing.ipynb` | 6K + custom QC and normalisation; exports counts for AUCell | AtoMx flat files → 6K h5ad | Supp 6D |
| `03_1k_processing_imputation_comparison.ipynb` | 1K QC and normalisation; comparison with scRNA-seq and the other panels; annotation of the GCN-imputed 1K; UMI and expression-range comparison; off-lineage expression of imputed genes | 1K and WTX objects, imputed 1K → plots, h5ad | Supp 6A-C, 6G-H, 6K |
| `04_validation_1k_processing.ipynb` | QC and annotation of the independent 1K validation set; marker dot plots before and after imputation | 1K validation run → h5ad, plots | Supp 6L |
| `05_aucell_annotation.R` | AUCell annotation of every UC panel with the shared 50-type gene sets (1K, imputed 1K, 6K, WTX) | normalised count CSVs → label CSVs | labels used by 01-04, 06-14 |
| `06_annotation_spatial_maps.ipynb` | Annotation maps of WTX, 6K and 1K sections; CD68 / PanCK protein vs RNA labels before and after imputation, Otsu thresholds | annotated h5ads → maps | Fig 3B-C, Supp 6D, 6G-I |
| `07_1k_differential_expression.R` | DE between inflamed rectal UC and other tissue in the CosMx 1K cohort (16 samples) | 1K Seurat objects → DE table | Fig 3A |
| `08_il1_cellchat_distance_sensitivity.R` | CellChat IL-1 signalling on the imputed 1K annotation; contact and range settings | imputed 1K labels → CellChat objects, plots | Supp 6M |
| `09_mac_fib_spatial_input.py` | Shared spatial input (CosMx 6K, middle section of each block) and the per-macrophage distance and receptor table | 6K h5ads → `uc_cosmx_colocal_input.h5ad`, distance table | input of 11-14 |
| `10_wtx_ligand_regulator_pairs.ipynb` | NicheNet-predicted regulators in OGN⁺RSPO3⁺ fibroblasts near vs far from Mac S+M+S+ (WTX), by expression and target score | WTX h5ad, NicheNet tables of `uc_scrna/14` → `ogn_spp1_mac_nichenet_tf_spatial_validation_combined_v2.csv` | input of 11 |
| `11_ligand_regulator_cooccurrence.py` | Ligand in the nearest Mac S+M+S+ vs regulator in each OGN⁺RSPO3⁺ fibroblast (random-effects pooled); co-occurrence of activated fibroblasts; neighbour models; NicheNet ligand activity plot | 09 input, 10 pair table, 1K validation, WTX → panels, supplementary tables | Fig 3D, 3E, 3I, 3J |
| `12_tnc_mmp1_spatial.py` | TNC in OGN⁺RSPO3⁺ fibroblasts near MMP1⁺ macrophages (logistic regression, maps); TNC vs its NicheNet regulators against random genes | 09 input, WTX, 11 per-piece table → panels, CSVs | Fig 3F, Supp 6E |
| `13_macrophage_distance_glmm.R` | Mixed-effects logistic regression of SPP1, MMP9 and receptor detection in SELENOP⁺ macrophages on distance to OGN⁺RSPO3⁺ fibroblasts | 09 distance table → `glmm_final.csv` | Fig 3H |
| `14_macrophage_distance_forest.py` | Forest plot and supplementary table of the models | 13 model table → `fig_distance_forest.png`, `supp_table_distance_models_cosmx.csv` | Fig 3H |

## Run order

`01, 02, 03, 04` up to their count exports → `05` → the rest of `01-04` → `06 → 08`. `07` is independent. Then `uc_scrna/14 → 10`, `09 → 11 → 12`, and `09 → 13 → 14`.

GCN imputation of the 1K panels (input of `03`, `04`, `06`, `08`, `11`) is run with `imputation_model/`.

## Not in this repository

- The QC that built the 16 CosMx 1K Seurat objects read by `07`.
