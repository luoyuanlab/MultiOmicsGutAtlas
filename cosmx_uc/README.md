# UC CosMx spatial transcriptomics

CosMx data of inflamed UC and control colon: whole-transcriptome (WTX), 6K with custom probes and 1K panels cut from the same blocks, an independent 1K validation set, and additional CosMx 6K sections of UC and control colon. The question is whether the macrophage-fibroblast signals predicted from the scRNA-seq atlas are seen in tissue: do SELENOP⁺MMP9⁺SPP1⁺ macrophages sit next to OGN⁺RSPO3⁺ and activated fibroblasts, and do ligand and target genes co-vary with distance?

Panel letters refer to the figures (Fig 3 = UC spatial; Supp 6 = panels and imputation).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_1k_processing_imputation_comparison.ipynb` | WTX QC, normalisation and AUCell labels; 1K QC and normalisation; genes per cell vs scRNA-seq; UMIs per detected gene and per-gene expression range across WTX, 6K and 1K; off-lineage expression of GCN-imputed lineage genes | AtoMx WTX flat files, 1K objects, imputed 1K, UC 6K h5ad → `qc_wtx_v7_merged_anno.h5ad`, plots | Fig 3B-C (WTX object), Supp 6A-C, 6K |
| `02_validation_1k_processing.ipynb` | QC and annotation of the independent 1K validation set; marker dot plots before and after imputation | 1K validation run → h5ad, plots | Supp 6L |
| `03_aucell_annotation.R` | AUCell annotation of every UC panel with the shared 50-type gene sets (1K, imputed 1K, 6K, WTX) | normalised count CSVs → label CSVs | labels used by 01, 02, 04-12 |
| `04_annotation_spatial_maps.ipynb` | Annotation maps of WTX, 6K and 1K sections; CD68 / PanCK protein vs RNA labels before and after imputation, Otsu thresholds | annotated h5ads, 1K label tables → `imp_ywl_anno.h5ad` (1K-only and GCN-imputed labels), `psuedo_1k.csv` (input of 03), maps | Fig 3B-C, Supp 6D, 6G-I |
| `05_1k_differential_expression.R` | DE between inflamed rectal UC and other tissue in the CosMx 1K cohort (16 samples) | 1K Seurat objects → DE table | Fig 3A |
| `06_il1_cellchat_distance_sensitivity.R` | CellChat IL-1 signalling on the imputed 1K annotation; contact and range settings | imputed 1K labels → CellChat objects, plots | Supp 6M |
| `07_mac_fib_spatial_input.py` | Shared spatial input (CosMx 6K, middle section of each block) and the per-macrophage distance and receptor table | 6K h5ads → `uc_cosmx_colocal_input.h5ad`, `cohort_summary.csv`, `ligand_programs.csv`, distance table | input of 09-12 |
| `08_wtx_ligand_regulator_pairs.ipynb` | NicheNet-predicted regulators in OGN⁺RSPO3⁺ fibroblasts near vs far from Mac S+M+S+ (WTX), by expression and target score | WTX h5ad, NicheNet tables of `uc_scrna/14` → `ogn_spp1_mac_nichenet_tf_spatial_validation_combined_v2.csv` | input of 09 |
| `09_ligand_regulator_cooccurrence.py` | Ligand in the nearest Mac S+M+S+ vs regulator in each OGN⁺RSPO3⁺ fibroblast (random-effects pooled); co-occurrence of activated fibroblasts; neighbour models; NicheNet ligand activity plot | 07 input, 08 pair table, 1K validation, WTX → panels, supplementary tables | Fig 3D, 3E, 3I, 3J |
| `10_tnc_mmp1_spatial.py` | TNC in OGN⁺RSPO3⁺ fibroblasts near MMP1⁺ macrophages (logistic regression, maps); TNC vs its NicheNet regulators against random genes | 07 input, WTX → panels, CSVs | Fig 3F, Supp 6E |
| `11_macrophage_distance_glmm.R` | Mixed-effects logistic regression of SPP1, MMP9 and receptor detection in SELENOP⁺ macrophages on distance to OGN⁺RSPO3⁺ fibroblasts | 07 distance table → `glmm_final.csv` | Fig 3H |
| `12_macrophage_distance_forest.py` | Forest plot and supplementary table of the models | 11 model table → `fig_distance_forest.png`, `supp_table_distance_models_cosmx.csv` | Fig 3H |

## Run order

`01, 02` up to their count exports → `03` → the rest of `01-02` → `04 → 06`. `05` is independent. Then `uc_scrna/14 → 08`, `07 → 09 → 10`, and `07 → 11 → 12`.

GCN imputation of the 1K panels (input of `01`, `02`, `04`, `06`, `09`) is run with `imputation_model/`.

## Not in this repository

- The QC that built the 16 CosMx 1K Seurat objects read by `05`.
- The QC and normalisation of the UC CosMx 6K + custom sections cut from the WTX/1K blocks (`uc_1k_6k_wtx/6k/Processed_merged_upd` and the per-batch `_upd` folders: `qc_6k_merged_log1p_*genes.h5ad`, `count_csv/norm_counts_6175genes.csv`, `qc_6k_*_log1p_anno_6375genes.h5ad`). `03` annotates its count export, and `01`, `04` and `07` read the objects (Supp 6B-D).
