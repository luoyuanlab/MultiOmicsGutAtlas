# Imputation benchmark

Benchmark of six gene-imputation methods on held-out CosMx WTX data with a pseudo-6K design: cells from held-out FOVs are reduced to the 6K gene panel, each method imputes the remaining genes, and the result is compared with the measured WTX expression. The GCN model itself is in `imputation_model/`.

Panel letters refer to Supp 6 (imputation) and Supp Table 7.

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_gcn_pseudo6k_postprocess.ipynb` | Post-processing of the GCN pseudo-6K run (test-cell split by FOV, merge of observed and imputed genes) | WTX h5ad, GCN output from `imputation_model/Tutorial.ipynb` → combined matrix | Supp 6J, Supp Table 7 |
| `02_tangram_imputation.ipynb` | Tangram imputation | WTX h5ad, 6K gene list, scRNA-seq reference → combined matrix | Supp 6J, Supp Table 7 |
| `03_gimvi_imputation.ipynb` | gimVI imputation | as above | Supp 6J, Supp Table 7 |
| `04_stplus_imputation.ipynb` | stPlus imputation | as above | Supp 6J, Supp Table 7 |
| `05_sprite_imputation.ipynb` | SPRITE imputation | as above | Supp 6J, Supp Table 7 |
| `06_cellplm_imputation.ipynb` | CellPLM imputation | as above | Supp 6J, Supp Table 7 |
| `07_imputation_comparison.ipynb` | Compares all methods: cell-type classification (Hamming loss, micro F1), clustering metrics, per-gene Pearson and Spearman correlation | outputs of 01-06 → heatmaps, tables | Supp 6J, Supp Table 7 |

## Run order

`01-06` in any order (each after its method has run), then `07`.
