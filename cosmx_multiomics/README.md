# CosMx spatial proteomics: UC vs CD

CosMx RNA + protein (64-antibody panel) on one UC rectum section (15 FOVs, 13,415 cells) and one inflamed CD ileum section (17 FOVs, 20,419 cells). The figure uses the protein data only, as a parallel check of the RNA-based findings: are myeloid cells more inflammatory in UC and more antigen-presenting in CD, and do Tregs in CD carry more co-stimulatory receptors? Disease and site are confounded, and the unit of analysis is the FOV.

Panel letters refer to the figures (Fig 7 = CosMx spatial proteomics).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_uc_rna_aucell_annotation.R` | AUCell RNA labels of the UC sections (sections `Cos4_UC_M_batch1_S2` = 385R_381R and `batch2` = 322R_367R) | normalised RNA CSVs of 03 → label CSVs | input of 03 |
| `02_cd_rna_aucell_annotation.R` | AUCell RNA labels of the CD section; `04` uses them to drop cycling cells | normalised RNA CSV of 04 → `aucell_anno_upd.csv` | input of 04 |
| `03_uc_rna_protein_processing.ipynb` | UC sections: RNA QC and normalisation, AUCell labels, protein matched to cells (per-protein centring and per-cell geometric-mean centring) | AtoMx objects, 01 labels → `RNA1_Cos4_385R_381R_S2_norm_wpr_anno_jd.h5ad` | input of 05 |
| `04_cd_rna_protein_processing.ipynb` | CD section: RNA QC and normalisation, protein matched to cells and centred on each cell's geometric mean | AtoMx flat files, 02 labels → `rna_norm_pr_geonorm_upd_anno.h5ad` | input of 05 |
| `05_protein_cell_types_myeloid_treg.py` | Protein gating of cell types and maps; composition by distance from epithelium; HLA-DR / iNOS myeloid maps; myeloid protein effect sizes (CD - UC); ICOS, GITR, CD27 on Tregs vs conventional CD4⁺ T cells | 03 and 04 objects → figure, `fig7_cell_types_depth_profiles.csv`, `fig7_myeloid_treg_markers.csv` | Fig 7A-E, Supp Table 12 (gates, markers) |

## Run order

`03` up to its count export → `01` → rest of `03`; `04` up to its count export → `02` → rest of `04`; then `05`.

The figure is drawn by one script (`05`) because all five parts share the same gates and layout.

## Not in this repository

- The AUCell label CSVs of the UC sections are no longer on disk; the processed UC object written by `03` is.
- The whole-panel protein marker screen behind the choice of supporting markers (Supp Table 12).
