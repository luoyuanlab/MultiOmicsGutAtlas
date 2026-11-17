# CosMx spatial proteomics (RNA + protein in the same cells)

CosMx RNA + protein (64-antibody panel) on one UC rectum section (15 FOVs, 13,415 cells) and one inflamed CD ileum section (17 FOVs, 20,419 cells). The protein data are a parallel check of the RNA-based findings: does each RNA-defined lineage carry its canonical protein, what do protein-defined lineages look like across each section, and do PD-L1⁺ macrophages sit near PD-1⁺ Tregs in CD ileum (protein-level check of the CosMx 6K ligand-receptor proximity in `cosmx_cd/18-19`)? The unit of analysis is the FOV.

Panel letters refer to the figures (Supp 6D-E = UC section; Supp 9I-J = CD section; Supp 11D = PD-L1 / PD-1 map).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_uc_rna_aucell_annotation.R` | AUCell RNA labels of the UC section 385R_381R (`Cos4_UC_M_batch1_S2`) | normalised RNA CSV of 03 → label CSV | input of 03 |
| `02_cd_rna_aucell_annotation.R` | AUCell RNA labels of the CD section; `04` uses them to drop cycling cells | normalised RNA CSV of 04 → `aucell_anno_upd.csv` | input of 04 (the same code is re-run on the matched RNA by 06) |
| `03_uc_rna_protein_processing.ipynb` | UC section 385R_381R: RNA QC and normalisation, AUCell labels, protein matched to cells (per-protein centring and per-cell geometric-mean centring) | AtoMx objects, 01 labels → `RNA1_Cos4_385R_381R_S2_norm_wpr_anno_jd.h5ad` | input of 07 |
| `04_cd_rna_protein_processing.ipynb` | CD section: RNA QC and normalisation, protein matched to cells and centred on each cell's geometric mean | AtoMx flat files, 02 labels → `rna_norm_pr_geonorm_upd_anno.h5ad` | input of 06 |
| `06_cd_rna_realign.py` | Matches the CD RNA counts to the right cells (the 04 object filled RNA rows by position, and AtoMx sorts the exprMat by FOV as text); recomputes QC and the `rna_qc_pass` flag; re-runs the AUCell annotation of `02` on the matched RNA; validation table | 04 object, AtoMx flat files, shared gene sets → `rna_norm_pr_geonorm_anno_aligned.h5ad`, `cd_multiomics_rna_realign_validation.csv` | input of 07 |
| `07_protein_lineage_maps_lr_proximity.py` | Protein gates (within-section percentiles); PD-L1⁺ macrophage → PD-1⁺ Treg proximity per FOV (30 µm; 50 µm sensitivity; permutation null within FOV, stratified pooled P); RNA lineage vs canonical protein per FOV (Wilcoxon, BH over 12 pairs per section); protein-defined lineage maps; PD-L1⁺ Mac / PD-1⁺ Treg map | 03 object, 06 object → `protein_lr_proximity.csv`, combined figure and CSV, one PNG per panel | Supp 6D-E, Supp 9I-J, Supp 11D |

No file is numbered 05. Section ids of `07`, in run order: `S11D_protein_lr_proximity`, `S6DE_S9IJ_S11D_figure`.

## Run order

`03` up to its count export → `01` → rest of `03`; `04` up to its count export → `02` → rest of `04` → `06`; then `07`.

## Not in this repository

- The AUCell label CSV of the UC section is no longer on disk; the processed UC object written by `03` is.
- `fig5_spatial_maps_v7.csv`, the per-gate cell counts of the earlier protein figure; `07` copies them into the `n_gate_v7` column of its CSV.
- The whole-panel protein marker screen behind the choice of supporting markers (Supp Table 12).
