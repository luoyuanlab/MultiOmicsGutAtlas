# TAURUS anti-TNF cohort

Public scRNA-seq of CD biopsies before and after adalimumab (TAURUS; Thomas et al., 2024): 16 patients (10 in remission, 6 not), 96 biopsies. The question is whether the Treg programs and the macrophage-derived signals found in the CD atlas and CD CosMx data change after anti-TNF, and whether they change differently in patients who reach remission.

Panel letters refer to the figures (Fig 6 = TAURUS anti-TNF; Supp 10).

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_taurus_preprocessing.ipynb` | CD biopsies only; QC, normalisation; attaches the AUCell labels | TAURUS h5ad → `xnorm_counts_for_aucell.h5ad.gz`, `raw_norm_counts_anno.h5ad.gz` | input of 02-05 |
| `02_taurus_aucell_annotation.R` | AUCell annotation with the shared 50-type gene sets | `xnorm_counts_for_aucell.h5ad.gz` → `anno/norm_aucell.csv` | labels used by 01 |
| `03_taurus_treg_programs.py` | Treg pseudotime; cell-level mixed models of the Treg programs; biopsy-level mixed models of the Treg Suppression Score, IL-2/STAT5, TNFR2 non-canonical NF-κB, the co-stimulatory/suppressive receptor module and TNFRSF4; Fig 6A | `raw_norm_counts_anno.h5ad.gz`, score genes (`cosmx_cd/10`), tables of `cd_scrna/12` and `cd_scrna/14` → CSVs, panel | Fig 6A, Supp Table 11 |
| `04_taurus_cytokine_sources.py` | Biopsy-level models of OSM, IL1B and IL6 by coarse cell group; split of the all-cell change into cell-group contributions; abundance contrasts; Fig 6B | `raw_norm_counts_anno.h5ad.gz` → cell cache, `cn10_cytokine_coarse_screen.csv`, `cn10_cytokine_contrib_bootstrap.csv`, panel | Fig 6B, Supp Table 11 |
| `05_taurus_ligand_senders.py` | Senders of the ligands of the four Treg activation programs, Remission vs Non-remission | h5ad, `04` cell cache → CSV, panel | Supp 10 |

Section ids, in run order: `03`: `treg_pseudotime`, `cell_level_lmm`, `F6A_programs_biopsy`, `F6A_suppression_score_biopsy`, `F6A_post_treatment_difference`, `F6A_pathways_biopsy`, `F6A_tnfrsf4_biopsy`, `F6A_treg_state_lines`. `04`: `cytokine_cell_cache`, `F6B_coarse_screen`, `F6B_contribution_se`, `F6B_cytokine_sources`. `05`: `S10_ligand_sources`, `S10_ligand_senders`.

## Run order

`01` up to the `xnorm_counts_for_aucell.h5ad.gz` export → `02` → rest of `01` → `03`, `04` → `05` (`05` reads the cell cache written by `04`).

Before `03`: `cd_scrna/12`, `cd_scrna/14` and `cosmx_cd/10` (Treg signatures and reference tables).

Model sections call `Rscript` (lme4, lmerTest) through temporary R scripts. Use `OPENBLAS_NUM_THREADS=32` for bit-identical model output.
