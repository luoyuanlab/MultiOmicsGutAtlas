# CITE-seq and TCR

CITE-seq (RNA + surface protein) with TCR sequencing of intestinal biopsies (GSE311585). It is used for three things: protein markers of SELENOP⁺ vs SELENOP⁻ macrophages, protein markers of CD9⁺ vs CD9⁻ macrophages, and the TCR repertoire of early- vs late-pseudotime Tregs.

## Files

| File | What it does | Inputs → outputs | Panels / tables |
|---|---|---|---|
| `01_aucell_annotation.R` | Builds the RNA + ADT object; AUCell annotation of myeloid cells and of T cells | count matrices, normalised exports of 02 and 04 → label CSVs | Supp 4A-B, Supp 7F, Supp 8I (Treg set) |
| `02_myeloid_protein_selenop.ipynb` | Quality filter, CLR protein normalisation, panel-wide effect sizes; SELENOP⁺ vs SELENOP⁻ macrophage protein comparison | count matrices, 01 myeloid labels → plots, tables | Supp 4A-B |
| `03_myeloid_protein_cd9.ipynb` | CD9⁺ vs CD9⁻ macrophages (CD9 transcript detected or not): Mann-Whitney U with BH and Cohen's d for every protein; heatmap and split violins of the top proteins | count matrices, 01 myeloid labels → plots | Supp 7F |
| `04_treg_subset_export.ipynb` | Normalises T/NK cells with TCR contigs for `01`; then subsets the AUCell Tregs and exports their raw counts for the label transfer | count matrices, 01 T cell labels → `good_anno_t_wtcr_norm_counts.csv`, `notstrict_treg_raw_counts.csv` | input of `cosmx_cd/09` |
| `05_treg_tcr_repertoire.py` | Clonotypes per lane, diversity, clonal sharing and TRBV usage of early (bins 0-1) vs late (bins 2-4) Tregs | TCR contigs, labels from `cosmx_cd/09` → panel, table | Supp 8I, Supp Table 9 |

## Run order

`01` section 1 → `02` (myeloid normalisation export) → `01` section 2 → rest of `02`, and `03`.
`04` (T/NK export) → `01` section 3 → rest of `04` → `cosmx_cd/09_treg_label_transfer.R` → `05`.
