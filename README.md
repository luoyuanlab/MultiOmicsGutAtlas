# MultiOmics Gut Atlas

Code accompanying **"An integrated single-cell and spatial atlas identifies OGN⁺RSPO3⁺ fibroblasts as a conserved disease mechanism across ulcerative colitis and Crohn's disease."**

This repository contains the analysis code used to construct and analyze an integrated single-cell RNA-seq, CITE-seq, and spatial transcriptomic/multi-omic atlas of ulcerative colitis (UC) and Crohn's disease (CD), including a graph neural network (GCN)-based imputation framework for lower-plex spatial transcriptomics panels.

## Overview

- Integration and annotation of public and newly generated scRNA-seq/CITE-seq datasets for UC and CD
- Cell-cell interaction, pseudotime trajectory, and ligand-receptor prioritization analyses
- CosMx spatial transcriptomics (1K, 6K, WTX) and spatial multi-omics (RNA + protein) processing and annotation
- Graph neural network-based gene expression imputation for lower-plex spatial panels, benchmarked against Tangram, gimVI, stPlus, SPRITE, and CellPLM
- Spatial cell-cell interaction, co-occurrence, and neighborhood enrichment analyses

## Repository Structure

```
MultiOmicsGutAtlas/
├── uc_scrna/               # UC single-cell RNA-seq analysis
├── cd_scrna/               # CD single-cell RNA-seq analysis
├── cosmx_uc/               # UC CosMx spatial transcriptomics
├── cosmx_cd/               # CD CosMx spatial transcriptomics
├── cosmx_multiomics/       # Multiomics (RNA + protein) spatial integration
├── treatment_validation/   # SAHA and TAURUS treatment cohort validation
├── imputation_benchmark/   # Benchmarking spatial gene imputation methods
├── imputation_model/       # GCN-based gene imputation model (CosMxGAE)
└── cite-tcr/               # CITE-seq and TCR clonotype analysis
```

### `uc_scrna/` — UC Single-Cell Atlas
Integration and annotation of UC single-cell RNA-seq data across multiple patients.

| Script | Description |
|--------|-------------|
| `00` | Dataset preparation |
| `01` | Atlas integration |
| `02–04` | Myeloid, fibroblast, and other lineage subclustering |
| `05` | Atlas assembly |
| `06–07` | Integration QC and UMAPs |
| `08–09` | Cell type abundance (Kruskal-Wallis, scCODA) |
| `10` | GSEA |
| `11` | Macrophage potency (CytoTRACE2) and pseudotime |
| `12` | Cell-cell interaction (CellChat) |
| `13` | Ligand-receptor inference (NicheNet) |
| `14` | Fibroblast pseudotime |

### `cd_scrna/` — CD Single-Cell Atlas
Parallel pipeline for CD single-cell RNA-seq, including T cell and epithelial subclustering.

| Script | Description |
|--------|-------------|
| `00–07` | Dataset preparation, integration, subclustering (myeloid, fibroblast, T cell, other, epithelial), atlas assembly |
| `08` | NicheNet Treg analysis |
| `09–10` | Atlas cleaning and UMAPs |
| `11` | scCODA compositional analysis |
| `12` | Treg pseudotime |

### `cosmx_uc/` — UC CosMx Spatial Transcriptomics
Processing, annotation, and spatial analysis of UC CosMx data across whole-transcriptome (WTX), 6K, and 1K panels.

| Script | Description |
|--------|-------------|
| `01` | WTX panel processing |
| `02` | 6K panel processing |
| `03` | 1K panel QC, normalization, annotation, and imputation comparison |
| `04` | AUCell cell type annotation (R) |
| `05` | Spatial plots, protein annotation, co-occurrence |
| `06` | CellChat CCI and distance sensitivity analysis (R) |
| `07` | Validation 1K processing |
| `08` | Fibroblast–macrophage colocalization (spatial mixed model) |
| `09` | OGN+RSPO3 fibroblast–macrophage spatial analysis (NicheNet + TNC) |
| `10` | Cross-disease fibroblast neighborhood comparison (UC vs CD) |

### `cosmx_cd/` — CD CosMx Spatial Transcriptomics
Spatial analysis of CD biopsies and resection samples.

| Script | Description |
|--------|-------------|
| `01` | Biopsy AUCell annotation (R) |
| `02` | Biopsy LIANA + NMF analysis |
| `03` | Resection sample processing (12 samples) |
| `04` | Resection AUCell annotation (R) |
| `05` | Per-FOV LIANA ligand-receptor inference |
| `06` | LR module analysis and LOWESS scatter |
| `07` | Treg proximity analysis (early vs late disease) |
| `08` | Treg label transfer — biopsy and resection (R) |

### `cosmx_multiomics/` — Multiomics Spatial Integration
Integration of RNA and protein data from CosMx multiomics panels in UC and CD.

| Script | Description |
|--------|-------------|
| `01` | UC multiomics AUCell annotation (R) |
| `02` | CD multiomics AUCell annotation (R) |
| `03` | Spatial RNA + protein integration and visualization |

### `treatment_validation/` — Treatment Cohort Validation
Spatial and single-cell validation of therapeutic interventions.

| Script | Description |
|--------|-------------|
| `01` | SAHA CosMx annotation (AUCell + Treg label transfer, R) |
| `02` | SAHA spatial analysis: Treg bins and fibroblast neighborhoods, pre/post-treatment |
| `03` | TAURUS CD scRNA annotation (AUCell + Treg label transfer, R) |
| `04` | TAURUS anti-TNF scRNA: Treg pseudotime delta pre/post-treatment |

### `imputation_benchmark/` — Gene Imputation Benchmarking
Evaluation of six spatial gene imputation methods on held-out CosMx WTX data using a pseudo-6K panel design.

**Methods benchmarked**: GCN (CosMxGAE), Tangram, gimVI, stPlus, SPRITE, CellPLM

**Metrics**: Cell-type classification (Hamming loss, Micro F1), clustering quality (Silhouette, Calinski-Harabasz, Davies-Bouldin), per-gene correlation (Pearson, Spearman)

| Script | Description |
|--------|-------------|
| `00–01` | GCN pseudo-6K post-processing |
| `02–06` | Per-method imputation (Tangram, gimVI, stPlus, SPRITE, CellPLM) |
| `07` | Benchmark comparison across all methods |

### `imputation_model/` — CosMxGAE Gene Imputation Model
A graph autoencoder (GAE) model for imputing gene expression from limited CosMx panels using a whole-transcriptome reference. The model uses mutual nearest neighbor (MNN) graphs to bridge reference (WTX) and query (panel) cells and trains to predict held-out gene expression via weighted MSE loss.

- `Tutorial.ipynb` — end-to-end usage example
- `Imputation_model/` — model source code (`GAE.py`, `adj.py`, `normalization.py`, `utils.py`)

### `cite-tcr/` — CITE-seq and TCR Analysis
Protein-level and clonotype analysis using CITE-seq and paired TCR sequencing.

| Script | Description |
|--------|-------------|
| `01` | CITE-seq annotation (R) |
| `02` | Macrophage CITE-seq analysis |
| `03` | Treg CITE-seq and TCR clonotype analysis |

---

## Requirements

**R (≥4.x):** Seurat (v5.1.0), CellChat (v2.2.0), NicheNet (v2.2.1), AUCell (v1.30.1), CytoTRACE2 (v1.1.0)

**Python (≥3.9):** Scanpy (v1.10.2), Squidpy (v1.6.2), decoupler (v2.1.1), Palantir (v1.4.2), LIANA (v1.6.1), pertpy (v0.10.0), scvi-tools (v1.2.0)

---

## Data

Newly generated sequencing and spatial data are deposited in the Gene Expression Omnibus (GEO):

- CITE-seq/TCR-seq: [GSE311585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311585)
- CosMx UC and CD 1K cohort: [GSE312415](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312415)
- CosMx WTX/6K/resection/multi-omics data: accession pending, will be added upon release

Public datasets used for atlas integration are cited in the manuscript's Supplementary Tables 1 and 8.

---

## Citation

If you use this code or data, please cite:

> [Full citation — add once accepted/assigned a DOI]

---

## License

> [Add a license (e.g., MIT, Apache 2.0) — currently unspecified]

---

## Contact

Questions can be directed to the corresponding authors:
Yuan Luo (yuan.luo@northwestern.edu) and Parambir S. Dulai (parambir.dulai@northwestern.edu).
