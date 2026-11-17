#!/usr/bin/env python
"""CD spatial multi-omics (run Cos4_315i_375i): RNA matched to the right cells, and the CD RNA labels redone on the matched RNA.

Panels: input of 07 (Supp 9I-J, Supp 11D); no panel of its own. Adapted from revision script cd_multiomics_rna_realign.py.
Why: in the object written by 04, obs (AtoMx metadata, coordinates, FOV, '_geo' protein) is consistent, but X / layers['raw'] were
  filled by row position from the RNA exprMat, which AtoMx sorts by FOV as text (10, 11, ..., 1, ...) while the metadata is sorted by
  FOV as a number, so every RNA row belonged to another cell.
Steps: 1. match the RNA exprMat (Updated/FlatFiles) to obs on the unique cell key (fov, cell_ID); 2. layers['raw'] = matched counts
  (same 6,195 genes), X = sc.pp.normalize_total (no log), as in 04; 3. scanpy QC columns from the matched counts and the 04 thresholds
  as the flag 'rna_qc_pass' (200 <= counts <= 3500, 100 <= genes <= 2000, Area <= 30000); 4. AUCell annotation exactly as
  02_cd_rna_aucell_annotation.R (default parameters, shared gsc gene sets, top-AUC label, set.seed(123)) on the matched normalised RNA,
  labels mapped with cell_type_category_map.csv as in 04; 5. obs_names = cell_id; 6. validation table (RNA row sums vs nCount_RNA,
  EPCAM vs EpCAM and PTPRC vs CD45 protein, before and after).
Inputs: rna_norm_pr_geonorm_upd_anno.h5ad (04); AtoMx flat files 315i_375i_exprMat_file.csv.gz, 315i_375i_metadata_file.csv.gz;
  shared gene sets cell_type_markers_gene_sets.RData; cell_type_category_map.csv.
Outputs: rna_norm_pr_geonorm_anno_aligned.h5ad (read by 07; the input object is not modified);
  SuppFig6D-E_9I-J_11D_spatial_proteomics/analysis/cd_multiomics_rna_realign_validation.csv.
Run: python 06_cd_rna_realign.py (calls Rscript for AUCell; R: Seurat, AUCell, GSEABase, data.table).
"""
import os
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse, stats

CD = '/share/fsmresfiles/UC/AtoMx/CD_multiOmics/Cos4_315i_375i'
FF = f'{CD}/Cos4_315i_375i_RNA_1_7.26/Updated/FlatFiles'
IN_H5 = f'{CD}/rna_norm_pr_geonorm_upd_anno.h5ad'
OUT_H5 = f'{CD}/rna_norm_pr_geonorm_anno_aligned.h5ad'
GSC = '/share/fsmresfiles/UC/scRNA-seq/merged_cd/cell_type_markers_gene_sets.RData'
MAPPER = '/share/fsmresfiles/UC/scRNA-seq/merged_cd/cell_type_category_map.csv'
TMP = '/data/ydn4687/tmp_cd_multiomics_realign'
OUT = '/share/fsmresfiles/UC/revision_jenny/SuppFig6D-E_9I-J_11D_spatial_proteomics/analysis'
os.makedirs(TMP, exist_ok=True)

# ── 1. input object and raw exports ──
a = sc.read_h5ad(IN_H5)
obs = a.obs.copy()
obs.index = obs['cell_id'].astype(str).values
genes = list(a.var_names)
key_obs = pd.MultiIndex.from_arrays([obs['fov'].astype(int).values, obs['cell_ID'].astype(int).values])

print('reading exprMat ...', flush=True)
expr = pd.read_csv(f'{FF}/315i_375i_exprMat_file.csv.gz', engine='pyarrow')
meta = pd.read_csv(f'{FF}/315i_375i_metadata_file.csv.gz')
expr = expr.set_index(['fov', 'cell_ID'])
meta = meta.set_index(['fov', 'cell_ID'], drop=False)

# root-cause check: the input raw row of a cell = the exprMat row at that cell's metadata row position
pos_meta = pd.Series(np.arange(len(meta)), index=meta.index).loc[key_obs].values
old_raw = a.layers['raw']
old_sum = np.asarray(old_raw.sum(1)).ravel()
e_sum = expr[genes].sum(1).values
rc = dict(frac_old_rows_equal_exprMat_row_at_metadata_position=float(np.mean(old_sum == e_sum[pos_meta])),
          exprMat_order_equals_metadata_order=bool((expr.index == meta.index).all()))
print(rc)

# ── 2. matched counts ──
X_raw = sparse.csr_matrix(expr.loc[key_obs, genes].values.astype(np.float32))
del expr
b = ad.AnnData(X=X_raw.copy(), obs=obs, var=pd.DataFrame(index=genes), obsm={k: a.obsm[k] for k in a.obsm})
b.layers['raw'] = X_raw

# ── 3. QC columns from matched counts; 04 thresholds as a flag ──
qc_cols = ['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts',
           'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes',
           'pct_counts_in_top_500_genes', 'n_counts', 'n_genes']
b.obs = b.obs.drop(columns=[c for c in qc_cols if c in b.obs])
sc.pp.calculate_qc_metrics(b, inplace=True, log1p=True)
b.obs['n_counts'] = b.obs['total_counts'].values
b.obs['n_genes'] = b.obs['n_genes_by_counts'].values
b.obs['rna_qc_pass'] = ((b.obs.total_counts >= 200) & (b.obs.total_counts <= 3500) & (b.obs.n_genes_by_counts >= 100)
                        & (b.obs.n_genes_by_counts <= 2000) & (b.obs['Area'] <= 30000)).values
sc.pp.normalize_total(b, inplace=True)

# ── 4. AUCell annotation (02_cd_rna_aucell_annotation.R, unchanged) on the matched normalised RNA ──
norm_csv, anno_csv = f'{TMP}/aligned_norm_nolog_rna.csv', f'{TMP}/aucell_anno_aligned.csv'
if not os.path.exists(anno_csv):
    print('writing normalised RNA for AUCell ...', flush=True)
    import pyarrow as pa, pyarrow.csv as pcsv
    Xd = b.X.toarray()
    tbl = pa.table({'': pa.array(b.obs_names.astype(str)), **{g: pa.array(Xd[:, j]) for j, g in enumerate(genes)}})
    pcsv.write_csv(tbl, norm_csv)
    del Xd, tbl
R_CODE = f'''
suppressMessages({{library(Seurat); library(AUCell); library(GSEABase); library(data.table)}})
set.seed(123)
load("{GSC}")
norm_1k <- fread("{norm_csv}")
norm_1k <- as.data.frame(norm_1k)
rownames(norm_1k) <- norm_1k[[1]]
norm_1k <- norm_1k[, -1, drop = FALSE]
norm_1k_so <- CreateSeuratObject(counts = t(norm_1k))
counts <- GetAssayData(object = norm_1k_so, layer = "counts")   # = slot "counts" of SeuratObject < 5
gut_sub <- subsetGeneSets(gsc, rownames(counts))
cells_AUC <- AUCell_run(counts, gut_sub)
cells_AUC_mat <- t(getAUC(cells_AUC))
max_col_indices <- max.col(cells_AUC_mat, ties.method = "first")
cell_labels <- colnames(cells_AUC_mat)[max_col_indices]
result <- data.frame(cells = rownames(cells_AUC_mat), label = cell_labels, stringsAsFactors = FALSE)
write.csv(result, "{anno_csv}", row.names = FALSE)
'''
if not os.path.exists(anno_csv):
    print('running AUCell ...', flush=True)
    subprocess.run(['Rscript', '-e', R_CODE], check=True)
anno = pd.read_csv(anno_csv).set_index('cells')

mapper = pd.read_csv(MAPPER)                                   # label mapping as in 04
mapper['cell_type_aucell'] = mapper['cell type'].str.replace(r'[+\-\/()]', ' ', regex=True)
mapper['cell_type_aucell'] = mapper['cell_type_aucell'].str.replace(r'\s+', ' ', regex=True).str.strip()
mapper['cell_type_aucell'] = (mapper.groupby('cell_type_aucell').cumcount().astype(str).replace('0', '', regex=False)
                              .radd(mapper['cell_type_aucell'] + ' ').str.strip())
mp = mapper.set_index('cell_type_aucell')
lab = anno.loc[b.obs_names, 'label'].str.rstrip()
b.obs['aucell_label'] = lab.values
for col, src in [('cell_type', 'cell type'), ('cell_type_short', 'cell type short'),
                 ('cell_category', 'cell category'), ('cell_type_general', 'cell type general')]:
    b.obs[col] = mp.loc[lab.values, src].values

# ── 5. var QC from matched counts; write ──
v = b.copy(); v.X = v.layers['raw']
sc.pp.calculate_qc_metrics(v, inplace=True, log1p=True)
b.var = v.var.copy(); del v
for c in b.obs.columns:
    if b.obs[c].dtype == object:
        b.obs[c] = b.obs[c].astype(str)
b.write_h5ad(OUT_H5)
print('wrote', OUT_H5, b.shape)

# ── 6. validation ──
raw_sum = np.asarray(b.layers['raw'].sum(1)).ravel()
X375 = (b.obs['fov'].values <= 29)
val = []
def add(check, subset, n, value, expected):
    val.append(dict(check=check, subset=subset, n_cells=n, value=value, expected=expected))
add('Spearman rho, raw RNA row sum vs nCount_RNA (input object)', 'all cells', len(b), stats.spearmanr(old_sum, b.obs.nCount_RNA)[0], '~0 (misaligned)')
add('Spearman rho, raw RNA row sum vs nCount_RNA (aligned)', 'all cells', len(b), stats.spearmanr(raw_sum, b.obs.nCount_RNA)[0], '1')
add('Fraction raw row sum == nCount_RNA (aligned)', 'all cells', len(b), float(np.mean(raw_sum == b.obs.nCount_RNA.values)), '1')
add('Duplicate obs names (input object)', 'all cells', len(a), int(a.obs_names.duplicated().sum()), '-')
add('Duplicate obs names (aligned)', 'all cells', len(b), int(b.obs_names.duplicated().sum()), '0')
add('Input raw row = exprMat row at metadata position', 'all cells', len(b), rc['frac_old_rows_equal_exprMat_row_at_metadata_position'], 'root cause')
add('Cells passing RNA QC thresholds of 04 (aligned)', 'all cells', len(b), int(b.obs.rna_qc_pass.sum()), '-')
lg = lambda M: np.log1p(np.asarray(M.toarray()).ravel())
for gene, prot in [('EPCAM', 'EpCAM'), ('PTPRC', 'CD45')]:
    for tag, Xr in [('input object', a.X), ('aligned', b.X)]:
        r = stats.spearmanr(lg(Xr[:, genes.index(gene)]), b.obs[f'{prot}_geo'].values)[0]
        add(f'Spearman rho, {gene} RNA (log1p norm) vs {prot} protein (_geo), {tag}', 'all cells', len(b), r, '>0 if aligned')
    for tag, Xr in [('input object', a.X), ('aligned', b.X)]:
        r = stats.spearmanr(lg(Xr[X375][:, genes.index(gene)]), b.obs.loc[X375, f'{prot}_geo'].values)[0]
        add(f'Spearman rho, {gene} RNA (log1p norm) vs {prot} protein (_geo), {tag}', 'CD ileum 375i (FOV <= 29)', int(X375.sum()), r, '>0 if aligned')
ep = b.obs['cell_category'] == 'Epithelial'
for prot in ['EpCAM', 'CD45']:
    add(f'Median {prot} protein, RNA-Epithelial minus other cells (aligned)', 'all cells', len(b),
        float(b.obs.loc[ep, f'{prot}_geo'].median() - b.obs.loc[~ep, f'{prot}_geo'].median()), 'EpCAM > 0, CD45 < 0')
old_gen = a.obs['cell_type_general'].astype(str).values
add('Fraction of cells with unchanged cell_type_general', 'all cells', len(b), float(np.mean(old_gen == b.obs.cell_type_general.values)), '-')
V = pd.DataFrame(val)
V['value'] = V['value'].map(lambda x: float(f'{x:.3g}'))
V.to_csv(f'{OUT}/cd_multiomics_rna_realign_validation.csv', index=False)
pd.set_option('display.width', 250)
print(V.to_string())
print(b.obs.loc[X375 & b.obs.rna_qc_pass].cell_type_general.value_counts().to_string())
