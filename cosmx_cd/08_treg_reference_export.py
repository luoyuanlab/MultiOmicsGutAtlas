#!/usr/bin/env python
"""CD scRNA-seq Tregs: label-transfer reference (raw counts and pseudotime bins) as CSV for 09_treg_label_transfer.R.

Input   objects/treg_pseudotime.h5ad (cd_scrna/12_treg_pseudotime.py)
Outputs objects/label_transfer/treg_raw_counts.csv, treg_pt_bin.csv (read by 09)
"""
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp

D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
OUT = f'{D}/objects/label_transfer'

t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
X = t.layers['counts']
X = X.toarray() if sp.issparse(X) else np.asarray(X)
X = np.round(X).astype(np.int32)
pd.DataFrame(X, index=t.obs_names, columns=t.var_names).to_csv(f'{OUT}/treg_raw_counts.csv')
pd.DataFrame({'pt_bin': t.obs['pt_bin'].astype(int).values}, index=t.obs_names).to_csv(f'{OUT}/treg_pt_bin.csv')
