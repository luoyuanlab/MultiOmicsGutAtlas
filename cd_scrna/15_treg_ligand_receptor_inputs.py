#!/usr/bin/env python
"""CD scRNA-seq: sender and Treg expression tables for the NicheNet step (16_treg_ligand_receptor_prioritization.R).

Inputs:  objects/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad, treg_pseudotime.h5ad (12_treg_pseudotime.py)
Outputs: objects/nichenet_inputs/strict_treg_imm_{meta,counts}.parquet (immune senders) and
         strict_treg_imm_general_{meta,counts}.parquet (immune, stromal and epithelial senders)
Panels:  none (inputs of Fig 4H-I, Supp 8G-H and the immune-sender ranking used by the Treg Suppression Score)
Run:     python <file>  or  RUN_SECTIONS=<id>[,<id>] python <file>  (sections: immune_sender_inputs, F4HI_nichenet_inputs)
"""
import os
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

warnings.filterwarnings('ignore')
_RUN_SECTIONS = [s.strip() for s in os.environ.get('RUN_SECTIONS', '').split(',') if s.strip()]
D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
OUT = f'{D}/objects/nichenet_inputs'


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


def load_labelled():
    """Atlas subset with pseudotime Tregs labelled 'Treg_B{pt_bin}' (other Tregs dropped) and sender labels in ct_lr_new."""
    adata = ad.read_h5ad(f'{D}/objects/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad')
    treg = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad', backed='r')
    tobs = treg.obs[['pt_bin', 'dpt_pseudotime']].copy()
    tregs_meta = pd.DataFrame({'index': tobs.index, 'ct_lr': ['Treg_B' + str(int(b)) for b in tobs['pt_bin']]})

    mask = (adata.obs['cell type general'] != 'Treg') | (adata.obs_names.isin(tregs_meta['index']))
    adata = adata[mask.values].copy()
    obs = adata.obs.copy()
    obs['index'] = obs.index.astype(str)
    obs = obs.merge(tregs_meta, on='index', how='left').set_index(obs.index)
    gen = obs['cell type general'].astype(str).replace({'Transitory cDC2/Mac': 'DC'})  # Transitory cDC2/Mac grouped with DC
    ct_lr_new = []
    for i, j, k in zip(obs['cell category'].astype(str), gen, obs['ct_lr'].tolist()):
        if j == 'Treg':
            ct_lr_new.append(k)
        elif i != 'Myeloid':
            ct_lr_new.append(i if i != 'T' else 'Other T')
        else:
            if 'DC' in j: ct_lr_new.append(j)
            elif j == 'Neutrophil': ct_lr_new.append('Granulocyte')
            elif 'Mac' in j: ct_lr_new.append('Macrophage')
            else: raise ValueError(j)
    obs['ct_lr_new'] = ct_lr_new
    adata.obs = obs
    return adata


def write_dense_parquet(a, path, chunk=20000):
    """Dense count parquet written in row chunks (avoids a 77 GB in-memory DataFrame)."""
    genes = a.var_names.tolist()
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    schema = pa.schema([pa.field(g, pa.float32()) for g in genes] + [pa.field('__index_level_0__', pa.string())])
    w = pq.ParquetWriter(path, schema, compression='snappy')
    for s in range(0, a.n_obs, chunk):
        blk = X[s:s + chunk].toarray().astype(np.float32)
        cols = [pa.array(blk[:, j]) for j in range(blk.shape[1])] + [pa.array(a.obs['index'].iloc[s:s + chunk].tolist())]
        w.write_table(pa.Table.from_arrays(cols, schema=schema))
    w.close()


def write_meta(a, path):
    m = a.obs.copy(); m['ct_lr'] = m['ct_lr_new'].tolist()
    for c in m.columns:
        if str(m[c].dtype) == 'category': m[c] = m[c].astype(str)
    m.to_parquet(path)


# ---- Immune senders (all labels except Epithelial, Endothelial, Connective); read by 16 section immune_sender_ranking
if _run('immune_sender_inputs'):
    adata = load_labelled()
    adata_imm = adata[(~adata.obs['ct_lr_new'].isin(['Epithelial', 'Endothelial', 'Connective'])).values].copy()
    print('immune subset', adata_imm.shape)
    write_meta(adata_imm, f'{OUT}/strict_treg_imm_meta.parquet')
    write_dense_parquet(adata_imm, f'{OUT}/strict_treg_imm_counts.parquet')


# ---- Immune, fibroblast, endothelial and epithelial senders (epithelial: random 30,000, seed 0); read by 16 section F4HI_nichenet_lr
if _run('F4HI_nichenet_inputs'):
    adata = load_labelled()
    adata.obs['ct_lr_new'] = adata.obs['ct_lr_new'].replace({'Connective': 'Fibroblast'})
    epi = np.where(adata.obs['ct_lr_new'].values == 'Epithelial')[0]
    drop = np.setdiff1d(epi, np.random.default_rng(0).choice(epi, min(30000, len(epi)), replace=False))
    keep = np.ones(adata.n_obs, bool); keep[drop] = False
    adata_imm = adata[keep].copy()
    print('immune + stromal/epithelial subset', adata_imm.shape)
    write_meta(adata_imm, f'{OUT}/strict_treg_imm_general_meta.parquet')
    write_dense_parquet(adata_imm, f'{OUT}/strict_treg_imm_general_counts.parquet')
