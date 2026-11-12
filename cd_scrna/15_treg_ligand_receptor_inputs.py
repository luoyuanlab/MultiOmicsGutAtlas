#!/usr/bin/env python
"""CD scRNA-seq: sender and Treg expression tables for the NicheNet ligand-receptor prioritization (16_treg_ligand_receptor_prioritization.R).

This step is in Python because the atlas and Treg objects are AnnData files; NicheNet itself runs in R (16).

Panels made
  none (inputs of Fig 4H-I and Supp 8G-H; inputs of the immune-sender ranking that defines the antigen-presentation LR module
  excluded from the Treg Suppression Score)

Inputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata_rowfix.h5ad,
    treg_pseudotime.h5ad and treg_pseudotime_bins.csv (12_treg_pseudotime.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/nichenet_inputs/ (strict_treg_imm_* and
    strict_treg_imm_general_* meta and count tables)

Run order of sections
  immune_sender_inputs -> F4HI_nichenet_inputs

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist. Command-line arguments and environment variables used by a section are listed in its own
  docstring.
  Paths are the lab's absolute paths (to be replaced by placeholders before public release).
  Environment: see the README section "Environment" (OPENBLAS_CORETYPE=SkylakeX on Intel Sapphire Rapids nodes).
"""
import os as _os
_THIS_DIR = _os.path.dirname(_os.path.abspath(__file__))   # repository folder of this file (used for utils/ imports)
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    """True if `section` should run: all sections run in file order unless RUN_SECTIONS=id1,id2,... is set."""
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# ======================================================================================================================
# Section immune_sender_inputs: NicheNet inputs with immune senders and pseudotime-binned Tregs.
# Panels: Treg Suppression Score prerequisite (Fig 5, Supp 9, Fig 6A): immune-sender NicheNet inputs for immune_sender_ranking (16_treg_ligand_receptor_prioritization.R)
# ======================================================================================================================
if _run('immune_sender_inputs'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/figure3/immune_sender_inputs.py'   # working folder of this step (relative paths resolve from here)
    #!/usr/bin/env python
    """NicheNet inputs with immune senders (read by section immune_sender_ranking of 16_treg_ligand_receptor_prioritization.R).

Cells: control/inactive/active object from 12_treg_pseudotime.py (703,284 cells). Tregs are restricted to the 2,941
pseudotime Tregs and labelled 'Treg_B{pt_bin}' (bins 0-4); other Tregs are dropped. Sender labels (ct_lr_new):
non-myeloid cells by 'cell category' (T cells other than Tregs = 'Other T'); myeloid cells by 'cell type general'
(DC types kept, with 'Transitory cDC2/Mac' grouped into 'DC'; Neutrophil -> Granulocyte; macrophage types -> Macrophage).
Dense count parquets are written in row chunks (avoids a 77 GB in-memory DataFrame).
Immune subset = all labels except Epithelial, Endothelial and Connective.
Outputs: objects/nichenet_inputs/strict_treg_imm_{meta,counts}_rowfix.parquet.
"""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, pyarrow as pa, pyarrow.parquet as pq, scipy.sparse as sp

    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
    OUT = f'{D}/objects/nichenet_inputs'
    log = open(f'{D}/logs/immune_sender_inputs.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()

    adata = ad.read_h5ad(f'{D}/objects/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata_rowfix.h5ad')
    treg = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad', backed='r')
    tobs = treg.obs[['pt_bin', 'dpt_pseudotime']].copy()
    bins = pd.read_csv(f'{D}/objects/treg_pseudotime_bins.csv', index_col=0)
    assert (bins.loc[tobs.index, 'pt_bin'].astype(int).values == tobs['pt_bin'].astype(int).values).all()
    P('control/inactive/active object', adata.shape, '| pseudotime Tregs', tobs.shape[0], '| in object', tobs.index.isin(adata.obs_names).sum())
    P('general label of strict Tregs in subset:', adata.obs.loc[tobs.index, 'cell type general'].value_counts().to_dict())
    assert (adata.obs.loc[tobs.index, 'cell type general'] == 'Treg').all()

    # ---- Tregs labelled 'Treg_B{pt_bin}' ----
    tregs_meta = pd.DataFrame({'index': tobs.index, 'ct_lr': ['Treg_B' + str(int(b)) for b in tobs['pt_bin']]})

    # ---- sender labels ----
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
    P('all_others ct_lr_new:\n' + adata.obs['ct_lr_new'].value_counts().to_string())

    # ---- immune subset ----
    imm_mask = ~adata.obs['ct_lr_new'].isin(['Epithelial', 'Endothelial', 'Connective'])
    adata_imm = adata[imm_mask.values].copy()
    P('imm subset', adata_imm.shape, '\n' + adata_imm.obs['ct_lr_new'].value_counts().to_string())

    def write_dense_parquet(a, path, chunk=20000):
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

    write_meta(adata_imm, f'{OUT}/strict_treg_imm_meta_rowfix.parquet')
    write_dense_parquet(adata_imm, f'{OUT}/strict_treg_imm_counts_rowfix.parquet')
    P('wrote imm parquet')


# ======================================================================================================================
# Section F4HI_nichenet_inputs: NicheNet inputs with immune, stromal and epithelial senders and pseudotime-binned Tregs.
# Panels: Fig 4H-I inputs
# ======================================================================================================================
if _run('F4HI_nichenet_inputs'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/figure3/F4HI_nichenet_inputs.py'   # working folder of this step (relative paths resolve from here)
    #!/usr/bin/env python
    """NicheNet inputs for Fig 4H-I (read by section F4HI_nichenet_lr of 16_treg_ligand_receptor_prioritization.R).
Senders = general cell types: the immune groups of immune_sender_inputs plus fibroblasts (Connective -> 'Fibroblast',
including myofibroblasts and pericytes), endothelial and epithelial cells (epithelial cells: random 30,000, seed 0, to
bound the dense parquet; ligand detection fractions are unaffected by random subsampling).

Cells: control/inactive/active object from 12_treg_pseudotime.py (703,284 cells). Tregs are restricted to the 2,941
pseudotime Tregs and labelled 'Treg_B{pt_bin}' (bins 0-4); other Tregs are dropped. Sender labels (ct_lr_new):
non-myeloid cells by 'cell category' (T cells other than Tregs = 'Other T'); myeloid cells by 'cell type general'
(DC types kept, with 'Transitory cDC2/Mac' grouped into 'DC'; Neutrophil -> Granulocyte; macrophage types -> Macrophage).
Dense count parquets are written in row chunks (avoids a 77 GB in-memory DataFrame).
Outputs: objects/nichenet_inputs/strict_treg_imm_general_{meta,counts}_rowfix.parquet.
"""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, pyarrow as pa, pyarrow.parquet as pq, scipy.sparse as sp

    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
    OUT = f'{D}/objects/nichenet_inputs'
    log = open(f'{D}/logs/F4HI_nichenet_inputs.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()

    adata = ad.read_h5ad(f'{D}/objects/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata_rowfix.h5ad')
    treg = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad', backed='r')
    tobs = treg.obs[['pt_bin', 'dpt_pseudotime']].copy()
    bins = pd.read_csv(f'{D}/objects/treg_pseudotime_bins.csv', index_col=0)
    assert (bins.loc[tobs.index, 'pt_bin'].astype(int).values == tobs['pt_bin'].astype(int).values).all()
    P('control/inactive/active object', adata.shape, '| pseudotime Tregs', tobs.shape[0], '| in object', tobs.index.isin(adata.obs_names).sum())
    P('general label of strict Tregs in subset:', adata.obs.loc[tobs.index, 'cell type general'].value_counts().to_dict())
    assert (adata.obs.loc[tobs.index, 'cell type general'] == 'Treg').all()

    # ---- Tregs labelled 'Treg_B{pt_bin}' ----
    tregs_meta = pd.DataFrame({'index': tobs.index, 'ct_lr': ['Treg_B' + str(int(b)) for b in tobs['pt_bin']]})

    # ---- sender labels ----
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
    P('all_others ct_lr_new:\n' + adata.obs['ct_lr_new'].value_counts().to_string())

    g2 = adata.obs['cell type general'].astype(str)
    P('Connective composition:', g2[adata.obs['ct_lr_new'] == 'Connective'].value_counts().to_dict())
    adata.obs['ct_lr_new'] = adata.obs['ct_lr_new'].replace({'Connective': 'Fibroblast'})
    epi = np.where(adata.obs['ct_lr_new'].values == 'Epithelial')[0]
    drop = np.setdiff1d(epi, np.random.default_rng(0).choice(epi, min(30000, len(epi)), replace=False))
    keep = np.ones(adata.n_obs, bool); keep[drop] = False
    adata_imm = adata[keep].copy()
    P('immune + general stromal/epithelial subset', adata_imm.shape, '| epithelial before / after', len(epi), (adata_imm.obs['ct_lr_new'] == 'Epithelial').sum(),
      '\n' + adata_imm.obs['ct_lr_new'].value_counts().to_string())

    def write_dense_parquet(a, path, chunk=20000):
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

    write_meta(adata_imm, f'{OUT}/strict_treg_imm_general_meta_rowfix.parquet')
    write_dense_parquet(adata_imm, f'{OUT}/strict_treg_imm_general_counts_rowfix.parquet')
    P('wrote imm + general parquet')
