#!/usr/bin/env python
"""CD scRNA-seq atlas: Treg selection, diffusion pseudotime and the Treg pseudotime panels.

Panels made
  Fig 4D  Treg subtype module scores on the Treg UMAP (2,941 Tregs)
  Fig 4E  diffusion pseudotime (left); Suppressive Core Score vs pseudotime, LMM and within-patient sign test (right)
  Fig 4F  marker genes across pseudotime bins

Inputs
  /share/fsmresfiles/UC/scRNA-seq/merged_cd/cleaned_annoed_all_cell_types_rowfix.h5ad (CD atlas, 09_atlas_h5ad_conversion.ipynb)
  /share/fsmresfiles/UC/scRNA-seq/merged_cd/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad (control, inactive and
    active CD cells without cycling cells; its obs is used as is)
  /share/fsmresfiles/UC/scRNA-seq/merged_cd/treg_ptbin_markers.xlsx (Fig 4F marker list)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/ (Treg input objects; treg_pseudotime.h5ad and
    treg_pseudotime_bins.csv, used by 13-15, cosmx_cd and treatment_validation)
  /share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/ (panels and CSVs)

Run order of sections
  treg_inputs -> treg_selection_dpt -> doublet_branch_dpt -> F4DF_module_scores_markers -> F4E_pseudotime_umap
  -> F4E_suppression_scatter

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
# Section treg_inputs: Treg input objects (control/inactive/active atlas without cycling cells and its Treg subset).
# Panels: prerequisite
# ======================================================================================================================
if _run('treg_inputs'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/figure3/treg_inputs.py'   # working folder of this step (relative paths resolve from here)
    """Treg input objects: control/inactive/active atlas without cycling cells and its Treg subset.

The cells, their order and the annotation columns come from ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad
(control, inactive and active CD samples without the P17 pooled-patient dataset and without cycling cells). Expression
is taken from the CD atlas (cleaned_annoed_all_cell_types_rowfix.h5ad) by cell name.
Arguments: --skip-full  write only the Treg subset.
Outputs (objects/):
  ctrl_inact_act_sccoda_no_cycle_correct_cat_adata_rowfix.h5ad   full object (703,284 cells)
  tregs_input_rowfix.h5ad                                         'cell type' == 'Treg' subset
"""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, time, sys
    D = '/share/fsmresfiles/UC/scRNA-seq/merged_cd/'
    OUT = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/'
    t0 = time.time()

    atlas = ad.read_h5ad(D + 'cleaned_annoed_all_cell_types_rowfix.h5ad')
    print('atlas', atlas.shape, f'{time.time()-t0:.0f}s', flush=True)
    der = ad.read_h5ad(D + 'ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad')
    print('control/inactive/active object', der.shape, f'{time.time()-t0:.0f}s', flush=True)
    assert np.array_equal(der.var_names.values, atlas.var_names.values)
    names = der.obs_names.values.astype(str)
    pos = atlas.obs_names.get_indexer(names)
    assert (pos >= 0).all()

    # same obs, X rows taken from the atlas by cell name
    X = atlas.X[pos]
    full = ad.AnnData(X=X, obs=der.obs.copy(), var=der.var.copy())
    del atlas

    # write outputs
    treg = np.where(der.obs['cell type'].astype(str).values == 'Treg')[0]
    tr = full[treg].copy()
    tr.write_h5ad(OUT + 'objects/tregs_input_rowfix.h5ad')
    print('wrote Treg inputs', tr.shape, f'{time.time()-t0:.0f}s', flush=True)
    if '--skip-full' not in sys.argv:
        full.write_h5ad(OUT + 'objects/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata_rowfix.h5ad')
        print('wrote full control/inactive/active object', full.shape, f'{time.time()-t0:.0f}s', flush=True)
    print('done')


# ======================================================================================================================
# Section treg_selection_dpt: Treg selection (doublet removal, two filtering rounds, strict Treg cluster) and DPT.
# Panels: prerequisite (Treg selection + DPT)
# ======================================================================================================================
if _run('treg_selection_dpt'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/figure3/treg_selection_dpt.py'   # working folder of this step (relative paths resolve from here)
    #!/usr/bin/env python
    """Treg selection and diffusion pseudotime (scanpy 1.10.2, anndata 0.11, harmonypy; run with OPENBLAS_CORETYPE=SkylakeX).

Steps
  1. Treg subset (objects/tregs_input_rowfix.h5ad); scrublet on raw counts, cells with doublet score > 0.2 removed.
  2. normalise, log1p, 3,000 HVGs (seurat_v3), scale, PCA (50), Harmony on 'dataset', kNN (k=20, 30 PCs), UMAP,
     Leiden; module scores (sc.tl.score_genes on the scaled X; adata.raw is cleared after the scrublet input is built).
  3. round 1: 0 < Treg_core < 2.5 and Cell_cycle < 2.5; re-embed.
  4. round 2: drop leiden_0.3 clusters whose median Treg_core is <= the 30th percentile of cells or whose median
     Cell_cycle is >= the 70th percentile; re-embed and re-score.
  5. strict Treg cluster: the leiden_0.1 cluster with the highest median Treg_core score whose median Cell_cycle score
     is below the 70th percentile of cells (high-confidence Treg cluster by core markers and low cycling).
  6. PAGA, PAGA-initialised UMAP; root = first cell of the leiden_0.6 sub-cluster (>= 20 cells) with the highest median
     Treg_core; DPT on the kNN graph from round 2 (X_pca_harmony, k=20, 30 PCs); a k=8 graph is built afterwards for
     display; five quantile bins of DPT.
Output: objects/treg_pseudotime_raw.h5ad (cluster statistics, doublet scores and obs tables in objects/intermediate/).
"""
    import warnings; warnings.filterwarnings('ignore')
    import sys, json, time
    import numpy as np, pandas as pd, scanpy as sc, scanpy.external as sce, anndata as ad
    BASE = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/'
    OUT = BASE
    log = open(f'{OUT}logs/treg_selection_dpt.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()
    t0 = time.time()

    # ---- Treg subset, raw copy, scrublet --------------------------------------------------------------------------
    tregs = ad.read_h5ad(f'{BASE}objects/tregs_input_rowfix.h5ad')
    adata = tregs.copy()
    adata.raw = adata.copy()
    adata_raw = adata.raw.to_adata()
    adata.raw = None   # module scores are computed on the scaled X (see docstring)
    adata_raw.X = pd.DataFrame(adata_raw.X.toarray()).astype(int).values
    sc.pp.scrublet(adata_raw, expected_doublet_rate=0.1, n_prin_comps=50)
    P('scrublet predicted_doublet:', adata_raw.obs['predicted_doublet'].value_counts().to_dict())
    adata_raw.obs['predicted_doublet_0.2'] = [0 if i <= 0.2 else 1 for i in adata_raw.obs['doublet_score'].tolist()]
    doublet = adata_raw[adata_raw.obs['predicted_doublet_0.2'] == 1].obs_names
    adata = adata[~adata.obs_names.isin(doublet)]
    adata_raw_nodoub = adata_raw[adata_raw.obs_names.isin(adata.obs_names.tolist())]
    adata.layers['counts'] = adata_raw_nodoub.X.copy()
    adata_raw.obs[['doublet_score', 'predicted_doublet_0.2']].to_csv(f'{OUT}objects/intermediate/doublet_scores.csv')
    P('cells after doublet removal:', adata.n_obs)

    # ---- normalise, HVG, scale, PCA, Harmony, neighbours, UMAP, Leiden, module scores ----------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3', layer='counts')
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    sce.pp.harmony_integrate(adata, key='dataset')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=20, n_pcs=30)
    sc.tl.umap(adata)
    for r in [0.1, 0.2, 0.3]:
        sc.tl.leiden(adata, resolution=r, key_added=f'leiden_{r}')
    panels = {
     'Activated_tissue': ['BATF', 'PRDM1', 'MAF', 'TNFRSF18', 'TNFRSF4'],
     'Tfr': ['BCL6', 'CXCR5', 'PDCD1', 'ICOS'],
     'Th17_like': ['RORC', 'CCR6', 'IL17A', 'IL17F'],
     'IFN_response': ['ISG15', 'IFI6', 'MX1', 'IFIT3'],
     'Cell_cycle': ['MKI67', 'TOP2A', 'STMN1'],
     'Treg_core': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT', 'ENTPD1', 'LRRC32', 'IL10']
    }
    for name, genes in panels.items():
        sc.tl.score_genes(adata, [g for g in genes if g in adata.var_names], score_name=name)
    adata.obs.to_csv(f'{OUT}objects/intermediate/stage1_tregs_nodoublet_obs.csv')

    # ---- round 1 = 0 < Treg_core < 2.5 & Cell_cycle < 2.5, re-embed ----------------------------------------------
    adata_clean = adata.copy()
    treg_strict = adata_clean[adata_clean.obs['Treg_core'] > 0].copy()
    treg_strict = treg_strict[treg_strict.obs['Cell_cycle'] < 2.5].copy()
    treg_strict = treg_strict[treg_strict.obs['Treg_core'] < 2.5].copy()
    adata_clean = treg_strict
    P('cells after round 1:', adata_clean.n_obs)
    sc.pp.highly_variable_genes(adata_clean, n_top_genes=3000, flavor='seurat_v3', layer='counts')
    sc.pp.scale(adata_clean, max_value=10)
    sc.tl.pca(adata_clean, n_comps=50, use_highly_variable=True)
    sce.pp.harmony_integrate(adata_clean, key='dataset')
    sc.pp.neighbors(adata_clean, use_rep='X_pca_harmony', n_neighbors=20, n_pcs=30)
    sc.tl.umap(adata_clean)
    for r in [0.1, 0.2, 0.3, 0.4]:
        sc.tl.leiden(adata_clean, resolution=r, key_added=f'leiden_{r}')
    adata_clean.obs.to_csv(f'{OUT}objects/intermediate/stage2_round1_obs.csv')

    # ---- round 2 = drop leiden_0.3 clusters with low Treg_core / high Cell_cycle -----------------------------------
    adata_clean2 = adata_clean.copy()
    cluster_key = 'leiden_0.3'
    cluster_stats = adata_clean2.obs.groupby(cluster_key, observed=False)[['Treg_core', 'Cell_cycle']].median()
    bad_clusters = cluster_stats.index[
        (cluster_stats['Treg_core'] <= adata_clean2.obs['Treg_core'].quantile(0.3)) |
        (cluster_stats['Cell_cycle'] >= adata_clean2.obs['Cell_cycle'].quantile(0.7))]
    P('round2 removing clusters:', list(bad_clusters), '\n', cluster_stats.round(3).to_string())
    adata_clean2 = adata_clean2[~adata_clean2.obs[cluster_key].isin(bad_clusters)].copy()
    P('cells after round 2:', adata_clean2.n_obs)
    sc.pp.highly_variable_genes(adata_clean2, n_top_genes=3000, flavor='seurat_v3', layer='counts')
    sc.pp.scale(adata_clean2, max_value=10)
    sc.tl.pca(adata_clean2, n_comps=50, use_highly_variable=True)
    sce.pp.harmony_integrate(adata_clean2, key='dataset')
    sc.pp.neighbors(adata_clean2, use_rep='X_pca_harmony', n_neighbors=20, n_pcs=30)
    sc.tl.umap(adata_clean2)
    for r in [0.1, 0.2, 0.3, 0.4]:
        sc.tl.leiden(adata_clean2, resolution=r, key_added=f'leiden_{r}')
    for r in [0.5, 0.6]:
        sc.tl.leiden(adata_clean2, resolution=r, key_added=f'leiden_{r}')
    panels2 = {
     'Suppresive': ['BATF', 'PRDM1', 'MAF', 'TNFRSF18', 'TNFRSF4'],
     'Tfr': ['BCL6', 'CXCR5', 'PDCD1', 'ICOS'],
     'Th17_like': ['RORC', 'CCR6', 'IL17A', 'IL17F'],
     'IFN_response': ['ISG15', 'IFI6', 'MX1', 'IFIT3'],
     'Cell_cycle': ['MKI67', 'TOP2A', 'STMN1'],
     'Treg_core': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT', 'ENTPD1', 'LRRC32', 'IL10']
    }
    for name, genes in panels2.items():
        sc.tl.score_genes(adata_clean2, [g for g in genes if g in adata_clean2.var_names], score_name=name)
    adata_tmp = adata_clean2.copy()
    adata_tmp.X = adata_clean2.layers['counts'].copy()
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    adata_clean2.layers['log1p'] = adata_tmp.X.copy()
    del adata_tmp
    adata_clean2.obs.to_csv(f'{OUT}objects/intermediate/stage3_round2_obs.csv')

    # ---- strict Treg cluster: highest median Treg_core among clusters with low median Cell_cycle -------------------
    o2 = adata_clean2.obs
    # cluster medians, as in the round-2 cluster rule
    cs = o2.groupby('leiden_0.1', observed=True)[['Treg_core', 'Cell_cycle']].median().add_suffix('_median')
    cs = cs.join(o2.groupby('leiden_0.1', observed=True)[['Treg_core', 'Cell_cycle']].mean().add_suffix('_mean'))
    cs['n'] = o2['leiden_0.1'].value_counts()
    cyc_cut = o2['Cell_cycle'].quantile(0.7)
    elig = cs[cs['Cell_cycle_median'] < cyc_cut]
    if len(elig) == 0:
        elig = cs
    strict_cl = elig['Treg_core_median'].idxmax()
    P('leiden_0.1 cluster stats (round2):\n', cs.round(3).to_string(), f'\ncycling cut (cell q0.7)={cyc_cut:.3f}; eligible={list(elig.index)}; strict cluster chosen = {strict_cl}')
    cs.to_csv(f'{OUT}objects/intermediate/strict_cluster_choice_leiden01.csv')
    clus0 = adata_clean2[adata_clean2.obs['leiden_0.1'] == strict_cl]
    P('cells in strict Treg cluster:', clus0.n_obs)

    # ---- PAGA, PAGA-initialised UMAP, root, DPT, kNN(8), 5 quantile bins -------------------------------------------
    sc.tl.paga(clus0, groups='leiden_0.6')
    sc.pl.paga(clus0, plot=False)
    sc.tl.umap(clus0, init_pos='paga')
    clus0.obs['treg_stability_score'] = clus0.obs['Activated_tissue'] - clus0.obs['Th17_like']
    clus0.obs['treg_stability_score2'] = clus0.obs['Treg_core'] - clus0.obs['Th17_like']
    cs6 = clus0.obs.groupby('leiden_0.6', observed=True)[['Treg_core', 'Cell_cycle']].median().add_suffix('_median')
    cs6 = cs6.join(clus0.obs.groupby('leiden_0.6', observed=True)[['Treg_core', 'Cell_cycle']].mean().add_suffix('_mean'))
    cs6['n'] = clus0.obs['leiden_0.6'].value_counts()
    root_cl = cs6[cs6['n'] >= 20]['Treg_core_median'].idxmax()   # ignore tiny (<20-cell) sub-clusters
    P('leiden_0.6 sub-cluster stats (strict cluster):\n', cs6.round(3).to_string(), f'\nroot sub-cluster = {root_cl}')
    cs6.to_csv(f'{OUT}objects/intermediate/root_subcluster_choice_leiden06.csv')
    clus0.uns['iroot'] = np.flatnonzero(clus0.obs['leiden_0.6'] == root_cl)[0]
    clus0.uns['root_subcluster'] = str(root_cl); clus0.uns['strict_cluster_leiden01'] = str(strict_cl)
    sc.tl.dpt(clus0)   # diffmap/DPT on the kNN graph from round 2 (Harmony, k=20, 30 PCs)
    sc.pp.neighbors(clus0, n_neighbors=8)   # display graph, built after the DPT
    pt = clus0.obs['dpt_pseudotime']
    clus0.obs['pt_bin'] = pd.qcut(pt, q=5, labels=False)
    P('iroot', clus0.uns['iroot'], clus0.obs_names[clus0.uns['iroot']])
    P('pt_bin x leiden_0.6\n', pd.crosstab(clus0.obs['pt_bin'], clus0.obs['leiden_0.6']).to_string())
    clus0.obs['pt_bin'] = clus0.obs['pt_bin'].astype(str)
    clus0.obs['Suppressive'] = clus0.obs['Activated_tissue'].tolist().copy()
    clus0 = clus0.copy()
    clus0.write_h5ad(f'{OUT}objects/treg_pseudotime_raw.h5ad')
    P(f'done {time.time()-t0:.0f}s', clus0.shape)


# ======================================================================================================================
# Section doublet_branch_dpt: Remove a small separated high-pseudotime branch (if any) and compute the Treg DPT.
# Panels: prerequisite (writes treg_pseudotime.h5ad)
# ======================================================================================================================
if _run('doublet_branch_dpt'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/figure3/doublet_branch_dpt.py'   # working folder of this step (relative paths resolve from here)
    """Remove a small, separated high-pseudotime branch (Treg-epithelial doublets or ambient RNA, e.g. PHGR1+/EPCAM+) if
present, then compute the Treg DPT with the parameters of treg_selection_dpt: Harmony kNN k=20, 30 PCs; root = first
cell of the leiden_0.6 sub-cluster (>=20 cells) with the highest median Treg_core; five quantile bins.
The branch is defined by the largest gap in raw DPT (no cells between the main body and the branch; gap > 0.1 and
< 5% of cells), not by a set percentile; with the current Treg object no cells meet this rule.
Inputs: objects/treg_pseudotime_raw.h5ad. Outputs: objects/treg_pseudotime.h5ad, objects/treg_pseudotime_bins.csv."""
    import os, anndata as ad, scanpy as sc, numpy as np, pandas as pd, scipy.sparse as sp
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
    IN = f'{D}/objects/treg_pseudotime_raw.h5ad'; OUTO = f'{D}/objects/treg_pseudotime.h5ad'
    t = ad.read_h5ad(IN)
    pt = np.sort(t.obs['dpt_pseudotime'].values)
    gaps = np.diff(pt); g = int(np.argmax(gaps)); cut = (pt[g] + pt[g + 1]) / 2
    print(f'largest gap in DPT: {pt[g]:.3f} -> {pt[g+1]:.3f}; cut at {cut:.3f}')
    branch = t.obs['dpt_pseudotime'] > cut
    if not (gaps[g] > 0.1 and branch.mean() < 0.05):
        branch[:] = False; print('no small separated branch -> nothing removed')
    X = sp.csc_matrix(t.layers['counts']); vn = np.asarray(t.var_names)
    def det(m, gene): return round(float((X[np.where(m)[0]][:, np.where(vn == gene)[0][0]].toarray() > 0).mean() * 100), 1)
    rep = {gene: (det(branch.values, gene), det(~branch.values, gene)) for gene in ['PHGR1', 'EPCAM', 'KRT8', 'KRT18', 'FOXP3', 'PTPRC']}
    print(f'removing {int(branch.sum())} branch cells; detection % (branch, rest):', rep)
    pd.DataFrame(rep, index=['branch_pct', 'rest_pct']).T.to_csv(f'{D}/objects/intermediate/removed_branch_markers.csv')
    t.obs.loc[branch].to_csv(f'{D}/objects/intermediate/removed_branch_cells.csv')
    c = t[~branch.values].copy()
    for k in ['X_diffmap']:
        c.obsm.pop(k, None)
    c.uns.pop('diffmap_evals', None); c.uns.pop('iroot', None)
    sc.pp.neighbors(c, use_rep='X_pca_harmony', n_neighbors=20, n_pcs=30)
    cs6 = c.obs.groupby('leiden_0.6', observed=True)['Treg_core'].median().to_frame('Treg_core_median')
    cs6['n'] = c.obs['leiden_0.6'].value_counts()
    root_cl = cs6[cs6['n'] >= 20]['Treg_core_median'].idxmax()
    c.uns['iroot'] = int(np.flatnonzero(c.obs['leiden_0.6'] == root_cl)[0]); c.uns['root_subcluster'] = str(root_cl)
    sc.tl.diffmap(c); sc.tl.dpt(c)
    c.obs['pt_bin'] = pd.qcut(c.obs['dpt_pseudotime'], q=5, labels=False)
    q = c.obs['dpt_pseudotime'].quantile([.05, .25, .5, .75, .95, 1]).round(3).to_dict()
    print('root sub-cluster', root_cl, '| new DPT quantiles', q, '| cells', c.n_obs, '| patients', c.obs['patient'].nunique())
    c.write_h5ad(OUTO)
    o = c.obs[['pt_bin', 'dpt_pseudotime', 'leiden_0.6', 'patient', 'sample', 'condition_summary', 'dataset']].copy(); o['dpt_pseudotime_rank'] = o['dpt_pseudotime'].rank(pct=True).round(4)
    o.to_csv(f'{D}/objects/treg_pseudotime_bins.csv')


# ======================================================================================================================
# Section F4DF_module_scores_markers: Treg subtype module scores on the Treg UMAP and marker genes across pseudotime bins.
# Panels: Fig 4D, 4F
# ======================================================================================================================
if _run('F4DF_module_scores_markers'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/F4DF_module_scores_markers.py'   # working folder of this step (relative paths resolve from here)
    """Fig 4D: Suppressive, Tfr and Th17-like module scores on the Treg UMAP (Reds, edges of a k=8 graph, shared vmax,
na_color #f2f2f2). Fig 4F: Wilcoxon test per pseudotime bin (non-Control Tregs, log1p); curated gene list from
treg_ptbin_markers.xlsx minus KLRD1 (not expressed) and ARID5B; colour scale -2..2.
Input: objects/treg_pseudotime.h5ad. Outputs: fig3D_module_scores.png, fig3F_marker_dotplot.png/.csv.
Run with OPENBLAS_CORETYPE=SkylakeX."""
    import os, anndata as ad, scanpy as sc, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns, matplotlib.colors as mcolors
    from matplotlib.lines import Line2D
    from scipy.stats import mannwhitneyu
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals'
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    panels = {'Suppressive': ['BATF', 'PRDM1', 'MAF', 'TNFRSF18', 'TNFRSF4'], 'Tfr': ['BCL6', 'CXCR5', 'PDCD1', 'ICOS'], 'Th17_like': ['RORC', 'CCR6', 'IL17A', 'IL17F']}
    for k, g in panels.items():
        sc.tl.score_genes(t, [x for x in g if x in t.var_names], score_name=k)
    p = t.copy(); sc.pp.neighbors(p, n_neighbors=8, use_rep='X_pca_harmony')   # graph used only to draw edges
    # ---- Fig 4D
    vmax = np.percentile(np.concatenate([p.obs[k].dropna() for k in panels]), 99)
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 4.0), gridspec_kw={'wspace': 0.0})
    for a, k in zip(ax, panels):
        sc.pl.umap(p, color=k, edges=True, edges_width=0.05, edges_color='#d9d9d9', cmap='Reds', vmin=0, vmax=vmax, frameon=False, na_color='#f2f2f2', ax=a, show=False, colorbar_loc=None, s=30)
        a.set_title(k, fontsize=16)
    fig.subplots_adjust(wspace=0.0, left=0.0, right=0.93)
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=mcolors.Normalize(0, vmax)); fig.colorbar(sm, cax=fig.add_axes([0.94, 0.15, 0.012, 0.65]))
    fig.savefig(f'{OUT}/fig3D_module_scores.png', dpi=300, bbox_inches='tight')
    # ---- Fig 4F
    f = t[t.obs['condition_summary'].astype(str) != 'Control'].copy(); f.X = f.layers['log1p']
    f.obs['pt_bin'] = f.obs['pt_bin'].astype(int).astype(str).astype('category')
    sc.tl.rank_genes_groups(f, groupby='pt_bin', method='wilcoxon', pts=True)
    res = sc.get.rank_genes_groups_df(f, group=None); res['group'] = res['group'].astype(str)
    def plot3f(genes, fname, order_by_peak=True, cap=1.5):
        sub = res[res.names.isin(genes)]
        lfc = sub.pivot(index='group', columns='names', values='logfoldchanges').reindex(['0', '1', '2', '3', '4'])
        padj = sub.pivot(index='group', columns='names', values='pvals_adj').reindex(['0', '1', '2', '3', '4'])
        if order_by_peak:
            pk = lfc.fillna(0).idxmax(axis=0); order = sorted(lfc.columns, key=lambda g: (int(pk[g]), -lfc[g].max()))
        else:
            order = [g for g in genes if g in lfc.columns]
        lfc, padj = lfc[order], padj[order]
        cmap = mcolors.LinearSegmentedColormap.from_list('logfc', ['#0000FF', '#FFFFFF', '#FF0000'])
        norm = mcolors.Normalize(-cap, cap)   # colour saturates at +/-cap so moderate effects are visible
        size = lambda z: 30 + min(abs(z), 3) * 120
        fig, ax = plt.subplots(figsize=(0.36 * len(order) + 1.2, 3.3))
        for i, b in enumerate(lfc.index):
            for j, g in enumerate(order):
                x, q = lfc.loc[b, g], padj.loc[b, g]
                if pd.notna(x):
                    ax.scatter(j, i, s=size(x), c=[cmap(norm(np.clip(x, -cap, cap)))], alpha=1.0 if q < 0.1 else 0.35, edgecolors='black', linewidth=0.6)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=13); ax.set_yticks(range(5)); ax.set_yticklabels(lfc.index, fontsize=13)
        ax.set_xlabel('Genes', fontsize=15); ax.set_ylabel('Pseudotime Bin', fontsize=15); ax.grid(True, alpha=0.3, linestyle='--'); ax.set_axisbelow(True)
        ax.set_xlim(-0.6, len(order) - 0.4); ax.set_ylim(-0.6, 4.6)
        h1 = [Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='black', markersize=11, alpha=a, label=l) for a, l in [(1.0, 'p_adj < 0.1'), (0.35, 'p_adj ≥ 0.1')]]
        h2 = [Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='black', markersize=np.sqrt(size(z)), label=f'|Log2FC| = {z}') for z in (2, 1, 0.5)]
        l1 = ax.legend(handles=h1, loc='upper left', bbox_to_anchor=(1.01, 1.02), frameon=False, fontsize=11)
        ax.add_artist(l1)
        ax.legend(handles=h2, loc='upper left', bbox_to_anchor=(1.01, 0.70), frameon=False, fontsize=11, labelspacing=1.6, borderpad=0.8, handletextpad=1.0)
        cax = ax.inset_axes([1.36, 0.08, 0.022, 0.62]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
        cb.set_label('Log2FC', fontsize=12, labelpad=6); cb.ax.tick_params(labelsize=10)
        fig.savefig(f'{OUT}/{fname}.png', dpi=300, bbox_inches='tight')
        out = lfc.T.round(3).add_prefix('log2FC_bin').join(padj.T.apply(lambda c: c.map(lambda z: float('%.3g' % z))).add_prefix('padj_bin'))
        out.to_csv(f'{OUT}/{fname}.csv'); return out

    genes = pd.read_excel('/share/fsmresfiles/UC/scRNA-seq/merged_cd/treg_ptbin_markers.xlsx').query("plot == 'plot'")['names'].unique()
    plot3f([g for g in genes if g not in ('KLRD1', 'ARID5B')], 'fig3F_marker_dotplot', cap=2.0)


# ======================================================================================================================
# Section F4E_pseudotime_umap: Treg UMAP coloured by diffusion pseudotime.
# Panels: Fig 4E left
# ======================================================================================================================
if _run('F4E_pseudotime_umap'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/F4E_pseudotime_umap.py'   # working folder of this step (relative paths resolve from here)
    """Fig 4E (left): Treg UMAP coloured by diffusion pseudotime, rank-normalized to [0, 1] for display (bins are quintiles,
so this is the scale the bins use). Input: objects/treg_pseudotime.h5ad. Output: fig3E_pseudotime_umap.png.
Run with OPENBLAS_CORETYPE=SkylakeX."""
    import anndata as ad, scanpy as sc, pandas as pd, numpy as np, matplotlib.pyplot as plt, matplotlib.colors as mcolors
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals'
    purple = mcolors.LinearSegmentedColormap.from_list('p', ['#e8d4cf', '#a16f8d', '#2d223a'])
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    # ---- left: pseudotime UMAP
    t.obs['Pseudotime'] = t.obs['dpt_pseudotime'].rank(pct=True).values
    fig, ax = plt.subplots(figsize=(6.1, 4.6))
    sc.pl.umap(t, color='Pseudotime', cmap=purple, vmin=0, vmax=1, frameon=False, s=22, ax=ax, show=False)
    ax.set_title('Diffusion Pseudotime', fontsize=15)
    fig.savefig(f'{OUT}/fig3E_pseudotime_umap.png', dpi=300, bbox_inches='tight'); plt.close(fig)


# ======================================================================================================================
# Section F4E_suppression_scatter: Suppressive score vs Treg pseudotime with linear and non-linear trend tests.
# Panels: Fig 4E right
# ======================================================================================================================
if _run('F4E_suppression_scatter'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/F4E_suppression_scatter.py'   # working folder of this step (relative paths resolve from here)
    """Fig 4E (right): suppressive score vs Treg pseudotime, with a test of whether the decline is linear.

Inputs and definitions:
  - objects/treg_pseudotime.h5ad (2,941 CD-atlas Tregs, 65 patients)
  - suppressive score = sc.tl.score_genes(FOXP3, CTLA4, IL2RA, BATF, CD27) on layers['log1p'] (computed here; the stored
    obs['Suppressive'] column holds the Activated_tissue module and is not used)
  - pseudotime = diffusion pseudotime, rank-normalized to [0, 1] (the scale of the UMAP and of the quintile bins 0-4;
    raw DPT is heavily skewed, 40% of its axis holds <1% of cells), bin midpoints 0.1, 0.3, ..., 0.9
  - patient-bin points = per-patient mean score in a bin with >=3 Tregs
Figure (single panel, PNG): per-cell hexbin with a LOWESS curve (patient-cluster bootstrap 95% band) and the cell-level
LMM linear fit; statistics printed above the axes so they do not cover data. The patient-bin analysis is in the CSV and
in the legend text.
Statistics (all in Fig3E_suppression_score_vs_pseudotime.csv):
  - LMM score ~ pt + (1 + pt | patient), at patient-bin and at cell level (REML estimates)
  - Spearman rho, per cell and per patient-bin (pooled)
  - non-linearity: ML likelihood-ratio tests of linear vs quadratic, vs cubic B-spline (cell level, df=4),
    vs bin-as-category (patient-bin level, lack-of-fit); LOWESS values at bin midpoints; slope with bin 0 removed
  - per-patient OLS slopes (patients with >=3 bins) with two-sided sign test
Outputs: Fig3E_suppression_score_vs_pseudotime.csv/.png and _legend.txt.
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=16 RUN_SECTIONS=F4E_suppression_scatter python <this file>
"""
    import os
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX', 'set OPENBLAS_CORETYPE=SkylakeX (wrong matmul on Xeon 8468)'
    import warnings
    warnings.filterwarnings('ignore')
    import anndata as ad, scanpy as sc, pandas as pd, numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt, matplotlib.colors as mcolors
    import statsmodels.formula.api as smf
    from statsmodels.nonparametric.smoothers_lowess import lowess
    from scipy.stats import spearmanr, binomtest, chi2, linregress

    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals'
    STEM = f'{OUT}/Fig3E_suppression_score_vs_pseudotime'
    GENES = ['FOXP3', 'CTLA4', 'IL2RA', 'BATF', 'CD27']
    MID = {b: 0.1 + 0.2 * b for b in range(5)}
    rng = np.random.default_rng(0)
    purple = mcolors.LinearSegmentedColormap.from_list('p', ['#e8d4cf', '#a16f8d', '#2d223a'])   # Fig 4E palette

    # ---------------- data ----------------
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    t.X = t.layers['log1p']
    sc.tl.score_genes(t, GENES, score_name='Suppressive')
    o = t.obs[['patient', 'pt_bin', 'dpt_pseudotime', 'Suppressive']].copy()
    o['patient'] = o['patient'].astype(str); o['pt_bin'] = o['pt_bin'].astype(int)
    o['pt'] = o['dpt_pseudotime'].rank(pct=True).values
    o = o.rename(columns={'Suppressive': 'score'})
    pb = o.groupby(['patient', 'pt_bin']).agg(score=('score', 'mean'), n=('score', 'size')).reset_index()
    pb = pb[pb.n >= 3].copy(); pb['pt'] = pb.pt_bin.map(MID)
    print(f'cells {len(o)}, patients {o.patient.nunique()}; patient-bin points {len(pb)} from {pb.patient.nunique()} patients')

    rows = []
    def add(analysis, level, term, est=np.nan, lo=np.nan, hi=np.nan, stat=np.nan, p=np.nan, n=np.nan, note=''):
        rows.append(dict(analysis=analysis, level=level, term=term, estimate=est, ci_low=lo, ci_high=hi,
                         statistic=stat, p_value=p, n=n, note=note))

    def mixed(formula, data, re='~pt', reml=True, strict=False):
        """MixedLM with random intercept + slope; falls back to random intercept only if the slope model fails (not when strict)."""
        for r in ([re, '~1'] if (re != '~1' and not strict) else [re]):
            try:
                m = smf.mixedlm(formula, data, groups=data['patient'], re_formula=r).fit(reml=reml, method=['lbfgs', 'powell'])
                if np.isfinite(m.llf) and (strict or np.all(np.isfinite(m.bse[m.fe_params.index]))):
                    return m, r   # strict (LRT) fits need only the likelihood
            except Exception as e:
                print(f'  MixedLM failed: {formula} re={r}: {e!r}')
        raise RuntimeError(formula)

    def lmm_block(data, level):
        m, r = mixed('score ~ pt', data)
        ci = m.conf_int().loc['pt']
        add('LMM linear slope', level, 'pt (score per full pseudotime range 0-1)', m.params['pt'], ci[0], ci[1],
            m.tvalues['pt'], m.pvalues['pt'], len(data), f'score ~ pt + (1 + pt | patient), REML; random effects {r}; '
        f'{data.patient.nunique()} patients')
        return m, r

    def lrt(data, level, f1, label, df_note, re):
        m0, _ = mixed('score ~ pt', data, re=re, reml=False, strict=True)
        m1, _ = mixed(f1, data, re=re, reml=False, strict=True)
        stat = 2 * (m1.llf - m0.llf); ddf = len(m1.fe_params) - len(m0.fe_params); p = chi2.sf(max(stat, 0), ddf)
        add('Non-linearity LRT', level, f'linear vs {label}', np.nan, np.nan, np.nan, stat, p, len(data),
            f'ML fits, same random effects ({re}); chi2 df = {ddf}; {df_note}')
        return m1, stat, p

    # ---------------- 1. linear trend ----------------
    m_pb, re_pb = lmm_block(pb, 'patient-bin')
    m_cell, re_cell = lmm_block(o, 'cell')
    rs, ps = spearmanr(o.pt, o.score); add('Spearman', 'cell', 'rho(pt, score)', rs, stat=rs, p=ps, n=len(o), note='pooled cells')
    rs2, ps2 = spearmanr(pb.pt, pb.score); add('Spearman', 'patient-bin', 'rho(bin midpoint, score)', rs2, stat=rs2, p=ps2,
                                                n=len(pb), note='pooled patient-bin means (ties: 5 x values)')
    # quadratic coefficient direction (convex = decelerating decline)
    m_q, _ = mixed('score ~ pt + I(pt**2)', o, re=re_cell)
    ciq = m_q.conf_int().loc['I(pt ** 2)']
    add('LMM quadratic term', 'cell', 'pt^2', m_q.params['I(pt ** 2)'], ciq[0], ciq[1], m_q.tvalues['I(pt ** 2)'],
        m_q.pvalues['I(pt ** 2)'], len(o), 'score ~ pt + pt^2 + (1 + pt | patient), REML; positive = decline flattens')
    m_q2, _ = mixed('score ~ pt + I(pt**2)', pb, re=re_pb)
    ciq2 = m_q2.conf_int().loc['I(pt ** 2)']
    add('LMM quadratic term', 'patient-bin', 'pt^2', m_q2.params['I(pt ** 2)'], ciq2[0], ciq2[1], m_q2.tvalues['I(pt ** 2)'],
        m_q2.pvalues['I(pt ** 2)'], len(pb), 'score ~ pt + pt^2 + (1 + pt | patient), REML; positive = decline flattens')

    # ---------------- 2. non-linearity ----------------
    _, lr_q_cell, p_q_cell = lrt(o, 'cell', 'score ~ pt + I(pt**2)', 'quadratic', '', re_cell)
    _, lr_s_cell, p_s_cell = lrt(o, 'cell', 'score ~ bs(pt, df=4)', 'cubic B-spline (df = 4)', 'patsy bs basis; spans the linear model', re_cell)
    _, lr_q_pb, p_q_pb = lrt(pb, 'patient-bin', 'score ~ pt + I(pt**2)', 'quadratic', '', re_pb)
    _, lr_c_pb, p_c_pb = lrt(pb, 'patient-bin', 'score ~ C(pt_bin)', 'bin as category (lack of fit)',
                             'C(pt_bin) = free mean per bin, contains the linear model; tests any departure from linearity', re_pb)
    # slope after bin 0 (does the decline continue past the first step?)
    for lab, sub_pb, sub_c in [('bins 0-1 only', pb[pb.pt_bin <= 1], o[o.pt_bin <= 1]), ('bins 1-4 only', pb[pb.pt_bin >= 1], o[o.pt_bin >= 1]),
                             ('bins 2-4 only', pb[pb.pt_bin >= 2], o[o.pt_bin >= 2])]:
        for level, dat in [('patient-bin', sub_pb), ('cell', sub_c)]:
            try:
                m, r = mixed('score ~ pt', dat)
                ci = m.conf_int().loc['pt']
                add(f'LMM slope, {lab}', level, 'pt (score per full pseudotime range 0-1)', m.params['pt'], ci[0], ci[1],
                    m.tvalues['pt'], m.pvalues['pt'], len(dat), f'random effects {r}; {dat.patient.nunique()} patients')
            except RuntimeError:
                add(f'LMM slope, {lab}', level, 'pt', note='did not converge')
    # LOWESS shape (cell level) + patient-cluster bootstrap
    FRAC = 0.3
    grid = np.linspace(0.02, 0.98, 97)
    def low(df):
        f = lowess(df.score.values, df.pt.values, frac=FRAC, it=1, return_sorted=True)
        xs, idx = np.unique(f[:, 0], return_index=True)
        return np.interp(grid, xs, f[idx, 1])
    lw = low(o)
    pats = o.patient.unique(); groups = {p: g for p, g in o.groupby('patient')}
    boot = []
    for _ in range(500):
        s = rng.choice(pats, len(pats), replace=True)
        df = pd.concat([groups[p] for p in s], ignore_index=True)
        df['pt'] = df.pt + rng.uniform(-1e-6, 1e-6, len(df))   # break duplicate-x ties from resampled patients
        boot.append(low(df))
    boot = np.array(boot); blo, bhi = np.percentile(boot, [2.5, 97.5], axis=0)
    for b, x in MID.items():
        i = np.argmin(np.abs(grid - x))
        add('LOWESS', 'cell', f'fitted score at bin {b} midpoint (pt = {x:.1f})', lw[i], blo[i], bhi[i], n=len(o),
            note=f'frac = {FRAC}; 95% CI = 500 patient-cluster bootstrap resamples')
    d = np.diff(lw)
    i_min = int(np.argmin(lw))
    add('LOWESS', 'cell', 'drop pt 0.1 -> 0.3 (bin 0 -> 1)', lw[np.argmin(np.abs(grid - .3))] - lw[np.argmin(np.abs(grid - .1))], n=len(o))
    add('LOWESS', 'cell', 'drop pt 0.3 -> 0.9 (bin 1 -> 4)', lw[np.argmin(np.abs(grid - .9))] - lw[np.argmin(np.abs(grid - .3))], n=len(o))
    add('LOWESS', 'cell', 'fraction of pseudotime grid where curve is decreasing', float((d < 0).mean()), n=len(o),
        note='97-point grid 0.02-0.98')
    add('LOWESS', 'cell', 'pseudotime of LOWESS minimum', grid[i_min], n=len(o),
        note=f'rise from minimum to end = {lw[-1] - lw[i_min]:.3g}')
    boot_mono = float(np.mean([(np.diff(b_) <= 0).mean() for b_ in boot]))
    add('LOWESS', 'cell', 'bootstrap: mean fraction of grid decreasing', boot_mono, n=500)

    # ---------------- 3. per-patient slopes ----------------
    sl = []
    for p, g in pb.groupby('patient'):
        if g.pt_bin.nunique() >= 3:
            r = linregress(g.pt, g.score); sl.append(dict(patient=p, slope=r.slope, n_bins=g.pt_bin.nunique(), n_cells=g.n.sum(),
                                                             rho=spearmanr(g.pt_bin, g.score)[0]))
    sl = pd.DataFrame(sl)
    neg = int((sl.slope < 0).sum()); sign_p = binomtest(neg, len(sl), 0.5).pvalue
    add('Per-patient OLS slope', 'patient', 'n patients with negative slope', neg, stat=neg, p=sign_p, n=len(sl),
        note='patient-bin means vs bin midpoint; patients with >=3 bins (>=3 Tregs each); two-sided sign test')
    add('Per-patient OLS slope', 'patient', 'median slope', sl.slope.median(), sl.slope.quantile(.25), sl.slope.quantile(.75), n=len(sl),
        note='ci_low/ci_high = IQR')
    add('Per-patient Spearman', 'patient', 'n patients with rho < 0', int((sl.rho < 0).sum()), n=len(sl))
    for _, r in sl.sort_values('slope').iterrows():
        add('Per-patient OLS slope', 'patient', r.patient, r.slope, n=int(r.n_cells), note=f'{int(r.n_bins)} bins; Spearman rho = {r.rho:.3g}')
    re_slopes = pd.Series({k: m_cell.fe_params['pt'] + v.get('pt', 0) for k, v in m_cell.random_effects.items()})
    add('LMM BLUP slope', 'cell', 'n patients with negative patient-specific slope', int((re_slopes < 0).sum()), n=len(re_slopes),
        note='fixed-effect + random slope from the cell-level LMM; all 65 patients (shrunken)')

    res = pd.DataFrame(rows)
    def sig3(v):
        if isinstance(v, (float, np.floating)) and np.isfinite(v): return float(f'{v:.3g}')
        return v
    for c in ['estimate', 'ci_low', 'ci_high', 'statistic', 'p_value']:
        res[c] = res[c].map(sig3)
    res.to_csv(f'{STEM}.csv', index=False)
    print(res.drop(columns='note').head(40).to_string())

    # ---------------- figure ----------------
    def fmt_p(p):
        if p < 1e-3:
            m, e = f'{p:.1e}'.split('e'); return f'P = {m} × 10$^{{{int(e)}}}$'
        return f'P = {p:.2g}'
    plt.rcParams.update({'font.size': 14, 'axes.labelsize': 16, 'axes.titlesize': 16, 'xtick.labelsize': 14,
                         'ytick.labelsize': 14, 'font.family': 'DejaVu Sans'})
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    xs = np.linspace(0.05, 0.95, 100)
    fig, bx = plt.subplots(figsize=(6.5, 5))
    light = mcolors.LinearSegmentedColormap.from_list('l', ['#f3eef1', '#b9a3b5', '#5e4a66'])
    hb = bx.hexbin(o.pt, o.score, gridsize=(22, 18), cmap=light, mincnt=1, bins='log', linewidths=0.2, edgecolors='white')
    bx.fill_between(grid, blo, bhi, color='#e07b39', alpha=0.45, lw=0, zorder=3)
    bx.plot(grid, lw, color='#e07b39', lw=3.2, zorder=4)
    fe_c = m_cell.fe_params[['Intercept', 'pt']].values
    bx.plot(xs, fe_c[0] + fe_c[1] * xs, color='black', lw=2.2, ls='--', zorder=4)
    bx.set_xlim(0, 1); bx.set_xticks([0, 0.25, 0.5, 0.75, 1]); bx.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
    bx.set_ylim(o.score.min() - 0.05, o.score.max() + 0.05)
    bx.set_xlabel('Treg Pseudotime', labelpad=2); bx.set_ylabel('Suppressive Score', labelpad=2)
    bx.spines[['top', 'right']].set_visible(False)
    # statistics above the axes (outside the data area); title above them
    bx.text(0.5, 1.02, f'Spearman ρ = {rs:.2f}, {fmt_p(ps)}\nNon-Linearity {fmt_p(p_s_cell)}  (n = {len(o):,} Tregs)', transform=bx.transAxes,
            ha='center', va='bottom', fontsize=13.5, linespacing=1.3)
    bx.set_title('Suppressive Score Along Treg Pseudotime', pad=56, fontsize=16)
    cb = fig.colorbar(hb, ax=bx, fraction=0.05, pad=0.02); cb.set_label('Cells per Hexagon', fontsize=14); cb.ax.tick_params(labelsize=12)
    hl = [(Line2D([], [], color='#e07b39', lw=3.2), Patch(color='#e07b39', alpha=0.45, lw=0)),
          Line2D([], [], color='black', lw=2.2, ls='--')]
    fig.legend(hl, ['LOWESS (95% CI)', 'Linear Mixed Model'], loc='upper center', bbox_to_anchor=(0.47, 0.03), ncol=2,
               fontsize=14, frameon=False, handlelength=1.8, columnspacing=1.4)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)

    # ---------------- legend ----------------
    g = lambda a, t_, lvl: res[(res.analysis == a) & (res.term.str.startswith(t_)) & (res.level == lvl)].iloc[0]
    s1 = g('LMM slope, bins 1-4 only', 'pt', 'patient-bin'); s2 = g('LMM slope, bins 2-4 only', 'pt', 'patient-bin')
    legend = f"""Fig 4E (right). Suppressive Score Along Treg Pseudotime.
Suppressive gene score (FOXP3, CTLA4, IL2RA, BATF, CD27; scanpy score_genes on log-normalized expression) versus diffusion
pseudotime, rank-normalized to 0-1, for {len(o):,} Tregs from {o.patient.nunique()} patients in the CD scRNA-seq atlas, shown as a hexagonal
density (log colour scale, cells per hexagon). Orange line and band: LOWESS fit (span {FRAC}) with 95% CI from 500 bootstrap
resamples of patients. Dashed line: fixed-effect fit of a linear mixed model with random intercept and slope per patient
(slope = {m_cell.params['pt']:.2f} per unit pseudotime, P = {m_cell.pvalues['pt']:.2g}). Spearman rho over all cells = {rs:.2f} (P = {ps:.2g}). Non-linearity P:
likelihood-ratio test of a cubic B-spline (df = 4) versus the linear mixed model (P = {p_s_cell:.2g}; quadratic, P = {p_q_cell:.2g}). At the
patient level (per-patient mean score in each of five equally sized pseudotime bins with >=3 Tregs; {len(pb)} means from
{pb.patient.nunique()} patients) the linear mixed-model slope was {m_pb.params['pt']:.2f} (P = {m_pb.pvalues['pt']:.2g}), and the per-patient slope was negative in
{neg} of {len(sl)} patients with Tregs in >=3 bins (two-sided sign test, P = {sign_p:.2g}). The decline is steepest between the first two
bins and flattens thereafter (bins 2-4 only: slope = {s2.estimate:.2f}, P = {s2.p_value:.2g}). Full statistics, including per-patient slopes:
Fig3E_suppression_score_vs_pseudotime.csv.
"""
    head_, body = legend.split('\n', 1)
    legend = head_ + '\n' + ' '.join(body.split()) + '\n'   # one paragraph
    with open(f'{STEM}_legend.txt', 'w') as fh:
        fh.write(legend)
    print(legend)
