#!/usr/bin/env python
"""CD scRNA-seq atlas: Treg selection, diffusion pseudotime and the Treg pseudotime panels (Fig 4D-F).

Inputs   CD atlas h5ad (09_atlas_h5ad_conversion.ipynb); control/inactive/active object without cycling cells
         (ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad, obs used as is); treg_ptbin_markers.xlsx (Fig 4F genes)
Outputs  cd_treg_rebuild_rowfix/objects/: Treg input objects, treg_pseudotime.h5ad and treg_pseudotime_bins.csv
         (read by 13-15, cosmx_cd and treatment_validation); Fig3D-I_treg_pseudotime_signals/: panels and CSVs
Sections treg_inputs -> treg_selection_dpt -> doublet_branch_dpt -> F4DF_module_scores_markers -> F4E_pseudotime_umap
         -> F4E_suppression_scatter   (RUN_SECTIONS=<id>,<id> runs a subset)
"""
import os as _os
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
FIG = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals'

# ---- Section treg_inputs: control/inactive/active atlas without cycling cells (X rows from the atlas) and its Treg subset
if _run('treg_inputs'):
    import sys, warnings; warnings.filterwarnings('ignore')
    import numpy as np, anndata as ad
    M = '/share/fsmresfiles/UC/scRNA-seq/merged_cd/'
    atlas = ad.read_h5ad(M + 'cd_atlas_annotated.h5ad')
    der = ad.read_h5ad(M + 'ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad')
    pos = atlas.obs_names.get_indexer(der.obs_names.values.astype(str))
    full = ad.AnnData(X=atlas.X[pos], obs=der.obs.copy(), var=der.var.copy())
    del atlas
    treg = np.where(der.obs['cell type'].astype(str).values == 'Treg')[0]
    full[treg].copy().write_h5ad(f'{D}/objects/tregs_input.h5ad')
    if '--skip-full' not in sys.argv:   # --skip-full writes only the Treg subset
        full.write_h5ad(f'{D}/objects/ctrl_inact_act_sccoda_no_cycle_correct_cat_adata.h5ad')


# ---- Section treg_selection_dpt: doublet removal, two filtering rounds, strict Treg cluster, PAGA, root and DPT
if _run('treg_selection_dpt'):
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, scanpy.external as sce, anndata as ad

    def embed(a, resolutions):
        sc.pp.highly_variable_genes(a, n_top_genes=3000, flavor='seurat_v3', layer='counts')
        sc.pp.scale(a, max_value=10)
        sc.tl.pca(a, n_comps=50, use_highly_variable=True)
        sce.pp.harmony_integrate(a, key='dataset')
        sc.pp.neighbors(a, use_rep='X_pca_harmony', n_neighbors=20, n_pcs=30)
        sc.tl.umap(a)
        for r in resolutions:
            sc.tl.leiden(a, resolution=r, key_added=f'leiden_{r}')

    def score_panels(a, first_name):
        panels = {first_name: ['BATF', 'PRDM1', 'MAF', 'TNFRSF18', 'TNFRSF4'],
                  'Tfr': ['BCL6', 'CXCR5', 'PDCD1', 'ICOS'],
                  'Th17_like': ['RORC', 'CCR6', 'IL17A', 'IL17F'],
                  'IFN_response': ['ISG15', 'IFI6', 'MX1', 'IFIT3'],
                  'Cell_cycle': ['MKI67', 'TOP2A', 'STMN1'],
                  'Treg_core': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT', 'ENTPD1', 'LRRC32', 'IL10']}
        for name, genes in panels.items():
            sc.tl.score_genes(a, [g for g in genes if g in a.var_names], score_name=name)

    # Treg subset; scrublet on raw counts, doublet score > 0.2 removed
    adata = ad.read_h5ad(f'{D}/objects/tregs_input.h5ad')
    adata.raw = adata.copy()
    adata_raw = adata.raw.to_adata()
    adata.raw = None   # module scores are computed on the scaled X
    adata_raw.X = pd.DataFrame(adata_raw.X.toarray()).astype(int).values
    sc.pp.scrublet(adata_raw, expected_doublet_rate=0.1, n_prin_comps=50)
    adata_raw.obs['predicted_doublet_0.2'] = [0 if i <= 0.2 else 1 for i in adata_raw.obs['doublet_score'].tolist()]
    doublet = adata_raw[adata_raw.obs['predicted_doublet_0.2'] == 1].obs_names
    adata = adata[~adata.obs_names.isin(doublet)]
    adata.layers['counts'] = adata_raw[adata_raw.obs_names.isin(adata.obs_names.tolist())].X.copy()

    # normalise, log1p, embed, module scores
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    embed(adata, [0.1, 0.2, 0.3])
    score_panels(adata, 'Activated_tissue')

    # round 1: 0 < Treg_core < 2.5 and Cell_cycle < 2.5
    adata_clean = adata[adata.obs['Treg_core'] > 0].copy()
    adata_clean = adata_clean[adata_clean.obs['Cell_cycle'] < 2.5].copy()
    adata_clean = adata_clean[adata_clean.obs['Treg_core'] < 2.5].copy()
    embed(adata_clean, [0.1, 0.2, 0.3, 0.4])

    # round 2: drop leiden_0.3 clusters with median Treg_core <= cell q0.3 or median Cell_cycle >= cell q0.7
    cluster_stats = adata_clean.obs.groupby('leiden_0.3', observed=False)[['Treg_core', 'Cell_cycle']].median()
    bad_clusters = cluster_stats.index[
        (cluster_stats['Treg_core'] <= adata_clean.obs['Treg_core'].quantile(0.3)) |
        (cluster_stats['Cell_cycle'] >= adata_clean.obs['Cell_cycle'].quantile(0.7))]
    adata_clean2 = adata_clean[~adata_clean.obs['leiden_0.3'].isin(bad_clusters)].copy()
    embed(adata_clean2, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    score_panels(adata_clean2, 'Suppresive')
    adata_tmp = adata_clean2.copy()
    adata_tmp.X = adata_clean2.layers['counts'].copy()
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    adata_clean2.layers['log1p'] = adata_tmp.X.copy()
    del adata_tmp

    # strict Treg cluster: highest median Treg_core among leiden_0.1 clusters with median Cell_cycle < cell q0.7
    o2 = adata_clean2.obs
    cs = o2.groupby('leiden_0.1', observed=True)[['Treg_core', 'Cell_cycle']].median().add_suffix('_median')
    elig = cs[cs['Cell_cycle_median'] < o2['Cell_cycle'].quantile(0.7)]
    if len(elig) == 0:
        elig = cs
    strict_cl = elig['Treg_core_median'].idxmax()
    clus0 = adata_clean2[adata_clean2.obs['leiden_0.1'] == strict_cl]
    print('strict cluster', strict_cl, '| cells', clus0.n_obs)

    # PAGA-initialised UMAP; root = first cell of the leiden_0.6 sub-cluster (>= 20 cells) with highest median Treg_core
    sc.tl.paga(clus0, groups='leiden_0.6')
    sc.pl.paga(clus0, plot=False)
    sc.tl.umap(clus0, init_pos='paga')
    clus0.obs['treg_stability_score'] = clus0.obs['Activated_tissue'] - clus0.obs['Th17_like']
    clus0.obs['treg_stability_score2'] = clus0.obs['Treg_core'] - clus0.obs['Th17_like']
    cs6 = clus0.obs.groupby('leiden_0.6', observed=True)['Treg_core'].median().to_frame('Treg_core_median')
    cs6['n'] = clus0.obs['leiden_0.6'].value_counts()
    root_cl = cs6[cs6['n'] >= 20]['Treg_core_median'].idxmax()
    clus0.uns['iroot'] = np.flatnonzero(clus0.obs['leiden_0.6'] == root_cl)[0]
    clus0.uns['root_subcluster'] = str(root_cl); clus0.uns['strict_cluster_leiden01'] = str(strict_cl)
    sc.tl.dpt(clus0)   # on the round-2 kNN graph (Harmony, k=20, 30 PCs)
    sc.pp.neighbors(clus0, n_neighbors=8)   # display graph, built after the DPT
    clus0.obs['pt_bin'] = pd.qcut(clus0.obs['dpt_pseudotime'], q=5, labels=False).astype(str)
    clus0.obs['Suppressive'] = clus0.obs['Activated_tissue'].tolist().copy()
    clus0.copy().write_h5ad(f'{D}/objects/treg_pseudotime_raw.h5ad')


# ---- Section doublet_branch_dpt: drop a small separated high-pseudotime branch (largest DPT gap > 0.1 and < 5% of
# cells; none with the current object), then recompute the Treg DPT with the same root rule
if _run('doublet_branch_dpt'):
    import anndata as ad, scanpy as sc, numpy as np, pandas as pd
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime_raw.h5ad')
    pt = np.sort(t.obs['dpt_pseudotime'].values)
    gaps = np.diff(pt); g = int(np.argmax(gaps)); cut = (pt[g] + pt[g + 1]) / 2
    branch = t.obs['dpt_pseudotime'] > cut
    if not (gaps[g] > 0.1 and branch.mean() < 0.05):
        branch[:] = False
    print(f'removing {int(branch.sum())} branch cells')
    c = t[~branch.values].copy()
    c.obsm.pop('X_diffmap', None)
    c.uns.pop('diffmap_evals', None); c.uns.pop('iroot', None)
    sc.pp.neighbors(c, use_rep='X_pca_harmony', n_neighbors=20, n_pcs=30)
    cs6 = c.obs.groupby('leiden_0.6', observed=True)['Treg_core'].median().to_frame('Treg_core_median')
    cs6['n'] = c.obs['leiden_0.6'].value_counts()
    root_cl = cs6[cs6['n'] >= 20]['Treg_core_median'].idxmax()
    c.uns['iroot'] = int(np.flatnonzero(c.obs['leiden_0.6'] == root_cl)[0]); c.uns['root_subcluster'] = str(root_cl)
    sc.tl.diffmap(c); sc.tl.dpt(c)
    c.obs['pt_bin'] = pd.qcut(c.obs['dpt_pseudotime'], q=5, labels=False)
    c.write_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    o = c.obs[['pt_bin', 'dpt_pseudotime', 'leiden_0.6', 'patient', 'sample', 'condition_summary', 'dataset']].copy()
    o['dpt_pseudotime_rank'] = o['dpt_pseudotime'].rank(pct=True).round(4)
    o.to_csv(f'{D}/objects/treg_pseudotime_bins.csv')


# ---- Section F4DF_module_scores_markers: module scores on the Treg UMAP (Fig 4D); markers by pseudotime bin (Fig 4F)
if _run('F4DF_module_scores_markers'):
    import anndata as ad, scanpy as sc, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.colors as mcolors
    from matplotlib.lines import Line2D
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    panels = {'Suppressive': ['BATF', 'PRDM1', 'MAF', 'TNFRSF18', 'TNFRSF4'], 'Tfr': ['BCL6', 'CXCR5', 'PDCD1', 'ICOS'],
              'Th17_like': ['RORC', 'CCR6', 'IL17A', 'IL17F']}
    for k, g in panels.items():
        sc.tl.score_genes(t, [x for x in g if x in t.var_names], score_name=k)
    p = t.copy(); sc.pp.neighbors(p, n_neighbors=8, use_rep='X_pca_harmony')   # graph used only to draw edges
    # Fig 4D
    vmax = np.percentile(np.concatenate([p.obs[k].dropna() for k in panels]), 99)
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 4.0), gridspec_kw={'wspace': 0.0})
    for a, k in zip(ax, panels):
        sc.pl.umap(p, color=k, edges=True, edges_width=0.05, edges_color='#d9d9d9', cmap='Reds', vmin=0, vmax=vmax,
                   frameon=False, na_color='#f2f2f2', ax=a, show=False, colorbar_loc=None, s=30)
        a.set_title(k, fontsize=16)
    fig.subplots_adjust(wspace=0.0, left=0.0, right=0.93)
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=mcolors.Normalize(0, vmax)); fig.colorbar(sm, cax=fig.add_axes([0.94, 0.15, 0.012, 0.65]))
    fig.savefig(f'{FIG}/fig3D_module_scores.png', dpi=300, bbox_inches='tight')
    # Fig 4F: Wilcoxon per pseudotime bin (non-Control Tregs, log1p); curated genes minus KLRD1 and ARID5B
    f = t[t.obs['condition_summary'].astype(str) != 'Control'].copy(); f.X = f.layers['log1p']
    f.obs['pt_bin'] = f.obs['pt_bin'].astype(int).astype(str).astype('category')
    sc.tl.rank_genes_groups(f, groupby='pt_bin', method='wilcoxon', pts=True)
    res = sc.get.rank_genes_groups_df(f, group=None); res['group'] = res['group'].astype(str)
    genes = pd.read_excel('/share/fsmresfiles/UC/scRNA-seq/merged_cd/treg_ptbin_markers.xlsx').query("plot == 'plot'")['names'].unique()
    genes = [g for g in genes if g not in ('KLRD1', 'ARID5B')]
    cap = 2.0   # colour saturates at +/-cap
    sub = res[res.names.isin(genes)]
    lfc = sub.pivot(index='group', columns='names', values='logfoldchanges').reindex(['0', '1', '2', '3', '4'])
    padj = sub.pivot(index='group', columns='names', values='pvals_adj').reindex(['0', '1', '2', '3', '4'])
    pk = lfc.fillna(0).idxmax(axis=0); order = sorted(lfc.columns, key=lambda g: (int(pk[g]), -lfc[g].max()))
    lfc, padj = lfc[order], padj[order]
    cmap = mcolors.LinearSegmentedColormap.from_list('logfc', ['#0000FF', '#FFFFFF', '#FF0000'])
    norm = mcolors.Normalize(-cap, cap)
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
    fig.savefig(f'{FIG}/fig3F_marker_dotplot.png', dpi=300, bbox_inches='tight')
    out = lfc.T.round(3).add_prefix('log2FC_bin').join(padj.T.apply(lambda c: c.map(lambda z: float('%.3g' % z))).add_prefix('padj_bin'))
    out.to_csv(f'{FIG}/fig3F_marker_dotplot.csv')


# ---- Section F4E_pseudotime_umap: Treg UMAP coloured by rank-normalised DPT (Fig 4E left)
if _run('F4E_pseudotime_umap'):
    import anndata as ad, scanpy as sc, matplotlib.pyplot as plt, matplotlib.colors as mcolors
    purple = mcolors.LinearSegmentedColormap.from_list('p', ['#e8d4cf', '#a16f8d', '#2d223a'])
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    t.obs['Pseudotime'] = t.obs['dpt_pseudotime'].rank(pct=True).values
    fig, ax = plt.subplots(figsize=(6.1, 4.6))
    sc.pl.umap(t, color='Pseudotime', cmap=purple, vmin=0, vmax=1, frameon=False, s=22, ax=ax, show=False)
    ax.set_title('Diffusion Pseudotime', fontsize=15)
    fig.savefig(f'{FIG}/fig3E_pseudotime_umap.png', dpi=300, bbox_inches='tight'); plt.close(fig)


# ---- Section F4E_suppression_scatter: suppressive score (FOXP3, CTLA4, IL2RA, BATF, CD27) vs rank-normalised DPT;
# LMMs, Spearman, non-linearity LRTs, LOWESS with patient bootstrap, per-patient slopes + sign test (Fig 4E right)
if _run('F4E_suppression_scatter'):
    import warnings; warnings.filterwarnings('ignore')
    import anndata as ad, scanpy as sc, pandas as pd, numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt, matplotlib.colors as mcolors
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import statsmodels.formula.api as smf
    from statsmodels.nonparametric.smoothers_lowess import lowess
    from scipy.stats import spearmanr, binomtest, chi2, linregress

    STEM = f'{FIG}/Fig3E_suppression_score_vs_pseudotime'
    MID = {b: 0.1 + 0.2 * b for b in range(5)}
    rng = np.random.default_rng(0)

    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    t.X = t.layers['log1p']
    sc.tl.score_genes(t, ['FOXP3', 'CTLA4', 'IL2RA', 'BATF', 'CD27'], score_name='Suppressive')
    o = t.obs[['patient', 'pt_bin', 'dpt_pseudotime', 'Suppressive']].copy()
    o['patient'] = o['patient'].astype(str); o['pt_bin'] = o['pt_bin'].astype(int)
    o['pt'] = o['dpt_pseudotime'].rank(pct=True).values
    o = o.rename(columns={'Suppressive': 'score'})
    pb = o.groupby(['patient', 'pt_bin']).agg(score=('score', 'mean'), n=('score', 'size')).reset_index()
    pb = pb[pb.n >= 3].copy(); pb['pt'] = pb.pt_bin.map(MID)   # patient-bin means, bins with >= 3 Tregs

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
                    return m, r
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

    # linear trend
    m_pb, re_pb = lmm_block(pb, 'patient-bin')
    m_cell, re_cell = lmm_block(o, 'cell')
    rs, ps = spearmanr(o.pt, o.score); add('Spearman', 'cell', 'rho(pt, score)', rs, stat=rs, p=ps, n=len(o), note='pooled cells')
    rs2, ps2 = spearmanr(pb.pt, pb.score); add('Spearman', 'patient-bin', 'rho(bin midpoint, score)', rs2, stat=rs2, p=ps2,
                                                n=len(pb), note='pooled patient-bin means (ties: 5 x values)')
    for lvl, dat, re_ in [('cell', o, re_cell), ('patient-bin', pb, re_pb)]:
        m_q, _ = mixed('score ~ pt + I(pt**2)', dat, re=re_)
        ciq = m_q.conf_int().loc['I(pt ** 2)']
        add('LMM quadratic term', lvl, 'pt^2', m_q.params['I(pt ** 2)'], ciq[0], ciq[1], m_q.tvalues['I(pt ** 2)'],
            m_q.pvalues['I(pt ** 2)'], len(dat), 'score ~ pt + pt^2 + (1 + pt | patient), REML; positive = decline flattens')

    # non-linearity
    lrt(o, 'cell', 'score ~ pt + I(pt**2)', 'quadratic', '', re_cell)
    _, _, p_s_cell = lrt(o, 'cell', 'score ~ bs(pt, df=4)', 'cubic B-spline (df = 4)', 'patsy bs basis; spans the linear model', re_cell)
    lrt(pb, 'patient-bin', 'score ~ pt + I(pt**2)', 'quadratic', '', re_pb)
    lrt(pb, 'patient-bin', 'score ~ C(pt_bin)', 'bin as category (lack of fit)',
        'C(pt_bin) = free mean per bin, contains the linear model; tests any departure from linearity', re_pb)
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
    # LOWESS (cell level) with 500 patient-cluster bootstrap resamples
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
    add('LOWESS', 'cell', 'bootstrap: mean fraction of grid decreasing', float(np.mean([(np.diff(b_) <= 0).mean() for b_ in boot])), n=500)

    # per-patient slopes (patients with >= 3 bins), two-sided sign test
    sl = []
    for p, g in pb.groupby('patient'):
        if g.pt_bin.nunique() >= 3:
            r = linregress(g.pt, g.score); sl.append(dict(patient=p, slope=r.slope, n_bins=g.pt_bin.nunique(), n_cells=g.n.sum(),
                                                             rho=spearmanr(g.pt_bin, g.score)[0]))
    sl = pd.DataFrame(sl)
    neg = int((sl.slope < 0).sum())
    add('Per-patient OLS slope', 'patient', 'n patients with negative slope', neg, stat=neg, p=binomtest(neg, len(sl), 0.5).pvalue, n=len(sl),
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
    for c in ['estimate', 'ci_low', 'ci_high', 'statistic', 'p_value']:
        res[c] = res[c].map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    res.to_csv(f'{STEM}.csv', index=False)

    # figure: hexbin, LOWESS + 95% band, cell-level LMM fit
    def fmt_p(p):
        if p < 1e-3:
            m, e = f'{p:.1e}'.split('e'); return f'P = {m} × 10$^{{{int(e)}}}$'
        return f'P = {p:.2g}'
    plt.rcParams.update({'font.size': 14, 'axes.labelsize': 16, 'axes.titlesize': 16, 'xtick.labelsize': 14,
                         'ytick.labelsize': 14, 'font.family': 'DejaVu Sans'})
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
