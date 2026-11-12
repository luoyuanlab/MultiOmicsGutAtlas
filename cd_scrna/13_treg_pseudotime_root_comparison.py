#!/usr/bin/env python
"""CD scRNA-seq Tregs: Th17-like Tregs vs conventional Th17 cells, and pseudotime root robustness (Supp 8A-E).

Inputs   CD atlas h5ad (09_atlas_h5ad_conversion.ipynb); T-cell sub-clusters tcell_iter3.RDS and tcell_iter3_umap.csv
         (05_t_cell_subclustering.R); treg_pseudotime.h5ad (12_treg_pseudotime.py)
Outputs  SuppFig8A-I_treg_pseudotime_robustness/: panels and CSVs
Sections S8AB_th17like_tregs -> S8CD_leiden_root -> S8E_palantir_main_root -> S8E_root_comparison
         -> S8E_root_comparison_plot   (RUN_SECTIONS=<id>,<id> runs a subset)
"""
import os as _os
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
OUT = '/share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness'


def w3(df, name):
    df.to_csv(f'{OUT}/{name}', index=False, float_format='%.3g')


# ---- Section S8AB_th17like_tregs: Th17-like Tregs (top 10% Th17_like score) vs conventional Th17 cells: position on the
# Fig 4B T-cell UMAP (kNN k = 30 in its Harmony space), markers (Fisher, BH), Scrublet doublets, Fig 4E trend without them
if _run('S8AB_th17like_tregs'):
    import os, subprocess, warnings
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
    from sklearn.neighbors import NearestNeighbors
    from scipy.stats import spearmanr, binomtest, friedmanchisquare, mannwhitneyu, fisher_exact
    from statsmodels.stats.multitest import multipletests
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    M = '/share/fsmresfiles/UC/scRNA-seq/merged_cd'
    CACHE = f'{D}/objects/intermediate/th17_check'
    os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
    GENES = ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TNFRSF18', 'RORC', 'IL17A', 'IL17F', 'IL23R', 'CCR6', 'KLRB1', 'IL22',
             'CD40LG', 'IL7R', 'ITGAE']
    EXTRA = ['CD4', 'CD8A', 'CD8B']
    TREG_LAB, TRM_LAB = 'Treg', 'CD4+ Trm'
    CD4_CLUSTERS = ['CD4+ Trm', 'CD4+ Trm(Prolif)', 'CD4+ Effector', 'Tfh', 'Exhausted CD4+ Tfh', 'Resting T']

    # Fig 4B embedding (Seurat Harmony, exported once to a cache) and UMAP
    harm_f = f'{CACHE}/tcell_iter3_harmony.csv.gz'
    if not os.path.exists(harm_f):
        rcode = f"""suppressPackageStartupMessages(library(Seurat))
s <- readRDS('{M}/t/tcell_iter3.RDS')
e <- Embeddings(s, 'harmony_tcell_iter3'); g <- gzfile('{harm_f}', 'w'); write.csv(e, g); close(g)"""
        subprocess.run(['/usr/bin/Rscript', '-e', rcode], check=True)
    H = pd.read_csv(harm_f, index_col=0)
    U = pd.read_csv(f'{M}/t/tcell_iter3_umap.csv').set_index('cell')
    U = U.rename(columns={'umapharmonytcelliter3_1': 'umap1', 'umapharmonytcelliter3_2': 'umap2'})
    lab = U['iter3_anno'].astype(str)

    # expression from the CD atlas, in Fig 4B row order
    A = ad.read_h5ad(f'{M}/cd_atlas_annotated.h5ad', backed='r')
    idx = A.obs_names.get_indexer(U.index)
    order = np.argsort(idx)
    T = A[idx[order]].to_memory(); T = T[np.argsort(order)].copy()
    del A
    T.layers['counts'] = T.X.copy()
    T.obs['n_umi'] = np.asarray(T.X.sum(1)).ravel(); T.obs['n_genes'] = np.diff(sp.csr_matrix(T.X).indptr)
    sc.pp.normalize_total(T, target_sum=1e4); sc.pp.log1p(T)
    G = GENES + EXTRA
    gi = [T.var_names.get_loc(g) for g in G]
    C = pd.DataFrame(T.layers['counts'][:, gi].toarray(), index=T.obs_names, columns=G)
    L = pd.DataFrame(T.X[:, gi].toarray(), index=T.obs_names, columns=G)

    # groups
    tp = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    thr = tp.obs['Th17_like'].quantile(0.9)
    th17_strict = tp.obs_names[tp.obs['Th17_like'] > thr]
    other_strict = tp.obs_names[tp.obs['Th17_like'] <= thr]
    rorc_treg = tp.obs_names[C.loc[tp.obs_names, 'RORC'].values > 0]
    cd4 = lab.isin(CD4_CLUSTERS)
    conv_th17 = U.index[cd4.values & (C['FOXP3'].values == 0) & ((C['IL17A'].values > 0) | (C['IL17F'].values > 0))]
    groups = {'Th17-like Treg (Fig 4D/E Set, Top 10%)': th17_strict,
              'Other Treg (Fig 4D/E Set)': other_strict,
              'RORC+ Treg (Fig 4D/E Set)': rorc_treg,
              'Conventional Th17 (FOXP3- IL17A/F+ CD4)': conv_th17}
    _fx = C.loc[th17_strict, 'FOXP3'].values > 0
    _il = ((C.loc[th17_strict, 'IL17A'] > 0) | (C.loc[th17_strict, 'IL17F'] > 0)).values
    groups['Th17-like Treg, FOXP3+ (Fig 4D/E Set)'] = th17_strict[_fx]
    groups['Th17-like Treg, FOXP3- IL17A/F+ (Fig 4D/E Set)'] = th17_strict[~_fx & _il]
    for k in sorted(lab.unique()):
        groups[f'Cluster: {k}'] = U.index[lab.values == k]

    # (1) position: centroids, Treg -> CD4+ Trm axis, kNN label composition (Harmony and UMAP)
    XY = U[['umap1', 'umap2']]
    cen = XY.groupby(lab).mean()
    c_treg, c_trm = cen.loc[TREG_LAB].values, cen.loc[TRM_LAB].values
    axis = c_trm - c_treg
    k = 30
    labs = np.array(sorted(lab.unique()))
    lab_code = pd.Categorical(lab, categories=labs).codes
    NB_LABS = (TREG_LAB, TRM_LAB, 'CD4+ Effector', 'CD4+ Trm(Prolif)', 'Tfh', 'Exhausted CD4+ Tfh')
    nn_h = NearestNeighbors(n_neighbors=k + 1).fit(H.values)
    nn_u = NearestNeighbors(n_neighbors=k + 1).fit(XY.values)
    pos_rows = []
    for g, cells in groups.items():
        ii = U.index.get_indexer(cells)
        xy = XY.values[ii]; m = xy.mean(0)
        proj = ((xy - c_treg) @ axis) / (axis @ axis)              # 0 = Treg centroid, 1 = CD4+ Trm centroid
        d_treg_cell = np.linalg.norm(xy - c_treg, axis=1); d_trm_cell = np.linalg.norm(xy - c_trm, axis=1)
        r = {'group': g, 'n_cells': len(ii), 'n_patients': T.obs.loc[cells, 'patient'].nunique(),
             'umap1_centroid': m[0], 'umap2_centroid': m[1],
             'dist_centroid_to_Treg_centroid': np.linalg.norm(m - c_treg), 'dist_centroid_to_CD4Trm_centroid': np.linalg.norm(m - c_trm),
             'median_cell_axis_position_Treg0_CD4Trm1': np.median(proj),
             'pct_cells_axis_position_0.25_to_0.75': 100 * np.mean((proj > 0.25) & (proj < 0.75)),
             'pct_cells_closer_to_CD4Trm_than_Treg_centroid': 100 * np.mean(d_trm_cell < d_treg_cell),
             'pct_cells_labelled_Treg_in_4B': 100 * np.mean(lab.values[ii] == TREG_LAB)}
        for space, nn, Z in [('harmony', nn_h, H.values), ('umap', nn_u, XY.values)]:
            nb = nn.kneighbors(Z[ii], return_distance=False)[:, 1:]
            comp = np.stack([(lab_code[nb] == j).mean(1) for j in range(len(labs))], 1).mean(0)
            for j, l in enumerate(labs):
                if l in NB_LABS:
                    r[f'knn{k}_{space}_pct_{l}'] = 100 * comp[j]
            r[f'knn{k}_{space}_pct_other'] = 100 * (1 - sum(comp[j] for j, l in enumerate(labs) if l in NB_LABS))
            if space == 'harmony':
                r[f'knn{k}_harmony_pct_cells_majority_Treg'] = 100 * np.mean((lab_code[nb] == list(labs).index(TREG_LAB)).mean(1) > 0.5)
        pos_rows.append(r)
    pos = pd.DataFrame(pos_rows)
    w3(pos, 'th17_check_position_neighbours.csv')

    # (2) markers: % expressing and mean log1p CP10k
    mk_rows = []
    for g, cells in groups.items():
        r = {'group': g, 'n_cells': len(cells)}
        for gene in G:
            r[f'pct_{gene}'] = 100 * float((C.loc[cells, gene] > 0).mean())
        for gene in G:
            r[f'mean_log1p_{gene}'] = float(L.loc[cells, gene].mean())
        mk_rows.append(r)
    w3(pd.DataFrame(mk_rows), 'th17_check_markers_by_group.csv')

    # (3) doublets: Scrublet per sample on all T cells (cached), UMI/gene counts, FOXP3 x IL17A/F co-detection
    sfile = f'{CACHE}/scrublet_tcells_per_sample.csv.gz'
    if not os.path.exists(sfile):
        rows = []
        cnt = sp.csr_matrix(T.layers['counts'])
        for smp, ii in T.obs.groupby('sample', observed=True).indices.items():
            if len(ii) < 100:
                rows.append(pd.DataFrame({'doublet_score': np.nan}, index=T.obs_names[ii])); continue
            a = ad.AnnData(cnt[ii].astype(np.float32)); a.obs_names = T.obs_names[ii]
            try:
                sc.pp.scrublet(a, expected_doublet_rate=0.06, n_prin_comps=min(30, len(ii) - 1), random_state=0, verbose=False)
                rows.append(a.obs[['doublet_score']])
            except Exception as e:
                print('scrublet failed', smp, len(ii), e); rows.append(pd.DataFrame({'doublet_score': np.nan}, index=a.obs_names))
        pd.concat(rows).to_csv(sfile)
    scr = pd.read_csv(sfile, index_col=0)['doublet_score'].reindex(T.obs_names)
    fox = C['FOXP3'] > 0; il17 = (C['IL17A'] > 0) | (C['IL17F'] > 0)
    db_rows = []
    for g, cells in groups.items():
        s = scr.loc[cells].dropna()
        db_rows.append({'group': g, 'n_cells': len(cells), 'n_scrublet_scored': len(s),
                        'median_scrublet_score': s.median(), 'pct_scrublet_gt_0.2': 100 * (s > 0.2).mean(),
                        'pct_scrublet_gt_0.3': 100 * (s > 0.3).mean(),
                        'median_n_umi': T.obs.loc[cells, 'n_umi'].median(), 'median_n_genes': T.obs.loc[cells, 'n_genes'].median(),
                        'pct_FOXP3pos': 100 * fox.loc[cells].mean(), 'pct_IL17AorFpos': 100 * il17.loc[cells].mean(),
                        'n_FOXP3pos_and_IL17AorFpos': int((fox & il17).loc[cells].sum()),
                        'n_FOXP3pos_and_IL17AorFpos_and_CD40LGpos': int((fox & il17 & (C['CD40LG'] > 0)).loc[cells].sum())})
    db = pd.DataFrame(db_rows)
    for col, ser in [('scrublet', scr), ('n_umi', T.obs['n_umi'])]:   # Mann-Whitney: Th17-like vs other Tregs
        x1, x0 = ser.loc[th17_strict].dropna(), ser.loc[other_strict].dropna()
        db.loc[db.group == 'Th17-like Treg (Fig 4D/E Set, Top 10%)', f'mwu_p_vs_other_treg_{col}'] = mannwhitneyu(x1, x0).pvalue
    w3(db, 'th17_check_doublets_by_group.csv')

    # (4) Fig 4E suppressive-score trend after excluding Th17-like Tregs (sign test, Friedman)
    tp.X = tp.layers['log1p']
    sc.tl.score_genes(tp, ['FOXP3', 'CTLA4', 'IL2RA', 'BATF', 'CD27'], score_name='Suppressive')   # as in 12 (Fig 4E)
    def trend(o):
        o = o[['patient', 'pt_bin', 'Suppressive']].copy(); o['pt_bin'] = o['pt_bin'].astype(int)
        pb = o.groupby(['patient', 'pt_bin']).agg(score=('Suppressive', 'mean'), n=('Suppressive', 'size')).reset_index(); pb = pb[pb.n >= 3]
        rho = pd.Series({p: spearmanr(g.pt_bin, g.score)[0] for p, g in pb.groupby('patient') if g.pt_bin.nunique() >= 3}).dropna()
        neg = int((rho < 0).sum()); w = pb.pivot(index='patient', columns='pt_bin', values='score').dropna()
        return dict(n_tregs=len(o), n_patients=len(rho), n_declining=neg, sign_test_p=binomtest(neg, len(rho), 0.5).pvalue,
                    friedman_n=len(w), friedman_p=friedmanchisquare(*[w[c] for c in w.columns]).pvalue,
                    **{f'mean_score_bin{b}': pb[pb.pt_bin == b].score.mean() for b in range(5)})
    foxneg_th17 = th17_strict[C.loc[th17_strict, 'FOXP3'].values == 0]
    foxneg_il17_th17 = th17_strict[(C.loc[th17_strict, 'FOXP3'].values == 0) & il17.loc[th17_strict].values]
    o = tp.obs
    ptd = pd.DataFrame([dict(analysis='All pseudotime Tregs (Fig 4E)', **trend(o)),
                        dict(analysis='Excluding FOXP3- Th17-like Tregs', **trend(o.drop(foxneg_th17))),
                        dict(analysis='Excluding FOXP3- IL17A/F+ Th17-like Tregs', **trend(o.drop(foxneg_il17_th17))),
                        dict(analysis='Excluding all Th17-like Tregs', **trend(o.drop(th17_strict)))])
    ptd['n_th17like_in_set'] = [len(th17_strict), len(th17_strict) - len(foxneg_th17), len(th17_strict) - len(foxneg_il17_th17), 0]
    w3(ptd, 'th17_check_fig3E_sensitivity.csv')
    bins = pd.crosstab(np.where(o.index.isin(th17_strict), 'Th17-like', 'Other'), o['pt_bin'].astype(int))
    bt = bins.T.reset_index(); bt['pct_th17like'] = 100 * bt['Th17-like'] / (bt['Th17-like'] + bt['Other'])
    w3(bt, 'th17_check_th17like_by_pseudotime_bin.csv')

    # (5) per-gene Fisher tests (BH within comparison): Th17-like Treg vs other Treg / vs conventional Th17
    st_rows = []
    for ref_name, ref in [('Other Treg (Fig 4D/E Set)', other_strict), ('Conventional Th17 (FOXP3- IL17A/F+ CD4)', conv_th17)]:
        for gene in GENES:
            a1 = int((C.loc[th17_strict, gene] > 0).sum()); n1 = len(th17_strict)
            a0 = int((C.loc[ref, gene] > 0).sum()); n0 = len(ref)
            orr, p = fisher_exact([[a1, n1 - a1], [a0, n0 - a0]])
            pos1 = L.loc[th17_strict, gene][C.loc[th17_strict, gene] > 0]; pos0 = L.loc[ref, gene][C.loc[ref, gene] > 0]
            st_rows.append(dict(comparison=f'Th17-like Treg (Fig 4D/E Set) vs {ref_name}', gene=gene,
                                pct_th17like_treg=100 * a1 / n1, pct_reference=100 * a0 / n0, odds_ratio=orr, fisher_p=p,
                                mean_log1p_in_expressing_th17like_treg=pos1.mean() if len(pos1) else np.nan,
                                mean_log1p_in_expressing_reference=pos0.mean() if len(pos0) else np.nan))
    stt = pd.DataFrame(st_rows)
    stt['fdr_bh'] = stt.groupby('comparison')['fisher_p'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
    stt['significant_fdr_0.05'] = stt['fdr_bh'] < 0.05
    w3(stt, 'th17_check_marker_tests.csv')

    # (6) figure: Fig 4B UMAP with the groups (labels auto-placed on empty space) + marker dotplot
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 15, 'axes.labelsize': 14, 'font.family': 'DejaVu Sans'})
    FIG_GENES = ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TNFRSF18', 'RORC', 'CCR6', 'IL17A', 'IL17F', 'CD40LG', 'IL7R']
    FW, FH = 15.0, 5.6
    fig = plt.figure(figsize=(FW, FH))
    ax = fig.add_axes([0.0, 0.02, FH * 0.86 / FW, 0.86])
    P_ = XY.values
    ax.scatter(P_[:, 0], P_[:, 1], s=1, c='#e6e6e6', lw=0, rasterized=True)
    for l, col in [(TREG_LAB, '#f6c9a0'), (TRM_LAB, '#c9dcef')]:
        z = P_[lab.values == l]; ax.scatter(z[:, 0], z[:, 1], s=1, c=col, lw=0, rasterized=True)
    z = XY.loc[conv_th17].values; ax.scatter(z[:, 0], z[:, 1], s=5, c='#2b6cb0', lw=0, alpha=0.8, rasterized=True)
    z = XY.loc[th17_strict].values; ax.scatter(z[:, 0], z[:, 1], s=12, c='#c0392b', lw=0, rasterized=True)
    x0, x1 = P_[:, 0].min(), P_[:, 0].max(); y0, y1 = P_[:, 1].min(), P_[:, 1].max()
    px, py = 0.06 * (x1 - x0), 0.30 * (y1 - y0)
    ax.set_xlim(x0 - px, x1 + px); ax.set_ylim(y0 - 0.04 * (y1 - y0), y1 + py)    # head-room for the legend
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('UMAP1', labelpad=2); ax.set_ylabel('UMAP2', labelpad=2)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title('Fig 4B T-Cell UMAP', loc='center', pad=4)
    pct_treg_nb = pos.set_index('group')['knn30_harmony_pct_Treg']
    r1n, r2n = pct_treg_nb['Th17-like Treg (Fig 4D/E Set, Top 10%)'], pct_treg_nb['Conventional Th17 (FOXP3- IL17A/F+ CD4)']
    leg = ax.legend(handles=[Line2D([], [], marker='o', ls='', color='#c0392b', ms=8,
                                    label=f'Th17-Like Treg (n = {len(th17_strict)}; {r1n:.0f}% Treg Neighbours)' if r1n else ''),
                             Line2D([], [], marker='o', ls='', color='#2b6cb0', ms=8,
                                    label=f'Conventional Th17 (n = {len(conv_th17):,}; {r2n:.0f}% Treg Neighbours)' if r2n else '')],
                    loc='upper left', bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=12.5, handletextpad=0.1, borderaxespad=0.1)
    fig.canvas.draw()
    inv = ax.transData.inverted()
    def bbox_data(artist, pad=0.0):
        bb = artist.get_window_extent(fig.canvas.get_renderer()).expanded(1 + pad, 1 + pad)
        (a, b), (c, d) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]]); return a, b, c, d
    def n_inside(bx):
        a, b, c, d = bx; return int(((P_[:, 0] > a) & (P_[:, 0] < c) & (P_[:, 1] > b) & (P_[:, 1] < d)).sum())
    def ovl(b1, b2): return not (b1[2] < b2[0] or b2[2] < b1[0] or b1[3] < b2[1] or b2[3] < b1[1])
    taken = [bbox_data(leg)]
    gx = np.linspace(x0 - px, x1 + px, 60); gy = np.linspace(y0, y1 + py, 60)
    for l, txt in [(TREG_LAB, 'Treg'), ('CD4+ Effector', 'CD4+ Effector'), (TRM_LAB, 'CD4+ Trm')]:
        cx, cy = cen.loc[l]
        cand = sorted(((x, y) for x in gx for y in gy), key=lambda q: (q[0] - cx) ** 2 + (q[1] - cy) ** 2)
        t = ax.text(0, 0, txt, ha='center', va='center', fontsize=14, fontweight='bold')
        for (x, y) in cand:
            t.set_position((x, y)); bx = bbox_data(t, pad=0.12)
            if bx[0] < x0 - px or bx[2] > x1 + px or bx[1] < y0 or bx[3] > y1 + py: continue
            if n_inside(bx) == 0 and not any(ovl(bx, o) for o in taken):
                break
        taken.append(bx)
        ax.annotate('', xy=(cx, cy), xytext=(x, y), arrowprops=dict(arrowstyle='-', color='black', lw=1.0, shrinkA=9, shrinkB=0))
    L2 = 0.495; ax2 = fig.add_axes([L2, 0.22, 0.86 - L2, 0.60])
    rows = [('Other Tregs', other_strict), ('Th17-Like Tregs', th17_strict), ('Conventional Th17', conv_th17),
            ('CD4+ Trm', groups['Cluster: CD4+ Trm'])]
    pct = np.array([[100 * (C.loc[c, g] > 0).mean() for g in FIG_GENES] for _, c in rows])
    mu = np.array([[L.loc[c, g].mean() for g in FIG_GENES] for _, c in rows])
    mus = mu / mu.max(0, keepdims=True)
    test_rows = [('Other Treg (Fig 4D/E Set)', 'Th17-Like vs Other Tregs'), ('Conventional Th17 (FOXP3- IL17A/F+ CD4)', 'Th17-Like vs Conv. Th17')]
    TY = [len(rows) - 0.25, len(rows) + 0.35]
    ax2.set_xlim(-0.5, len(FIG_GENES) - 0.5); ax2.set_ylim(TY[1] + 0.35, -0.5)
    fig.canvas.draw()   # dot size: 100% dot diameter = 76% of the smaller grid-cell side (points)
    bb = ax2.get_window_extent(); cw = bb.width / len(FIG_GENES) * 72 / fig.dpi; ch = bb.height / (TY[1] + 0.85) * 72 / fig.dpi
    SZ = (0.76 * min(cw, ch)) ** 2
    xx, yy = np.meshgrid(np.arange(len(FIG_GENES)), np.arange(len(rows)))
    sca = ax2.scatter(xx.ravel(), yy.ravel(), s=(pct.ravel() / 100) * SZ + 2, c=mus.ravel(), cmap='Reds', vmin=0, vmax=1,
                      edgecolor='#555555', lw=0.4)
    for j, (cmp, _) in enumerate(test_rows):
        sub = stt[stt.comparison.str.endswith(cmp)].set_index('gene')
        for i, g in enumerate(FIG_GENES):
            if sub.loc[g, 'significant_fdr_0.05']:
                up = sub.loc[g, 'pct_th17like_treg'] > sub.loc[g, 'pct_reference']
                ax2.text(i, TY[j], '▲' if up else '▼', ha='center', va='center', fontsize=10, color='#c0392b' if up else '#2b6cb0')
    ax2.set_yticks(list(range(len(rows))) + TY); ax2.set_yticklabels([r[0] for r in rows] + [t[1] for t in test_rows])
    for tl in ax2.get_yticklabels()[len(rows):]: tl.set_fontsize(12.5); tl.set_color('#444444')
    ax2.set_xticks(range(len(FIG_GENES))); ax2.set_xticklabels(FIG_GENES, rotation=45, ha='right', rotation_mode='anchor')
    ax2.tick_params(length=0, pad=3)
    for xv in (4.5, 8.5): ax2.axvline(xv, color='#aaaaaa', lw=0.8, ls=':')
    ax2.axhline(len(rows) - 0.6, color='#cccccc', lw=0.8)
    for xc, txt in [(2, 'Treg'), (6.5, 'Th17'), (9.5, 'Effector')]:
        ax2.text(xc, -0.55, txt, ha='center', va='bottom', fontsize=14, fontweight='bold')
    for s in ax2.spines.values(): s.set_visible(False)
    ax2.set_title('Marker Expression (▲▼ Fisher FDR < 0.05)', loc='center', pad=26)
    cax = fig.add_axes([0.875, 0.56, 0.012, 0.26])
    cb = fig.colorbar(sca, cax=cax); cb.set_label('Scaled Mean', fontsize=12); cb.ax.tick_params(labelsize=11)
    cb.set_ticks([0, 0.5, 1])
    hs = [Line2D([], [], marker='o', ls='', mfc='white', mec='#555555', ms=np.sqrt((p_ / 100) * SZ + 2)) for p_ in (25, 50, 100)]
    fig.legend(hs, ['25%', '50%', '100%'], title='% Expressing', loc='upper left', bbox_to_anchor=(0.866, 0.50), frameon=False,
               fontsize=12, title_fontsize=12, labelspacing=1.25, borderaxespad=0, handletextpad=0.9, handlelength=1.6)
    fig.suptitle('Th17-Like Tregs Are Tregs, Not Conventional Th17', fontsize=16, fontweight='bold', x=0.5, y=0.985)
    fig.savefig(f'{OUT}/th17_check_umap_dotplot.png', dpi=300, bbox_inches='tight', pad_inches=0.06); plt.close(fig)


# ---- Section S8CD_leiden_root: Leiden sub-clusters and pseudotime bins on the Treg UMAP (Supp 8C); root selection =
# leiden_0.6 sub-cluster with >= 20 cells and the highest median Treg_core (Supp 8D)
if _run('S8CD_leiden_root'):
    import os, warnings
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, scanpy as sc
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.patches import Patch

    os.makedirs(OUT, exist_ok=True)
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 15, 'savefig.dpi': 300})

    def save(fig, name):
        fig.savefig(f'{OUT}/{name}', dpi=300, bbox_inches='tight'); plt.close(fig)

    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    o = t.obs
    o['pt_bin'] = pd.Categorical(o['pt_bin'].astype(int).astype(str), categories=['0', '1', '2', '3', '4'])
    o['leiden_0.6'] = o['leiden_0.6'].astype(str)
    o['leiden_0.6'] = pd.Categorical(o['leiden_0.6'], categories=sorted(o['leiden_0.6'].unique(), key=int))
    ROOT = str(t.uns['root_subcluster']); IROOT = int(t.uns['iroot'])
    U = t.obsm['X_umap']
    base = pd.DataFrame({'cell': o.index, 'umap1': U[:, 0], 'umap2': U[:, 1]})

    # Leiden + pseudotime-bin UMAPs
    fig, ax = plt.subplots(figsize=(4.8, 4.4))
    sc.pl.umap(t, color='leiden_0.6', s=18, legend_loc='on data', legend_fontsize=11, legend_fontoutline=2,
               frameon=False, ax=ax, show=False)
    ax.set_title('Leiden Sub-Clusters', loc='center')
    save(fig, 'ext8A_leiden_umap.png')
    w3(base.assign(leiden_0_6=o['leiden_0.6'].values), 'ext8A_leiden_umap.csv')

    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    sc.pl.umap(t, color='pt_bin', s=18, alpha=0.8, frameon=False, ax=ax, show=False)
    ax.set_title('Pseudotime Bins', loc='center')
    ax.legend(title='Pseudotime Bin', loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=10,
              title_fontsize=10, markerscale=1.2)
    save(fig, 'ext8A_pseudotime_bins_umap.png')
    w3(base.assign(pt_bin=o['pt_bin'].astype(str).values, dpt_pseudotime=o['dpt_pseudotime'].values),
       'ext8A_pseudotime_bins_umap.csv')
    ct = pd.crosstab(o['leiden_0.6'], o['pt_bin'])
    ct.columns = [f'n_bin{c}' for c in ct.columns]
    ct['n_cells'] = ct.sum(1); ct['frac_in_bins0_1'] = (ct['n_bin0'] + ct['n_bin1']) / ct['n_cells']
    w3(ct.reset_index().rename(columns={'leiden_0.6': 'leiden_0_6'}), 'ext8A_leiden_x_ptbin_counts.csv')

    # root selection table
    cs = o.groupby('leiden_0.6', observed=True).agg(n_cells=('Treg_core', 'size'), Treg_core_median=('Treg_core', 'median'),
                                                    Treg_core_mean=('Treg_core', 'mean'), Cell_cycle_median=('Cell_cycle', 'median'),
                                                    dpt_median=('dpt_pseudotime', 'median'))
    cs['eligible_n_ge_20'] = cs['n_cells'] >= 20
    cs['chosen_root'] = cs.index == ROOT
    cs['frac_in_bins0_1'] = ct['frac_in_bins0_1']
    cs = cs.sort_values('Treg_core_median', ascending=False)

    def centroid_cell(cl):   # cell closest to the UMAP centroid of the sub-cluster
        cells = t.obs_names[o['leiden_0.6'].astype(str) == cl]
        idx = t.obs_names.get_indexer(cells); c = U[idx].mean(axis=0)
        return cells[np.argmin(((U[idx] - c) ** 2).sum(axis=1))]

    ALTS = [c for c in o['leiden_0.6'].astype(str).value_counts().index if c != ROOT][:2]   # two largest non-root sub-clusters
    roots = [('DPT root cell', o.index[IROOT], ROOT)] + [(f'Clus {ROOT} centroid', centroid_cell(ROOT), ROOT)] + \
            [(f'Clus {c} centroid', centroid_cell(c), c) for c in ALTS]
    cs['root_cells'] = [';'.join(f'{lab}={cell}' for lab, cell, cl in roots if cl == i) for i in cs.index]
    w3(cs.reset_index().rename(columns={'leiden_0.6': 'leiden_0_6'}), 'ext8C_root_selection_treg_core.csv')

    # Treg core per sub-cluster ordered by median; root highlighted, ineligible (< 20 cells) grey
    order = list(cs.index)
    pal = {c: ('#a16f8d' if c == ROOT else ('#d9d9d9' if not cs.loc[c, 'eligible_n_ge_20'] else '#9ecae1')) for c in order}
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    sns.boxplot(data=o, x='leiden_0.6', y='Treg_core', order=order, hue='leiden_0.6', hue_order=order, palette=pal,
                legend=False, fliersize=0, width=0.65, ax=ax)
    sns.stripplot(data=o, x='leiden_0.6', y='Treg_core', order=order, color='k', size=1.2, alpha=0.25, jitter=0.25, ax=ax)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{c}\n(n={cs.loc[c, 'n_cells']})" for c in order], fontsize=9.5)
    ax.set_ylim(o['Treg_core'].quantile(0.002), o['Treg_core'].quantile(0.998))
    ax.set_xlabel('Leiden Sub-Cluster (Ordered by Median)', fontsize=11)
    ax.set_ylabel('Treg Core Score', fontsize=11)
    ax.set_title('Root Selection', loc='center')
    ax.legend(handles=[Patch(color='#a16f8d', label=f'Root (Clus {ROOT})'), Patch(color='#9ecae1', label='Eligible (>=20 Cells)'),
                       Patch(color='#d9d9d9', label='Excluded (<20 Cells)')], loc='upper center', bbox_to_anchor=(0.5, -0.17),
              ncol=3, frameon=False, fontsize=8, handlelength=1)
    sns.despine(ax=ax)
    save(fig, 'ext8C_root_selection_treg_core.png')


# ---- Section S8E_palantir_main_root: Palantir from the main DPT root cell on the Harmony PCA (5 and 10 diffusion
# components); Spearman rho vs the main DPT and early/late agreement (value quoted in the Methods)
if _run('S8E_palantir_main_root'):
    import warnings
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, palantir
    from scipy import stats
    treg = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad'); o = treg.obs; o['leiden_0.6'] = o['leiden_0.6'].astype(str)
    ROOT = str(treg.uns['root_subcluster']); rc = o.index[int(treg.uns['iroot'])]; coords = treg.obsm['X_umap']
    ref = o['dpt_pseudotime'].values; ref_early = pd.qcut(pd.Series(ref).rank(method='first'), 5, labels=False).values <= 1
    isroot = (o['leiden_0.6'] == ROOT).values; rows = []; PT = {}
    for nc in (5, 10):
        p = treg.copy(); palantir.utils.run_diffusion_maps(p, pca_key='X_pca_harmony', n_components=nc); palantir.utils.determine_multiscale_space(p)
        palantir.core.run_palantir(p, rc, num_waypoints=500, use_early_cell_as_start=True)
        v = p.obs['palantir_pseudotime'].values.copy(); PT[f'PalantirHarmony_exactroot_nc{nc}'] = v
        early = pd.qcut(pd.Series(v).rank(method='first'), 5, labels=False).values <= 1; rho, pv = stats.spearmanr(v, ref)
        rows.append(dict(pseudotime=f'PalantirHarmony_exactroot_nc{nc}', root_cell=rc, root_subcluster=ROOT,
                         start_cell_is_root=bool(o.index[int(np.argmin(v))] == rc),
                         spearman_rho_vs_main_DPT=rho, spearman_p=pv, early_late_agreement=float(np.mean(early == ref_early)),
                         root_subcluster_frac_in_bins0_1=float(early[isroot].mean())))
    pd.DataFrame(rows).to_csv(f'{OUT}/ext8E_palantir_main_root_concordance.csv', index=False, float_format='%.3g')
    pd.DataFrame({'cell': treg.obs_names, 'leiden_0_6': o['leiden_0.6'].values, 'umap1': coords[:, 0], 'umap2': coords[:, 1],
                  'DPT_main': ref, **PT}).to_csv(f'{OUT}/ext8E_palantir_main_root_per_cell.csv', index=False, float_format='%.3g')


# ---- Section S8E_root_comparison: DPT (stored Harmony k=20 graph) and Palantir (10 components) from ten alternative
# root cells in the root sub-cluster: UMAP centroid cell, highest Treg_core cell, 8 random cells (seed 0)
if _run('S8E_root_comparison'):
    import warnings
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, scanpy as sc, palantir
    from scipy import stats
    treg = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad'); o = treg.obs; o['leiden_0.6'] = o['leiden_0.6'].astype(str)
    ROOT = str(treg.uns['root_subcluster']); main_root = o.index[int(treg.uns['iroot'])]; ref = o['dpt_pseudotime'].values
    cells = o.index[o['leiden_0.6'] == ROOT]; idx = treg.obs_names.get_indexer(cells); um = treg.obsm['X_umap'][idx]
    cent = cells[np.argmin(((um - um.mean(0)) ** 2).sum(1))]; topcore = o.loc[cells, 'Treg_core'].idxmax()
    rng = np.random.default_rng(0); rand = list(rng.choice([c for c in cells if c not in (main_root, cent, topcore)], 8, replace=False))
    ROOTS = [('main root (first cell)', main_root), ('UMAP centroid', cent), ('highest Treg-core', topcore)] + [(f'random {k + 1}', c) for k, c in enumerate(rand)]
    p = treg.copy(); palantir.utils.run_diffusion_maps(p, pca_key='X_pca_harmony', n_components=10); palantir.utils.determine_multiscale_space(p)
    rows = []; PERCELL = {'main_DPT': ref}
    for lab, rc in ROOTS:
        t = treg.copy(); t.uns['iroot'] = int(treg.obs_names.get_loc(rc)); sc.tl.dpt(t); v_dpt = t.obs['dpt_pseudotime'].values
        palantir.core.run_palantir(p, rc, num_waypoints=500, use_early_cell_as_start=True); v_pal = p.obs['palantir_pseudotime'].values
        PERCELL[f'DPT|{lab}'] = v_dpt; PERCELL[f'Palantir|{lab}'] = v_pal
        rows.append(dict(root=lab, cell=rc, Treg_core=round(float(o.loc[rc, 'Treg_core']), 3),
                         rho_DPT_vs_main=stats.spearmanr(v_dpt, ref)[0], rho_Palantir_vs_main=stats.spearmanr(v_pal, ref)[0]))
    pd.DataFrame(PERCELL, index=treg.obs_names).to_csv(f'{OUT}/ext8DE_root_robustness_within_clus{ROOT}_per_cell.csv.gz', float_format='%.4g')
    pd.DataFrame(rows).to_csv(f'{OUT}/ext8DE_root_robustness_within_clus{ROOT}.csv', index=False, float_format='%.3g')


# ---- Section S8E_root_comparison_plot: root cells tested on the Treg UMAP; Spearman rho vs the main DPT (Supp 8E)
if _run('S8E_root_comparison_plot'):
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    treg = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad', backed='r'); o = treg.obs; um = np.asarray(treg.obsm['X_umap'])
    ROOT = str(treg.uns['root_subcluster']); R = pd.read_csv(f'{OUT}/ext8DE_root_robustness_within_clus{ROOT}.csv')
    inr = (o['leiden_0.6'].astype(str) == ROOT).values
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={'width_ratios': [1.1, 1]})
    a1.scatter(um[~inr, 0], um[~inr, 1], s=4, c='#d9d9d9', lw=0, rasterized=True)
    a1.scatter(um[inr, 0], um[inr, 1], s=6, c='#6baed6', lw=0, rasterized=True, label=f'Root sub-cluster (clus {ROOT})')
    for _, r in R.iterrows():
        i = treg.obs_names.get_loc(r['cell']); main = r['root'].startswith('main')
        a1.scatter(um[i, 0], um[i, 1], marker='*' if main else 'o', s=260 if main else 55, c='k' if main else '#e6550d', edgecolor='w', lw=0.8, zorder=5)
    a1.scatter([], [], marker='*', s=160, c='k', label='Main root cell'); a1.scatter([], [], marker='o', s=45, c='#e6550d', label='Alternative root cells (n = 10)')
    a1.legend(loc='lower left', fontsize=9.5, frameon=False); a1.set_axis_off(); a1.set_title('Root Cells Tested', fontsize=15, loc='center')
    alt = R[~R['root'].str.startswith('main')]
    for k, (col, lab, colr) in enumerate([('rho_DPT_vs_main', 'Diffusion\nPseudotime', '#3182bd'), ('rho_Palantir_vs_main', 'Palantir\n(10 components)', '#e6550d')]):
        v = alt[col].values; x = k + np.random.default_rng(1).uniform(-0.12, 0.12, len(v))
        a2.scatter(x, v, s=55, c=colr, edgecolor='w', lw=0.6, zorder=3)
        m = np.median(v); a2.plot([k - 0.25, k + 0.25], [m, m], c='k', lw=2.2, zorder=4)
        a2.text(k + 0.29, m, f'median {m:.2f}', va='center', fontsize=10.5)
    a2.axhline(1, c='grey', lw=0.8, ls='--'); a2.set_xticks([0, 1]); a2.set_xticklabels(['Diffusion\nPseudotime', 'Palantir\n(10 components)'], fontsize=12)
    a2.set_xlim(-0.6, 1.75); a2.set_ylim(0.4, 1.03); a2.set_ylabel('Spearman ρ vs Main Pseudotime', fontsize=12)
    a2.set_title('Alternative Roots Within Clus 4', fontsize=15, loc='center'); a2.spines[['top', 'right']].set_visible(False)
    plt.tight_layout(); fig.savefig(f'{OUT}/ext8DE_root_robustness_within_clus{ROOT}.png', dpi=300, bbox_inches='tight')
