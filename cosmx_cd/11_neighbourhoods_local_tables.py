#!/usr/bin/env python
"""CD CosMx 6K resections: cellular neighbourhoods and per-Treg local-environment tables.

Panels made
  none directly; writes the shared tables used by 12-17 (Fig 5A-F, Supp 9A-F)

Inputs
  /share/fsmresfiles/UC/AtoMx/CD_surgical_resection/combined/ (all_samples.parquet, resection cell metadata with global
    coordinates and AUCell labels; spot_lr_df.parquet from 07_resection_liana_merge.ipynb; raw_counts_strict_tregs.csv from
    05_resection_treg_counts.ipynb); per-sample Processed/norm_log1p.h5ad (03_resection_processing.ipynb)
  /share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged/h5ad_obj/norm_anno_tregbin.h5ad (WTX biopsies,
    02_wtx_biopsy_treg_labels_liana.ipynb; discovery neighbourhoods)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/D_shared_inputs/ (10_treg_suppression_score.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/B_resection_analyses/
    cellular_neighborhoods_v2/ (neighbourhoods per cell, composition and Treg-state tables, biopsy discovery / resection validation)
    ligand_receptor_state/ (per-Treg local ligand, receptor, cell-type and CXCL12 tables; CXCL12 models)

Run order of sections
  cellular_neighbourhoods -> treg_local_ligand_receptor -> treg_local_composition -> cxcl12_expected_producer_field

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist. Command-line arguments and environment variables used by a section are listed in its
  own docstring.
  Paths are the lab's absolute paths; edit them to match the local data location.
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
# Section cellular_neighbourhoods: cellular neighbourhoods (20 nearest cells, k-means k = 10) and Treg states per neighbourhood
# Panels: Fig 5B-C, 5F, Supp 9A-C (cellular neighbourhoods, k = 10)
# ======================================================================================================================
if _run('cellular_neighbourhoods'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/cellular_neighbourhoods.py'   # working folder of this step (relative paths resolve from here)
    """Cellular neighbourhoods (CN) in the CD resections (CosMx 6K + custom) and the spatial niches of activated (early-like)
vs low-suppressive (late-like) Tregs.
Method (Schurch et al., Cell 2020):
  cells      all cells of the 11 resection sections (combined/all_samples.parquet; global coordinates CenterX/Y_global_px,
             um = px x 0.12028); cell types = the cell-type grouping of the LR modules (network_anno) with all Tregs as one 'Treg'
             type, neutrophils split from other granulocytes, Tfh split from exhausted CD4 Tfh, CD8 Trm split by GZMK, and SELENOP+
             Fib cells removed (not part of any window or CN), so the CNs do not depend on Treg state.
  windows    composition of the 20 nearest cells (excluding the cell itself), neighbours within the same section only.
  CNs        MiniBatchKMeans on the composition fractions, k = 10 (primary), k = 8 and 12 (sensitivity; adjusted Rand
             index vs k = 10), random_state 0. CN names = cell types with the highest log2 enrichment (CN vs all cells).
  Treg state D6_tertile (scRNA-seq bin-4-vs-bin-0 signature, top vs bottom tertile within section;
             treg_definition_comparison_per_cell.csv). Early = activated (early-like), Late = low-suppressive (late-like).
  test       per CN: 2 x 2 per section (state x in CN / not in CN), Cochran-Mantel-Haenszel odds ratio across sections
             (OR > 1 = activated Tregs enriched in the CN relative to low-suppressive), BH across CNs; per-section sign (activated
             fraction in CN > low-suppressive fraction in CN, sections with >= 5 Tregs of each state). Treg density per CN.
  modules    per-cell LIANA bivariate LR scores (combined/spot_lr_df.parquet; module columns only; NaN = 0, no interaction),
             per-cell module score = mean over the module's LR pairs (4 modules; Ag Presentation/Treg = the 13-pair module of
             10_treg_suppression_score.py, fig4_ag_treg_module_lrs.csv), averaged per CN.
Discovery / validation:
  discovery  CD biopsies, CosMx WTX (norm_anno_tregbin.h5ad; 4 patients I0368, I0476 inflamed CD, I0453 non-inflamed CD,
             I0477 healthy control; 10 tissue pieces = obs 'slide', one patient each; obsm['spatial_fov'] global px, neighbours
             within piece). CNs are built on the biopsies (same window / k-means settings), Treg-state enrichment per CN by CMH
             stratified by tissue piece; descriptive by condition.
  shared vocabulary  both cohorts carry the same three-level annotation (cell_category / cell_type_general / cell_type_short),
             so the same network_anno rule maps both to one vocabulary (mapping table written to CSV; unmapped types reported).
  validation resection cells: 20-NN windows on the shared vocabulary assigned to the nearest biopsy CN centroid (k-means
             predict); Treg-state enrichment per CN by CMH stratified by section; per-CN concordance (direction, q) and
             correlation of CN composition profiles between cohorts.
  secondary  the resection de-novo CNs (above) are kept as a sensitivity result; ARI between resection de-novo CNs and
             projected biopsy CNs.
Outputs: panels/fig4_and_related/B_resection_analyses/cellular_neighborhoods_v2/ (resection de novo; cn_per_cell.parquet,
cn_k10_celltype_log2_enrichment.csv, cn_treg_state_enrichment_all.csv, ...), .../discovery_biopsy/, .../validation_resection/
Run: OPENBLAS_CORETYPE=SkylakeX /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, gc, warnings
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX', 'run with OPENBLAS_CORETYPE=SkylakeX'
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    from scipy.spatial import cKDTree
    from scipy import sparse
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.metrics import adjusted_rand_score
    from statsmodels.stats.contingency_tables import StratifiedTable
    from statsmodels.stats.multitest import multipletests
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt, seaborn as sns

    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
    P = f'{D}/panels/fig4_and_related'; OUT = f'{P}/B_resection_analyses/cellular_neighborhoods_v2'; os.makedirs(OUT, exist_ok=True)
    UM, KNN, KS, K0 = 0.12028, 20, [8, 10, 12], 10
    LOG = []
    def log(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.append(s)


    def stars(q):
        return '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''


    def network_anno(df):                                              # cell-type grouping of the LR modules; Tregs = one type; neutrophils split
        out = []
        for i, j, k in zip(df['cell_category'].astype(str), df['cell_type_general'].astype(str), df['cell_type_short'].astype(str)):
            if i == 'Myeloid':
                if 'DC' in j: out.append(j)
                elif k in ['Mac S+XS-', 'Mac S+SG+']: out.append('Mac S')
                elif k == 'Mac S+M+P+': out.append('Mac S_M_P')
                elif k == 'Mac S+M+S+': out.append('Mac S_M_S')
                elif 'CD9' in j: out.append('CD9 Mac')
                else: out.append(k)
            elif i == 'Connective':
                if 'OGN' in j: out.append('OGN_RSPO3 Fib')
                elif 'ADAMDEC1' in j: out.append('ADAMDEC1 Fib')
                elif 'PLA2G2A' in j: out.append('PLA2G2A Fib')
                elif 'SELENOP' in j: out.append('DROP')   # SELENOP+ Fib not part of any window or CN
                else: out.append(j)
            elif i == 'T':
                if j == 'Treg': out.append('Treg')
                elif j == 'Exhausted CD4 Tfh': out.append('Exhausted CD4 Tfh')   # split from Tfh
                elif 'Tfh' in j: out.append('Tfh')
                elif j == 'CD8 Trm': out.append('CD8 GZMK+ Trm' if 'GZMK+' in k else 'CD8 GZMK- Trm')   # split by GZMK
                else: out.append(j)
            elif i == 'Granulocyte': out.append('Neutrophil' if j == 'Neutrophil' else 'Granulocyte')
            elif i.lower() == 'cycling': out.append('Cycling')
            else: out.append(i)
        return out


    # ------------------------------------------------------------------------------------------------ cells + windows
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['fov', 'CenterX_global_px', 'CenterY_global_px', 'cell_category',
                                                                        'cell_type_general', 'cell_type_short', 'patient_sec',
                                                                        '__index_level_0__']).to_pandas()
    sa.index = [p + '_' + str(i) for p, i in zip(sa['patient_sec'], sa.index)]
    sa['section'] = sa['patient_sec'].astype(str); sa['fov'] = sa['fov'].astype(str)
    sa['x'], sa['y'] = sa.CenterX_global_px.values * UM, sa.CenterY_global_px.values * UM
    sa['ct'] = network_anno(sa); sa = sa[sa.ct != 'DROP'].copy()
    CT = sorted(sa['ct'].unique()); ct_code = pd.Categorical(sa['ct'], categories=CT).codes
    log('cells', len(sa), '| sections', sa.section.nunique(), '| cell types', len(CT))
    rows, cols = [], []
    pos = np.arange(len(sa))
    for s, idx in sa.groupby('section').indices.items():
        xy = sa[['x', 'y']].values[idx]
        _, nn = cKDTree(xy).query(xy, k=KNN + 1)
        nn = idx[nn[:, 1:]]                                               # drop self (first neighbour)
        rows.append(np.repeat(idx, KNN)); cols.append(ct_code[nn.ravel()])
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    W = sparse.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(sa), len(CT))).toarray() / KNN
    del rows, cols; gc.collect()
    log('window matrix', W.shape, '| row sums ~1:', float(W.sum(1).min()), float(W.sum(1).max()))

    # ------------------------------------------------------------------------------------------------ CNs
    LAB = {}
    for k in KS:
        km = MiniBatchKMeans(n_clusters=k, random_state=0, batch_size=20000, n_init=5, max_iter=300).fit(W)
        LAB[k] = km.labels_
    ari = {k: adjusted_rand_score(LAB[K0], LAB[k]) for k in KS if k != K0}
    log('adjusted Rand index vs k=10:', {k: round(v, 3) for k, v in ari.items()})
    overall = pd.Series(np.bincount(ct_code, minlength=len(CT)) / len(sa), index=CT)


    def cn_tables(lab, k):
        comp = pd.DataFrame(W).groupby(lab).mean(); comp.columns = CT                   # mean window composition per CN
        cell = pd.crosstab(lab, sa['ct'].values, normalize='index').reindex(columns=CT, fill_value=0)   # cells of each type in the CN
        enr = np.log2((cell + 1e-4) / (overall + 1e-4))
        names = {}
        for c in range(k):
            top = enr.loc[c].sort_values(ascending=False)
            top = [t for t in top.index if top[t] > 0.5][:3] or [top.index[0]]
            names[c] = f'CN{c + 1}: ' + ' / '.join(top)
        return comp, cell, enr, names


    TAB = {k: cn_tables(LAB[k], k) for k in KS}
    comp, cell, enr, NAMES = TAB[K0]
    sa['CN'] = LAB[K0]
    cn_size = sa['CN'].value_counts().sort_index()
    log('CNs (k = 10):'); [log(f'  {NAMES[c]}  (n = {cn_size[c]:,}, {100 * cn_size[c] / len(sa):.1f}% of cells)') for c in range(K0)]
    enr.rename(index=NAMES).to_csv(f'{OUT}/cn_k10_celltype_log2_enrichment.csv', float_format='%.3g')
    cell.rename(index=NAMES).to_csv(f'{OUT}/cn_k10_celltype_fraction.csv', float_format='%.4g')
    for k in [8, 12]:
        TAB[k][2].rename(index=TAB[k][3]).to_csv(f'{OUT}/cn_k{k}_celltype_log2_enrichment.csv', float_format='%.3g')

    # ------------------------------------------------------------------------------------------------ Treg states
    PC = pd.read_csv(f'{P}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0)
    PC = PC[PC.cohort == 'resection']
    assert PC.index.isin(sa.index).all(), 'Treg not in all_samples'
    treg_ct = sa.loc[PC.index, 'ct'].value_counts()
    log('labelled Tregs:', len(PC), '| their network_anno type:', treg_ct.to_dict())
    DEFS = {'D6_tertile': 'signature (bin 4 vs bin 0), top vs bottom tertile'}


    def treg_cn_test(dname, k):
        t = pd.DataFrame({'section': sa.loc[PC.index, 'section'].values, 'state': PC[dname].values, 'CN': LAB[k][sa.index.get_indexer(PC.index)]})
        t = t[t.state.isin(['Early', 'Late'])]
        out = []
        for c in range(k):
            tabs, sgn = [], []
            for s, g in t.groupby('section'):
                e, l = g[g.state == 'Early'], g[g.state == 'Late']
                a, b = int((e.CN == c).sum()), int((e.CN != c).sum()); cc, d = int((l.CN == c).sum()), int((l.CN != c).sum())
                if a + cc == 0 or b + d == 0: continue                      # no Tregs in / out of the CN in this section
                tabs.append(np.array([[a, b], [cc, d]], float))
                if len(e) >= 5 and len(l) >= 5: sgn.append(a / len(e) > cc / len(l))
            st = StratifiedTable([x + 0.5 if (x == 0).any() else x for x in tabs]) if tabs else None
            n_e, n_l = int(((t.state == 'Early') & (t.CN == c)).sum()), int(((t.state == 'Late') & (t.CN == c)).sum())
            out.append(dict(definition=dname, k=k, CN=c, n_activated_in_CN=n_e, n_low_suppressive_in_CN=n_l,
                            frac_activated_in_CN=n_e / (t.state == 'Early').sum(), frac_low_suppressive_in_CN=n_l / (t.state == 'Late').sum(),
                            cmh_OR_activated_vs_low_suppressive=st.oddsratio_pooled if st else np.nan,
                            cmh_p=st.test_null_odds(correction=True).pvalue if st else np.nan,
                            sections_activated_higher=int(np.sum(sgn)), sections_tested=len(sgn)))
        o = pd.DataFrame(out); ok = o.cmh_p.notna(); o['cmh_q'] = np.nan
        o.loc[ok, 'cmh_q'] = multipletests(o.loc[ok, 'cmh_p'], method='fdr_bh')[1]
        return o


    RES = pd.concat([treg_cn_test(d, k) for d in DEFS for k in KS], ignore_index=True)
    RES['CN_name'] = [TAB[k][3][c] for k, c in zip(RES.k, RES.CN)]
    # Treg density per CN (all Tregs, state-independent)
    dens = sa.groupby('CN')['ct'].apply(lambda x: (x == 'Treg').mean())
    RES = RES.merge(pd.DataFrame({'k': K0, 'CN': dens.index, 'treg_fraction_of_cells': dens.values}), on=['k', 'CN'], how='left')
    RES.to_csv(f'{OUT}/cn_treg_state_enrichment_all.csv', index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    log(RES[RES.k == K0].drop(columns=['k']).round(3).to_string(index=False))

    # key CNs for robustness: per k, the CN most enriched for OGN_RSPO3 Fib, DCs (all DC types), B, Neutrophil, Epithelial
    DCS = [c for c in CT if 'DC' in c]
    def key_cn(k, kind):
        e = TAB[k][2]
        v = e[DCS].max(1) if kind == 'DC' else e[kind]
        return int(v.idxmax())
    rob = []
    for kind in ['OGN_RSPO3 Fib', 'DC', 'B', 'Neutrophil', 'Epithelial', 'Treg']:
        for d in DEFS:
            for k in KS:
                c = key_cn(k, kind); r = RES[(RES.definition == d) & (RES.k == k) & (RES.CN == c)].iloc[0]
                rob.append(dict(key_celltype=kind, definition=d, k=k, CN=r.CN_name, OR=r.cmh_OR_activated_vs_low_suppressive, p=r.cmh_p, q=r.cmh_q,
                                sections=f'{r.sections_activated_higher}/{r.sections_tested}'))
    ROB = pd.DataFrame(rob); ROB.to_csv(f'{OUT}/cn_robustness_key_niches.csv', index=False, float_format='%.4g')
    log(ROB.round(3).to_string(index=False))

    # ------------------------------------------------------------------------------------------------ module scores per CN
    ag = pd.read_csv(f'{P}/D_shared_inputs/fig4_ag_treg_module_lrs.csv')['lr_pair'].tolist()
    MODS = {'Ag Presentation/Treg': ag,
            'Immunomodulatory': ['LGALS1^PTPRC', 'IL12A^CD28', 'CD58^CD2', 'PTPRC^CD22', 'CD22^PTPRC'],
            'Pro-Fibrotic/ECM': ['TIMP2^MMP2', 'COL3A1^DDR2', 'COL17A1^ITGB1', 'COL9A1^ITGA2', 'COL27A1^ITGB1', 'TGM2^SDC4'],
            'Pro-Inflammatory': ['TNFSF13^TNFRSF1A', 'TNFSF13^TNFRSF14', 'IL6^IL6ST', 'IL1A^IL1RAP', 'CCL2^CCR3', 'CCL19^CCR7',
                                 'CCL25^CCR9', 'CCL25^ACKR4', 'CCL11^ACKR4']}
    names = pq.read_schema(f'{R}/combined/spot_lr_df.parquet').names
    cols_ = [c for c in dict.fromkeys(sum(MODS.values(), [])) if c in names]
    log('module LR columns present:', len(cols_), '| absent:', [c for c in dict.fromkeys(sum(MODS.values(), [])) if c not in names])
    lr = pq.read_table(f'{R}/combined/spot_lr_df.parquet', columns=cols_ + ['__index_level_0__']).to_pandas()
    if '__index_level_0__' in lr.columns: lr.index = lr.pop('__index_level_0__')
    lr = lr.astype(np.float32).fillna(0.0)
    inlr = sa.index.isin(lr.index); log('cells with LIANA scores:', int(inlr.sum()), 'of', len(sa))
    MS = pd.DataFrame({m: lr[[c for c in v if c in lr.columns]].mean(1) for m, v in MODS.items()})
    MS['CN'] = sa.loc[MS.index, 'CN'].values if MS.index.isin(sa.index).all() else sa['CN'].reindex(MS.index).values
    mod_cn = MS.groupby('CN')[list(MODS)].mean()
    mod_cn_z = (mod_cn - mod_cn.mean()) / mod_cn.std(ddof=0)
    mod_cn.rename(index=NAMES).assign(n_cells_with_LR=MS.groupby('CN').size().values).to_csv(f'{OUT}/cn_k10_module_scores.csv', float_format='%.4g')
    log(mod_cn.rename(index=NAMES).round(4).to_string())
    del lr; gc.collect()

    # ------------------------------------------------------------------------------------------------ figures
    order = list(range(K0))
    show_ct = [c for c in CT if c not in ('Cycling', 'Neuron')]
    # 1. composition heatmap
    fig, ax = plt.subplots(figsize=(13, 5.2))
    hm = enr.loc[order, show_ct].clip(-3, 3); hm.index = [NAMES[c] for c in order]
    sns.heatmap(hm, cmap='RdBu_r', center=0, vmin=-3, vmax=3, ax=ax, cbar_kws={'label': 'log2 Enrichment vs All Cells', 'shrink': 0.8},
                linewidths=0.3, linecolor='white')
    ax.set_xlabel('Cell Type', fontsize=11); ax.set_ylabel('Cellular Neighbourhood (k = 10)', fontsize=11)
    ax.tick_params(axis='x', labelsize=8.5); ax.tick_params(axis='y', labelsize=9)
    ax.set_title('Cellular Neighbourhood Composition, CD Resections (11 Sections, 4 Patients)', fontsize=12)
    fig.savefig(f'{OUT}/fig_cn_composition_heatmap.png', dpi=300, bbox_inches='tight'); plt.close(fig)
    # 2. forest plot: CMH OR activated vs low-suppressive per CN
    r10 = RES[RES.k == K0].copy()
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    cols_d = {'D6_tertile': '#7a3b69'}; off = {'D6_tertile': -0.15}
    for d in DEFS:
        g = r10[r10.definition == d].set_index('CN').loc[order]
        y = np.arange(K0) + off[d]; orv = np.log2(g.cmh_OR_activated_vs_low_suppressive.values)
        ax.scatter(orv, y, color=cols_d[d], s=45, zorder=3, label='Signature (top vs bottom tertile)')
        for yi, o, q, sgn in zip(y, orv, g.cmh_q.values, g.sections_activated_higher.astype(str) + '/' + g.sections_tested.astype(str)):
            ax.text(o + (0.06 if o >= 0 else -0.06), yi, stars(q), ha='left' if o >= 0 else 'right', va='center', fontsize=10, color=cols_d[d])
    ax.axvline(0, color='0.4', lw=0.8, ls='--')
    ax.set_yticks(range(K0)); ax.set_yticklabels([NAMES[c] for c in order], fontsize=8.5); ax.invert_yaxis()
    ax.set_xlabel('log2 CMH Odds Ratio (Activated vs Low-suppressive Tregs in CN)\n← Low-suppressive Enriched     Activated Enriched →', fontsize=10)
    ax.set_title('Treg State Enrichment per Cellular Neighbourhood', fontsize=12); ax.legend(fontsize=8.5, loc='lower right', frameon=False)
    sns.despine(ax=ax); fig.savefig(f'{OUT}/fig_cn_treg_state_forest.png', dpi=300, bbox_inches='tight'); plt.close(fig)
    # 3. module score heatmap per CN (z-scored across CNs)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    mz = mod_cn_z.loc[order]; mz.index = [NAMES[c] for c in order]
    sns.heatmap(mz, cmap='RdBu_r', center=0, annot=mod_cn.loc[order].values, fmt='.3f', annot_kws={'fontsize': 7.5}, ax=ax,
                cbar_kws={'label': 'Module Score (z across CNs)', 'shrink': 0.8}, linewidths=0.3, linecolor='white')
    ax.set_xlabel('LR Module', fontsize=10); ax.set_ylabel(''); ax.tick_params(axis='y', labelsize=8.5); ax.tick_params(axis='x', labelsize=9, rotation=30)
    ax.set_title('LR Module Scores per Cellular Neighbourhood\n(Mean Per-Cell LIANA Score; Numbers = Raw Means)', fontsize=11)
    fig.savefig(f'{OUT}/fig_cn_module_scores.png', dpi=300, bbox_inches='tight'); plt.close(fig)
    # 4. spatial maps
    pal = dict(zip(order, sns.color_palette('tab10', K0)))
    for sec, fovs, fn in [('11808_A4_Cos3', [str(i) for i in range(42, 54)], 'fig_cn_map_11808_A4_Cos3_fov42-53.png'),
                          ('54033_B7', None, 'fig_cn_map_54033_B7.png')]:
        s = sa[(sa.section == sec) & (sa.fov.isin(fovs) if fovs else True)]
        tr = PC[(PC.section == sec) & PC.D6_tertile.isin(['Early', 'Late'])]
        tr = tr[tr.index.isin(s.index)]
        fig, ax = plt.subplots(figsize=(8.5, 8.5) if fovs is None else (6.5, 8.5))
        for c in order:
            m = s.CN.values == c
            ax.scatter(s.x.values[m], s.y.values[m], s=1.2 if fovs is None else 4, color=pal[c], lw=0, rasterized=True, label=NAMES[c])
        for st, mk, col in [('Early', '^', 'black'), ('Late', 'o', 'white')]:
            t = s.loc[tr.index[tr.D6_tertile == st]]
            ax.scatter(t.x, t.y, s=14 if fovs is None else 26, marker=mk, facecolor=col, edgecolor='black', lw=0.5, zorder=4,
                       label=f"{'Activated' if st == 'Early' else 'Low-suppressive'} Treg (n = {len(t)})")
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"Cellular Neighbourhoods, {sec}{' FOVs 42-53' if fovs else ''}", fontsize=11)
        ax.legend(fontsize=7, markerscale=3, bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
        sb = 500 if fovs is None else 200
        x0, y0 = s.x.min(), s.y.min() - (s.y.max() - s.y.min()) * 0.03
        ax.plot([x0, x0 + sb], [y0, y0], 'k-', lw=2); ax.text(x0 + sb / 2, y0, f'{sb} µm', ha='center', va='top', fontsize=8)
        fig.savefig(f'{OUT}/{fn}', dpi=300, bbox_inches='tight'); plt.close(fig)
    pd.DataFrame({'section': sa.section.values, 'fov': sa.fov.values, 'x_um': sa.x.values, 'y_um': sa.y.values, 'cell_type': sa.ct.values,
                  'CN_k10': sa.CN.values + 1, 'CN_k8': LAB[8] + 1, 'CN_k12': LAB[12] + 1}, index=sa.index).to_parquet(f'{OUT}/cn_per_cell.parquet')
    pd.DataFrame({'k': list(ari), 'ARI_vs_k10': list(ari.values())}).to_csv(f'{OUT}/cn_k_sensitivity_ari.csv', index=False, float_format='%.3g')
    open(f'{OUT}/cn_run_log.txt', 'w').write('\n'.join(LOG))
    print('done')

    # =================================================================================================================
    # DISCOVERY (CD biopsies, WTX) and VALIDATION (resections projected onto biopsy CNs)
    # =================================================================================================================
    import h5py
    from anndata.io import read_elem
    WTX = '/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged'
    OD, OV = f'{OUT}/discovery_biopsy', f'{OUT}/validation_resection'; os.makedirs(OD, exist_ok=True); os.makedirs(OV, exist_ok=True)
    LOG2 = []
    def log2_(*a):
        s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG2.append(s_)

    with h5py.File(f'{WTX}/h5ad_obj/norm_anno_tregbin.h5ad', 'r') as h:
        bo = read_elem(h['obs']); bxy = h['obsm/spatial_fov'][:]
    bo = bo.set_index('cell', drop=False); bo.index = bo.index.astype(str)
    bo['section'] = bo['slide'].astype(str); bo['fov'] = bo['fov'].astype(str)
    bo['x'], bo['y'] = bxy[:, 0] * UM, bxy[:, 1] * UM
    bo['ct'] = network_anno(bo); bo = bo[bo.ct != 'DROP'].copy()
    COND = {'I0368': 'Inflamed CD', 'I0476': 'Inflamed CD', 'I0453': 'Non-inflamed CD', 'I0477': 'Healthy control'}
    bo['condition_lab'] = bo['patient'].astype(str).map(COND)
    assert (bo.groupby('section')['patient'].nunique() == 1).all(), 'a tissue piece spans >1 patient'

    # ---- shared vocabulary + mapping table
    mp = pd.concat([sa.assign(cohort='resection')[['cohort', 'cell_category', 'cell_type_general', 'cell_type_short', 'ct']],
                    bo.assign(cohort='biopsy')[['cohort', 'cell_category', 'cell_type_general', 'cell_type_short', 'ct']]])
    mp = mp.astype(str).groupby(['cohort', 'cell_category', 'cell_type_general', 'cell_type_short', 'ct']).size().reset_index(name='n_cells')
    mp = mp.rename(columns={'ct': 'shared_type'})
    VB, VR = set(bo['ct']), set(sa['ct'])
    mp['shared_type_in_both_cohorts'] = mp['shared_type'].isin(VB & VR)
    mp.to_csv(f'{OUT}/shared_celltype_vocabulary_mapping.csv', index=False)
    log2_('shared vocabulary: biopsy types', len(VB), '| resection types', len(VR), '| in both', len(VB & VR),
          '| only biopsy', sorted(VB - VR), '| only resection', sorted(VR - VB))
    CTS = sorted(VB | VR)


    def windows(df):
        code = pd.Categorical(df['ct'], categories=CTS).codes
        assert (code >= 0).all()
        rr, cc = [], []
        for s_, idx in df.groupby('section').indices.items():
            xy = df[['x', 'y']].values[idx]
            _, nn = cKDTree(xy).query(xy, k=KNN + 1)
            rr.append(np.repeat(idx, KNN)); cc.append(code[idx[nn[:, 1:]].ravel()])
        rr, cc = np.concatenate(rr), np.concatenate(cc)
        return sparse.csr_matrix((np.ones(len(rr), np.float32), (rr, cc)), shape=(len(df), len(CTS))).toarray() / KNN


    WB = windows(bo); WR = windows(sa)
    log2_('biopsy cells', len(bo), '| pieces', bo.section.nunique(), '| resection cells', len(sa))

    # ---- discovery CNs
    KMB, LB, LRP = {}, {}, {}
    for k in KS:
        KMB[k] = MiniBatchKMeans(n_clusters=k, random_state=0, batch_size=4096, n_init=10, max_iter=500).fit(WB)
        LB[k] = KMB[k].labels_; LRP[k] = KMB[k].predict(WR)
    ariB = {k: adjusted_rand_score(LB[K0], LB[k]) for k in KS if k != K0}
    log2_('biopsy ARI vs k=10:', {k: round(v, 3) for k, v in ariB.items()})
    overB = pd.Series(np.bincount(pd.Categorical(bo['ct'], categories=CTS).codes, minlength=len(CTS)) / len(bo), index=CTS)


    def name_cns(ct_series, lab, k, overall_, prefix):
        cellf = pd.crosstab(lab, ct_series.values, normalize='index').reindex(index=range(k), columns=CTS, fill_value=0)
        e = np.log2((cellf + 1e-4) / (overall_ + 1e-4))
        nm = {}
        for c in range(k):
            t = e.loc[c].sort_values(ascending=False); t = [x for x in t.index if t[x] > 0.5][:3] or [t.index[0]]
            nm[c] = f'{prefix}{c + 1}: ' + ' / '.join(t)
        return cellf, e, nm


    cellB, enrB, NB = name_cns(bo['ct'], LB[K0], K0, overB, 'CN')
    bo['CN'] = LB[K0]; sa['CNp'] = LRP[K0]
    log2_('biopsy CNs (k = 10):'); [log2_(f'  {NB[c]} (n = {(bo.CN == c).sum():,})') for c in range(K0)]
    enrB.rename(index=NB).to_csv(f'{OD}/cn_k10_celltype_log2_enrichment.csv', float_format='%.3g')
    cellB.rename(index=NB).to_csv(f'{OD}/cn_k10_celltype_fraction.csv', float_format='%.4g')
    for k in [8, 12]:
        name_cns(bo['ct'], LB[k], k, overB, 'CN')[1].rename(index=name_cns(bo['ct'], LB[k], k, overB, 'CN')[2]).to_csv(
            f'{OD}/cn_k{k}_celltype_log2_enrichment.csv', float_format='%.3g')

    # ---- Treg states (biopsy: all labelled Tregs of 10_treg_suppression_score.py)
    PCa = pd.read_csv(f'{P}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0)
    PB = PCa[PCa.cohort == 'biopsy'].copy(); PB.index = PB.index.astype(str)
    assert PB.index.isin(bo.index).all(); assert (bo.loc[PB.index, 'ct'] == 'Treg').all()
    log2_('biopsy labelled Tregs', len(PB), '| D6_tertile', PB.D6_tertile.value_counts().to_dict())


    def cn_test(df, lab_all, labels, k, dname, cohort):
        t = pd.DataFrame({'section': df.loc[labels.index, 'section'].values, 'state': labels[dname].values,
                          'CN': lab_all[df.index.get_indexer(labels.index)]})
        t = t[t.state.isin(['Early', 'Late'])]
        out = []
        for c in range(k):
            tabs, sgn = [], []
            for s_, g in t.groupby('section'):
                e, l = g[g.state == 'Early'], g[g.state == 'Late']
                a, b = int((e.CN == c).sum()), int((e.CN != c).sum()); cc, d = int((l.CN == c).sum()), int((l.CN != c).sum())
                if a + cc == 0 or b + d == 0: continue
                tabs.append(np.array([[a, b], [cc, d]], float))
                if len(e) >= 5 and len(l) >= 5: sgn.append(a / len(e) > cc / len(l))
            st = StratifiedTable([x + 0.5 if (x == 0).any() else x for x in tabs]) if tabs else None
            ne, nl = int(((t.state == 'Early') & (t.CN == c)).sum()), int(((t.state == 'Late') & (t.CN == c)).sum())
            out.append(dict(cohort=cohort, definition=dname, k=k, CN=c, n_activated_in_CN=ne, n_low_suppressive_in_CN=nl,
                            frac_activated_in_CN=ne / (t.state == 'Early').sum(), frac_low_suppressive_in_CN=nl / (t.state == 'Late').sum(),
                            cmh_OR_activated_vs_low_suppressive=st.oddsratio_pooled if st else np.nan,
                            cmh_p=st.test_null_odds(correction=True).pvalue if st else np.nan,
                            strata_activated_higher=int(np.sum(sgn)), strata_tested=len(sgn)))
        o = pd.DataFrame(out); ok = o.cmh_p.notna(); o['cmh_q'] = np.nan
        o.loc[ok, 'cmh_q'] = multipletests(o.loc[ok, 'cmh_p'], method='fdr_bh')[1]
        return o


    PR = PC.copy()   # resection Treg labels (D6_tertile), as above
    RB = pd.concat([cn_test(bo, LB[k], PB, k, d, 'biopsy (discovery)') for d in DEFS for k in KS], ignore_index=True)
    RV = pd.concat([cn_test(sa, LRP[k], PR, k, d, 'resection (validation, projected)') for d in DEFS for k in KS], ignore_index=True)
    for R_ in (RB, RV):
        R_['CN_name'] = [name_cns(bo['ct'], LB[k], k, overB, 'CN')[2][c] for k, c in zip(R_.k, R_.CN)]
    RB.to_csv(f'{OD}/cn_treg_state_enrichment_biopsy.csv', index=False, float_format='%.4g')
    RV.to_csv(f'{OV}/cn_treg_state_enrichment_resection_projected.csv', index=False, float_format='%.4g')
    log2_(RB[RB.k == K0].drop(columns=['k', 'cohort']).round(3).to_string(index=False))
    log2_(RV[RV.k == K0].drop(columns=['k', 'cohort']).round(3).to_string(index=False))

    # ---- descriptive by condition (biopsy)
    tb = pd.DataFrame({'condition': bo.loc[PB.index, 'condition_lab'].values, 'state': PB['D6_tertile'].values, 'CN': bo.loc[PB.index, 'CN'].values})
    tb = tb[tb.state.isin(['Early', 'Late'])]
    dc_ = (tb.groupby(['condition', 'state', 'CN']).size() / tb.groupby(['condition', 'state']).size()).rename('fraction').reset_index()
    dc_['n'] = tb.groupby(['condition', 'state', 'CN']).size().values
    dc_['CN_name'] = dc_['CN'].map(NB); dc_.to_csv(f'{OD}/cn_treg_state_by_condition_D6_tertile.csv', index=False, float_format='%.4g')
    log2_('Treg states by condition (n):', tb.groupby(['condition', 'state']).size().to_dict())

    # ---- concordance + composition correlation
    cellRp = pd.crosstab(LRP[K0], sa['ct'].values, normalize='index').reindex(index=range(K0), columns=CTS, fill_value=0)
    compB = pd.DataFrame(WB).groupby(LB[K0]).mean(); compB.columns = CTS
    compR = pd.DataFrame(WR).groupby(LRP[K0]).mean().reindex(range(K0)); compR.columns = CTS
    cc_rows = []
    for c in range(K0):
        from scipy.stats import pearsonr, spearmanr
        cc_rows.append(dict(CN=NB[c], n_biopsy_cells=int((LB[K0] == c).sum()), n_resection_cells=int((LRP[K0] == c).sum()),
                            frac_resection_cells=float((LRP[K0] == c).mean()), frac_biopsy_cells=float((LB[K0] == c).mean()),
                            pearson_window_composition=pearsonr(compB.loc[c], compR.loc[c].fillna(0))[0],
                            pearson_celltype_profile=pearsonr(cellB.loc[c], cellRp.loc[c])[0]))
    CC = pd.DataFrame(cc_rows)
    con = []
    for d in DEFS:
        b_, v_ = RB[(RB.definition == d) & (RB.k == K0)].set_index('CN'), RV[(RV.definition == d) & (RV.k == K0)].set_index('CN')
        for c in range(K0):
            ob, ov = b_.loc[c, 'cmh_OR_activated_vs_low_suppressive'], v_.loc[c, 'cmh_OR_activated_vs_low_suppressive']
            qb, qv = b_.loc[c, 'cmh_q'], v_.loc[c, 'cmh_q']
            same = np.sign(np.log(ob)) == np.sign(np.log(ov))
            con.append(dict(definition=d, CN=NB[c], biopsy_OR=ob, biopsy_q=qb, biopsy_pieces=f"{b_.loc[c, 'strata_activated_higher']}/{b_.loc[c, 'strata_tested']}",
                            resection_OR=ov, resection_q=qv, resection_sections=f"{v_.loc[c, 'strata_activated_higher']}/{v_.loc[c, 'strata_tested']}",
                            same_direction=bool(same), replicates=bool(same and qb < 0.05 and qv < 0.05),
                            direction_in_validation_if_discovery_sig=('same' if same else 'opposite') if qb < 0.05 else ''))
    CON = pd.DataFrame(con).merge(CC, on='CN', how='left')
    CON.to_csv(f'{OV}/cn_concordance_biopsy_vs_resection.csv', index=False, float_format='%.4g')
    log2_(CON.round(3).to_string(index=False))
    ari_proj = {k: adjusted_rand_score(LAB[k], LRP[k]) for k in KS}
    log2_('ARI resection de-novo CNs vs projected biopsy CNs (same k):', {k: round(v, 3) for k, v in ari_proj.items()})
    xt = pd.crosstab(pd.Series(LAB[K0]).map(NAMES).values, pd.Series(LRP[K0]).map(NB).values, normalize='index')
    xt.to_csv(f'{OV}/resection_denovo_vs_projected_CN_crosstab.csv', float_format='%.3g')
    pd.DataFrame({'k': list(ari_proj), 'ARI_resection_denovo_vs_projected': list(ari_proj.values()),
                  'ARI_biopsy_vs_k10': [ariB.get(k, 1.0) for k in KS]}).to_csv(f'{OV}/cn_ari_summary.csv', index=False, float_format='%.3g')
    modR = MS.assign(CNp=pd.Series(LRP[K0], index=sa.index).reindex(MS.index).values).groupby('CNp')[list(MODS)].mean()
    modR.rename(index=NB).to_csv(f'{OV}/cn_projected_module_scores.csv', float_format='%.4g')
    log2_(modR.rename(index=NB).round(4).to_string())


    # ---- key questions: (a) low-suppressive in OGN+RSPO3+ fibroblast CN, (b) activated in lymphoid B/Tfh or DC/CD9 Mac CN
    def key(e, kinds):
        return int(e[[x for x in kinds if x in e.columns]].mean(1).idxmax())
    kq = []
    for kind, kinds in [('OGN+RSPO3+ fibroblast CN', ['OGN_RSPO3 Fib']), ('B/Tfh lymphoid CN', ['B', 'Tfh']),
                        ('DC/CD9 Mac CN', [x for x in CTS if 'DC' in x] + ['CD9 Mac'])]:
        for k in KS:
            e = name_cns(bo['ct'], LB[k], k, overB, 'CN')[1]; c = key(e, kinds)
            for d in DEFS:
                b_ = RB[(RB.definition == d) & (RB.k == k) & (RB.CN == c)].iloc[0]; v_ = RV[(RV.definition == d) & (RV.k == k) & (RV.CN == c)].iloc[0]
                kq.append(dict(question=kind, k=k, definition=d, CN=b_.CN_name, biopsy_OR=b_.cmh_OR_activated_vs_low_suppressive, biopsy_q=b_.cmh_q,
                               biopsy_pieces=f'{b_.strata_activated_higher}/{b_.strata_tested}', resection_OR=v_.cmh_OR_activated_vs_low_suppressive,
                               resection_q=v_.cmh_q, resection_sections=f'{v_.strata_activated_higher}/{v_.strata_tested}'))
    KQ = pd.DataFrame(kq); KQ.to_csv(f'{OV}/cn_key_questions.csv', index=False, float_format='%.4g')
    log2_(KQ.round(3).to_string(index=False))

    # ---- figures
    show2 = [c for c in CTS if c not in ('Cycling', 'Neuron')]
    fig, ax = plt.subplots(figsize=(13, 5.2))
    hm = enrB.loc[range(K0), show2].clip(-3, 3); hm.index = [NB[c] for c in range(K0)]
    sns.heatmap(hm, cmap='RdBu_r', center=0, vmin=-3, vmax=3, ax=ax, cbar_kws={'label': 'log2 Enrichment vs All Cells', 'shrink': 0.8},
                linewidths=0.3, linecolor='white')
    ax.set_xlabel('Cell Type', fontsize=11); ax.set_ylabel('Cellular Neighbourhood (k = 10)', fontsize=11)
    ax.tick_params(axis='x', labelsize=8.5); ax.tick_params(axis='y', labelsize=9)
    ax.set_title('Cellular Neighbourhoods, CD Biopsies (Discovery; 4 Patients, 10 Tissue Pieces)', fontsize=12)
    fig.savefig(f'{OD}/fig_cn_composition_heatmap_biopsy.png', dpi=300, bbox_inches='tight'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
    for ax, (R_, ttl) in zip(axes, [(RB, 'Discovery: CD Biopsies (CMH by Tissue Piece)'), (RV, 'Validation: CD Resections (CMH by Section)')]):
        for d in DEFS:
            g = R_[(R_.definition == d) & (R_.k == K0)].set_index('CN').loc[range(K0)]
            y = np.arange(K0) + off[d]; ov_ = np.log2(g.cmh_OR_activated_vs_low_suppressive.values)
            ax.scatter(ov_, y, color=cols_d[d], s=45, zorder=3, label='Signature (top vs bottom tertile)')
            for yi, o_, q_ in zip(y, ov_, g.cmh_q.values):
                if np.isfinite(o_): ax.text(o_ + (0.08 if o_ >= 0 else -0.08), yi, stars(q_), ha='left' if o_ >= 0 else 'right', va='center', fontsize=10, color=cols_d[d])
        ax.axvline(0, color='0.4', lw=0.8, ls='--'); ax.set_title(ttl, fontsize=11); sns.despine(ax=ax)
        ax.set_xlabel('log2 CMH Odds Ratio, Activated vs Low-suppressive Tregs\n← Low-suppressive Enriched     Activated Enriched →', fontsize=9.5)
    axes[0].set_yticks(range(K0)); axes[0].set_yticklabels([NB[c] for c in range(K0)], fontsize=8.5); axes[0].invert_yaxis()
    axes[1].legend(fontsize=8.5, loc='lower right', frameon=False)
    fig.suptitle('Treg State Enrichment per Biopsy-Defined Cellular Neighbourhood (FDR: * <0.05, ** <0.01, *** <0.001)', fontsize=12, y=1.0)
    plt.tight_layout(); fig.savefig(f'{OV}/fig_cn_treg_state_forest_biopsy_vs_resection.png', dpi=300, bbox_inches='tight'); plt.close(fig)

    palB = dict(zip(range(K0), sns.color_palette('tab10', K0)))
    def cn_map(df, labs_, states, title, fn, s_cell, s_treg, sb):
        fig, ax = plt.subplots(figsize=(8, 8))
        for c in range(K0):
            m = labs_ == c
            ax.scatter(df.x.values[m], df.y.values[m], s=s_cell, color=palB[c], lw=0, rasterized=True, label=NB[c])
        for st, mk, col in [('Early', '^', 'black'), ('Late', 'o', 'white')]:
            t = df.loc[states.index[states == st].intersection(df.index)]
            ax.scatter(t.x, t.y, s=s_treg, marker=mk, facecolor=col, edgecolor='black', lw=0.5, zorder=4,
                       label=f"{'Activated' if st == 'Early' else 'Low-suppressive'} Treg (n = {len(t)})")
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=11)
        ax.legend(fontsize=7, markerscale=3, bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
        x0, y0 = df.x.min(), df.y.min() - (df.y.max() - df.y.min()) * 0.03
        ax.plot([x0, x0 + sb], [y0, y0], 'k-', lw=2); ax.text(x0 + sb / 2, y0, f'{sb} µm', ha='center', va='top', fontsize=8)
        fig.savefig(fn, dpi=300, bbox_inches='tight'); plt.close(fig)
    pc_ = bo[bo.section == 'i0476_slide1']
    cn_map(pc_, pc_['CN'].values, PB['D6_tertile'], 'Biopsy CNs, Inflamed CD (I0476, Piece i0476_slide1)', f'{OD}/fig_cn_map_biopsy_i0476_slide1.png', 3, 24, 200)
    rs = sa[(sa.section == '11808_A4_Cos3') & sa.fov.isin([str(i) for i in range(42, 54)])]
    cn_map(rs, rs['CNp'].values, PR['D6_tertile'], 'Projected Biopsy CNs, Resection 11808_A4_Cos3 FOVs 42-53', f'{OV}/fig_cn_map_resection_11808_A4_Cos3_fov42-53_projected.png', 4, 26, 200)
    bo[['section', 'patient', 'condition_lab', 'fov', 'x', 'y', 'ct']].assign(CN_k10=LB[K0] + 1, CN_k8=LB[8] + 1, CN_k12=LB[12] + 1).to_csv(
        f'{OD}/cn_per_cell_biopsy.csv.gz', float_format='%.5g')
    pd.DataFrame({'CNp_k10': LRP[K0] + 1, 'CNp_k8': LRP[8] + 1, 'CNp_k12': LRP[12] + 1}, index=sa.index).to_parquet(f'{OV}/cn_projected_per_cell_resection.parquet')
    open(f'{OV}/cn_discovery_validation_log.txt', 'w').write('\n'.join(LOG2))
    print('discovery/validation done')


# ======================================================================================================================
# Section treg_local_ligand_receptor: per-Treg local ligand availability, own receptor expression and local APC density
# Panels: Fig 5D, Supp 9D-F (writes treg_local_ligand_receptor.csv.gz)
# ======================================================================================================================
if _run('treg_local_ligand_receptor'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/treg_local_ligand_receptor.py'   # working folder of this step (relative paths resolve from here)
    """Local ligand -> Treg receptor -> Treg state, all 11 CD resection sections.
For every strict Treg:
  local ligand availability  mean within-section z of each ligand over ALL cells within 50 um (>= 5 neighbours): CXCL12, CXCL14 (stromal),
                             IL6 (inflammatory), TNF, CCL22 (myeloid / DC) - ligand side only, so no LIANA receptor-density confound.
  own receptor expression    CXCR4, CXCR3 (CXCL12 / CXCL14 receptors), IL6R, IL6ST, TNFRSF1B, CCR7 - within-section z.
  local APC density          DC, SELENOP+ Mac, CD9+ M2 Mac, Mo-Mac, B cell within 50 um (count and fraction of neighbours).
  act_z                      Treg Suppression Score = -signature_score_D6, z within section.
Tests (linear mixed models, section random intercept; BH within block):
  A correspondence  receptor z ~ matched local ligand (CXCR4 and CXCR3 ~ CXCL12 and CXCL14; IL6R and IL6ST ~ IL6; TNFRSF1B ~ TNF; CCR7 ~ CCL22),
                    with every unmatched ligand-receptor combination reported as the specificity background.
  B state           act_z ~ each local ligand.
  C co-occurrence   Tregs split by local CXCL12 tertile within section: local APC density and act_z per tertile, and
                    act_z ~ local CXCL12 adjusted for local APC density.
Outputs panels/fig4_and_related/B_resection_analyses/ligand_receptor_state/: treg_local_ligand_receptor.csv.gz (per-Treg table read
by the later sections and by 13-17), ligand_receptor_state_stats.csv, by_cxcl12_tertile.csv."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
    from anndata.io import read_elem
    from scipy import sparse
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/B_resection_analyses/ligand_receptor_state'; os.makedirs(OUT, exist_ok=True)
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50
    LIG = ['CXCL12', 'CXCL14', 'IL6', 'TNF', 'CCL22']
    REC = ['CXCR4', 'CXCR3', 'IL6R', 'IL6ST', 'TNFRSF1B', 'CCR7']
    MATCH = {'CXCR4': ['CXCL12', 'CXCL14'], 'CXCR3': ['CXCL12', 'CXCL14'], 'IL6R': ['IL6'], 'IL6ST': ['IL6'], 'TNFRSF1B': ['TNF'], 'CCR7': ['CCL22']}
    APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str); sa['apc'] = sa.ct.isin(APC)
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); pcr = PC[PC.cohort == 'resection']
    T = sa.loc[pcr.index, ['section', 'x', 'y']].copy()
    T['act_z'] = -pcr.signature_score_D6.values; T['act_z'] = T.groupby('section').act_z.transform(lambda v: (v - v.mean()) / v.std())
    for g in LIG: T['lig_' + g] = np.nan
    for g in REC: T['rec_' + g] = np.nan
    T['apc_n'] = np.nan; T['nb_n'] = np.nan
    for s_, f in FOLDER.items():
        sub = sa[sa.section == s_]
        with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
            ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
        nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
        gl = [g for g in LIG + REC if g in v]; keep = nm.isin(sub.index)
        E = pd.DataFrame(X[np.where(keep)[0]][:, [v.get_loc(g) for g in gl]].toarray(), index=nm[keep], columns=gl); del X; gc.collect()
        E = E.loc[E.index.intersection(sub.index)]; cs = sub.loc[E.index]
        Z = (E - E.mean()) / E.std().replace(0, np.nan)
        tg = T[T.section == s_]
        if not len(tg): continue
        tree = cKDTree(cs[['x', 'y']].values); nb = tree.query_ball_point(tg[['x', 'y']].values, r=RAD)
        apcv = cs.apc.values.astype(float)
        T.loc[tg.index, 'nb_n'] = [len(i) for i in nb]; T.loc[tg.index, 'apc_n'] = [apcv[i].sum() if len(i) else np.nan for i in nb]
        for g in [x for x in LIG if x in Z.columns]:
            vals = Z[g].values; T.loc[tg.index, 'lig_' + g] = [vals[i].mean() if len(i) >= 5 else np.nan for i in nb]
        common = tg.index.intersection(Z.index)
        for g in [x for x in REC if x in Z.columns]: T.loc[common, 'rec_' + g] = Z.loc[common, g].values
        print(s_, 'Tregs', len(tg), 'cells', len(cs), flush=True)
    T['apc_frac'] = T.apc_n / T.nb_n.replace(0, np.nan)
    T.to_csv(f'{OUT}/treg_local_ligand_receptor.csv.gz', float_format='%.4g')
    rows = []
    def mm(d, formula, term, block, outcome, predictor, extra=''):
        d = d.dropna(subset=[c for c in d.columns if c in ['y', 'x1', 'x2']])
        if len(d) < 100 or d.x1.std() == 0: return
        m = smf.mixedlm(formula, d, groups=d.section).fit(reml=True, method='lbfgs')
        rows.append(dict(block=block, outcome=outcome, predictor=predictor, adjusted=extra, coef=m.params[term], se=m.bse[term], p=m.pvalues[term], n=len(d)))
    for r_ in REC:
        if 'rec_' + r_ not in T or T['rec_' + r_].notna().sum() < 100: continue
        for l_ in LIG:
            if 'lig_' + l_ not in T or T['lig_' + l_].notna().sum() < 100: continue
            d = T.assign(y=T['rec_' + r_], x1=T['lig_' + l_])[['section', 'y', 'x1']]
            mm(d, 'y ~ x1', 'x1', 'A ligand -> receptor', r_, l_, 'matched' if l_ in MATCH.get(r_, []) else 'unmatched')
    for l_ in LIG:
        if 'lig_' + l_ not in T or T['lig_' + l_].notna().sum() < 100: continue
        d = T.assign(y=T.act_z, x1=T['lig_' + l_])[['section', 'y', 'x1']]
        mm(d, 'y ~ x1', 'x1', 'B ligand -> Treg activation', 'activation z', l_)
        d = T.assign(y=T.act_z, x1=T['lig_' + l_], x2=T.apc_frac)[['section', 'y', 'x1', 'x2']]
        mm(d, 'y ~ x1 + x2', 'x1', 'C adjusted for local APC fraction', 'activation z', l_, 'local APC fraction')
    S = pd.DataFrame(rows); S['q'] = S.groupby('block').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
    S.to_csv(f'{OUT}/ligand_receptor_state_stats.csv', index=False, float_format='%.3g')
    pd.set_option('display.width', 240); pd.set_option('display.max_rows', 80)
    for b, g in S.groupby('block'):
        print('==', b); print(g.sort_values('p')[['outcome', 'predictor', 'adjusted', 'coef', 'p', 'q', 'n']].round(4).to_string(index=False))
    T['cxcl12_tertile'] = T.groupby('section')['lig_CXCL12'].transform(lambda v: pd.qcut(v, 3, labels=['low', 'mid', 'high'], duplicates='drop'))
    tab = T.groupby('cxcl12_tertile').agg(n_tregs=('act_z', 'size'), mean_activation=('act_z', 'mean'), mean_apc_frac=('apc_frac', 'mean'),
                                          mean_apc_n=('apc_n', 'mean'), mean_neighbours=('nb_n', 'mean'))
    tab.to_csv(f'{OUT}/by_cxcl12_tertile.csv', float_format='%.4g'); print('== Tregs by local CXCL12 tertile (within section)'); print(tab.round(3).to_string())


# ======================================================================================================================
# Section treg_local_composition: per-Treg local cell-type composition within 50 um and region-specificity models
# Panels: Fig 5D, Supp 9D (writes treg_local_composition.csv.gz)
# ======================================================================================================================
if _run('treg_local_composition'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/treg_local_composition.py'   # working folder of this step (relative paths resolve from here)
    """Local cell-type composition of every strict Treg and the specificity of the APC x stromal-region interaction.
For every strict Treg (all 11 sections): fraction of each cell_type_general among ALL cells within 50 um (cKDTree per section, same
neighbourhoods as treg_local_ligand_receptor). Per stromal region (one stromal cell type, plus the 'CXCL12 ligand' reference using
lig_CXCL12): is the loss of the APC-associated rise of the Treg Suppression Score specific to OGN+RSPO3+ / CXCL12, or shared by
every stromal region?
  hi       upper vs lower within-section tertile of the local fraction (rank(method='first') then qcut, so ties at 0 don't collapse bins);
           middle tertile dropped.
  apcz     within-section quintile of apc_frac, z-scored over the analysed Tregs.
  cz       log1p(total raw counts), z within section (library size).
  model    act_z ~ hi * apcz + cz + (1|section), mixedlm REML lbfgs; BH across regions for hi and hi:apcz.
  robust   per-section OLS of the same model; count sections with negative hi:apcz (sections with < 60 Tregs in high+low skipped).
Outputs panels/fig4_and_related/B_resection_analyses/ligand_receptor_state/treg_local_composition.csv.gz (read by 15 and 17) and
territory_specificity.csv."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; LR = f'{D}/panels/fig4_and_related/B_resection_analyses/ligand_receptor_state'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50; MIN_SEC = 60
    TERR = ['OGN+RSPO3+ Fib', 'Myofibroblast', 'ADAMDEC1+ Fib', 'Crypt Top Fib', 'Activated Fib', 'PLA2G2A+ ECM Fib', 'SELENOP+ Fib', 'Pericyte', 'Endothelial']
    REF = 'CXCL12 ligand'
    # ---- per-Treg table and covariates (library size cz, APC-fraction quintile) ----
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    T['cz'] = T.groupby('section').cz.transform(lambda v: (v - v.mean()) / v.std())
    T['apc_q'] = T.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    # ---- local composition within 50 um (all cells, neighbourhoods of treg_local_ligand_receptor) ----
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
    assert T.index.isin(sa.index).all()
    CT = sorted(sa.ct.unique()); assert all(t in CT for t in TERR)
    comp = []
    for s_, tg in T.groupby('section'):
        cs = sa[sa.section == s_]; code = pd.Categorical(cs.ct, categories=CT).codes
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
        M = np.zeros((len(tg), len(CT)))
        for k, i in enumerate(nb): M[k] = np.bincount(code[i], minlength=len(CT))
        n = M.sum(1); C_ = pd.DataFrame(M / n[:, None], index=tg.index, columns=CT); C_.insert(0, 'nb_n', n); C_.insert(0, 'section', s_)
        assert np.array_equal(n, tg.nb_n.values), s_  # same neighbourhoods as treg_local_ligand_receptor.csv.gz
        comp.append(C_)
    C = pd.concat(comp).loc[T.index]; del sa; gc.collect()
    C.to_csv(f'{LR}/treg_local_composition.csv.gz', float_format='%.4g')
    # ---- stromal-region models ----
    def tert(v):  # within-section tertile; rank first so that ties (many 0 fractions) do not collapse bins
        return pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high'])
    rows, per_sec = [], []
    for terr in [REF] + TERR:
        frac = T.lig_CXCL12 if terr == REF else C[terr]
        t3 = frac.groupby(T.section).transform(tert)
        d = T.assign(t3=t3, frac=frac)
        d = d[d.t3.isin(['low', 'high']) & d.apc_q.notna()].copy(); d['hi'] = (d.t3 == 'high').astype(float)
        d['apcz'] = (d.apc_q - d.apc_q.mean()) / d.apc_q.std()
        d = d[['section', 'act_z', 'hi', 'apcz', 'cz', 'frac']].dropna().rename(columns={'act_z': 'y'})
        m = smf.mixedlm('y ~ hi * apcz + cz', d, groups=d.section).fit(reml=True, method='lbfgs')
        neg = tested = 0
        for s_, g in d.groupby('section'):
            if len(g) < MIN_SEC: continue
            f = smf.ols('y ~ hi * apcz + cz', g).fit(); tested += 1; neg += int(f.params['hi:apcz'] < 0)
            per_sec.append(dict(territory=terr, section=s_, n=len(g), coef_hi_x_apcz=f.params['hi:apcz'], p_hi_x_apcz=f.pvalues['hi:apcz']))
        hi_, lo_ = d[d.hi == 1].frac, d[d.hi == 0].frac
        rows.append(dict(territory=terr, n=len(d), mean_frac_high=hi_.mean(), mean_frac_low=lo_.mean(), pct_high_zero=100 * (hi_ == 0).mean(),
                         coef_hi=m.params['hi'], se_hi=m.bse['hi'], p_hi=m.pvalues['hi'],
                         coef_apcz=m.params['apcz'], se_apcz=m.bse['apcz'], p_apcz=m.pvalues['apcz'],
                         coef_hi_x_apcz=m.params['hi:apcz'], se_hi_x_apcz=m.bse['hi:apcz'], p_hi_x_apcz=m.pvalues['hi:apcz'],
                         sections_neg_hi_x_apcz=neg, sections_tested=tested))
        print(terr, 'n', len(d), 'hi:apcz {:.4f} p {:.3g}'.format(m.params['hi:apcz'], m.pvalues['hi:apcz']), f'neg {neg}/{tested}', flush=True)
    S = pd.DataFrame(rows)
    S['q_hi'] = multipletests(S.p_hi, method='fdr_bh')[1]; S['q_hi_x_apcz'] = multipletests(S.p_hi_x_apcz, method='fdr_bh')[1]
    S['sections_neg_of_tested'] = S.sections_neg_hi_x_apcz.astype(str) + '/' + S.sections_tested.astype(str)
    S = S[['territory', 'n', 'mean_frac_high', 'mean_frac_low', 'pct_high_zero', 'coef_hi', 'se_hi', 'p_hi', 'q_hi', 'coef_apcz', 'se_apcz', 'p_apcz',
           'coef_hi_x_apcz', 'se_hi_x_apcz', 'p_hi_x_apcz', 'q_hi_x_apcz', 'sections_neg_hi_x_apcz', 'sections_tested', 'sections_neg_of_tested']]
    S.to_csv(f'{LR}/territory_specificity.csv', index=False, float_format='%.3g')
    pd.set_option('display.width', 250)
    print(S[['territory', 'n', 'pct_high_zero', 'coef_hi', 'q_hi', 'coef_apcz', 'coef_hi_x_apcz', 'p_hi_x_apcz', 'q_hi_x_apcz', 'sections_neg_of_tested']].round(4).to_string(index=False))


# ======================================================================================================================
# Section cxcl12_expected_producer_field: local CXCL12 from expected stromal producers per Treg, and the CXCL12 models
# Panels: Fig 5E, Supp 9F (expected-producer CXCL12 field per Treg and the CXCL12-High Region statistics)
# ======================================================================================================================
if _run('cxcl12_expected_producer_field'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/cxcl12_expected_producer_field.py'   # working folder of this step (relative paths resolve from here)
    """Local CXCL12 restricted to expected producers, and the CXCL12 statistics of Fig 5E / Supp 9F.
Averaging CXCL12 over all neighbouring cells includes CosMx background in immune cells, so local CXCL12 is taken over the stromal
cell types expected to express CXCL12:
  EXPECTED_PRODUCERS (cell_type_general) = OGN+RSPO3+ Fib, ADAMDEC1+ Fib, Activated Fib, Crypt Top Fib, PLA2G2A+ ECM Fib, Myofibroblast,
  Pericyte, Endothelial.  Excluded: all SELENOP+ cells (SELENOP+ Mac and SELENOP+ Fib), all immune cells, all epithelium.
  Producer table  per cell type CXCL12 % expressing / mean (all 11 CD resection sections, cell set of treg_local_ligand_receptor),
          with the expected-producer flag; descriptive -> cxcl12_producer_table_exp.csv.
  Local field  lig_CXCL12_exp = mean of the within-section CXCL12 z over EXPECTED_PRODUCERS cells within 50 um, >= 3 required (else NaN);
          n_exp_producers_50um. The all-cell lig_CXCL12 and the neighbour counts of treg_local_ligand_receptor are recomputed
          as a guard. -> treg_local_cxcl12_producer_exp.csv.gz. Within-section tertiles of lig_CXCL12_exp use rank(method='first')
          because Tregs whose nearby producers all have zero CXCL12 share one value.
  Models  linear mixed models, section random intercept, REML, lbfgs; act_z = Treg Suppression Score (-signature_score_D6, z within
          section); cz = log1p(total raw counts), z within section; slopes also per within-section SD of x (coef_per_sd_x).
     (a) CXCR4 z ~ x;  ~ x + act_z;  ~ x + act_z + cz.     (b) act_z ~ x;  ~ x + cz; per-section OLS sign tests.
     (c) CXCL12-High vs CXCL12-Low Region (upper vs lower within-section tertile) x within-section APC-fraction quintile: per-quintile table
         (Mann-Whitney, BH over 5), pooled act_z ~ hi + C(apc_q), library-size-adjusted interaction act_z ~ hi * apcz + cz.
     (d) single-sender interactions act_z ~ hi * s_z + cz for DC, Macrophage, B cell, Plasma (BH over the 4), per-section negative counts.
     Sensitivity: Tregs with < 3 expected producers within 50 um as a separate group (S1) and merged into CXCL12-Low (S2).
Outputs (panels/fig4_and_related/B_resection_analyses/ligand_receptor_state/, suffix _exp): cxcl12_producer_table_exp.csv,
  treg_local_cxcl12_producer_exp.csv.gz, a_cxcr4_cxcl12_exp.csv, b_suppression_cxcl12_exp.csv, b_per_section_sign_tests_exp.csv,
  c_apc_matched_cxcl12_high_vs_low_exp.csv, c_pooled_apc_stratified_exp.csv, c_interaction_apc_by_cxcl12_countadj_exp.csv,
  which_sender_lost_exp.csv, sensitivity_no_producer_exp.csv."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
    from anndata.io import read_elem
    from scipy import sparse, stats
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; LR = f'{D}/panels/fig4_and_related/B_resection_analyses/ligand_receptor_state'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50; MIN_PROD = 3; MIN_SEC = 60
    EXPECTED_PRODUCERS = ['OGN+RSPO3+ Fib', 'ADAMDEC1+ Fib', 'Activated Fib', 'Crypt Top Fib', 'PLA2G2A+ ECM Fib', 'Myofibroblast', 'Pericyte', 'Endothelial']
    APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']
    SEND = {'dc': ['DC'], 'mac': ['SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac'], 'b': ['B cell'], 'plasma': ['Plasma']}
    NAME = {'dc': 'DC', 'mac': 'Macrophage', 'b': 'B cell', 'plasma': 'Plasma'}
    NOPROD = 'Fewer than 3 Expected Producers within 50 um'
    SS = 'Treg Suppression Score'
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    sig3 = lambda v: float(f'{v:.3g}') if pd.notna(v) else v
    def r3(df):
        df = df.copy()
        for c in df.columns:
            if df[c].dtype.kind == 'f' and not (df[c].dropna() % 1 == 0).all(): df[c] = df[c].map(sig3)   # counts are kept exact
        return df
    def wcsv(df, path, **kw):   # NFS share: plain python open('w')
        with open(path, 'w') as fh: df.to_csv(fh, **kw)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 150); pd.set_option('display.max_columns', 40)
    # ---------------- cells (cell set of treg_local_ligand_receptor) ----------------
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
    sa = sa[['section', 'x', 'y', 'ct']]
    missing = [c for c in EXPECTED_PRODUCERS if c not in set(sa.ct)]; assert not missing, missing
    T0 = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)   # per-Treg table (read-only)
    print('Tregs', len(T0), 'sections', T0.section.nunique(), flush=True)
    # ---------------- per-cell CXCL12, all 11 sections ----------------
    CELLS = {}
    for s_, f in FOLDER.items():
        sub = sa[sa.section == s_]
        with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
            ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
        nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsc() if sparse.issparse(X) else sparse.csc_matrix(X)
        keep = nm.isin(sub.index)
        e = pd.Series(np.asarray(X[:, v.get_loc('CXCL12')].toarray()).ravel()[keep], index=nm[keep]); del X, ob; gc.collect()
        e = e.loc[e.index.intersection(sub.index)]
        cs = sub.loc[e.index].copy(); cs['expr'] = e.values
        cs['z'] = (e.values - e.mean()) / e.std()                        # within-section z, as in treg_local_ligand_receptor
        CELLS[s_] = cs; print(s_, 'cells', len(cs), flush=True)
    allc = pd.concat(CELLS.values()); allc['pos'] = allc.expr > 0
    tot_pct = 100 * allc.pos.mean()
    P = allc.groupby('ct').agg(n_cells=('pos', 'size'), n_expressing=('pos', 'sum'), mean_expr=('expr', 'mean'))
    P['pct_expressing'] = 100 * P.n_expressing / P.n_cells
    P['mean_expr_in_expressing'] = allc[allc.pos].groupby('ct').expr.mean().reindex(P.index)
    P['fold_over_all_cells'] = P.pct_expressing / tot_pct
    P['expected_producer'] = np.where(P.index.isin(EXPECTED_PRODUCERS), 'True', 'False')
    P = P.sort_values(['expected_producer', 'pct_expressing'], ascending=[False, False])
    ex = allc[allc.ct.isin(EXPECTED_PRODUCERS)]
    ROWS = pd.DataFrame([dict(n_cells=len(allc), n_expressing=int(allc.pos.sum()), mean_expr=allc.expr.mean(), pct_expressing=tot_pct,
                              mean_expr_in_expressing=allc.expr[allc.pos].mean(), fold_over_all_cells=1.0, expected_producer=''),
                         dict(n_cells=len(ex), n_expressing=int(ex.pos.sum()), mean_expr=ex.expr.mean(), pct_expressing=100 * ex.pos.mean(),
                              mean_expr_in_expressing=ex.expr[ex.pos].mean(), fold_over_all_cells=ex.pos.mean() * 100 / tot_pct, expected_producer='')],
                        index=pd.Index(['ALL CELLS', 'ALL EXPECTED PRODUCERS'], name='ct'))
    Pout = pd.concat([ROWS, P])[['n_cells', 'n_expressing', 'pct_expressing', 'fold_over_all_cells', 'mean_expr', 'mean_expr_in_expressing', 'expected_producer']]
    Pout.index.name = 'cell_type_general'
    wcsv(r3(Pout), f'{LR}/cxcl12_producer_table_exp.csv')
    print('== CXCL12 by cell type (expected producers flagged)'); print(Pout.round(3).to_string(), flush=True)
    # ---------------- expected-producer local CXCL12 ----------------
    T = T0.copy()
    for c in ['lig_CXCL12_exp', 'n_exp_producers_50um', 'lig_all_check', 'nb_check'] + [k + '_n' for k in SEND]: T[c] = np.nan
    for s_, cs in CELLS.items():
        tg = T[T.section == s_]
        if not len(tg): continue
        nb = [np.asarray(i, dtype=int) for i in cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)]
        z = cs.z.values; pm = cs.ct.isin(EXPECTED_PRODUCERS).values; ctv = cs.ct.values
        npd = np.array([pm[i].sum() for i in nb])
        T.loc[tg.index, 'n_exp_producers_50um'] = npd
        T.loc[tg.index, 'lig_CXCL12_exp'] = [z[i[pm[i]]].mean() if k >= MIN_PROD else np.nan for i, k in zip(nb, npd)]
        T.loc[tg.index, 'lig_all_check'] = [z[i].mean() if len(i) >= 5 else np.nan for i in nb]
        T.loc[tg.index, 'nb_check'] = [len(i) for i in nb]
        for k, types in SEND.items():
            m_ = np.isin(ctv, types); T.loc[tg.index, k + '_n'] = [m_[i].sum() for i in nb]
    del CELLS, allc, ex; gc.collect()
    assert (T.nb_check == T.nb_n).all(), 'neighbourhood differs from treg_local_ligand_receptor'
    ok = T[['lig_CXCL12', 'lig_all_check']].dropna()
    assert len(ok) == T.lig_CXCL12.notna().sum() and np.allclose(ok.lig_CXCL12, ok.lig_all_check, rtol=1e-3, atol=1e-4), 'lig_CXCL12 not reproduced'
    print('guard: all-cell lig_CXCL12 of treg_local_ligand_receptor reproduced for {:,} Tregs'.format(len(ok)), flush=True)
    NEW = T[['lig_CXCL12_exp', 'n_exp_producers_50um']].copy(); NEW['n_exp_producers_50um'] = NEW.n_exp_producers_50um.astype(int)
    with open(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', 'wb') as fh: NEW.to_csv(fh, float_format='%.4g', compression='gzip')
    T['lig_CXCL12_exp'] = pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0).lig_CXCL12_exp.reindex(T.index)   # analyse as stored
    n_nan = int(T.lig_CXCL12_exp.isna().sum()); n_ok = len(T) - n_nan
    print('== Tregs with >= {} expected producers within 50 um: {:,} / {:,}; fewer (NaN): {:,} ({:.1f}%); 0 producers: {:,}'.format(
        MIN_PROD, n_ok, len(T), n_nan, 100 * n_nan / len(T), int((T.n_exp_producers_50um == 0).sum())))
    print('   median expected producers within 50 um: {:.0f}; per section NaN:'.format(T.n_exp_producers_50um.median()),
          T.groupby('section').lig_CXCL12_exp.apply(lambda v: f'{v.isna().sum()}/{len(v)}').to_dict())
    print('   Spearman lig_CXCL12 vs _exp {:.3f}'.format(T[['lig_CXCL12', 'lig_CXCL12_exp']].corr('spearman').iloc[0, 1]), flush=True)
    # ---------------- covariates ----------------
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    assert len(cnt) == len(T)
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    zs = lambda v: (v - v.mean()) / v.std()
    T['cz'] = T.groupby('section').cz.transform(zs)
    for k in SEND: T[k + '_frac'] = T[k + '_n'] / T.nb_n.replace(0, np.nan)
    for k in SEND: T[k + '_z'] = T.groupby('section')[k + '_frac'].transform(zs)
    T['apc_q'] = T.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    XV = {'exp': 'lig_CXCL12_exp'}
    for k in ['exp']:
        T['t_' + k] = T.groupby('section')[XV[k]].transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))
    def lmm(d, f, terms):
        m = smf.mixedlm(f, d, groups=d.section).fit(reml=True, method='lbfgs')
        return {t: dict(coef=m.params[t], se=m.bse[t], p=m.pvalues[t], n=len(d)) for t in terms}
    def perSD(d, coef):
        sd = (d.x - d.groupby('section').x.transform('mean')).std(); return dict(sd_x_within_section=sd, coef_per_sd_x=coef * sd)
    # ---------------- (a) CXCR4 ~ local CXCL12 ----------------
    rows = []
    for newlab, f in [('CXCR4 ~ CXCL12', 'y ~ x'), (f'CXCR4 ~ CXCL12 + {SS}', 'y ~ x + a'), (f'CXCR4 ~ CXCL12 + {SS} + library size', 'y ~ x + a + cz')]:
        res = {}
        for k, x in XV.items():
            d = T.assign(y=T.rec_CXCR4, x=T[x], a=T.act_z)[['section', 'y', 'x', 'a', 'cz']].dropna()
            res[k] = lmm(d, f, ['x'])['x']; res[k].update(perSD(d, res[k]['coef']))
            rows.append(dict(model=newlab, cxcl12=x, **res[k]))
    wcsv(r3(pd.DataFrame(rows)), f'{LR}/a_cxcr4_cxcl12_exp.csv', index=False)
    # ---------------- (b) suppression score ~ local CXCL12 ----------------
    rows = []
    for newlab, f in [(f'{SS} ~ CXCL12', 'y ~ x'), (f'{SS} ~ CXCL12 + library size', 'y ~ x + cz')]:
        res = {}
        for k, x in XV.items():
            d = T.assign(y=T.act_z, x=T[x])[['section', 'y', 'x', 'cz']].dropna()
            res[k] = lmm(d, f, ['x'])['x']; res[k].update(perSD(d, res[k]['coef'])); rows.append(dict(model=newlab, cxcl12=x, **res[k]))
    wcsv(r3(pd.DataFrame(rows)), f'{LR}/b_suppression_cxcl12_exp.csv', index=False)
    rows = []
    for lab, y, cov in [('CXCR4 ~ CXCL12', 'rec_CXCR4', None), (f'CXCR4 ~ CXCL12 + {SS} + library size', 'rec_CXCR4', ['act_z', 'cz']),
                        (f'{SS} ~ CXCL12', 'act_z', None), (f'{SS} ~ CXCL12 + library size', 'act_z', ['cz'])]:
        for k, x in XV.items():
            pe = []
            for s_, g in T.groupby('section'):
                c = [y, x] + (cov or []); g = g[c].dropna()
                if len(g) < 50 or g[x].std() == 0: continue
                dd = pd.DataFrame({'y': g[y].values, 'x': g[x].values}); f = 'y ~ x'
                for i, cv in enumerate(cov or []): dd[f'c{i}'] = g[cv].values; f += f' + c{i}'
                pe.append(smf.ols(f, dd).fit().params['x'])
            pe = np.array(pe); sg = np.sign(np.median(pe)); kk = int((np.sign(pe) == sg).sum())
            rows.append(dict(test=lab, cxcl12=x, n_sections=len(pe), concordant=f'{kk}/{len(pe)}', median_slope=np.median(pe),
                             sign_test_p=stats.binomtest(kk, len(pe), 0.5).pvalue))
    SGN = pd.DataFrame(rows); wcsv(r3(SGN), f'{LR}/b_per_section_sign_tests_exp.csv', index=False)
    # ---------------- (c) APC quintile x CXCL12-High / Low Region ----------------
    def per_quintile(tcol, lowset=('low',)):
        sub = T[T[tcol].isin(list(lowset) + ['high']) & T.apc_q.notna()]
        out = []
        for q, g in sub.groupby('apc_q'):
            hi, lo = g[g[tcol] == 'high'].act_z.dropna(), g[g[tcol].isin(lowset)].act_z.dropna()
            if len(hi) < 20 or len(lo) < 20: continue
            u = stats.mannwhitneyu(hi, lo); sd = np.sqrt((hi.var() * (len(hi) - 1) + lo.var() * (len(lo) - 1)) / (len(hi) + len(lo) - 2))
            out.append(dict(apc_quintile=int(q) + 1, mean_apc_frac=g.apc_frac.mean(), n_high=len(hi), n_low=len(lo),
                            suppression_cxcl12_high=hi.mean(), suppression_cxcl12_low=lo.mean(), delta_sd=(hi.mean() - lo.mean()) / sd, p=u.pvalue))
        out = pd.DataFrame(out); out['q'] = multipletests(out.p, method='fdr_bh')[1]; return out
    def pooled(tcol, lowset=('low',)):
        sub = T[T[tcol].isin(list(lowset) + ['high']) & T.apc_q.notna()]
        d = sub.assign(y=sub.act_z, hi=(sub[tcol] == 'high').astype(float))[['section', 'apc_q', 'y', 'hi']].dropna()
        return lmm(d, 'y ~ hi + C(apc_q)', ['hi'])['hi']
    def interaction(tcol, lowset=('low',)):
        d = T[T[tcol].isin(list(lowset) + ['high']) & T.apc_q.notna()].copy(); d['hi'] = (d[tcol] == 'high').astype(float)
        d['apcz'] = (d.apc_q - d.apc_q.mean()) / d.apc_q.std(); d = d[['section', 'act_z', 'hi', 'apcz', 'cz']].dropna().rename(columns={'act_z': 'y'})
        return lmm(d, 'y ~ hi * apcz + cz', ['hi', 'apcz', 'hi:apcz', 'cz'])
    Cq = {k: per_quintile('t_' + k) for k in XV}
    wcsv(r3(Cq['exp']), f'{LR}/c_apc_matched_cxcl12_high_vs_low_exp.csv', index=False)
    Po = {k: pooled('t_' + k) for k in XV}
    wcsv(pd.DataFrame([dict(estimate_sd=Po['exp']['coef'], se=Po['exp']['se'], p=Po['exp']['p'], n=Po['exp']['n'],
                            model=f'{SS} ~ CXCL12-High Region (expected producers) + C(APC quintile) + (1|section)')]).pipe(r3), f'{LR}/c_pooled_apc_stratified_exp.csv', index=False)
    IXr = {k: interaction('t_' + k) for k in XV}
    wcsv(pd.DataFrame([dict(term=t, **IXr['exp'][t]) for t in ['hi', 'apcz', 'hi:apcz', 'cz']]).pipe(r3), f'{LR}/c_interaction_apc_by_cxcl12_countadj_exp.csv', index=False)
    # ---------------- (d) single-sender interactions ----------------
    def senders(tcol, lowset=('low',)):
        d = T[T[tcol].isin(list(lowset) + ['high'])].copy(); d['hi'] = (d[tcol] == 'high').astype(float); d = d.rename(columns={'act_z': 'yy'})
        out = []
        for k in SEND:
            dd = d[['section', 'yy', 'hi', k + '_z', 'cz']].dropna().rename(columns={k + '_z': 's_z', 'yy': 'y'})
            m = lmm(dd, 'y ~ hi * s_z + cz', ['s_z', 'hi:s_z', 'hi'])
            neg = tested = 0
            for s_, g in dd.groupby('section'):
                if len(g) < MIN_SEC or g.s_z.std() == 0 or g.hi.nunique() < 2: continue
                tested += 1; neg += int(smf.ols('y ~ hi * s_z + cz', g).fit().params['hi:s_z'] < 0)
            out.append(dict(sender=NAME[k], m=m, neg=neg, tested=tested))
        qs = multipletests([o['m']['s_z']['p'] for o in out], method='fdr_bh')[1]; qi = multipletests([o['m']['hi:s_z']['p'] for o in out], method='fdr_bh')[1]
        for o, a, b in zip(out, qs, qi): o['s_q'], o['int_q'] = a, b
        return out
    SD = {k: senders('t_' + k) for k in XV}
    rows = []
    for n in SD['exp']:
        M = f'single sender: {SS} ~ hi * s_z + cz'
        rows += [dict(model=M, sender=n['sender'], term='s_z (sender effect in CXCL12-Low Region)', **n['m']['s_z'], q=n['s_q']),
                 dict(model=M, sender=n['sender'], term='hi:s_z (change in CXCL12-High Region)', **n['m']['hi:s_z'], q=n['int_q'],
                      sections_negative=n['neg'], sections_tested=n['tested']),
                 dict(model=M, sender=n['sender'], term='hi', **n['m']['hi'])]
    WSE = pd.DataFrame(rows)[['model', 'sender', 'term', 'coef', 'se', 'p', 'q', 'n', 'sections_negative', 'sections_tested']]
    wcsv(r3(WSE), f'{LR}/which_sender_lost_exp.csv', index=False)
    # ---------------- sensitivity: Tregs with < 3 expected producers ----------------
    T['grp'] = T.t_exp.astype(object).where(T.t_exp.notna(), NOPROD)
    T['grp_low'] = T.grp.replace({NOPROD: 'low'})
    SEN = []
    desc = T.groupby('grp').agg(n_tregs=('act_z', 'size'), mean_suppression_z=('act_z', 'mean'), mean_CXCR4_z=('rec_CXCR4', 'mean'), mean_apc_frac=('apc_frac', 'mean'),
                                mean_neighbours=('nb_n', 'mean'), mean_exp_producers=('n_exp_producers_50um', 'mean'), mean_lig_CXCL12_allcell=('lig_CXCL12', 'mean'))
    desc = desc.reindex([NOPROD, 'low', 'mid', 'high'])
    for g, r in desc.iterrows(): SEN.append(dict(analysis='S1 descriptives by group', group=g, **r.to_dict()))
    dg = T.assign(g=T.grp.map({NOPROD: 'a_noprod', 'low': 'b_low', 'mid': 'c_mid', 'high': 'd_high'}))
    for lab, y, f in [(f'S1 {SS} ~ group + library size', 'act_z', "y ~ C(g, Treatment('b_low')) + cz"),
                      (f'S1 CXCR4 ~ group + {SS} + library size', 'rec_CXCR4', "y ~ C(g, Treatment('b_low')) + act_z + cz")]:
        d = dg.assign(y=dg[y])[['section', 'y', 'g', 'act_z', 'cz']].dropna()
        m = smf.mixedlm(f, d, groups=d.section).fit(reml=True, method='lbfgs')
        for lv, nm_ in [('a_noprod', NOPROD), ('c_mid', 'mid'), ('d_high', 'high')]:
            t = f"C(g, Treatment('b_low'))[T.{lv}]"
            SEN.append(dict(analysis=lab, group=f'{nm_} vs low', coef=m.params[t], se=m.bse[t], p=m.pvalues[t], n=len(d)))
    s2 = interaction('grp_low')
    for t in ['hi', 'apcz', 'hi:apcz']: SEN.append(dict(analysis='S2 interaction, <3-producer Tregs merged into CXCL12-Low', group=t, **s2[t]))
    s2p = pooled('grp_low'); SEN.append(dict(analysis=f'S2 pooled {SS} ~ CXCL12-High + C(APC quintile), <3-producer merged into low', group='hi', **s2p))
    s2q = per_quintile('grp_low')
    for _, r in s2q.iterrows():
        SEN.append(dict(analysis='S2 per APC quintile, <3-producer merged into low (delta SD, BH over 5)', group=f'APC quintile {int(r.apc_quintile)}',
                        coef=r.delta_sd, p=r.p, q=r.q, n=int(r.n_high + r.n_low)))
    for o in senders('grp_low'):
        SEN.append(dict(analysis='S2 single sender, <3-producer merged into low', group=f"{o['sender']} hi:s_z", **o['m']['hi:s_z'], q=o['int_q'],
                        sections_negative=f"{o['neg']}/{o['tested']}"))
    SEN = pd.DataFrame(SEN); wcsv(r3(SEN), f'{LR}/sensitivity_no_producer_exp.csv', index=False)
    print('== per-section sign tests'); print(SGN.to_string(index=False))
    print('== _exp per-quintile table'); print(Cq['exp'].round(4).to_string(index=False))
    print('== sensitivity'); print(SEN.round(4).to_string(index=False))
