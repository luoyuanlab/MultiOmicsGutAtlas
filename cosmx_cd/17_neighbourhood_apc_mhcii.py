#!/usr/bin/env python
"""CD CosMx 6K resections: antigen-presenting cells (APCs) and MHC-II per cellular neighbourhood (CN).

Panels: Supp 9C (APC % of cells, adjusted per-APC MHC-II score and inside-minus-outside Treg Suppression Score per CN).
Inputs: resection combined/ tables and per-section norm_log1p.h5ad; cellular_neighborhoods_v2/cn_per_cell.parquet (11);
        lr_modules_cellchat/treg_receiver_module_scores_per_treg.csv.gz (13); A_fig4_panels/fig4_cn_treg_state_forest.csv (12);
        archive/supp_apc_abundance_by_neighbourhood.csv (16, section S9E_cxcl12_by_celltype).
Outputs: archive/ext_cn10_antigen_presentation_abundance.csv, A_fig4_panels/supp_cn10_apc_abundance_mhcii.csv,
         SuppFig9A-F_niche_cxcl12_apc/supp_cn10_apc_one_panel.{png,csv}.
Run: python 17_neighbourhood_apc_mhcii.py, or RUN_SECTIONS=<id>[,<id>...] (file order). APC_CACHE=<path.pkl> stores / re-uses the per-APC table.
"""
import os, gc, warnings
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from anndata.io import read_elem
from scipy import sparse, stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
warnings.filterwarnings('ignore'); warnings.simplefilter('ignore', ConvergenceWarning)

_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


RJ = '/share/fsmresfiles/UC/revision_jenny'
P4 = f'{RJ}/cd_treg_rebuild_rowfix/panels/fig4_and_related'; OUTA = f'{P4}/A_fig4_panels'
CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'
R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM, RAD = 0.12028, 50
FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
          '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
          '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
          '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
          '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
          '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
APC = ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac']           # SELENOP+ Mac is an APC, SELENOP+ Fib is not
APC_GRP = {'DC': 'DC', 'B cell': 'B cell', 'CD9+ M2 Mac': 'Macrophage', 'Mo-Mac': 'Macrophage', 'SELENOP+ Mac': 'Macrophage'}
APCG = ['DC', 'B cell', 'Macrophage']
HLA = ['HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DRA']   # antigen-presentation ligands (-> CD4)
CN78 = [1, 2, 3, 4, 5, 7, 8, 10]; OTH78 = [c for c in CN78 if c != 10]    # CN6 / CN9 omitted (too few Tregs)
zs = lambda v: (v - v.mean()) / v.std()
sig3 = lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
APCC = os.environ.get('APC_CACHE')


def load_apcs():
    """Per APC: section, type, x, y (um), log library size, MHC-II score (mean of the within-section z, over all cells, of the 9 HLA genes), CN."""
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', 'total_counts', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
    ab = []
    for s_, f in FOLDER.items():
        sub = sa[sa.section == s_]
        with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
            ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
        nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
        keep = np.where(nm.isin(sub.index))[0]
        E = X[keep][:, [v.get_loc(g) for g in HLA]].toarray().astype(np.float32); del X; gc.collect()
        idx = nm[keep]; cs = sub.loc[idx]
        sd = E.std(0); Z = np.where(sd > 0, (E - E.mean(0)) / np.where(sd > 0, sd, 1), np.nan)   # z over all cells of the section
        am = np.isin(cs.ct.values, APC)
        a = pd.DataFrame({'section': s_, 'apc_type': cs.ct.values[am], 'x': cs.x.values[am], 'y': cs.y.values[am],
                          'log_total_counts': np.log1p(cs.total_counts.values[am].astype(float)), 'mhc2': np.nanmean(Z[am], 1)}, index=idx[am])
        a['CN'] = cn.loc[a.index, 'CN_k10'].astype(int).values
        ab.append(a); print(s_, 'APCs', len(a), flush=True); del E, Z; gc.collect()
    return pd.concat(ab), sa[['x', 'y']]


# ======================================================================================================================
# Section cn10_apc_abundance: APCs near Tregs by CN, per-APC MHC-II score by CN (mixed models, BH over 7 contrasts), and robustness of the
# CN10 antigen-presentation ligand-supply model (APC threshold, leave-one-section-out, abundance-weighted supply)
# ======================================================================================================================
if _run('cn10_apc_abundance'):
    STEM = f'{P4}/archive/ext_cn10_antigen_presentation_abundance'
    MIN_TREG_CN = 10; MIN_APC_SEC = 10
    M = 'Antigen presentation'; CNS = CN78; OTH = OTH78; THRS = [1, 3, 5, 10]
    FIGCN = {1: 'CN1: Epithelial', 2: 'CN2: Plasma', 3: 'CN3: Myofib / Neuron', 4: 'CN4: DC / Mac', 5: 'CN5: InfMo-Mac',
             7: 'CN7: Activ Fib / Neut', 8: 'CN8: B / Tfh', 10: 'CN10: OGN⁺RSPO3⁺ Fib'}
    W = pd.read_csv(f'{OUTB}/treg_receiver_module_scores_per_treg.csv.gz', index_col=0)   # per-Treg table of 13 (CN, library size, ligand supply)
    W['lib_z'] = W.groupby('section').log_total_counts.transform(zs)
    W['n74'] = W['n|MHCII_APC'].astype(int)
    if APCC and os.path.exists(APCC):
        A, TT = pd.read_pickle(APCC)
    else:
        A0, XY = load_apcs(); ab, tb = [], []
        for s_ in FOLDER:
            a = A0[A0.section == s_].copy(); tg = W[W.section == s_]; tx = XY.loc[tg.index, ['x', 'y']].values
            nb = cKDTree(a[['x', 'y']].values).query_ball_point(tx, r=RAD)          # Treg -> APCs within 50 um
            t = pd.DataFrame(index=tg.index)
            t['n_chk'] = [len(i) for i in nb]
            t['L_chk'] = [a.mhc2.values[i].mean() if len(i) else np.nan for i in nb]
            t['sum_mhc2'] = [a.mhc2.values[i].sum() if len(i) else 0.0 for i in nb]
            for ty in APC: t[f'n_{ty}'] = [int((a.apc_type.values[i] == ty).sum()) for i in nb]
            tb.append(t)
            tcn = tg.CN.values; rev = cKDTree(tx).query_ball_point(a[['x', 'y']].values, r=RAD)   # APC -> Tregs within 50 um
            dist, nn = cKDTree(tx).query(a[['x', 'y']].values, k=1, distance_upper_bound=RAD)
            a['near_treg'] = np.isfinite(dist); a['nearest_treg_CN'] = np.where(np.isfinite(dist), tcn[np.minimum(nn, len(tcn) - 1)], -1)
            a['n_tregs_50um'] = [len(i) for i in rev]; a['n_cn10_tregs_50um'] = [int((tcn[i] == 10).sum()) if len(i) else 0 for i in rev]
            ab.append(a)
        A = pd.concat(ab); TT = pd.concat(tb).loc[W.index]; del A0, XY; gc.collect()
        if APCC: pd.to_pickle((A, TT), APCC)
    SAME74 = bool((TT.n_chk.values == W.n74.values).all())   # True if the per-Treg table of 13 used the same five APC types
    W['nAPC'] = TT.n_chk.values.astype(int); W['L_apc'] = TT.L_chk.values; W['S_sum'] = TT.sum_mhc2.values
    for g_ in APCG: W[f'n_{g_}'] = TT[[f'n_{ty}' for ty in APC if APC_GRP[ty] == g_]].sum(1).values
    A['cell_type'] = A.apc_type; A['apc_type'] = A.cell_type.map(APC_GRP)
    A['patient'] = A.section.str.split('_').str[0]
    A['lib_z'] = A.groupby('section').log_total_counts.transform(zs)
    A['y'] = A.groupby('section').mhc2.transform(zs)
    A['high'] = A.mhc2 >= A.groupby('section').mhc2.transform(lambda v: v.quantile(0.75))
    ROWS = []
    # ---- 1. APCs within 50 um of Tregs, by the Treg CN
    def dist_row(v):
        v = np.asarray(v)
        return dict(n_tregs=len(v), median_apcs_50um=np.median(v), q1_apcs_50um=np.percentile(v, 25), q3_apcs_50um=np.percentile(v, 75),
                    mean_apcs_50um=v.mean(), pct_ge1=100 * (v >= 1).mean(), pct_ge3=100 * (v >= 3).mean(), pct_ge5=100 * (v >= 5).mean(), pct_ge10=100 * (v >= 10).mean())
    for c in CNS:
        g = W[W.CN == c]; r = dict(analysis='1_apc_abundance', level='CN', CN=c, CN_label=FIGCN[c], section='all', **dist_row(g.nAPC))
        tot = g[[f'n_{t}' for t in APCG]].sum(); r.update({f'pct_{t}_of_nearby_apcs': 100 * tot[f'n_{t}'] / max(tot.sum(), 1) for t in APCG})
        r['n_sections_with_ge10_tregs'] = int((g.groupby('section').size() >= MIN_TREG_CN).sum())
        ROWS.append(r)
        for s_, gs in g.groupby('section'):
            ROWS.append(dict(analysis='1_apc_abundance', level='CN x section', CN=c, CN_label=FIGCN[c], section=s_, **dist_row(gs.nAPC)))
    f10 = A[A.n_cn10_tregs_50um > 0]                                          # unique APCs feeding CN10 Tregs, concentration of the links
    links = np.sort(f10.n_cn10_tregs_50um.values)[::-1]; top10 = links[:max(1, int(np.ceil(0.1 * len(links))))].sum() / links.sum()
    ROWS.append(dict(analysis='1_apc_abundance', level='unique APCs within 50 um of >= 1 CN10 Treg', CN=10, CN_label=FIGCN[10], section='all',
                     n_unique_apcs=len(f10), n_sections_with_such_apcs=f10.section.nunique(), median_unique_apcs_per_section=f10.groupby('section').size().median(),
                     min_unique_apcs_per_section=f10.groupby('section').size().min(), median_cn10_tregs_per_apc=np.median(links),
                     pct_links_from_top10pct_apcs=100 * top10, **{f'pct_{t}_of_nearby_apcs': 100 * (f10.apc_type == t).mean() for t in APCG}))

    # ---- 2. per-APC MHC-II score: y ~ C(CN, ref CN10) + lib_z (+ APC group) + (1 | section), APCs grouped by own CN or by nearest Treg's CN
    def fit(formula, d, term_prefix='C(grp, Treatment(10))'):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm']:
            try:
                f = smf.mixedlm(formula, d, groups=d.section).fit(reml=True, method=meth)
                ks = [k for k in f.params.index if k.startswith(term_prefix)]
                if all(np.isfinite(f.bse[k]) and f.bse[k] > 0 for k in ks): return f, meth
            except Exception: pass
        return None, 'failed'
    for gname, gcol in {'own_CN': 'CN', 'nearest_treg_CN': 'nearest_treg_CN'}.items():
        d0 = A[A[gcol].isin(CNS)].copy(); d0['grp'] = d0[gcol].astype(int)
        for tset in ['all APC types (type-adjusted)'] + APCG:
            alltypes = tset.startswith('all')
            d = d0 if alltypes else d0[d0.apc_type == tset]
            f, meth = fit('y ~ C(grp, Treatment(10)) + lib_z' + (' + C(apc_type)' if alltypes else ''), d)
            res = []                                                          # section dots: OLS residual per section, mean in CN minus mean in CN10
            for s_, gs in d.groupby('section'):
                rf = 'y ~ lib_z + C(apc_type)' if alltypes and gs.apc_type.nunique() > 1 else 'y ~ lib_z'
                res.append(gs[['grp']].assign(r=gs.y - smf.ols(rf, gs).fit().fittedvalues, section=s_))
            res = pd.concat(res); out = []
            for c in OTH:
                k = f'C(grp, Treatment(10))[T.{c}]'; dsec = []
                for s_, gs in res.groupby('section'):
                    a_, b_ = gs[gs.grp == c].r, gs[gs.grp == 10].r
                    if len(a_) >= MIN_APC_SEC and len(b_) >= MIN_APC_SEC: dsec.append((s_, a_.mean() - b_.mean()))
                b, se, p = (f.params[k], f.bse[k], f.pvalues[k]) if f is not None and k in f.params else (np.nan,) * 3
                out.append(dict(analysis='2_per_apc_mhc2', grouping=gname, apc_types=tset, CN=c, CN_label=FIGCN[c], section='all (mixed model)',
                                contrast='CN minus CN10 (negative = CN10 APCs higher)', estimate=b, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se, p=p,
                                n_apcs_in_CN=int((d.grp == c).sum()), n_apcs_CN10=int((d.grp == 10).sum()),
                                n_sections_ge10_both=len(dsec), sections_CN10_higher=f'{sum(x < 0 for _, x in dsec)}/{len(dsec)}',
                                optimiser=meth, lib_size_coef=f.params['lib_z'] if f is not None else np.nan))
                for s_, x in dsec:
                    out.append(dict(analysis='2_per_apc_mhc2', grouping=gname, apc_types=tset, CN=c, CN_label=FIGCN[c], section=s_, row_type='section',
                                    contrast='CN minus CN10 (type/lib-adjusted residual means)', delta=x,
                                    n_apcs_in_CN=int(((res.section == s_) & (res.grp == c)).sum()), n_apcs_CN10=int(((res.section == s_) & (res.grp == 10)).sum())))
            o = pd.DataFrame(out); mm = o.section == 'all (mixed model)'; okp = mm & o.p.notna()
            o['q'] = np.nan; o.loc[okp, 'q'] = multipletests(o.loc[okp, 'p'], method='fdr_bh')[1]
            o['row_type'] = np.where(mm, 'mixed_model', 'section'); ROWS += o.to_dict('records')
    for gname, gcol in {'own_CN': 'CN', 'nearest_treg_CN': 'nearest_treg_CN'}.items():
        d0 = A[A[gcol].isin(CNS)]
        for c in CNS:                                                         # descriptive: n APCs, mean score, % MHC-II-high (top quartile)
            g = d0[d0[gcol] == c]
            ROWS.append(dict(analysis='2_per_apc_descriptive', grouping=gname, apc_types='all', CN=c, CN_label=FIGCN[c], section='all', n_apcs=len(g),
                             mean_mhc2_score=g.mhc2.mean(), pct_mhc2_high_top_quartile=100 * g.high.mean(), n_sections_ge10_apcs=int((g.groupby('section').size() >= MIN_APC_SEC).sum()),
                             **{f'pct_{t}': 100 * (g.apc_type == t).mean() for t in APCG}))
            for t in APCG:
                gt = g[g.apc_type == t]
                ROWS.append(dict(analysis='2_per_apc_descriptive', grouping=gname, apc_types=t, CN=c, CN_label=FIGCN[c], section='all', n_apcs=len(gt),
                                 mean_mhc2_score=gt.mhc2.mean(), pct_mhc2_high_top_quartile=100 * gt.high.mean() if len(gt) else np.nan))
            for s_, gs in g.groupby('section'):
                ROWS.append(dict(analysis='2_per_apc_descriptive', grouping=gname, apc_types='all', CN=c, CN_label=FIGCN[c], section=s_, n_apcs=len(gs),
                                 mean_mhc2_score=gs.mhc2.mean(), pct_mhc2_high_top_quartile=100 * gs.high.mean()))
    dB = A[A.near_treg].copy(); dB['near_cn10'] = (dB.n_cn10_tregs_50um > 0).astype(float)   # APCs near >= 1 CN10 Treg vs near other-CN Tregs only
    for tset in ['all APC types (type-adjusted)'] + APCG:
        d = dB if tset.startswith('all') else dB[dB.apc_type == tset]
        f, meth = fit('y ~ near_cn10 + lib_z' + (' + C(apc_type)' if tset.startswith('all') else ''), d, 'near_cn10')
        b, se, p = (f.params['near_cn10'], f.bse['near_cn10'], f.pvalues['near_cn10']) if f is not None else (np.nan,) * 3
        sec = d.groupby('section').apply(lambda g: g[g.near_cn10 == 1].y.mean() - g[g.near_cn10 == 0].y.mean() if (g.near_cn10 == 1).sum() >= MIN_APC_SEC else np.nan).dropna()
        ROWS.append(dict(analysis='2_per_apc_mhc2', grouping='near CN10 Treg vs near other-CN Tregs only', apc_types=tset, CN=10, CN_label=FIGCN[10], row_type='mixed_model',
                         section='all (mixed model)', contrast='APCs within 50 um of >= 1 CN10 Treg minus APCs within 50 um of other-CN Tregs only (positive = CN10-adjacent higher)',
                         estimate=b, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se, p=p, q=p, n_apcs_in_CN=int(d.near_cn10.sum()), n_apcs_CN10=int(d.near_cn10.sum()),
                         n_apcs_other=int((d.near_cn10 == 0).sum()), pct_mhc2_high_near_cn10=100 * d[d.near_cn10 == 1].high.mean(),
                         pct_mhc2_high_near_other=100 * d[d.near_cn10 == 0].high.mean(), n_sections_ge10_both=len(sec),
                         sections_CN10_higher=f'{int((sec > 0).sum())}/{len(sec)}', optimiser=meth))

    # ---- 3-5. Treg-level ligand-supply model y ~ in_CN + lib_z + (1 | section), BH over CNs
    def treg_fit(col, d_all, c=10, drop=None):
        d = d_all[['section', 'CN', 'lib_z', col]].dropna().copy()
        if drop: d = d[d.section != drop]
        d['y'] = d.groupby('section')[col].transform(zs)
        secs, dls = [], []
        for s_, g in d.groupby('section'):
            inn, oth = g[g.CN == c], g[g.CN != c]
            if len(inn) >= MIN_TREG_CN and len(oth) >= MIN_TREG_CN: secs.append(s_); dls.append(inn.y.mean() - oth.y.mean())
        o = dict(n_sections=len(secs), n_tregs_in_CN=int(((d.CN == c) & d.section.isin(secs)).sum()), n_tregs_total=int(d.section.isin(secs).sum()),
                 n_tregs_scored_in_CN_all_sections=int((d.CN == c).sum()), n_tregs_scored_all=len(d))
        if len(secs) < 3: return {**o, 'estimate': np.nan}
        dd = d[d.section.isin(secs)].assign(in_cn=lambda x: (x.CN == c).astype(float))
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm('y ~ in_cn + lib_z', dd, groups=dd.section).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: break
            except Exception: f = None
        b, se, p = f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn']
        return {**o, 'estimate': b, 'ci_low': b - 1.96 * se, 'ci_high': b + 1.96 * se, 'p': p, 'optimiser': meth,
                'sections_agreeing': f'{int((np.sign(dls) == np.sign(b)).sum())}/{len(dls)}', 'section_deltas': dls, 'secs': secs}

    def all_cn(col, d_all, extra):
        rs = []
        for c in range(1, 11):
            r = treg_fit(col, d_all, c); r.update(CN=c, **extra); rs.append(r)
        o = pd.DataFrame(rs); okp = o.p.notna() if 'p' in o else pd.Series(False, index=o.index)
        o['q'] = np.nan; o.loc[okp, 'q'] = multipletests(o.loc[okp, 'p'], method='fdr_bh')[1]
        return o
    TH = []
    for thr in THRS:                                                          # 3. only Tregs with >= thr APCs within 50 um
        col = f'L_thr{thr}'
        W[col] = np.where(W.nAPC >= thr, W.L_apc.values, np.nan)
        if SAME74 and thr in (1, 3): W[col] = W[f'{M}|L' if thr == 1 else f'{M}|L min3']   # module columns of the per-Treg table of 13
        TH.append(all_cn(col, W, dict(analysis='3_threshold', min_apcs_50um=thr)))
    TH = pd.concat(TH, ignore_index=True)
    LO = []                                                                   # 4. leave one section out (threshold 1)
    for s_ in sorted(W.section.unique()):
        r = treg_fit('L_thr1', W, 10, drop=s_); r.update(analysis='4_leave_one_section_out', CN=10, dropped_section=s_, min_apcs_50um=1); LO.append(r)
    LO = pd.DataFrame(LO)
    SU = all_cn('S_sum', W, dict(analysis='5_abundance_weighted_sum', min_apcs_50um=0))   # 5. sum of APC MHC-II scores within 50 um, all Tregs
    for df in (TH, LO, SU):
        df['CN_label'] = df.CN.map(lambda c: FIGCN.get(c, f'CN{c} (omitted in Fig 4D)'))
        df['section'] = 'all (mixed model)'; df['row_type'] = 'mixed_model'
        df['section_deltas_z'] = df.get('section_deltas', pd.Series(index=df.index, dtype=object)).map(lambda v: ';'.join(f'{x:.3g}' for x in v) if isinstance(v, list) else '')
        df['sections_used'] = df.get('secs', pd.Series(index=df.index, dtype=object)).map(lambda v: ';'.join(v) if isinstance(v, list) else '')
        for k in ('section_deltas', 'secs'):
            if k in df: del df[k]
    OUT = pd.concat([pd.DataFrame(ROWS), TH, LO, SU], ignore_index=True)
    lead = ['analysis', 'row_type', 'grouping', 'apc_types', 'min_apcs_50um', 'dropped_section', 'CN', 'CN_label', 'section', 'level', 'contrast',
            'estimate', 'ci_low', 'ci_high', 'p', 'q', 'sections_agreeing', 'sections_CN10_higher', 'delta']
    OUT[[c for c in lead if c in OUT] + [c for c in OUT if c not in lead]].apply(sig3).to_csv(f'{STEM}.csv', index=False)


# ======================================================================================================================
# Section cn_apc_mhcii_table: APC share of cells per CN (table of 16) and model-adjusted per-APC MHC-II score per CN (read by S9C_cn_apc_mhcii)
# ======================================================================================================================
if _run('cn_apc_mhcii_table'):
    SRC78 = f'{P4}/archive/ext_cn10_antigen_presentation_abundance.csv'; SRC80 = f'{P4}/archive/supp_apc_abundance_by_neighbourhood.csv'
    STEM = f'{OUTA}/supp_cn10_apc_abundance_mhcii'
    A80 = pd.read_csv(SRC80)
    SALL = A80[A80.row_type == 'section']
    LAB = A80[A80.row_type == 'summary'].set_index('CN').CN_label.to_dict()
    # all-cell APC %: per section APCs / all cells, mean (95% CI) over the 11 sections; pooled for reference
    sec_all = SALL.groupby('section')[['n_APC', 'n_cells']].sum(); sec_all['apc_pct'] = 100 * sec_all.n_APC / sec_all.n_cells
    ALL_M = sec_all.apc_pct.mean(); ALL_SE = sec_all.apc_pct.std(ddof=1) / np.sqrt(len(sec_all))
    ALL_LO, ALL_HI = ALL_M + np.array([-1, 1]) * stats.t.ppf(0.975, len(sec_all) - 1) * ALL_SE
    ALL_POOL = 100 * sec_all.n_APC.sum() / sec_all.n_cells.sum()
    A78 = pd.read_csv(SRC78, low_memory=False)
    MM = A78[(A78.analysis == '2_per_apc_mhc2') & (A78.grouping == 'own_CN') & (A78.apc_types == 'all APC types (type-adjusted)')]
    MMo = MM[MM.row_type == 'mixed_model'].set_index('CN').loc[OTH78]
    if APCC and os.path.exists(APCC):
        A = pd.read_pickle(APCC); A = A[0] if isinstance(A, tuple) else A
    else:
        A = load_apcs()[0]
    A = A[['section', 'apc_type', 'log_total_counts', 'mhc2', 'CN']].copy()
    A['cell_type'] = A.apc_type; A['apc_type'] = A.cell_type.map(APC_GRP)
    A['lib_z'] = A.groupby('section').log_total_counts.transform(zs)
    A['y'] = A.groupby('section').mhc2.transform(zs)                               # z within section over all APCs
    dA = A[A.CN.isin(CN78)].copy(); dA['grp'] = dA.CN.astype(int)
    FRM = 'y ~ C(grp, Treatment(10)) + lib_z + C(apc_type)'
    # a fit is accepted only with a finite likelihood and finite fixed-effect SEs (lbfgs stops at the boundary here; its intercept is unusable)
    FIT, METH = None, None
    for meth in ['lbfgs', 'bfgs', 'powell', 'nm']:
        try:
            f_ = smf.mixedlm(FRM, dA, groups=dA.section).fit(reml=True, method=meth)
            if np.isfinite(f_.llf) and all(np.isfinite(f_.bse[k]) and f_.bse[k] > 0 for k in f_.fe_params.index): FIT, METH = f_, meth; break
        except Exception: pass
    KC = {c: f'C(grp, Treatment(10))[T.{c}]' for c in OTH78}
    qv = multipletests([FIT.pvalues[KC[c]] for c in OTH78], method='fdr_bh')[1]
    # model-adjusted mean per CN: intercept (CN10) + CN effect, lib_z = 0, APC-group effects weighted by the overall APC-group mix
    FE = FIT.fe_params; COV = FIT.cov_params().loc[FE.index, FE.index]
    WMIX = dA.apc_type.value_counts(normalize=True)
    KT = {g: f'C(apc_type)[T.{g}]' for g in WMIX.index if f'C(apc_type)[T.{g}]' in FE.index}
    ADJ = {}
    for c in CN78:
        L = pd.Series(0.0, index=FE.index); L['Intercept'] = 1.0
        if c != 10: L[KC[c]] = 1.0
        for g, k in KT.items(): L[k] = WMIX[g]
        m = float(L @ FE); se = float(np.sqrt(L @ COV @ L)); ADJ[c] = (m, m - 1.96 * se, m + 1.96 * se)
    p1 = A80[A80.row_type.isin(['summary', 'section'])][['row_type', 'CN', 'CN_label', 'section', 'n_cells', 'n_APC', 'mean_apc_pct', 'ci_low', 'ci_high',
          'mean_pct_DC', 'mean_pct_B cell', 'mean_pct_Macrophage', 'share_of_APCs_DC', 'share_of_APCs_B cell', 'share_of_APCs_Macrophage',
          'estimate_vs_CN10', 'estimate_ci_low', 'estimate_ci_high', 'p', 'q', 'used_in_mean']].copy()
    p1 = p1.rename(columns={'mean_apc_pct': 'apc_pct_of_cells'}).assign(panel='APCs as % of Cells')
    p1.loc[p1.row_type == 'summary', 'section'] = 'mean of sections'
    pa = pd.concat([pd.DataFrame(dict(row_type='section', CN='all', CN_label='All Cells', section=s, n_cells=r.n_cells, n_APC=r.n_APC, apc_pct_of_cells=r.apc_pct), index=[0])
                    for s, r in sec_all.iterrows()] + [pd.DataFrame(dict(row_type='summary', CN='all', CN_label='All Cells', section='mean of sections',
                    n_cells=sec_all.n_cells.sum(), n_APC=sec_all.n_APC.sum(), apc_pct_of_cells=ALL_M, ci_low=ALL_LO, ci_high=ALL_HI, pooled_apc_pct=ALL_POOL), index=[0])],
                   ignore_index=True).assign(panel='APCs as % of Cells')
    p2 = pd.DataFrame([dict(row_type='model_adjusted_mean', CN=c, section='all (mixed model)', adjusted_mean=ADJ[c][0], ci_low=ADJ[c][1], ci_high=ADJ[c][2],
                            estimate_vs_CN10=FIT.params[KC[c]] if c != 10 else 0.0, p=FIT.pvalues[KC[c]] if c != 10 else np.nan,
                            q=qv[OTH78.index(c)] if c != 10 else np.nan, q_78=MMo.loc[c, 'q'] if c != 10 else np.nan, n_apcs_in_CN=int((dA.grp == c).sum()),
                            optimiser=METH, weight_DC=WMIX['DC'], weight_B_cell=WMIX['B cell'], weight_Macrophage=WMIX['Macrophage']) for c in CN78])
    p2 = p2.assign(CN_label=p2.CN.map(LAB), panel='Per-APC MHC-II Score')
    OUT = pd.concat([p1, pa, p2], ignore_index=True)
    lead = ['panel', 'row_type', 'CN', 'CN_label', 'section']
    OUT[lead + [c for c in OUT if c not in lead]].apply(sig3).to_csv(f'{STEM}.csv', index=False)


# ======================================================================================================================
# Section S9C_cn_apc_mhcii: Supp 9C, one point per CN: APCs as % of cells (x), adjusted per-APC MHC-II score (y),
# inside-minus-outside Treg Suppression Score (fill, forest model of 12)
# ======================================================================================================================
if _run('S9C_cn_apc_mhcii'):
    SRC = f'{OUTA}/supp_cn10_apc_abundance_mhcii.csv'; SRCF = f'{OUTA}/fig4_cn_treg_state_forest.csv'
    STEM = f'{RJ}/SuppFig9A-F_niche_cxcl12_apc/supp_cn10_apc_one_panel'
    CNS = CN78
    VLIM = 0.4; CMAP = plt.get_cmap('RdBu_r'); NORM = matplotlib.colors.Normalize(-VLIM, VLIM)
    GREY, GREYT, GREYS = '#8c8c8c', '#666666', '#9a9a9a'
    EMPH_CN = (4, 5, 10)                                                            # highlighted: black-edged larger point, bold label
    S = pd.read_csv(SRC, low_memory=False)
    AP = S[(S.panel == 'APCs as % of Cells') & (S.row_type == 'summary')].copy()
    ALLR = AP[AP.CN.astype(str) == 'all'].iloc[0]; AP = AP[AP.CN.astype(str) != 'all']; AP['CN'] = AP.CN.astype(int); AP = AP.set_index('CN')
    MH = S[(S.panel == 'Per-APC MHC-II Score') & (S.row_type == 'model_adjusted_mean')].copy(); MH['CN'] = MH.CN.astype(int); MH = MH.set_index('CN')
    T = pd.DataFrame({'CN_label': AP.loc[CNS, 'CN_label'], 'apc_pct': AP.loc[CNS, 'apc_pct_of_cells'], 'apc_ci_low': AP.loc[CNS, 'ci_low'],
                      'apc_ci_high': AP.loc[CNS, 'ci_high'], 'q_apc': AP.loc[CNS, 'q'], 'mhc2': MH.loc[CNS, 'adjusted_mean'],
                      'mhc2_ci_low': MH.loc[CNS, 'ci_low'], 'mhc2_ci_high': MH.loc[CNS, 'ci_high'], 'q_mhc2': MH.loc[CNS, 'q_78'],
                      'mhc2_diff_vs_CN10': MH.loc[CNS, 'estimate_vs_CN10']})
    FO = pd.read_csv(SRCF).set_index('CN')
    T['treg_score'] = FO.loc[CNS, 'treg_suppression_score_z_in_minus_out'].values; T['treg_score_ci_low'] = FO.loc[CNS, 'ci_lo'].values
    T['treg_score_ci_high'] = FO.loc[CNS, 'ci_hi'].values; T['q_treg'] = FO.loc[CNS, 'mixed_q'].values
    T['sig_line'] = [f'Treg {stars(q)}' if q < 0.05 else '' for q in T.q_treg]
    T['fill'] = [matplotlib.colors.to_hex(CMAP(NORM(v))) for v in T.treg_score]
    ALL_PCT = ALLR.apc_pct_of_cells
    plt.rcParams.update({'font.size': 14, 'font.family': 'DejaVu Sans'})
    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=14)
    ax.axhline(0, color='#9e9e9e', lw=0.9, zorder=1)
    ax.axvline(ALL_PCT, color='#8c8c8c', lw=1.3, ls='--', zorder=1)
    for c in CNS:
        r = T.loc[c]; em = c in EMPH_CN
        ax.plot([r.apc_ci_low, r.apc_ci_high], [r.mhc2] * 2, color=GREY, lw=1.2, zorder=3, solid_capstyle='butt')
        ax.plot([r.apc_pct] * 2, [r.mhc2_ci_low, r.mhc2_ci_high], color=GREY, lw=1.2, zorder=3, solid_capstyle='butt')
        ax.plot([r.apc_pct], [r.mhc2], 'o', ms=14 if em else 12, mfc=r.fill, mec='#000000' if em else '#8c8c8c', mew=2.0 if em else 1.0, zorder=4)
    ax.set_xlim(0, 40); ax.set_ylim(-0.315, 0.225)
    ax.set_xlabel('APCs (% of All Cells in the CN)', fontsize=16); ax.set_ylabel('MHC-II Score per APC (z)', fontsize=16)
    ax.set_title('APC Abundance, MHC-II per APC and Treg Suppression Score', fontsize=16, loc='left', pad=6)
    # label anchors (x, top y, alignment, leader line to the point), placed by hand in free space
    POS = {10: (14.7, 0.212, 'center', False), 8: (31.8, 0.205, 'center', False), 3: (0.6, 0.060, 'left', False), 4: (32.6, 0.057, 'center', False),
           7: (19.0, -0.058, 'left', False), 5: (18.0, -0.150, 'left', False), 2: (9.0, -0.214, 'left', True), 1: (0.6, -0.262, 'left', True)}
    BGW = dict(facecolor='white', edgecolor='none', pad=0.8)
    for c in CNS:
        x, ytop, ha, lead = POS[c]
        ax.text(x, ytop, T.loc[c, 'CN_label'], ha=ha, va='top', fontsize=12, color=GREYT, zorder=7, bbox=BGW, fontweight='bold' if c in EMPH_CN else 'normal')
        if T.loc[c, 'sig_line']: ax.text(x, ytop - 0.0235, T.loc[c, 'sig_line'], ha=ha, va='top', fontsize=11, color=GREYS, zorder=7, bbox=BGW)
        if lead: ax.plot([T.loc[c, 'apc_pct']] * 2, [T.loc[c, 'mhc2_ci_low'] - 0.006, ytop + 0.004], color='#b0b0b0', lw=0.8, zorder=2)
    ax.text(ALL_PCT + 0.4, -0.307, 'All Cells', ha='left', va='bottom', fontsize=12, color='#6b6b6b', zorder=6)
    cax = ax.inset_axes([0.60, 0.075, 0.36, 0.032])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation='horizontal', ticks=[-0.4, -0.2, 0, 0.2, 0.4])
    cb.ax.tick_params(labelsize=12); cb.outline.set_linewidth(0.6)
    cb.ax.set_xticklabels(['−0.4', '−0.2', '0', '0.2', '0.4'])
    cb.ax.set_title('Treg Suppression Score (z)', fontsize=12.5, pad=4, color='#252525')
    fig.subplots_adjust(left=0.115, right=0.985, top=0.935, bottom=0.11)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    OUT = T.reset_index().rename(columns={'index': 'CN'}).assign(all_cells_apc_pct=ALL_PCT,
           metric='x = APCs as % of all cells in the CN (mean over sections, 95% CI t; q_apc = mixed model vs CN10, BH over 9); y = per-APC MHC-II score, '
              'model-adjusted mean (mixed model; 95% CI; q_mhc2 = q vs CN10, BH over 7); fill = treg_score: Treg Suppression Score, Tregs in minus outside the CN '
              '(CN forest LMM, q_treg BH over 8 tested CNs), RdBu_r +-0.4')
    OUT['stars_treg'] = [stars(q) for q in OUT.q_treg]
    OUT['stars_apc'] = [stars(q) if c != 10 else '' for c, q in zip(OUT.CN, OUT.q_apc)]; OUT['stars_mhc2'] = [stars(q) if c != 10 else '' for c, q in zip(OUT.CN, OUT.q_mhc2)]
    OUT.apply(sig3).to_csv(f'{STEM}.csv', index=False)
