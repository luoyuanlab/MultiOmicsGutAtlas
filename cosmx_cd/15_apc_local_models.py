#!/usr/bin/env python
"""CD CosMx 6K resections: Treg Suppression Score vs local cell types, with APC and DC models (11 sections, 4 patients, 6,554 strict Tregs).

Panels   Fig 5D joint APC-type model (left) and DC effect by within-section tertile of local OGN+RSPO3+ fibroblasts (right);
         Supp 9D score per SD of each cell type's local fraction; Results numbers (cd_spatial_results_all.csv, cd_interaction_selenop_tests.csv)
Inputs   combined/all_samples.parquet, combined/raw_counts_strict_tregs.csv; ligand_receptor_state/treg_local_*.csv.gz (11);
         A_fig4_panels/fig4_cn_treg_state_forest.csv (12)
Outputs  updated_figures/ (cd_spatial_results_all.csv, cd_interaction_selenop_tests.csv); Fig4A-F_treg_suppression_niches/ (Fig 5D panel,
         fig4d_dc_by_ogn_tertile_checks.csv, fig4d_dc_slope_by_ogn_tertile_per_section.csv); SuppFig9A-F_niche_cxcl12_apc/ (Supp 9D panel)
Sections spatial_results_table -> interaction_tests -> F5D_dc_by_ogn_tertile -> F5D_apc_types_dc_ogn -> S9D_celltype_forest
Run      python <file>, or RUN_SECTIONS=<id>[,<id>] python <file>
Models   act_z = Treg Suppression Score (z within section); x_z = local fraction among all cells within 50 um, z within section;
         cz = log1p(library size), z within section. Section-clustered SEs use t(G - 1) inference (G = 11 sections).
"""
import os, gc, warnings
import numpy as np, pandas as pd, pyarrow.parquet as pq
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')
_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
LR = f'{P4}/B_resection_analyses/ligand_receptor_state'; R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
UPD = '/share/fsmresfiles/UC/revision_jenny/updated_figures'; FIG = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches'
UM, RAD = 0.12028, 50
SUB = ['SELENOP+MMP9+PLA2G2D+ Mac', 'SELENOP+MMP9+SPP1+ Mac', 'SELENOP+SIGLEC1+ Mac', 'SELENOP+XIST- Mac']   # the 4 SELENOP+ Mac labels
# APC types of Fig 5D: key -> (source table, label(s), display name)
APC_TYPES = {'dc': ('comp', 'DC', 'DC'),
             's_mp': ('fine', 'SELENOP+MMP9+PLA2G2D+ Mac', 'SELENOP⁺MMP9⁺PLA2G2D⁺ Mac'),
             's_ms': ('fine', 'SELENOP+MMP9+SPP1+ Mac', 'SELENOP⁺MMP9⁺SPP1⁺ Mac'),
             's_sgxs': ('fine', ['SELENOP+SIGLEC1+ Mac', 'SELENOP+XIST- Mac'], 'SELENOP⁺ Mac (SIGLEC1⁺/XIST⁻)'),
             'cd9': ('comp', 'CD9+ M2 Mac', 'CD9⁺ M2 Mac'),
             'momac': ('comp', 'Mo-Mac', 'Mo-Mac'),
             'bc': ('comp', 'B cell', 'B cell')}
K = list(APC_TYPES); TT = ['bottom', 'middle', 'top']; TL = dict(bottom='Low', middle='Mid', top='High')


def zs(v):
    s = v.std()
    return (v - v.mean()) / s if s > 0 else v * 0.0


def load_tregs():
    """Per-Treg table and local composition (11), in the Treg order of the raw-count table, with library size cz."""
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    C = pd.read_csv(f'{LR}/treg_local_composition.csv.gz', index_col=0).loc[T.index]
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    T['cz'] = T.groupby('section').cz.transform(zs)
    return T, C.loc[T.index]


def selenop_subtypes(T):
    """Fraction of each SELENOP+ Mac subtype among all cells within 50 um (neighbour definition of the composition table)."""
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['cell_type'] = sa.cell_type.astype(str)
    F = pd.DataFrame(0.0, index=T.index, columns=SUB)
    for s_, tg in T.groupby('section'):
        cs = sa[sa.section == s_]; code = pd.Categorical(cs.cell_type, categories=SUB).codes      # -1 = not a SELENOP+ Mac subtype
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
        M = np.zeros((len(tg), len(SUB)))
        for k, i in enumerate(nb):
            cc = code[i]; M[k] = np.bincount(cc[cc >= 0], minlength=len(SUB))
        F.loc[tg.index] = M / np.array([len(i) for i in nb])[:, None]
    return F


def apc_frame():
    """Analysis frame of Fig 5D: local fraction (_f) and within-section z of the 7 APC types and of OGN+RSPO3+ Fib."""
    T, COMP = load_tregs(); F = selenop_subtypes(T)
    for k in K:
        src, lab, _ = APC_TYPES[k]
        T[k + '_f'] = COMP[lab].values if src == 'comp' else F[lab].sum(1).values if isinstance(lab, list) else F[lab].values
        T[k] = T.groupby('section')[k + '_f'].transform(zs)
    T['ogn_f'] = COMP['OGN+RSPO3+ Fib'].values; T['ogn'] = T.groupby('section').ogn_f.transform(zs)
    return T[['section', 'act_z', 'cz', 'ogn', 'ogn_f'] + K + [k + '_f' for k in K]].dropna()


def clfit(dd, rhs, w=None):
    """OLS (or WLS) with section fixed effects and section-clustered SEs."""
    g = pd.factorize(dd.section)[0]
    if w is None: return smf.ols(f'act_z ~ {rhs} + C(section)', dd).fit(cov_type='cluster', cov_kwds=dict(groups=g))
    return smf.wls(f'act_z ~ {rhs} + C(section)', dd, weights=w).fit(cov_type='cluster', cov_kwds=dict(groups=g))


def summ(o, t, g):
    b, se = o.params[t], o.bse[t]; tq = stats.t.ppf(0.975, g - 1)
    return dict(estimate=b, lo=b - tq * se, hi=b + tq * se, p=2 * stats.t.sf(abs(b / se), g - 1))


def agree(dd, rhs, t, b):
    """Sections whose per-section OLS coefficient has the sign of the pooled estimate b."""
    a = n = 0
    for _, gd in dd.groupby('section'):
        if len(gd) < 10 or any(gd[c].std() == 0 for c in t.split(':')): continue
        n += 1; a += int(np.sign(smf.ols(f'act_z ~ {rhs}', gd).fit().params[t]) == np.sign(b))
    return f'{a}/{n}'


def tert(v):  # within-section tertile; rank first so that ties (many 0 fractions) do not collapse bins
    return pd.qcut(v.rank(method='first'), 3, labels=['bottom', 'middle', 'top'])


def dc_by_ogn(d):
    """DC x OGN+RSPO3+ interaction and the DC slope within each within-section OGN+RSPO3+ tertile (adds column t to d)."""
    n_s = d.section.nunique()
    IX = summ(clfit(d, 'dc * ogn + cz'), 'dc:ogn', n_s); IX['agree'] = agree(d, 'dc * ogn + cz', 'dc:ogn', IX['estimate'])
    d['t'] = d.groupby('section').ogn_f.transform(tert); TER = {}
    for tt in TT:
        dt = d[d.t == tt]; g = dt.section.nunique()
        TER[tt] = summ(clfit(dt, 'dc + cz'), 'dc', g); TER[tt]['n'] = len(dt); TER[tt]['agree'] = agree(dt, 'dc + cz', 'dc', TER[tt]['estimate'])
    return IX, TER


# ======================================================================================================================
# Section spatial_results_table: all-feature table of the resection spatial results (neighbourhoods, local cell types, local ligands, APC x CXCL12)
#   separate  act_z ~ x_z + cz + (1|section) (MixedLM REML); BH across the block; clustered P and per-section sign agreement
#   joint     act_z ~ all x_z of the block + cz + (1|section), fitted only if max |Pearson r| between the x_z < 0.6; BH across the model
# ======================================================================================================================
if _run('spatial_results_table'):
    MIN_PREV, MAX_R = 5.0, 0.6
    T, COMP = load_tregs()
    T = T.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left')
    N_ALL = len(T)
    CT_ALL = sorted([c for c in COMP.columns if c not in ('section', 'nb_n', 'Treg')], key=str.lower)
    PREV = {c: 100 * (COMP[c] > 0).mean() for c in CT_ALL}
    CTS = [c for c in CT_ALL if PREV[c] >= MIN_PREV]
    LIGS = sorted([c for c in T.columns if c.startswith('lig_') and c != 'lig_CXCL12_exp'], key=str.lower)
    LIGS_ALL = LIGS + ['lig_CXCL12_exp']
    X = pd.DataFrame(index=T.index); X['section'] = T.section; X['act_z'] = T.act_z; X['cz'] = T.cz
    key = {}
    for i, c in enumerate(CT_ALL): key[c] = f'ct{i}'; X[f'ct{i}'] = COMP[c].values
    for i, c in enumerate(LIGS_ALL): key[c] = f'lg{i}'; X[f'lg{i}'] = T[c].values
    for c in CT_ALL + LIGS_ALL: X[key[c]] = X.groupby('section')[key[c]].transform(zs)

    def lmm(d, f, terms):
        m = smf.mixedlm(f, d, groups=d.section).fit(reml=True, method='lbfgs'); ci = m.conf_int()
        return {t: dict(coef=m.params[t], lo=ci.loc[t, 0], hi=ci.loc[t, 1], p=m.pvalues[t], n=len(d)) for t in terms}

    def clustered_p(d, f, t):
        o = smf.ols(f, d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]), use_t=True)
        return 2 * stats.t.sf(abs(o.params[t] / o.bse[t]), d.section.nunique() - 1)

    def sec_agree(d, f, t, sign):
        a, n = 0, 0
        for _, g in d.groupby('section'):
            if g[t].std() == 0 or len(g) < 5: continue
            n += 1; a += int(np.sign(smf.ols(f, g).fit().params[t]) == sign)
        return a, n

    def screen(feats):
        out = {}
        for c in feats:
            k = key[c]; d = X[['section', 'act_z', 'cz', k]].dropna(); f = f'act_z ~ {k} + cz'
            r = lmm(d, f, [k])[k]; r['agree'], r['nsec'] = sec_agree(d, f, k, np.sign(r['coef'])); r['cl_p'] = clustered_p(d, f, k)
            out[c] = r
        for c, qq in zip(feats, multipletests([out[c]['p'] for c in feats], method='fdr_bh')[1]): out[c]['q'] = qq
        return out

    def joint(feats):
        ks = [key[c] for c in feats]; d = X[['section', 'act_z', 'cz'] + ks].dropna()
        Rz = d[ks].corr().values
        if np.abs(Rz[np.triu_indices(len(ks), 1)]).max() >= MAX_R: return None
        r = lmm(d, 'act_z ~ ' + ' + '.join(ks) + ' + cz', ks)
        q = multipletests([r[k]['p'] for k in ks], method='fdr_bh')[1]
        return {c: dict(r[key[c]], q=qq) for c, qq in zip(feats, q)}

    RT = screen(CTS); JT = joint(CTS)
    RL = screen(LIGS_ALL); JL = joint(LIGS)   # lig_CXCL12_exp redefines lig_CXCL12 (r ~ 0.8): separate model only
    # APC x CXCL12: region x APC-quintile interaction and difference-in-differences of the rise from the lowest to the highest APC quintile
    A = T[['section', 'act_z', 'cz', 'apc_frac', 'lig_CXCL12_exp']].copy()
    A['apc_q'] = A.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    A['t_exp'] = A.groupby('section').lig_CXCL12_exp.transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))
    di = A[A.t_exp.isin(['low', 'high']) & A.apc_q.notna()].copy(); di['hi'] = (di.t_exp == 'high').astype(float)
    di['qz'] = (di.apc_q - di.apc_q.mean()) / di.apc_q.std(); di = di[['section', 'act_z', 'hi', 'qz', 'cz']].dropna()
    ri = lmm(di, 'act_z ~ hi * qz + cz', ['hi:qz'])['hi:qz']
    ri['agree'], ri['nsec'] = sec_agree(di.assign(hixqz=di.hi * di.qz), 'act_z ~ hi + qz + hixqz + cz', 'hixqz', np.sign(ri['coef']))
    ri['cl_p'] = clustered_p(di.assign(hixqz=di.hi * di.qz), 'act_z ~ hi + qz + hixqz + cz', 'hixqz')
    dd = A[A.t_exp.isin(['low', 'high']) & A.apc_q.isin([0, 4])].copy()
    dd['hi'] = (dd.t_exp == 'high').astype(float); dd['top'] = (dd.apc_q == 4).astype(float); dd = dd[['section', 'act_z', 'hi', 'top', 'cz']].dropna()
    rd = lmm(dd, 'act_z ~ hi * top + cz', ['hi:top'])['hi:top']; dd['hixtop'] = dd.hi * dd.top
    rd['agree'], rd['nsec'] = sec_agree(dd, 'act_z ~ hi + top + hixtop + cz', 'hixtop', np.sign(rd['coef']))
    rd['cl_p'] = clustered_p(dd, 'act_z ~ hi + top + hixtop + cz', 'hixtop')
    rise = {}
    for g in ['low', 'high']:
        s = A[(A.t_exp == g)]; q1, q5 = s[s.apc_q == 0].act_z, s[s.apc_q == 4].act_z
        rise[g] = (q5.mean() - q1.mean(), np.sqrt(q1.sem() ** 2 + q5.sem() ** 2))
    # assemble
    COLS = ['block', 'feature', 'model', 'n_tregs', 'effect_sd', 'ci_low', 'ci_high', 'p', 'q', 'sections_same_sign', 'clustered_p', 'joint_effect', 'joint_q', 'note']
    rows = []
    for r in pd.read_csv(f'{P4}/A_fig4_panels/fig4_cn_treg_state_forest.csv').sort_values('CN').itertuples():
        tested = bool(r.tested)
        rows.append(dict(block='neighbourhood', feature=r.CN_name, model='LMM: score ~ in_CN + (1|section) (Fig 4 forest)',
                         n_tregs=N_ALL, effect_sd=r.treg_suppression_score_z_in_minus_out, ci_low=r.ci_lo, ci_high=r.ci_hi, p=r.mixed_p,
                         q=r.mixed_q if tested else np.nan, sections_same_sign=f'{int(r.sections_agreeing)}/{int(r.sections)}', clustered_p=np.nan,
                         joint_effect=np.nan, joint_q=np.nan,
                         note=(f'{int(r.n_tregs_in_CN):,} Tregs in CN; effect = in-minus-out difference in score (SD of score); sections = those with >= 3 Tregs in CN'
                               + ('' if tested else f'; not tested (< 100 Tregs or < 5 sections with >= 10 Tregs: {int(r.n_sections_ge10_tregs)}); excluded from BH'))))

    def add(block, feats, RES, JNT, label, model):
        for c in feats:
            v = RES[c]; j = JNT.get(c) if JNT is not None else None; note = []
            if label == 'type':
                note.append(f'present near {PREV[c]:.3g}% of Tregs')
                if c == 'Absoprtive Epi': note.append("label spelled as in the annotation; a separate label from 'Absorptive Epi' (both co-occur in every section)")
            if c == 'lig_CXCL12_exp': note.append('CXCL12 over expected stromal producers (fibroblasts, myofibroblasts, pericytes, endothelial; >= 3 cells), the definition behind the Fig 4F regions')
            elif label == 'lig': note.append('mean within-section z of the ligand over all cells within 50 um (>= 5 cells)')
            if j is not None: note.append(f'joint 95% CI {j["lo"]:.3g} to {j["hi"]:.3g}')
            elif c == 'lig_CXCL12_exp': note.append('not in the joint model (r = 0.80 with lig_CXCL12)')
            rows.append(dict(block=block, feature=c, model=model, n_tregs=v['n'], effect_sd=v['coef'], ci_low=v['lo'], ci_high=v['hi'], p=v['p'], q=v['q'],
                             sections_same_sign=f'{v["agree"]}/{v["nsec"]}', clustered_p=v['cl_p'],
                             joint_effect=j['coef'] if j is not None else np.nan, joint_q=j['q'] if j is not None else np.nan, note='; '.join(note)))

    add('local cell type', CTS, RT, JT, 'type', 'LMM: score ~ type_z + cz + (1|section); type_z = fraction within 50 um, z within section')
    for c in [c for c in CT_ALL if c not in CTS]:
        rows.append(dict(block='local cell type', feature=c, model='not tested', n_tregs=N_ALL, note=f'present near {PREV[c]:.3g}% of Tregs (< {MIN_PREV:g}%)'))
    add('local ligand availability', LIGS_ALL, RL, JL, 'lig', 'LMM: score ~ ligand_z + cz + (1|section); ligand_z = local ligand, z within section')
    rows.append(dict(block='APC x CXCL12', feature='CXCL12 region (High vs Low) x local APC quintile (z)', model='LMM: score ~ CXCL12_high * APC_quintile_z + cz + (1|section) (Fig 4F)',
                     n_tregs=ri['n'], effect_sd=ri['coef'], ci_low=ri['lo'], ci_high=ri['hi'], p=ri['p'], sections_same_sign=f'{ri["agree"]}/{ri["nsec"]}',
                     clustered_p=ri['cl_p'], note='Fig 4F in-panel interaction P (equal to c_interaction_apc_by_cxcl12_countadj_exp.csv); single test, no BH family; CXCL12 regions = lower and upper within-section tertiles of lig_CXCL12_exp'))
    rows.append(dict(block='APC x CXCL12', feature='Rise lowest -> highest APC quintile, CXCL12-High minus CXCL12-Low', model='LMM: score ~ CXCL12_high * top_quintile + cz + (1|section), Q1 + Q5 Tregs only (Fig 4F)',
                     n_tregs=rd['n'], effect_sd=rd['coef'], ci_low=rd['lo'], ci_high=rd['hi'], p=rd['p'], sections_same_sign=f'{rd["agree"]}/{rd["nsec"]}',
                     clustered_p=rd['cl_p'], note=(f'Fig 4F difference P; rise {rise["low"][0]:+.2f} +- {rise["low"][1]:.2f} SD in CXCL12-Low vs '
                                               f'{rise["high"][0]:+.2f} +- {rise["high"][1]:.2f} SD in CXCL12-High Regions; single test, no BH family')))
    O = pd.DataFrame(rows)[COLS]

    def g3(x):
        if isinstance(x, str) or x is None or (isinstance(x, float) and np.isnan(x)): return '' if not isinstance(x, str) else x
        return f'{x:.3g}'

    for c in ['effect_sd', 'ci_low', 'ci_high', 'p', 'q', 'clustered_p', 'joint_effect', 'joint_q']: O[c] = O[c].map(g3)
    O['n_tregs'] = O.n_tregs.map(lambda v: '' if pd.isna(v) else str(int(v)))
    with open(f'{UPD}/cd_spatial_results_all.csv', 'w') as fh: O.to_csv(fh, index=False)


# ======================================================================================================================
# Section interaction_tests: DC x OGN+RSPO3+ Fib, DC x B cell, and SELENOP+ Mac subtypes on the Treg Suppression Score
#   primary: OLS with section FE, section-clustered SE, t(G - 1), always + cz; reference: MixedLM (section random intercept) + cz
# ======================================================================================================================
if _run('interaction_tests'):
    T, COMP = load_tregs()
    SK = {'SELENOP+MMP9+PLA2G2D+ Mac': 's_mp', 'SELENOP+MMP9+SPP1+ Mac': 's_ms', 'SELENOP+SIGLEC1+ Mac': 's_sg', 'SELENOP+XIST- Mac': 's_xs'}
    F = selenop_subtypes(T).rename(columns=SK)
    PREV = {k: 100 * (F[k] > 0).mean() for k in F.columns}
    X = pd.DataFrame(index=T.index); X['section'] = T.section; X['act_z'] = T.act_z; X['cz'] = T.cz
    RAW = {'dc': 'DC', 'ogn': 'OGN+RSPO3+ Fib', 'bc': 'B cell', 'cd9': 'CD9+ M2 Mac', 'selmac': 'SELENOP+ Mac'}
    for k, c in RAW.items(): X[k + '_f'] = COMP[c].values
    for k in SK.values(): X[k + '_f'] = F[k].values
    for k in list(RAW) + list(SK.values()): X[k] = X.groupby('section')[k + '_f'].transform(zs)
    X = X.dropna()

    def fit(d, rhs, terms):
        o = smf.ols(f'act_z ~ {rhs} + C(section)', d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]))
        g = d.section.nunique(); tc = stats.t.ppf(0.975, g - 1)
        m = smf.mixedlm(f'act_z ~ {rhs}', d, groups=d.section).fit(reml=True, method='lbfgs')
        sec = {}
        for s_, gd in d.groupby('section'):
            if len(gd) < 10 or any(gd[t.split(':')[0]].std() == 0 or gd[t.split(':')[-1]].std() == 0 for t in terms): continue
            sec[s_] = smf.ols(f'act_z ~ {rhs}', gd).fit().params[terms]
        out = {}
        for t in terms:
            b, se = o.params[t], o.bse[t]; a = sum(int(np.sign(v[t]) == np.sign(b)) for v in sec.values())
            out[t] = dict(estimate=b, lo=b - tc * se, hi=b + tc * se, clustered_p=2 * stats.t.sf(abs(b / se), g - 1), mixedlm_p=m.pvalues[t],
                          mixedlm_est=m.params[t], sections_agree=f'{a}/{len(sec)}', n=len(d))
        return out

    rows = []

    def add(block, term, r, note=''):
        rows.append(dict(block=block, term=term, **{k: r[k] for k in ['estimate', 'lo', 'hi', 'clustered_p', 'mixedlm_p', 'sections_agree', 'n']}, note=note))

    # 1-2: DC x OGN+RSPO3+ Fib and DC x B cell; DC slope in the bottom vs top within-section tertile of the second type
    for blk, mk, lab in [('1 DC x OGN+RSPO3+ Fib', 'ogn', 'OGN+RSPO3+ Fib'), ('2 DC x B cell', 'bc', 'B cell')]:
        r = fit(X, f'dc * {mk} + cz', ['dc', mk, f'dc:{mk}'])
        add(blk, 'DC_z (at average local ' + lab + ')', r['dc'], 'act_z ~ DC_z * x_z + cz')
        add(blk, lab + '_z (at average local DC)', r[mk])
        add(blk, f'DC_z x {lab}_z interaction', r[f'dc:{mk}'], 'negative = DC effect weaker where ' + lab + ' is locally enriched')
        X['t'] = X.groupby('section')[mk + '_f'].transform(tert)
        for tt in ['bottom', 'top']:
            d = X[X.t == tt]; z0 = 100 * (d[mk + '_f'] == 0).mean()
            add(blk, f'DC_z slope, {tt} tertile of local {lab}', fit(d, 'dc + cz', ['dc'])['dc'],
                f'act_z ~ DC_z + cz in Tregs of the {tt} within-section tertile of local {lab} fraction; {z0:.3g}% of these Tregs have 0 {lab} within 50 um')
        d = X[X.t.isin(['bottom', 'top'])].copy(); d['top'] = (d.t == 'top').astype(float); d['dctop'] = d.dc * d.top
        add(blk, f'DC_z slope difference, top minus bottom tertile of {lab}', fit(d, 'dc + top + dctop + cz', ['dctop'])['dctop'],
            'act_z ~ DC_z * top_tertile + cz, bottom + top tertile Tregs')
    X = X.drop(columns='t')
    # 3: SELENOP+ Mac subtypes, separate and joint (with DC_z and CD9+ M2 Mac_z)
    LAB = {v: k for k, v in SK.items()}; blk = '3 SELENOP+ Mac subtypes'
    add(blk, 'SELENOP+ Mac (merged)_z, separate', fit(X, 'selmac + cz', ['selmac'])['selmac'], f'present near {100 * (X.selmac_f > 0).mean():.3g}% of Tregs; act_z ~ x_z + cz')
    for k in SK.values():
        add(blk, f'{LAB[k]}_z, separate', fit(X, f'{k} + cz', [k])[k], f'present near {PREV[k]:.3g}% of Tregs; act_z ~ x_z + cz')
    JK = list(SK.values()) + ['dc', 'cd9']
    Rz = X[JK].corr().values; mr = np.abs(Rz[np.triu_indices(len(JK), 1)]).max()
    rj = fit(X, ' + '.join(JK) + ' + cz', JK)
    for k in JK:
        add(blk, f'{LAB.get(k, RAW.get(k))}_z, joint', rj[k], f'act_z ~ 4 SELENOP+ Mac subtype_z + DC_z + CD9+ M2 Mac_z + cz; max |r| between the 6 predictors = {mr:.2f}')
    O = pd.DataFrame(rows)

    def g3(x): return '' if x is None or (isinstance(x, float) and np.isnan(x)) else f'{x:.3g}'

    O['ci'] = [f'{g3(a)} to {g3(b)}' if pd.notna(a) else '' for a, b in zip(O.lo, O.hi)]
    for c in ['estimate', 'clustered_p', 'mixedlm_p']: O[c] = O[c].map(g3)
    O['n'] = O.n.map(lambda v: '' if pd.isna(v) else str(int(v))); O['sections_agree'] = O.sections_agree.fillna('')
    with open(f'{UPD}/cd_interaction_selenop_tests.csv', 'w') as fh:
        O[['block', 'term', 'estimate', 'ci', 'clustered_p', 'mixedlm_p', 'sections_agree', 'n', 'note']].to_csv(fh, index=False)


# ======================================================================================================================
# Section F5D_dc_by_ogn_tertile: local DC amount and DC slope by OGN+RSPO3+ tertile (DC % row of Fig 5D right, plus sensitivity analyses)
#   1 DC descriptives per tertile (pooled, per section); 2 tertile slopes unweighted, DC-matched (IPW to the section's pooled DC-decile
#   distribution) and on common support; Low - High slope difference from one fully interacted model; 3 non-linearity (dc^2 terms)
# ======================================================================================================================
if _run('F5D_dc_by_ogn_tertile'):
    d = apc_frame(); n_t, n_s = len(d), d.section.nunique()
    IX, _ = dc_by_ogn(d)
    d['dcp'] = 100 * d.dc_f   # DC fraction in % of all cells within 50 um
    ROWS = []

    def wq(x, w, q):
        o = np.argsort(x); x, w = x[o], w[o]; c = np.cumsum(w) / w.sum()
        return x[np.searchsorted(c, q)]

    def desc(dd, w=None):
        """DC descriptives (weighted if w); sd_within = SD of dc after removing the section means of this subset."""
        w = np.ones(len(dd)) if w is None else np.asarray(w, float); r = {}
        for v, lab in [('dcp', 'dcpct'), ('dc', 'dcz')]:
            x = dd[v].values; m = np.average(x, weights=w)
            r[f'{lab}_mean'] = m; r[f'{lab}_median'] = wq(x, w, 0.5); r[f'{lab}_sd'] = np.sqrt(np.average((x - m) ** 2, weights=w))
            r[f'{lab}_q1'] = wq(x, w, 0.25); r[f'{lab}_q3'] = wq(x, w, 0.75); r[f'{lab}_iqr'] = r[f'{lab}_q3'] - r[f'{lab}_q1']
        sm = pd.Series(w * dd.dc.values, index=dd.index).groupby(dd.section).transform('sum') / pd.Series(w, index=dd.index).groupby(dd.section).transform('sum')
        r['dcz_sd_within'] = np.sqrt(np.average((dd.dc.values - sm.values) ** 2, weights=w))
        r['pct_zero_dc'] = 100 * np.average(dd.dc_f.values == 0, weights=w); r['ogn_pct_mean'] = 100 * np.average(dd.ogn_f.values, weights=w)
        return r

    # 1. descriptives
    for tt in TT:
        ROWS.append(dict(block='1_desc_pooled', section='all', tertile=TL[tt], n=int((d.t == tt).sum()), **desc(d[d.t == tt])))
    for s_, gs in d.groupby('section'):
        kw = stats.kruskal(*[gs.loc[gs.t == tt, 'dcp'] for tt in TT]).pvalue
        for tt in TT:
            ROWS.append(dict(block='1_desc_section', section=s_, tertile=TL[tt], n=int((gs.t == tt).sum()), **desc(gs[gs.t == tt]),
                             note=f'Kruskal-Wallis DC fraction across tertiles P = {kw:.3g}' if tt == 'bottom' else ''))
    o = smf.ols('dc ~ C(t, Treatment("bottom")) + C(section)', d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]))
    for tt in ['middle', 'top']:
        ROWS.append(dict(block='1_desc_test', section='all', tertile=f'{TL[tt]} - Low', term='mean dc (z) difference, section FE, clustered', n=n_t,
                         **summ(o, f'C(t, Treatment("bottom"))[T.{tt}]', n_s)))

    # 2. tertile slopes
    def slopes(dd, label, w=None, note=''):
        for tt in TT:
            m = (dd.t == tt).values; dt = dd[m]
            r = summ(clfit(dt, 'dc + cz', None if w is None else w[m]), 'dc', dt.section.nunique())
            r['agree'] = agree(dt, 'dc + cz', 'dc', r['estimate']) if w is None else ''
            ROWS.append(dict(block=label, section='all', tertile=TL[tt], term='dc', n=int(m.sum()), **r, **desc(dt, None if w is None else w[m]), note=note))
        # Low - High slope difference: fully interacted model (section x tertile FE, dc and cz slopes by tertile)
        dd = dd.assign(st=dd.section.astype(str) + '|' + dd.t.astype(str))
        f = 'act_z ~ 0 + C(st) + dc:C(t, Treatment("bottom")) + cz:C(t, Treatment("bottom"))'
        oo = (smf.ols(f, dd) if w is None else smf.wls(f, dd, weights=w)).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(dd.section)[0]))
        nm = {tt: f'dc:C(t, Treatment("bottom"))[{tt}]' for tt in TT}
        for a, b in [('bottom', 'top'), ('bottom', 'middle')]:
            L = np.zeros(len(oo.params)); L[list(oo.params.index).index(nm[a])] = 1; L[list(oo.params.index).index(nm[b])] = -1
            est = L @ oo.params.values; se = np.sqrt(L @ oo.cov_params().values @ L); tq = stats.t.ppf(0.975, n_s - 1)
            ROWS.append(dict(block=label, section='all', tertile=f'{TL[a]} - {TL[b]}', term='dc slope difference', n=len(dd),
                             estimate=est, lo=est - tq * se, hi=est + tq * se, p=2 * stats.t.sf(abs(est / se), n_s - 1), note=note))

    slopes(d, '2_slopes_unweighted')
    d['bin'] = d.groupby('section').dc_f.transform(lambda v: pd.qcut(v, 10, labels=False, duplicates='drop'))
    ps = d.groupby(['section', 'bin']).size() / d.groupby('section').size()
    pst = d.groupby(['section', 't', 'bin'], observed=True).size() / d.groupby(['section', 't'], observed=True).size()
    W = ps.reindex(pd.MultiIndex.from_arrays([d.section, d.bin])).values / pst.reindex(pd.MultiIndex.from_arrays([d.section, d.t, d.bin])).values
    full = pd.MultiIndex.from_product([sorted(d.section.unique()), TT, sorted(d.bin.unique())], names=['section', 't', 'bin'])
    cnt_st = d.groupby(['section', 't', 'bin'], observed=True).size().reindex(full, fill_value=0)
    miss_mass = 100 * sum(ps[(s_, b)] for (s_, tt, b), c in cnt_st.items() if c == 0 and (s_, b) in ps.index) / (n_s * 3)   # pooled bin mass unmatched
    nbins = d.groupby('section').bin.nunique()
    slopes(d, '2A_slopes_dc_matched_ipw', w=W, note=(f'IPW within section to pooled DC-fraction decile distribution (tied edges dropped; {nbins.min()}-{nbins.max()} bins per section); '
                                                     f'weight range {W.min():.3g}-{W.max():.3g}; unmatched pooled bin mass {miss_mass:.3g}%'))
    for lab, qlo, qhi, txt in [('2B_slopes_common_support_minmax', 0, 1, 'overlap of tertile min-max DC-fraction ranges'),
                               ('2C_slopes_common_support_p5_p95', 0.05, 0.95, 'overlap of tertile 5th-95th percentile DC-fraction ranges')]:
        keep = pd.Series(False, index=d.index)
        for s_, gs in d.groupby('section'):
            lo = max(gs.loc[gs.t == tt, 'dc_f'].quantile(qlo) for tt in TT); hi = min(gs.loc[gs.t == tt, 'dc_f'].quantile(qhi) for tt in TT)
            keep.loc[gs.index] = (gs.dc_f >= lo) & (gs.dc_f <= hi)
        dd = d[keep].copy(); slopes(dd, lab, note=f'within section, {txt}; kept {len(dd):,}/{n_t:,} Tregs')

    # 3. non-linearity
    for lab, rhs in [('3_nonlinear_dc2', 'dc * ogn + I(dc**2) + cz'), ('3_nonlinear_dc2_x_ogn', 'dc * ogn + I(dc**2) + I(dc**2):ogn + cz')]:
        o = clfit(d, rhs)
        for tm in ['dc', 'ogn', 'dc:ogn', 'I(dc ** 2)'] + (['I(dc ** 2):ogn'] if 'I(dc**2):ogn' in rhs else []):
            r = summ(o, tm, n_s)
            ag = agree(d.assign(dc2=d.dc ** 2), rhs.replace('I(dc**2)', 'dc2'), tm.replace('I(dc ** 2)', 'dc2'), r['estimate'])
            ROWS.append(dict(block=lab, section='all', tertile='continuous', term=tm, n=n_t, **r, agree=ag, note=f'act_z ~ {rhs} + C(section)'))
    ROWS.append(dict(block='3_reference_linear', section='all', tertile='continuous', term='dc:ogn', n=n_t, **IX, note='act_z ~ dc * ogn + cz + C(section) (main model)'))
    cols = ['block', 'section', 'tertile', 'term', 'n', 'estimate', 'lo', 'hi', 'p', 'agree', 'ogn_pct_mean', 'dcpct_mean', 'dcpct_median', 'dcpct_sd',
            'dcpct_q1', 'dcpct_q3', 'dcpct_iqr', 'pct_zero_dc', 'dcz_mean', 'dcz_median', 'dcz_sd', 'dcz_q1', 'dcz_q3', 'dcz_iqr', 'dcz_sd_within', 'note']
    R_ = pd.DataFrame(ROWS).reindex(columns=cols)
    for c in cols:
        if R_[c].dtype.kind == 'f': R_[c] = R_[c].map(lambda x: x if not isinstance(x, (float, np.floating)) or not np.isfinite(x) else float(f'{x:.3g}'))
    with open(f'{FIG}/fig4d_dc_by_ogn_tertile_checks.csv', 'w') as fh: R_.to_csv(fh, index=False)


# ======================================================================================================================
# Section F5D_apc_types_dc_ogn: Fig 5D. Left, one regression of the score on the local fractions of the 7 APC types (BH over 7; dots = per-section
# OLS coefficients; strip = mean local %). Right, DC slope per OGN+RSPO3+ tertile (filled = P < 0.05) over per-section slopes; DC % row from
# fig4d_dc_by_ogn_tertile_checks.csv
# ======================================================================================================================
if _run('F5D_apc_types_dc_ogn'):
    d = apc_frame(); n_t, n_s = len(d), d.section.nunique()
    RHS = ' + '.join(K) + ' + cz'
    o = clfit(d, RHS); RES = {k: summ(o, k, n_s) for k in K}
    for k, qq in zip(K, multipletests([RES[k]['p'] for k in K], method='fdr_bh')[1]): RES[k]['q'] = qq
    SEC = pd.DataFrame({s_: smf.ols(f'act_z ~ {RHS}', g).fit().params[K] for s_, g in d.groupby('section')}).T
    ABUND = {k: 100 * d[k + '_f'].mean() for k in K}
    IX, TER = dc_by_ogn(d)
    # per-section DC slope in each tertile
    PAT = lambda sec: sec.split('_')[0]
    PS = pd.DataFrame([dict(section=s_, patient=PAT(s_), tertile=TL[tt], tertile_key=tt, slope=f_.params['dc'], se=f_.bse['dc'], n_tregs=len(gd))
                       for (s_, tt), gd in d.groupby(['section', 't'], observed=True) for f_ in [smf.ols('act_z ~ dc + cz', gd).fit()]])
    WIDE = PS.pivot(index='section', columns='tertile_key', values='slope')
    SECS = sorted(WIDE.index, key=lambda x: (PAT(x), x))
    SPAL = ['#1f77b4', '#ff7f0e', '#2ca02c', '#CC3311', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#DAA520', '#17255A']
    SCOL = dict(zip(SECS, SPAL)); PS['color'] = PS.section.map(SCOL)
    OUTCSV = PS[['section', 'patient', 'tertile', 'slope', 'se', 'n_tregs', 'color']].copy()
    OUTCSV['tertile'] = pd.Categorical(OUTCSV.tertile, ['Low', 'Mid', 'High']); OUTCSV = OUTCSV.sort_values(['patient', 'section', 'tertile'])
    for c in ['slope', 'se']: OUTCSV[c] = OUTCSV[c].map(lambda v: float(f'{v:.3g}'))
    with open(f'{FIG}/fig4d_dc_slope_by_ogn_tertile_per_section.csv', 'w') as fh: OUTCSV.to_csv(fh, index=False)
    # plot
    YLO, YHI = -0.3, 0.4; YTK_R = [-0.2, 0, 0.2]
    FT, FL, FSm, FKEY = 15, 15, 13, 10; DARK = '#555555'
    plt.rcParams.update({'font.size': FT, 'font.family': 'DejaVu Sans', 'font.weight': 'normal', 'axes.titleweight': 'normal',
                         'axes.labelweight': 'normal', 'axes.linewidth': 1.0})
    fig = plt.figure(figsize=(7.6, 4.75))   # columns: abundance strip, forest, gap for the right y label, DC slopes
    gs = fig.add_gridspec(1, 4, width_ratios=[0.7, 1.52, 1.06, 2.0], wspace=0.076, left=0.253, right=0.987, bottom=0.33, top=0.903)
    axa, axf, axr = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[3])

    def est(ax, ci, pt, filled):   # pooled estimate and 95% CI, with a white halo over the section dots
        ax.plot(*ci, color='white', lw=3.4, zorder=1.9, solid_capstyle='butt'); ax.plot(*ci, color='black', lw=1.8, zorder=2, solid_capstyle='butt')
        ax.plot(*pt, 'o', ms=8.5 + 2.6, mfc='white', mec='white', mew=0, zorder=1.95)
        ax.plot(*pt, 'o', ms=8.5, mec='black', mew=1.4, mfc='black' if filled else 'white', zorder=3)

    WRAP = {'s_mp': 'SELENOP⁺MMP9⁺\nPLA2G2D⁺ Mac', 's_ms': 'SELENOP⁺MMP9⁺\nSPP1⁺ Mac', 's_sgxs': 'SELENOP⁺ Mac\n(SIGLEC1⁺/XIST⁻)'}
    LBL = {k: WRAP.get(k, APC_TYPES[k][2]) for k in K}; NL = {k: LBL[k].count('\n') + 1 for k in K}
    ypos = dict(zip(K, -np.cumsum([0] + [0.115 * (NL[a] + NL[b]) + 0.045 for a, b in zip(K[:-1], K[1:])])))   # 2-line rows get more room
    YL = (ypos[K[-1]] - 0.2, 0.2)
    amax = max(ABUND.values())   # abundance strip
    for k in K:
        axa.barh(ypos[k], ABUND[k], height=0.11, color='#9E9E9E', lw=0)
        axa.text(ABUND[k] + 0.06 * amax, ypos[k], f'{ABUND[k]:.1f}', ha='left', va='center', fontsize=FSm, color=DARK)
    axa.set_xlim(0, amax * 2.9); axa.set_ylim(YL)
    axa.set_yticks([ypos[k] for k in K]); axa.set_yticklabels([LBL[k] for k in K], fontsize=FT, linespacing=0.95)
    axa.set_xticks([]); axa.spines[['top', 'right', 'bottom']].set_visible(False)
    axa.text(0, YL[1] + 0.02, 'Mean %\nwithin 50 µm', ha='left', va='bottom', fontsize=FSm, color=DARK, linespacing=1.0)
    rng = np.random.default_rng(0)   # forest
    axf.axvline(0, color='#777777', lw=1.1, ls='--', zorder=1)
    for k in K:
        v = RES[k]; y = ypos[k]
        axf.scatter(SEC[k].values, y + rng.uniform(-0.075, 0.075, len(SEC)), s=13, color=[SCOL[x] for x in SEC.index], alpha=0.5, lw=0, zorder=1.5)
        est(axf, ([v['lo'], v['hi']], [y, y]), (v['estimate'], y), v['q'] < 0.05)
    lo_all = min(SEC.values.min(), min(RES[k]['lo'] for k in K)); hi_all = max(SEC.values.max(), max(RES[k]['hi'] for k in K))
    pad = 0.05 * (hi_all - lo_all); axf.set_xlim(lo_all - pad, hi_all + pad); axf.set_ylim(YL)
    axf.set_xticks([-0.2, 0, 0.2]); axf.set_xticklabels(['−0.2', '0', '0.2'])
    axf.set_yticks([]); axf.spines[['top', 'right', 'left']].set_visible(False); axf.tick_params(labelsize=FT)
    axf.set_xlabel('Effect on Treg\nSuppression Score\n(SD per SD of\nlocal fraction)', fontsize=FL, linespacing=1.0)
    axr.axhline(0, color='#777777', lw=1.1, ls='--', zorder=1)   # right: DC slope by OGN+RSPO3+ tertile
    OFF = dict(zip(SECS, np.linspace(-0.13, 0.13, len(SECS))))
    for s_ in SECS:
        xs = [i + OFF[s_] for i in range(3)]; ys = [WIDE.loc[s_, tt] for tt in TT]
        axr.plot(xs, ys, '-', color=SCOL[s_], lw=1.0, alpha=0.4, zorder=1.5)
        axr.scatter(xs, ys, s=14, color=SCOL[s_], lw=0, alpha=0.5, zorder=1.6)
    for i, tt in enumerate(TT):
        v = TER[tt]; est(axr, ([i, i], [v['lo'], v['hi']]), (i, v['estimate']), v['p'] < 0.05)
    axr.set_xlim(-0.4, 2.4); axr.set_xticks(range(3)); axr.set_xticklabels([TL[tt] for tt in TT], fontsize=FT)
    CHK = pd.read_csv(f'{FIG}/fig4d_dc_by_ogn_tertile_checks.csv'); CHK = CHK[CHK.block == '1_desc_pooled'].set_index('tertile')
    for lab, vals, dy in [('n =', [f'{TER[tt]["n"]:,}' for tt in TT], 27),   # n and DC % rows under the tick labels
                          ('DC %', [f'{float(CHK.loc[TL[tt], "dcpct_mean"]):.1f}' for tt in TT], 27 + FSm * 1.25 + 2)]:
        for x, txt, ha, dx in [(-0.4, lab, 'right', -4)] + [(i, t, 'center', 0) for i, t in enumerate(vals)]:
            axr.annotate(txt, (x, 0), xycoords=('data', 'axes fraction'), xytext=(dx, -dy), textcoords='offset points', ha=ha, va='top', fontsize=FSm, color=DARK)
    axr.set_ylim(YLO, YHI); axr.set_yticks(YTK_R); axr.set_yticklabels([('−' if t < 0 else '') + f'{abs(t):g}' for t in YTK_R])
    axr.spines[['top', 'right']].set_visible(False); axr.tick_params(axis='y', labelsize=FT)
    axr.set_ylabel('DC effect on Treg\nSuppression Score\n(SD per SD)', fontsize=FL, linespacing=1.0)
    axr.set_xlabel('Local OGN⁺RSPO3⁺\nfibroblasts\n(within-section tertile)', fontsize=FL, linespacing=1.0, loc='right', labelpad=42, multialignment='center')
    axr.annotate(f'DC × OGN⁺RSPO3⁺ interaction\nP = {IX["p"]:.3g}; {IX["agree"]} sections', (1.0, 1.0), xycoords='axes fraction',
                 xytext=(0, 3), textcoords='offset points', ha='right', va='bottom', fontsize=FSm, color='black', linespacing=1.05)
    axr.text(0.02, 0.995, 'Thin lines: sections\n(one color each)\nBold: pooled ± 95% CI', transform=axr.transAxes, ha='left', va='top',
             fontsize=FKEY, color=DARK, linespacing=1.05, zorder=4)
    with open(f'{FIG}/fig4_apc_types_dc_ogn_combined.png', 'wb') as fh: fig.savefig(fh, format='png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)


# ======================================================================================================================
# Section S9D_celltype_forest: every local cell type vs the score, from cd_spatial_results_all.csv (nothing refitted); SELENOP+ Fib not shown;
# 'Absoprtive Epi' (annotation spelling) = BEST4+ cells; filled = q < 0.05, clustered P < 0.05 and >= 9/11 sections same sign; DC row in red
# ======================================================================================================================
if _run('S9D_celltype_forest'):
    PNG = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_all_celltypes_treg_score_forest.png'
    DISP = {'OGN+RSPO3+ Fib': 'OGN⁺RSPO3⁺ Fib', 'ADAMDEC1+ Fib': 'ADAMDEC1⁺ Fib', 'PLA2G2A+ ECM Fib': 'PLA2G2A⁺ ECM Fib',
            'CD9+ M2 Mac': 'CD9⁺ M2 Mac', 'Absoprtive Epi': 'BEST4⁺ Epi', 'Gd T': 'γδ T', 'SELENOP+ Mac': 'SELENOP⁺ Mac'}
    T = pd.read_csv(f'{UPD}/cd_spatial_results_all.csv')
    T = T[(T.block == 'local cell type') & (T.feature != 'SELENOP+ Fib')].copy()
    T['n_same'] = T.sections_same_sign.str.split('/').str[0].astype(int)
    T['filled'] = (T.q < 0.05) & (T.clustered_p < 0.05) & (T.n_same >= 9)
    T['label'] = T.feature.map(lambda s: DISP.get(s, s))
    T = T.sort_values('effect_sd', ascending=True).reset_index(drop=True)   # most negative at bottom
    y = np.arange(len(T))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13, 'axes.linewidth': 1.0})
    FS_LAB, FS_SEC, FS_KEY = 15, 13, 13
    fig, ax = plt.subplots(figsize=(6.0, 10.0)); fig.subplots_adjust(left=0.335, right=0.855, bottom=0.171, top=0.95)
    ax.axvline(0, color='0.45', ls='--', lw=1.1, zorder=0)
    RED = '#c0392b'; col = np.where(T.feature.values == 'DC', RED, 'black')
    ax.hlines(y, T.ci_low, T.ci_high, colors=col, lw=1.6, zorder=1)
    r = T.filled.values
    ax.scatter(T.effect_sd[r], y[r], s=70, facecolor=col[r], edgecolor=col[r], lw=1.3, zorder=3)
    ax.scatter(T.effect_sd[~r], y[~r], s=70, facecolor='white', edgecolor='black', lw=1.3, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(T.label, fontsize=FS_LAB)
    for t in ax.get_yticklabels():
        if t.get_text() == 'DC': t.set_color(RED)
    ax.set_ylim(-0.7, len(T) - 0.3); ax.set_xlim(-0.125, 0.135)
    ax.set_xticks([-0.1, 0, 0.1]); ax.set_xticklabels(['−0.1', '0', '0.1'], fontsize=FS_LAB); ax.set_xticks([-0.05, 0.05], minor=True)
    ax.tick_params(axis='y', length=0, pad=3)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    ax.set_xlabel('Effect on Treg Suppression Score\n(SD per SD of local fraction within 50 µm)', fontsize=FS_LAB, labelpad=4)
    tr = ax.get_yaxis_transform()   # sections with the same sign, right of the axes
    for yi, t in zip(y, T.sections_same_sign): ax.text(1.02, yi, t, transform=tr, ha='left', va='center', fontsize=FS_SEC, color='0.25', clip_on=False)
    ax.text(1.02, len(T) - 0.1, 'Same\nsign', transform=tr, ha='left', va='bottom', fontsize=FS_SEC, color='0.25', clip_on=False)
    h = [Line2D([], [], marker='o', ls='', mfc='black', mec='black', ms=8.5, label='q < 0.05\nSection-clustered P < 0.05\n≥ 9/11 sections same sign'),
         Line2D([], [], marker='o', ls='', mfc='white', mec='black', ms=8.5, label='Otherwise')]
    fig.legend(handles=h, loc='lower center', bbox_to_anchor=(0.5, 0.004), ncol=2, frameon=False, fontsize=FS_KEY,
               handletextpad=0.2, columnspacing=1.0, borderpad=0.1, labelspacing=0.4)
    fig.savefig(PNG, dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
