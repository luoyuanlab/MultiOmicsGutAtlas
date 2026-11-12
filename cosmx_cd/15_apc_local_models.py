#!/usr/bin/env python
"""CD CosMx 6K resections: Treg Suppression Score vs local cell types, with APC and DC models.

Panels made
  Fig 5D   joint model of local APC-type fractions (left); DC effect by within-section tertile of local OGN+RSPO3+ fibroblasts, with per-section
           slopes and the DC % row (right)
  Supp 9D  Treg Suppression Score per SD of each cell type's local fraction (separate LMMs)
  Numbers in the Results: all-feature summary table (cd_spatial_results_all.csv) and interaction tests (cd_interaction_selenop_tests.csv)

Inputs
  /share/fsmresfiles/UC/AtoMx/CD_surgical_resection/combined/ (all_samples.parquet, raw_counts_strict_tregs.csv)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/
    B_resection_analyses/ligand_receptor_state/ (treg_local_*.csv.gz, c_interaction_apc_by_cxcl12_countadj_exp.csv; 11_neighbourhoods_local_tables.py)
    A_fig4_panels/fig4_cn_treg_state_forest.csv (12_neighbourhood_treg_score.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/updated_figures/ (cd_spatial_results_all.csv + _README.txt, cd_interaction_selenop_tests.csv)
  /share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/ (Fig 5D panel, legend and CSVs)
  /share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/ (Supp 9D panel and legend)

Run order of sections
  spatial_results_table -> interaction_tests -> F5D_dc_by_ogn_tertile -> F5D_apc_types_dc_ogn -> S9D_celltype_forest
  F5D_apc_types_dc_ogn reads the table of F5D_dc_by_ogn_tertile; S9D_celltype_forest reads cd_spatial_results_all.csv of spatial_results_table
  (also read by section S9F_apc_cxcl12_quintiles of 16_cxcl12_regions.py).

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its upstream
  outputs exist. Command-line arguments and environment variables used by a section are listed in its own docstring.
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
# Section spatial_results_table: all-feature summary table of the CD resection spatial results (neighbourhoods, local cell types, local ligands, APC x CXCL12)
# Panels: Supp 9D values, Results text (cd_spatial_results_all.csv)
# ======================================================================================================================
if _run('spatial_results_table'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/updated_figures/cd_spatial_results_all.py'   # working folder of this step (relative paths resolve from here)
    """All-feature summary table of the CD resection spatial results: every feature of each block is tested the same way and reported.
CD resection CosMx 6K, 11 sections from 4 patients, 6,554 strict Tregs. Read-only inputs.
Blocks (rows in neutral order within each block: CN number or alphabetical):
  neighbourhood      fig4_cn_treg_state_forest.csv (all 10 CNs, in-vs-out LMM of the CN forest plot of 12_neighbourhood_treg_score.py; copied as is).
  local cell type    every cell_type_general column of treg_local_composition.csv.gz (per-Treg fraction among all cells within 50 um), except Treg;
                     tested if present (fraction > 0) near >= 5% of Tregs. type_z = fraction z-scored within section.
  local ligand       every lig_* column of treg_local_ligand_receptor.csv.gz (mean within-section z of the ligand over all cells within 50 um),
                     plus lig_CXCL12_exp (CXCL12 over expected stromal producers, the definition behind the CXCL12-High Regions), re-z-scored
                     within section (tables of 11_neighbourhoods_local_tables.py).
  APC x CXCL12       CXCL12-High vs -Low Region x APC quintile interaction P (refitted and asserted equal to c_interaction_apc_by_cxcl12_countadj_exp.csv)
                     and the difference-in-differences P of the rise from the lowest to the highest APC quintile (16_cxcl12_regions.py).
Models (cell type and ligand blocks), cz = within-section z of log1p(library size):
  separate  act_z ~ x_z + cz + (1|section)      MixedLM REML lbfgs
  joint     act_z ~ all x_z of the block + cz + (1|section), fitted only if max |Pearson r| between the x_z < 0.6
  robustness OLS act_z ~ x_z + cz with section-clustered SE, t(G - 1) inference (G = 11)
  per-section OLS act_z ~ x_z + cz -> sections whose coefficient has the sign of the pooled LMM estimate.
BH q: separate across all tested features of the block; joint across all features of the joint model.
Outputs: cd_spatial_results_all.csv and cd_spatial_results_all_README.txt in CDALL_OUT (default /share/fsmresfiles/UC/revision_jenny/updated_figures)."""
    import os, sys, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    import numpy as np, pandas as pd
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
    from scipy import stats
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    LR = f'{P4}/B_resection_analyses/ligand_receptor_state'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
    FOREST = f'{P4}/A_fig4_panels/fig4_cn_treg_state_forest.csv'
    OUT = os.environ.get('CDALL_OUT', '/share/fsmresfiles/UC/revision_jenny/updated_figures')
    CSV, README = f'{OUT}/cd_spatial_results_all.csv', f'{OUT}/cd_spatial_results_all_README.txt'
    MIN_PREV = 5.0; MAX_R = 0.6
    def zs(v):
        s = v.std()
        return (v - v.mean()) / s if s > 0 else v * 0.0
    def close(a, b, rtol=5e-3, atol=1e-4): return bool(np.isclose(a, b, rtol=rtol, atol=atol))
    # ---------------- data ----------------
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    T = T.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left'); assert T.n_exp_producers_50um.notna().all()
    C = pd.read_csv(f'{LR}/treg_local_composition.csv.gz', index_col=0)
    assert set(C.index) == set(T.index); C = C.loc[T.index]
    assert (C.nb_n.values == T.nb_n.values).all() and (C.section.values == T.section.values).all()
    rs = C.drop(columns=['section', 'nb_n']).sum(1)
    assert np.allclose(rs[C.nb_n > 0], 1, atol=5e-3), 'composition is not stored as fractions'
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    assert len(cnt) == len(T)
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()          # library size (Treg order of the raw-count table)
    T['cz'] = T.groupby('section').cz.transform(zs)
    C = C.loc[T.index]
    N_ALL, N_SEC = len(T), T.section.nunique()
    assert N_ALL == 6554 and N_SEC == 11, (N_ALL, N_SEC)
    TYPES_ALL = sorted([c for c in C.columns if c not in ('section', 'nb_n', 'Treg')], key=str.lower)
    PREV = {c: 100 * (C[c] > 0).mean() for c in TYPES_ALL}
    TYPES = [c for c in TYPES_ALL if PREV[c] >= MIN_PREV]
    LIGS = sorted([c for c in T.columns if c.startswith('lig_') and c != 'lig_CXCL12_exp'], key=str.lower)
    assert LIGS == ['lig_CCL22', 'lig_CXCL12', 'lig_CXCL14', 'lig_IL6', 'lig_TNF'], LIGS
    LIGS_ALL = LIGS + ['lig_CXCL12_exp']
    X = pd.DataFrame(index=T.index); X['section'] = T.section; X['act_z'] = T.act_z; X['cz'] = T.cz
    key = {}
    for i, c in enumerate(TYPES_ALL): key[c] = f'ct{i}'; X[f'ct{i}'] = C[c].values
    for i, c in enumerate(LIGS_ALL): key[c] = f'lg{i}'; X[f'lg{i}'] = T[c].values
    for c in TYPES_ALL + LIGS_ALL:                                    # within-section z (NaN kept for Tregs without a value)
        X[key[c]] = X.groupby('section')[key[c]].transform(zs)
    # ---------------- model helpers ----------------
    def lmm(d, f, terms):
        m = smf.mixedlm(f, d, groups=d.section).fit(reml=True, method='lbfgs')
        ci = m.conf_int()
        return {t: dict(coef=m.params[t], lo=ci.loc[t, 0], hi=ci.loc[t, 1], p=m.pvalues[t], n=len(d)) for t in terms}
    def clustered_p(d, f, t):
        o = smf.ols(f, d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]), use_t=True)
        G = d.section.nunique(); tv = o.params[t] / o.bse[t]
        return 2 * stats.t.sf(abs(tv), G - 1)                         # t(G - 1) inference
    def sec_agree(d, f, t, sign):
        a, n = 0, 0
        for _, g in d.groupby('section'):
            if g[t].std() == 0 or len(g) < 5: continue                   # coefficient not estimable in this section
            b = smf.ols(f, g).fit().params[t]; n += 1; a += int(np.sign(b) == sign)
        return a, n
    def screen(feats):
        out = {}
        for c in feats:
            k = key[c]; d = X[['section', 'act_z', 'cz', k]].dropna()
            f = f'act_z ~ {k} + cz'
            r = lmm(d, f, [k])[k]
            r['agree'], r['nsec'] = sec_agree(d, f, k, np.sign(r['coef']))
            r['cl_p'] = clustered_p(d, f, k)
            out[c] = r
        q = multipletests([out[c]['p'] for c in feats], method='fdr_bh')[1]
        for c, qq in zip(feats, q): out[c]['q'] = qq
        return out
    def joint(feats, label):
        ks = [key[c] for c in feats]; d = X[['section', 'act_z', 'cz'] + ks].dropna()
        Rz = d[ks].corr().values; mr = np.abs(Rz[np.triu_indices(len(ks), 1)]).max()
        ia, ib = np.unravel_index(np.argmax(np.abs(Rz - np.eye(len(ks)))), Rz.shape)
        pair = (feats[ia], feats[ib], Rz[ia, ib])
        print(f'{label}: max |r| between features = {mr:.3f} ({pair[0]} vs {pair[1]})', flush=True)
        if mr >= MAX_R: return None, mr, pair, len(d)
        r = lmm(d, 'act_z ~ ' + ' + '.join(ks) + ' + cz', ks)
        q = multipletests([r[k]['p'] for k in ks], method='fdr_bh')[1]
        return {c: dict(r[key[c]], q=qq) for c, qq in zip(feats, q)}, mr, pair, len(d)
    # ---------------- cell types and ligands ----------------
    RT = screen(TYPES); JT, mrT, prT, njT = joint(TYPES, 'cell types')
    RL = screen(LIGS_ALL); JL, mrL, prL, njL = joint(LIGS, 'ligands (lig_* columns of the table)')   # lig_CXCL12_exp is a redefinition of lig_CXCL12 (r ~ 0.8): separate model only
    _d = X[['section', 'act_z', 'cz', key['DC']]]; _m = smf.mixedlm(f"act_z ~ {key['DC']} + cz", _d, groups=_d.section).fit(reml=True, method='lbfgs')
    RE_VAR = float(_m.cov_re.iloc[0, 0]); print(f'section random-intercept variance (DC model) = {RE_VAR:.2e}', flush=True)
    # ---------------- APC x CXCL12: interaction (asserted equal to the table of 11_neighbourhoods_local_tables.py) and difference-in-differences ----------------
    IX = pd.read_csv(f'{LR}/c_interaction_apc_by_cxcl12_countadj_exp.csv').set_index('term')
    A = T[['section', 'act_z', 'cz', 'apc_frac', 'lig_CXCL12_exp']].copy()
    A['apc_q'] = A.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    A['t_exp'] = A.groupby('section').lig_CXCL12_exp.transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))
    di = A[A.t_exp.isin(['low', 'high']) & A.apc_q.notna()].copy(); di['hi'] = (di.t_exp == 'high').astype(float)
    di['qz'] = (di.apc_q - di.apc_q.mean()) / di.apc_q.std(); di = di[['section', 'act_z', 'hi', 'qz', 'cz']].dropna()
    FI = 'act_z ~ hi * qz + cz'; ri = lmm(di, FI, ['hi:qz'])['hi:qz']
    assert close(ri['coef'], IX.coef['hi:apcz']) and close(ri['p'], IX.p['hi:apcz']) and ri['n'] == int(IX.n['hi:apcz']), ri
    ri['agree'], ri['nsec'] = sec_agree(di.assign(**{'hixqz': di.hi * di.qz}), 'act_z ~ hi + qz + hixqz + cz', 'hixqz', np.sign(ri['coef']))
    ri['cl_p'] = clustered_p(di.assign(hixqz=di.hi * di.qz), 'act_z ~ hi + qz + hixqz + cz', 'hixqz')
    dd = A[A.t_exp.isin(['low', 'high']) & A.apc_q.isin([0, 4])].copy()
    dd['hi'] = (dd.t_exp == 'high').astype(float); dd['top'] = (dd.apc_q == 4).astype(float); dd = dd[['section', 'act_z', 'hi', 'top', 'cz']].dropna()
    FD = 'act_z ~ hi * top + cz'; rd = lmm(dd, FD, ['hi:top'])['hi:top']
    dd['hixtop'] = dd.hi * dd.top
    rd['agree'], rd['nsec'] = sec_agree(dd, 'act_z ~ hi + top + hixtop + cz', 'hixtop', np.sign(rd['coef']))
    rd['cl_p'] = clustered_p(dd, 'act_z ~ hi + top + hixtop + cz', 'hixtop')
    rise = {}
    for g in ['low', 'high']:
        s = A[(A.t_exp == g)]; q1, q5 = s[s.apc_q == 0].act_z, s[s.apc_q == 4].act_z
        rise[g] = (q5.mean() - q1.mean(), np.sqrt(q1.sem() ** 2 + q5.sem() ** 2))
    print('APC x CXCL12: interaction P = {:.3g}, difference P = {:.3g}'.format(ri['p'], rd['p']), flush=True)
    # ---------------- assemble ----------------
    COLS = ['block', 'feature', 'model', 'n_tregs', 'effect_sd', 'ci_low', 'ci_high', 'p', 'q', 'sections_same_sign', 'clustered_p',
            'joint_effect', 'joint_q', 'note']
    rows = []
    FO = pd.read_csv(FOREST).sort_values('CN')
    assert len(FO) == 10 and list(FO.CN) == list(range(1, 11))
    for r in FO.itertuples():
        tested = bool(r.tested)
        rows.append(dict(block='neighbourhood', feature=r.CN_name, model='LMM: score ~ in_CN + (1|section) (Fig 4 forest)',
                         n_tregs=N_ALL, effect_sd=r.treg_suppression_score_z_in_minus_out, ci_low=r.ci_lo, ci_high=r.ci_hi, p=r.mixed_p,
                         q=r.mixed_q if tested else np.nan, sections_same_sign=f'{int(r.sections_agreeing)}/{int(r.sections)}', clustered_p=np.nan,
                         joint_effect=np.nan, joint_q=np.nan,
                         note=(f'{int(r.n_tregs_in_CN):,} Tregs in CN; effect = in-minus-out difference in score (SD of score); sections = those with >= 3 Tregs in CN'
                               + ('' if tested else f'; not tested (< 100 Tregs or < 5 sections with >= 10 Tregs: {int(r.n_sections_ge10_tregs)}); excluded from BH'))))
    def add(block, feats, RES, JNT, label, model):
        for c in feats:
            v = RES[c]; j = JNT.get(c) if JNT is not None else None
            note = []
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
    add('local cell type', TYPES, RT, JT, 'type', 'LMM: score ~ type_z + cz + (1|section); type_z = fraction within 50 um, z within section')
    for c in [c for c in TYPES_ALL if c not in TYPES]:
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
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_colwidth', 40)
    print(O.drop(columns=['model', 'note']).to_string(), flush=True)
    with open(CSV, 'w') as fh: O.to_csv(fh, index=False)
    print('wrote', CSV, flush=True)
    # ---------------- strict-criteria summary (printed only) ----------------
    def strict(block, feats, RES):
        for c in feats:
            v = RES[c]
            if v['q'] < 0.05 and v['agree'] >= 9 and v['cl_p'] < 0.05:
                print(f'STRICT {block}: {c} {v["coef"]:+.3f} q {v["q"]:.2g} {v["agree"]}/{v["nsec"]} clustered P {v["cl_p"]:.2g}', flush=True)
    strict('cell type', TYPES, RT); strict('ligand', LIGS_ALL, RL)
    # ---------------- README ----------------
    jt = (f'fitted (max |Pearson r| between the {len(TYPES)} within-section z = {mrT:.2f}, {prT[0]} vs {prT[1]}; {njT:,} Tregs)' if JT is not None
          else f'NOT fitted (max |r| = {mrT:.2f} >= {MAX_R}, {prT[0]} vs {prT[1]})')
    jl = (f'fitted on the 5 lig_* columns of the table (max |r| = {mrL:.2f}, {prL[0]} vs {prL[1]}; {njL:,} Tregs); lig_CXCL12_exp left out of the joint model because it is an alternative definition of lig_CXCL12 (r = 0.80)' if JL is not None else f'NOT fitted (max |r| = {mrL:.2f} >= {MAX_R})')
    excl = [c for c in TYPES_ALL if c not in TYPES]
    TXT = f"""cd_spatial_results_all.csv - all CD resection spatial results, every feature of each block tested the same way
Cohort: Crohn's disease resections, CosMx 6K, {N_SEC} sections from 4 patients, {N_ALL:,} strict Tregs. Script: cd_spatial_results_all.py (this folder).
Outcome: Treg Suppression Score (act_z; negated D6 late-vs-early pseudotime signature, z-scored within section).
Rows are in neutral order within each block (neighbourhood: CN number; cell type and ligand: alphabetical). No ranking, no selection.

Blocks
1. neighbourhood (10 rows): all 10 cellular neighbourhoods (CN1-CN10; k-means k = 10 on 20-nearest-neighbour composition), copied as is from
   fig4_cn_treg_state_forest.csv (Fig 4). Model: linear mixed model, score ~ in-CN indicator + (1|section), all {N_ALL:,} Tregs; effect = in-minus-out
   difference (SD of score). BH across the 8 tested CNs; CN6 and CN9 (< 100 Tregs or < 5 sections with >= 10 Tregs) were pre-specified as not tested
   and have P but no q. sections_same_sign = sections (with >= 3 Tregs in the CN) where the in-CN mean is on the pooled side. No library-size
   covariate and no clustered P in that analysis (blank).
2. local cell type ({len(TYPES)} rows tested{', ' + str(len(excl)) + ' not tested' if excl else ''}): every cell type of the composition table except Treg
   (fraction > 0 near >= {MIN_PREV:g}% of Tregs; {'all passed, lowest ' + min(TYPES, key=lambda c: PREV[c]) + ' ' + format(min(PREV[c] for c in TYPES), '.3g') + '%' if not excl else 'excluded: ' + ', '.join(excl)}).
   type_z = the Treg's local fraction of that type among all cells within a 50 um radius, z-scored within section.
3. local ligand availability ({len(LIGS_ALL)} rows): every lig_* column of treg_local_ligand_receptor.csv.gz (CCL22, CXCL12, CXCL14, IL6, TNF: mean
   within-section z of the ligand over all cells within 50 um, >= 5 cells) plus lig_CXCL12_exp (CXCL12 over expected stromal producers, the
   definition behind the Fig 4F CXCL12 regions). Re-z-scored within section across Tregs. These are ligand expression near the Treg, not
   ligand-receptor signalling.
4. APC x CXCL12 (2 rows): Fig 4F, re-fitted here (interaction asserted equal to the CXCL12 interaction table; interaction P = {ri['p']:.2g}; difference-in-differences
   P = {rd['p']:.2g}). Single pre-specified tests, not part of a BH family (q blank).

Models (blocks 2 and 3; each feature separately)
  separate: act_z ~ x_z + cz + (1|section); linear mixed model, REML; effect_sd = SD change in score per SD of the feature; 95% Wald CI.
  cz = log1p(library size of the Treg) z-scored within section (library-size covariate).
  clustered_p: OLS act_z ~ x_z + cz with section-clustered SE, t distribution with G - 1 = {N_SEC - 1} df (robustness).
  sections_same_sign: per-section OLS act_z ~ x_z + cz; sections whose coefficient has the sign of the pooled estimate (x / sections estimable).
  joint: act_z ~ all x_z of the block + cz + (1|section); fitted only if max |Pearson r| between the within-section z < {MAX_R}.
         Cell types: {jt}. Ligands: {jl}. Joint 95% CI in the note column.
BH families: q = Benjamini-Hochberg across all tested features of the block (cell types: {len(TYPES)}; ligands: {len(LIGS_ALL)}; neighbourhoods: 8);
  joint_q = BH across all features of the joint model of that block (cell types: {len(TYPES)}; ligands: {len(LIGS)}).
Note on P: score and features are z-scored within section, so the section random-intercept variance is estimated at or near zero
  (DC model: {RE_VAR:.1e}; boundary warnings from statsmodels) and the LMM P is close to a pooled OLS P that treats Tregs as independent.
  clustered_p allows for correlation of Tregs within a section and is the more conservative column.
Note on the joint cell-type model: local fractions are compositional (they sum to 1 with Treg), so joint estimates are conditional on all other
  types and can differ in sign from the separate estimates.
Radius: 50 um around each Treg for cell-type fractions and ligand availability (the Treg itself is counted among the cells within 50 um).
Every feature in each block was tested and is reported, significant or not. Caveat: sections are nested in 4 patients; section is the
  random effect and the cluster unit, so patient-level dependence is not modelled separately.
Numbers: 3 significant figures.
"""
    with open(README, 'w') as fh: fh.write(TXT)
    print(TXT, flush=True)


# ======================================================================================================================
# Section interaction_tests: DC x OGN+RSPO3+ Fib and DC x B cell interactions, and SELENOP+ Mac subtypes, on the Treg Suppression Score
# Panels: Fig 5D interaction P, Results text (cd_interaction_selenop_tests.csv)
# ======================================================================================================================
if _run('interaction_tests'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/updated_figures/cd_interaction_selenop_tests.py'   # working folder of this step (relative paths resolve from here)
    """Interaction and subtype tests on the CD resection Treg Suppression Score (numbers only, one CSV).
CD resection CosMx 6K, 11 sections from 4 patients, 6,554 strict Tregs. Loading as in spatial_results_table:
  act_z = Treg Suppression Score (z within section); <type>_z = local fraction among all cells within 50 um (treg_local_composition.csv.gz of
  11_neighbourhoods_local_tables.py), z-scored within section; cz = log1p(library size) z within section.
Primary model: OLS with section fixed effects, section-clustered SE, t(G - 1) inference (G = 11), always + cz.
Reference: MixedLM with section random intercept (REML, lbfgs), + cz.  sections_agree: per-section OLS (same formula, no FE), sign of the pooled OLS estimate.
  1. DC x OGN+RSPO3+ Fib: act_z ~ DC_z * OGN_z + cz; DC slope in the bottom vs top within-section tertile of local OGN+RSPO3+ Fib fraction.
  2. DC x B cell: same with B cell.
  3. SELENOP+ Mac subtypes (cell_type of all_samples.parquet; the 4 labels whose cell_type_general is 'SELENOP+ Mac'): per-Treg fraction within 50 um
     with the neighbour definition of the composition table (cKDTree per section, all cells, global px * 0.12028 um), asserted against the merged
     column; separate and joint (with DC_z and CD9+ M2 Mac_z) models.
Output: cd_interaction_selenop_tests.csv in CDIX_OUT (default /share/fsmresfiles/UC/revision_jenny/updated_figures; 3 sig figs)."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    import statsmodels.formula.api as smf
    from scipy import stats
    from scipy.spatial import cKDTree
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    LR = f'{P4}/B_resection_analyses/ligand_receptor_state'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50
    OUT = os.environ.get('CDIX_OUT', '/share/fsmresfiles/UC/revision_jenny/updated_figures')
    CSV = f'{OUT}/cd_interaction_selenop_tests.csv'
    def zs(v):
        s = v.std()
        return (v - v.mean()) / s if s > 0 else v * 0.0
    # ---------------- data (as cd_spatial_results_all.py) ----------------
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    COMP = pd.read_csv(f'{LR}/treg_local_composition.csv.gz', index_col=0)
    assert set(COMP.index) == set(T.index); COMP = COMP.loc[T.index]
    assert (COMP.nb_n.values == T.nb_n.values).all() and (COMP.section.values == T.section.values).all()
    rs = COMP.drop(columns=['section', 'nb_n']).sum(1)
    assert np.allclose(rs[COMP.nb_n > 0], 1, atol=5e-3), 'composition is not stored as fractions'
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    assert len(cnt) == len(T)
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    T['cz'] = T.groupby('section').cz.transform(zs)
    COMP = COMP.loc[T.index]
    assert len(T) == 6554 and T.section.nunique() == 11
    G = T.section.nunique()
    # ---------------- SELENOP+ Mac fine subtypes within 50 um (neighbour definition of the composition table) ----------------
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type', 'cell_type_short',
                                                                      'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM
    for c in ['cell_type', 'cell_type_short', 'cell_type_general']: sa[c] = sa[c].astype(str)
    assert T.index.isin(sa.index).all()
    MAC = sa[sa.cell_type_general.str.contains('Mac')].groupby(['cell_type_general', 'cell_type', 'cell_type_short']).size()
    print('macrophage labels (cell_type_general / cell_type / cell_type_short): n cells'); print(MAC.to_string(), flush=True)
    SUB = sorted(sa.loc[sa.cell_type_general == 'SELENOP+ Mac', 'cell_type'].unique())
    assert SUB == ['SELENOP+MMP9+PLA2G2D+ Mac', 'SELENOP+MMP9+SPP1+ Mac', 'SELENOP+SIGLEC1+ Mac', 'SELENOP+XIST- Mac'], SUB
    SK = {'SELENOP+MMP9+PLA2G2D+ Mac': 's_mp', 'SELENOP+MMP9+SPP1+ Mac': 's_ms', 'SELENOP+SIGLEC1+ Mac': 's_sg', 'SELENOP+XIST- Mac': 's_xs'}
    F = pd.DataFrame(0.0, index=T.index, columns=list(SK.values())); NB = pd.Series(0, index=T.index)
    for s_, tg in T.groupby('section'):
        cs = sa[sa.section == s_]; code = pd.Categorical(cs.cell_type, categories=SUB).codes      # -1 = not a SELENOP+ Mac subtype
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
        M = np.zeros((len(tg), len(SUB)))
        for k, i in enumerate(nb):
            cc = code[i]; M[k] = np.bincount(cc[cc >= 0], minlength=len(SUB))
        n = np.array([len(i) for i in nb]); NB.loc[tg.index] = n
        F.loc[tg.index] = M / n[:, None]
    del sa; gc.collect()
    assert (NB.values == COMP.nb_n.values).all(), 'neighbour counts differ from the composition table'
    dif = (F.sum(1) - COMP['SELENOP+ Mac']).abs()
    AGREE_TXT = (f'neighbour counts identical to the composition table for all {len(T):,} Tregs; sum of the 4 subtype fractions vs merged SELENOP+ Mac: '
             f'max |diff| = {dif.max():.2g}, Pearson r = {np.corrcoef(F.sum(1), COMP["SELENOP+ Mac"])[0, 1]:.6f} (table stored at 4 sig figs)')
    assert dif.max() < 5e-4, dif.max()
    print('CHECK', AGREE_TXT, flush=True)
    PREV = {k: 100 * (F[k] > 0).mean() for k in F.columns}
    # ---------------- analysis frame ----------------
    X = pd.DataFrame(index=T.index); X['section'] = T.section; X['act_z'] = T.act_z; X['cz'] = T.cz
    RAW = {'dc': 'DC', 'ogn': 'OGN+RSPO3+ Fib', 'bc': 'B cell', 'cd9': 'CD9+ M2 Mac', 'selmac': 'SELENOP+ Mac'}
    for k, c in RAW.items(): X[k + '_f'] = COMP[c].values
    for k in SK.values(): X[k + '_f'] = F[k].values
    for k in list(RAW) + list(SK.values()): X[k] = X.groupby('section')[k + '_f'].transform(zs)
    X = X.dropna(); assert len(X) == 6554
    # ---------------- model helpers ----------------
    def fit(d, rhs, terms, agree=True):
        """clustered OLS with section FE (primary) + MixedLM (reference) + per-section sign agreement."""
        o = smf.ols(f'act_z ~ {rhs} + C(section)', d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]))
        g = d.section.nunique(); tc = stats.t.ppf(0.975, g - 1)
        m = smf.mixedlm(f'act_z ~ {rhs}', d, groups=d.section).fit(reml=True, method='lbfgs')
        sec = {}
        if agree:
            for s_, gd in d.groupby('section'):
                if len(gd) < 10 or any(gd[t.split(':')[0]].std() == 0 or gd[t.split(':')[-1]].std() == 0 for t in terms): continue
                sec[s_] = smf.ols(f'act_z ~ {rhs}', gd).fit().params[terms]
        out = {}
        for t in terms:
            b, se = o.params[t], o.bse[t]
            a = sum(int(np.sign(v[t]) == np.sign(b)) for v in sec.values())
            out[t] = dict(estimate=b, lo=b - tc * se, hi=b + tc * se, clustered_p=2 * stats.t.sf(abs(b / se), g - 1), mixedlm_p=m.pvalues[t],
                          mixedlm_est=m.params[t], sections_agree=f'{a}/{len(sec)}' if agree else '', n=len(d))
        return out
    def tert(v):  # within-section tertile, rank first so that ties (many 0 fractions) do not collapse bins
        return pd.qcut(v.rank(method='first'), 3, labels=['bottom', 'middle', 'top'])
    rows = []
    def add(block, term, r, note=''):
        rows.append(dict(block=block, term=term, **{k: r[k] for k in ['estimate', 'lo', 'hi', 'clustered_p', 'mixedlm_p', 'sections_agree', 'n']}, note=note))
    # ---------------- tests 1 and 2: DC x OGN+RSPO3+ Fib, DC x B cell ----------------
    for blk, mk, lab in [('1 DC x OGN+RSPO3+ Fib', 'ogn', 'OGN+RSPO3+ Fib'), ('2 DC x B cell', 'bc', 'B cell')]:
        r = fit(X, f'dc * {mk} + cz', ['dc', mk, f'dc:{mk}'])
        add(blk, 'DC_z (at average local ' + lab + ')', r['dc'], 'act_z ~ DC_z * x_z + cz')
        add(blk, lab + '_z (at average local DC)', r[mk])
        add(blk, f'DC_z x {lab}_z interaction', r[f'dc:{mk}'], 'negative = DC effect weaker where ' + lab + ' is locally enriched')
        X['t'] = X.groupby('section')[mk + '_f'].transform(tert)
        for tt in ['bottom', 'top']:
            d = X[X.t == tt]
            z0 = 100 * (d[mk + '_f'] == 0).mean()
            rr = fit(d, 'dc + cz', ['dc'])['dc']
            add(blk, f'DC_z slope, {tt} tertile of local {lab}', rr,
                f'act_z ~ DC_z + cz in Tregs of the {tt} within-section tertile of local {lab} fraction; {z0:.3g}% of these Tregs have 0 {lab} within 50 um')
        d = X[X.t.isin(['bottom', 'top'])].copy(); d['top'] = (d.t == 'top').astype(float); d['dctop'] = d.dc * d.top
        rr = fit(d, 'dc + top + dctop + cz', ['dctop'])['dctop']
        add(blk, f'DC_z slope difference, top minus bottom tertile of {lab}', rr, 'act_z ~ DC_z * top_tertile + cz, bottom + top tertile Tregs')
    X = X.drop(columns='t')
    # ---------------- test 3: SELENOP+ Mac fine subtypes ----------------
    LAB = {v: k for k, v in SK.items()}
    blk = '3 SELENOP+ Mac subtypes'
    r = fit(X, 'selmac + cz', ['selmac'])['selmac']
    add(blk, 'SELENOP+ Mac (merged)_z, separate', r, f'present near {100 * (X.selmac_f > 0).mean():.3g}% of Tregs; act_z ~ x_z + cz')
    for k in SK.values():
        r = fit(X, f'{k} + cz', [k])[k]
        add(blk, f'{LAB[k]}_z, separate', r, f'present near {PREV[k]:.3g}% of Tregs; act_z ~ x_z + cz')
    JK = list(SK.values()) + ['dc', 'cd9']
    Rz = X[JK].corr().values; mr = np.abs(Rz[np.triu_indices(len(JK), 1)]).max()
    rj = fit(X, ' + '.join(JK) + ' + cz', JK)
    for k in JK:
        add(blk, f'{LAB.get(k, RAW.get(k))}_z, joint', rj[k],
            f'act_z ~ 4 SELENOP+ Mac subtype_z + DC_z + CD9+ M2 Mac_z + cz; max |r| between the 6 predictors = {mr:.2f}')
    rows.append(dict(block=blk, term='check: subtype fractions sum to merged SELENOP+ Mac', n=len(X), note=AGREE_TXT))
    # ---------------- output ----------------
    O = pd.DataFrame(rows)
    def g3(x): return '' if x is None or (isinstance(x, float) and np.isnan(x)) else f'{x:.3g}'
    O['ci'] = [f'{g3(a)} to {g3(b)}' if pd.notna(a) else '' for a, b in zip(O.lo, O.hi)]
    for c in ['estimate', 'clustered_p', 'mixedlm_p']: O[c] = O[c].map(g3)
    O['n'] = O.n.map(lambda v: '' if pd.isna(v) else str(int(v)))
    O['sections_agree'] = O.sections_agree.fillna('')
    O = O[['block', 'term', 'estimate', 'ci', 'clustered_p', 'mixedlm_p', 'sections_agree', 'n', 'note']]
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
    print(O.drop(columns='note').to_string(), flush=True)
    with open(CSV, 'w') as fh: O.to_csv(fh, index=False)
    print('wrote', CSV, flush=True)


# ======================================================================================================================
# Section F5D_dc_by_ogn_tertile: local DC amount and DC slope by within-section tertile of local OGN+RSPO3+ fibroblasts (numbers only)
# Panels: Fig 5D right (DC % row; DC effect by OGN+RSPO3+ tertile)
# ======================================================================================================================
if _run('F5D_dc_by_ogn_tertile'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/104c_dc_ogn_tertile_checks.py'   # working folder of this step (relative paths resolve from here)
    """Local DC amount and DC slope by within-section tertile of local OGN+RSPO3+ Fib fraction (numbers only, no figure). The pooled mean DC fraction
per tertile is the 'DC %' row of the Fig 5D right panel (F5D_apc_types_dc_ogn).
Question: is the local DC amount comparable across tertiles? If the high-OGN tertile has fewer DCs / a narrower DC range, the weaker DC slope there could
reflect range restriction or non-linearity rather than OGN+RSPO3+ Fib modifying the DC association.
Data loading, z-scoring, tertiles and models as in F5D_apc_types_dc_ogn and interaction_tests (test 1).
  1. Descriptives of local DC fraction (% of all cells within 50 um; and dc = z within section) per OGN tertile, pooled and per section:
     mean, median, SD, IQR, % of Tregs with zero DCs; sd_within = SD of dc after removing section means (the variation the FE slope uses).
  2. Sensitivity A (DC-matched): within each section, DC fraction binned into deciles (pd.qcut, tied edges dropped, so zeros form one bin);
     inverse-probability weights w = p_section(bin) / p_section,tertile(bin) give every tertile the section's pooled DC-bin distribution;
     weighted OLS act_z ~ dc + cz + C(section), section-clustered SE, t(G - 1). Balance is reported as weighted DC descriptives.
     A2 (common support): within each section keep Tregs whose DC fraction lies within the overlap of the three tertiles' DC ranges
     [max of tertile minima, min of tertile maxima]; A3 = same with the tertiles' 5th-95th percentile ranges. Unweighted refit.
     For each: Low - High slope difference from one fully interacted model (section x tertile FE; dc, cz slopes by tertile), clustered, t(G - 1).
  3. Sensitivity B (non-linearity): act_z ~ dc * ogn + I(dc**2) + cz + C(section), and + I(dc**2):ogn; clustered, t(G - 1); per-section sign agreement.
Output: Fig4A-F_treg_suppression_niches/fig4d_dc_by_ogn_tertile_checks.csv (3 sig figs)."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    import statsmodels.formula.api as smf
    from scipy import stats
    from scipy.spatial import cKDTree
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; LR = f'{D}/panels/fig4_and_related/B_resection_analyses/ligand_receptor_state'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50
    OUT = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/fig4d_dc_by_ogn_tertile_checks.csv'
    # ---------------- APC types: key -> (source table, label(s), display name), as in F5D_apc_types_dc_ogn ----------------
    TYPES = {'dc': ('comp', 'DC', 'DC'),
             's_mp': ('fine', 'SELENOP+MMP9+PLA2G2D+ Mac', 'SELENOP⁺MMP9⁺PLA2G2D⁺ Mac'),
             's_ms': ('fine', 'SELENOP+MMP9+SPP1+ Mac', 'SELENOP⁺MMP9⁺SPP1⁺ Mac'),
             's_sgxs': ('fine', ['SELENOP+SIGLEC1+ Mac', 'SELENOP+XIST- Mac'], 'SELENOP⁺ Mac (SIGLEC1⁺/XIST⁻)'),
             'cd9': ('comp', 'CD9+ M2 Mac', 'CD9⁺ M2 Mac'),
             'momac': ('comp', 'Mo-Mac', 'Mo-Mac'),
             'bc': ('comp', 'B cell', 'B cell')}
    K = list(TYPES); SUB = ['SELENOP+MMP9+PLA2G2D+ Mac', 'SELENOP+MMP9+SPP1+ Mac', 'SELENOP+SIGLEC1+ Mac', 'SELENOP+XIST- Mac']
    def zs(v):
        s = v.std()
        return (v - v.mean()) / s if s > 0 else v * 0.0
    # ---------------- data (as in F5D_apc_types_dc_ogn) ----------------
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)[['section', 'act_z', 'nb_n', 'apc_frac']]
    COMP = pd.read_csv(f'{LR}/treg_local_composition.csv.gz', index_col=0)
    assert set(COMP.index) == set(T.index); COMP = COMP.loc[T.index]
    assert (COMP.nb_n.values == T.nb_n.values).all() and (COMP.section.values == T.section.values).all()
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    assert len(cnt) == len(T)
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    T['cz'] = T.groupby('section').cz.transform(zs); COMP = COMP.loc[T.index]
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type', 'cell_type_general', 'patient_sec',
                                                                      '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['cell_type'] = sa.cell_type.astype(str)
    assert sorted(sa.loc[sa.cell_type_general.astype(str) == 'SELENOP+ Mac', 'cell_type'].unique()) == sorted(SUB)
    F = pd.DataFrame(0.0, index=T.index, columns=SUB); NB = pd.Series(0, index=T.index)
    for s_, tg in T.groupby('section'):
        cs = sa[sa.section == s_]; code = pd.Categorical(cs.cell_type, categories=SUB).codes
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
        M = np.zeros((len(tg), len(SUB)))
        for k, i in enumerate(nb):
            cc = code[i]; M[k] = np.bincount(cc[cc >= 0], minlength=len(SUB))
        n = np.array([len(i) for i in nb]); NB.loc[tg.index] = n; F.loc[tg.index] = M / n[:, None]
    del sa; gc.collect()
    assert (NB.values == COMP.nb_n.values).all() and (F.sum(1) - COMP['SELENOP+ Mac']).abs().max() < 5e-4
    for k in K:
        T[k + '_f'] = COMP[TYPES[k][1]].values if TYPES[k][0] == 'comp' else F[TYPES[k][1]].sum(1).values if isinstance(TYPES[k][1], list) else F[TYPES[k][1]].values
        T[k] = T.groupby('section')[k + '_f'].transform(zs)
    assert np.allclose(T[[k + '_f' for k in K]].sum(1), T.apc_frac, atol=2e-3), 'rows do not resolve the APC fraction (apc_frac)'
    T['ogn_f'] = COMP['OGN+RSPO3+ Fib'].values; T['ogn'] = T.groupby('section').ogn_f.transform(zs)
    d = T[['section', 'act_z', 'cz', 'ogn', 'ogn_f'] + K + [k + '_f' for k in K]].dropna(); n_t, n_s = len(d), d.section.nunique()
    assert n_t == 6554 and n_s == 11 and len(K) == 7
    def clfit(dd, rhs, w=None):
        g = pd.factorize(dd.section)[0]
        if w is None: return smf.ols(f'act_z ~ {rhs} + C(section)', dd).fit(cov_type='cluster', cov_kwds=dict(groups=g))
        return smf.wls(f'act_z ~ {rhs} + C(section)', dd, weights=w).fit(cov_type='cluster', cov_kwds=dict(groups=g))
    def summ(o, t, g):
        b, se = o.params[t], o.bse[t]; tq = stats.t.ppf(0.975, g - 1)
        return dict(estimate=b, lo=b - tq * se, hi=b + tq * se, p=2 * stats.t.sf(abs(b / se), g - 1))
    def agree(dd, rhs, t, b):
        a = n = 0
        for s_, gd in dd.groupby('section'):
            if len(gd) < 10 or any(gd[c].std() == 0 for c in t.split(':')): continue
            n += 1; a += int(np.sign(smf.ols(f'act_z ~ {rhs}', gd).fit().params[t]) == np.sign(b))
        return f'{a}/{n}'
    IX = summ(clfit(d, 'dc * ogn + cz'), 'dc:ogn', n_s); IX['agree'] = agree(d, 'dc * ogn + cz', 'dc:ogn', IX['estimate'])
    def tert(v):  # as in interaction_tests
        return pd.qcut(v.rank(method='first'), 3, labels=['bottom', 'middle', 'top'])
    d['t'] = d.groupby('section').ogn_f.transform(tert)
    TT = ['bottom', 'middle', 'top']; TL = dict(bottom='Low', middle='Mid', top='High')
    TER = {}
    for tt in TT:
        dt = d[d.t == tt]; g = dt.section.nunique()
        TER[tt] = summ(clfit(dt, 'dc + cz'), 'dc', g); TER[tt]['n'] = len(dt); TER[tt]['g'] = g
        TER[tt]['agree'] = agree(dt, 'dc + cz', 'dc', TER[tt]['estimate'])
    d['dcp'] = 100 * d.dc_f   # DC fraction in % of all cells within 50 um
    d['dc_w'] = d.dc - d.groupby('section').dc.transform('mean')   # within-section deviation (pooled over tertiles)
    ROWS = []
    def wq(x, w, q):
        o = np.argsort(x); x, w = x[o], w[o]; c = np.cumsum(w) / w.sum()
        return x[np.searchsorted(c, q)]
    def desc(dd, w=None):
        """DC descriptives (weighted if w). sd_within = SD of dc within section x tertile (section means of this subset removed)."""
        w = np.ones(len(dd)) if w is None else np.asarray(w, float); r = {}
        for v, lab in [('dcp', 'dcpct'), ('dc', 'dcz')]:
            x = dd[v].values; m = np.average(x, weights=w)
            r[f'{lab}_mean'] = m; r[f'{lab}_median'] = wq(x, w, 0.5); r[f'{lab}_sd'] = np.sqrt(np.average((x - m) ** 2, weights=w))
            r[f'{lab}_q1'] = wq(x, w, 0.25); r[f'{lab}_q3'] = wq(x, w, 0.75); r[f'{lab}_iqr'] = r[f'{lab}_q3'] - r[f'{lab}_q1']
        sm = pd.Series(w * dd.dc.values, index=dd.index).groupby(dd.section).transform('sum') / pd.Series(w, index=dd.index).groupby(dd.section).transform('sum')
        r['dcz_sd_within'] = np.sqrt(np.average((dd.dc.values - sm.values) ** 2, weights=w))
        r['pct_zero_dc'] = 100 * np.average(dd.dc_f.values == 0, weights=w); r['ogn_pct_mean'] = 100 * np.average(dd.ogn_f.values, weights=w)
        return r
    # ---------------- 1. descriptives ----------------
    for tt in TT:
        ROWS.append(dict(block='1_desc_pooled', section='all', tertile=TL[tt], n=int((d.t == tt).sum()), **desc(d[d.t == tt])))
    for tt in TT:   # Kruskal-Wallis across tertiles is reported per section on the first row of each section block (in p)
        pass
    for s_, gs in d.groupby('section'):
        kw = stats.kruskal(*[gs.loc[gs.t == tt, 'dcp'] for tt in TT]).pvalue
        for tt in TT:
            ROWS.append(dict(block='1_desc_section', section=s_, tertile=TL[tt], n=int((gs.t == tt).sum()), **desc(gs[gs.t == tt]),
                             note=f'Kruskal-Wallis DC fraction across tertiles P = {kw:.3g}' if tt == 'bottom' else ''))
    # section-clustered test of whether mean DC (z) differs by tertile (High vs Low)
    o = smf.ols('dc ~ C(t, Treatment("bottom")) + C(section)', d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]))
    for tt in ['middle', 'top']:
        r = summ(o, f'C(t, Treatment("bottom"))[T.{tt}]', n_s)
        ROWS.append(dict(block='1_desc_test', section='all', tertile=f'{TL[tt]} - Low', term='mean dc (z) difference, section FE, clustered', n=n_t, **r))
    # ---------------- 2. tertile slopes: unweighted, DC-matched, common support ----------------
    def slopes(dd, label, w=None, note=''):
        out = {}
        for tt in TT:
            m = (dd.t == tt).values; dt = dd[m]; g = dt.section.nunique()
            r = summ(clfit(dt, 'dc + cz', None if w is None else w[m]), 'dc', g)
            r['agree'] = agree(dt, 'dc + cz', 'dc', r['estimate']) if w is None else ''
            out[tt] = r
            ROWS.append(dict(block=label, section='all', tertile=TL[tt], term='dc', n=int(m.sum()), **r,
                             **desc(dt, None if w is None else w[m]), note=note))
        # Low - High slope difference: fully interacted model (section x tertile FE, dc and cz slopes by tertile)
        dd = dd.assign(st=dd.section.astype(str) + '|' + dd.t.astype(str))
        f = 'act_z ~ 0 + C(st) + dc:C(t, Treatment("bottom")) + cz:C(t, Treatment("bottom"))'
        g_ = pd.factorize(dd.section)[0]
        oo = (smf.ols(f, dd) if w is None else smf.wls(f, dd, weights=w)).fit(cov_type='cluster', cov_kwds=dict(groups=g_))
        nm = {tt: f'dc:C(t, Treatment("bottom"))[{tt}]' for tt in TT}
        for tt in TT: assert abs(oo.params[nm[tt]] - out[tt]['estimate']) < 1e-8
        for a, b in [('bottom', 'top'), ('bottom', 'middle')]:
            L = np.zeros(len(oo.params)); L[list(oo.params.index).index(nm[a])] = 1; L[list(oo.params.index).index(nm[b])] = -1
            est = L @ oo.params.values; se = np.sqrt(L @ oo.cov_params().values @ L); tq = stats.t.ppf(0.975, n_s - 1)
            ROWS.append(dict(block=label, section='all', tertile=f'{TL[a]} - {TL[b]}', term='dc slope difference', n=len(dd),
                             estimate=est, lo=est - tq * se, hi=est + tq * se, p=2 * stats.t.sf(abs(est / se), n_s - 1), note=note))
        return out
    S0 = slopes(d, '2_slopes_unweighted')
    for tt in TT: assert abs(S0[tt]['estimate'] - TER[tt]['estimate']) < 1e-12
    # A: IPW to the section's pooled DC-decile distribution
    d['bin'] = d.groupby('section').dc_f.transform(lambda v: pd.qcut(v, 10, labels=False, duplicates='drop'))
    ps = d.groupby(['section', 'bin']).size() / d.groupby('section').size()
    pst = d.groupby(['section', 't', 'bin'], observed=True).size() / d.groupby(['section', 't'], observed=True).size()
    idx_s = pd.MultiIndex.from_arrays([d.section, d.bin]); idx_st = pd.MultiIndex.from_arrays([d.section, d.t, d.bin])
    W = ps.reindex(idx_s).values / pst.reindex(idx_st).values
    assert np.isfinite(W).all() and (W > 0).all()
    # pooled bin mass with no Tregs in some tertile (unmatched)
    full = pd.MultiIndex.from_product([sorted(d.section.unique()), TT, sorted(d.bin.unique())], names=['section', 't', 'bin'])
    cnt_st = d.groupby(['section', 't', 'bin'], observed=True).size().reindex(full, fill_value=0)
    miss = [(s_, tt, b) for (s_, tt, b), c in cnt_st.items() if c == 0 and (s_, b) in ps.index]
    miss_mass = 100 * sum(ps[(s_, b)] for s_, tt, b in miss) / (n_s * 3)
    nbins = d.groupby('section').bin.nunique()
    noteA = (f'IPW within section to pooled DC-fraction decile distribution (tied edges dropped; {nbins.min()}-{nbins.max()} bins per section); '
         f'weight range {W.min():.3g}-{W.max():.3g}; unmatched pooled bin mass {miss_mass:.3g}%')
    SA = slopes(d, '2A_slopes_dc_matched_ipw', w=W, note=noteA)
    # A2 / A3: common support within section
    def support(dd, qlo, qhi):
        keep = pd.Series(False, index=dd.index); info = []
        for s_, gs in dd.groupby('section'):
            lo = max(gs.loc[gs.t == tt, 'dc_f'].quantile(qlo) for tt in TT); hi = min(gs.loc[gs.t == tt, 'dc_f'].quantile(qhi) for tt in TT)
            keep.loc[gs.index] = (gs.dc_f >= lo) & (gs.dc_f <= hi)
        return dd[keep]
    for lab, qlo, qhi, txt in [('2B_slopes_common_support_minmax', 0, 1, 'overlap of tertile min-max DC-fraction ranges'),
                               ('2C_slopes_common_support_p5_p95', 0.05, 0.95, 'overlap of tertile 5th-95th percentile DC-fraction ranges')]:
        dd = support(d, qlo, qhi).copy()
        ROWS_note = f'within section, {txt}; kept {len(dd):,}/{n_t:,} Tregs'
        out = slopes(dd, lab, note=ROWS_note)
        if lab.startswith('2B'): SB = out; nB = len(dd)
        else: SC = out; nC = len(dd)
    # ---------------- 3. non-linearity ----------------
    NL = {}
    for lab, rhs in [('3_nonlinear_dc2', 'dc * ogn + I(dc**2) + cz'), ('3_nonlinear_dc2_x_ogn', 'dc * ogn + I(dc**2) + I(dc**2):ogn + cz')]:
        o = clfit(d, rhs); NL[lab] = {}
        terms = ['dc', 'ogn', 'dc:ogn', 'I(dc ** 2)'] + (['I(dc ** 2):ogn'] if 'I(dc**2):ogn' in rhs else [])
        for tm in terms:
            r = summ(o, tm, n_s); NL[lab][tm] = r
            ag = agree(d.assign(dc2=d.dc ** 2), rhs.replace('I(dc**2)', 'dc2'), tm.replace('I(dc ** 2)', 'dc2'), r['estimate'])
            ROWS.append(dict(block=lab, section='all', tertile='continuous', term=tm, n=n_t, **r, agree=ag, note=f'act_z ~ {rhs} + C(section)'))
    r = dict(IX); ROWS.append(dict(block='3_reference_linear', section='all', tertile='continuous', term='dc:ogn', n=n_t, **r,
                                   note='act_z ~ dc * ogn + cz + C(section) (main model)'))
    # ---------------- write ----------------
    R_ = pd.DataFrame(ROWS)
    cols = ['block', 'section', 'tertile', 'term', 'n', 'estimate', 'lo', 'hi', 'p', 'agree', 'ogn_pct_mean', 'dcpct_mean', 'dcpct_median', 'dcpct_sd',
            'dcpct_q1', 'dcpct_q3', 'dcpct_iqr', 'pct_zero_dc', 'dcz_mean', 'dcz_median', 'dcz_sd', 'dcz_q1', 'dcz_q3', 'dcz_iqr', 'dcz_sd_within', 'note']
    R_ = R_.reindex(columns=cols)
    def sig3(x):
        return x if not isinstance(x, (float, np.floating)) or not np.isfinite(x) else float(f'{x:.3g}')
    for c in cols:
        if R_[c].dtype.kind == 'f': R_[c] = R_[c].map(sig3)
    with open(OUT, 'w') as fh: R_.to_csv(fh, index=False)
    print('wrote', OUT, len(R_), 'rows', flush=True)
    # ---------------- summary ----------------
    P = R_[R_.block == '1_desc_pooled'].set_index('tertile')
    print('DC per tertile:', '; '.join(f'{t} mean {P.loc[t, "dcpct_mean"]}% median {P.loc[t, "dcpct_median"]}% SD {P.loc[t, "dcpct_sd"]} '
                                   f'zero {P.loc[t, "pct_zero_dc"]}% sd_within(z) {P.loc[t, "dcz_sd_within"]}' for t in ['Low', 'Mid', 'High']))
    print(R_[R_.block == '1_desc_test'][['tertile', 'estimate', 'lo', 'hi', 'p']].to_string(index=False))
    for b in ['2_slopes_unweighted', '2A_slopes_dc_matched_ipw', '2B_slopes_common_support_minmax', '2C_slopes_common_support_p5_p95']:
        x = R_[R_.block == b]
        print(b, ' | '.join(f'{r.tertile} {r.estimate:+} [{r.lo}, {r.hi}] P {r.p} n {r.n} {r.agree if isinstance(r.agree, str) else ""} dc% {r.dcpct_mean} zero {r.pct_zero_dc}'
                            for r in x.itertuples()), '|', x.note.iloc[0] if isinstance(x.note.iloc[0], str) else '')
    print(R_[R_.block.str.startswith('3')][['block', 'term', 'estimate', 'lo', 'hi', 'p', 'agree']].to_string(index=False))


# ======================================================================================================================
# Section F5D_apc_types_dc_ogn: Fig 5D, joint APC-type model (left) and DC slope by OGN+RSPO3+ tertile with per-section slopes (right)
# Panels: Fig 5D
# ======================================================================================================================
if _run('F5D_apc_types_dc_ogn'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/104b_apc_types_dc_ogn_lines.py'   # working folder of this step (relative paths resolve from here)
    """Fig 5D: local APC types and the Treg Suppression Score (CD resection CosMx 6K, 11 sections from 4 patients, 6,554 strict Tregs).
  Left  = one regression of the score on the local fractions of seven APC types (DC; SELENOP+MMP9+PLA2G2D+ Mac; SELENOP+MMP9+SPP1+ Mac; SELENOP+ Mac
          SIGLEC1+/XIST-; CD9+ M2 Mac; Mo-Mac; B cell): act_z ~ types + cz + C(section), OLS with section-clustered SE, t(G - 1), BH over the 7 types;
          section dots = per-section OLS coefficients, coloured by section as in the right panel (alpha 0.5); pooled markers and CIs with a thin white
          edge; strip = mean local fraction (%).
  Right = pooled DC slope per within-section tertile of local OGN+RSPO3+ fibroblast fraction (act_z ~ dc + cz + C(section), OLS, section-clustered SE,
          t(G - 1); filled = P < 0.05), overlaid on per-section DC slopes: OLS act_z ~ dc + cz within that section's Tregs in that tertile (dc, cz
          z-scored within section; plain OLS SE). Thin lines connect each section's Low / Mid / High slopes, one colour per section (11 colours; mapping
          in the CSV 'color' column). Patient (CSV only) = case ID = patient_sec prefix: 11808, 18888, 43563, 54033. Per-section signs are asserted
          against the pooled sign-agreement counts. 'DC %' row under the n row = pooled mean local DC fraction per tertile, read from
          fig4d_dc_by_ogn_tertile_checks.csv (F5D_dc_by_ogn_tertile). The DC x OGN+RSPO3+ interaction (estimate, clustered P, sections) is printed in the panel.
Outputs (APC_OUT, default Fig4A-F_treg_suppression_niches/): fig4_apc_types_dc_ogn_combined.png + _legend.txt, fig4d_dc_slope_by_ogn_tertile_per_section.csv.
Environment: YLO / YHI set the right-panel y range (default -0.3 / 0.4)."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
    from scipy import stats
    from scipy.spatial import cKDTree
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; LR = f'{D}/panels/fig4_and_related/B_resection_analyses/ligand_receptor_state'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50
    OUT = os.environ.get('APC_OUT', '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches') + '/fig4_apc_types_dc_ogn_combined'
    CSV = os.path.dirname(OUT) + '/fig4d_dc_slope_by_ogn_tertile_per_section.csv'
    SS = 'Treg Suppression Score'
    # ---------------- APC types: key -> (source table, label(s), display name) ----------------
    TYPES = {'dc': ('comp', 'DC', 'DC'),
             's_mp': ('fine', 'SELENOP+MMP9+PLA2G2D+ Mac', 'SELENOP⁺MMP9⁺PLA2G2D⁺ Mac'),
             's_ms': ('fine', 'SELENOP+MMP9+SPP1+ Mac', 'SELENOP⁺MMP9⁺SPP1⁺ Mac'),
             's_sgxs': ('fine', ['SELENOP+SIGLEC1+ Mac', 'SELENOP+XIST- Mac'], 'SELENOP⁺ Mac (SIGLEC1⁺/XIST⁻)'),
             'cd9': ('comp', 'CD9+ M2 Mac', 'CD9⁺ M2 Mac'),
             'momac': ('comp', 'Mo-Mac', 'Mo-Mac'),
             'bc': ('comp', 'B cell', 'B cell')}
    K = list(TYPES); SUB = ['SELENOP+MMP9+PLA2G2D+ Mac', 'SELENOP+MMP9+SPP1+ Mac', 'SELENOP+SIGLEC1+ Mac', 'SELENOP+XIST- Mac']
    def zs(v):
        s = v.std()
        return (v - v.mean()) / s if s > 0 else v * 0.0
    # ---------------- data (as in interaction_tests) ----------------
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)[['section', 'act_z', 'nb_n', 'apc_frac']]
    COMP = pd.read_csv(f'{LR}/treg_local_composition.csv.gz', index_col=0)
    assert set(COMP.index) == set(T.index); COMP = COMP.loc[T.index]
    assert (COMP.nb_n.values == T.nb_n.values).all() and (COMP.section.values == T.section.values).all()
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    assert len(cnt) == len(T)
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    T['cz'] = T.groupby('section').cz.transform(zs); COMP = COMP.loc[T.index]
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type', 'cell_type_general', 'patient_sec',
                                                                      '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['cell_type'] = sa.cell_type.astype(str)
    assert sorted(sa.loc[sa.cell_type_general.astype(str) == 'SELENOP+ Mac', 'cell_type'].unique()) == sorted(SUB)
    F = pd.DataFrame(0.0, index=T.index, columns=SUB); NB = pd.Series(0, index=T.index)
    for s_, tg in T.groupby('section'):
        cs = sa[sa.section == s_]; code = pd.Categorical(cs.cell_type, categories=SUB).codes
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
        M = np.zeros((len(tg), len(SUB)))
        for k, i in enumerate(nb):
            cc = code[i]; M[k] = np.bincount(cc[cc >= 0], minlength=len(SUB))
        n = np.array([len(i) for i in nb]); NB.loc[tg.index] = n; F.loc[tg.index] = M / n[:, None]
    del sa; gc.collect()
    assert (NB.values == COMP.nb_n.values).all() and (F.sum(1) - COMP['SELENOP+ Mac']).abs().max() < 5e-4
    for k in K:
        T[k + '_f'] = COMP[TYPES[k][1]].values if TYPES[k][0] == 'comp' else F[TYPES[k][1]].sum(1).values if isinstance(TYPES[k][1], list) else F[TYPES[k][1]].values
        T[k] = T.groupby('section')[k + '_f'].transform(zs)
    assert np.allclose(T[[k + '_f' for k in K]].sum(1), T.apc_frac, atol=2e-3), 'rows do not resolve the APC fraction (apc_frac)'
    T['ogn_f'] = COMP['OGN+RSPO3+ Fib'].values; T['ogn'] = T.groupby('section').ogn_f.transform(zs)
    d = T[['section', 'act_z', 'cz', 'ogn', 'ogn_f'] + K + [k + '_f' for k in K]].dropna(); n_t, n_s = len(d), d.section.nunique()
    assert n_t == 6554 and n_s == 11 and len(K) == 7
    tc = stats.t.ppf(0.975, n_s - 1)
    def clfit(dd, rhs):
        return smf.ols(f'act_z ~ {rhs} + C(section)', dd).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(dd.section)[0]))
    def summ(o, t, g):
        b, se = o.params[t], o.bse[t]; tq = stats.t.ppf(0.975, g - 1)
        return dict(estimate=b, lo=b - tq * se, hi=b + tq * se, p=2 * stats.t.sf(abs(b / se), g - 1))
    # ---------------- left: joint APC-type model ----------------
    RHS = ' + '.join(K) + ' + cz'
    o = clfit(d, RHS)
    RES = {k: summ(o, k, n_s) for k in K}
    q = multipletests([RES[k]['p'] for k in K], method='fdr_bh')[1]
    for k, qq in zip(K, q): RES[k]['q'] = qq
    SEC = pd.DataFrame({s_: smf.ols(f'act_z ~ {RHS}', g).fit().params[K] for s_, g in d.groupby('section')}).T
    assert len(SEC) == n_s
    for k in K: RES[k]['agree'] = int((np.sign(SEC[k]) == np.sign(RES[k]['estimate'])).sum())
    ABUND = {k: 100 * d[k + '_f'].mean() for k in K}
    for k in K:
        v = RES[k]
        print(f'{TYPES[k][2]:30s} {v["estimate"]:+.4f} [{v["lo"]:+.4f}, {v["hi"]:+.4f}] q {v["q"]:.3g} {v["agree"]}/{n_s}  '
          f'per-section {SEC[k].min():+.3f} to {SEC[k].max():+.3f}  mean local {ABUND[k]:.3g}%', flush=True)
    # ---------------- right: cd_interaction_selenop_tests.py test 1 ----------------
    def agree(dd, rhs, t, b):
        a = n = 0
        for s_, gd in dd.groupby('section'):
            if len(gd) < 10 or any(gd[c].std() == 0 for c in t.split(':')): continue
            n += 1; a += int(np.sign(smf.ols(f'act_z ~ {rhs}', gd).fit().params[t]) == np.sign(b))
        return f'{a}/{n}'
    IX = summ(clfit(d, 'dc * ogn + cz'), 'dc:ogn', n_s); IX['agree'] = agree(d, 'dc * ogn + cz', 'dc:ogn', IX['estimate'])
    def tert(v):  # as in interaction_tests
        return pd.qcut(v.rank(method='first'), 3, labels=['bottom', 'middle', 'top'])
    d['t'] = d.groupby('section').ogn_f.transform(tert)
    TER = {}
    for tt in ['bottom', 'middle', 'top']:
        dt = d[d.t == tt]; g = dt.section.nunique()
        TER[tt] = summ(clfit(dt, 'dc + cz'), 'dc', g); TER[tt]['n'] = len(dt); TER[tt]['g'] = g
        TER[tt]['agree'] = agree(dt, 'dc + cz', 'dc', TER[tt]['estimate']); TER[tt]['ogn_pct'] = 100 * dt.ogn_f.mean()
    assert sum(TER[t]['n'] for t in TER) == n_t and all(TER[t]['g'] == n_s for t in TER)
    print(f'interaction {IX["estimate"]:+.4f} [{IX["lo"]:+.4f}, {IX["hi"]:+.4f}] P {IX["p"]:.3g} {IX["agree"]}', flush=True)
    for tt in TER:
        v = TER[tt]
        print(f'DC slope {tt:6s} n={v["n"]} {v["estimate"]:+.4f} [{v["lo"]:+.4f}, {v["hi"]:+.4f}] P {v["p"]:.3g} {v["agree"]}  '
          f'mean local OGN+RSPO3+ Fib {v["ogn_pct"]:.3g}%', flush=True)
    # ---------------- per-section DC slope in each tertile ----------------
    PAT = lambda sec: sec.split('_')[0]
    rows = []
    for (s_, tt), gd in d.groupby(['section', 't'], observed=True):
        f_ = smf.ols('act_z ~ dc + cz', gd).fit()
        rows.append(dict(section=s_, patient=PAT(s_), tertile=dict(bottom='Low', middle='Mid', top='High')[tt], tertile_key=tt,
                         slope=f_.params['dc'], se=f_.bse['dc'], n_tregs=len(gd)))
    PS = pd.DataFrame(rows)
    assert len(PS) == 3 * n_s and PS.n_tregs.sum() == n_t and PS.patient.nunique() == 4
    for tt in ['bottom', 'middle', 'top']:   # signs must match the agree() counts of the pooled slopes
        p_ = PS[PS.tertile_key == tt]; a_ = int((np.sign(p_.slope) == np.sign(TER[tt]['estimate'])).sum())
        assert f'{a_}/{len(p_)}' == TER[tt]['agree'], (tt, a_, TER[tt]['agree'])
    WIDE = PS.pivot(index='section', columns='tertile_key', values='slope')
    N_LH = int((WIDE.bottom > WIDE.top).sum())
    print(f'per-section DC slope range {PS.slope.min():+.3f} to {PS.slope.max():+.3f}; Low > High in {N_LH}/{n_s} sections', flush=True)
    for tt in ['bottom', 'middle', 'top']:
        p_ = PS[PS.tertile_key == tt]
        print(f'  {tt:6s} per-section {p_.slope.min():+.3f} to {p_.slope.max():+.3f}, median {p_.slope.median():+.3f}, n Tregs {p_.n_tregs.min()}-{p_.n_tregs.max()}', flush=True)
    print(WIDE[['bottom', 'middle', 'top']].round(3).assign(patient=[PAT(x) for x in WIDE.index]).to_string(), flush=True)
    SECS = sorted(WIDE.index, key=lambda x: (PAT(x), x))
    SPAL = ['#1f77b4', '#ff7f0e', '#2ca02c', '#CC3311', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#DAA520', '#17255A']   # tab10 minus grey (red swapped for #CC3311) + gold, navy
    assert len(SECS) == len(SPAL) == len(set(SPAL)); SCOL = dict(zip(SECS, SPAL))
    PS['color'] = PS.section.map(SCOL)
    OUTCSV = PS[['section', 'patient', 'tertile', 'slope', 'se', 'n_tregs', 'color']].copy()
    OUTCSV['tertile'] = pd.Categorical(OUTCSV.tertile, ['Low', 'Mid', 'High']); OUTCSV = OUTCSV.sort_values(['patient', 'section', 'tertile'])
    for c in ['slope', 'se']: OUTCSV[c] = OUTCSV[c].map(lambda v: float(f'{v:.3g}'))
    with open(CSV, 'w') as fh: OUTCSV.to_csv(fh, index=False)
    print('wrote', CSV, flush=True)
    # ---------------- plot ----------------
    YLIM_R = (float(os.environ.get('YLO', -0.3)), float(os.environ.get('YHI', 0.4))); YTK_R = [-0.2, 0, 0.2]
    FKEY = 10; KEY_Y = 0.995
    FT, FL, FSm = 15, 15, 13          # tick/row labels, axis labels, small annotations (% / n / header / interaction)
    GREY = '#BDBDBD'; DARK = '#555555'
    plt.rcParams.update({'font.size': FT, 'font.family': 'DejaVu Sans', 'font.weight': 'normal', 'axes.titleweight': 'normal',
                         'axes.labelweight': 'normal', 'axes.linewidth': 1.0})
    W, H = 7.6, 4.75; fig = plt.figure(figsize=(W, H))   # includes 0.25 in at the bottom for the DC % row
    AH = 2.72   # axes height (in); row positions below are in inches
    def ax_in(x0, w, y0=1.57, h=AH): return fig.add_axes([x0 / W, y0 / H, w / W, h / H])
    axa, axf, axr = ax_in(1.92, 0.7), ax_in(2.72, 1.52), ax_in(5.5, 2.0)
    WRAP = {'s_mp': 'SELENOP⁺MMP9⁺\nPLA2G2D⁺ Mac', 's_ms': 'SELENOP⁺MMP9⁺\nSPP1⁺ Mac', 's_sgxs': 'SELENOP⁺ Mac\n(SIGLEC1⁺/XIST⁻)'}
    NL = {k: WRAP.get(k, TYPES[k][2]).count('\n') + 1 for k in K}; HALF = {1: 0.115, 2: 0.23}   # half label height (in) at 15 pt
    ypos = {K[0]: 0.0}
    for a, b in zip(K[:-1], K[1:]): ypos[b] = ypos[a] - (HALF[NL[a]] + HALF[NL[b]] + 0.045)   # uneven spacing: 2-line rows get more room
    YL = (ypos[K[-1]] - 0.2, 0.2)
    # abundance strip (header gives the unit, so values are printed without '%')
    amax = max(ABUND.values())
    for k in K:
        axa.barh(ypos[k], ABUND[k], height=0.11, color='#9E9E9E', lw=0)
        axa.text(ABUND[k] + 0.06 * amax, ypos[k], f'{ABUND[k]:.1f}', ha='left', va='center', fontsize=FSm, color=DARK)
    axa.set_xlim(0, amax * 2.9); axa.set_ylim(YL)
    axa.set_yticks([ypos[k] for k in K]); axa.set_yticklabels([WRAP.get(k, TYPES[k][2]) for k in K], fontsize=FT, linespacing=0.95)
    axa.set_xticks([]); axa.spines[['top', 'right', 'bottom']].set_visible(False)
    axa.text(0, YL[1] + 0.02, 'Mean %\nwithin 50 µm', ha='left', va='bottom', fontsize=FSm, color=DARK, linespacing=1.0)
    # forest
    rng = np.random.default_rng(0)
    axf.axvline(0, color='#777777', lw=1.1, ls='--', zorder=1)
    for k in K:
        v = RES[k]; y = ypos[k]
        axf.scatter(SEC[k].values, y + rng.uniform(-0.075, 0.075, len(SEC)), s=13, color=[SCOL[x] for x in SEC.index], alpha=0.5, lw=0, zorder=1.5)   # section colours as in the right panel
        axf.plot([v['lo'], v['hi']], [y, y], color='white', lw=3.4, zorder=1.9, solid_capstyle='butt')   # thin white edge
        axf.plot(v['estimate'], y, 'o', ms=8.5 + 2.6, mfc='white', mec='white', mew=0, zorder=1.95)   # halo under the black CI line, above section dots
        axf.plot([v['lo'], v['hi']], [y, y], color='black', lw=1.8, zorder=2, solid_capstyle='butt')
        axf.plot(v['estimate'], y, 'o', ms=8.5, mec='black', mew=1.4, mfc='black' if v['q'] < 0.05 else 'white', zorder=3)
    lo_all = min(SEC.values.min(), min(RES[k]['lo'] for k in K)); hi_all = max(SEC.values.max(), max(RES[k]['hi'] for k in K))
    pad = 0.05 * (hi_all - lo_all); axf.set_xlim(lo_all - pad, hi_all + pad); axf.set_ylim(YL)
    axf.set_xticks([-0.2, 0, 0.2]); axf.set_xticklabels(['−0.2', '0', '0.2'])
    axf.set_yticks([]); axf.spines[['top', 'right', 'left']].set_visible(False); axf.tick_params(labelsize=FT)
    axf.set_xlabel('Effect on Treg\nSuppression Score\n(SD per SD of\nlocal fraction)', fontsize=FL, linespacing=1.0)
    # right: DC slope by OGN+RSPO3+ tertile
    TT = ['bottom', 'middle', 'top']; XL = ['Low', 'Mid', 'High']
    axr.axhline(0, color='#777777', lw=1.1, ls='--', zorder=1)
    OFF = dict(zip(SECS, np.linspace(-0.13, 0.13, len(SECS))))
    for s_ in SECS:   # thin per-section lines (one color per section)
        xs = [i + OFF[s_] for i in range(3)]; ys = [WIDE.loc[s_, tt] for tt in TT]
        axr.plot(xs, ys, '-', color=SCOL[s_], lw=1.0, alpha=0.4, zorder=1.5)
        axr.scatter(xs, ys, s=14, color=SCOL[s_], lw=0, alpha=0.5, zorder=1.6)
    for i, tt in enumerate(TT):   # pooled estimates on top
        v = TER[tt]
        axr.plot([i, i], [v['lo'], v['hi']], color='white', lw=3.4, zorder=1.9, solid_capstyle='butt')   # thin white edge
        axr.plot(i, v['estimate'], 'o', ms=8.5 + 2.6, mfc='white', mec='white', mew=0, zorder=1.95)   # halo under the black CI line, above section dots
        axr.plot([i, i], [v['lo'], v['hi']], color='black', lw=1.8, zorder=2, solid_capstyle='butt')
        axr.plot(i, v['estimate'], 'o', ms=8.5, mec='black', mew=1.4, mfc='black' if v['p'] < 0.05 else 'white', zorder=3)
    axr.set_xlim(-0.4, 2.4); axr.set_xticks(range(3)); axr.set_xticklabels(XL, fontsize=FT)
    for i, tt in enumerate(TT):   # n per tertile, own row under the tick labels
        axr.annotate(f'{TER[tt]["n"]:,}', (i, 0), xycoords=('data', 'axes fraction'), xytext=(0, -27), textcoords='offset points',
                     ha='center', va='top', fontsize=FSm, color=DARK)
    axr.annotate('n =', (-0.4, 0), xycoords=('data', 'axes fraction'), xytext=(-4, -27), textcoords='offset points',
                 ha='right', va='top', fontsize=FSm, color=DARK)
    # DC % row: pooled mean local DC fraction (% of cells within 50 um) per tertile, read from fig4d_dc_by_ogn_tertile_checks.csv (F5D_dc_by_ogn_tertile)
    CHK = pd.read_csv(os.path.dirname(OUT) + '/fig4d_dc_by_ogn_tertile_checks.csv')
    CHK = CHK[CHK.block == '1_desc_pooled'].set_index('tertile')
    DCP = {tt: float(CHK.loc[lab, 'dcpct_mean']) for tt, lab in zip(TT, XL)}
    assert all(int(CHK.loc[lab, 'n']) == TER[tt]['n'] for tt, lab in zip(TT, XL)), 'tertile n differ from fig4d_dc_by_ogn_tertile_checks.csv'
    print('DC % row (from checks CSV):', {XL[i]: DCP[tt] for i, tt in enumerate(TT)}, flush=True)
    DY2 = 27 + FSm * 1.25 + 2   # second row, one n-row height below
    for i, tt in enumerate(TT):
        axr.annotate(f'{DCP[tt]:.1f}', (i, 0), xycoords=('data', 'axes fraction'), xytext=(0, -DY2), textcoords='offset points',
                     ha='center', va='top', fontsize=FSm, color=DARK)
    axr.annotate('DC %', (-0.4, 0), xycoords=('data', 'axes fraction'), xytext=(-4, -DY2), textcoords='offset points',
                 ha='right', va='top', fontsize=FSm, color=DARK)
    YLO, YHI = YLIM_R
    assert PS.slope.min() > YLO and PS.slope.max() < YHI and min(TER[t]['lo'] for t in TT) > YLO and max(TER[t]['hi'] for t in TT) < YHI
    axr.set_ylim(YLO, YHI); axr.set_yticks(YTK_R); axr.set_yticklabels([('−' if t < 0 else '') + f'{abs(t):g}' for t in YTK_R])
    axr.spines[['top', 'right']].set_visible(False); axr.tick_params(axis='y', labelsize=FT)
    axr.set_ylabel('DC effect on Treg\nSuppression Score\n(SD per SD)', fontsize=FL, linespacing=1.0)
    axr.set_xlabel('Local OGN⁺RSPO3⁺\nfibroblasts\n(within-section tertile)', fontsize=FL, linespacing=1.0, ha='right', va='top', multialignment='center')
    axr.xaxis.set_label_coords(1.0, -(DY2 + FSm * 1.25 + 4) / 72 / AH)   # right-aligned under the DC % row
    IXT = axr.annotate(f'DC × OGN⁺RSPO3⁺ interaction\nP = {IX["p"]:.3g}; {IX["agree"]} sections', (1.0, 1.0), xycoords='axes fraction',
                       xytext=(0, 3), textcoords='offset points', ha='right', va='bottom',
                   fontsize=FSm, color='black', linespacing=1.05)
    KEY = axr.text(0.02, KEY_Y, 'Thin lines: sections\n(one color each)\nBold: pooled ± 95% CI', transform=axr.transAxes, ha='left', va='top',
                   fontsize=FKEY, color=DARK, linespacing=1.05, zorder=4)
    # overlap check: every text box against every other, and against the data area of axes it does not belong to
    fig.canvas.draw(); rr = fig.canvas.get_renderer()
    TX = [t for a in (axa, axf, axr) for t in a.get_xticklabels() + a.get_yticklabels() + [a.xaxis.label, a.yaxis.label] + a.texts
          if t.get_visible() and t.get_text()]
    BB = [t.get_window_extent(rr) for t in TX]
    for i in range(len(TX)):
        for j in range(i + 1, len(TX)):
            assert not BB[i].overlaps(BB[j]), ('text overlap', TX[i].get_text(), TX[j].get_text())
    for t, b in zip(TX, BB):
        for a in (axa, axf, axr):
            if a is not t.axes: assert not b.overlaps(a.get_window_extent(rr)), ('text inside another axes', t.get_text())
    # interaction text must sit above every tertile CI
    assert IXT.get_window_extent(rr).y0 > axr.get_window_extent(rr).y1, 'interaction text inside the axes'
    FB = fig.bbox
    for t, b in zip(TX, BB): assert b.x0 >= FB.x0 and b.x1 <= FB.x1 and b.y0 >= FB.y0 and b.y1 <= FB.y1, ('text outside figure', t.get_text(), b, FB)
    print('CHECK layout: no text overlaps', len(TX), 'text boxes', flush=True)
    with open(OUT + '.png', 'wb') as fh: fig.savefig(fh, format='png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig); print('wrote', OUT + '.png', flush=True)
    # ---------------- legend ----------------
    L = (f'Local APC types and {SS} in CD resections (CosMx 6K; {n_s} sections, 4 patients; {n_t:,} Tregs). '
     f'Local fraction = each cell type\'s share of all cells within 50 µm of a Treg, z-scored within section; strip, mean %. '
     f'Left, one regression of the {SS} on all seven APC types (section fixed effects, library-size covariate): dot, estimate; bar, 95% CI '
     f'(section-clustered SEs, t {n_s - 1} df); filled, BH q < 0.05; small dots, sections (colors as in right panel). '
     f'Right, DC slope within within-section tertiles of local OGN⁺RSPO3⁺ fibroblast fraction (same SEs; filled, P < 0.05); '
     f'thin lines, per-section slopes, one color per section; DC %, mean local DC fraction within 50 µm; '
     f'DC × OGN⁺RSPO3⁺ interaction {IX["estimate"]:+.3f}, clustered P = {IX["p"]:.3g}, same sign in {IX["agree"]} sections.')
    nw = len(L.split()); print('legend words', nw); assert nw <= 130   # legend word limit
    with open(OUT + '_legend.txt', 'w') as fh: fh.write(L + '\n')
    print(L)


# ======================================================================================================================
# Section S9D_celltype_forest: forest plot of every local cell type vs the Treg Suppression Score (separate LMMs)
# Panels: Supp 9D
# ======================================================================================================================
if _run('S9D_celltype_forest'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/105_all_celltypes_forest.py'   # working folder of this step (relative paths resolve from here)
    """Supp 9D: every local cell type vs the Treg Suppression Score (CD resections, CosMx 6K; 11 sections, 4 patients, 6,554 strict Tregs), drawn from
the table of spatial_results_table; nothing is refitted.
Input: updated_figures/cd_spatial_results_all.csv, block 'local cell type' (33 rows; see cd_spatial_results_all_README.txt).
  Model per type (separate): act_z ~ type_z + cz + (1|section); type_z = local fraction of the type among all cells within 50 um of the Treg,
  z-scored within section; cz = log1p(Treg library size) z-scored within section. effect_sd with 95% Wald CI; q = BH across the 33 types;
  clustered_p = OLS with section-clustered SE (10 df); sections_same_sign = per-section OLS sign agreement.
Display: SELENOP+ Fib is not shown (SELENOP+ Mac is shown), so 32 rows. 'Absoprtive Epi' is the annotation's spelling for BEST4+ cells (a different
  type from 'Absorptive Epi' = absorptive progenitors; separate rows of the table), so it is not merged and is displayed as 'BEST4+ Epi'.
  The DC row (label, dot, CI) is red.
Filled marker: q < 0.05 AND section-clustered P < 0.05 AND >= 9/11 sections with the same sign; open otherwise; the criteria are written in an
  on-plot key. Figure size is fixed at 6.0 x 10.0 in (saved without bbox 'tight').
Outputs (PNG only, no CSV): SuppFig9A-F_niche_cxcl12_apc/supp_all_celltypes_treg_score_forest.png + _legend.txt"""
    import os, sys
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    import numpy as np, pandas as pd
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    IN = '/share/fsmresfiles/UC/revision_jenny/updated_figures/cd_spatial_results_all.csv'
    OUT = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc'
    PNG = f'{OUT}/supp_all_celltypes_treg_score_forest.png'; LEG = f'{OUT}/supp_all_celltypes_treg_score_forest_legend.txt'

    DISP = {'OGN+RSPO3+ Fib': 'OGN⁺RSPO3⁺ Fib', 'ADAMDEC1+ Fib': 'ADAMDEC1⁺ Fib', 'PLA2G2A+ ECM Fib': 'PLA2G2A⁺ ECM Fib',
            'CD9+ M2 Mac': 'CD9⁺ M2 Mac', 'Absoprtive Epi': 'BEST4⁺ Epi', 'Gd T': 'γδ T', 'SELENOP+ Mac': 'SELENOP⁺ Mac'}
    EXCL = ['SELENOP+ Fib']                                   # not shown (SELENOP+ Mac is shown)

    T = pd.read_csv(IN)
    T = T[T.block == 'local cell type'].copy()
    assert len(T) == 33 and T.feature.is_unique and (T.n_tregs == 6554).all()
    assert {'Absoprtive Epi', 'Absorptive Epi'} <= set(T.feature)       # two separate rows -> different types (BEST4+ vs progenitor); not merged
    T = T[~T.feature.isin(EXCL)].copy(); assert len(T) == 32
    T['n_same'] = T.sections_same_sign.str.split('/').str[0].astype(int)
    T['n_sec'] = T.sections_same_sign.str.split('/').str[1].astype(int); assert (T.n_sec == 11).all()
    T['filled'] = (T.q < 0.05) & (T.clustered_p < 0.05) & (T.n_same >= 9)
    T['label'] = T.feature.map(lambda s: DISP.get(s, s))
    assert not T.label.str.contains('Absoprtive').any() and T.label.is_unique
    T = T.sort_values('effect_sd', ascending=True).reset_index(drop=True)   # most negative at bottom
    y = np.arange(len(T))

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13, 'axes.linewidth': 1.0})
    FS_LAB, FS_SEC, FS_KEY = 15, 13, 13                 # labels/ticks/axis label 15 pt; section text and key 13 pt
    W, H = 6.0, 10.0                                      # PNG size in inches; saved without bbox 'tight' so the size is exact
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.4, 0.2, 0.4, 0.75])              # placeholder; set from measured text extents below
    ax.axvline(0, color='0.45', ls='--', lw=1.1, zorder=0)
    RED = '#c0392b'                                       # DC row (label, dot, CI) in red
    col = np.where(T.feature.values == 'DC', RED, 'black'); assert (col == RED).sum() == 1 and T.filled[T.feature == 'DC'].all()
    ax.hlines(y, T.ci_low, T.ci_high, colors=col, lw=1.6, zorder=1)
    r = T.filled.values
    ax.scatter(T.effect_sd[r], y[r], s=70, facecolor=col[r], edgecolor=col[r], lw=1.3, zorder=3)
    ax.scatter(T.effect_sd[~r], y[~r], s=70, facecolor='white', edgecolor='black', lw=1.3, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(T.label, fontsize=FS_LAB)
    for t in ax.get_yticklabels():                       # 'DC' row label in red, regular weight
        if t.get_text() == 'DC': t.set_color(RED)
    assert sum(t.get_color() == '#c0392b' for t in ax.get_yticklabels()) == 1 and all(t.get_fontweight() == 'normal' for t in ax.get_yticklabels())
    ax.set_ylim(-0.7, len(T) - 0.3)
    xl = (-0.125, 0.135); ax.set_xlim(*xl)
    assert T.ci_low.min() > xl[0] and T.ci_high.max() < xl[1]
    ax.set_xticks([-0.1, 0, 0.1]); ax.set_xticklabels(['−0.1', '0', '0.1'], fontsize=FS_LAB)
    ax.set_xticks([-0.05, 0.05], minor=True)
    ax.tick_params(axis='y', length=0, pad=3)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    XLAB = 'Effect on Treg Suppression Score\n(SD per SD of local fraction within 50 µm)'   # two lines (one line is ~7.5 in at 15 pt: does not fit 6 in)
    ax.set_xlabel(XLAB, fontsize=FS_LAB, labelpad=4)
    # section agreement as small text at the right (axes-fraction x, data y)
    tr = ax.get_yaxis_transform()
    sec = [ax.text(1.02, yi, t, transform=tr, ha='left', va='center', fontsize=FS_SEC, color='0.25', clip_on=False)
           for yi, t in zip(y, T.sections_same_sign)]
    sec.append(ax.text(1.02, len(T) - 0.1, 'Same\nsign', transform=tr, ha='left', va='bottom', fontsize=FS_SEC, color='0.25', clip_on=False))
    # criteria key under the x-axis label (no data there)
    h = [Line2D([], [], marker='o', ls='', mfc='black', mec='black', ms=8.5,
                label='q < 0.05\nSection-clustered P < 0.05\n≥ 9/11 sections same sign'),
         Line2D([], [], marker='o', ls='', mfc='white', mec='black', ms=8.5, label='Otherwise')]
    lg = fig.legend(handles=h, loc='lower center', bbox_to_anchor=(0.5, 0.004), ncol=2, frameon=False, fontsize=FS_KEY,
                    handletextpad=0.2, columnspacing=1.0, borderpad=0.1, labelspacing=0.4)
    # layout from measured extents (inches): left = widest y label, right = section column, bottom = x label + key
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); dpi = fig.dpi
    ylab_w = max(t.get_window_extent(rd).width for t in ax.get_yticklabels()) / dpi
    sec_w = max(t.get_window_extent(rd).width for t in sec) / dpi
    lg_h = lg.get_window_extent(rd).height / dpi
    ax_bot_in = lg_h + 0.06 + 0.12 + (XLAB.count('\n') + 1) * FS_LAB / 72 * 1.2 + FS_LAB / 72 * 1.3 + 0.12   # key + gap + 3-line x label + tick labels + pads
    left_in = 0.04 + ylab_w + 0.06; right_in = 0.03 + sec_w + 0.06; top_in = 0.12 + sec[-1].get_window_extent(rd).height / dpi   # room for the 2-line section header
    ax.set_position([left_in / W, ax_bot_in / H, 1 - (left_in + right_in) / W, 1 - (ax_bot_in + top_in) / H])
    fig.canvas.draw()
    # the 3-line x label is wider than the axes: keep it centred on the axes unless it would leave the figure, then shift it just enough
    xb = ax.xaxis.label.get_window_extent(rd); axb = ax.get_window_extent(rd); over = xb.x1 - (fig.bbox.x1 - 0.03 * dpi)
    if over > 0:
        ax.xaxis.label.set_x(0.5 - over / axb.width)          # label x is in axes fraction; y stays automatic
    fig.canvas.draw()
    # checks: every text inside the figure; no two text boxes overlap; key covers no data and sits below the x label
    F = fig.bbox
    arts = list(ax.get_yticklabels()) + [t for t in ax.get_xticklabels() if t.get_text()] + [ax.xaxis.label] + sec + [lg]
    boxes = [a.get_window_extent(rd) for a in arts]
    for a_, bx in zip(arts, boxes):
        assert bx.x0 >= F.x0 - 0.5 and bx.x1 <= F.x1 + 0.5 and bx.y0 >= F.y0 - 0.5 and bx.y1 <= F.y1 + 0.5, ('outside figure', a_, bx)
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            ov = (min(boxes[i].x1, boxes[j].x1) - max(boxes[i].x0, boxes[j].x0) > 0.5) and (min(boxes[i].y1, boxes[j].y1) - max(boxes[i].y0, boxes[j].y0) > 0.5)
            assert not ov, ('text overlap', arts[i], arts[j])
    assert boxes[-1].y1 < ax.xaxis.label.get_window_extent(rd).y0 and boxes[-1].y1 < ax.get_window_extent(rd).y0
    from PIL import Image
    print(f'figure {W} x {H} in; axes width {ax.get_position().width * W:.2f} in; y-label column {ylab_w:.2f} in; section column {sec_w:.2f} in')
    fig.savefig(PNG, dpi=300)
    wpx, hpx = Image.open(PNG).size; assert (wpx, hpx) == (round(W * 300), round(H * 300)), (wpx, hpx)
    print(f'saved PNG {wpx} x {hpx} px = {wpx / 300:.2f} x {hpx / 300:.2f} in')
    plt.close(fig)

    nr = int(T.filled.sum())
    leg = ('Local cell types and the Treg Suppression Score in Crohn\'s disease resections (CosMx 6K; 6,554 Tregs, 11 sections, 4 patients). '
       'Each cell type was tested separately: linear mixed model, score ~ local fraction of the type among all cells within 50 µm of the Treg '
       '+ Treg library size + (1 | section); score, fraction and library size z-scored within section. Dots, estimate; bars, 95% CI. '
       'Filled, q < 0.05 (BH, 33 types) and section-clustered P < 0.05 and same sign in ≥ 9/11 sections; open, otherwise. '
       'Right, sections with the same sign. SELENOP⁺ fibroblasts not shown.')
    nw = len(leg.split()); assert nw <= 100, nw
    for bad in ('quiescent', 'colocaliz', 'SELENOP⁺ Fib', 'SELENOP+ Fib', 'territor', 'blunt', 'ctivation score'):
        assert bad not in leg
    with open(LEG, 'w') as fh: fh.write(leg + '\n')
    print(f'rows {len(T)}; filled {nr}: {list(T.feature[T.filled])}; legend words {nw}')
    print(T[['label', 'effect_sd', 'ci_low', 'ci_high', 'q', 'clustered_p', 'sections_same_sign', 'filled']].iloc[::-1].to_string())
