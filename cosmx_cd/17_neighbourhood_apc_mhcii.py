#!/usr/bin/env python
"""CD CosMx 6K resections: antigen-presenting cells (APCs) and MHC-II per cellular neighbourhood.

Panels made
  Supp 9C  APC percentage, adjusted per-APC MHC-II score and inside-minus-outside Treg Suppression Score per CN

Inputs
  /share/fsmresfiles/UC/AtoMx/CD_surgical_resection/combined/ (all_samples.parquet, raw_counts_strict_tregs.csv) and <sample>/Processed/norm_log1p.h5ad
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/
    B_resection_analyses/lr_modules_cellchat/ (treg_receiver_module_scores_per_treg.csv.gz, treg_receiver_pairs_kept.csv: 13_neighbourhood_signalling.py,
      section treg_receiver_lr_modules)
    B_resection_analyses/cellular_neighborhoods_v2/ and ligand_receptor_state/ (11_neighbourhoods_local_tables.py)
    D_shared_inputs/treg_definition_comparison_per_cell.csv (10_treg_suppression_score.py)
    A_fig4_panels/fig4_cn_treg_state_forest.csv and fig4_cn_composition_heatmap.csv (12_neighbourhood_treg_score.py)
    archive/supp_apc_abundance_by_neighbourhood.csv (APC share of cells per CN: 16_cxcl12_regions.py, section S9E_cxcl12_by_celltype)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/ (APC tables; the per-APC table of cn10_apc_abundance is written
  to archive/ext_cn10_antigen_presentation_abundance.csv) and /share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/ (Supp 9C panel)

Run order of sections
  cn10_apc_abundance -> cn_apc_mhcii_table -> cn10_apc_matched_score -> S9C_cn_apc_mhcii
  cn_apc_mhcii_table and cn10_apc_matched_score also read the APC table of section S9E_cxcl12_by_celltype of 16_cxcl12_regions.py
  (run first with RUN_FIG2=1 RUN_SECTIONS=S9E_cxcl12_by_celltype).

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
# Section cn10_apc_abundance: APCs around Tregs and per-APC MHC-II expression by cellular neighbourhood
# Panels: Supp 9C inputs (per-APC MHC-II table, per-Treg APC counts)
# ======================================================================================================================
if _run('cn10_apc_abundance'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/78_cn10_antigen_presentation_abundance.py'   # working folder of this step (relative paths resolve from here)
    """APCs around Tregs and per-APC MHC-II expression by cellular neighbourhood (CD resections, CosMx 6K; inputs of Supp 9C).
Definitions follow section treg_receiver_lr_modules of 13_neighbourhood_signalling.py: Antigen Presentation = NicheNet pairs grouped by the CellChatDB
MHC-II pathway; ligands on the 6K panel = HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPB1, HLA-DQA1, HLA-DQA2, HLA-DQB1, HLA-DRA (-> CD4).
APCs = cell_type_general DC, B cell, CD9+ M2 Mac, Mo-Mac and SELENOP+ Mac (SELENOP+ Fib is not an APC), reported in three groups: DC, B cell and
Macrophage (= SELENOP+ Mac + CD9+ M2 Mac + Mo-Mac), also for the per-group models and the group adjustment. The per-Treg table of
treg_receiver_lr_modules (treg_receiver_module_scores_per_treg.csv.gz; DC, B cell, CD9+ M2 Mac and Mo-Mac as senders) supplies the CN, library size
and the per-ligand ligand supply; APC counts and ligand supply are recomputed here with the five APC types, and the four-type recomputation is
asserted to agree with the per-Treg table. Receivers = the 6,554 strict Tregs of 11 CD resection sections (4 patients); ligand z = within-section z
over all cells of the section (norm_log1p.h5ad); 50 um radius (px x 0.12028); neighbourhood = CN_k10 of cellular_neighborhoods_v2/cn_per_cell.parquet
(1-based); CN6 / CN9 omitted (too few Tregs).
Analyses (Antigen Presentation only):
 1. APC abundance: per Treg, number of APCs within 50 um; median, IQR, % >= 1/3/5/10 per CN and per CN x section; APC group composition near Tregs
    per CN; unique APCs near CN10 Tregs and how concentrated the Treg-APC links are.
 2. Per-APC MHC-II score (APC level, not via Tregs): APC score = mean within-section z of the 9 HLA ligands (the quantity averaged over the APCs near
    a Treg in the ligand supply). y = APC score z-standardised within section over APCs. Grouping (A) the APC's own CN; (B) the CN of the APC's
    nearest strict Treg within 50 um (APCs that feed Tregs; each APC once); (C) binary: APC within 50 um of >= 1 CN10 Treg vs APC within 50 um of
    Tregs of other CNs only.
    Model: y ~ C(CN, ref = CN10) + lib_z + C(APC group) + (1 | section) (MixedLM, REML), lib_z = within-section z of log1p total counts; per APC
    group the same model without the group term. Contrast reported = CN - CN10 (negative = CN10 APCs higher); 95% Wald CI; BH over the 7 contrasts
    within each grouping x APC set. Section dots = per section, mean residual of OLS y ~ C(group) + lib_z in the CN minus that in CN10 (sections with
    >= 10 APCs in both). MHC-II-high = top quartile of the APC score within section over all APCs.
 3. Threshold sensitivity: Treg-level ligand-supply model (y = within-section z of the Treg's ligand supply L, y ~ in_CN + lib_z + (1 | section);
    sections with >= 10 scored Tregs in and outside the CN; BH over CNs) refitted with L computed only for Tregs with >= 1 / 3 / 5 / 10 APCs within 50 um.
 4. Leave-one-section-out of the CN10 fit (threshold 1).
 5. Abundance-weighted supply: S = SUM over the APCs within 50 um of their MHC-II score (= L x n; 0 for Tregs with no APC), all Tregs; same model and
    BH as 3.
Outputs: panels/fig4_and_related/archive/ext_cn10_antigen_presentation_abundance.{csv,png} + _legend.txt (figure: panels a, b; analyses 3-5 in the
CSV and one legend sentence) and B_resection_analyses/lr_modules_cellchat/cn10_apc_per_treg_78.csv.gz (per-Treg APC counts by APC group). Both
tables are read by cn_apc_mhcii_table and cn10_apc_matched_score.
Environment: APC_CACHE=<path.pkl> stores / re-uses the per-APC and per-Treg tables of the compute stage."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
    from anndata.io import read_elem
    from scipy import sparse, stats as sst
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    from statsmodels.tools.sm_exceptions import ConvergenceWarning; warnings.simplefilter('ignore', ConvergenceWarning)
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    STEM = f'{P4}/archive/ext_cn10_antigen_presentation_abundance'; PERTREG = f'{OUTB}/cn10_apc_per_treg_78.csv.gz'; CACHE74 = f'{OUTB}/treg_receiver_module_scores_per_treg.csv.gz'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50; MIN_TREG_CN = 10; MIN_APC_SEC = 10
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    APC = ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac']   # SELENOP+ Mac is an APC
    APC74 = ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac']                 # APC set of the per-Treg table (consistency check only)
    APC_GRP = {'DC': 'DC', 'B cell': 'B cell', 'CD9+ M2 Mac': 'Macrophage', 'Mo-Mac': 'Macrophage', 'SELENOP+ Mac': 'Macrophage'}
    APCG = ['DC', 'B cell', 'Macrophage']                             # all macrophages in one group
    HLA = ['HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DRA']
    M = 'Antigen presentation'; EXCL = [6, 9]; CNS = [1, 2, 3, 4, 5, 7, 8, 10]; OTH = [c for c in CNS if c != 10]; THRS = [1, 3, 5, 10]
    FIGCN = {1: 'CN1: Epithelial', 2: 'CN2: Plasma', 3: 'CN3: Myofib / Neuron', 4: 'CN4: DC / Mac', 5: 'CN5: InfMo-Mac',
             7: 'CN7: Activ Fib / Neut', 8: 'CN8: B / Tfh', 10: 'CN10: OGN⁺RSPO3⁺ Fib'}
    # ---------------- per-Treg table of treg_receiver_lr_modules and the module's ligand set
    W = pd.read_csv(CACHE74, index_col=0); assert len(W) == 6554
    KP = pd.read_csv(f'{OUTB}/treg_receiver_pairs_kept.csv'); KP = KP[(KP.module == M) & KP.kept]
    assert sorted(KP.ligand) == sorted(HLA) and set(KP.receptor) == {'CD4'} and set(KP.sender_set) == {'MHCII_APC'}
    LC = [f'L|MHCII_APC|{g}' for g in HLA]; assert all(c in W for c in LC)
    zs = lambda v: (v - v.mean()) / v.std()
    W['lib_z'] = W.groupby('section').log_total_counts.transform(zs)
    W['n74'] = W['n|MHCII_APC'].astype(int)   # four-type APC count of the per-Treg table; the analyses use the five-type count below
    APCC = os.environ.get('APC_CACHE')
    if APCC and os.path.exists(APCC):
        A, TT = pd.read_pickle(APCC)
    else:
        sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', 'total_counts', '__index_level_0__']).to_pandas()
        sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
        sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
        assert set(APC) <= set(sa.ct.unique())
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
        ab, tb = [], []
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
                              'log_total_counts': np.log1p(cs.total_counts.values[am].astype(float)),
                              'mhc2': np.nanmean(Z[am], 1)}, index=idx[am])
            a['CN'] = cn.loc[a.index, 'CN_k10'].astype(int).values; assert (cn.loc[a.index, 'section'].astype(str).values == s_).all()
            tg = W[W.section == s_]; tx = sa.loc[tg.index, ['x', 'y']].values
            # Treg -> APCs within 50 um (consistency with n and L of the per-Treg table; per-type composition)
            nb = cKDTree(a[['x', 'y']].values).query_ball_point(tx, r=RAD)
            t = pd.DataFrame(index=tg.index)
            t['n_chk'] = [len(i) for i in nb]
            t['L_chk'] = [a.mhc2.values[i].mean() if len(i) else np.nan for i in nb]
            o74 = np.isin(a.apc_type.values, APC74)
            t['n_chk74'] = [int(o74[i].sum()) for i in nb]
            t['L_chk74'] = [a.mhc2.values[i][o74[i]].mean() if o74[i].sum() else np.nan for i in nb]
            t['sum_mhc2'] = [a.mhc2.values[i].sum() if len(i) else 0.0 for i in nb]
            for ty in APC: t[f'n_{ty}'] = [int((a.apc_type.values[i] == ty).sum()) for i in nb]
            tb.append(t)
            # APC -> Tregs within 50 um: nearest Treg's CN, any CN10 Treg, number of (CN10) Tregs fed
            tcn = tg.CN.values; rev = cKDTree(tx).query_ball_point(a[['x', 'y']].values, r=RAD)
            dist, nn = cKDTree(tx).query(a[['x', 'y']].values, k=1, distance_upper_bound=RAD)
            a['near_treg'] = np.isfinite(dist); a['nearest_treg_CN'] = np.where(np.isfinite(dist), tcn[np.minimum(nn, len(tcn) - 1)], -1)
            a['n_tregs_50um'] = [len(i) for i in rev]; a['n_cn10_tregs_50um'] = [int((tcn[i] == 10).sum()) if len(i) else 0 for i in rev]
            ab.append(a); print(s_, 'APCs', len(a), 'Tregs', len(tg), flush=True)
            del E, Z; gc.collect()
        A = pd.concat(ab); TT = pd.concat(tb).loc[W.index]; del sa, cn; gc.collect()
        if APCC: pd.to_pickle((A, TT), APCC)
    # ---------------- consistency with the per-Treg table of treg_receiver_lr_modules
    SAME74 = bool((TT.n_chk.values == W.n74.values).all())   # True only if the per-Treg table was built with the five APC types
    nv, lv = (TT.n_chk, TT.L_chk) if SAME74 else (TT.n_chk74, TT.L_chk74)
    print('per-Treg table APC set == five types:', SAME74, flush=True)
    assert (nv.values == W.n74.values).all(), 'APC counts differ from the per-Treg table'
    ok = W[f'{M}|L'].notna()
    dl = np.abs(lv[ok] - W.loc[ok, f'{M}|L']); print('max |L recomputed - per-Treg table L| =', dl.max(), flush=True)
    assert dl.max() < 5e-3 * max(1, W.loc[ok, f'{M}|L'].abs().max()), 'per-Treg L differs from the per-Treg table'
    assert (nv.values > 0).sum() == ok.sum()
    W['nAPC'] = TT.n_chk.values.astype(int); W['L_apc'] = TT.L_chk.values      # five APC types
    print('Tregs with >= 1 APC within 50 um: four types', int((TT.n_chk74 > 0).sum()), '| five types', int((W.nAPC > 0).sum()), flush=True)
    W['S_sum'] = TT.sum_mhc2.values
    for g_ in APCG: W[f'n_{g_}'] = TT[[f'n_{ty}' for ty in APC if APC_GRP[ty] == g_]].sum(1).values
    assert (W[[f'n_{g_}' for g_ in APCG]].sum(1) == W.nAPC).all()
    with open(PERTREG, 'wb') as fh: W[['section', 'CN', 'nAPC', 'nb_n'] + [f'n_{g_}' for g_ in APCG]].to_csv(fh, compression='gzip')
    print('wrote', PERTREG, flush=True)
    A['cell_type'] = A.apc_type; A['apc_type'] = A.cell_type.map(APC_GRP); assert A.apc_type.notna().all()
    A['patient'] = A.section.str.split('_').str[0]
    A['lib_z'] = A.groupby('section').log_total_counts.transform(zs)
    A['y'] = A.groupby('section').mhc2.transform(zs)
    A['high'] = A.mhc2 >= A.groupby('section').mhc2.transform(lambda v: v.quantile(0.75))
    print('APCs:', len(A), A.apc_type.value_counts().to_dict(), flush=True)
    sig3 = lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    ROWS = []
    # ================= 1. APC abundance near Tregs
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
    # unique APCs feeding CN10 Tregs and concentration of the Treg-APC links
    f10 = A[A.n_cn10_tregs_50um > 0]
    links = np.sort(f10.n_cn10_tregs_50um.values)[::-1]; top10 = links[:max(1, int(np.ceil(0.1 * len(links))))].sum() / links.sum()
    CONC = dict(analysis='1_apc_abundance', level='unique APCs within 50 um of >= 1 CN10 Treg', CN=10, CN_label=FIGCN[10], section='all',
                n_unique_apcs=len(f10), n_sections_with_such_apcs=f10.section.nunique(), median_unique_apcs_per_section=f10.groupby('section').size().median(),
                min_unique_apcs_per_section=f10.groupby('section').size().min(), median_cn10_tregs_per_apc=np.median(links),
                pct_links_from_top10pct_apcs=100 * top10, **{f'pct_{t}_of_nearby_apcs': 100 * (f10.apc_type == t).mean() for t in APCG})
    ROWS.append(CONC)
    print(pd.DataFrame([r for r in ROWS if r['level'] == 'CN'])[['CN', 'n_tregs', 'median_apcs_50um', 'q1_apcs_50um', 'q3_apcs_50um', 'pct_ge1', 'pct_ge3', 'pct_ge5', 'pct_ge10'] + [f'pct_{t}_of_nearby_apcs' for t in APCG]].round(1).to_string(index=False))
    print(CONC, flush=True)
    # ================= 2. per-APC MHC-II score
    def fit(formula, d, term_prefix='C(grp, Treatment(10))'):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm']:
            try:
                f = smf.mixedlm(formula, d, groups=d.section).fit(reml=True, method=meth)
                ks = [k for k in f.params.index if k.startswith(term_prefix)]
                if all(np.isfinite(f.bse[k]) and f.bse[k] > 0 for k in ks): return f, meth
            except Exception: pass
        return None, 'failed'
    GROUPINGS = {'own_CN': 'CN', 'nearest_treg_CN': 'nearest_treg_CN'}
    APC2 = {}
    for gname, gcol in GROUPINGS.items():
        d0 = A[A[gcol].isin(CNS)].copy(); d0['grp'] = d0[gcol].astype(int)
        for tset in ['all APC types (type-adjusted)'] + APCG:
            d = d0 if tset.startswith('all') else d0[d0.apc_type == tset]
            frm = 'y ~ C(grp, Treatment(10)) + lib_z' + (' + C(apc_type)' if tset.startswith('all') else '')
            f, meth = fit(frm, d)
            # section dots: OLS residual (type + lib) per section, mean in CN minus mean in CN10
            rf = 'y ~ lib_z' + (' + C(apc_type)' if tset.startswith('all') and True else '')
            res = []
            for s_, gs in d.groupby('section'):
                rr = gs.y - smf.ols(rf if gs.apc_type.nunique() > 1 or not tset.startswith('all') else 'y ~ lib_z', gs).fit().fittedvalues
                res.append(gs[['grp']].assign(r=rr, section=s_))
            res = pd.concat(res)
            out = []
            for c in OTH:
                k = f'C(grp, Treatment(10))[T.{c}]'
                dsec = []
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
            o['row_type'] = np.where(mm, 'mixed_model', 'section')
            APC2[(gname, tset)] = o; ROWS += o.to_dict('records')
            print(gname, '|', tset, '|', meth); print(o[mm][['CN', 'estimate', 'ci_low', 'ci_high', 'q', 'n_apcs_in_CN', 'n_apcs_CN10', 'sections_CN10_higher']].round(3).to_string(index=False), flush=True)
    # descriptive per-APC table: n APCs, mean score, % MHC-II-high per CN (own / nearest-Treg), per type, per section
    for gname, gcol in GROUPINGS.items():
        d0 = A[A[gcol].isin(CNS)]
        for c in CNS:
            g = d0[d0[gcol] == c]
            r = dict(analysis='2_per_apc_descriptive', grouping=gname, apc_types='all', CN=c, CN_label=FIGCN[c], section='all', n_apcs=len(g),
                     mean_mhc2_score=g.mhc2.mean(), pct_mhc2_high_top_quartile=100 * g.high.mean(), n_sections_ge10_apcs=int((g.groupby('section').size() >= MIN_APC_SEC).sum()),
                     **{f'pct_{t}': 100 * (g.apc_type == t).mean() for t in APCG})
            ROWS.append(r)
            for t in APCG:
                gt = g[g.apc_type == t]
                ROWS.append(dict(analysis='2_per_apc_descriptive', grouping=gname, apc_types=t, CN=c, CN_label=FIGCN[c], section='all', n_apcs=len(gt),
                                 mean_mhc2_score=gt.mhc2.mean(), pct_mhc2_high_top_quartile=100 * gt.high.mean() if len(gt) else np.nan))
            for s_, gs in g.groupby('section'):
                ROWS.append(dict(analysis='2_per_apc_descriptive', grouping=gname, apc_types='all', CN=c, CN_label=FIGCN[c], section=s_, n_apcs=len(gs),
                                 mean_mhc2_score=gs.mhc2.mean(), pct_mhc2_high_top_quartile=100 * gs.high.mean()))
    # (C) binary: APCs within 50 um of >= 1 CN10 Treg vs APCs near Tregs of other CNs only
    dB = A[A.near_treg].copy(); dB['near_cn10'] = (dB.n_cn10_tregs_50um > 0).astype(float)
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
        print('near-CN10 binary', tset, round(b, 3), (round(b - 1.96 * se, 3), round(b + 1.96 * se, 3)), f'p={p:.2g}', f'{int((sec > 0).sum())}/{len(sec)}', flush=True)
    # ================= 3-5. Treg-level ligand-supply model
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
    def all_cn(col, d_all, tag, extra):
        rs = []
        for c in range(1, 11):
            r = treg_fit(col, d_all, c); r.update(CN=c, **extra); rs.append(r)
        o = pd.DataFrame(rs); okp = o.p.notna() if 'p' in o else pd.Series(False, index=o.index)
        o['q'] = np.nan; o.loc[okp, 'q'] = multipletests(o.loc[okp, 'p'], method='fdr_bh')[1]
        return o
    TH = []
    for thr in THRS:
        col = f'L_thr{thr}'
        W[col] = np.where(W.nAPC >= thr, W.L_apc.values, np.nan)
        if SAME74 and thr in (1, 3):   # use the module columns of the per-Treg table (recomputation from the rounded per-ligand columns agrees to < 1e-3)
            src = f'{M}|L' if thr == 1 else f'{M}|L min3'; assert (W[col].isna() == W[src].isna()).all() and np.nanmax(np.abs(W[col] - W[src])) < 2e-3
            W[col] = W[src]
        o = all_cn(col, W, 'thr', dict(analysis='3_threshold', min_apcs_50um=thr)); TH.append(o)
        for meth in ['bfgs', 'powell']:   # optimiser sensitivity of the CN10 REML fit (lbfgs = main fit)
            d_ = W[['section', 'CN', 'lib_z', col]].dropna().copy(); d_['y'] = d_.groupby('section')[col].transform(zs)
            secs_ = o.loc[o.CN == 10, 'secs'].iloc[0] if isinstance(o.loc[o.CN == 10, 'secs'].iloc[0], list) else []
            if len(secs_) >= 3:
                dd_ = d_[d_.section.isin(secs_)].assign(in_cn=lambda x: (x.CN == 10).astype(float))
                f_ = smf.mixedlm('y ~ in_cn + lib_z', dd_, groups=dd_.section).fit(reml=True, method=meth)
                o.loc[o.CN == 10, f'estimate_{meth}'] = f_.params['in_cn']; o.loc[o.CN == 10, f'p_{meth}'] = f_.pvalues['in_cn']
        r10 = o[o.CN == 10].iloc[0]
        print(f'thr {thr}: CN10 est {r10.estimate:.3f} [{r10.get("ci_low", np.nan):.3f}, {r10.get("ci_high", np.nan):.3f}] q={r10.q:.2g} n_in={r10.n_tregs_in_CN} sections={r10.n_sections} agree={r10.get("sections_agreeing")}', flush=True)
        print(o[['CN', 'n_sections', 'n_tregs_in_CN', 'estimate', 'q']].round(3).to_string(index=False), flush=True)
    TH = pd.concat(TH, ignore_index=True)
    # 4. leave-one-section-out (threshold 1)
    LO = []
    for s_ in sorted(W.section.unique()):
        r = treg_fit('L_thr1', W, 10, drop=s_); r.update(analysis='4_leave_one_section_out', CN=10, dropped_section=s_, min_apcs_50um=1); LO.append(r)
    LO = pd.DataFrame(LO); print(LO[['dropped_section', 'n_sections', 'n_tregs_in_CN', 'estimate', 'ci_low', 'ci_high', 'p']].round(4).to_string(index=False), flush=True)
    # 5. abundance-weighted (sum) supply, all Tregs (0 when no APC within 50 um)
    SU = all_cn('S_sum', W, 'sum', dict(analysis='5_abundance_weighted_sum', min_apcs_50um=0))
    print(SU[['CN', 'n_sections', 'n_tregs_in_CN', 'estimate', 'ci_low', 'ci_high', 'q', 'sections_agreeing']].round(3).to_string(index=False), flush=True)
    # ================= CSV
    for df in (TH, LO, SU):
        df['CN_label'] = df.CN.map(lambda c: FIGCN.get(c, f'CN{c} (omitted in Fig 4D)'))
        df['section'] = 'all (mixed model)'; df['row_type'] = 'mixed_model'
        df['section_deltas_z'] = df.get('section_deltas', pd.Series(index=df.index, dtype=object)).map(lambda v: ';'.join(f'{x:.3g}' for x in v) if isinstance(v, list) else '')
        df['sections_used'] = df.get('secs', pd.Series(index=df.index, dtype=object)).map(lambda v: ';'.join(v) if isinstance(v, list) else '')
        for k in ('section_deltas', 'secs'):
            if k in df: del df[k]
    OUT = pd.concat([pd.DataFrame(ROWS), TH, LO, SU], ignore_index=True)
    OUT['metric_note'] = OUT.analysis.map({
        '1_apc_abundance': 'APCs (DC, B cell, CD9+ M2 Mac, Mo-Mac, SELENOP+ Mac; groups DC / B cell / Macrophage) within 50 um of each strict Treg, by the Treg CN',
        '2_per_apc_mhc2': 'per-APC MHC-II score (mean within-section z of 9 HLA ligands), z within section over APCs; MixedLM y ~ CN + lib_z (+ APC type) + (1|section)',
        '2_per_apc_descriptive': 'per-APC MHC-II score (mean within-section z of 9 HLA ligands, all cells); MHC-II-high = top quartile within section over APCs',
        '3_threshold': 'Fig 4D model on Antigen Presentation ligand supply L, Tregs with >= min_apcs_50um APCs; BH over CNs per threshold',
        '4_leave_one_section_out': 'Fig 4D CN10 model (threshold 1) without the dropped section',
        '5_abundance_weighted_sum': 'SUM of APC MHC-II scores within 50 um (0 if none), all Tregs; Fig 4D model; BH over CNs'})
    lead = ['analysis', 'row_type', 'grouping', 'apc_types', 'min_apcs_50um', 'dropped_section', 'CN', 'CN_label', 'section', 'level', 'contrast',
            'estimate', 'ci_low', 'ci_high', 'p', 'q', 'sections_agreeing', 'sections_CN10_higher', 'delta']
    OUT = OUT[[c for c in lead if c in OUT] + [c for c in OUT if c not in lead]].apply(sig3)
    OUT.to_csv(f'{STEM}.csv', index=False); print('wrote', f'{STEM}.csv', OUT.shape, flush=True)
    # ================= figure (PNG only)
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    LC_, DC_ = '#9ecae1', '#08519c'; HI, HIL = '#d94801', '#fdae6b'
    fig = plt.figure(figsize=(11, 10.6))   # panels a over b, shared x (threshold, leave-one-out and total-supply statistics in the CSV)
    gsp = fig.add_gridspec(2, 1, height_ratios=[0.85, 1], hspace=0.13, top=0.89)
    axa = fig.add_subplot(gsp[0]); axb = fig.add_subplot(gsp[1], sharex=axa)
    rng = np.random.default_rng(1); X = {c: i for i, c in enumerate(CNS)}
    # (a) APCs within 50 um of each Treg, by the Treg's CN (box = all Tregs, whiskers 5-95th pct; dots = section medians, sections with >= 10 Tregs)
    YCAP = 10 * np.ceil(max(np.percentile(W.loc[W.CN == c, 'nAPC'], 95) for c in CNS) / 10)
    for c in CNS:
        v = W.loc[W.CN == c, 'nAPC'].values; col = HI if c == 10 else DC_
        axa.boxplot([v], positions=[X[c]], widths=0.55, whis=(5, 95), showfliers=False, patch_artist=True, medianprops=dict(color=col, lw=2.6),
                    boxprops=dict(facecolor='white', edgecolor=col, lw=1.8), whiskerprops=dict(color=col, lw=1.6), capprops=dict(color=col, lw=1.6))
        sm = W[W.CN == c].groupby('section').nAPC.agg(['median', 'size']); sm = sm[sm['size'] >= MIN_TREG_CN]['median'].values
        axa.scatter(X[c] + rng.uniform(-0.17, 0.17, len(sm)), sm, s=45, color=HIL if c == 10 else LC_, edgecolor='white', linewidth=0.5, zorder=3)
        axa.text(X[c], YCAP * 1.03, f'{100 * (v >= 3).mean():.0f}%', ha='center', va='bottom', fontsize=14, color=col, fontweight='bold' if c == 10 else 'normal')
    assert all(np.percentile(W.loc[W.CN == c, 'nAPC'], 95) <= YCAP for c in CNS)
    axa.text(-0.5, YCAP * 1.13, 'Tregs with ≥3 APCs within 50 µm:', ha='left', va='bottom', fontsize=14, color='#444444')
    axa.set_ylim(-0.02 * YCAP, YCAP * 1.22)
    axa.set_ylabel('APCs within 50 µm\nof Each Treg', fontsize=16)
    axa.set_title('a   CN10 Tregs Have Few APCs Nearby', loc='left', fontsize=17, fontweight='bold', pad=6)
    plt.setp(axa.get_xticklabels(), visible=False); axa.tick_params(axis='x', length=0)
    # (b) per-APC MHC-II score by the APC's own CN, type- and library-size-adjusted, vs CN10 APCs
    o = APC2[('own_CN', 'all APC types (type-adjusted)')]; mmo = o[o.row_type == 'mixed_model'].set_index('CN'); so = o[o.row_type == 'section']
    axb.axhline(0, color=HI, lw=1.6, ls='--', zorder=0)
    for c in OTH:
        v = so[so.CN == c].delta.values
        axb.scatter(X[c] + rng.uniform(-0.15, 0.15, len(v)), v, s=45, color=LC_, edgecolor='white', linewidth=0.5, zorder=2)
        e, lo, hi, q = mmo.loc[c, ['estimate', 'ci_low', 'ci_high', 'q']]
        axb.errorbar(X[c], e, yerr=[[e - lo], [hi - e]], fmt='o', ms=10, color=DC_, mfc=DC_ if q < 0.05 else 'white', mew=2.2, elinewidth=2.6, capsize=0, zorder=3)
    vals = np.concatenate([so.delta.values, mmo.ci_low.values, mmo.ci_high.values]); span = vals.max() - vals.min()
    for c in OTH:
        q = mmo.loc[c, 'q']; top = max(so[so.CN == c].delta.max() if (so.CN == c).any() else -np.inf, mmo.loc[c, 'ci_high'])
        axb.text(X[c], top + 0.025 * span, stars(q), ha='center', va='bottom', fontsize=16 if q < 0.05 else 14, color=DC_, fontweight='bold' if q < 0.05 else 'normal')
    axb.text(X[10], 0, 'Ref.', ha='center', va='center', fontsize=15, color=HI, fontweight='bold', bbox=dict(facecolor='white', edgecolor='none', pad=1))
    axb.set_ylim(vals.min() - 0.05 * span, vals.max() + 0.16 * span)
    axb.set_ylabel('Per-APC MHC-II Score\nvs CN10 APCs (z)', fontsize=16)
    axb.set_xticks(range(len(CNS))); axb.set_xticklabels([FIGCN[c] for c in CNS], rotation=35, ha='right', rotation_mode='anchor', fontsize=14)
    axb.set_xlabel('Neighbourhood of the Treg (a) or of the APC (b)', fontsize=16)
    axb.set_xlim(-0.55, len(CNS) - 0.45)
    axb.set_title('b   CN10 APCs Are Among the Most MHC-II-High', loc='left', fontsize=17, fontweight='bold', pad=6)
    t10 = TH[TH.CN == 10].set_index('min_apcs_50um'); s10 = SU[SU.CN == 10].iloc[0]   # panel-c statistics: CSV + legend only
    for ax in (axa, axb): ax.spines[['top', 'right']].set_visible(False); ax.tick_params(axis='y', labelsize=14)
    # title: statement matching panels a-b, asserted on the data
    assert np.median(W.loc[W.CN == 10, 'nAPC']) < np.median(W.loc[W.CN.isin(OTH), 'nAPC'])   # 'few APCs'
    assert all(mmo.loc[c, 'estimate'] < 0.1 for c in OTH)   # no CN's APCs exceed CN10's by > 0.1 SD
    TITLE = 'CN10 Has Few APCs, but Those APCs Are MHC-II-High'
    fig.suptitle(TITLE, fontsize=19, fontweight='bold', x=0.02, ha='left', y=0.985)
    h = [Line2D([], [], marker='o', ls='', color='#bdbdbd', ms=9, label='One Section'),
         Line2D([], [], marker='o', ls='-', color='#252525', ms=10, lw=2.4, label='Estimate ± 95% CI, q < 0.05'),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=10, mew=2, lw=2.4, label='q ≥ 0.05')]
    fig.legend(handles=h, loc='upper left', bbox_to_anchor=(0.02, 0.962), ncol=3, frameon=False, fontsize=14, handletextpad=0.3, columnspacing=1.2)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    print('wrote', f'{STEM}.png', flush=True)
    # ================= legend
    own = APC2[('own_CN', 'all APC types (type-adjusted)')]; own = own[own.row_type == 'mixed_model'].set_index('CN')
    fmt = lambda e, lo, hi: f'{e:+.2f} (95% CI {lo:+.2f} to {hi:+.2f})'
    hi_lo = ', '.join(f'CN{c}' for c in OTH if own.loc[c, 'q'] < 0.05 and own.loc[c, 'estimate'] < 0)
    hi_hi = ', '.join(f'CN{c}' for c in OTH if own.loc[c, 'q'] < 0.05 and own.loc[c, 'estimate'] > 0)
    ns_ = ', '.join(f'CN{c}' for c in OTH if not own.loc[c, 'q'] < 0.05)
    ab = {r['CN']: r for r in ROWS if r.get('analysis') == '1_apc_abundance' and r.get('level') == 'CN'}
    leg = ('Antigen-presentation (MHC-II) ligand supply to Tregs in the OGN⁺RSPO3⁺ fibroblast neighbourhood (CN10): abundance of the sender cells '
       'versus their per-cell MHC-II expression. CD resections (CosMx 6K; 11 sections, 4 patients; 6,554 strict Tregs). As in Fig 4D, the antigen-'
       'presentation module comprises the NicheNet v2 pairs grouped by the CellChatDB MHC-II pathway with both genes on the 6K panel (HLA-DMA, HLA-DMB, '
       'HLA-DOA, HLA-DOB, HLA-DPB1, HLA-DQA1, HLA-DQA2, HLA-DQB1 and HLA-DRA, each to CD4). The sender cells (antigen-presenting cells, APCs) are DCs, B cells, '
       'CD9+ M2 macrophages, Mo-Macs and SELENOP+ macrophages (types grouped as DC, B cell and Macrophage); SELENOP+ fibroblasts are not APCs. Ligand expression is the within-section z over all cells of the section. '
       'A cell\'s neighbourhood is its own CN_k10. CN6 and CN9 are omitted because they contain too few Tregs. '
       f'(a) The number of APCs within 50 µm of each Treg, by the Treg\'s neighbourhood. Boxes show the median and IQR over all Tregs, and whiskers the 5th–95th percentiles. '
       'Dots are section medians (sections with at least 10 Tregs in the neighbourhood). Numbers above the boxes give the percentage of Tregs with at least 3 APCs within 50 µm. '
       f'CN10 Tregs: median {ab[10]["median_apcs_50um"]:.0f} APCs (IQR {ab[10]["q1_apcs_50um"]:.0f}–{ab[10]["q3_apcs_50um"]:.0f}); CN4 Tregs: median {ab[4]["median_apcs_50um"]:.0f} '
       f'(IQR {ab[4]["q1_apcs_50um"]:.0f}–{ab[4]["q3_apcs_50um"]:.0f}). '
       f'{CONC["n_unique_apcs"]:,} distinct APCs lie within 50 µm of at least one CN10 Treg ({CONC["n_sections_with_such_apcs"]} sections; median {CONC["median_unique_apcs_per_section"]:.0f} per section), '
       f'and the 10% most-connected of them account for {CONC["pct_links_from_top10pct_apcs"]:.0f}% of the Treg–APC links. '
       '(b) The per-APC MHC-II score, analysed at the level of APCs and not through Tregs. Each APC\'s score is the mean of its within-section z over the 9 HLA ligand genes, '
       f'z-standardised within section over all APCs ({len(A):,} APCs). This is the quantity that Fig 4D averages over the APCs near each Treg. '
       'Each APC is grouped by its own neighbourhood. Linear mixed model: score ~ neighbourhood (reference CN10) + APC type + APC library size '
       '(log total counts, within-section z) + (1 | section). Each point is the difference of a neighbourhood from CN10 APCs, with its 95% CI; '
       'negative values mean that CN10 APCs score higher. Light dots: one section each, the mean type- and library-size-adjusted residual in the neighbourhood minus that in CN10 '
       '(sections with at least 10 APCs in both). Stars give the BH q over the 7 contrasts: * q < 0.05, ** q < 0.01, *** q < 0.001; ns, q ≥ 0.05. '
       f'CN10 APCs score higher than APCs in {hi_lo} and do not differ from APCs in {ns_}' +
           (f', and score slightly lower than APCs in {hi_hi} ({fmt(*own.loc[8, ["estimate", "ci_low", "ci_high"]])}, where B cells predominate). ' if hi_hi == 'CN8' else
            ' (CN8, where B cells predominate: ' + fmt(*own.loc[8, ["estimate", "ci_low", "ci_high"]]) + f', q = {own.loc[8, "q"]:.2g}); no neighbourhood\'s APCs score higher than CN10 APCs. ' if not hi_hi else
            f', and score lower than APCs in {hi_hi}. ') +
           'Section-level agreement is moderate; the per-section counts are in the source data. '
       'The source data also give: the same model for each APC type; the model with APCs grouped by the neighbourhood of their nearest Treg within 50 µm; '
       'APCs within 50 µm of a CN10 Treg compared with APCs near Tregs of other neighbourhoods only; and the number and percentage of MHC-II-high APCs '
       '(top quartile within section) per neighbourhood and per section. '
       'Robustness of the Fig 4D CN10 call (Treg-level ligand-supply model, with the five APC types above; statistics in the source data): the CN10 estimate stays positive and significant when Tregs are '
       'scored only with at least 1, 3, 5 or 10 APCs within 50 µm (' + ', '.join(f'{e:+.2f}' for e in t10.loc[THRS, 'estimate']) + f'; q ≤ {t10.q.max():.1g}) and when any one section is left out '
       f'({LO.estimate.min():+.2f} to {LO.estimate.max():+.2f}), whereas the abundance-weighted total supply (sum over the APCs within 50 µm) is lower in CN10 than in other Tregs '
       f'({fmt(s10.estimate, s10.ci_low, s10.ci_high)}, q = {s10.q:.1g}, {s10.sections_agreeing} sections).')
    for bad in ('activation', 'quiescent', 'serial'): assert bad not in leg.lower()
    open(f'{STEM}_legend.txt', 'w').write(leg + '\n'); print(leg)


# ======================================================================================================================
# Section cn_apc_mhcii_table: APC share of cells and model-adjusted per-APC MHC-II score per cellular neighbourhood
# Panels: Supp 9C inputs (supp_cn10_apc_abundance_mhcii.csv and legend)
# ======================================================================================================================
if _run('cn_apc_mhcii_table'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/81b_supp_cn10_apc_two_panel.py'   # working folder of this step (relative paths resolve from here)
    """APC abundance and per-APC MHC-II score per cellular neighbourhood (CD resections, CosMx 6K; 11 sections, 4 patients); the CSV and legend are
read by S9C_cn_apc_mhcii. The log2 enrichment of the CN composition heatmap (fig4_cn_composition_heatmap.csv, 12_neighbourhood_treg_score.py) is
relative abundance; this step gives absolute levels for every CN including CN10.
  Left   'APCs as % of Cells'    <- supp_apc_abundance_by_neighbourhood.csv (16_cxcl12_regions.py, section S9E_cxcl12_by_celltype), values asserted:
                                    stacked DC / B Cell / Macrophage = mean % of all cells; section dots, mean +- 95% CI and mixed-model stars for the total;
                                    'All Cells' line = APCs / all cells per section, averaged over sections, from the same section rows.
  Right  'Per-APC MHC-II Score'  <- per-APC table rebuilt with the definitions of cn10_apc_abundance (APC_CACHE=<pkl> re-uses it); the model
                                    y ~ C(CN, ref CN10) + lib_z + C(APC group) + (1|section) is refitted and asserted to reproduce the CSV of
                                    cn10_apc_abundance (estimates, n, BH q / stars). Bars = model-adjusted means per CN (CN10 intercept + CN effect, at
                                    lib_z = 0 and the overall APC-group mix), 95% CI from the fixed-effect covariance; bar differences = the CN10
                                    contrasts, so bars and stars agree. Dots = per-section mean of the group- and library-adjusted score (sections with
                                    >= 10 APCs in the CN).
The per-Treg APC counts within 50 um (cn10_apc_per_treg_78.csv.gz) are not drawn; they are in the CSV and one legend sentence.
APCs = DC, B cell, Macrophage (= SELENOP+ Mac + CD9+ M2 Mac + Mo-Mac); SELENOP+ Fib is not an APC. CN labels: CN_label column of the APC table.
CN6 / CN9 are not in the per-APC MHC-II model and are left out of the right panel.
Outputs: panels/fig4_and_related/A_fig4_panels/supp_cn10_apc_abundance_mhcii.{png,csv} + _legend.txt"""
    import os, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.legend_handler import HandlerTuple
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUTA = f'{P4}/A_fig4_panels'
    SRC78 = f'{P4}/archive/ext_cn10_antigen_presentation_abundance.csv'; SRC80 = f'{P4}/archive/supp_apc_abundance_by_neighbourhood.csv'
    CACHE78 = f'{P4}/B_resection_analyses/lr_modules_cellchat/cn10_apc_per_treg_78.csv.gz'
    HEAT = f'{OUTA}/fig4_cn_composition_heatmap.csv'
    STEM = f'{OUTA}/supp_cn10_apc_abundance_mhcii'
    APC = ['DC', 'B cell', 'Macrophage']
    APC_DISP = {'DC': 'DC', 'B cell': 'B Cell', 'Macrophage': 'Macrophage'}
    APC_COL = {'DC': '#08519c', 'B cell': '#6baed6', 'Macrophage': '#31a354'}     # APC group colours
    HI, HIL, BAND = '#b2182b', '#e7a3a3', '#fbeaea'      # CN10 highlight
    DC_, LC_ = '#08519c', '#9ecae1'                      # other CNs
    DOT, DOTE = '#525252', '#252525'                     # section dots / total mean ± CI on the stacked bars
    CNO = list(range(1, 11)); CN78 = [1, 2, 3, 4, 5, 7, 8, 10]; OTH78 = [c for c in CN78 if c != 10]
    X1 = {c: i for i, c in enumerate(CNO)}; X2 = {c: i for i, c in enumerate(CN78)}
    S_DOT = 34; MS_EST = 9                                # uniform dot sizes (section dots / estimate markers)
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    sig3 = lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    # ======================================================= inputs
    A80 = pd.read_csv(SRC80)
    T80 = A80[A80.row_type == 'summary'].set_index('CN').loc[CNO]
    SALL = A80[A80.row_type == 'section']                                  # all section x CN rows (for the all-cell reference)
    S80 = SALL[SALL.used_in_mean == True]                                  # rows that enter the CN means (>= 100 cells)
    LAB = T80.CN_label.to_dict()
    assert LAB[10] == 'CN10: OGN⁺RSPO3⁺/PLA2G2A⁺ Fib' and LAB[4] == 'CN4: DC / Mac / CD9 Mac'
    assert len(SALL) == 110 and SALL.section.nunique() == 11
    for c in CNO:
        assert abs(T80.loc[c, [f'mean_pct_{a}' for a in APC]].sum() - T80.loc[c, 'mean_apc_pct']) < 0.02 * max(1, T80.loc[c, 'mean_apc_pct'])
        assert abs(S80[S80.CN == c].mean_apc_pct.mean() - T80.loc[c, 'mean_apc_pct']) < 0.01 * max(1, T80.loc[c, 'mean_apc_pct']) + 0.01
        assert SALL[SALL.CN == c].n_APC.sum() == T80.loc[c, 'n_APC'] and SALL[SALL.CN == c].n_cells.sum() == T80.loc[c, 'n_cells'], c
    # all-cell APC %: per section APCs / all cells (all neighbourhoods), mean (95% CI) over the 11 sections; pooled for reference
    from scipy import stats
    sec_all = SALL.groupby('section')[['n_APC', 'n_cells']].sum(); sec_all['apc_pct'] = 100 * sec_all.n_APC / sec_all.n_cells
    ALL_M = sec_all.apc_pct.mean(); ALL_SE = sec_all.apc_pct.std(ddof=1) / np.sqrt(len(sec_all))
    ALL_LO, ALL_HI = ALL_M + np.array([-1, 1]) * stats.t.ppf(0.975, len(sec_all) - 1) * ALL_SE
    ALL_POOL = 100 * sec_all.n_APC.sum() / sec_all.n_cells.sum()
    A78 = pd.read_csv(SRC78, low_memory=False)
    AB = A78[(A78.analysis == '1_apc_abundance') & (A78.level == 'CN')].set_index('CN')
    MM = A78[(A78.analysis == '2_per_apc_mhc2') & (A78.grouping == 'own_CN') & (A78.apc_types == 'all APC types (type-adjusted)')]
    MMo = MM[MM.row_type == 'mixed_model'].set_index('CN').loc[OTH78]; MMs = MM[MM.row_type == 'section']
    assert set(MMs.CN) == set(OTH78)
    # ---- per-APC MHC-II table with the definitions of cn10_apc_abundance (its compute stage, APC part only)
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
    from statsmodels.tools.sm_exceptions import ConvergenceWarning; warnings.simplefilter('ignore', ConvergenceWarning)
    RAW = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'; MIN_APC_SEC = 10
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    APC5 = ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac']
    APC_GRP = {'DC': 'DC', 'B cell': 'B cell', 'CD9+ M2 Mac': 'Macrophage', 'Mo-Mac': 'Macrophage', 'SELENOP+ Mac': 'Macrophage'}
    HLA = ['HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DRA']
    zs = lambda v: (v - v.mean()) / v.std()
    APCC = os.environ.get('APC_CACHE')
    if APCC and os.path.exists(APCC):
        A = pd.read_pickle(APCC); A = A[0] if isinstance(A, tuple) else A          # (A, TT) cache of cn10_apc_abundance or the A cache of this section
        A = A[['section', 'apc_type', 'log_total_counts', 'mhc2', 'CN']].copy()
    else:
        import gc, h5py, pyarrow.parquet as pq
        from anndata.io import read_elem
        from scipy import sparse
        sa = pq.read_table(f'{RAW}/combined/all_samples.parquet', columns=['cell_type_general', 'patient_sec', 'total_counts', '__index_level_0__']).to_pandas()
        sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]; sa['ct'] = sa.cell_type_general.astype(str)
        assert set(APC5) <= set(sa.ct.unique())
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
        ab = []
        for s_, f in FOLDER.items():
            sub = sa[sa.section == s_]
            with h5py.File(f'{RAW}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
                ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
            nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
            keep = np.where(nm.isin(sub.index))[0]
            E = X[keep][:, [v.get_loc(g) for g in HLA]].toarray().astype(np.float32); del X; gc.collect()
            idx = nm[keep]; cs = sub.loc[idx]
            sd = E.std(0); Z = np.where(sd > 0, (E - E.mean(0)) / np.where(sd > 0, sd, 1), np.nan)   # z over all cells of the section
            am = np.isin(cs.ct.values, APC5)
            a = pd.DataFrame({'section': s_, 'apc_type': cs.ct.values[am], 'log_total_counts': np.log1p(cs.total_counts.values[am].astype(float)),
                              'mhc2': np.nanmean(Z[am], 1)}, index=idx[am])
            a['CN'] = cn.loc[a.index, 'CN_k10'].astype(int).values; assert (cn.loc[a.index, 'section'].astype(str).values == s_).all()
            ab.append(a); print(s_, 'APCs', len(a), flush=True); del E, Z; gc.collect()
        A = pd.concat(ab); del sa, cn; gc.collect()
        if APCC: pd.to_pickle(A, APCC)
    A['cell_type'] = A.apc_type; A['apc_type'] = A.cell_type.map(APC_GRP); assert A.apc_type.notna().all()
    A['lib_z'] = A.groupby('section').log_total_counts.transform(zs)
    A['y'] = A.groupby('section').mhc2.transform(zs)                               # z within section over all APCs
    dA = A[A.CN.isin(CN78)].copy(); dA['grp'] = dA.CN.astype(int)
    FRM = 'y ~ C(grp, Treatment(10)) + lib_z + C(apc_type)'                        # own-CN model with all APC groups, as in cn10_apc_abundance
    # optimisers in the order of cn10_apc_abundance; a fit is accepted only with a finite likelihood and finite intercept / fixed-effect SEs, because the adjusted
    # means use the intercept. (lbfgs stops at the boundary for this model: section variance 0, infinite log-likelihood, intercept stuck at 0. Its CN contrasts
    # are unaffected - they equal the section-fixed-effects OLS contrasts - but its intercept is unusable; bfgs / powell converge to the same contrasts.)
    FIT, METH = None, None
    for meth in ['lbfgs', 'bfgs', 'powell', 'nm']:
        try:
            f_ = smf.mixedlm(FRM, dA, groups=dA.section).fit(reml=True, method=meth)
            if np.isfinite(f_.llf) and all(np.isfinite(f_.bse[k]) and f_.bse[k] > 0 for k in f_.fe_params.index): FIT, METH = f_, meth; break
            print('optimiser', meth, 'rejected (llf', f_.llf, ')', flush=True)
        except Exception as e: print('optimiser', meth, 'failed', e, flush=True)
    assert FIT is not None; print('MHC-II model optimiser:', METH, '| section variance', float(FIT.cov_re.iloc[0, 0]), flush=True)
    KC = {c: f'C(grp, Treatment(10))[T.{c}]' for c in OTH78}
    qv = multipletests([FIT.pvalues[KC[c]] for c in OTH78], method='fdr_bh')[1]
    for i, c in enumerate(OTH78):   # the refit reproduces the CSV of cn10_apc_abundance (3 significant figures) and its stars
        assert np.isclose(FIT.params[KC[c]], MMo.loc[c, 'estimate'], rtol=5e-3, atol=5e-4), (c, FIT.params[KC[c]], MMo.loc[c, 'estimate'])
        assert np.isclose(FIT.params[KC[c]] - 1.96 * FIT.bse[KC[c]], MMo.loc[c, 'ci_low'], rtol=5e-3, atol=5e-4), c
        assert float(f"{(dA.grp == c).sum():.3g}") == MMo.loc[c, 'n_apcs_in_CN'] and stars(qv[i]) == stars(MMo.loc[c, 'q']), c
    assert float(f"{(dA.grp == 10).sum():.3g}") == MMo.n_apcs_CN10.iloc[0]
    # model-adjusted mean per CN: intercept (CN10) + CN effect + lib_z = 0 + APC-group effects weighted by the overall APC-group mix of the modelled APCs
    FE = FIT.fe_params; COV = FIT.cov_params().loc[FE.index, FE.index]
    WMIX = dA.apc_type.value_counts(normalize=True)
    KT = {g: f'C(apc_type)[T.{g}]' for g in WMIX.index if f'C(apc_type)[T.{g}]' in FE.index}
    assert len(KT) == 2                                                            # 3 groups, one reference
    ADJ = {}
    for c in CN78:
        L = pd.Series(0.0, index=FE.index); L['Intercept'] = 1.0
        if c != 10: L[KC[c]] = 1.0
        for g, k in KT.items(): L[k] = WMIX[g]
        m = float(L @ FE); se = float(np.sqrt(L @ COV @ L)); assert np.isfinite(m) and np.isfinite(se); ADJ[c] = (m, m - 1.96 * se, m + 1.96 * se)
    for c in OTH78: assert np.isclose(ADJ[c][0] - ADJ[10][0], FIT.params[KC[c]])  # bar differences = the CN10 contrasts
    # section dots: score adjusted to the same reference (lib_z = 0, overall APC-group mix) with the model's fixed effects; mean per section x CN
    teff = dA.apc_type.map(lambda g: FE[KT[g]] if g in KT else 0.0)
    dA['y_adj'] = dA.y - FE['lib_z'] * dA.lib_z - teff + sum(WMIX[g] * FE[k] for g, k in KT.items())
    SEC2 = dA.groupby(['grp', 'section']).agg(y_adj=('y_adj', 'mean'), n=('y_adj', 'size')).reset_index(); SEC2 = SEC2[SEC2.n >= MIN_APC_SEC]
    print('adjusted means', {c: round(v[0], 3) for c, v in ADJ.items()}, '| mean of section dots', SEC2.groupby('grp').y_adj.mean().round(3).to_dict(), flush=True)
    W = pd.read_csv(CACHE78, index_col=0, usecols=['Unnamed: 0', 'section', 'CN', 'nAPC', 'nb_n'])
    assert len(W) == 6554; W['nAPC'] = W['nAPC'].astype(int)
    for c in CN78:   # the per-Treg numbers (legend sentence + CSV) agree with the CSV of cn10_apc_abundance
        v = W.loc[W.CN == c, 'nAPC'].values
        assert np.median(v) == AB.loc[c, 'median_apcs_50um'] and np.percentile(v, 25) == AB.loc[c, 'q1_apcs_50um'], c
        assert np.isclose(np.percentile(v, 75), AB.loc[c, 'q3_apcs_50um'], rtol=5e-3) and np.isclose(100 * (v >= 3).mean(), AB.loc[c, 'pct_ge3'], rtol=5e-3), c
    HM = pd.read_csv(HEAT, index_col=0); assert HM.index[9] == LAB[10]
    HM_COLS = ['DC', 'B', 'CD9 Mac', 'Mac S', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac']
    H10 = HM.loc[LAB[10], HM_COLS]
    # ======================================================= figure
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw=dict(wspace=0.34, width_ratios=[len(CNO), len(CN78) + 0.4]))
    rng = np.random.default_rng(2)
    for ax, X in ((a1, X1), (a2, X2)):
        ax.axvspan(X[10] - 0.48, X[10] + 0.48, color=BAND, zorder=0, lw=0)
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(axis='y', labelsize=14)
    TFS = 16
    # ---- left: stacked APC composition (% of all cells) + total section dots / mean ± CI / stars
    bot = np.zeros(len(CNO))
    for a in APC:
        w = T80.loc[CNO, f'mean_pct_{a}'].values
        a1.bar(range(len(CNO)), w, bottom=bot, width=0.74, color=APC_COL[a], edgecolor='white', linewidth=0.6, zorder=2)
        bot += w
    a1.bar([X1[10]], [bot[X1[10]]], width=0.74, facecolor='none', edgecolor=HI, linewidth=2.2, zorder=3)
    a1.axhline(ALL_M, color='#737373', lw=1.2, ls='--', zorder=1.5)
    ymax = max(S80.mean_apc_pct.max(), T80.ci_high.max())
    for c in CNO:
        i = X1[c]; r = T80.loc[c]; v = S80[S80.CN == c].mean_apc_pct.values
        a1.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=S_DOT, color=DOT, edgecolor='white', linewidth=0.5, zorder=4, alpha=0.85)
        a1.plot([i, i], [r.ci_low, r.ci_high], color=DOTE, lw=2.4, zorder=5, solid_capstyle='butt')
        a1.plot([i - 0.17, i + 0.17], [r.mean_apc_pct] * 2, color=DOTE, lw=2.4, zorder=5, solid_capstyle='butt')
        top = max(v.max(), r.ci_high)
        if c != 10 and r.q < 0.05: a1.text(i, top + 0.012 * ymax, stars(r.q), ha='center', va='bottom', fontsize=15, color='#252525')
    # 'All Cells' dashed line is keyed in the figure legend
    a1.set_ylim(0, ymax * 1.10); a1.set_ylabel('APCs (% of All Cells in the CN)', fontsize=15)
    a1.set_title('APCs as % of Cells', loc='left', fontsize=TFS, pad=5)
    a1.set_xlim(-0.6, len(CNO) - 0.4)
    # ---- right: per-APC MHC-II score by the APC's own CN vs CN10
    BARC = '#c6c6c6'                                     # single-quantity bars (neutral; CN10 outlined in red as in the left panel)
    a2.axhline(0, color='#252525', lw=0.9, zorder=1)
    for c in CN78:
        i = X2[c]; m, lo, hi = ADJ[c]; v = SEC2[SEC2.grp == c].y_adj.values
        a2.bar([i], [m], width=0.74, color=BARC, edgecolor=HI if c == 10 else 'none', linewidth=2.2 if c == 10 else 0, zorder=2)
        a2.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=S_DOT, color=DOT, edgecolor='white', linewidth=0.5, zorder=4, alpha=0.85)
        a2.plot([i, i], [lo, hi], color=DOTE, lw=2.4, zorder=5, solid_capstyle='butt')
        a2.plot([i - 0.17, i + 0.17], [m] * 2, color=DOTE, lw=2.4, zorder=5, solid_capstyle='butt')
    vals = np.concatenate([SEC2.y_adj.values, [v_ for c in CN78 for v_ in ADJ[c]], [0]]); span = vals.max() - vals.min()
    for c in CN78:
        top = max(SEC2[SEC2.grp == c].y_adj.max(), ADJ[c][2], ADJ[c][0], 0)
        if c == 10: a2.text(X2[c], top + 0.012 * span, 'Ref.', ha='center', va='bottom', fontsize=13, color=HI)
        elif MMo.loc[c, 'q'] < 0.05: a2.text(X2[c], top + 0.012 * span, stars(MMo.loc[c, 'q']), ha='center', va='bottom', fontsize=15, color='#252525')
    a2.set_ylim(vals.min() - 0.04 * span, vals.max() + 0.11 * span)
    a2.set_ylabel('MHC-II Score per APC (z)', fontsize=15)
    a2.set_title('Per-APC MHC-II Score', loc='left', fontsize=TFS, pad=5)
    a2.set_xlim(-0.6, len(CN78) - 0.4)
    for ax, cns in ((a1, CNO), (a2, CN78)):
        ax.set_xticks(range(len(cns))); ax.set_xticklabels([LAB[c] for c in cns], rotation=40, ha='right', rotation_mode='anchor', fontsize=13)
        for t_, c in zip(ax.get_xticklabels(), cns):
            if c == 10: t_.set_color(HI)
    # figure key, two rows: row 1 = APC groups, All Cells line, CN10 band; row 2 = section dots, mean ± CI, stars (same markers in both panels)
    hs = Line2D([], [], marker='o', ls='', color=DOT, mec='white', ms=7.5)
    he = Line2D([], [], marker='+', ls='', color=DOTE, ms=12, mew=2.4)
    R1 = [(Patch(facecolor=APC_COL[a], edgecolor='none'), APC_DISP[a]) for a in APC] + [
          (Line2D([], [], color='#737373', lw=1.2, ls='--'), 'All Cells'), (Patch(facecolor=BAND, edgecolor=HI, lw=1.5), 'CN10 (Ref.)')]
    R2 = [(hs, 'One Section'), (he, 'Mean ± 95% CI'), (Line2D([], [], ls=''), '*, **, *** q < 0.05, 0.01, 0.001 vs CN10 (mixed model, BH)')]
    for k, R in enumerate((R1, R2)):   # two independent single-row legends (rows not forced into shared columns)
        fig.add_artist(fig.legend([h_ for h_, _ in R], [l_ for _, l_ in R], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.4)}, loc='lower left',
                                  bbox_to_anchor=(0.0, 1.075 - 0.075 * k), ncol=len(R), frameon=False, fontsize=13.5, handletextpad=0.35,
                                  columnspacing=1.3, handlelength=1.6, borderaxespad=0.1))
    fig.subplots_adjust(top=0.92, left=0.08, right=0.995, bottom=0.08)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig); print('wrote', f'{STEM}.png', flush=True)
    # ======================================================= CSV (plotted data + legend numbers, one 'panel' column)
    p1 = A80[A80.row_type.isin(['summary', 'section'])][['row_type', 'CN', 'CN_label', 'section', 'n_cells', 'n_APC', 'mean_apc_pct', 'ci_low', 'ci_high',
          'mean_pct_DC', 'mean_pct_B cell', 'mean_pct_Macrophage', 'share_of_APCs_DC', 'share_of_APCs_B cell', 'share_of_APCs_Macrophage',
          'estimate_vs_CN10', 'estimate_ci_low', 'estimate_ci_high', 'p', 'q', 'used_in_mean']].copy()
    p1 = p1.rename(columns={'mean_apc_pct': 'apc_pct_of_cells'}).assign(panel='APCs as % of Cells',
          metric='APCs (DC, B cell, Macrophage = SELENOP+ Mac + CD9+ M2 Mac + Mo-Mac) as % of all cells in the CN; summary: mean over sections (used_in_mean, >= 100 cells) '
             'with 95% CI (t), stacked mean_pct_<group>, share_of_APCs pooled; mixed model APC % ~ CN (ref CN10) + (1|section), BH over 9')
    p1.loc[p1.row_type == 'summary', 'section'] = 'mean of sections'
    pa = pd.concat([pd.DataFrame(dict(row_type='section', CN='all', CN_label='All Cells', section=s, n_cells=r.n_cells, n_APC=r.n_APC, apc_pct_of_cells=r.apc_pct), index=[0])
                    for s, r in sec_all.iterrows()] + [pd.DataFrame(dict(row_type='summary', CN='all', CN_label='All Cells', section='mean of sections',
                    n_cells=sec_all.n_cells.sum(), n_APC=sec_all.n_APC.sum(), apc_pct_of_cells=ALL_M, ci_low=ALL_LO, ci_high=ALL_HI, pooled_apc_pct=ALL_POOL), index=[0])],
                   ignore_index=True).assign(panel='APCs as % of Cells', metric='All Cells line: APCs / all cells of the section (all CNs), mean over 11 sections; pooled_apc_pct = pooled over sections')
    p2 = pd.DataFrame([dict(row_type='model_adjusted_mean', CN=c, section='all (mixed model)', adjusted_mean=ADJ[c][0], ci_low=ADJ[c][1], ci_high=ADJ[c][2],
                            estimate_vs_CN10=FIT.params[KC[c]] if c != 10 else 0.0, p=FIT.pvalues[KC[c]] if c != 10 else np.nan,
                            q=qv[OTH78.index(c)] if c != 10 else np.nan, q_78=MMo.loc[c, 'q'] if c != 10 else np.nan, n_apcs_in_CN=int((dA.grp == c).sum()),
                            optimiser=METH, weight_DC=WMIX['DC'], weight_B_cell=WMIX['B cell'], weight_Macrophage=WMIX['Macrophage']) for c in CN78] +
                      [dict(row_type='section', CN=r.grp, section=r.section, adjusted_mean=r.y_adj, n_apcs_in_CN=r.n) for r in SEC2.itertuples()])
    p2 = p2.assign(CN_label=p2.CN.map(LAB), panel='Per-APC MHC-II Score',
                   metric='per-APC MHC-II score (mean within-section z of 9 HLA ligands, z within section over APCs; 78 definitions). Mixed model (78) y ~ CN (ref CN10) + APC group '
                      '(DC / B cell / Macrophage) + lib_z + (1|section); adjusted_mean = intercept + CN effect at lib_z = 0 and the overall APC-group mix (weight_*), 95% CI from '
                      'the fixed-effect covariance; estimate_vs_CN10 / q = contrast vs CN10, BH over 7 (q_78 = value in 78 CSV); section = mean adjusted score (>= 10 APCs)')
    p2b = MM[['row_type', 'CN', 'section', 'estimate', 'ci_low', 'ci_high', 'p', 'q', 'delta', 'n_apcs_in_CN', 'n_apcs_CN10', 'sections_CN10_higher']].copy()
    p2b = p2b.assign(CN_label=p2b.CN.map(LAB), panel='Per-APC MHC-II Score: contrasts vs CN10 from 78 (stars; not plotted as points)',
                     metric='78 CSV rows: CN minus CN10; section = adjusted residual mean difference')
    bs = []
    for c in CN78:
        v = W.loc[W.CN == c, 'nAPC'].values
        bs.append(dict(row_type='summary', CN=c, CN_label=LAB[c], section='all Tregs', n_tregs=len(v), median_apcs_50um=np.median(v), q1_apcs_50um=np.percentile(v, 25),
                       q3_apcs_50um=np.percentile(v, 75), p5_apcs_50um=np.percentile(v, 5), p95_apcs_50um=np.percentile(v, 95), pct_tregs_ge3_apcs=100 * (v >= 3).mean(),
                       median_all_cells_50um=W.loc[W.CN == c, 'nb_n'].median()))
    p3 = pd.DataFrame(bs).assign(panel='APCs Within 50 µm of Tregs (not plotted; legend sentence)',
         metric='APCs within 50 um of each strict Treg, by the Treg CN (from 78); median_all_cells_50um = cells of any type within 50 um')
    p4 = pd.DataFrame(dict(row_type='summary', CN=10, CN_label=LAB[10], section='pooled', heatmap_cell_type=HM_COLS, log2_enrichment_vs_all_cells=H10.values)).assign(
         panel='Composition heatmap CN10 APC cells (not plotted; legend reference)', metric='log2 enrichment of the cell type in CN10 vs all cells, read from fig4_cn_composition_heatmap.csv')
    OUT = pd.concat([p1, pa, p2, p2b, p3, p4], ignore_index=True)
    lead = ['panel', 'row_type', 'CN', 'CN_label', 'section']
    OUT = OUT[lead + [c for c in OUT if c not in lead and c != 'metric'] + ['metric']].apply(sig3)
    with open(f'{STEM}.csv', 'w') as fh: OUT.to_csv(fh, index=False)
    print('wrote', f'{STEM}.csv', OUT.shape, flush=True)
    # ======================================================= legend (all numbers from the data)
    r10 = T80.loc[10]; c4, c8 = T80.loc[4], T80.loc[8]
    ratio = r10.mean_apc_pct / T80.loc[[4, 8], 'mean_apc_pct'].mean()
    FRAC = 'about a third of' if 0.25 < ratio < 0.42 else 'about half of' if 0.42 <= ratio < 0.6 else f'{ratio:.2f} times'
    SH = sorted(((r10[f'share_of_APCs_{a}'], a) for a in APC), reverse=True)
    NAME = {'Macrophage': 'macrophages', 'DC': 'DCs', 'B cell': 'B cells'}
    COMP = ', '.join(f'{NAME[a]} {v:.0f}%' for v, a in SH)
    assert SH[0][1] == 'Macrophage' and SH[1][1] == 'DC'                     # 'mostly macrophages and DCs'
    more = [c for c in CNO if c != 10 and T80.loc[c, 'q'] < 0.05 and T80.loc[c, 'estimate_vs_CN10'] < 0]   # CN10 has more APCs than these
    fewer = [c for c in CNO if c != 10 and T80.loc[c, 'q'] < 0.05 and T80.loc[c, 'estimate_vs_CN10'] > 0]
    ns1 = [c for c in CNO if c != 10 and not T80.loc[c, 'q'] < 0.05]
    assert fewer == [4, 8] and r10.mean_apc_pct > ALL_M
    cl = lambda L: ', '.join(f'CN{c}' for c in L[:-1]) + (' and ' if len(L) > 1 else '') + f'CN{L[-1]}'
    lo2 = [c for c in OTH78 if MMo.loc[c, 'q'] < 0.05 and MMo.loc[c, 'estimate'] < 0]
    hi2 = [c for c in OTH78 if MMo.loc[c, 'q'] < 0.05 and MMo.loc[c, 'estimate'] > 0]
    ns2 = [c for c in OTH78 if not MMo.loc[c, 'q'] < 0.05]
    assert 4 in ns2 and 8 in ns2 and not hi2                                 # 'matches CN4 / CN8'
    f2 = lambda v: '0.00' if abs(v) < 0.005 else f'{v:+.2f}'
    fa = lambda c: f'{f2(ADJ[c][0])}, 95% CI {f2(ADJ[c][1])} to {f2(ADJ[c][2])}'
    MIXT = ', '.join(f'{NAME[g]} {100 * WMIX[g]:.0f}%' for g in ['DC', 'B cell', 'Macrophage'])
    lo2_txt = 'differences ' + ' to '.join(f2(v) for v in sorted([MMo.loc[c, 'estimate'] for c in lo2])[::len(lo2) - 1]) if len(lo2) > 1 else ''
    fm = lambda c: f'{f2(MMo.loc[c, "estimate"])}, 95% CI {f2(MMo.loc[c, "ci_low"])} to {f2(MMo.loc[c, "ci_high"])}'
    nb10, nb4 = W.loc[W.CN == 10, 'nb_n'].median(), W.loc[W.CN == 4, 'nb_n'].median(); assert nb10 < 0.5 * nb4
    hm_txt = ', '.join(f'{k} {f2(v)}' for k, v in zip(['DC', 'B cell'], H10[['DC', 'B']].values))
    mac_rng = H10[['CD9 Mac', 'Mac S', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac']]
    leg = ('Antigen-presenting cell (APC) abundance and per-APC MHC-II expression across the cellular neighbourhoods (CN1–CN10, Fig 4) of the CD resections '
       f'(CosMx 6K; 11 sections, 4 patients; {int(T80.n_cells.sum()):,} cells). CN10 is the OGN⁺RSPO3⁺/PLA2G2A⁺ fibroblast neighbourhood (shaded; red). '
       'APCs are DCs, B cells and macrophages (SELENOP⁺ macrophages, CD9⁺ M2 macrophages and Mo-Macs); SELENOP⁺ fibroblasts are not counted. '
       'APCs as % of Cells: stacked bars give the mean percentage of all cells of the neighbourhood made up by each APC group; dots are the total APC percentage per section '
       '(section × neighbourhood combinations with at least 100 cells), black lines the mean and 95% CI (t distribution) across sections, and the dashed line the '
       f'all-cell APC percentage (mean of sections). Stars: linear mixed model APC % ~ neighbourhood (reference CN10) + (1 | section), Benjamini–Hochberg (BH) q over 9 contrasts. '
       'Per-APC MHC-II Score: each APC\'s mean within-section z over the 9 HLA ligand genes of the Fig 4D antigen-presentation module, z-standardised within section '
       'over APCs (0 = average APC of the section) and grouped by the APC\'s own neighbourhood. Bars are model-adjusted means from the linear mixed model '
       'score ~ neighbourhood (reference CN10) + APC group + library size + (1 | section): the CN10 intercept plus each neighbourhood\'s effect, at average library size and '
       f'the overall APC-group mix ({MIXT}); black lines are 95% CIs. Differences between bars are the model contrasts with CN10 (stars; BH q over 7 contrasts). '
       'Dots: per-section means of the group- and library-size-adjusted score (sections with at least 10 APCs in the neighbourhood); the model weights sections by their APC numbers, so a bar need not sit at the centre of its dots. '
       'CN6 and CN9 are not in this model (too few Tregs, as in Fig 4D). * q < 0.05, ** q < 0.01, *** q < 0.001. '
       f'CN10 contains APCs: {r10.mean_apc_pct:.1f}% of its cells (95% CI {r10.ci_low:.1f}–{r10.ci_high:.1f}%; {COMP} of its APCs), {FRAC} the APC-rich '
       f'neighbourhoods CN4 ({c4.mean_apc_pct:.0f}%) and CN8 ({c8.mean_apc_pct:.0f}%), significantly more than {cl(more)}, and not different from {cl(ns1)}. '
       f'CN10 is slightly above the all-cell APC percentage ({ALL_M:.1f}%, 95% CI {ALL_LO:.1f}–{ALL_HI:.1f}%). The log2 enrichment in the Fig 4 composition heatmap (CN10: {hm_txt}, macrophage subsets '
       f'{mac_rng.min():+.2f} to {mac_rng.max():+.2f}) is abundance relative to all cells, not absolute abundance; near-zero or negative values mean average or '
       'lower proportions, not absence of APCs. '
       f'Per APC, the adjusted MHC-II score of CN10 ({fa(10)}) matches CN4 ({fa(4)}) and CN8 ({fa(8)}; differences from CN10 {f2(MMo.loc[4, "estimate"])} and '
       f'{f2(MMo.loc[8, "estimate"])}, q ≥ 0.05) and is higher than in {cl(lo2)} ({lo2_txt}). '
       f'Around CN10 Tregs, APCs within 50 µm are fewer (median {AB.loc[10, "median_apcs_50um"]:.0f}, IQR {AB.loc[10, "q1_apcs_50um"]:.0f}–{AB.loc[10, "q3_apcs_50um"]:.0f}; '
       f'{AB.loc[10, "pct_ge3"]:.0f}% of Tregs with at least 3) than around CN4 Tregs (median {AB.loc[4, "median_apcs_50um"]:.0f}), partly because stromal tissue is less '
       f'cell-dense (median {nb10:.0f} versus {nb4:.0f} cells within 50 µm). Per-section values and model estimates are in the accompanying CSV.')
    for bad in ('activation', 'territory', 'quiescent', 'serial'): assert bad not in leg.lower()
    with open(f'{STEM}_legend.txt', 'w') as fh: fh.write(leg + '\n')
    print(leg)


# ======================================================================================================================
# Section cn10_apc_matched_score: APCs near CN10, CN4 and CN8 Tregs and the Treg Suppression Score at matched APC exposure
# Panels: Supp 9C inputs (APC counts near Tregs quoted in the legend)
# ======================================================================================================================
if _run('cn10_apc_matched_score'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/81d_cn10_apc_matched_treg.py'   # working folder of this step (relative paths resolve from here)
    """Treg-centred, APC-matched comparison of CN10 Tregs with CN4 / CN8 Tregs (CD resections, CosMx 6K; 11 sections, 4 patients; 6,554 strict Tregs).
The CSV rows '1 APC exposure per niche: mean of section means' give the APC counts near Tregs quoted in the Supp 9C legend.
  Per-Treg table      = ligand_receptor_state/treg_local_ligand_receptor.csv.gz (11_neighbourhoods_local_tables.py): apc_n = APCs within 50 um,
                        nb_n = all cells within 50 um, apc_frac = apc_n / nb_n, act_z; APC = DC, B cell, SELENOP+ Mac, CD9+ M2 Mac, Mo-Mac; Tregs in
                        the order of combined/raw_counts_strict_tregs.csv.
  APC quintile        = within-section quintile of apc_frac (the APC FRACTION within 50 um, not the count), as in the CXCL12 analysis; the per-quintile
                        CXCL12-High vs -Low table is recomputed and asserted equal to c_apc_matched_cxcl12_high_vs_low_exp.csv (11_neighbourhoods_local_tables.py).
  Treg Suppression Score = -signature_score_D6 z within section (D_shared_inputs/treg_definition_comparison_per_cell.csv); asserted = act_z.
  CN                  = cellular_neighborhoods_v2/cn_per_cell.parquet CN_k10. Tested CNs = the rule of the CN forest plot (>= 100 Tregs and >= 10 Tregs
                        in >= 5 sections), asserted against fig4_cn_treg_state_forest.csv (12_neighbourhood_treg_score.py).
  Library size (cz)   = within-section z of log1p(total raw counts) of the Treg (raw_counts_strict_tregs.csv).
  DC count            = treg_local_composition.csv.gz 'DC' fraction x nb_n (APC columns asserted to sum to apc_n).
  MHC-II-weighted APC exposure = S_sum of cn10_apc_abundance: sum over the APCs within 50 um of each APC's MHC-II score (mean within-section z, over
                        all cells, of the 9 HLA ligand genes); APC_CACHE=<the (A, TT) pickle of cn10_apc_abundance> re-uses it; APC counts asserted
                        against cn10_apc_per_treg_78.csv.gz.
Analyses
  1 APC exposure per niche: per-section mean APC count within 50 um of the Tregs of each tested CN (sections with >= 10 Tregs in the CN),
    mean +- t 95% CI; mixed model log1p(apc_n) ~ C(CN, ref CN10) + (1|section) on Tregs of the tested CNs, BH over the 7 contrasts.
  2 APC-matched score: per APC quintile, per-section mean score of CN4 / CN8 / CN10 Tregs (>= 3 Tregs per section x CN x quintile), mean +- t CI;
    model on CN4 / CN8 / CN10 Tregs: score ~ C(CN, ref CN10) + log1p(apc_n) + cz + (1|section), BH over the 2 contrasts.
    Sensitivities (CSV): no exposure term; log1p(DC count); apc_frac; APC quintile as factor; log1p(nb_n) added; MHC-II-weighted exposure (S_sum).
CN labels: CN_label column of supp_apc_abundance_by_neighbourhood.csv (16_cxcl12_regions.py, section S9E_cxcl12_by_celltype).
Outputs: panels/fig4_and_related/A_fig4_panels/fig4_cn10_apc_matched_treg_score.{png,csv} + _legend.txt"""
    import os, sys, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    from statsmodels.tools.sm_exceptions import ConvergenceWarning; warnings.simplefilter('ignore', ConvergenceWarning)
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUTA = f'{P4}/A_fig4_panels'
    LR = f'{P4}/B_resection_analyses/ligand_receptor_state'; CNV = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; RAW = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
    STEM = f'{OUTA}/fig4_cn10_apc_matched_treg_score'
    UM, RAD = 0.12028, 50
    APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']
    HLA = ['HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DRA']
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    MIN_TREG_TOTAL, MIN_TREG_SEC, MIN_SEC = 100, 10, 5        # testing rule of the CN forest plot (12_neighbourhood_treg_score.py)
    MIN_Q = 3                                                  # Tregs per section x CN x quintile for a section mean (right panel)
    MIN_QSEC = 3                                               # sections needed to draw a quintile point (mean +- CI)
    HI, HIE, DK, LT = '#b2182b', '#4a0a14', '#3d3d3d', '#c4c4c4'
    zs = lambda v: (v - v.mean()) / v.std()
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    sig3 = lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    def tci(v):
        v = np.asarray(v, float); n = len(v); m = v.mean()
        if n < 2: return m, np.nan, np.nan, n
        h = stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n); return m, m - h, m + h, n
    def lmm(formula, d, terms):
        """MixedLM (REML), optimisers lbfgs, bfgs, powell, nm; a fit is accepted only with a finite likelihood and finite SEs (lbfgs can stop at the boundary)."""
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm']:
            try:
                f = smf.mixedlm(formula, d, groups=d.section).fit(reml=True, method=meth)
                if np.isfinite(f.llf) and all(np.isfinite(f.bse[k]) and f.bse[k] > 0 for k in f.fe_params.index): return f, meth
            except Exception: pass
        raise RuntimeError(formula)
    # ======================================================= per-Treg table
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)[['section', 'x', 'y', 'act_z', 'apc_frac', 'nb_n', 'apc_n']]
    T = T.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left')
    cnt = pd.read_csv(f'{RAW}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    assert len(cnt) == len(T) == 6554
    T = T.loc[cnt.index]                                       # Treg order of the raw-count table (rank ties)
    T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    T['cz'] = T.groupby('section').cz.transform(zs)            # library size
    # APC quintile of the CXCL12 analysis, recomputed
    T['apc_q'] = T.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    T['t_exp'] = T.groupby('section').lig_CXCL12_exp.transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))
    assert (T.groupby('section').apc_q.max() == 4).all() and T.apc_q.notna().all()
    CQ = pd.read_csv(f'{LR}/c_apc_matched_cxcl12_high_vs_low_exp.csv').set_index('apc_quintile')
    rows = []
    for q, g in T[T.t_exp.isin(['low', 'high']) & T.apc_q.notna()].groupby('apc_q'):
        hi_, lo_ = g[g.t_exp == 'high'].act_z.dropna(), g[g.t_exp == 'low'].act_z.dropna()
        rows.append(dict(apc_quintile=int(q) + 1, n_high=len(hi_), n_low=len(lo_), hi=hi_.mean(), lo=lo_.mean(), mf=g.apc_frac.mean(), p=stats.mannwhitneyu(hi_, lo_).pvalue))
    RQ = pd.DataFrame(rows).set_index('apc_quintile'); RQ['q'] = multipletests(RQ.p, method='fdr_bh')[1]
    for qq in RQ.index:
        a_, b_ = RQ.loc[qq], CQ.loc[qq]
        assert a_.n_high == b_.n_high and a_.n_low == b_.n_low, qq
        assert np.isclose(a_.hi, b_.suppression_cxcl12_high, rtol=2e-3, atol=1e-4) and np.isclose(a_.lo, b_.suppression_cxcl12_low, rtol=2e-3, atol=1e-4), qq
        assert np.isclose(a_.mf, b_.mean_apc_frac, rtol=2e-3) and np.isclose(a_.q, b_.q, rtol=5e-3), qq
    print('APC quintiles agree with c_apc_matched_cxcl12_high_vs_low_exp.csv (5 quintiles, n and means asserted)', flush=True)
    T['apc_q'] = T.apc_q.astype(int) + 1                       # 1 = lowest APC fraction ... 5 = highest
    # score and CN
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); PC = PC[PC.cohort == 'resection']
    assert set(PC.index) == set(T.index)
    sc = (-PC.signature_score_D6).groupby(PC.section).transform(zs)
    assert (PC.section.loc[T.index].values == T.section.values).all() and np.allclose(sc.loc[T.index].values, T.act_z.values, atol=2e-3)
    T['score'] = sc.loc[T.index].values                        # Treg Suppression Score (full precision)
    cn = pq.read_table(f'{CNV}/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas()
    T['CN'] = cn.loc[T.index, 'CN_k10'].astype(int).values
    # DC count (composition table) and APC count consistency
    CP = pd.read_csv(f'{LR}/treg_local_composition.csv.gz', index_col=0).loc[T.index]   # cell-type FRACTIONS within 50 um (4 significant figures)
    assert (CP.nb_n.values == T.nb_n.values).all() and (np.round(CP[APC].sum(1).values * T.nb_n.values) == T.apc_n.values).all()
    W78 = pd.read_csv(f'{OUTB}/cn10_apc_per_treg_78.csv.gz', index_col=0).loc[T.index]
    assert (W78.nAPC.values == T.apc_n.values).all() and (W78.CN.values == T.CN.values).all()     # per-Treg APC count and CN of cn10_apc_abundance
    # ======================================================= MHC-II-weighted APC exposure (S_sum of cn10_apc_abundance)
    APCC = os.environ.get('APC_CACHE')
    if APCC and os.path.exists(APCC):
        _, TT = pd.read_pickle(APCC); TT = TT.loc[T.index, ['n_chk', 'sum_mhc2', 'n_DC']]
    else:
        import h5py
        from anndata.io import read_elem
        from scipy import sparse
        from scipy.spatial import cKDTree
        sa = pq.read_table(f'{RAW}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
        sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
        sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
        tb = []
        for s_, f in FOLDER.items():                           # compute stage of cn10_apc_abundance (APC MHC-II score, Treg -> APCs within 50 um)
            sub = sa[sa.section == s_]
            with h5py.File(f'{RAW}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
                ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
            nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
            keep = np.where(nm.isin(sub.index))[0]
            E = X[keep][:, [v.get_loc(g) for g in HLA]].toarray().astype(np.float32); del X; gc.collect()
            cs = sub.loc[nm[keep]]
            sd = E.std(0); Z = np.where(sd > 0, (E - E.mean(0)) / np.where(sd > 0, sd, 1), np.nan)
            am = np.isin(cs.ct.values, APC); mh = np.nanmean(Z[am], 1); axy = cs[['x', 'y']].values[am]
            tg = T[T.section == s_]; tx = sa.loc[tg.index, ['x', 'y']].values
            nb = cKDTree(axy).query_ball_point(tx, r=RAD)
            isdc = cs.ct.values[am] == 'DC'
            tb.append(pd.DataFrame({'n_chk': [len(i) for i in nb], 'sum_mhc2': [mh[i].sum() if len(i) else 0.0 for i in nb],
                                    'n_DC': [int(isdc[i].sum()) for i in nb]}, index=tg.index))
            print(s_, 'Tregs', len(tg), flush=True); del E, Z; gc.collect()
        TT = pd.concat(tb).loc[T.index]; del sa; gc.collect()
    assert (TT.n_chk.values == T.apc_n.values).all()
    T['mhc2_exp'] = TT.sum_mhc2.values
    T['dc_n'] = TT.n_DC.values.astype(int)                     # DC count within 50 um; = DC fraction of the composition table x nb_n
    assert (np.round(CP['DC'].values * T.nb_n.values) == T.dc_n.values).all()
    # ======================================================= tested CNs (rule of the CN forest plot) and labels
    nT = T.groupby('CN').size(); nS = T.groupby(['CN', 'section']).size().unstack(fill_value=0)
    TESTED = [c for c in range(1, 11) if nT.get(c, 0) >= MIN_TREG_TOTAL and (nS.loc[c] >= MIN_TREG_SEC).sum() >= MIN_SEC]
    FO = pd.read_csv(f'{OUTA}/fig4_cn_treg_state_forest.csv').set_index('CN')
    assert TESTED == [c for c in FO.index if FO.loc[c, 'tested']]          # = tested set of the CN forest plot
    assert 10 in TESTED and 4 in TESTED and 8 in TESTED
    L80 = pd.read_csv(f'{P4}/archive/supp_apc_abundance_by_neighbourhood.csv'); LAB = L80[L80.row_type == 'summary'].set_index('CN').CN_label.to_dict()
    assert LAB[10] == 'CN10: OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'
    OTH = [c for c in TESTED if c != 10]
    ROWS = []
    # ======================================================= 1. APC exposure per niche
    E1 = {}
    for c in TESTED:
        g = T[T.CN == c]; sm = g.groupby('section').agg(m=('apc_n', 'mean'), mf=('apc_frac', 'mean'), mn=('nb_n', 'mean'), md=('dc_n', 'mean'), n=('apc_n', 'size'))
        sm = sm[sm.n >= MIN_TREG_SEC]; E1[c] = sm
        m, lo, hi, ns = tci(sm.m)
        ROWS.append(dict(part='1 APC exposure per niche: mean of section means', CN=c, CN_label=LAB[c], n_tregs=len(g), n_sections=ns, mean_apc_n=m, ci_low=lo, ci_high=hi,
                         median_apc_n_pooled=g.apc_n.median(), mean_apc_frac=tci(sm.mf)[0], mean_nb_n=tci(sm.mn)[0], mean_dc_n=tci(sm.md)[0]))
        for s_, r in sm.iterrows():
            ROWS.append(dict(part='1 APC exposure per niche: section mean', CN=c, CN_label=LAB[c], section=s_, n_tregs=int(r.n), mean_apc_n=r.m, mean_apc_frac=r.mf, mean_nb_n=r.mn, mean_dc_n=r.md))
    d1 = T[T.CN.isin(TESTED)].assign(y=lambda x: np.log1p(x.apc_n), grp=lambda x: x.CN)
    f1, m1 = lmm('y ~ C(grp, Treatment(10))', d1, None)
    K = {c: f'C(grp, Treatment(10))[T.{c}]' for c in OTH}
    q1 = dict(zip(OTH, multipletests([f1.pvalues[K[c]] for c in OTH], method='fdr_bh')[1]))
    d1f = d1.assign(y=d1.apc_frac); f1f, _ = lmm('y ~ C(grp, Treatment(10))', d1f, None)
    q1f = dict(zip(OTH, multipletests([f1f.pvalues[K[c]] for c in OTH], method='fdr_bh')[1]))
    for c in OTH:
        b, se = f1.params[K[c]], f1.bse[K[c]]
        ROWS.append(dict(part='1 APC exposure per niche: mixed model log1p(APC count) ~ CN (ref CN10) + (1|section)', CN=c, CN_label=LAB[c], estimate=b, ci_low=b - 1.96 * se,
                         ci_high=b + 1.96 * se, ratio_vs_CN10=np.exp(b), p=f1.pvalues[K[c]], q=q1[c], optimiser=m1, n_tregs=len(d1)))
        b, se = f1f.params[K[c]], f1f.bse[K[c]]
        ROWS.append(dict(part='1 sensitivity: mixed model APC fraction within 50 um ~ CN (ref CN10) + (1|section)', CN=c, CN_label=LAB[c], estimate=b, ci_low=b - 1.96 * se,
                         ci_high=b + 1.96 * se, p=f1f.pvalues[K[c]], q=q1f[c], n_tregs=len(d1f)))
    print('1. APC count near Tregs (mean of section means):', {c: round(tci(E1[c].m)[0], 1) for c in TESTED}, '| q vs CN10:', {c: f'{q1[c]:.2g}' for c in OTH},
          '| ratio vs CN10:', {c: round(np.exp(f1.params[K[c]]), 2) for c in OTH}, flush=True)
    # ======================================================= 2. APC-matched Treg Suppression Score
    FOC = [10, 4, 8]
    d2 = T[T.CN.isin(FOC)].assign(grp=lambda x: x.CN, la=lambda x: np.log1p(x.apc_n), ld=lambda x: np.log1p(x.dc_n), lnb=lambda x: np.log1p(x.nb_n),
                                  lm=lambda x: np.sign(x.mhc2_exp) * np.log1p(np.abs(x.mhc2_exp)), y=lambda x: x.score)
    K2 = {c: f'C(grp, Treatment(10))[T.{c}]' for c in (4, 8)}
    MODELS = [('primary: + log1p(APC count) + library size', 'y ~ C(grp, Treatment(10)) + la + cz'),
              ('no exposure term: + library size only', 'y ~ C(grp, Treatment(10)) + cz'),
              ('sensitivity: + log1p(DC count) + library size', 'y ~ C(grp, Treatment(10)) + ld + cz'),
              ('sensitivity: + APC fraction within 50 um + library size', 'y ~ C(grp, Treatment(10)) + apc_frac + cz'),
              ('sensitivity: + APC quintile (87, factor) + library size', 'y ~ C(grp, Treatment(10)) + C(apc_q) + cz'),
              ('sensitivity: + log1p(APC count) + log1p(all cells within 50 um) + library size', 'y ~ C(grp, Treatment(10)) + la + lnb + cz'),
              ('sensitivity: + MHC-II-weighted APC exposure (78 S_sum) + library size', 'y ~ C(grp, Treatment(10)) + mhc2_exp + cz'),
              ('sensitivity: + signed log1p(MHC-II-weighted APC exposure) + library size', 'y ~ C(grp, Treatment(10)) + lm + cz')]
    P5, P95 = np.percentile(T.loc[T.CN == 10, 'apc_n'], [5, 95])
    MODELS.append((f'sensitivity: common support, APC count within CN10 5th-95th percentile ({P5:.0f}-{P95:.0f}), + log1p(APC count) + library size', 'y ~ C(grp, Treatment(10)) + la + cz'))
    M2 = {}
    for name, frm in MODELS:
        dm = d2[(d2.apc_n >= P5) & (d2.apc_n <= P95)] if name.startswith('sensitivity: common support') else d2
        f, meth = lmm(frm, dm, None); qq = multipletests([f.pvalues[K2[c]] for c in (4, 8)], method='fdr_bh')[1]
        M2[name] = {c: (f.params[K2[c]], f.bse[K2[c]], f.pvalues[K2[c]], qq[i]) for i, c in enumerate((4, 8))}
        expo = [t for t in f.fe_params.index if t not in ('Intercept', 'cz') and not t.startswith('C(grp')]
        for c in (4, 8):
            b, se, p, q = M2[name][c]
            ROWS.append(dict(part=f'2 model on CN4 / CN8 / CN10 Tregs: {name}', CN=c, CN_label=LAB[c], contrast=f'CN{c} minus CN10', estimate=b, ci_low=b - 1.96 * se,
                             ci_high=b + 1.96 * se, p=p, q=q, optimiser=meth, n_tregs=len(dm), n_tregs_CN10=int((dm.CN == 10).sum()), formula=frm.replace('y ~', 'score ~') + ' + (1|section)',
                             exposure_coef=';'.join(f'{t}={f.params[t]:.3g}' for t in expo)))
        print(f'2. {name}:', {c: (round(M2[name][c][0], 3), f'q={M2[name][c][3]:.2g}') for c in (4, 8)}, flush=True)
    # quintile x CN section means (focus CNs + others for faint lines)
    QM = {}
    for c in TESTED:
        for qn in range(1, 6):
            g = T[(T.CN == c) & (T.apc_q == qn)]; sm = g.groupby('section').score.agg(['mean', 'size']); sm = sm[sm['size'] >= MIN_Q]
            m, lo, hi, ns = tci(sm['mean']) if len(sm) else (np.nan, np.nan, np.nan, 0)
            QM[(c, qn)] = (m, lo, hi, ns, len(g))
            ROWS.append(dict(part='2 APC quintile x CN: mean of section means (>= 3 Tregs per section)', CN=c, CN_label=LAB[c], apc_quintile=qn, n_tregs=len(g), n_sections=ns,
                             mean_score=m, ci_low=lo, ci_high=hi, pooled_mean_score=g.score.mean(), mean_apc_n=g.apc_n.mean(), mean_apc_frac=g.apc_frac.mean()))
            if c in FOC:
                for s_, r in sm.iterrows():
                    ROWS.append(dict(part='2 APC quintile x CN: section mean', CN=c, CN_label=LAB[c], apc_quintile=qn, section=s_, n_tregs=int(r['size']), mean_score=r['mean']))
    below = {c: [qn for qn in range(1, 6) if QM[(10, qn)][0] < QM[(c, qn)][0]] for c in (4, 8)}
    print('2. quintile means CN10 / CN4 / CN8:', {qn: tuple(round(QM[(c, qn)][0], 2) for c in FOC) for qn in range(1, 6)},
          '| n sections:', {qn: tuple(QM[(c, qn)][3] for c in FOC) for qn in range(1, 6)}, '| CN10 below at quintiles', below, flush=True)
    QN10 = {qn: QM[(10, qn)][4] for qn in range(1, 6)}; print('   CN10 Tregs per quintile:', QN10, flush=True)
    # ======================================================= figure
    FS_T, FS_L, FS_TK, FS_X = 16, 15, 14, 13
    plt.rcParams.update({'font.size': FS_TK, 'font.family': 'DejaVu Sans'})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw=dict(width_ratios=[len(TESTED), 6.2], wspace=0.24))
    for ax in (a1, a2): ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=FS_TK)
    rng = np.random.default_rng(3)
    STY = {10: dict(c=HI, ms=11, lw=3.2, mk='o', mec=HIE), 4: dict(c=DK, ms=8, lw=2.0, mk='s', mec='white'), 8: dict(c=DK, ms=9, lw=2.0, mk='^', mec='white')}
    # ---- left: APC count near each Treg, per-section means
    X1 = {c: i for i, c in enumerate(TESTED)}
    a1.axvspan(X1[10] - 0.45, X1[10] + 0.45, color='#fbeaea', zorder=0, lw=0)
    ymax1 = max(E1[c].m.max() for c in TESTED)
    for c in TESTED:
        v = E1[c].m.values; m, lo, hi, ns = tci(v); col = HI if c == 10 else DK if c in (4, 8) else '#9a9a9a'
        a1.scatter(X1[c] + rng.uniform(-0.18, 0.18, len(v)), v, s=30, color=col, alpha=0.45, edgecolor='white', linewidth=0.5, zorder=2)
        a1.plot([X1[c]] * 2, [lo, hi], color=col, lw=2.6, zorder=3, solid_capstyle='butt')
        a1.plot([X1[c]], [m], STY.get(c, dict(mk='o'))['mk'] if c in STY else 'o', ms=11 if c == 10 else 9, color=col, mec=HIE if c == 10 else 'white', mew=1.2, zorder=4)
        if c != 10 and q1[c] < 0.05: a1.text(X1[c], max(v.max(), hi) + 0.02 * ymax1, stars(q1[c]), ha='center', va='bottom', fontsize=15, color='#252525')
        if c == 10: a1.text(X1[c], max(v.max(), hi) + 0.02 * ymax1, 'Ref.', ha='center', va='bottom', fontsize=13, color=HI)
    a1.set_ylim(0, ymax1 * 1.12); a1.set_xlim(-0.6, len(TESTED) - 0.4)
    a1.set_xticks(range(len(TESTED))); a1.set_xticklabels([LAB[c] for c in TESTED], rotation=40, ha='right', rotation_mode='anchor', fontsize=FS_X)
    for t_, c in zip(a1.get_xticklabels(), TESTED):
        t_.set_color(HI if c == 10 else '#1a1a1a' if c in (4, 8) else '#7a7a7a')
    a1.set_ylabel('APCs Within 50 µm of a Treg\n(Section Mean)', fontsize=FS_L)
    a1.set_title('APCs Near Each Treg', loc='left', fontsize=FS_T, pad=6)
    # ---- right: score by APC quintile
    a2.axhline(0, color='#9e9e9e', lw=0.9, zorder=1)
    for c in [c for c in TESTED if c not in FOC]:
        a2.plot(range(1, 6), [QM[(c, qn)][0] if QM[(c, qn)][3] >= MIN_QSEC else np.nan for qn in range(1, 6)], color=LT, lw=1.2, zorder=1)
    DX = {10: 0.0, 4: -0.13, 8: 0.13}
    for c in (4, 8, 10):
        st = STY[c]; xs = np.arange(1, 6) + DX[c]; ms_ = [QM[(c, qn)][0] if QM[(c, qn)][3] >= MIN_QSEC else np.nan for qn in range(1, 6)]
        a2.plot(xs, ms_, color=st['c'], lw=st['lw'], zorder=4 if c == 10 else 3)
        for x_, qn in zip(xs, range(1, 6)):
            m, lo, hi, ns, _ = QM[(c, qn)]
            if ns < MIN_QSEC: continue
            if np.isfinite(lo): a2.plot([x_, x_], [lo, hi], color=st['c'], lw=st['lw'] * 0.6, zorder=4 if c == 10 else 3, solid_capstyle='butt')
            a2.plot([x_], [m], st['mk'], ms=st['ms'], color=st['c'], mec=st['mec'], mew=1.2, zorder=5 if c == 10 else 4)
    a2.set_xticks(range(1, 6)); a2.set_xticklabels(['1\nFewest', '2', '3', '4', '5\nMost'], fontsize=FS_TK)
    a2.set_xlim(0.6, 5.4)
    a2.set_xlabel('Local APC Fraction Quintile (Within Section)', fontsize=FS_L)
    a2.set_ylabel('Treg Suppression Score\n(Section Mean, z)', fontsize=FS_L)
    a2.set_title('Treg Suppression Score by APC Exposure', loc='left', fontsize=FS_T, pad=6)
    f2 = lambda v: f'{v:+.2f}'.replace('-', '−')
    P = M2['primary: + log1p(APC count) + library size']
    txt = 'Adjusted for APC count and library size:\n' + '\n'.join(f'CN{c} − CN10  {f2(P[c][0])}  {stars(P[c][3]) or "n.s."}' for c in (4, 8))
    a2.text(0.03, 0.97, txt, transform=a2.transAxes, ha='left', va='top', fontsize=12.5, color='#252525')
    h = [Line2D([], [], color=HI, lw=3.2, marker='o', ms=10, mec=HIE, label=LAB[10]),
         Line2D([], [], color=DK, lw=2.0, marker='s', ms=8, mec='white', label=LAB[4]),
         Line2D([], [], color=DK, lw=2.0, marker='^', ms=9, mec='white', label=LAB[8]),
         Line2D([], [], color=LT, lw=1.2, label='Other CNs')]
    a2.legend(handles=h, loc='upper left', bbox_to_anchor=(-0.02, -0.25), frameon=False, fontsize=13, handlelength=2.2, borderaxespad=0.0, labelspacing=0.35, ncol=1)   # below the axis (free space beside the rotated CN labels)
    fig.subplots_adjust(left=0.07, right=0.995, top=0.93, bottom=0.12)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig); print('wrote', f'{STEM}.png', flush=True)
    # ======================================================= CSV
    OUT = pd.DataFrame(ROWS)
    lead = ['part', 'CN', 'CN_label', 'section', 'apc_quintile', 'contrast']
    OUT = OUT[[c for c in lead if c in OUT] + [c for c in OUT if c not in lead]].apply(sig3)
    with open(f'{STEM}.csv', 'w') as fh: OUT.to_csv(fh, index=False)
    print('wrote', f'{STEM}.csv', OUT.shape, flush=True)
    # ======================================================= legend (numbers from the data above)
    m10, m4, m8 = (tci(E1[c].m) for c in (10, 4, 8))
    r4, r8 = np.exp(f1.params[K[4]]), np.exp(f1.params[K[8]])
    U = M2['no exposure term: + library size only']; Dn = M2['sensitivity: + log1p(DC count) + library size']; Mx = M2['sensitivity: + MHC-II-weighted APC exposure (78 S_sum) + library size']
    qf = lambda q: '< 1e-300' if q < 1e-300 else f'= {q:.1g}'
    NB = {c: tci(E1[c].mn)[0] for c in FOC}
    TQ = {c: T[(T.CN == c) & (T.apc_q == 5)].apc_n.mean() for c in FOC}
    fe = lambda R, c: f'{f2(R[c][0])} (95% CI {f2(R[c][0] - 1.96 * R[c][1])} to {f2(R[c][0] + 1.96 * R[c][1])}; q = {R[c][3]:.2g})'
    leg = ('Treg-centred, APC-matched comparison of CN10 (OGN⁺RSPO3⁺/PLA2G2A⁺ fibroblast neighbourhood) Tregs with CN4 and CN8 Tregs, CD resections (CosMx 6K; '
       f'11 sections, 4 patients; {len(T):,} strict Tregs). APCs are DCs, B cells, SELENOP⁺ macrophages, CD9⁺ M2 macrophages and Mo-Macs. '
       'Treg Suppression Score: −D6 signature score, z-standardised within section. Neighbourhoods tested as in the Fig 4 forest plot '
       f'(at least {MIN_TREG_TOTAL} Tregs and at least {MIN_TREG_SEC} Tregs in at least {MIN_SEC} sections; CN6 and CN9 are not tested). '
       'APCs Near Each Treg: number of APCs within 50 µm of each Treg, averaged per section (dots; sections with at least 10 Tregs in the neighbourhood); '
       'large markers and lines, mean and 95% CI (t distribution) across sections. Stars: linear mixed model log(1 + APC count) ~ neighbourhood (reference CN10) '
       '+ (1 | section) on the Tregs of the tested neighbourhoods, Benjamini–Hochberg (BH) q over 7 contrasts. '
       'Treg Suppression Score by APC Exposure: Tregs grouped by the within-section quintile of the local APC fraction (APCs among all cells within 50 µm; '
       f'the APC quintiles of the CXCL12 analysis); mean and 95% CI across section means (at least {MIN_Q} Tregs per section, neighbourhood and quintile; points shown when at least {MIN_QSEC} sections qualify); '
       'grey lines, other tested neighbourhoods (means only). In-panel estimates: linear mixed model on CN4, CN8 and CN10 Tregs, score ~ neighbourhood '
       '(reference CN10) + log(1 + APC count) + library size + (1 | section), BH over 2 contrasts. * q < 0.05, ** q < 0.01, *** q < 0.001. '
       f'CN10 Tregs have fewer APCs nearby than CN4 and CN8 Tregs ({m10[0]:.1f} versus {m4[0]:.1f} and {m8[0]:.1f} APCs within 50 µm, means of section means; '
       f'model ratios {r4:.1f}× and {r8:.1f}×; q {qf(q1[4])} and {qf(q1[8])}). '
       f'Within an APC-fraction quintile, CN10 Tregs still have fewer APCs in absolute number (top quintile: {TQ[10]:.0f} versus {TQ[4]:.0f} and {TQ[8]:.0f} APCs), '
       f'because CN10 tissue is less cell-dense ({NB[10]:.0f} versus {NB[4]:.0f} and {NB[8]:.0f} cells of any type within 50 µm of a Treg, means of section means). '
       f'Adjusted for APC count and library size, the score difference from CN10 is {fe(P, 4)} for CN4 Tregs and {fe(P, 8)} for CN8 Tregs '
       f'(without the APC term: {f2(U[4][0])} and {f2(U[8][0])}; with DC count instead of APC count: {f2(Dn[4][0])} and {f2(Dn[8][0])}; with MHC-II-weighted '
       f'APC exposure: {f2(Mx[4][0])} and {f2(Mx[8][0])}). Per-section values, all model terms and sensitivity analyses are in the accompanying CSV.')
    for bad_ in ('activation score', 'quiescent', 'territory', 'serial', 'blunt'): assert bad_ not in leg.lower()
    with open(f'{STEM}_legend.txt', 'w') as fh: fh.write(leg + '\n')
    print(leg)


# ======================================================================================================================
# Section S9C_cn_apc_mhcii: APC abundance, per-APC MHC-II score and Treg Suppression Score per cellular neighbourhood (one panel)
# Panels: Supp 9C
# ======================================================================================================================
if _run('S9C_cn_apc_mhcii'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/81f_cn10_apc_one_panel.py'   # working folder of this step (relative paths resolve from here)
    """Supp 9C: APC abundance, per-APC MHC-II score and Treg Suppression Score of the cellular neighbourhoods in one panel (CD resections, CosMx 6K;
11 sections, 4 patients). One point per neighbourhood: x = APCs as % of all cells (mean across sections +- 95% CI), y = per-APC MHC-II score
(model-adjusted mean +- 95% CI), both read from the CSV of cn_apc_mhcii_table; point fill = the neighbourhood's Treg Suppression Score (Tregs inside
minus outside the CN, fig4_cn_treg_state_forest.csv of 12_neighbourhood_treg_score.py) on RdBu_r, +-0.4, with a colourbar. Second label line: forest
q-stars only ('Treg ***'); the APC % / MHC-II q values go to the CSV and legend.
Highlight: CN4, CN5 and CN10 in one neutral style (EMPH_CN = (4, 5, 10), EMPH_MARKER: bold label and a black-edged, slightly larger point; fill
unchanged). HILITE_CN10 = True switches to a CN10-only style (larger marker, thick black edge, red label) and CN10_BAND = True draws the CN10 MHC-II
95% CI as a grey band; both are off, and the legend gives the CN10 MHC-II CI instead (CN10 is the reference level of both models).
The APC counts near CN10 / CN4 / CN8 Tregs quoted in the legend are read from the CSV of cn10_apc_matched_score.
Labels: hand-placed offsets, checked programmatically against points, CI bars, leaders, band, colourbar and each other.
Outputs: /share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn10_apc_one_panel.{png,csv} + _legend.txt"""
    import os, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    SRC = f'{P4}/A_fig4_panels/supp_cn10_apc_abundance_mhcii.csv'
    STEM = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn10_apc_one_panel'
    SRCF = f'{P4}/A_fig4_panels/fig4_cn_treg_state_forest.csv'; SRC81D = f'{P4}/A_fig4_panels/fig4_cn10_apc_matched_treg_score.csv'
    CNS = [1, 2, 3, 4, 5, 7, 8, 10]
    HI, HIE, BAND = '#b2182b', '#000000', '#ececec'                              # CN10 label red / CN10 marker edge black / CN10 MHC-II CI band (grey: no clash with RdBu_r)
    VLIM = 0.4; CMAP = plt.get_cmap('RdBu_r'); NORM = matplotlib.colors.Normalize(-VLIM, VLIM)
    GREY, GREYT, GREYS = '#8c8c8c', '#666666', '#9a9a9a'           # other CNs: points + bars / label / significance line
    # highlight: CN4, CN5 and CN10 in one neutral style; CN10-only styling (HILITE_CN10) and the CN10 MHC-II CI band (CN10_BAND) off
    EMPH_CN = (4, 5, 10); HILITE_CN10 = False; CN10_BAND = False
    EMPH_MARKER = True; MS_EMPH, MEW_EMPH, MEC_EMPH = 14, 2.0, '#000000'           # highlighted point = black edge, slightly larger
    EM = lambda c: EMPH_MARKER and c in EMPH_CN and not (HILITE_CN10 and c == 10)
    H10 = lambda c: HILITE_CN10 and c == 10                                         # CN10-specific styling only when HILITE_CN10
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    # ======================================================= inputs (CSV of cn_apc_mhcii_table)
    S = pd.read_csv(SRC, low_memory=False)
    AP = S[(S.panel == 'APCs as % of Cells') & (S.row_type == 'summary')].copy()
    ALLR = AP[AP.CN.astype(str) == 'all'].iloc[0]; AP = AP[AP.CN.astype(str) != 'all']; AP['CN'] = AP.CN.astype(int); AP = AP.set_index('CN')
    MH = S[(S.panel == 'Per-APC MHC-II Score') & (S.row_type == 'model_adjusted_mean')].copy(); MH['CN'] = MH.CN.astype(int); MH = MH.set_index('CN')
    assert sorted(MH.index) == CNS and set(CNS) <= set(AP.index) and not {6, 9} & set(MH.index)
    T = pd.DataFrame({'CN_label': AP.loc[CNS, 'CN_label'], 'apc_pct': AP.loc[CNS, 'apc_pct_of_cells'], 'apc_ci_low': AP.loc[CNS, 'ci_low'],
                      'apc_ci_high': AP.loc[CNS, 'ci_high'], 'q_apc': AP.loc[CNS, 'q'], 'mhc2': MH.loc[CNS, 'adjusted_mean'],
                      'mhc2_ci_low': MH.loc[CNS, 'ci_low'], 'mhc2_ci_high': MH.loc[CNS, 'ci_high'], 'q_mhc2': MH.loc[CNS, 'q_78'],
                      'mhc2_diff_vs_CN10': MH.loc[CNS, 'estimate_vs_CN10']})
    assert (T.CN_label == MH.loc[CNS, 'CN_label']).all() and T.loc[10, 'CN_label'] == 'CN10: OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'
    for c in CNS:
        assert T.loc[c, 'apc_ci_low'] < T.loc[c, 'apc_pct'] < T.loc[c, 'apc_ci_high'] and T.loc[c, 'mhc2_ci_low'] < T.loc[c, 'mhc2'] < T.loc[c, 'mhc2_ci_high']
        if c != 10:
            assert stars(MH.loc[c, 'q']) == stars(MH.loc[c, 'q_78'])                         # refit stars == stars of cn10_apc_abundance
            assert abs((T.loc[c, 'mhc2'] - T.loc[10, 'mhc2']) - T.loc[c, 'mhc2_diff_vs_CN10']) < 2e-3   # bar differences = contrasts (3 sig figs)
    FO = pd.read_csv(SRCF).set_index('CN'); assert all(FO.loc[c, 'tested'] for c in CNS)
    T['treg_score'] = FO.loc[CNS, 'treg_suppression_score_z_in_minus_out'].values; T['treg_score_ci_low'] = FO.loc[CNS, 'ci_lo'].values
    T['treg_score_ci_high'] = FO.loc[CNS, 'ci_hi'].values; T['q_treg'] = FO.loc[CNS, 'mixed_q'].values
    assert all(FO.loc[c, 'CN_name'].startswith(f'CN{c}:') for c in CNS) and T.treg_score.abs().max() <= VLIM   # colour scale +-0.4 not clipped
    T['sig_line'] = [f'Treg {stars(q)}' if q < 0.05 else '' for q in T.q_treg]
    T['fill'] = [matplotlib.colors.to_hex(CMAP(NORM(v))) for v in T.treg_score]
    ALL_PCT = ALLR.apc_pct_of_cells
    # ======================================================= figure
    FS_T, FS_L, FS_TK, FS_P, FS_S, FS_P10 = 16, 16, 14, 12, 11, 15
    plt.rcParams.update({'font.size': FS_TK, 'font.family': 'DejaVu Sans'})
    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=FS_TK)
    ax.axhline(0, color='#9e9e9e', lw=0.9, zorder=1)
    ax.axvline(ALL_PCT, color='#8c8c8c', lw=1.3, ls='--', zorder=1)
    MS = {c: 17 if H10(c) else MS_EMPH if EM(c) else 12 for c in CNS}; COL = {c: '#262626' if H10(c) else GREY for c in CNS}   # CI bar colour
    LW = {c: 2.6 if H10(c) else 1.2 for c in CNS}
    band = ax.axhspan(T.loc[10, 'mhc2_ci_low'], T.loc[10, 'mhc2_ci_high'], color=BAND, lw=0, zorder=0) if CN10_BAND else None   # CN10 MHC-II 95% CI, full x range
    OBST = []                                                                       # artists text must not touch
    for c in CNS:
        r = T.loc[c]
        zc = 5 if H10(c) else 3
        hl, = ax.plot([r.apc_ci_low, r.apc_ci_high], [r.mhc2] * 2, color=COL[c], lw=LW[c], zorder=zc, solid_capstyle='butt')
        vl, = ax.plot([r.apc_pct] * 2, [r.mhc2_ci_low, r.mhc2_ci_high], color=COL[c], lw=LW[c], zorder=zc, solid_capstyle='butt')
        pt, = ax.plot([r.apc_pct], [r.mhc2], 'o', ms=MS[c], mfc=r.fill, mec=HIE if H10(c) else MEC_EMPH if EM(c) else '#8c8c8c', mew=3.0 if H10(c) else MEW_EMPH if EM(c) else 1.0,
                       zorder=zc + 1)   # EM = neutral highlight
        OBST += [(f'CN{c} x-CI', hl), (f'CN{c} y-CI', vl), (f'CN{c} point', pt)]
    ax.set_xlim(0, 40); ax.set_ylim(-0.315, 0.225)
    ax.set_xlabel('APCs (% of All Cells in the CN)', fontsize=FS_L); ax.set_ylabel('MHC-II Score per APC (z)', fontsize=FS_L)
    ax.set_title('APC Abundance, MHC-II per APC and Treg Suppression Score', fontsize=FS_T, loc='left', pad=6)
    # label anchors (data units): (x, y of the label's top edge, horizontal alignment, leader line?) - chosen by hand for free space, then checked below
    POS = {10: (14.7, 0.212, 'center', False),    # above CN10
           8: (31.8, 0.205, 'center', False),     # above CN8
           3: (0.6, 0.060, 'left', False),        # below CN3 (band between y = 0 and the upper CIs)
           4: (32.6, 0.057, 'center', False),     # below CN4
           7: (19.0, -0.058, 'left', False),      # right of CN7
           5: (18.0, -0.150, 'left', False),      # right of CN5
           2: (9.0, -0.214, 'left', True),        # bottom band, leader to CN2
           1: (0.6, -0.262, 'left', True)}        # bottom band, leader to CN1
    TX = {}
    BGW = dict(facecolor='white', edgecolor='none', pad=0.8)                     # text background: the dashed / zero lines pass behind labels
    LH = {c: 0.032 if H10(c) else 0.0235 for c in CNS}                               # data-units drop from label line to significance line
    for c in CNS:
        x, ytop, ha, lead = POS[c]
        t1 = ax.text(x, ytop, T.loc[c, 'CN_label'], ha=ha, va='top', fontsize=FS_P10 if H10(c) else FS_P, color=HI if H10(c) else GREYT, zorder=7, bbox=BGW,
                     fontweight='bold' if c in EMPH_CN else 'normal')   # bold labels for EMPH_CN
        t2 = ax.text(x, ytop - LH[c], T.loc[c, 'sig_line'], ha=ha, va='top', fontsize=13 if H10(c) else FS_S, color=HI if H10(c) else GREYS, zorder=7, bbox=BGW) if T.loc[c, 'sig_line'] else None
        TX[c] = (t1, t2)
        if lead:   # vertical leader from the lower end of the point's y-CI down to the label
            ax.plot([T.loc[c, 'apc_pct']] * 2, [T.loc[c, 'mhc2_ci_low'] - 0.006, ytop + 0.004], color='#b0b0b0', lw=0.8, zorder=2)
            OBST.append((f'CN{c} leader', ax.lines[-1]))
    at = ax.text(ALL_PCT + 0.4, -0.307, 'All Cells', ha='left', va='bottom', fontsize=FS_P, color='#6b6b6b', zorder=6)
    cax = ax.inset_axes([0.60, 0.075, 0.36, 0.032])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation='horizontal', ticks=[-0.4, -0.2, 0, 0.2, 0.4])
    cb.ax.tick_params(labelsize=12); cb.outline.set_linewidth(0.6)
    cb.ax.set_xticklabels(['−0.4', '−0.2', '0', '0.2', '0.4'])
    cb.ax.set_title('Treg Suppression Score (z)', fontsize=12.5, pad=4, color='#252525')
    fig.subplots_adjust(left=0.115, right=0.985, top=0.935, bottom=0.11)
    # ---- programmatic overlap check (display coordinates, 2 px margin)
    fig.canvas.draw(); rd = fig.canvas.get_renderer()
    def tb(t): return t.get_window_extent(rd).expanded(1.0, 1.0).padded(2)
    TEXTS = [(f'CN{c} label', TX[c][0]) for c in CNS] + [(f'CN{c} sig', TX[c][1]) for c in CNS if TX[c][1] is not None] + [('All Cells', at)]
    def line_bb(l):
        xy = ax.transData.transform(np.column_stack(l.get_data())); lw = l.get_linewidth() * fig.dpi / 72 / 2
        if l.get_marker() not in (None, 'None', ''): lw = max(lw, l.get_markersize() * fig.dpi / 72 / 2 + 1)
        return matplotlib.transforms.Bbox([xy.min(0) - lw, xy.max(0) + lw])
    bad = []
    for i, (n1, t1) in enumerate(TEXTS):
        b1 = tb(t1)
        assert ax.bbox.contains(b1.x0, b1.y0) and ax.bbox.contains(b1.x1, b1.y1), f'{n1} outside the axes'
        for n2, t2 in TEXTS[i + 1:]:
            if b1.overlaps(tb(t2)): bad.append((n1, n2))
        for n2, l in OBST:
            if b1.overlaps(line_bb(l)): bad.append((n1, n2))
        if band is not None and n1 != 'All Cells' and b1.overlaps(band.get_window_extent(rd)): bad.append((n1, 'CN10 band'))
        cbb = cax.get_tightbbox(rd).padded(2)
        if b1.overlaps(cbb): bad.append((n1, 'colourbar'))
    for n2, l in OBST:
        if line_bb(l).overlaps(cax.get_tightbbox(rd).padded(2)): bad.append(('colourbar', n2))
    for n1, l1 in [o for o in OBST if 'leader' in o[0]]:   # leaders must not cross other CNs' bars / points
        for n2, l2 in OBST:
            if n2.split()[0] != n1.split()[0] and line_bb(l1).overlaps(line_bb(l2)): bad.append((n1, n2))
    print('overlaps:', bad, flush=True)
    if bad: fig.savefig(os.environ.get('DEBUG_PNG', '/dev/null'), dpi=150, format='png')
    assert not bad, bad
    # style check - only EMPH_CN labels bold; highlighted points share one style, the others another; CN10 not special unless HILITE_CN10
    assert all((TX[c][0].get_fontweight() == 'bold') == (c in EMPH_CN) for c in CNS)
    assert (band is None) == (not CN10_BAND)
    if not HILITE_CN10:
        P_ = {c: [l for n, l in OBST if n == f'CN{c} point'][0] for c in CNS}
        for grp in ([c for c in CNS if EM(c)], [c for c in CNS if not EM(c)]):
            assert len({(P_[c].get_markersize(), P_[c].get_markeredgecolor(), P_[c].get_markeredgewidth()) for c in grp}) <= 1, grp
        assert all(P_[c].get_markersize() == MS_EMPH and P_[c].get_markeredgecolor() == MEC_EMPH for c in CNS if EM(c))
        assert len({(TX[c][0].get_color(), TX[c][0].get_fontsize()) for c in CNS}) == 1 and all(TX[c][1] is None or TX[c][1].get_color() == GREYS for c in CNS)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig); print('wrote', f'{STEM}.png', flush=True)
    # ======================================================= CSV (plotted values)
    OUT = T.reset_index().rename(columns={'index': 'CN'}).assign(all_cells_apc_pct=ALL_PCT,
           metric='x = APCs as % of all cells in the CN (mean over sections, 95% CI t; q_apc = mixed model vs CN10, BH over 9); y = per-APC MHC-II score, '
              'model-adjusted mean (mixed model; 95% CI; q_mhc2 = q vs CN10, BH over 7); fill = treg_score: Treg Suppression Score, Tregs in minus outside the CN '
              '(CN forest LMM, q_treg BH over 8 tested CNs), RdBu_r +-0.4')
    OUT['stars_treg'] = [stars(q) for q in OUT.q_treg]
    OUT['stars_apc'] = [stars(q) if c != 10 else '' for c, q in zip(OUT.CN, OUT.q_apc)]; OUT['stars_mhc2'] = [stars(q) if c != 10 else '' for c, q in zip(OUT.CN, OUT.q_mhc2)]
    OUT = OUT.apply(lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v))
    with open(f'{STEM}.csv', 'w') as fh: OUT.to_csv(fh, index=False)
    print('wrote', f'{STEM}.csv', OUT.shape, flush=True)
    # ======================================================= legend (numbers from the cn_apc_mhcii_table CSV)
    r10 = AP.loc[10]
    NAME = {'Macrophage': 'macrophages', 'DC': 'DCs', 'B cell': 'B cells'}
    MIX10 = ', '.join(f'{v:.0f}% {NAME[a]}' for v, a in sorted(((r10[f'share_of_APCs_{a}'], a) for a in NAME), reverse=True))
    f2 = lambda v: '0.00' if abs(v) < 0.005 else f'{v:+.2f}'
    w = MH.iloc[0]; WMIXT = f'DCs {100 * w.weight_DC:.0f}%, B cells {100 * w.weight_B_cell:.0f}%, macrophages {100 * w.weight_Macrophage:.0f}%'
    cl = lambda L: ', '.join(f'CN{c}' for c in L[:-1]) + (' and ' if len(L) > 1 else '') + f'CN{L[-1]}'
    t10, t4, t8 = T.loc[10], T.loc[4], T.loc[8]
    assert T.loc[4, 'q_mhc2'] >= 0.05 and T.loc[8, 'q_mhc2'] >= 0.05 and T.loc[4, 'q_apc'] < 0.05 and T.loc[8, 'q_apc'] < 0.05
    assert 0.42 <= t10.apc_pct / ((t4.apc_pct + t8.apc_pct) / 2) < 0.6                  # 'about half'
    INB = [c for c in CNS if c != 10 and t10.mhc2_ci_low <= T.loc[c, 'mhc2'] <= t10.mhc2_ci_high]
    low_m = [c for c in CNS if c != 10 and T.loc[c, 'q_mhc2'] < 0.05 and T.loc[c, 'mhc2_diff_vs_CN10'] < 0]
    low_a = [c for c in CNS if c != 10 and T.loc[c, 'q_apc'] < 0.05 and T.loc[c, 'apc_pct'] < t10.apc_pct]
    A1 = pd.read_csv(SRC81D); A1 = A1[A1.part == '1 APC exposure per niche: mean of section means'].set_index('CN')
    na = {c: A1.loc[c, 'mean_apc_n'] for c in (10, 4, 8)}; assert na[10] < na[4] and na[10] < na[8]
    qs = lambda c, k: ('q < 1e-300' if T.loc[c, k] < 1e-300 else f'q = {T.loc[c, k]:.2g}')
    APCQ = '; '.join(f'CN{c} {qs(c, "q_apc")}' for c in CNS if c != 10); MHCQ = '; '.join(f'CN{c} {qs(c, "q_mhc2")}' for c in CNS if c != 10)
    TRQ = ', '.join(f'CN{c} {f2(T.loc[c, "treg_score"])}{(" " + stars(T.loc[c, "q_treg"])) if T.loc[c, "q_treg"] < 0.05 else ""}' for c in CNS)
    high_m = [c for c in CNS if c != 10 and T.loc[c, 'mhc2'] >= t10.mhc2_ci_low]
    assert T.loc[10, 'treg_score'] == T.loc[high_m + [10], 'treg_score'].min()     # CN10 = lowest Treg score among the high-MHC-II CNs
    leg = ('APC abundance, per-APC MHC-II expression and Treg Suppression Score of the cellular neighbourhoods (Fig 4) of the CD resections (CosMx 6K; 11 sections, '
       '4 patients). Each point is one neighbourhood'
           + ('; CN10, the OGN⁺RSPO3⁺/PLA2G2A⁺ fibroblast neighbourhood, has the larger marker with a black edge and a red label' if HILITE_CN10 else '')
           + ('; ' + ' and '.join(f'CN{c}' for c in EMPH_CN) + ' labels are shown in bold' if EMPH_CN and not EMPH_MARKER else '') + '. '
           + ('CN4, CN5 and CN10 (highlighted) differ in APC abundance and per-APC MHC-II, in line with their Treg Suppression Scores. '
              if EMPH_MARKER and tuple(EMPH_CN) == (4, 5, 10) else '')
           +               'APCs are DCs, B cells and macrophages (SELENOP⁺ macrophages, CD9⁺ M2 macrophages and Mo-Macs). '
       'x: APCs as a percentage of all cells of the neighbourhood, mean and 95% CI (t distribution) across sections (section × neighbourhood combinations with at least '
       f'100 cells); the dashed line is the all-cell APC percentage ({ALL_PCT:.1f}%, mean of sections). '
       'y: per-APC MHC-II score (each APC\'s mean within-section z over the 9 HLA ligand genes of the Fig 4D antigen-presentation module, z-standardised within section '
       'over APCs, so 0 = average APC of the section; APCs grouped by their own neighbourhood), shown as model-adjusted means with 95% CI from the linear mixed model '
       'score ~ neighbourhood (reference CN10) + APC group + library size + (1 | section): the CN10 intercept plus each neighbourhood\'s effect, at average library '
       f'size and the overall APC-group mix ({WMIXT}). '
           + (f'The grey band spans the 95% CI of the CN10 MHC-II score ({f2(t10.mhc2_ci_low)} to {f2(t10.mhc2_ci_high)}). ' if CN10_BAND else
              f'CN10 is the reference level of both models (95% CI of its MHC-II score {f2(t10.mhc2_ci_low)} to {f2(t10.mhc2_ci_high)}). ') + f''
       'Point colour: Treg Suppression Score (−D6 signature score, z-standardised within section) of the Tregs inside minus outside the neighbourhood, from the '
       'linear mixed model of the Fig 4 forest plot (section random intercept); colour scale ±0.4. Second label line: its Benjamini–Hochberg q over the 8 tested '
       'neighbourhoods, shown only when significant (* q < 0.05, ** q < 0.01, *** q < 0.001). CN6 and CN9 are omitted because they are not in the per-APC '
       'MHC-II model, which, as in Fig 4D, excludes neighbourhoods with too few Tregs. '
       f'Comparisons with CN10 (not marked in the panel): APC %, linear mixed model APC % ~ neighbourhood (reference CN10) + (1 | section), BH over 9 contrasts '
       f'({APCQ}); MHC-II, the model above, BH over 7 contrasts ({MHCQ}). '
       f'CN10 APCs make up {t10.apc_pct:.1f}% of its cells (95% CI {t10.apc_ci_low:.1f}–{t10.apc_ci_high:.1f}%; {MIX10}), about half of CN4 ({t4.apc_pct:.1f}%) and '
       f'CN8 ({t8.apc_pct:.1f}%). Their adjusted MHC-II score ({f2(t10.mhc2)}) equals that of CN4 ({f2(t4.mhc2)}) and CN8 ({f2(t8.mhc2)}), yet CN10 Tregs have the '
       f'lowest Treg Suppression Score of the high-MHC-II neighbourhoods (Treg Suppression Score, in minus outside: {TRQ}). '
       'The panel describes neighbourhoods, not matched APC exposure: around the Tregs themselves, CN10 Tregs have fewer APCs within 50 µm than CN4 and CN8 Tregs '
       f'({na[10]:.0f} versus {na[4]:.0f} and {na[8]:.0f}, means of section means; fig4_cn10_apc_matched_treg_score). Values are in the accompanying CSV.')
    assert 4 in INB and 8 in INB and all(T.loc[c, 'mhc2_ci_high'] < t10.mhc2_ci_low for c in low_m), (INB, low_m)
    for bad_ in ('activation', 'territory', 'quiescent', 'serial', 'blunt'): assert bad_ not in leg.lower()
    with open(f'{STEM}_legend.txt', 'w') as fh: fh.write(leg + '\n')
    print(leg)
