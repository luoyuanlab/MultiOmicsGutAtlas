#!/usr/bin/env python
"""CD CosMx 6K resections: cellular neighbourhoods (CN) and per-Treg local-environment tables, read by 12-17.

Inputs   /share/fsmresfiles/UC/AtoMx/CD_surgical_resection/ (combined/all_samples.parquet, combined/raw_counts_strict_tregs.csv,
         per-run Processed/norm_log1p.h5ad); WTX norm_anno_tregbin.h5ad (02); D_shared_inputs/treg_definition_comparison_per_cell.csv (10)
Outputs  .../B_resection_analyses/cellular_neighborhoods_v2/ and .../B_resection_analyses/ligand_receptor_state/
Sections cellular_neighbourhoods -> treg_local_ligand_receptor -> treg_local_composition -> cxcl12_expected_producer_field
Run      python <file>, or RUN_SECTIONS=<id>[,<id>] python <file>
"""
import os, gc, warnings
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
from anndata.io import read_elem
from scipy import sparse, stats
from scipy.spatial import cKDTree
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
warnings.filterwarnings('ignore')
_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'; LR = f'{P4}/B_resection_analyses/ligand_receptor_state'
R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; WTX = '/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged'
UM, RAD = 0.12028, 50
APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']
FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
          '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
          '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
          '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
          '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
          '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}


def load_cells(cols=('cell_type_general',)):
    """Resection cells with section, global coordinates in um and the requested annotation columns."""
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'patient_sec', *cols, '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px.values * UM, sa.CenterY_global_px.values * UM
    return sa


def resection_tregs():
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0)
    return PC[PC.cohort == 'resection']


def norm_expr(sec, genes):
    """log1p-normalised expression of `genes` for every cell of one resection run."""
    with h5py.File(f'{R}/{FOLDER[sec]}/Processed/norm_log1p.h5ad', 'r') as h:
        ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
    X = X.tocsc() if sparse.issparse(X) else sparse.csc_matrix(X)
    gl = [g for g in genes if g in v]
    return pd.DataFrame(X[:, [v.get_loc(g) for g in gl]].toarray(), index=pd.Index(sec + '_' + ob.index.astype(str)), columns=gl)


def network_anno(df):
    """Cell-type grouping of the LR modules: Tregs one type, neutrophils, Tfh and CD8 Trm split, SELENOP+ Fib dropped."""
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
            elif 'SELENOP' in j: out.append('DROP')
            else: out.append(j)
        elif i == 'T':
            if j == 'Treg': out.append('Treg')
            elif j == 'Exhausted CD4 Tfh': out.append('Exhausted CD4 Tfh')
            elif 'Tfh' in j: out.append('Tfh')
            elif j == 'CD8 Trm': out.append('CD8 GZMK+ Trm' if 'GZMK+' in k else 'CD8 GZMK- Trm')
            else: out.append(j)
        elif i == 'Granulocyte': out.append('Neutrophil' if j == 'Neutrophil' else 'Granulocyte')
        elif i.lower() == 'cycling': out.append('Cycling')
        else: out.append(i)
    return out


# ======================================================================================================================
# Section cellular_neighbourhoods: 20-nearest-cell windows, MiniBatchKMeans k = 10 (k = 8, 12 sensitivity), Treg state per CN (CMH, BH)
# ======================================================================================================================
if _run('cellular_neighbourhoods'):
    os.makedirs(f'{CND}/discovery_biopsy', exist_ok=True)
    KNN, KS, K0 = 20, [8, 10, 12], 10

    def windows(df, cats):
        """Cell-type composition of the KNN nearest cells (self excluded), neighbours within the same section."""
        code = pd.Categorical(df['ct'], categories=cats).codes; rr, cc = [], []
        for _, idx in df.groupby('section').indices.items():
            xy = df[['x', 'y']].values[idx]
            _, nn = cKDTree(xy).query(xy, k=KNN + 1)
            rr.append(np.repeat(idx, KNN)); cc.append(code[idx[nn[:, 1:]].ravel()])
        rr, cc = np.concatenate(rr), np.concatenate(cc)
        return sparse.csr_matrix((np.ones(len(rr), np.float32), (rr, cc)), shape=(len(df), len(cats))).toarray() / KNN

    def name_cns(cts, lab, k, cats, overall):
        """log2 enrichment of each cell type in each CN vs all cells; CN name = up to 3 types with log2 enrichment > 0.5."""
        cell = pd.crosstab(lab, cts.values, normalize='index').reindex(index=range(k), columns=cats, fill_value=0)
        enr = np.log2((cell + 1e-4) / (overall + 1e-4)); names = {}
        for c in range(k):
            top = enr.loc[c].sort_values(ascending=False); top = [t for t in top.index if top[t] > 0.5][:3] or [top.index[0]]
            names[c] = f'CN{c + 1}: ' + ' / '.join(top)
        return enr, names

    # resection CNs
    sa = load_cells(['fov', 'cell_category', 'cell_type_general', 'cell_type_short'])
    sa['fov'] = sa['fov'].astype(str); sa['ct'] = network_anno(sa); sa = sa[sa.ct != 'DROP'].copy()
    CT = sorted(sa['ct'].unique()); W = windows(sa, CT)
    LAB = {k: MiniBatchKMeans(n_clusters=k, random_state=0, batch_size=20000, n_init=5, max_iter=300).fit(W).labels_ for k in KS}
    ari = {k: adjusted_rand_score(LAB[K0], LAB[k]) for k in KS if k != K0}
    overall = pd.Series(np.bincount(pd.Categorical(sa['ct'], categories=CT).codes, minlength=len(CT)) / len(sa), index=CT)
    TAB = {k: name_cns(sa['ct'], LAB[k], k, CT, overall) for k in KS}
    enr, NAMES = TAB[K0]; sa['CN'] = LAB[K0]
    enr.rename(index=NAMES).to_csv(f'{CND}/cn_k10_celltype_log2_enrichment.csv', float_format='%.3g')

    # Treg state (D6 tertile) per CN: 2 x 2 per section, CMH odds ratio across sections, BH across CNs
    PC = resection_tregs()

    def treg_cn_test(dname, k):
        t = pd.DataFrame({'section': sa.loc[PC.index, 'section'].values, 'state': PC[dname].values, 'CN': LAB[k][sa.index.get_indexer(PC.index)]})
        t = t[t.state.isin(['Early', 'Late'])]; out = []
        for c in range(k):
            tabs, sgn = [], []
            for _, g in t.groupby('section'):
                e, l = g[g.state == 'Early'], g[g.state == 'Late']
                a, b = int((e.CN == c).sum()), int((e.CN != c).sum()); cc, d = int((l.CN == c).sum()), int((l.CN != c).sum())
                if a + cc == 0 or b + d == 0: continue
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

    RES = pd.concat([treg_cn_test('D6_tertile', k) for k in KS], ignore_index=True)
    RES['CN_name'] = [TAB[k][1][c] for k, c in zip(RES.k, RES.CN)]
    dens = sa.groupby('CN')['ct'].apply(lambda x: (x == 'Treg').mean())
    RES = RES.merge(pd.DataFrame({'k': K0, 'CN': dens.index, 'treg_fraction_of_cells': dens.values}), on=['k', 'CN'], how='left')
    RES.to_csv(f'{CND}/cn_treg_state_enrichment_all.csv', index=False, float_format='%.4g')
    pd.DataFrame({'section': sa.section.values, 'fov': sa.fov.values, 'x_um': sa.x.values, 'y_um': sa.y.values, 'cell_type': sa.ct.values,
                  'CN_k10': sa.CN.values + 1, 'CN_k8': LAB[8] + 1, 'CN_k12': LAB[12] + 1}, index=sa.index).to_parquet(f'{CND}/cn_per_cell.parquet')
    pd.DataFrame({'k': list(ari), 'ARI_vs_k10': list(ari.values())}).to_csv(f'{CND}/cn_k_sensitivity_ari.csv', index=False, float_format='%.3g')
    del W; gc.collect()

    # WTX biopsy CNs on the shared vocabulary (per-cell table used by the Fig 5A maps of 14)
    with h5py.File(f'{WTX}/h5ad_obj/norm_anno_tregbin.h5ad', 'r') as h:
        bo = read_elem(h['obs']); bxy = h['obsm/spatial_fov'][:]
    bo = bo.set_index('cell', drop=False); bo.index = bo.index.astype(str)
    bo['section'] = bo['slide'].astype(str); bo['fov'] = bo['fov'].astype(str)
    bo['x'], bo['y'] = bxy[:, 0] * UM, bxy[:, 1] * UM
    bo['ct'] = network_anno(bo); bo = bo[bo.ct != 'DROP'].copy()
    bo['condition_lab'] = bo['patient'].astype(str).map({'I0368': 'Inflamed CD', 'I0476': 'Inflamed CD', 'I0453': 'Non-inflamed CD', 'I0477': 'Healthy control'})
    WB = windows(bo, sorted(set(bo['ct']) | set(sa['ct'])))
    LB = {k: MiniBatchKMeans(n_clusters=k, random_state=0, batch_size=4096, n_init=10, max_iter=500).fit(WB).labels_ for k in KS}
    bo[['section', 'patient', 'condition_lab', 'fov', 'x', 'y', 'ct']].assign(CN_k10=LB[K0] + 1, CN_k8=LB[8] + 1, CN_k12=LB[12] + 1).to_csv(
        f'{CND}/discovery_biopsy/cn_per_cell_biopsy.csv.gz', float_format='%.5g')


# ======================================================================================================================
# Section treg_local_ligand_receptor: per-Treg local ligand availability (50 um), own receptor expression, local APC density, score
# ======================================================================================================================
if _run('treg_local_ligand_receptor'):
    os.makedirs(LR, exist_ok=True)
    LIG = ['CXCL12', 'CXCL14', 'IL6', 'TNF', 'CCL22']
    REC = ['CXCR4', 'CXCR3', 'IL6R', 'IL6ST', 'TNFRSF1B', 'CCR7']
    sa = load_cells(); sa['ct'] = sa.cell_type_general.astype(str); sa['apc'] = sa.ct.isin(APC)
    pcr = resection_tregs()
    T = sa.loc[pcr.index, ['section', 'x', 'y']].copy()
    T['act_z'] = -pcr.signature_score_D6.values; T['act_z'] = T.groupby('section').act_z.transform(lambda v: (v - v.mean()) / v.std())   # Treg Suppression Score
    for g in LIG: T['lig_' + g] = np.nan
    for g in REC: T['rec_' + g] = np.nan
    T['apc_n'] = np.nan; T['nb_n'] = np.nan
    for s_ in FOLDER:
        sub = sa[sa.section == s_]
        E = norm_expr(s_, LIG + REC); E = E[E.index.isin(sub.index)]; cs = sub.loc[E.index]
        Z = (E - E.mean()) / E.std().replace(0, np.nan)                     # within-section z
        tg = T[T.section == s_]
        if not len(tg): continue
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(tg[['x', 'y']].values, r=RAD)
        apcv = cs.apc.values.astype(float)
        T.loc[tg.index, 'nb_n'] = [len(i) for i in nb]; T.loc[tg.index, 'apc_n'] = [apcv[i].sum() if len(i) else np.nan for i in nb]
        for g in [x for x in LIG if x in Z.columns]:
            vals = Z[g].values; T.loc[tg.index, 'lig_' + g] = [vals[i].mean() if len(i) >= 5 else np.nan for i in nb]
        common = tg.index.intersection(Z.index)
        for g in [x for x in REC if x in Z.columns]: T.loc[common, 'rec_' + g] = Z.loc[common, g].values
    T['apc_frac'] = T.apc_n / T.nb_n.replace(0, np.nan)
    T.to_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', float_format='%.4g')


# ======================================================================================================================
# Section treg_local_composition: per-Treg fraction of each cell type among all cells within 50 um
# ======================================================================================================================
if _run('treg_local_composition'):
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    idx = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0, usecols=[0]).index
    T = T.loc[idx.intersection(T.index)]                                 # Treg order of the raw-count table
    sa = load_cells(); sa['ct'] = sa.cell_type_general.astype(str)
    CT = sorted(sa.ct.unique()); comp = []
    for s_, tg in T.groupby('section'):
        cs = sa[sa.section == s_]; code = pd.Categorical(cs.ct, categories=CT).codes
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
        M = np.zeros((len(tg), len(CT)))
        for k, i in enumerate(nb): M[k] = np.bincount(code[i], minlength=len(CT))
        n = M.sum(1); C_ = pd.DataFrame(M / n[:, None], index=tg.index, columns=CT); C_.insert(0, 'nb_n', n); C_.insert(0, 'section', s_)
        comp.append(C_)
    pd.concat(comp).loc[T.index].to_csv(f'{LR}/treg_local_composition.csv.gz', float_format='%.4g')


# ======================================================================================================================
# Section cxcl12_expected_producer_field: local CXCL12 over expected stromal producers per Treg; CXCR4 and CXCL12-High Region models
# ======================================================================================================================
if _run('cxcl12_expected_producer_field'):
    MIN_PROD, MIN_SEC = 3, 60
    EXPECTED_PRODUCERS = ['OGN+RSPO3+ Fib', 'ADAMDEC1+ Fib', 'Activated Fib', 'Crypt Top Fib', 'PLA2G2A+ ECM Fib', 'Myofibroblast', 'Pericyte', 'Endothelial']
    SEND = {'dc': ['DC'], 'mac': ['SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac'], 'b': ['B cell'], 'plasma': ['Plasma']}
    NAME = {'dc': 'DC', 'mac': 'Macrophage', 'b': 'B cell', 'plasma': 'Plasma'}
    NOPROD = 'Fewer than 3 Expected Producers within 50 um'
    SS = 'Treg Suppression Score'; XC = 'lig_CXCL12_exp'
    sig3 = lambda v: float(f'{v:.3g}') if pd.notna(v) else v

    def r3(df):
        df = df.copy()
        for c in df.columns:
            if df[c].dtype.kind == 'f' and not (df[c].dropna() % 1 == 0).all(): df[c] = df[c].map(sig3)   # counts are kept exact
        return df

    def wcsv(df, path, **kw):   # NFS share: plain python open('w')
        with open(path, 'w') as fh: df.to_csv(fh, **kw)

    sa = load_cells(); sa['ct'] = sa.cell_type_general.astype(str); sa = sa[['section', 'x', 'y', 'ct']]
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    for c in [XC, 'n_exp_producers_50um'] + [k + '_n' for k in SEND]: T[c] = np.nan
    for s_ in FOLDER:
        sub = sa[sa.section == s_]
        e = norm_expr(s_, ['CXCL12'])['CXCL12']; e = e[e.index.isin(sub.index)]
        cs = sub.loc[e.index]; z = ((e.values - e.mean()) / e.std())       # within-section z
        tg = T[T.section == s_]
        if not len(tg): continue
        nb = [np.asarray(i, dtype=int) for i in cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)]
        pm = cs.ct.isin(EXPECTED_PRODUCERS).values; ctv = cs.ct.values
        npd = np.array([pm[i].sum() for i in nb])
        T.loc[tg.index, 'n_exp_producers_50um'] = npd
        T.loc[tg.index, XC] = [z[i[pm[i]]].mean() if k >= MIN_PROD else np.nan for i, k in zip(nb, npd)]
        for k, types in SEND.items():
            m_ = np.isin(ctv, types); T.loc[tg.index, k + '_n'] = [m_[i].sum() for i in nb]
    NEW = T[[XC, 'n_exp_producers_50um']].copy(); NEW['n_exp_producers_50um'] = NEW.n_exp_producers_50um.astype(int)
    with open(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', 'wb') as fh: NEW.to_csv(fh, float_format='%.4g', compression='gzip')
    T[XC] = pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0)[XC].reindex(T.index)   # analyse the stored (4 sig fig) values

    # covariates: library size (cz), sender fractions, APC-fraction quintile, CXCL12 tertile (rank first: ties at zero)
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()
    zs = lambda v: (v - v.mean()) / v.std()
    T['cz'] = T.groupby('section').cz.transform(zs)
    for k in SEND: T[k + '_frac'] = T[k + '_n'] / T.nb_n.replace(0, np.nan)
    for k in SEND: T[k + '_z'] = T.groupby('section')[k + '_frac'].transform(zs)
    T['apc_q'] = T.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    T['t_exp'] = T.groupby('section')[XC].transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))

    def lmm(d, f, terms):
        m = smf.mixedlm(f, d, groups=d.section).fit(reml=True, method='lbfgs')
        return {t: dict(coef=m.params[t], se=m.bse[t], p=m.pvalues[t], n=len(d)) for t in terms}

    # (a) CXCR4 ~ local CXCL12
    rows = []
    for lab, f in [('CXCR4 ~ CXCL12', 'y ~ x'), (f'CXCR4 ~ CXCL12 + {SS}', 'y ~ x + a'), (f'CXCR4 ~ CXCL12 + {SS} + library size', 'y ~ x + a + cz')]:
        d = T.assign(y=T.rec_CXCR4, x=T[XC], a=T.act_z)[['section', 'y', 'x', 'a', 'cz']].dropna()
        r = lmm(d, f, ['x'])['x']; sd = (d.x - d.groupby('section').x.transform('mean')).std()
        rows.append(dict(model=lab, cxcl12=XC, **r, sd_x_within_section=sd, coef_per_sd_x=r['coef'] * sd))
    wcsv(r3(pd.DataFrame(rows)), f'{LR}/a_cxcr4_cxcl12_exp.csv', index=False)

    # (c) CXCL12-High vs -Low Region (within-section tertiles) x APC quintile
    def per_quintile(tcol, lowset=('low',)):
        sub = T[T[tcol].isin(list(lowset) + ['high']) & T.apc_q.notna()]; out = []
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

    def senders(tcol, lowset=('low',)):
        """Single-sender interactions act_z ~ hi * s_z + cz (BH over the 4 senders), per-section negative counts."""
        d = T[T[tcol].isin(list(lowset) + ['high'])].copy(); d['hi'] = (d[tcol] == 'high').astype(float); d = d.rename(columns={'act_z': 'yy'})
        out = []
        for k in SEND:
            dd = d[['section', 'yy', 'hi', k + '_z', 'cz']].dropna().rename(columns={k + '_z': 's_z', 'yy': 'y'})
            m = lmm(dd, 'y ~ hi * s_z + cz', ['s_z', 'hi:s_z', 'hi'])
            neg = tested = 0
            for _, g in dd.groupby('section'):
                if len(g) < MIN_SEC or g.s_z.std() == 0 or g.hi.nunique() < 2: continue
                tested += 1; neg += int(smf.ols('y ~ hi * s_z + cz', g).fit().params['hi:s_z'] < 0)
            out.append(dict(sender=NAME[k], m=m, neg=neg, tested=tested))
        qi = multipletests([o['m']['hi:s_z']['p'] for o in out], method='fdr_bh')[1]
        for o, b in zip(out, qi): o['int_q'] = b
        return out

    wcsv(r3(per_quintile('t_exp')), f'{LR}/c_apc_matched_cxcl12_high_vs_low_exp.csv', index=False)
    Po = pooled('t_exp')
    wcsv(pd.DataFrame([dict(estimate_sd=Po['coef'], se=Po['se'], p=Po['p'], n=Po['n'],
                            model=f'{SS} ~ CXCL12-High Region (expected producers) + C(APC quintile) + (1|section)')]).pipe(r3), f'{LR}/c_pooled_apc_stratified_exp.csv', index=False)
    IXr = interaction('t_exp')
    wcsv(pd.DataFrame([dict(term=t, **IXr[t]) for t in ['hi', 'apcz', 'hi:apcz', 'cz']]).pipe(r3), f'{LR}/c_interaction_apc_by_cxcl12_countadj_exp.csv', index=False)

    # sensitivity: Tregs with < 3 expected producers as their own group (S1) or merged into CXCL12-Low (S2)
    T['grp'] = T.t_exp.astype(object).where(T.t_exp.notna(), NOPROD)
    T['grp_low'] = T.grp.replace({NOPROD: 'low'})
    SEN = []
    desc = T.groupby('grp').agg(n_tregs=('act_z', 'size'), mean_suppression_z=('act_z', 'mean'), mean_CXCR4_z=('rec_CXCR4', 'mean'), mean_apc_frac=('apc_frac', 'mean'),
                                mean_neighbours=('nb_n', 'mean'), mean_exp_producers=('n_exp_producers_50um', 'mean'), mean_lig_CXCL12_allcell=('lig_CXCL12', 'mean'))
    for g, r in desc.reindex([NOPROD, 'low', 'mid', 'high']).iterrows(): SEN.append(dict(analysis='S1 descriptives by group', group=g, **r.to_dict()))
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
    SEN.append(dict(analysis=f'S2 pooled {SS} ~ CXCL12-High + C(APC quintile), <3-producer merged into low', group='hi', **pooled('grp_low')))
    for _, r in per_quintile('grp_low').iterrows():
        SEN.append(dict(analysis='S2 per APC quintile, <3-producer merged into low (delta SD, BH over 5)', group=f'APC quintile {int(r.apc_quintile)}',
                        coef=r.delta_sd, p=r.p, q=r.q, n=int(r.n_high + r.n_low)))
    for o in senders('grp_low'):
        SEN.append(dict(analysis='S2 single sender, <3-producer merged into low', group=f"{o['sender']} hi:s_z", **o['m']['hi:s_z'], q=o['int_q'],
                        sections_negative=f"{o['neg']}/{o['tested']}"))
    wcsv(r3(pd.DataFrame(SEN)), f'{LR}/sensitivity_no_producer_exp.csv', index=False)
