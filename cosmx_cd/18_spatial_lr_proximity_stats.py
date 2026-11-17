#!/usr/bin/env python
"""CD CosMx 6K resections: proximity of ligand+ sender cells to receptor+ cells for the Treg activation ligand-receptor pairs.

Panels: Fig 6C (O/E, permutation P and BH q per pair), Supp 11B (sections with O/E > 1, observed pairs); drawn by 19.
Adapted from revision scripts SuppFig10B_spatial_lr_pairspecific_stats.py (statistics) and spatial_lr_check.py (gene and
negative-probe cache only).
Senders: for each ligand, the TAURUS sender group with the largest contribution to the Remission - Non-remission difference in
  change (treatment_validation/05, treg_activation_ligand_sources.csv part '3 decomposition'; Stromal = Endothelial + Fibroblast +
  Pericyte); any other group with >= 25% of the summed |contributions| is a second sender (tested alone and combined).
  CosMx mapping: Mac = SELENOP+ Mac + CD9+ M2 Mac; Mo-Mac; Stromal = cell_category Connective + Endothelial; strict Tregs never send.
Test: per section, sum over receptor+ cells of ligand+ senders within 30 um; null = ligand+ status permuted among the section's cells
  of that sender group (1,000 permutations, seed 0); O/E = observed / permutation mean; two-sided permutation P. Pooled over the 11
  sections: sum observed / sum expected, stratified permutation P, sign test of sections with O/E > 1; BH across the 14 tested rows.
  The same draws give receptor- cells, receptor+ strict Tregs (the Fig 6C receivers) and receptor- strict Tregs; 50 um = sensitivity.
  Also: receptor enrichment in strict Tregs vs conventional CD4 T, other T and all other cells (per-section paired Wilcoxon, BH).
Inputs: resection combined/all_samples.parquet and per-run Processed/norm_log1p.h5ad (03); strict Tregs from
  treg_definition_comparison_per_cell.csv (10); treg_activation_ligand_sources.csv (treatment_validation/05).
Outputs (SuppFig11_treg_ligand_senders_spatial_lr/): cache/spatial_lr/lr_genes_lognorm_allcells.parquet,
  cache/spatial_lr/negprobe_counts_allcells.parquet, SuppFig10B_spatial_lr_pairspecific_stats.csv (read by 19).
Sections (run order, after 10, 11 and treatment_validation/05): spatial_lr_cache -> F6C_lr_proximity_stats.
Run: python 18_spatial_lr_proximity_stats.py, or RUN_SECTIONS=<id>[,<id>...] for single sections (file order).
"""
import os
import warnings
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
from anndata.io import read_elem
from multiprocessing import Pool
from scipy import sparse, stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


OUT = '/share/fsmresfiles/UC/revision_jenny/SuppFig11_treg_ligand_senders_spatial_lr'
CACHE = f'{OUT}/cache/spatial_lr'
STEM = f'{OUT}/SuppFig10B_spatial_lr_pairspecific_stats'
R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
P4 = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related'
S10 = f'{OUT}/treg_activation_ligand_sources.csv'
UM, RADS, NPERM, SEED, DRIVER = 0.12028, (30, 50), 1000, 0, 25.0
sig3 = lambda v: float(f'{v:.3g}') if pd.notna(v) and np.isfinite(v) else v
FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
          '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
          '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
          '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
          '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
          '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
# 16 genes of the cache: ligands and receptors of the Treg activation pairs (incl. TNFSF4-TNFRSF4 and CD70-CD27, not tested here)
CACHE_GENES = sorted({'CD86', 'CD28', 'CD80', 'CTLA4', 'CD274', 'PDCD1', 'IL10', 'IL10RA', 'IL10RB', 'TNFSF18', 'TNFRSF18', 'PDCD1LG2',
                      'TNFSF4', 'TNFRSF4', 'CD70', 'CD27'})


def extract(sec):
    """Log-normalised expression of CACHE_GENES for every cell of one run (CSR X read once, columns picked)."""
    with h5py.File(f'{R}/{FOLDER[sec]}/Processed/norm_log1p.h5ad', 'r') as h:
        obs_i = read_elem(h['obs'][h['obs'].attrs['_index']]).astype(str); var = pd.Index(read_elem(h['var'][h['var'].attrs['_index']]).astype(str))
        ip, ind, dat = h['X/indptr'][:], h['X/indices'][:], h['X/data'][:]
    col = np.full(len(var), -1); col[[var.get_loc(g) for g in CACHE_GENES]] = np.arange(len(CACHE_GENES))
    k = col[ind]; m = k >= 0; rows = np.repeat(np.arange(len(obs_i)), np.diff(ip))[m]
    M = np.zeros((len(obs_i), len(CACHE_GENES)), np.float32); M[rows, k[m]] = dat[m]
    return pd.DataFrame(M, index=sec + '_' + pd.Index(obs_i), columns=CACHE_GENES)


def negcount(sec):
    """Negative-probe counts per cell (sum over the 21 negative probes, AtoMx nCount_negprobes)."""
    with h5py.File(f'{R}/{FOLDER[sec]}/Processed/norm_log1p.h5ad', 'r') as h:
        obs_i = read_elem(h['obs'][h['obs'].attrs['_index']]).astype(str); nc = read_elem(h['obs/nCount_negprobes'])
    return pd.Series(np.asarray(nc, float), index=sec + '_' + pd.Index(obs_i))


# ======================================================================================================================
# Section spatial_lr_cache: 16-gene log-norm extract and negative-probe counts of all 2,959,056 resection cells (built once)
# ======================================================================================================================
if _run('spatial_lr_cache'):
    os.makedirs(CACHE, exist_ok=True)
    NCACHE = f'{CACHE}/negprobe_counts_allcells.parquet'
    if not os.path.exists(NCACHE):
        with Pool(len(FOLDER)) as pool: pd.concat(pool.map(negcount, list(FOLDER))).rename('nCount_negprobes').to_frame().to_parquet(NCACHE)
    GCACHE = f'{CACHE}/lr_genes_lognorm_allcells.parquet'
    if not os.path.exists(GCACHE):
        with Pool(len(FOLDER)) as pool: E = pd.concat(pool.map(extract, list(FOLDER)))
        E.to_parquet(GCACHE)


# ======================================================================================================================
# Section F6C_lr_proximity_stats: senders from the TAURUS decomposition, proximity permutation tests per section and pooled,
# receptor enrichment in strict Tregs. Output: SuppFig10B_spatial_lr_pairspecific_stats.csv
# ======================================================================================================================
if _run('F6C_lr_proximity_stats'):
    # ---------------- senders (dot colour of Supp 11A = contribution to the difference in change, % of pre-treatment mean) ----------------
    LIGREC = [('CD86', ['CD28', 'CTLA4']), ('CD80', ['CD28', 'CTLA4']), ('IL10', ['IL10RA']), ('CD274', ['PDCD1']), ('TNFSF18', ['TNFRSF18']),
              ('PDCD1LG2', ['PDCD1'])]                                                   # Supp 11A row order; Diff P < 0.05 plus GITRL (0.050)
    STROMAL, OTHER = ['Endothelial', 'Fibroblast', 'Pericyte'], ['Granulocyte/Neutrophil', 'Other T/NK', 'Epithelial', 'Neuron (no coarse group)']
    d = pd.read_csv(S10); dc = d[d.part == '3 decomposition'].copy()
    dc['col'] = dc.cell_group.where(~dc.cell_group.isin(STROMAL), 'Stromal').where(~dc.cell_group.isin(OTHER), 'Other')

    def contrib(rg):
        x = dc[dc.remission_group == rg]; return x.groupby(['feature', 'col']).contribution_to_change.sum().unstack().div(x.groupby('feature').all_cell_mean_pre.first(), axis=0)
    VAL = 100 * (contrib('remission') - contrib('non-remission'))
    DIFFP = d[(d.part == '2 all-cell model') & (d.contrast == 'Differential change (remission − non-remission)')].set_index('feature').p
    SHARE = 100 * VAL.abs().div(VAL.abs().sum(1), axis=0)
    SND = {}
    for lig, _ in LIGREC:
        v = VAL.loc[lig].drop('Other'); top = v.abs().idxmax(); sec = [g for g in SHARE.loc[lig].drop(['Other', top]).index if SHARE.loc[lig, g] >= DRIVER]
        SND[lig] = [top] + sec
        print(f'{lig}: Diff P {DIFFP[lig]:.3g}; senders {SND[lig]}; values ' + ', '.join(f'{g} {VAL.loc[lig, g]:.1f} ({SHARE.loc[lig, g]:.1f}%)' for g in SND[lig]), flush=True)
    # result: CD86 Macrophage + Monocyte-derived Mac; CD80, IL10, CD274, TNFSF18 Macrophage; PDCD1LG2 Stromal + Macrophage
    SHORT = {'Macrophage': 'Mac', 'Monocyte-derived Mac': 'Mo-Mac', 'Stromal': 'Stromal'}
    PAIRS = []                                                                            # (ligand, receptor, sender key, role)
    for lig, recs in LIGREC:
        for rec in recs:
            g = [SHORT[x] for x in SND[lig]]
            PAIRS.append((lig, rec, g[0], 'primary'))
            if len(g) > 1: PAIRS += [(lig, rec, g[1], 'second sender'), (lig, rec, ' + '.join(g), 'combined')]
    CT = {'Mac': ('ct', ['SELENOP+ Mac', 'CD9+ M2 Mac']), 'Mo-Mac': ('ct', ['Mo-Mac']), 'Stromal': ('cat', ['Connective', 'Endothelial'])}
    LIGS = sorted({p[0] for p in PAIRS}); RECS = sorted({p[1] for p in PAIRS})

    # ---------------- cells ----------------
    E = pd.read_parquet(f'{CACHE}/lr_genes_lognorm_allcells.parquet')
    NEGF = 1 - np.exp(-pd.read_parquet(f'{CACHE}/negprobe_counts_allcells.parquet')['nCount_negprobes'] / 21)   # % cells with >= 1 count of one of 21 probes (Poisson)
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'patient_sec', 'cell_type_general',
                                                                      'cell_category', 'total_counts', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px.values * UM, sa.CenterY_global_px.values * UM
    sa['ct'], sa['cat'] = sa.cell_type_general.astype(str), sa.cell_category.astype(str)
    E = E.loc[sa.index]; NEGF = NEGF.loc[sa.index].values
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); TREG = PC.index[PC.cohort == 'resection']   # 6,554 strict Tregs
    sa['is_treg'] = sa.index.isin(TREG)
    for k, (c, v) in CT.items(): sa['snd_' + k] = sa[c].isin(v) & ~sa.is_treg
    for p in PAIRS:
        if ' + ' in p[2]: sa['snd_' + p[2]] = np.logical_or.reduce([sa['snd_' + g].values for g in p[2].split(' + ')])
    POS = E > 0
    SECS = sorted(sa.section.unique(), key=lambda s: (s.split('_')[0], s))                # 11 sections
    PAT = lambda s: s.split('_')[0]

    def run_section(si):
        """Proximity counts and permutation nulls for one section, every pair row and radius."""
        s_ = SECS[si]; cs = sa[sa.section == s_]; PS = POS.loc[cs.index]; xy = cs[['x', 'y']].values; tr = cs.is_treg.values
        snd_any = np.logical_or.reduce([cs['snd_' + k].values for k in CT]); iu = np.flatnonzero(snd_any)
        tree = cKDTree(xy); A = {}
        for r in RADS:
            nb = tree.query_ball_point(xy[iu], r=r)
            A[r] = sparse.csr_matrix((np.ones(sum(map(len, nb)), np.float32), (np.repeat(np.arange(len(iu)), [len(i) for i in nb]), np.concatenate(nb))),
                                     shape=(len(iu), len(cs)))
        pos_in_u = np.full(len(cs), -1); pos_in_u[iu] = np.arange(len(iu))
        out, nulls = [], {}
        for pi, (lig, rec, snd, role) in enumerate(PAIRS):
            sm = cs['snd_' + snd].values; ui = pos_in_u[sm]; Ls = PS[lig].values[sm]; k, nS = int(Ls.sum()), int(sm.sum())
            rp = PS[rec].values; rn = ~rp
            rec_masks = {'all': rp, 'all_rneg': rn, 'treg': rp & tr, 'treg_rneg': rn & tr}
            base = dict(pair=f'{lig}->{rec}', ligand=lig, receptor=rec, sender=snd, role=role, row='section', section=s_, patient=PAT(s_),
                        n_cells=len(cs), n_receptor_pos_cells=int(rp.sum()), n_strict_tregs=int(tr.sum()), n_receptor_pos_tregs=int((rp & tr).sum()),
                        n_senders=nS, n_ligand_pos_senders=k)
            for r in RADS:
                Ar = A[r][ui]
                W = np.column_stack([np.asarray(Ar @ m.astype(np.float32)).ravel() - (m[sm] if key.startswith('all') else 0) for key, m in rec_masks.items()])
                o = {**base, 'radius_um': r}
                if k == 0 or rp.sum() == 0:
                    out.append({**o, 'note': 'no ligand+ senders or receptor+ cells'}); continue
                obs = Ls.astype(np.float32) @ W; rng = np.random.default_rng([SEED, pi, si, r]); perm = np.zeros((NPERM, 4))
                for b in range(NPERM): perm[b] = W[rng.permutation(nS)[:k]].sum(0)
                ex = perm.mean(0); nulls[(pi, s_, r)] = (obs - ex, perm - ex)
                for j, key in enumerate(rec_masks):
                    sfx = '' if key == 'all' else '_' + key.replace('all_', '')
                    nrec = rec_masks[key].sum()
                    o[f'obs_total{sfx}'], o[f'exp_total{sfx}'] = obs[j], ex[j]
                    o[f'obs_exp{sfx}'] = obs[j] / ex[j] if ex[j] > 0 else np.nan
                    o[f'obs_mean{sfx}'] = obs[j] / nrec if nrec else np.nan
                    ph, pl = (1 + (perm[:, j] >= obs[j]).sum()) / (NPERM + 1), (1 + (perm[:, j] <= obs[j]).sum()) / (NPERM + 1)
                    o[f'perm_p{sfx}'] = min(1.0, 2 * min(ph, pl))
                out.append(o)
        print(f'{s_}: {len(cs)} cells, {int(tr.sum())} strict Tregs', flush=True)
        return out, nulls

    with Pool(len(SECS)) as pool: res = pool.map(run_section, range(len(SECS)))
    S = pd.DataFrame([o for r_ in res for o in r_[0]]); NULL = {k: v for r_ in res for k, v in r_[1].items()}
    SFX = {'': 'receptor+ cells', '_rneg': 'receptor- cells', '_treg': 'receptor+ strict Tregs', '_treg_rneg': 'receptor- strict Tregs'}
    # ---------------- detection ----------------
    DET = {}
    for lig, rec, snd, role in PAIRS:
        sm = sa['snd_' + snd].values
        DET[(lig, rec, snd)] = dict(pct_ligand_pos_senders=100 * POS[lig].values[sm].mean(), pct_negprobe_floor_senders=100 * NEGF[sm].mean(),
                                    pct_receptor_pos_all_cells=100 * POS[rec].mean(), pct_negprobe_floor_all_cells=100 * NEGF.mean(),
                                    pct_receptor_pos_strict_tregs=100 * POS.loc[TREG, rec].mean())
    # ---------------- pooled ----------------
    POOL = []
    for pi, (lig, rec, snd, role) in enumerate(PAIRS):
        for r in RADS:
            t = S[(S.pair == f'{lig}->{rec}') & (S.sender == snd) & (S.radius_um == r)].dropna(subset=['obs_exp'])
            o = dict(pair=f'{lig}->{rec}', ligand=lig, receptor=rec, sender=snd, role=role, row='pooled', section='pooled (11 sections)', patient='4 patients',
                     radius_um=r, n_cells=t.n_cells.sum(), n_receptor_pos_cells=t.n_receptor_pos_cells.sum(), n_strict_tregs=t.n_strict_tregs.sum(),
                     n_receptor_pos_tregs=t.n_receptor_pos_tregs.sum(), n_senders=t.n_senders.sum(), n_ligand_pos_senders=t.n_ligand_pos_senders.sum(),
                     **DET[(lig, rec, snd)], sections_tested=len(t))
            for j, sfx in enumerate(SFX):
                tt = t.dropna(subset=[f'obs_exp{sfx}'])
                st = sum(NULL[(pi, s_, r)][0][j] for s_ in tt.section); nl = sum(NULL[(pi, s_, r)][1][:, j] for s_ in tt.section)
                up, nt = int((tt[f'obs_exp{sfx}'] > 1).sum()), int((tt[f'obs_exp{sfx}'] != 1).sum())
                o.update({f'obs_exp{sfx}': tt[f'obs_total{sfx}'].sum() / tt[f'exp_total{sfx}'].sum(),
                          f'perm_p{sfx}': min(1.0, 2 * min((1 + (nl >= st).sum()) / (NPERM + 1), (1 + (nl <= st).sum()) / (NPERM + 1))),
                          f'sections_obs_exp_gt1{sfx}': up, f'sections_tested{sfx}': len(tt),
                          f'sign_test_p{sfx}': stats.binomtest(up, nt, 0.5).pvalue if nt else np.nan})
            o['obs_mean'] = t.obs_total.sum() / t.n_receptor_pos_cells.sum(); o['exp_mean'] = t.exp_total.sum() / t.n_receptor_pos_cells.sum()
            POOL.append(o)
    P = pd.DataFrame(POOL)
    for r in RADS:
        m = P.radius_um == r
        for sfx in SFX: P.loc[m, f'perm_q_bh{sfx}'] = multipletests(P.loc[m, f'perm_p{sfx}'], method='fdr_bh')[1]
    # ---------------- receptor enrichment in strict Tregs ----------------
    GRP = {'Strict Treg': sa.is_treg.values,
           'Conventional CD4 T': sa.ct.isin(['CD4 Effector', 'CD4 Trm', 'Tfh', 'Exhausted CD4 Tfh']).values & ~sa.is_treg.values,
           'Other T': sa.ct.isin(['CD8 Trm', 'Gd T', 'Resting T']).values,
           'All other cells': (sa.ct != 'Treg').values}                                # annotated Tregs outside the strict set are left out
    ENR = []
    for rec in RECS:
        x, xp = E[rec].values, POS[rec].values
        per = {g: pd.DataFrame({'pct': pd.Series(100 * xp[m]).groupby(sa.section.values[m]).mean(), 'mean': pd.Series(x[m]).groupby(sa.section.values[m]).mean()})
               for g, m in GRP.items()}
        T_ = GRP['Strict Treg']
        for g in list(GRP)[1:]:
            m = GRP[g]; a_, b_ = per['Strict Treg'].reindex(SECS), per[g].reindex(SECS)
            o = dict(receptor=rec, comparison_group=g, n_tregs=int(T_.sum()), n_comparison_cells=int(m.sum()),
                     pct_pos_treg=100 * xp[T_].mean(), pct_pos_comparison=100 * xp[m].mean(), mean_lognorm_treg=x[T_].mean(), mean_lognorm_comparison=x[m].mean(),
                     mean_total_counts_treg=sa.total_counts.values[T_].mean(), mean_total_counts_comparison=sa.total_counts.values[m].mean())
            o['pct_pos_ratio'] = o['pct_pos_treg'] / o['pct_pos_comparison']; o['pct_pos_difference'] = o['pct_pos_treg'] - o['pct_pos_comparison']
            o['mean_lognorm_ratio'] = o['mean_lognorm_treg'] / o['mean_lognorm_comparison']
            for met in ['pct', 'mean']:
                dd = (a_[met] - b_[met]).dropna()
                o[f'sections_treg_higher_{met}'] = int((dd > 0).sum()); o[f'sections_{met}'] = len(dd)
                o[f'wilcoxon_p_{met}'] = stats.wilcoxon(dd).pvalue
            ENR.append(o)
    ENR = pd.DataFrame(ENR)
    for met in ['pct', 'mean']: ENR[f'wilcoxon_q_bh_{met}'] = multipletests(ENR[f'wilcoxon_p_{met}'], method='fdr_bh')[1]

    # ---------------- CSV ----------------
    def tidy(df):
        for c in df.columns:
            if df[c].dtype.kind == 'f' and not (c.startswith('n_') or c.startswith('sections_') or c == 'radius_um'): df[c] = df[c].map(sig3)
            elif df[c].dtype.kind == 'f': df[c] = df[c].round().astype('Int64')
        return df
    ALL = pd.concat([P, S], ignore_index=True)
    ALL['pi'] = [next(i for i, p in enumerate(PAIRS) if f'{p[0]}->{p[1]}' == a and p[2] == b) for a, b in zip(ALL.pair, ALL.sender)]
    ALL['si'] = ALL.section.map({s: i for i, s in enumerate(SECS)}).fillna(-1)
    ALL = ALL.sort_values(['pi', 'radius_um', 'si']).drop(columns=['pi', 'si'])
    lead = ['pair', 'ligand', 'receptor', 'sender', 'role', 'row', 'section', 'patient', 'radius_um', 'obs_exp', 'perm_p', 'perm_q_bh', 'sections_obs_exp_gt1',
            'sections_tested', 'sign_test_p', 'obs_exp_rneg', 'obs_exp_treg', 'perm_p_treg', 'perm_q_bh_treg', 'sections_obs_exp_gt1_treg', 'sign_test_p_treg',
            'obs_exp_treg_rneg']
    ALL = ALL[lead + [c for c in ALL.columns if c not in lead and c != 'note'] + (['note'] if 'note' in ALL else [])]
    ALL.insert(0, 'table', 'proximity'); ENR.insert(0, 'table', 'receptor enrichment (strict Treg vs comparison group)')
    with open(f'{STEM}.csv', 'w') as fh: tidy(pd.concat([ALL, ENR], ignore_index=True)).to_csv(fh, index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    print(P[['pair', 'sender', 'radius_um', 'pct_ligand_pos_senders', 'pct_negprobe_floor_senders', 'pct_receptor_pos_all_cells', 'obs_exp', 'perm_p', 'perm_q_bh',
             'sections_obs_exp_gt1', 'sign_test_p', 'obs_exp_rneg', 'obs_exp_treg', 'perm_q_bh_treg', 'sections_obs_exp_gt1_treg', 'obs_exp_treg_rneg']].to_string())
    print(ENR[['receptor', 'comparison_group', 'pct_pos_treg', 'pct_pos_comparison', 'pct_pos_ratio', 'mean_lognorm_ratio', 'sections_treg_higher_pct',
               'wilcoxon_q_bh_pct', 'wilcoxon_q_bh_mean', 'mean_total_counts_treg', 'mean_total_counts_comparison']].to_string(), flush=True)
    print('wrote', f'{STEM}.csv', flush=True)
