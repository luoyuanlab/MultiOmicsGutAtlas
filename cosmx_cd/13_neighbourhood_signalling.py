#!/usr/bin/env python
"""CD CosMx 6K resections: ligand-receptor (LR) signalling per cellular neighbourhood (CN).

Panels: Fig 5C (F5C_score_signalling), Supp 9A (S9A_cn_characterisation), Supp 9B (S9B_cn_signalling); Supp Table 10 (LR module membership).
Also writes the Treg receiver pair table and per-Treg module scores read by 17_neighbourhood_apc_mhcii.py.
Inputs: combined/spot_lr_df.parquet (07), all_samples.parquet, raw_counts_strict_tregs.csv (05), per-section norm_log1p.h5ad (03), CellChatDB snapshot,
  NicheNet network, and the Treg score (10), CN tables (11) and CN composition / CN forest tables (12).
Run: python 13_neighbourhood_signalling.py, or RUN_SECTIONS=<id>[,<id>] python 13_neighbourhood_signalling.py (sections run in file order).
"""
import os, gc, warnings
import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.ticker, seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.transforms import blended_transform_factory, ScaledTranslation
from seaborn.utils import relative_luminance
from scipy import sparse, stats as sst
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
warnings.filterwarnings('ignore')

_RUN_SECTIONS = [s.strip() for s in os.environ.get('RUN_SECTIONS', '').split(',') if s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# ---------------- paths
D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; SPOT = f'{R}/combined/spot_lr_df.parquet'
CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'; CCDB = f'{D}/objects/cellchatdb_human_interactions_snapshot.csv'
C45 = f'{P4}/A_fig4_panels/fig4_cn_composition_heatmap.csv'; C26 = f'{P4}/A_fig4_panels/fig4_cn_treg_state_forest.csv'   # from 12_neighbourhood_treg_score.py
OUTA = f'{P4}/A_fig4_panels'; LRC = f'{P4}/B_resection_analyses/lr_modules_cellchat'; INF = f'{P4}/B_resection_analyses/inflammatory_axes'
FIG5 = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches'; SUPP9 = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc'
# files passed between sections
MEMB82 = f'{LRC}/lr_modules_allcells_membership.csv'; CELL82 = f'{LRC}/lr_modules_allcells_per_cell.parquet'
MEMB84 = f'{LRC}/lr_modules_v2_membership.csv'; PAIR84 = f'{LRC}/lr_pairs_v2_per_section_cn.csv'; DET84 = f'{LRC}/lr_pairs_v2_detection.csv'
C82D = f'{OUTA}/fig4_lr_modules_allcells_costim_dotplot.csv'; C86 = f'{OUTA}/supp_inflammatory_axes_by_neighbourhood.csv'
C89 = f'{SUPP9}/supp_cn_niche_characterisation.csv'

# ---------------- shared helpers
CNS = list(range(1, 11)); MIN_CELLS_CN = 100
stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
zs = lambda x: (x - x.mean()) / x.std()
gs = lambda x: frozenset(y.strip() for y in str(x).replace(';', ',').split(',') if y.strip() and y.strip() != 'nan')
GROUPS = [('T-Cell Signals', ['Antigen Presentation', 'Co-Stimulatory']), ('Inflammatory Cytokines', ['IL-6', 'IL-1', 'OSM', 'IL-11', 'TNF']),
          ('Stroma', ['Matrix–Integrin'])]
COLS = [c for _, cs in GROUPS for c in cs]; GOF = {c: g for g, cs in GROUPS for c in cs}
SRC = {c: ('Ligand–Receptor Modules' if c in ('Antigen Presentation', 'Co-Stimulatory', 'Matrix–Integrin') else 'Inflammatory Signalling') for c in COLS}
PAL = dict(zip(CNS, sns.color_palette('tab10', 10).as_hex()))
SQX, SQW = 0.05, 0.13   # CN colour key: gap to the y axis and width (inches)
umf = lambda: matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−'))   # tick labels with a true minus sign


def style(fs):
    plt.rcParams.update({'font.size': fs, 'font.family': 'DejaVu Sans', 'font.weight': 'normal', 'axes.titleweight': 'normal',
                         'axes.labelweight': 'normal', 'figure.titleweight': 'normal'})


def cn_key(ax, y, h):
    """CN colour key between the y axis and its row labels; y = row centres (data), h = key height (inches)."""
    tr = blended_transform_factory(ax.transAxes, ax.transData) + ScaledTranslation(-(SQX + SQW / 2), 0, ax.figure.dpi_scale_trans)
    ax.scatter(np.zeros(len(y)), y, s=(72 * h) ** 2, marker=[(-SQW, -h), (SQW, -h), (SQW, h), (-SQW, h)], c=[PAL[c] for c in CNS], lw=0,
               transform=tr, clip_on=False)
    ax.tick_params(axis='y', length=0, pad=(SQX + SQW + 0.07) * 72)


def cellchat():
    """CellChatDB snapshot with ligand / receptor gene sets, and the four CellChatDB modules (pathway_name lists)."""
    cc = pd.read_csv(CCDB); cc['L'] = cc['ligand.symbol'].map(gs); cc['R'] = cc['receptor.symbol'].map(gs)
    MOD = {'Antigen presentation': ['MHC-II'], 'Co-stimulation / suppressive Treg': ['CD80', 'CD86', 'PD-L1', 'PDL2', 'ICOS', 'OX40', 'GITRL', 'CD70', 'CD40', 'IL2', 'IL10'],
           'Inflammatory cytokine': ['IL1', 'IL6', 'TNF', 'IFN-II', 'IL12', 'IL17', 'OSM'],
           'Matrix-integrin': sorted(cc.loc[cc.annotation == 'ECM-Receptor', 'pathway_name'].unique())}
    return cc, MOD


def cn_labels():
    """Full CN names and the compact labels used in the figures."""
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [(' Fib / PLA2G2A⁺ Fib / Endothelial', '/PLA2G2A⁺ Fib'), (' / Crypt Top Fib / Cycling', ' / CryptTop / Cycl'),
            (' / ADAMDEC1 Fib', ' / ADAMDEC1'), (' / Neuron / Pericyte', ' / Neuron / Peri'), (' / Mac S_M_P / CD9 Mac', ' / Mac / CD9 Mac'),
            ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'), ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'),
            ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    LAB = {}
    for c in CNS:
        x = NM[c]
        for a, b in TRIM: x = x.replace(a, b)
        LAB[c] = x
    return NM, LAB


def fit_in_cn(dd):
    """y ~ in_cn + (1 | section); optimiser fallback if the Hessian fails."""
    for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
        try:
            f = smf.mixedlm('y ~ in_cn', dd, groups=dd.section).fit(reml=True, method=meth)
            if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], meth
        except Exception: pass
    return np.nan, np.nan, np.nan, 'failed'


def cn_stats(d, c):
    """One CN vs the other CNs on per-section within-section z (column y): mean + t 95% CI, Wilcoxon, mixed model."""
    dd = d.assign(in_cn=(d.CN == c).astype(float)); zc = dd[dd.CN == c].y.values
    b, se, p, meth = fit_in_cn(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
    t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
    return dict(n_sections=len(zc), mean_within_section_z=zc.mean(), mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)),
                mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)), sections_above_0=int((zc > 0).sum()),
                wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan, lmm_estimate=b, lmm_ci_low=b - 1.96 * se, lmm_ci_high=b + 1.96 * se,
                lmm_p=p, optimiser=meth)


def add_bh(ST, by):
    """BH over the 10 CNs within each module / axis, for the mixed model and the Wilcoxon test."""
    ST['lmm_q'] = np.nan; ST['wilcoxon_q'] = np.nan
    for _, g in ST.groupby(by):
        for k_ in ['lmm', 'wilcoxon']:
            ok = g[f'{k_}_p'].notna(); ST.loc[g.index[ok], f'{k_}_q'] = multipletests(g[f'{k_}_p'][ok], method='fdr_bh')[1]
    return ST


def cn_section_table():
    """Factorised CN x section groups of the LIANA cells with a CN (shared by the per-pair computations)."""
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
    rows_ = pd.Index(pq.read_table(SPOT, columns=['__index_level_0__']).column(0).to_pylist()); keep = rows_.isin(cn.index)
    sec = cn.section.reindex(rows_[keep]).values; cnn = cn.CN_k10.reindex(rows_[keep]).values.astype(int)
    codes, uniq = pd.factorize(pd.MultiIndex.from_arrays([sec, cnn], names=['section', 'CN']))
    M = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    return M, np.asarray(M.sum(1)).ravel(), np.where(keep)[0], uniq


# ======================================================================================================================
# Section treg_receiver_lr_modules: per-Treg ligand supply from the canonical sender cells (input of 17_neighbourhood_apc_mhcii.py)
# ======================================================================================================================
if _run('treg_receiver_lr_modules'):
    import h5py, pyreadr
    from anndata.io import read_elem
    from scipy.spatial import cKDTree
    UM = 0.12028; RAD = 50; MIN_SEND = 1; MIN_SEND_SENS = 3
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    cc, MOD = cellchat(); MODS = list(MOD)
    SETS = {'MHCII_APC': ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac'], 'DC_B': ['DC', 'B cell'], 'CD4_T': ['CD4 Effector', 'CD4 Trm', 'Tfh'],
            'IL10_myeloid': ['DC', 'CD9+ M2 Mac', 'Mo-Mac'], 'Inflam_myeloid': ['Mo-Mac', 'CD9+ M2 Mac', 'DC', 'Neutrophil'],
            'Inflam_T_NK': ['CD4 Effector', 'CD8 Trm', 'NK', 'Gd T'],
            'Fibroblast': ['OGN+RSPO3+ Fib', 'ADAMDEC1+ Fib', 'Activated Fib', 'Crypt Top Fib', 'PLA2G2A+ ECM Fib', 'Myofibroblast']}

    def sender_set(module, pathway):   # main canonical senders per pathway (Tregs are never senders)
        if module == 'Antigen presentation': return 'MHCII_APC'
        if module == 'Co-stimulation / suppressive Treg': return 'CD4_T' if pathway in ('IL2', 'CD40') else 'IL10_myeloid' if pathway == 'IL10' else 'DC_B'
        if module == 'Inflammatory cytokine': return 'Inflam_T_NK' if pathway in ('IFN-II', 'IL17') else 'Inflam_myeloid'
        return 'Fibroblast'
    with h5py.File(f'{R}/{FOLDER["43563_A8"]}/Processed/norm_log1p.h5ad', 'r') as h: PANEL = set(read_elem(h['var']).index.astype(str))
    # NicheNet pairs; a pair joins a module if a CellChatDB interaction of the module has the ligand among its ligand genes and the receptor among its receptor genes
    NN = list(pyreadr.read_r('/share/fsmresfiles/ydn4687/nichenet/lr_network_human_21122021.rds').values())[0][['from', 'to']].drop_duplicates()
    PR = []
    for m in MODS:
        sub = cc.pathway_name.isin(MOD[m])
        for l, r in zip(NN['from'], NN['to']):
            hit = sorted({pw for pw, a, b in zip(cc.pathway_name[sub], cc.L[sub], cc.R[sub]) if l in a and r in b})
            if not hit: continue
            miss = [g for g in (l, r) if g not in PANEL]
            PR.append(dict(module=m, pathway='/'.join(hit), interaction=f'{l}_{r}', ligand=l, receptor=r, sender_set=sender_set(m, hit[0]), kept=not miss,
                           missing_genes=','.join(miss), pair_source='NicheNet v2 lr_network_human_21122021'))
    PR = pd.DataFrame(PR); PR['sender_cell_types'] = PR.sender_set.map(lambda s: ', '.join(SETS[s]))
    PR.to_csv(f'{LRC}/treg_receiver_pairs_kept.csv', index=False)
    KP = PR[PR.kept].reset_index(drop=True)
    LUN = sorted(set(KP.ligand)); RUN = sorted(set(KP.receptor)); GENES = sorted({g for u in LUN + RUN for g in u.split('_')})
    LSET = KP[['ligand', 'sender_set']].drop_duplicates()
    NEED = {s: sorted(LSET[LSET.sender_set == s].ligand) for s in LSET.sender_set.unique()}
    DENS = {m: sorted({c for s in KP[KP.module == m].sender_set.unique() for c in SETS[s]}) for m in MODS}
    # per-Treg ligand (mean within-section z over the pair's senders within 50 um), own receptor z and sender counts
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); pcr = PC[PC.cohort == 'resection']
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
    T = sa.loc[pcr.index, ['section']].copy()
    LCOL = {(s, u): f'L|{s}|{u}' for s in NEED for u in NEED[s]}; RCOL = {u: f'R|{u}' for u in RUN}
    blocks = []
    for s_, f in FOLDER.items():
        sub = sa[sa.section == s_]
        with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
            ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
        nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
        keep = np.where(nm.isin(sub.index))[0]
        E = pd.DataFrame(X[keep][:, [v.get_loc(g) for g in GENES]].toarray().astype(np.float32), index=nm[keep], columns=GENES); del X; gc.collect()
        cs = sub.loc[E.index]
        U = {u: E[u.split('_')].min(1).values for u in sorted(set(LUN) | set(RUN))}       # complex = minimum of subunits
        Z = {}
        for u, x in U.items():
            sd = x.std(); Z[u] = (x - x.mean()) / sd if sd > 0 else np.full_like(x, np.nan)
        tg = T[T.section == s_]
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
        ctv = cs.ct.values; out = pd.DataFrame(index=tg.index)
        out['nb_n'] = [len(i) for i in nb]
        rows_ = np.repeat(np.arange(len(nb)), [len(i) for i in nb]); cols_ = np.concatenate([np.asarray(i, dtype=int) for i in nb])
        for sname, types in SETS.items():
            sel = np.isin(ctv, types)[cols_]; Nm = sparse.csr_matrix((np.ones(sel.sum(), np.float32), (rows_[sel], cols_[sel])), shape=(len(nb), len(cs)))
            n = np.asarray(Nm.sum(1)).ravel(); out[f'n|{sname}'] = n
            if sname in NEED:
                ZL = np.nan_to_num(np.column_stack([Z[u] for u in NEED[sname]]), nan=0.0)
                Lm = np.asarray(Nm @ ZL) / np.where(n > 0, n, 1)[:, None]; Lm[n < 1] = np.nan
                for j, u in enumerate(NEED[sname]): out[LCOL[(sname, u)]] = Lm[:, j] if np.isfinite(Z[u]).any() else np.nan
        for m_, types in DENS.items(): out[f'dens|{m_}'] = np.asarray([np.isin(ctv[np.asarray(i, dtype=int)], types).sum() for i in nb])
        pos = E.index.get_indexer(tg.index)
        for u in RUN: out[RCOL[u]] = Z[u][pos]
        blocks.append(out); print(s_, 'Tregs', len(tg), 'cells', len(cs), flush=True)
        del E, U, Z; gc.collect()
    W = pd.concat(blocks).loc[T.index]; del sa; gc.collect()
    W.insert(0, 'section', T.section)
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
    W.insert(1, 'CN', cn.loc[W.index, 'CN_k10'].astype(int).values)
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0)
    W.insert(2, 'log_total_counts', np.log1p(cnt.loc[W.index].sum(1).values)); del cnt; gc.collect()
    # module scores: each pair from its own sender group, counted if >= MIN_SEND senders within 50 um (>= MIN_SEND_SENS: sensitivity)
    def Lmat(units, thr):
        return np.column_stack([np.where(W[f'n|{s}'].values >= thr, W[f'L|{s}|{l}'].values, np.nan) for l, s in units])
    for m in MODS:
        kp = KP[KP.module == m]; lu = kp[['ligand', 'sender_set']].drop_duplicates()
        Rp = np.column_stack([W[f'R|{r}'].values for r in kp.receptor])
        for thr, sfx in [(MIN_SEND, ''), (MIN_SEND_SENS, ' min3')]:
            Lp = Lmat(zip(kp.ligand, kp.sender_set), thr); cnt_ = np.isfinite(Lp)
            W[f'{m}|n_pairs_counting{sfx}'] = cnt_.sum(1)
            W[f'{m}|LR{sfx}'] = np.nanmean(np.where(cnt_, np.clip(Lp, 0, None) * np.clip(Rp, 0, None), np.nan), 1)
            W[f'{m}|LxR raw{sfx}'] = np.nanmean(np.where(cnt_, Lp * Rp, np.nan), 1)
            W[f'{m}|L{sfx}'] = np.nanmean(Lmat(zip(lu.ligand, lu.sender_set), thr), 1)
        W[f'{m}|Sender density'] = np.log1p(W[f'dens|{m}'].values)
    W.to_csv(f'{LRC}/treg_receiver_module_scores_per_treg.csv.gz', float_format='%.4g')


# ======================================================================================================================
# Section cn_lr_modules_allcells: CellChatDB module membership (Supp Table 10) and per-cell raw module scores over all cells
# ======================================================================================================================
if _run('cn_lr_modules_allcells'):
    cc, MOD = cellchat()
    cols = [c for c in pq.read_schema(SPOT).names if '^' in c]

    REFS = {m: list(zip(cc[cc.pathway_name.isin(pws)].L, cc[cc.pathway_name.isin(pws)].R)) for m, pws in MOD.items()}

    def module_of(c):   # same ligand set and an overlapping receptor set
        l, r = c.split('^'); L, Rr = frozenset(l.split('_')), frozenset(r.split('_'))
        return [m for m, refs in REFS.items() if any((L == a) and (Rr == b or Rr <= b or b <= Rr) for a, b in refs)]
    MEM = {c: m for c, m in ((c, module_of(c)) for c in cols) if m}
    pd.DataFrame([(c, m) for c, ms in MEM.items() for m in ms], columns=['lr_pair', 'module']).to_csv(MEMB82, index=False)
    idx_all = pq.read_table(SPOT, columns=['__index_level_0__']).column(0).to_pylist()
    acc = {m: np.zeros(len(idx_all)) for m in MOD}; npair = {m: 0 for m in MOD}
    for i in range(0, len(cols), 250):
        ch = [c for c in cols[i:i + 250] if c in MEM]
        if not ch: continue
        A = pq.read_table(SPOT, columns=ch).to_pandas().astype(np.float32).fillna(0.0)
        for m in MOD:
            ps = [c for c in ch if m in MEM[c]]
            if ps: acc[m] += A[ps].values.sum(1); npair[m] += len(ps)
        del A; gc.collect()
    C = pd.DataFrame({f'raw|{m}': acc[m] / npair[m] for m in MOD}, index=pd.Index(idx_all, name='cell'))   # mean raw local score over the module's pairs
    C.astype(np.float32).to_parquet(CELL82)


# ======================================================================================================================
# Section lr_costim_inhibitory_modules: pre-specified co-stimulatory / immunoregulatory / gamma-c modules (Supp Table 10) and their per-pair CN x section means
# ======================================================================================================================
if _run('lr_costim_inhibitory_modules'):
    CS, IR, GC = 'Co-Stimulatory (Activating)', 'Immunoregulatory (Inhibitory / Suppressive)', 'Treg-Supportive γc Cytokines'
    AP, IC, MI, BROAD = 'Antigen presentation', 'Inflammatory cytokine', 'Matrix-integrin', 'Co-stimulation / suppressive Treg'
    cc, _ = cellchat(); F = frozenset

    def ccref(pw, lig=None, rec_any=None, rec_eq=None):
        s = cc[cc.pathway_name == pw]; out = []
        for a, b in zip(s.L, s.R):
            if lig is not None and not (a <= F(lig)): continue
            if rec_any is not None and not (b & F(rec_any)) or rec_any is not None and not (b <= F(rec_any)): continue
            if rec_eq is not None and b != F(rec_eq): continue
            out.append((a, b))
        return out
    gn = lambda ligs, recs: [(F([l]), F([r])) for l in ligs for r in recs]
    SPEC = {CS: [('CD80 -> CD28', 'CellChatDB CD80 (receptor CD28 only)', ccref('CD80', rec_eq=['CD28'])),
                 ('CD86 -> CD28', 'CellChatDB CD86 (receptor CD28 only)', ccref('CD86', rec_eq=['CD28'])),
                 ('ICOS', 'CellChatDB ICOS (receptors ICOS, CD28)', ccref('ICOS', rec_any=['ICOS', 'CD28'])),
                 ('OX40', 'CellChatDB OX40', ccref('OX40')), ('GITRL', 'CellChatDB GITRL', ccref('GITRL')),
                 ('4-1BB', 'CellChatDB CD137', ccref('CD137', lig=['TNFSF9'], rec_eq=['TNFRSF9'])),
                 ('CD70', 'CellChatDB CD70', ccref('CD70')), ('CD40', 'CellChatDB CD40 (CD40LG -> CD40 only)', ccref('CD40', rec_eq=['CD40']))],
            IR: [('PD-L1', 'CellChatDB PD-L1', ccref('PD-L1')), ('PDL2', 'CellChatDB PDL2', ccref('PDL2')),
                 ('CD80 / CD86 -> CTLA4', 'CellChatDB CD80, CD86 (receptor CTLA4 only)', ccref('CD80', rec_eq=['CTLA4']) + ccref('CD86', rec_eq=['CTLA4'])),
                 ('NECTIN / PVR -> TIGIT, CD96', 'gene names (CellChatDB lists CD96 as ligand)', gn(['NECTIN2', 'PVR'], ['TIGIT', 'CD96'])),
                 ('GALECTIN-9 -> TIM-3', 'CellChatDB GALECTIN (receptor HAVCR2 only)', ccref('GALECTIN', rec_eq=['HAVCR2'])),
                 ('HVEM -> BTLA, CD160', 'gene names (CellChatDB lists BTLA / CD160 as ligand)', gn(['TNFRSF14'], ['BTLA', 'CD160'])),
                 ('CD200 -> CD200R1', 'CellChatDB CD200', ccref('CD200')),
                 ('HLA-E -> NKG2A', 'gene names', gn(['HLA-E'], ['KLRC1'])),
                 ('IL10', 'CellChatDB IL10 (ligand IL10 only)', ccref('IL10', lig=['IL10'])),
                 ('TGFb', 'CellChatDB TGFb (receptors TGFBR1, TGFBR2)', ccref('TGFb', rec_eq=['TGFBR1', 'TGFBR2']))],
            GC: [('IL2', 'CellChatDB IL2 (ligand IL2)', ccref('IL2', lig=['IL2'])), ('IL7', 'CellChatDB IL2 (ligand IL7)', ccref('IL2', lig=['IL7'])),
                 ('IL15', 'CellChatDB IL2 (ligand IL15)', ccref('IL2', lig=['IL15']))]}
    cols = [c for c in pq.read_schema(SPOT).names if '^' in c]

    def match(c, refs):
        l, r = c.split('^'); L, Rr = F(l.split('_')), F(r.split('_'))
        return any((L == a) and (Rr == b or Rr <= b or b <= Rr) for a, b in refs)
    NEW = pd.DataFrame([dict(module=m, subfamily=fam, source=src_, lr_pair=c) for m, fams in SPEC.items() for fam, src_, refs in fams for c in cols if match(c, refs)])
    m82 = pd.read_csv(MEMB82)
    MEMB = pd.concat([m82[m82.module.isin([AP, IC, MI])].assign(subfamily='', source='cn_lr_modules_allcells'), NEW], ignore_index=True)[['module', 'subfamily', 'source', 'lr_pair']]
    MEMB['in_costim_suppressive_treg_module'] = MEMB.lr_pair.isin(m82[m82.module == BROAD].lr_pair)
    MEMB.to_csv(MEMB84, index=False)
    # per pair x section x CN means over all cells, and the fraction of cells with a non-zero score
    NP = NEW.lr_pair.tolist(); M, ncell, kidx, uniq = cn_section_table(); out = []; det = []
    for i in range(0, len(NP), 64):
        ch = NP[i:i + 64]; A = pq.read_table(SPOT, columns=ch).to_pandas().astype(np.float64).fillna(0.0).values[kidx]
        out.append(pd.DataFrame(np.asarray(M @ A) / ncell[:, None], columns=ch))
        det.append(pd.DataFrame({'lr_pair': ch, 'frac_cells_nonzero': (A > 0).mean(0), 'mean_raw_all_cells': A.mean(0)})); del A; gc.collect()
    W = pd.concat(out, axis=1); W.index = pd.MultiIndex.from_tuples(list(uniq), names=['section', 'CN']); W['n_cells'] = ncell
    W.reset_index().melt(id_vars=['section', 'CN', 'n_cells'], var_name='lr_pair', value_name='raw').to_csv(PAIR84, index=False, float_format='%.8g')
    pd.concat(det).to_csv(DET84, index=False, float_format='%.4g')


# ======================================================================================================================
# Section cn_lr_modules_costim: module statistics per CN (Antigen Presentation, Inflammatory Cytokine, 7-pair Co-Stimulatory, Matrix-Integrin; reference
# modules in the CSV). CSV = module source of S9A_cn_characterisation, F5C_score_signalling and S9B_cn_signalling.
# ======================================================================================================================
if _run('cn_lr_modules_costim'):
    AP, BROAD, IC, MI = 'Antigen presentation', 'Co-stimulation / suppressive Treg', 'Inflammatory cytokine', 'Matrix-integrin'
    CS9, IR, IMM, CS7 = 'Co-Stimulatory (Activating)', 'Immunoregulatory (Inhibitory / Suppressive)', 'Immunomodulatory', 'Co-Stimulatory (APC -> T cell)'
    MODS = [AP, IC, CS7, MI]; REF = [IR, CS9, IMM, BROAD]; MOD4 = [AP, BROAD, IC, MI]
    vmem = pd.read_csv(MEMB84)
    CS_P = vmem[vmem.module == CS9].lr_pair.tolist(); IR_P = vmem[vmem.module == IR].lr_pair.tolist(); IMM_P = CS_P + IR_P
    CS7_P = ['CD80^CD28', 'CD86^CD28', 'ICOSLG^ICOS', 'TNFSF4^TNFRSF4', 'TNFSF18^TNFRSF18', 'TNFSF9^TNFRSF9', 'CD70^CD27']   # ligand on the APC, receptor on the T cell
    CS_DROP = {'CD40LG^CD40': 'T cell -> APC direction (CD40L licenses APCs; does not co-stimulate T cells)', 'ICOSLG^CD28': 'non-canonical interaction'}
    PAIRS = {IMM: IMM_P, CS9: CS_P, IR: IR_P, CS7: CS7_P}
    NPAIR = {AP: 9, BROAD: 50, IC: 34, MI: 227, IMM: 28, CS9: 9, IR: 19, CS7: 7}
    # CellChatDB modules: per-cell scores -> CN x section means (float32 path)
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
    C = pd.read_parquet(CELL82, columns=[f'raw|{m}' for m in MOD4])
    X = C.join(cn[['CN_k10', 'section']], how='inner').rename(columns={'CN_k10': 'CN'}); G = X.groupby(['section', 'CN'])
    V82 = (G[[f'raw|{m}' for m in MOD4]].mean().set_axis(MOD4, axis=1).reset_index().melt(id_vars=['section', 'CN'], var_name='module', value_name='raw')
           .merge(G.size().rename('n_cells').reset_index(), on=['section', 'CN']))
    del C, X
    # pre-specified modules: mean over the per-pair CN x section means
    L = pd.read_csv(PAIR84)
    VN = pd.concat([L[L.lr_pair.isin(ps)].groupby(['section', 'CN']).agg(raw=('raw', 'mean'), n_cells=('n_cells', 'first')).reset_index().assign(module=m)
                    for m, ps in PAIRS.items()], ignore_index=True)
    VN['raw'] = VN.raw.astype(np.float32)
    V = pd.concat([V82[['section', 'CN', 'module', 'raw', 'n_cells']], VN[['section', 'CN', 'module', 'raw', 'n_cells']]], ignore_index=True)
    V = V[V.n_cells >= MIN_CELLS_CN].copy()
    V['z_raw'] = V.groupby(['section', 'module']).raw.transform(zs); V.loc[~np.isfinite(V.z_raw), 'z_raw'] = np.nan
    NM, LAB = cn_labels()
    ST = []
    for m in MODS + REF:
        d = V[V.module == m].dropna(subset=['z_raw']); d = d.assign(y=d.z_raw)
        for c in CNS: ST.append(dict(module=m, CN=c, CN_label=LAB[c], CN_name=NM[c], n_pairs=NPAIR[m], **cn_stats(d, c)))
    ST = add_bh(pd.DataFrame(ST), 'module')
    SHORTM = {AP: 'Antigen Presentation', CS7: 'Co-Stimulatory', IC: 'Inflammatory Cytokine', MI: 'Matrix–Integrin', IR: 'Inhibitory / Suppressive',
              CS9: 'Co-Stimulatory (9 pairs, reference)', IMM: 'Immunomodulatory (28 pairs, reference)', BROAD: 'Co-Stimulation / Suppressive Treg (50 pairs, reference)'}
    # per-pair co-stimulatory values (9 pairs, CSV only)
    DET = pd.read_csv(DET84).set_index('lr_pair')
    LP = L[L.lr_pair.isin(CS_P)].copy(); LP = LP[LP.n_cells >= MIN_CELLS_CN].copy()
    LP['z'] = LP.groupby(['section', 'lr_pair']).raw.transform(zs); LP.loc[~np.isfinite(LP.z), 'z'] = np.nan
    PP = []
    for p in CS_P:
        d2 = LP[LP.lr_pair == p].dropna(subset=['z'])
        for c in CNS:
            zc = d2[d2.CN == c].z.values
            PP.append(dict(module='Co-Stimulatory pair (not plotted)', lr_pair=p, pair_label=p.replace('^', ' → '),
                           pair_status='kept (in the 7-pair module)' if p in CS7_P else 'removed: ' + CS_DROP[p].replace('->', '→'),
                           frac_cells_nonzero=DET.loc[p, 'frac_cells_nonzero'], CN=c, CN_label=LAB[c], CN_name=NM[c], n_sections=len(zc),
                           mean_within_section_z=zc.mean() if len(zc) else np.nan, sections_above_0=int((zc > 0).sum())))
    PP = pd.DataFrame(PP)
    METRIC = 'LIANA bivariate local score (raw), mean over module pairs per cell, mean over ALL cells of the CN per section; z across CNs within section, mean over sections'
    MODEL = 'z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within module'
    PAIRS_TXT = {m: '; '.join(p.replace('^', ' → ') for p in ps) for m, ps in [(CS7, CS7_P), (CS9, CS_P), (IR, IR_P), (IMM, IMM_P)]}
    SUM = ST.assign(shown_in_figure=ST.module.isin(MODS), stars=ST.lmm_q.map(stars), lr_pairs=ST.module.map(PAIRS_TXT), module=ST.module.map(SHORTM), row_type='summary',
                    cell_text=[f'{e:+.2f}{stars(q)}' for e, q in zip(ST.mean_within_section_z, ST.lmm_q)])   # exact 2-dp heatmap text
    SEC = V.assign(module=V.module.map(SHORTM), shown_in_figure=V.module.isin(MODS), row_type='section').rename(columns={'z_raw': 'within_section_z'})
    PPO = PP.assign(shown_in_figure=False, row_type='pair_summary', n_pairs=1,
                    metric='per pair: LIANA bivariate local score (raw), mean over ALL cells of the CN per section; z across CNs within section, mean over sections; '
                           'frac_cells_nonzero = fraction of all cells with a non-zero local score', model='')
    OUTS = pd.concat([SUM.assign(metric=METRIC, model=MODEL), SEC.assign(metric=METRIC, model=MODEL), PPO], ignore_index=True)
    for c_ in ['n_cells', 'n_pairs', 'n_sections', 'sections_above_0']: OUTS[c_] = OUTS[c_].round().astype('Int64')
    OUTS.apply(sig3).to_csv(C82D, index=False)


# ======================================================================================================================
# Section cn_inflammatory_axes: pre-specified inflammatory axes (IL-6, IL-1, OSM, IL-11; TNF contrast) per CN, computed like an LR module
# ======================================================================================================================
if _run('cn_inflammatory_axes'):
    AXES = {'IL-6': ['IL6^IL6R', 'IL6^IL6ST'], 'IL-1': ['IL1A^IL1R1', 'IL1A^IL1RAP', 'IL1B^IL1R1', 'IL1B^IL1RAP'],
            'OSM': ['OSM^OSMR', 'OSM^LIFR', 'OSM^IL6ST'], 'IL-11': ['IL11^IL11RA', 'IL11^IL6ST'], 'TNF': ['TNF^TNFRSF1A', 'TNF^TNFRSF1B']}
    AXN = list(AXES); PAIRS = [p for ps in AXES.values() for p in ps]; C1 = f'{INF}/axes_per_section_cn.csv'
    os.makedirs(INF, exist_ok=True)
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
    A = pq.read_table(SPOT, columns=PAIRS + ['__index_level_0__']).to_pandas()
    A = A[PAIRS].astype(np.float64).fillna(0.0)
    Cax = pd.DataFrame({a: A[ps].values.mean(1) for a, ps in AXES.items()}, index=A.index); del A; gc.collect()
    X = Cax.join(cn[['CN_k10', 'section']], how='inner').rename(columns={'CN_k10': 'CN'}); G = X.groupby(['section', 'CN'])
    VD = G[AXN].mean().reset_index().melt(id_vars=['section', 'CN'], var_name='axis', value_name='raw').merge(G.size().rename('n_cells').reset_index(), on=['section', 'CN'])
    VD.to_csv(C1, index=False, float_format='%.8g'); del X, Cax, cn; gc.collect()
    VD = pd.read_csv(C1)
    V = VD[VD.n_cells >= MIN_CELLS_CN].copy()
    V['z'] = V.groupby(['section', 'axis']).raw.transform(zs); V.loc[~np.isfinite(V.z), 'z'] = np.nan
    NM, LAB = cn_labels()
    ST = []
    for a in AXN:
        d = V[V.axis == a].dropna(subset=['z']); d = d.assign(y=d.z)
        for c in CNS:
            ST.append(dict(axis=a, lr_pairs='; '.join(p_.replace('^', ' -> ') for p_ in AXES[a]), CN=c, CN_label=LAB[c], CN_name=NM[c], **cn_stats(d, c),
                           mean_raw_score=V[(V.axis == a) & (V.CN == c)].raw.mean()))
    ST = add_bh(pd.DataFrame(ST), 'axis')
    OUT1 = ST.assign(cell_text=[f'{e:+.2f}{stars(q)}' for e, q in zip(ST.mean_within_section_z, ST.lmm_q)],
                     metric='LIANA bivariate local score (raw), mean over the axis pairs per cell, mean over ALL cells of the CN per section; z across CNs within section, mean over sections',
                     model='z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within axis')
    OUT1 = pd.concat([OUT1.assign(row_type='summary'), V.rename(columns={'z': 'within_section_z', 'raw': 'raw_score'}).assign(row_type='section')], ignore_index=True)
    OUT1.apply(sig3).to_csv(C86, index=False)


# ======================================================================================================================
# Section S9A_cn_characterisation: CN cell-type composition, LR modules and inflammatory axes (Supp 9A)
# ======================================================================================================================
if _run('S9A_cn_characterisation'):
    STEM = f'{SUPP9}/supp_cn_niche_characterisation'
    MINUS = lambda t: t.replace('-', '−')
    # composition (12_neighbourhood_treg_score.py)
    E = pd.read_csv(C45, index_col=0)
    LIN = {'Epithelial': ['Epithelial'],
           'Stromal': ['OGN_RSPO3 Fib', 'PLA2G2A Fib', 'ADAMDEC1 Fib', 'Activated Fib', 'Crypt Top Fib', 'Myofibroblast', 'Pericyte', 'Endothelial', 'Neuron'],
           'Myeloid': ['DC', 'CD9 Mac', 'Mac S', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac', 'Neutrophil', 'Granulocyte'],
           'Lymphoid': ['Treg', 'CD4 Trm', 'CD4 Effector', 'Tfh', 'Resting T', 'CD8 Trm', 'Gd T', 'NK', 'B', 'Plasma'],
           'Other': ['Cycling']}
    cols = [c for v in LIN.values() for c in v if c in E.columns]
    H = E[cols].clip(-3, 3)
    DISP = {'OGN_RSPO3 Fib': 'OGN⁺RSPO3⁺ Fib', 'PLA2G2A Fib': 'PLA2G2A⁺ Fib', 'Mac S': 'Mac S+SG+/XS-', 'Mac S_M_P': 'Mac S+M+P+', 'Mac S_M_S': 'Mac S+M+S+', 'Gd T': 'γδ T'}
    XL = [DISP.get(c, c) for c in cols]
    # inflammatory axes and LR modules
    A = pd.read_csv(C86, low_memory=False); S = A[A.row_type == 'summary'].copy(); S['CN'] = S.CN.astype(int)
    AXN = ['IL-6', 'IL-1', 'OSM', 'IL-11', 'TNF']; CONTRAST = 'TNF'
    Z = S.pivot(index='CN', columns='axis', values='mean_within_section_z').loc[CNS, AXN]
    TX = S.pivot(index='CN', columns='axis', values='cell_text').loc[CNS, AXN]
    LAB = S.drop_duplicates('CN').set_index('CN').CN_label.loc[CNS].tolist(); NMF = S.drop_duplicates('CN').set_index('CN').CN_name.loc[CNS]
    ANN = TX.apply(lambda s: s.map(MINUS)).values
    MODS = ['Antigen Presentation', 'Inflammatory Cytokine', 'Co-Stimulatory', 'Matrix–Integrin']
    MODL = ['Antigen\nPresentation', 'Inflammatory\nCytokine', 'Co-\nStimulatory', 'Matrix–\nIntegrin']
    LR = pd.read_csv(C82D, low_memory=False); LS = LR[(LR.row_type == 'summary') & LR.module.isin(MODS)].copy(); LS['CN'] = LS.CN.astype(int)
    ZL = LS.pivot(index='CN', columns='module', values='mean_within_section_z').loc[CNS, MODS]
    TXL = LS.pivot(index='CN', columns='module', values='cell_text').loc[CNS, MODS]
    ANNL = TXL.apply(lambda s: s.map(MINUS)).values
    VM = float(np.ceil(max(np.abs(Z.values).max(), np.abs(ZL.values).max()) * 10) / 10)   # one shared range for both bottom heatmaps
    # layout: composition on top; LR modules | inflammatory axes | TNF contrast below, sharing one colour scale
    FS_TICK, FS_COL, FS_GRP, FS_NUM, FS_CB, FS_TITLE, FS_YLAB = 14, 14, 15, 13, 13, 16, 15
    style(14)
    fig = plt.figure(figsize=(16, 11), layout='constrained')
    gsp = fig.add_gridspec(2, 3, width_ratios=[4 * 1.45, 4 * 0.98, 0.98])
    axc = fig.add_subplot(gsp[0, :]); axl, ax1, ax2 = (fig.add_subplot(gsp[1, i]) for i in range(3))
    sns.heatmap(H, cmap='RdBu_r', center=0, vmin=-3, vmax=3, ax=axc, cbar=False, linewidths=0.4, linecolor='white', xticklabels=True, yticklabels=True)
    axc.set_xticks(np.arange(len(cols)) + 0.5); axc.set_xticklabels(XL, rotation=45, ha='right', va='center', rotation_mode='anchor', fontsize=FS_COL)
    x = 0
    for lin, cs in LIN.items():
        n = len([c for c in cs if c in cols])
        if x: axc.axvline(x, color='k', lw=1.5)
        axc.text(x + n / 2, -0.25, lin, ha='center', va='bottom', fontsize=FS_GRP, fontweight='normal'); x += n
    axc.set_title('Cell-Type Composition', fontsize=FS_TITLE, pad=30)
    kw = dict(cmap='RdBu_r', center=0, vmin=-VM, vmax=VM, fmt='', annot_kws={'fontsize': FS_NUM}, linewidths=0.5, linecolor='white', cbar=False, xticklabels=True)
    sns.heatmap(ZL.set_axis(LAB), annot=ANNL, ax=axl, yticklabels=True, **kw)
    sns.heatmap(Z[AXN[:4]].set_axis(LAB), annot=ANN[:, :4], ax=ax1, yticklabels=False, **kw)
    sns.heatmap(Z[[CONTRAST]].set_axis(LAB), annot=ANN[:, 4:], ax=ax2, yticklabels=False, **kw)
    axl.set_xticks(np.arange(4) + 0.5); axl.set_xticklabels(MODL, rotation=0, fontsize=FS_COL, linespacing=1.0)
    ax1.set_xticks(np.arange(4) + 0.5); ax1.set_xticklabels(AXN[:4], rotation=0, fontsize=FS_COL)
    ax2.set_xticks([0.5]); ax2.set_xticklabels([CONTRAST], rotation=0, fontsize=FS_COL)
    axl.set_title('Ligand–Receptor Modules', fontsize=FS_TITLE, pad=8)
    ax1.set_title('Inflammatory Signalling', fontsize=FS_TITLE, pad=8, x=0.64)   # centred over the IL-6 to TNF columns
    for a_ in (ax1, ax2): a_.set_yticks([]); a_.set_ylabel('')
    for ax in (axc, axl):
        ax.set_yticks(np.arange(10) + 0.5); ax.set_yticklabels(LAB, rotation=0, fontsize=FS_TICK)
        ax.set_ylabel('Spatial Neighborhood', fontsize=FS_YLAB); cn_key(ax, np.arange(10) + 0.5, 0.32)
    for ax in (axc, axl, ax1, ax2): ax.tick_params(axis='x', length=3, width=0.8, pad=2); ax.set_xlabel('')
    for src_ax, par, lab_ in ((axc, axc, 'log2 Enrichment vs All Cells'), (ax2, [axl, ax1, ax2], 'Mean Within-Section z')):
        cb = fig.colorbar(src_ax.collections[0], ax=par, shrink=0.85, aspect=15)
        cb.set_label(lab_, fontsize=FS_CB); cb.ax.tick_params(labelsize=FS_CB); cb.outline.set_visible(False)
    cb.ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1.0 if VM > 1.5 else 0.5)); cb.ax.yaxis.set_major_formatter(umf())
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    # source data: every plotted cell
    rows = []
    for c, lab in zip(CNS, LAB):
        for k, col in enumerate(cols):
            rows.append(dict(heatmap='Cell-Type Composition', CN=c, CN_label=lab, CN_name=NMF.loc[c], column=XL[k], source_column=col,
                             value=E.loc[lab, col], value_plotted=H.loc[lab, col], cell_text='',
                             metric='log2((fraction of the cell type among the CN cells + 1e-4) / (fraction among all cells + 1e-4)); colour clipped to +/-3',
                             model='descriptive', source_csv=os.path.basename(C45)))
    LSi = LS.set_index(['module', 'CN']); Si = S.set_index(['axis', 'CN'])
    KEEP = ['n_sections', 'mean_ci_low', 'mean_ci_high', 'sections_above_0', 'lmm_estimate', 'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q', 'wilcoxon_q']
    for hm, src, keys, tx, srcf, extra in (('Ligand–Receptor Modules', LSi, MODS, TXL, C82D, 'n_pairs'), ('Inflammatory Signalling', Si, AXN, TX, C86, 'lr_pairs')):
        for c, lab in zip(CNS, LAB):
            for k in keys:
                r = src.loc[(k, c)]
                rows.append(dict(heatmap=hm, CN=c, CN_label=lab, CN_name=NMF.loc[c], column=k, source_column=k, value=r.mean_within_section_z,
                                 value_plotted=r.mean_within_section_z, cell_text=tx.loc[c, k], stars=stars(r.lmm_q), **{q_: r[q_] for q_ in KEEP},
                                 n_pairs=r['n_pairs'] if extra == 'n_pairs' else len(r['lr_pairs'].split('; ')),
                                 lr_pairs=(r['lr_pairs'] if isinstance(r.get('lr_pairs'), str) else ''), metric=r.metric,
                                 model=r.model.replace('within module', 'within each module').replace('within axis', 'within each axis'),
                                 source_csv=os.path.basename(srcf)))
    OUTC = pd.DataFrame(rows)
    for c_ in ['CN', 'n_sections', 'sections_above_0', 'n_pairs']: OUTC[c_] = OUTC[c_].round().astype('Int64')
    OUTC['colour_scale'] = np.where(OUTC.heatmap == 'Cell-Type Composition', 'RdBu_r, -3 to 3 (own colourbar)', f'RdBu_r, -{VM:g} to {VM:g} (shared by both bottom heatmaps)')
    OUTC.apply(sig3).to_csv(f'{STEM}.csv', index=False)


def niche_signalling_table():
    """The 8 signalling columns (modules + axes) of the S9A table, used by Fig 5C and Supp 9B."""
    N89 = pd.read_csv(C89); N89 = N89[N89.heatmap.isin(['Ligand–Receptor Modules', 'Inflammatory Signalling'])]
    X = pd.concat([N89[(N89.heatmap == SRC[c]) & (N89.column == c)] for c in COLS], ignore_index=True); X['CN'] = X.CN.astype(int)
    return X, X.drop_duplicates('CN').set_index('CN').CN_label.loc[CNS].tolist()


# ======================================================================================================================
# Section F5C_score_signalling: Treg Suppression Score per CN (forest, from 12) beside the niche signalling heatmap (Fig 5C; colours saturate at +/-1.2)
# ======================================================================================================================
if _run('F5C_score_signalling'):
    STEM = f'{FIG5}/fig4_cn_treg_score_and_signaling_vm12_wide'; VM = 1.2
    EMPH_CN = (4, 5, 10)   # CN4, CN5 and CN10 row labels in bold
    X, LAB = niche_signalling_table()
    V = X.pivot(index='CN', columns='column', values='value').loc[CNS, COLS]
    ANN = X.pivot(index='CN', columns='column', values='cell_text').loc[CNS, COLS].apply(lambda s: s.str.replace('-', '−'))
    F = pd.read_csv(C26).set_index('CN').loc[CNS]; EST = 'treg_suppression_score_z_in_minus_out'
    FCOL = {cn: ('#9e9e9e' if not (F.loc[cn, 'mixed_q'] < 0.05) else '#c26a77' if F.loc[cn, EST] > 0 else '#5b7fa6') for cn in CNS if F.loc[cn, 'tested']}
    FST = {cn: stars(F.loc[cn, 'mixed_q']) for cn in CNS if F.loc[cn, 'tested']}
    # layout: forest | heatmap (column groups separated by gaps; x data unit = cell width scale) | colourbar
    FS_LAB, FS_TITLE, FS_CB, FS_TICK, FS_STAR, FS_NT, FS_NUM = 14, 16, 13, 13, 18, 11, 11
    style(FS_LAB)
    CW, GAPG, GAPT, ROTA, XPAD = 0.88, 0.20, 0.13, 28, 0.12   # cell width, gaps between groups / before TNF, label angle, forest room for stars
    XLAB = {'Antigen Presentation': 'Antigen\nPresentation', 'Co-Stimulatory': 'Co-\nStimulatory', 'Matrix–Integrin': 'Matrix–\nIntegrin'}
    XL = {}; x = 0.0
    for gi, (g, cs) in enumerate(GROUPS):
        if gi: x += GAPG
        for c in cs:
            if c == 'TNF': x += GAPT
            XL[c] = x; x += CW
    WT = x
    fig, (axf, axh) = plt.subplots(1, 2, figsize=(17, 6.2), width_ratios=[2.9, WT], layout='constrained')
    # forest panel
    yc = np.arange(10) + 0.5
    axf.set_ylim(10, 0); axf.axvline(0, color='0.4', lw=1, ls='--', zorder=1)
    for i, cn in enumerate(CNS):
        if not F.loc[cn, 'tested']: continue
        lo, hi, e = F.loc[cn, 'ci_lo'], F.loc[cn, 'ci_hi'], F.loc[cn, EST]
        axf.hlines(yc[i], lo, hi, color='0.35', lw=2.0, zorder=2)
        axf.scatter([e], [yc[i]], c=[FCOL[cn]], s=150, zorder=3, edgecolor='k', lw=0.7)
        if FST[cn]: axf.text(hi + 0.02 if e >= 0 else lo - 0.02, yc[i], FST[cn], ha='left' if e >= 0 else 'right', va='center', fontsize=FS_STAR)
    axf.set_yticks(yc); axf.set_yticklabels(LAB, fontsize=FS_LAB); cn_key(axf, yc, 0.32)
    for i, cn in enumerate(CNS):
        if cn in EMPH_CN: axf.get_yticklabels()[i].set_fontweight('bold')
    axf.set_ylabel('Spatial Neighborhood', fontsize=FS_LAB)
    axf.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([-0.4, -0.2, 0.0, 0.2, 0.4])); axf.xaxis.set_major_formatter(umf())
    axf.tick_params(axis='x', labelsize=FS_TICK)
    axf.set_xlabel('Treg Suppression Score (z)\n← Low-suppressive   Activated →', fontsize=FS_TICK, linespacing=1.3)
    for k_ in ('top', 'right', 'left'): axf.spines[k_].set_visible(False)
    axf.spines['bottom'].set_position(('outward', 6))
    x0, x1 = min(-0.45, F.loc[F.tested, 'ci_lo'].min() - 0.02) - XPAD, max(0.45, F.loc[F.tested, 'ci_hi'].max() + 0.02) + XPAD   # x range: ticks, CIs and stars
    axf.set_xlim(x0, x1)
    for i, cn in enumerate(CNS):
        if not F.loc[cn, 'tested']:
            axf.text((x0 + x1) / 2, yc[i], 'Not tested (<100 Tregs)', ha='center', va='center', fontsize=FS_NT, color='0.45', style='italic',
                     bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'), zorder=4)
    # heatmap
    axh.set_xlim(0, WT); axh.set_ylim(10, 0)
    for sp in axh.spines.values(): sp.set_visible(False)
    cmap = plt.get_cmap('RdBu_r'); norm = matplotlib.colors.Normalize(-VM, VM)
    for c in COLS:
        for i, cn in enumerate(CNS):
            col = cmap(norm(V.loc[cn, c]))
            axh.add_patch(Rectangle((XL[c], i), CW, 1, facecolor=col, edgecolor='white', lw=0.5))
            axh.text(XL[c] + CW / 2, i + 0.5, ANN.loc[cn, c], ha='center', va='center', fontsize=FS_NUM, color='.15' if relative_luminance(col) > .408 else 'w')
    axh.set_yticks([]); axh.set_xticks([XL[c] + CW / 2 for c in COLS])
    axh.set_xticklabels([XLAB.get(c, c) for c in COLS], fontsize=FS_LAB, rotation=ROTA, ha='right', va='top', rotation_mode='anchor', linespacing=1.0)
    axh.tick_params(axis='x', length=0, pad=5)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axh, shrink=0.85, aspect=15, extend='both')
    cb.set_label('Mean Within-Section z', fontsize=FS_CB); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=FS_CB)
    cb.ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5)); cb.ax.yaxis.set_major_formatter(umf())

    def header(ax, xcrd, a, b, texts):   # rule from a to b just above the axes, then (text, fontsize, offset in pt) above it
        tr = blended_transform_factory(xcrd, ax.transAxes)
        ax.plot([a, b], [1, 1], color='k', lw=1.0, solid_capstyle='butt', clip_on=False, transform=tr + ScaledTranslation(0, 0.084, fig.dpi_scale_trans))
        for t_, fs, dy in texts: ax.annotate(t_, ((a + b) / 2, 1), xycoords=tr, xytext=(0, dy), textcoords='offset points', ha='center', va='bottom', fontsize=fs)
    for g, cs in GROUPS: header(axh, axh.transData, XL[cs[0]] + 0.04, XL[cs[-1]] + CW - 0.04, [(g, FS_LAB, 9)])
    header(axf, axf.transAxes, 0.014, 0.986, [('In vs Outside Neighborhood', FS_TICK, 9), ('Treg Suppression Score', FS_LAB, 27)])
    fig.suptitle('Treg Suppression Score and Niche Signaling', fontsize=FS_TITLE)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    # source data: forest rows + heatmap cells
    FO = F.reset_index().rename(columns={EST: 'estimate', 'ci_lo': 'ci_low', 'ci_hi': 'ci_high', 'mixed_p': 'p', 'mixed_q': 'q'})
    FO = FO.assign(panel='Treg Suppression Score (forest)', column='Treg Suppression Score, in vs outside neighborhood', CN_label=LAB,
                   stars=[FST.get(c, '') for c in CNS], color=[FCOL.get(c, 'not tested') for c in CNS],
                   model='Treg Suppression Score (z within section) ~ in_CN + (1 | section); BH over the 8 tested neighborhoods', source_csv=os.path.basename(C26))
    FO = FO[['panel', 'CN', 'CN_label', 'CN_name', 'column', 'estimate', 'ci_low', 'ci_high', 'p', 'q', 'stars', 'tested', 'n_tregs_in_CN', 'n_sections_ge10_tregs',
             'sections_agreeing', 'sections', 'color', 'model', 'source_csv']]
    KEEP = ['CN', 'CN_label', 'CN_name', 'column', 'value', 'cell_text', 'stars', 'n_sections', 'mean_ci_low', 'mean_ci_high', 'sections_above_0', 'lmm_estimate',
            'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q', 'wilcoxon_q', 'n_pairs', 'lr_pairs', 'metric', 'model', 'source_csv']
    HO = X[KEEP].rename(columns={'value': 'mean_within_section_z'}).copy(); HO.insert(0, 'panel', 'Niche Signaling (heatmap)'); HO.insert(1, 'group', HO.column.map(GOF))
    HO['column'] = pd.Categorical(HO.column, COLS); HO = HO.sort_values(['column', 'CN']); HO['column'] = HO.column.astype(str)
    HO['source_csv'] = HO.source_csv + ' (via ' + os.path.basename(C89) + ')'
    OUTC = pd.concat([FO, HO], ignore_index=True); OUTC['highlighted'] = OUTC.CN == 10
    for c_ in ['CN', 'n_tregs_in_CN', 'n_sections_ge10_tregs', 'sections_agreeing', 'sections', 'n_sections', 'sections_above_0', 'n_pairs']: OUTC[c_] = OUTC[c_].round().astype('Int64')
    OUTC.apply(sig3).to_csv(f'{STEM}.csv', index=False)


# ======================================================================================================================
# Section S9B_cn_signalling: per-niche mixed-model estimate and 95% CI for every column of the niche signalling heatmap (Supp 9B)
# ======================================================================================================================
if _run('S9B_cn_signalling'):
    STEM = f'{SUPP9}/supp_cn_signaling_errorbars'
    HI = '#b2182b'; GREY = '#6e6e6e'; PTC = '#252525'
    PANEL_HI = {'Antigen Presentation': 4, 'Co-Stimulatory': 4, 'IL-6': 10, 'IL-1': 10, 'OSM': 10, 'IL-11': 10, 'TNF': 4, 'Matrix–Integrin': 10}   # niche discussed in the text
    LMMC = ['lmm_estimate', 'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q']
    X, LAB = niche_signalling_table()
    SUMR = []
    for c in COLS:
        hc = X[X.column == c].set_index('CN').loc[CNS]
        for cn in CNS:
            h = hc.loc[cn]; st = h.stars if isinstance(h.stars, str) else ''
            e, lo, hi, p, q = (float(h[k_]) for k_ in LMMC)
            ratio = e / h.value if abs(h.value) > 0.01 else np.nan
            SUMR.append(dict(group=GOF[c], column=c, CN=cn, CN_label=LAB[cn - 1], CN_name=h.CN_name, n_sections=int(h.n_sections),
                             lmm_estimate=e, lmm_ci_low=lo, lmm_ci_high=hi, lmm_p=p, lmm_q=q, stars=st, filled_point=bool(q < 0.05), ci_excludes_0=bool(lo > 0 or hi < 0),
                             heatmap_mean_within_section_z=h.value, estimate_to_heatmap_ratio=ratio, heatmap_cell_text=h.cell_text,
                             source_csv=os.path.basename(C89) + ' = ' + h.source_csv))
    SUMR = pd.DataFrame(SUMR)
    # layout: 4 x 2 panels, row labels and CN colour key on the left column, stars in one column right of all CIs
    FS_LAB, FS_TICK, FS_TITLE, FS_GROUP, FS_STAR, FS_KEY = 13, 13, 14, 13, 14, 13
    style(FS_TICK)
    MS, MEW, ELW, MS_HI = 7.5, 1.6, 1.9, 9.5
    FP = FontProperties(family='DejaVu Sans', size=FS_STAR)
    GLY = {s: TextPath((0, 0), s, size=FS_STAR, prop=FP).get_extents() for s in ['*', '**', '***']}   # glyph extents, to centre the stars on the row
    XLO = SUMR.lmm_ci_low.min() - 0.15
    xc = max(SUMR.lmm_ci_high.max(), SUMR.lmm_estimate.max()) + 0.35; XHI = xc + 0.3
    fig, AXS = plt.subplots(4, 2, figsize=(11, 12.3), sharex=True, layout='constrained')
    for k, (c, ax) in enumerate(zip(COLS, AXS.ravel())):
        left = k % 2 == 0
        ax.set_xlim(XLO, XHI); ax.set_ylim(9.55, -0.55)
        ax.axvline(0, ls=(0, (4, 3)), color='#555555', lw=1.0, zorder=1)
        h_ = SUMR[SUMR.column == c].set_index('CN')
        for i, cn in enumerate(CNS):
            r = h_.loc[cn]
            ax.hlines(i, r.lmm_ci_low, r.lmm_ci_high, color=PTC, lw=ELW, zorder=3)
            hc_ = PTC if PANEL_HI[c] != cn else HI
            ax.plot(r.lmm_estimate, i, 'o', ms=MS if hc_ == PTC else MS_HI, color=hc_, mfc=hc_ if r.filled_point else 'white', mew=MEW, zorder=4 if hc_ == PTC else 5)
            if r.stars:
                gb = GLY[r.stars]
                ax.annotate(r.stars, xy=(xc, i), xytext=(-(gb.x0 + gb.x1) / 2, -(gb.y0 + gb.y1) / 2), textcoords='offset points',
                            ha='left', va='baseline', fontsize=FS_STAR, color='k', annotation_clip=False)
        ax.set_title(c, loc='left', fontsize=FS_TITLE, pad=4)
        ax.set_title(GOF[c], loc='right', fontsize=FS_GROUP, color=GREY, pad=4)
        ax.set_yticks(range(10)); ax.tick_params(axis='y', length=0)
        if left:   # shared row labels; red for the niche highlighted in the row's left panel
            ax.set_yticklabels(LAB, fontsize=FS_LAB); cn_key(ax, range(10), 0.19)
            ax.get_yticklabels()[CNS.index(PANEL_HI[c])].set_color(HI)
        else:
            ax.set_yticklabels([])
        ax.set_xticks([-2, -1, 0, 1, 2]); ax.tick_params(axis='x', labelsize=FS_TICK, length=3)
        ax.spines[['top', 'right', 'left']].set_visible(False)
    fig.supxlabel('Niche vs Other Niches (Mixed-Model Estimate, z)', fontsize=FS_TITLE)
    key = [Line2D([], [], marker='o', ls='-', color=PTC, mfc=PTC, ms=MS, mew=MEW, lw=ELW, label='Estimate ± 95% CI, q < 0.05'),
           Line2D([], [], marker='o', ls='-', color=PTC, mfc='white', ms=MS, mew=MEW, lw=ELW, label='q ≥ 0.05'),
           Line2D([], [], ls=(0, (4, 3)), color='#555555', lw=1.0, label='0 = No Difference')]
    fig.legend(handles=key, loc='outside upper right', ncol=4, frameon=False, fontsize=FS_KEY, handletextpad=0.4, columnspacing=1.2, handlelength=1.8)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    OUTC = SUMR.copy(); OUTC['highlighted'] = OUTC.CN == 10
    OUTC = OUTC.assign(metric='per section: mean LIANA bivariate local score over the column\'s pairs and over ALL cells of the niche (>= 100 cells), z across the '
                              '10 niches within the section; PLOTTED: lmm_estimate (niche minus the other niches, within-section z) with its Wald 95% CI',
                       model='z ~ in_CN + (1 | section) on per-section z; BH over the 10 niches within each column (values from the heatmap; not recomputed)')
    for c_ in ['CN', 'n_sections']: OUTC[c_] = OUTC[c_].round().astype('Int64')
    OUTC.apply(sig3).to_csv(f'{STEM}.csv', index=False)
