#!/usr/bin/env python
"""CosMx spatial proteomics (RNA and a 64-antibody protein panel in the same cells): protein-defined lineage maps, RNA lineage vs
canonical protein, and PD-L1+ macrophage / PD-1+ Treg proximity in CD ileum.

Panels: Supp 6D (UC rectum protein-defined lineage map), Supp 6E (UC RNA lineage vs canonical protein dot plot),
        Supp 9I (CD ileum protein-defined lineage map), Supp 9J (CD RNA lineage vs canonical protein dot plot),
        Supp 11D (CD ileum PD-L1+ macrophage / PD-1+ Treg map with the proximity statistics).
Adapted from revision scripts analysis/protein_lr_proximity.py (statistics) and supp_spatial_multiomics.py (figure and panels).
Sections:
  S11D_protein_lr_proximity: CD ileum 375i, 17 FOVs. Macrophages = RNA SELENOP+ Mac or CD9+ M2 Mac, Tregs = RNA Treg (rna_qc_pass);
    PD-L1+ / PD-1+ = protein > p75 over all cells of the section. Per FOV, observed = sum over PD-L1+ macrophages of PD-1+ Tregs within
    r um (same FOV); null = PD-L1+ status permuted among the FOV's macrophages (1,000 permutations, seed 0); O/E, two-sided permutation
    P; pooled O/E with FOV bootstrap CI and stratified permutation P; 30 um primary, 50 um sensitivity. Also protein-gate sensitivity
    rows, PD-L1+ Mac near CTLA4+ Tregs and HLA-DR-high Mac -> CTLA4+ / CTLA4- Tregs (CSV only).
  S6DE_S9IJ_S11D_figure: per-FOV dot plot of the RNA lineage minus all other cells (protein z) for 12 lineage-canonical protein pairs
    (two-sided Wilcoxon across FOVs, BH over 12 per section); full-section maps coloured by protein-defined lineage (within-section
    percentile gates, one label per cell by z-order); PD-L1+ macrophage / PD-1+ Treg map (colour = within-section protein percentile,
    near = partner within 30 um in the same FOV) with the median per-FOV O/E and the
    permutation P of S11D_protein_lr_proximity.
Inputs: UC rectum 385R RNA1_Cos4_385R_381R_S2_norm_wpr_anno_jd.h5ad (03); CD ileum rna_norm_pr_geonorm_anno_aligned.h5ad (06);
  fig5_spatial_maps_v7.csv (gate counts of the earlier protein figure; fills the n_gate_v7 column only).
Outputs (/share/fsmresfiles/UC/revision_jenny/SuppFig6D-E_9I-J_11D_spatial_proteomics/): analysis/protein_lr_proximity.csv; supp_spatial_multiomics.png, .csv;
  panels/uc_rectum_protein_annotation_map.png (Supp 6D), panels/uc_rectum_rna_lineage_protein_dotplot.png (Supp 6E),
  panels/cd_ileum_protein_annotation_map.png (Supp 9I), panels/cd_ileum_rna_lineage_protein_dotplot.png (Supp 9J),
  panels/cd_ileum_pdl1_mac_pd1_treg_spatial_map.png (Supp 11D), panels/protein_lineage_key_scale_bar.png (key of the lineage maps).
Run: python 07_protein_lineage_maps_lr_proximity.py, or RUN_SECTIONS=<id>[,<id>...] for single sections (file order).
"""
import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
warnings.filterwarnings('ignore')

_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


OUT = '/share/fsmresfiles/UC/revision_jenny/SuppFig6D-E_9I-J_11D_spatial_proteomics'
PANELS = f'{OUT}/panels'
UC_PATH = '/share/fsmresfiles/UC/scOmics_model/data/Spatial_Protein_CosMx'
CD_PATH = '/share/fsmresfiles/UC/AtoMx/CD_multiOmics/Cos4_315i_375i'
REF7 = '/share/fsmresfiles/UC/revision_jenny/_shared_inputs/Fig5_spatial_protein/fig5_spatial_maps_v7.csv'
LRCSV = f'{OUT}/analysis/protein_lr_proximity.csv'
FOVS_CD = [6, 7, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29]
px_per_um = np.sqrt(3266 / 47.25)
RADS, NPERM, NBOOT, SEED, MIN_PAIRS = (30, 50), 1000, 10000, 0, 15       # proximity test
N_BOOT, BOOT_SEED = 2000, 20261116                                        # bootstrap of the per-FOV median O/E (CSV only)
sig3 = lambda v: float(f'{v:.3g}') if pd.notna(v) and np.isfinite(v) else v

# ── protein gates (within-section percentiles; shared by both sections) ──
GATE_MARKERS = ['CD45', 'CD3', 'CD20', 'CD19', 'CD68', 'CD163', 'CD11c', 'CD14', 'CD15', 'CD16', 'CD123']


def gate_myeloid(adata):
    v = {m: adata.obs[f'{m}_geo'].values.astype(float) for m in GATE_MARKERS}
    thr = {m: {'p75': np.percentile(v[m], 75), 'p90': np.percentile(v[m], 90)} for m in GATE_MARKERS}
    pos = {m: v[m] > thr[m]['p75'] for m in GATE_MARKERS}
    high = {m: v[m] > thr[m]['p90'] for m in GATE_MARKERS}
    parent = (pos['CD45'] & ~pos['CD3'] & ~pos['CD20'] & ~pos['CD19']
              & (pos['CD68'] | pos['CD163'] | pos['CD11c'] | pos['CD14']))
    lab = np.full(adata.n_obs, 'non-myeloid', dtype=object)
    free = parent.copy()

    def assign(mask, name):
        nonlocal free
        m = free & mask
        lab[m] = name
        free = free & ~m
    assign(pos['CD15'] & pos['CD16'], 'Neutrophil')
    assign(high['CD123'] & ~pos['CD11c'], 'pDC')
    assign((pos['CD68'] | pos['CD163']) & ~high['CD11c'], 'Macrophage')
    assign(high['CD11c'] & ~pos['CD68'] & ~pos['CD14'] & ~pos['CD163'], 'DC')
    assign(np.ones(adata.n_obs, bool), 'Other myeloid')
    return lab, pos


def masks_of(adata):
    lab, pos = gate_myeloid(adata)
    g = lambda p: adata.obs[f'{p}_geo'].values.astype(float)
    above = lambda p, q: g(p) > np.percentile(g(p), q)
    notabove = lambda p, q: g(p) <= np.percentile(g(p), q)
    epi = above('EpCAM', 75)
    linneg = notabove('CD45', 75) & notabove('EpCAM', 75) & notabove('CD31', 75)
    stroma = above('Fibronectin', 60) & above('Vimentin', 50) & notabove('SMA', 85) & linneg
    mac = lab == 'Macrophage'
    dc_strict = lab == 'DC'
    dc_relax = (above('CD11c', 75) & above('HLA-DR', 75) & notabove('CD68', 75) & notabove('CD14', 75)
                & ~pos['CD3'] & ~pos['CD20'] & ~pos['CD19'] & ~mac)
    bcell = (pos['CD20'] | pos['CD19']) & ~pos['CD3']
    treg = pos['CD3'] & above('CD4', 75) & above('FOXP3', 75) & above('CTLA4', 75)
    oth_t = pos['CD3'] & ~treg
    return dict(Epithelium=epi, Stroma=stroma, Macrophage=mac, **{'DC strict': dc_strict,
                'DC relaxed': dc_relax, 'B cell': bcell, 'Treg': treg, 'Other T': oth_t})


# ======================================================================================================================
# Section S11D_protein_lr_proximity: PD-L1+ macrophage -> PD-1+ Treg proximity (per FOV, pooled), CD ileum 375i, 17 FOVs.
# Output: analysis/protein_lr_proximity.csv (read by the figure section)
# ======================================================================================================================
if _run('S11D_protein_lr_proximity'):
    STEM = f'{OUT}/analysis/protein_lr_proximity'
    # ── Load: CD ileum 375i, 17-FOV subset (verbatim from cd_verification.py) ──
    A = sc.read_h5ad(f'{CD_PATH}/rna_norm_pr_geonorm_anno_aligned.h5ad')
    A = A[A.obs['fov'] <= 29].copy()
    fovs = [6, 7, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29]
    A = A[A.obs['fov'].isin(fovs)].copy()

    M = masks_of(A)

    ob = A.obs
    qc = ob['rna_qc_pass'].values.astype(bool)
    fov = ob['fov'].values.astype(int)
    ctg = ob['cell_type_general'].astype(str).values
    xy = np.asarray(A.obsm['spatial_fov'])[:, :2].astype(float) / px_per_um
    above75 = lambda p: (lambda v: v > np.percentile(v, 75))(ob[f'{p}_geo'].values.astype(float))
    PDL1, PD1, CTLA4, HLADR = above75('PD-L1'), above75('PD-1'), above75('CTLA4'), above75('HLA-DR')
    LIG = {'PD-L1': PDL1, 'HLA-DR': HLADR}

    MAC = {'RNA': qc & np.isin(ctg, ['SELENOP+ Mac', 'CD9+ M2 Mac']), 'protein': M['Macrophage']}
    TREG = {'RNA': qc & (ctg == 'Treg'), 'protein': M['Treg']}
    REC = {'PD-1': PD1, 'CTLA4': CTLA4, 'CTLA4-neg': ~CTLA4}
    # (pair label, ligand, receptor, mac def, treg def, definition role, readout note)
    LR = 'ligand-receptor pair'
    EXP = 'exploratory receptor readout; CTLA4 does not bind PD-L1, not a ligand-receptor pair'
    HD = 'HLA-DR-high Macrophage (Activated APC Proxy)'
    HDN = 'activated-APC proxy for the CTLA4 axis (CD80/CD86 not on the protein panel)'
    HDC = 'specificity control for the activated-APC proxy: CTLA4- RNA-labelled Tregs'
    TESTS = [('PD-L1 (Mac) -> PD-1 (Treg)', 'PD-L1', 'PD-1', 'RNA', 'RNA', 'primary', LR),
             ('PD-L1 (Mac) -> PD-1 (Treg)', 'PD-L1', 'PD-1', 'protein', 'protein', 'sensitivity (protein gates)', LR),
             ('PD-L1 (Mac) near CTLA4+ Treg', 'PD-L1', 'CTLA4', 'RNA', 'RNA', 'primary', EXP),
             ('PD-L1 (Mac) near CTLA4+ Treg', 'PD-L1', 'CTLA4', 'protein', 'RNA', 'sensitivity (protein macrophage gate)', EXP),
             (f'{HD} -> CTLA4+ Treg', 'HLA-DR', 'CTLA4', 'RNA', 'RNA', 'primary', HDN),
             (f'{HD} -> CTLA4+ Treg', 'HLA-DR', 'CTLA4', 'protein', 'RNA', 'sensitivity (protein macrophage gate)', HDN),
             (f'{HD} near CTLA4- Treg (control)', 'HLA-DR', 'CTLA4-neg', 'RNA', 'RNA', 'primary', HDC),
             (f'{HD} near CTLA4- Treg (control)', 'HLA-DR', 'CTLA4-neg', 'protein', 'RNA', 'sensitivity (protein macrophage gate)', HDC)]

    ROWS, POOL = [], []
    for ti, (pair, lig, rec, md, td, role, note) in enumerate(TESTS):
        mac, recT, LP = MAC[md], TREG[td] & REC[rec], LIG[lig]
        for r in RADS:
            per = []
            for fi, f in enumerate(fovs):
                inf = fov == f
                mi, ri = np.flatnonzero(mac & inf), np.flatnonzero(recT & inf)
                L = LP[mi]; k, nM = int(L.sum()), len(mi)
                o = dict(pair=pair, ligand=lig, receptor=rec, readout=note, definition=role, macrophage_def=md, treg_def=td, radius_um=r, row='FOV', fov=f,
                         n_macrophages=nM, n_ligand_pos_macrophages=k, n_tregs=int((TREG[td] & inf).sum()), n_receptor_pos_tregs=len(ri))
                if k == 0 or len(ri) == 0:
                    per.append({**o, 'obs_pairs': 0 if nM else np.nan, 'note': 'not tested: no ligand+ macrophages or receptor+ Tregs'}); continue
                W = np.array([len(x) for x in cKDTree(xy[ri]).query_ball_point(xy[mi], r=r)], float)   # receptor+ Tregs within r of each mac
                obs = float(W[L].sum()); rng = np.random.default_rng([SEED, ti, fi, r])
                perm = np.array([W[rng.permutation(nM)[:k]].sum() for _ in range(NPERM)])
                ex = perm.mean()
                ph, pl = (1 + (perm >= obs).sum()) / (NPERM + 1), (1 + (perm <= obs).sum()) / (NPERM + 1)
                o.update(obs_pairs=obs, exp_pairs=ex, obs_exp=obs / ex if ex > 0 else np.nan, perm_p=min(1.0, 2 * min(ph, pl)) if ex > 0 else np.nan,
                         n_ligand_pos_macs_near_receptor_pos_treg=int((W[L] > 0).sum()), _dev=obs - ex, _null=perm - ex)
                if ex == 0: o['note'] = 'not tested: no receptor+ Treg within radius of any macrophage'
                per.append(o)
            T = pd.DataFrame(per); t = T[T.obs_exp.notna()]
            allm = mac; allr = recT
            p = dict(pair=pair, ligand=lig, receptor=rec, readout=note, definition=role, macrophage_def=md, treg_def=td, radius_um=r, row='pooled',
                     fov=f'pooled ({len(fovs)} FOVs)', n_macrophages=int(allm.sum()), n_ligand_pos_macrophages=int((allm & LP).sum()),
                     n_tregs=int(TREG[td].sum()), n_receptor_pos_tregs=int(allr.sum()), fovs_tested=len(t))
            if len(t):
                obs_t, exp_t = t.obs_pairs.sum(), t.exp_pairs.sum()
                st = t._dev.sum(); nl = np.sum(np.vstack(t._null.values), axis=0)
                brng = np.random.default_rng([SEED, ti, r, 99]); oo, ee = t.obs_pairs.values, t.exp_pairs.values
                bi = brng.integers(0, len(t), (NBOOT, len(t))); boot = oo[bi].sum(1) / ee[bi].sum(1)
                up, nt = int((t.obs_exp > 1).sum()), int((t.obs_exp != 1).sum())
                p.update(obs_pairs=obs_t, exp_pairs=exp_t, obs_exp=obs_t / exp_t, ci95_lo=np.percentile(boot, 2.5), ci95_hi=np.percentile(boot, 97.5),
                         perm_p=min(1.0, 2 * min((1 + (nl >= st).sum()) / (NPERM + 1), (1 + (nl <= st).sum()) / (NPERM + 1))),
                         fovs_obs_exp_gt1=up, sign_test_p=stats.binomtest(up, nt, 0.5).pvalue if nt else np.nan,
                         n_ligand_pos_macs_near_receptor_pos_treg=int(t.n_ligand_pos_macs_near_receptor_pos_treg.sum()))
            else:
                p.update(obs_pairs=0, note='no FOV testable')
            p['underpowered'] = bool(p.get('obs_pairs', 0) < MIN_PAIRS)
            p['supports_proximity'] = bool(len(t) and p['perm_p'] < 0.05 and p['obs_exp'] > 1 and p['fovs_obs_exp_gt1'] > len(t) / 2)
            POOL.append(p); ROWS += per

    P = pd.DataFrame(POOL); S = pd.DataFrame(ROWS).drop(columns=['_dev', '_null'])
    ALL = pd.concat([P, S], ignore_index=True)
    lead = ['pair', 'readout', 'definition', 'macrophage_def', 'treg_def', 'radius_um', 'row', 'fov', 'obs_exp', 'ci95_lo', 'ci95_hi', 'perm_p',
            'fovs_obs_exp_gt1', 'fovs_tested', 'sign_test_p', 'obs_pairs', 'exp_pairs', 'n_ligand_pos_macrophages', 'n_receptor_pos_tregs',
            'n_ligand_pos_macs_near_receptor_pos_treg', 'n_macrophages', 'n_tregs', 'underpowered', 'supports_proximity']
    ALL = ALL[lead + [c for c in ALL.columns if c not in lead]]
    for c in ALL.columns:
        if c.startswith('n_') or c in ('fovs_obs_exp_gt1', 'fovs_tested', 'radius_um', 'obs_pairs'):
            ALL[c] = pd.to_numeric(ALL[c]).round().astype('Int64')
        elif ALL[c].dtype.kind == 'f':
            ALL[c] = ALL[c].map(sig3)
    ALL['notes'] = ('CD ileum 375i, 17 FOVs; protein + = > p75 over all section cells; RNA cells = rna_qc_pass; observed = ligand+ macrophage - '
                    'receptor+ Treg pairs within radius (same FOV); null = ligand+ (PD-L1 or HLA-DR > p75) status permuted among macrophages within FOV, 1,000 permutations; '
                    'pooled O/E = sum obs / sum exp over tested FOVs; CI = FOV bootstrap 10,000; perm P two-sided, stratified; underpowered = < 15 pairs')
    with open(f'{STEM}.csv', 'w') as fh:
        ALL.to_csv(fh, index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    print(P[['pair', 'definition', 'radius_um', 'n_ligand_pos_macrophages', 'n_receptor_pos_tregs', 'obs_pairs', 'exp_pairs', 'obs_exp', 'ci95_lo',
             'ci95_hi', 'perm_p', 'fovs_obs_exp_gt1', 'fovs_tested', 'n_ligand_pos_macs_near_receptor_pos_treg', 'underpowered',
             'supports_proximity']].to_string())
    print(S[S.radius_um == 30][['pair', 'definition', 'fov', 'n_ligand_pos_macrophages', 'n_receptor_pos_tregs', 'obs_pairs', 'exp_pairs',
                                'obs_exp']].to_string())

# ======================================================================================================================
# Section S6DE_S9IJ_S11D_figure: RNA lineage vs canonical protein dot plots, protein-defined lineage maps (UC, CD) and the
# PD-L1+ Mac / PD-1+ Treg map; combined PNG + CSV and one PNG per panel in panels/
# ======================================================================================================================
if _run('S6DE_S9IJ_S11D_figure'):
    STEM = f'{OUT}/supp_spatial_multiomics'
    # ── Load ──
    U = sc.read_h5ad(f'{UC_PATH}/RNA1_Cos4_385R_381R_S2_norm_wpr_anno_jd.h5ad')
    C = sc.read_h5ad(f'{CD_PATH}/rna_norm_pr_geonorm_anno_aligned.h5ad')
    C = C[(C.obs['fov'] <= 29) & C.obs['fov'].isin(FOVS_CD)].copy()
    n_ab = sum(c.endswith('_geo') and not c.startswith('Channel-') for c in U.obs.columns)
    print('antibodies (incl. 2 isotype controls):', n_ab)

    zs = lambda v: (v - v.mean()) / v.std()
    CANON = [('Epithelium', 'EpCAM'), ('Fibroblast/Stroma', 'Vimentin'), ('Fibroblast/Stroma', 'SMA'),
             ('Endothelium', 'CD31'), ('Macrophage/Myeloid', 'CD68'), ('Macrophage/Myeloid', 'CD163'), ('DC', 'HLA-DR'),
             ('DC', 'CD11c'), ('T Cells', 'CD3'), ('T Cells', 'CD4'), ('B/Plasma Cells', 'CD20'), ('B/Plasma Cells', 'CD138')]
    MARK = sorted({m for _, m in CANON})
    LIN = ['Epithelium', 'Fibroblast/Stroma', 'Endothelium', 'Macrophage/Myeloid', 'DC', 'T Cells', 'B/Plasma Cells']
    LIN6 = [x for x in LIN if x != 'Endothelium']          # lineages with a protein gate
    UC_DC = ['cDC1', 'cDC2', 'CCR7+ DC']


    def lineages(cat, dc, neut, keep):
        return {'Epithelium': keep & (cat == 'Epithelial'), 'Fibroblast/Stroma': keep & (cat == 'Connective'),
                'Endothelium': keep & (cat == 'Endothelial'),
                'Macrophage/Myeloid': keep & (cat == 'Myeloid') & ~dc & ~neut, 'DC': keep & dc,
                'T Cells': keep & (cat == 'T'), 'B/Plasma Cells': keep & np.isin(cat, ['B', 'Plasma'])}


    S = {}
    o = U.obs
    ct, cat = o['aucell_cell_type_short'].astype(str).values, o['aucell_cell_category_from_type'].astype(str).values
    keep = np.ones(U.n_obs, bool)
    S['UC'] = dict(ad=U, obs=o, fov=o['fov'].values.astype(int), keep=keep, cat=cat,
                   L=lineages(cat, np.isin(ct, UC_DC), ct == 'Neutrophil', keep))
    o = C.obs
    ctg, cat = o['cell_type_general'].astype(str).values, o['cell_category'].astype(str).values
    keep = o['rna_qc_pass'].values.astype(bool)
    S['CD'] = dict(ad=C, obs=o, fov=o['fov'].values.astype(int), keep=keep, cat=cat,
                   L=lineages(cat, ctg == 'DC', ctg == 'Neutrophil', keep))
    for s in S.values():
        s['PZ'] = {m: zs(s['obs'][f'{m}_geo'].values.astype(float)) for m in MARK}
        s['xy'] = np.asarray(s['ad'].obsm['spatial_fov'])[:, :2].astype(float) / px_per_um   # global slide coordinates, um


    def fov_diffs(v, g1, g0, fov):
        rows = [(f, v[g1 & (fov == f)].mean(), v[g0 & (fov == f)].mean()) for f in np.unique(fov)
                if (g1 & (fov == f)).sum() >= 3 and (g0 & (fov == f)).sum() >= 3]
        return pd.DataFrame(rows, columns=['fov', 'mean_group', 'mean_ref'])


    def wtest(d):
        x = (d.mean_group - d.mean_ref).values
        p = stats.wilcoxon(x).pvalue if len(x) >= 3 and np.any(x != 0) else np.nan
        return dict(n_FOV=len(x), n_FOV_higher=int((x > 0).sum()), estimate=np.median(x) if len(x) else np.nan, P=p)


    sig_of = lambda q: '' if not np.isfinite(q) else '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'

    # ── top: canonical lineage-marker pairs (computation unchanged from v2 panel A) ──
    HM, PA = [], []
    for sec, s_ in S.items():
        rows = []
        for lin, m in CANON:
            g1, g0 = s_['L'][lin], s_['keep'] & ~s_['L'][lin]
            v = s_['PZ'][m]
            d = fov_diffs(v, g1, g0, s_['fov'])
            rows.append(dict(record='top canonical pair', section=sec, group=lin, reference='All other cells', marker=m,
                             n_cells=int(g1.sum()), n_cells_ref=int(g0.sum()), mean_z=v[g1].mean(), mean_z_ref=v[g0].mean(),
                             **wtest(d)))
            for _, x in d.iterrows():
                PA.append(dict(record='top per FOV mean', section=sec, group=lin, reference='All other cells', marker=m,
                               fov=int(x.fov), mean_z=x.mean_group, mean_z_ref=x.mean_ref))
        r = pd.DataFrame(rows)
        ok = r.P.notna()
        r.loc[ok, 'q_BH'] = multipletests(r.loc[ok, 'P'], method='fdr_bh')[1]
        r['sig'] = [sig_of(q) if e > 0 else ('ns' if np.isfinite(q) else '') for q, e in zip(r.q_BH, r.estimate)]
        r['bh_family'] = f'{sec}: {len(CANON)} canonical lineage-marker pairs'
        HM.append(r)
    HM = pd.concat(HM, ignore_index=True)

    GATE_TXT = {'Epithelium': 'EpCAM > p75',
                'Fibroblast/Stroma': 'Fibronectin > p60 & Vimentin > p50 & SMA <= p85 & CD45, EpCAM, CD31 each <= p75',
                'Macrophage/Myeloid': ('myeloid parent [CD45 > p75 & CD3, CD20, CD19 <= p75 & (CD68|CD163|CD11c|CD14) > p75, '
                                       'minus neutrophil and pDC] & (CD68 > p75 | CD163 > p75) & CD11c <= p90'),
                'DC': 'CD11c > p75 & HLA-DR > p75 & CD68, CD14, CD3, CD20, CD19 <= p75 & not macrophage (Fig 5 relaxed gate)',
                'T Cells': 'CD3 > p75 (Treg = CD3, CD4, FOXP3, CTLA4 each > p75; Other T = CD3 > p75 & not Treg)',
                'B/Plasma Cells': '(CD20 > p75 | CD19 > p75) & CD3 <= p75 (protein B-cell gate; no CD138 gate in Fig 5)'}
    ZORDER = [('Epithelium', 'Epithelium'), ('B/Plasma Cells', 'B cell'), ('T Cells', 'Other T'),
              ('Fibroblast/Stroma', 'Stroma'), ('Macrophage/Myeloid', 'Macrophage'), ('DC', 'DC relaxed'), ('T Cells', 'Treg')]
    GATE_OF = {'Epithelium': ['Epithelium'], 'B/Plasma Cells': ['B cell'], 'T Cells': ['Other T', 'Treg'],
               'Fibroblast/Stroma': ['Stroma'], 'Macrophage/Myeloid': ['Macrophage'], 'DC': ['DC relaxed']}
    ref7 = pd.read_csv(REF7, low_memory=False)
    ref7 = ref7[ref7.record == 'n_cells_class']
    MAPROWS = []
    for sec, s in S.items():
        M = masks_of(s['ad'])
        n7 = {}
        for k in ['Epithelium', 'Stroma', 'Macrophage', 'DC relaxed', 'B cell', 'Treg', 'Other T']:
            n7[k] = int(ref7[(ref7.disease == sec) & (ref7.feature == k)].n_cells.iloc[0])
        prot = np.full(s['ad'].n_obs, 'Other', dtype=object)
        for lin, k in ZORDER:
            prot[M[k]] = lin
        s['prot'] = prot
        for lin in LIN6:
            MAPROWS.append(dict(record='maps section protein lineage counts', section=sec, group=lin,
                                n_cells=int(s['ad'].n_obs), n_FOV=int(len(np.unique(s['fov']))),
                                n_gate=int(sum(M[k].sum() for k in GATE_OF[lin])),
                                n_gate_v7=int(sum(n7[k] for k in GATE_OF[lin])),
                                n_protein=int((prot == lin).sum()),
                                gate=GATE_TXT[lin]))
        MAPROWS.append(dict(record='maps section protein lineage counts', section=sec, group='Other or Unassigned',
                            n_cells=int(s['ad'].n_obs), n_FOV=int(len(np.unique(s['fov']))),
                            n_protein=int((prot == 'Other').sum())))
    MR = pd.DataFrame(MAPROWS)

    # map geometry: global coordinates, every FOV, one um-per-inch for both maps, no rotation
    EXT = {sec: (S[sec]['xy'].min(0), S[sec]['xy'].max(0)) for sec in S}
    SPAN = {sec: EXT[sec][1] - EXT[sec][0] for sec in S}
    for sec in S:
        print(sec, 'section span (um) x, y =', SPAN[sec].round(0))
    MPAD = 25                                                   # um margin around each section
    H_UM = max(SPAN[sec][1] for sec in S) + 2 * MPAD            # both maps share this height (um)
    W_UM = {sec: SPAN[sec][0] + 2 * MPAD for sec in S}
    SB_UM = 500                                                 # scale bar (um)

    # ── bottom: PD-L1+ macrophage near PD-1+ Treg (read; median of per-FOV O/E recomputed from the CSV) ──
    lr = pd.read_csv(LRCSV, low_memory=False)
    lr = lr[(lr.pair == 'PD-L1 (Mac) -> PD-1 (Treg)') & (lr.definition == 'primary')].copy()
    lr['fov'] = lr['fov'].astype(str)
    LP = lr[lr.row == 'pooled'].set_index('radius_um')
    LF = lr[(lr.row != 'pooled') & lr.obs_exp.notna()].copy()          # the tested FOVs; O/E = 0 kept as 0
    MED = {}
    for r_ in [30, 50]:
        v = LF[LF.radius_um == r_].obs_exp.values.astype(float)
        rng_b = np.random.default_rng(BOOT_SEED)
        boot = np.array([np.median(v[rng_b.integers(0, len(v), len(v))]) for _ in range(N_BOOT)])
        MED[r_] = dict(median=float(np.median(v)), lo=float(np.percentile(boot, 2.5)), hi=float(np.percentile(boot, 97.5)),
                       n_fov=len(v), n_zero=int((v == 0).sum()))
    print('per-FOV median O/E', MED)
    BR = []
    for r_ in [30, 50]:
        x, m = LP.loc[r_], MED[r_]
        BR.append(dict(record='bottom LR proximity median of FOVs', section='CD',
                       group='PD-L1+ macrophage (RNA) -> PD-1+ Treg (RNA)', radius_um=r_,
                       median_fov_obs_exp=m['median'], median_ci95_lo=m['lo'], median_ci95_hi=m['hi'],
                       n_FOV=m['n_fov'], n_FOV_obs_exp_0=m['n_zero'], n_boot=N_BOOT, boot_seed=BOOT_SEED,
                       obs_exp_pooled=x.obs_exp, pooled_ci95_lo=x.ci95_lo, pooled_ci95_hi=x.ci95_hi,
                       perm_p=x.perm_p, sig=sig_of(x.perm_p), fovs_obs_exp_gt1=x.fovs_obs_exp_gt1, fovs_tested=x.fovs_tested,
                       obs_pairs=x.obs_pairs, exp_pairs=x.exp_pairs,
                       note=('median of the per-FOV O/E over the tested FOVs (O/E = 0 counted as 0), printed on the '
                             'spatial map with perm_p and stars; no O/E plot. median_ci95 = percentile 95% CI of the median, FOVs resampled with '
                             f'replacement ({N_BOOT} resamples, seed {BOOT_SEED}): NOT PLOTTED (CSV only; at 30 um the 3 FOVs '
                             'with O/E = 0 put the lower bound at 0). obs_exp_pooled / pooled_ci95 = pooled O/E with '
                             'FOV-bootstrap CI (10,000): NOT PLOTTED. Stars = perm_p (two-sided stratified permutation of '
                             'PD-L1+ status among macrophages within FOV, 1,000; one pre-specified test per radius, 30 um '
                             'primary, 50 um sensitivity); protein_lr_proximity.csv, definition primary: macrophages = RNA '
                             'SELENOP+ Mac or CD9+ M2 Mac, Tregs = RNA Treg (rna_qc_pass); PD-L1 / PD-1 protein > p75')))
    for x in LF.itertuples():
        BR.append(dict(record='bottom LR proximity per FOV', section='CD', group='PD-L1+ macrophage (RNA) -> PD-1+ Treg (RNA)',
                       radius_um=int(x.radius_um), fov=int(x.fov), obs_exp=x.obs_exp, perm_p=x.perm_p, obs_pairs=x.obs_pairs,
                       exp_pairs=x.exp_pairs, note='per-FOV O/E from protein_lr_proximity.csv (input to the median; not plotted)'))
    BR = pd.DataFrame(BR)

    # ── bottom map: PD-L1+ macrophage / PD-1+ Treg pairs (cells, gates and pair rule verbatim from the primary analysis of
    #    the S11D_protein_lr_proximity section) ──
    ob = C.obs
    qc = ob['rna_qc_pass'].values.astype(bool)
    fovC = ob['fov'].values.astype(int)
    ctgC = ob['cell_type_general'].astype(str).values
    xyC = S['CD']['xy']                                         # spatial_fov / px_per_um, as protein_lr_proximity.py
    above75 = lambda p: (lambda v: v > np.percentile(v, 75))(ob[f'{p}_geo'].values.astype(float))
    PDL1, PD1 = above75('PD-L1'), above75('PD-1')
    MAC_R = qc & np.isin(ctgC, ['SELENOP+ Mac', 'CD9+ M2 Mac'])
    TREG_R = qc & (ctgC == 'Treg')
    LIGM, RECT = MAC_R & PDL1, TREG_R & PD1                     # PD-L1+ macrophages, PD-1+ Tregs
    PAIRS = {30: [], 50: []}
    for f in FOVS_CD:
        inf = fovC == f
        mi, ri = np.flatnonzero(LIGM & inf), np.flatnonzero(RECT & inf)
        if len(mi) == 0 or len(ri) == 0:
            continue
        for r_ in PAIRS:
            nb = cKDTree(xyC[ri]).query_ball_point(xyC[mi], r=r_)
            PAIRS[r_] += [(f, mi[a], ri[b]) for a, lst in enumerate(nb) for b in lst]
    PAIRS = {r_: pd.DataFrame(v, columns=['fov', 'mac', 'treg']) for r_, v in PAIRS.items()}
    PMAC = np.unique(PAIRS[30].mac.values)                     # PD-L1+ macrophages in a 30 um pair
    PTREG = np.unique(PAIRS[30].treg.values)                   # PD-1+ Tregs in a 30 um pair
    UMAC = np.setdiff1d(np.flatnonzero(LIGM), PMAC)
    UTREG = np.setdiff1d(np.flatnonzero(RECT), PTREG)
    print('PD-L1+ mac', int(LIGM.sum()), 'PD-1+ Treg', int(RECT.sum()), '| 30 um pairs', len(PAIRS[30]), 'paired mac', len(PMAC),
          'paired Treg', len(PTREG), '| 50 um pairs', len(PAIRS[50]))
    STAT = [f"{r_} µm: Median O/E {MED[r_]['median']:.2f}, P = {LP.loc[r_].perm_p:.3g} {sig_of(LP.loc[r_].perm_p)}" for r_ in [30, 50]]
    PDROWS = [dict(record='bottom map PD-L1+ Mac / PD-1+ Treg', section='CD', group='PD-L1+ macrophage (RNA) -> PD-1+ Treg (RNA)',
                   radius_um=r_, n_FOV=len(FOVS_CD), obs_pairs=len(PAIRS[r_]), n_pdl1_mac=int(LIGM.sum()), n_pd1_treg=int(RECT.sum()),
                   n_pdl1_mac_paired=int(len(np.unique(PAIRS[r_].mac))), n_pd1_treg_paired=int(len(np.unique(PAIRS[r_].treg))),
                   median_fov_obs_exp=MED[r_]['median'], perm_p=LP.loc[r_].perm_p, sig=sig_of(LP.loc[r_].perm_p),
                   note=('observed pairs (obs_pairs) = for each PD-L1+ macrophage, the number of PD-1+ Tregs within r um in the same FOV, summed over PD-L1+ macrophages (protein_lr_proximity.py: W = receptor+ Tregs within r of each macrophage, obs = W[ligand+].sum()); many-to-many, so a macrophage near two Tregs counts twice and a Treg near two macrophages counts twice; cells in different FOVs are never paired. Cells and gates of the primary analysis of protein_lr_proximity.py; counts asserted equal to '
                         'its obs_pairs, per FOV and pooled. n_pdl1_mac_paired / n_pd1_treg_paired = cells with >= 1 partner of the '
                         'other type within radius ("near" on the map, 30 um only; 50 um counts are CSV only). Stat line printed on '
                         'the map = median_fov_obs_exp and perm_p.')) for r_ in [30, 50]]


    def pct_of(p):
        """Within-section percentile: 100 x fraction of all section cells (17 FOVs) with a value <= the cell's value."""
        v = ob[f'{p}_geo'].values.astype(float)
        return 100 * np.searchsorted(np.sort(v), v, side='right') / len(v)


    PCT = {'PD-L1': pct_of('PD-L1'), 'PD-1': pct_of('PD-1')}
    for x in PAIRS[30].itertuples():
        PDROWS.append(dict(record='bottom map pair (30 um)', section='CD', fov=int(x.fov), radius_um=30,
                           pair_distance_um=float(np.linalg.norm(xyC[x.mac] - xyC[x.treg])),
                           note=f'cell index mac {x.mac}, treg {x.treg}'))

    # ── Figure setup ──
    FT = 15
    plt.rcParams.update({'font.size': FT, 'axes.labelsize': FT, 'xtick.labelsize': 14, 'ytick.labelsize': 14,
                         'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
    LCOL = {'Epithelium': '#A3A3A3', 'Fibroblast/Stroma': '#3A86A8', 'Endothelium': '#222222',
            'Macrophage/Myeloid': '#F28E2B', 'DC': '#C2185B', 'T Cells': '#4DAF4A', 'B/Plasma Cells': '#7B4FA6'}
    C_OTHER = '#E3E3E3'
    TXTCOL = dict(LCOL, Epithelium='#6E6E6E', Endothelium='#222222')
    SECT = {'UC': 'UC Rectum', 'CD': 'CD Ileum'}

    DOT_X, DOT_W, ROW_H, ROW_GAP, TOP = 1.5, 6.8, 4.8, 0.3, 0.3
    UM_PER_IN = H_UM / ROW_H                       # same for both maps; map height = dot-plot row height
    MAP_W = {sec: W_UM[sec] / UM_PER_IN for sec in S}
    MAP_X = DOT_X + DOT_W + 0.35

    PAD = pd.DataFrame(PA)
    PAD['diff'] = PAD.mean_z - PAD.mean_z_ref
    rng = np.random.default_rng(0)                 # jitter drawn in the v4 order (UC pairs, then CD pairs)
    JIT = {}
    for sec in ['UC', 'CD']:
        for lin, m in CANON:
            n = int(((PAD.section == sec) & (PAD.group == lin) & (PAD.marker == m)).sum())
            JIT[(sec, lin, m)] = rng.uniform(-0.2, 0.2, n)
    ylims = {'UC': (-1.15, 2.95), 'CD': (-0.85, 0.95)}
    MAP_ORDER = ['Epithelium', 'B/Plasma Cells', 'T Cells', 'Fibroblast/Stroma', 'Macrophage/Myeloid', 'DC']
    MAP_S = {'Other': 0.7, 'Epithelium': 1.6, 'B/Plasma Cells': 2.4, 'T Cells': 2.4, 'Fibroblast/Stroma': 2.6,
             'Macrophage/Myeloid': 7, 'DC': 8}


    def draw_dot(ax, sec, xlabels):
        h = HM[HM.section == sec].reset_index(drop=True)
        nx = len(CANON)
        for i, (lin, m) in enumerate(CANON):
            d = PAD[(PAD.section == sec) & (PAD.group == lin) & (PAD.marker == m)]['diff'].values
            ax.scatter(i + JIT[(sec, lin, m)], d, s=75, color=LCOL[lin], alpha=0.85, edgecolor='white', lw=0.5, zorder=2)
            ax.plot([i - 0.34, i + 0.34], [h.estimate[i]] * 2, color='black', lw=3, zorder=3)
            if h.sig[i].startswith('*'):
                ax.text(i, 0.985, h.sig[i], transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=17,
                        fontweight='bold')
            if i and CANON[i - 1][0] != lin:
                ax.axvline(i - 0.5, color='#E2E2E2', lw=0.8, zorder=0)
        lo_, hi_ = ylims[sec]
        ax.set_ylim(lo_, hi_)
        ax.axhline(0, color='#888888', ls='--', lw=1, zorder=1)
        ax.set_xlim(-0.6, nx - 0.4)
        ax.set_xticks(range(nx))
        if xlabels:
            ax.set_xticklabels([f'{lin}: {m}' for lin, m in CANON], rotation=55, ha='right', rotation_mode='anchor')
            for t, (lin, _) in zip(ax.get_xticklabels(), CANON):
                t.set_color(TXTCOL[lin])
            ax.set_xlabel('RNA-Defined Lineage: Canonical Protein', labelpad=6)
        else:
            ax.set_xticklabels([])
        ax.set_ylabel(f'{SECT[sec]}\nProtein z-Score (Lineage − Other)')
        ax.spines[['top', 'right']].set_visible(False)


    def scale_bar(ax, x0, y0):
        """500 um bar in data units (um) of a map axis, label underneath."""
        ax.plot([x0, x0 + SB_UM], [y0, y0], color='black', lw=4, solid_capstyle='butt', clip_on=False)
        ax.text(x0 + SB_UM / 2, y0 + 0.06 * UM_PER_IN, f'{SB_UM} µm', ha='center', va='top', fontsize=14, clip_on=False)


    def draw_map(ax, sec, bar):
        s = S[sec]
        xy = s['xy'] - EXT[sec][0] + MPAD
        L_ = s['prot']
        ax.scatter(*xy[L_ == 'Other'].T, s=MAP_S['Other'], c=C_OTHER, lw=0, zorder=1, rasterized=True)
        for z, lin in enumerate(MAP_ORDER):
            mk = L_ == lin
            ax.scatter(*xy[mk].T, s=MAP_S[lin], c=LCOL[lin], lw=0, zorder=2 + z, rasterized=True)
        ax.set_xlim(0, W_UM[sec])
        ax.set_ylim(H_UM, 0)                        # image orientation, as Fig 5
        ax.set_aspect('equal')
        ax.axis('off')
        ax.text(0.02 if sec == 'CD' else 0.98, 0.995, SECT[sec], transform=ax.transAxes,
                ha='left' if sec == 'CD' else 'right', va='top', fontsize=15)
        if bar:
            scale_bar(ax, W_UM[sec] - SB_UM - MPAD, H_UM + 0.12 * UM_PER_IN)
        return {lin: int((L_ == lin).sum()) for lin in MAP_ORDER + ['Other']}


    KEY = ['Epithelium', 'Fibroblast/Stroma', 'Macrophage/Myeloid', 'DC', 'T Cells', 'B/Plasma Cells', 'Other']
    KLAB = {'Other': 'Other or Unassigned'}


    def draw_key(kax):
        kax.axis('off')
        hd = [Line2D([], [], marker='o', ls='', color=(C_OTHER if l == 'Other' else LCOL[l]), ms=12, label=KLAB.get(l, l))
              for l in KEY]
        kax.legend(handles=hd, ncol=1, loc='upper left', bbox_to_anchor=(0, 1.0), frameon=False, fontsize=14,
                   handletextpad=0.1, labelspacing=0.3, borderaxespad=0, borderpad=0, title='Protein-Defined Lineage',
                   title_fontsize=15, alignment='left')


    KEY_W, KEY_H = 2.6, 2.1                         # key box (inches): title + 7 rows


    def draw_keybar(fig, x, ytop, FW, FH):
        """Shared key and, under it, one scale bar at the map um-per-inch."""
        kax = fig.add_axes([x / FW, 1 - (ytop + KEY_H) / FH, KEY_W / FW, KEY_H / FH])
        draw_key(kax)
        sb_w = 1.2 * SB_UM / UM_PER_IN
        sax = fig.add_axes([(x + 0.05) / FW, 1 - (ytop + KEY_H + 0.2 + 0.4) / FH, sb_w / FW, 0.4 / FH])
        sax.set_xlim(0, 1.2 * SB_UM)
        sax.set_ylim(1, 0)
        sax.axis('off')
        sax.plot([0, SB_UM], [0.15, 0.15], color='black', lw=4, solid_capstyle='butt')
        sax.text(SB_UM / 2, 0.35, f'{SB_UM} µm', ha='center', va='top', fontsize=14)
        return kax, sax


    C_BG, S_BG, A_BG = '#D9D9D9', 0.5, 1.0                        # background as v10/v11 (v12 fainter one reverted, author)
    S_NEAR, LW_NEAR, S_FAR, A_FAR = 16, 0.8, 6, 0.85             # v12: near = larger, 0.8 pt black edge, on top; not near =
                                                                 # smaller, thin dark-grey edge (v13), 85% alpha
    CM_L = matplotlib.colors.ListedColormap(plt.get_cmap('OrRd')(np.linspace(0.12, 1, 256)))
    CM_R = matplotlib.colors.ListedColormap(plt.get_cmap('Blues')(np.linspace(0.3, 1, 256)))   # v14: rings, palest still a line
    S_RNEAR, LW_RNEAR, S_RFAR, LW_RFAR = 60, 1.5, 22, 0.8         # v14: PD-1+ Treg rings (no fill), above the macrophage dots
    NORM = matplotlib.colors.Normalize(75, 100)                  # gate-positive range uses the full colour range; linear in
    # percentile (gate-positive percentiles are spread evenly over 75-100: PD-L1 macrophages 33/36/33/26/34 and PD-1 Tregs
    # 9/8/21/13/16 per 5-point bin, so no rank transform is needed)
    C_EDGE, LW_EDGE = '#444444', 0.3                            # v13: thin dark-grey edge on every coloured dot
    PD_UPI, PD_H, PD_W = UM_PER_IN, ROW_H, MAP_W['CD']          # same um per inch and frame as the CD annotation map
    COL_W = 3.75                                                 # stat-line / key block width (in)
    CB_W = 2.7                                                   # colourbar width (in)


    def pd_key(ax):
        cm_, cr_ = CM_L(0.65), CM_R(0.65)
        dot = lambda ms, ec, mew, al=1: Line2D([], [], marker='o', ls='', mfc=cm_, mec=ec, mew=mew, ms=ms, alpha=al)
        ring = lambda ms, mew: Line2D([], [], marker='o', ls='', mfc='none', mec=cr_, mew=mew, ms=ms)
        hk = [dot(1.5 * np.sqrt(S_NEAR), C_EDGE, LW_EDGE * 1.3), ring(np.sqrt(S_RNEAR), LW_RNEAR),
              (dot(1.5 * np.sqrt(S_NEAR), 'black', LW_NEAR * 1.3), ring(np.sqrt(S_RNEAR), LW_RNEAR)),
              (dot(1.5 * np.sqrt(S_FAR), C_EDGE, LW_EDGE * 1.3, A_FAR), ring(np.sqrt(S_RFAR), LW_RFAR))]
        ax.legend(handles=hk, labels=['PD-L1⁺ Macrophage', 'PD-1⁺ Treg', 'Within 30 µm of a Partner Cell',
                                      'No Partner Within 30 µm'], ncol=1,
                  handler_map={tuple: HandlerTuple(ndivide=None, pad=0.5)}, loc='upper left', bbox_to_anchor=(0, 1),
                  frameon=False, fontsize=14, handletextpad=0.5, labelspacing=0.5, borderaxespad=0, borderpad=0,
                  handlelength=2.2)


    def pd_cells(a, D):
        """162 PD-L1+ macrophages = filled dots (OrRd by PD-L1 percentile): not near small, thin dark-grey edge, 85% alpha;
        near larger, 0.8 pt black edge. 67 PD-1+ Tregs = hollow rings (Blues by PD-1 percentile), drawn above all macrophage
        dots: not near small, 0.8 pt line; near larger, 1.5 pt line. Within each layer high values on top."""
        lm, lt = np.flatnonzero(LIGM), np.flatnonzero(RECT)
        near_m, near_t = np.isin(lm, PMAC), np.isin(lt, PTREG)
        srt = lambda ix, p_: ix[np.argsort(PCT[p_][ix], kind='stable')]
        for near in (False, True):
            ix = srt(lm[near_m == near], 'PD-L1')
            a.scatter(*D[ix].T, s=(S_NEAR if near else S_FAR), c=PCT['PD-L1'][ix], cmap=CM_L, norm=NORM,
                      ec=('black' if near else C_EDGE), lw=(LW_NEAR if near else LW_EDGE), alpha=(1 if near else A_FAR),
                      zorder=(3 if near else 2), rasterized=True)
        for near in (False, True):
            ix = srt(lt[near_t == near], 'PD-1')
            a.scatter(*D[ix].T, s=(S_RNEAR if near else S_RFAR), facecolors='none', edgecolors=CM_R(NORM(PCT['PD-1'][ix])),
                      lw=(LW_RNEAR if near else LW_RFAR), zorder=(5 if near else 4), rasterized=True)
        return dict(near_mac=int(near_m.sum()), far_mac=int((~near_m).sum()), near_treg=int(near_t.sum()),
                    far_treg=int((~near_t).sum()))


    def pd_cbars(fig, add, x, y):
        """Two compact colourbars (PD-L1 on macrophages, PD-1 on Tregs), percentile 75-100."""
        for j, (cm, lab) in enumerate([(CM_L, 'PD-L1, Macrophages (Percentile)'), (CM_R, 'PD-1, Tregs (Ring Colour, Percentile)')]):
            c = add(x, y + 0.32 + j * 0.85, CB_W, 0.16)
            cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=cm), cax=c, orientation='horizontal')
            cb.set_ticks([75, 80, 85, 90, 95, 100])
            cb.ax.tick_params(labelsize=12, length=3)
            cb.outline.set_linewidth(0.8)
            c.set_title(lab, fontsize=13, loc='left', pad=4)


    def draw_pd(fig, FW, FH, MXY, SXY, KXY, CXY):
        """Full-section PD-L1+ Mac / PD-1+ Treg map (top-left at MXY), stat lines (SXY), near / not-near key (KXY) and the
        two colourbars (CXY); inches from the figure's top-left."""
        add = lambda x, y, w, h: fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH])
        D = xyC - EXT['CD'][0] + MPAD
        other = ~(LIGM | RECT)
        a = add(MXY[0], MXY[1], PD_W, PD_H)
        a.scatter(*D[other].T, s=S_BG, c=C_BG, alpha=A_BG, lw=0, zorder=1, rasterized=True)
        nn = pd_cells(a, D)
        a.set_xlim(0, W_UM['CD'])
        a.set_ylim(H_UM, 0)
        a.set_aspect('equal')
        a.axis('off')
        a.text(0.02, 0.995, 'CD Ileum', transform=a.transAxes, ha='left', va='top', fontsize=15)
        sbx = MPAD
        a.plot([sbx, sbx + SB_UM], [H_UM + 0.12 * PD_UPI] * 2, color='black', lw=4, solid_capstyle='butt', clip_on=False)
        a.text(sbx + SB_UM / 2, H_UM + 0.18 * PD_UPI, f'{SB_UM} µm', ha='center', va='top', fontsize=14, clip_on=False)
        t = add(SXY[0], SXY[1], COL_W, 0.6)
        t.axis('off')
        t.text(0, 1, '\n'.join(STAT), transform=t.transAxes, ha='left', va='top', fontsize=14, linespacing=1.5)
        k = add(KXY[0], KXY[1], COL_W, 1.4)
        k.axis('off')
        pd_key(k)
        pd_cbars(fig, add, *CXY)
        n_draw = dict(pmac=nn['near_mac'], ptreg=nn['near_treg'], umac=nn['far_mac'], utreg=nn['far_treg'],
                      other=int(other.sum()))
        return a, n_draw


    # ── Combined figure ──
    MAPS_W = max(MAP_W.values())
    KEY_X = MAP_X + MAPS_W + 0.25
    FW = KEY_X + KEY_W + 0.1
    MAP_BOT = TOP + 2 * ROW_H + ROW_GAP
    PD_X = KEY_X                                                 # PD-L1 / PD-1 map: right column, beside the CD annotation map
    FW = max(FW, PD_X + PD_W + 0.05)                            # v6b: PD stat lines + key under the CD annotation map
    FH = MAP_BOT + 3.1                                           # v11: stat lines + key under the CD annotation map; colourbars under the PD map
    fig = plt.figure(figsize=(FW, FH))


    def put(x, ytop, w, h):
        return fig.add_axes([x / FW, 1 - (ytop + h) / FH, w / FW, h / FH])


    AX = {}
    for j, sec in enumerate(['UC', 'CD']):
        AX[f'dot_{sec}'] = put(DOT_X, TOP + j * (ROW_H + ROW_GAP), DOT_W, ROW_H)
        draw_dot(AX[f'dot_{sec}'], sec, xlabels=(sec == 'CD'))
    MAPN = []
    for j, sec in enumerate(['UC', 'CD']):
        AX[f'map_{sec}'] = put(MAP_X, TOP + j * (ROW_H + ROW_GAP), MAP_W[sec], ROW_H)
        n_ = draw_map(AX[f'map_{sec}'], sec, bar=False)
        MAPN += [(sec, lin, n) for lin, n in n_.items()]
    MAPN = pd.DataFrame(MAPN, columns=['section', 'lineage', 'n_cells_on_map'])
    print(MAPN.to_string())
    KEY_TOP = TOP
    AX['key'], AX['bar'] = draw_keybar(fig, KEY_X, KEY_TOP, FW, FH)
    PD_TOP = TOP + ROW_H + ROW_GAP                               # aligned with the CD row
    AX['pd'], PD_N = draw_pd(fig, FW, FH, (PD_X, PD_TOP), (MAP_X, MAP_BOT + 0.62), (MAP_X, MAP_BOT + 1.37),
                             (MAP_X + COL_W + 0.1, MAP_BOT + 0.6))          # colourbars right of the stat lines
    print('PD map drawn', PD_N)

    # ── CSV ──
    LAYOUT = [dict(record='maps layout', section=sec, n_cells=int(S[sec]['ad'].n_obs),
                   n_FOV=int(len(np.unique(S[sec]['fov']))), width_um=SPAN[sec][0], height_um=SPAN[sec][1],
                   um_per_inch=UM_PER_IN, map_width_in=MAP_W[sec], map_height_in=ROW_H, scale_bar_um=SB_UM, rotation_deg=0,
                   note=(f"full section, every FOV ({'UC rectum 385R, all 15 FOVs' if sec == 'UC' else 'CD ileum 375i, the 17-FOV section of Fig 5 (FOVs ' + ', '.join(map(str, FOVS_CD)) + ')'}), "
                         "global coordinates (obsm spatial_fov / px_per_um, as fig5_combined_v1.py); all cells drawn, one "
                         'protein label per cell; same um-per-inch for both maps; one shared 500 um scale bar under the '
                         'key; no rotation (both sections are taller than wide, and the maps fill the dot-plot row height). '
                         'Section sizes differ in width (UC narrower than CD); scale kept equal, so the UC map is narrower. '
                         'Protein label = Fig 5 gates (fig5_combined_v1.py compute_v7; per-gate section counts asserted '
                         'equal to fig5_spatial_maps_v7.csv n_cells_class), one label per cell by Fig 5 map z-order (later '
                         'wins): Epithelium < B cell < Other T < Stroma < Macrophage < DC (relaxed) < Treg; T Cells = Treg + '
                         'Other T. n_gate = cells passing the gate (= n_gate_v7); n_protein = cells drawn in that colour '
                         'after z-order.')) for sec in ['UC', 'CD']]
    PDROWS.append(dict(record='bottom map layout', section='CD', radius_um=30, n_FOV=len(FOVS_CD),
                       um_per_inch=PD_UPI, map_width_in=PD_W, map_height_in=PD_H, scale_bar_um=SB_UM, rotation_deg=0,
                       note=(f'full section, no rotation, {PD_UPI:.0f} um per inch (same as the annotation maps), 500 um scale '
                             'bar; no zoom inset (v8). v14 markers: only the 162 PD-L1+ macrophages and 67 PD-1+ Tregs; colour '
                             '= within-section protein percentile (100 x fraction of all section cells <= value; > 75 = the gate, '
                             'asserted), linear over 75-100 (gate-positive percentiles are spread evenly). Macrophages = filled dots, '
                             'OrRd 0.12-1.0; PD-1+ Tregs = hollow rings (no fill), Blues 0.3-1.0, drawn above the macrophage dots. '
                             'Near = >= 1 partner of the other type within 30 um in the same FOV (24 macrophages, 22 Tregs): '
                             f'macrophage dot {S_NEAR} pt^2 with a {LW_NEAR} pt black edge, Treg ring {S_RNEAR} pt^2 with a '
                             f'{LW_RNEAR} pt line; not near (138 macrophages, 45 Tregs): dot {S_FAR} pt^2, {LW_EDGE} pt dark-grey '
                             f'edge, {A_FAR:.0%} alpha; ring {S_RFAR} pt^2, {LW_RFAR} pt line. Other cells light grey (#D9D9D9, '
                             '0.5 pt^2, as v10); no connecting lines (pairs are many-to-many). '
                             f'observed pairs (obs_pairs) = for each PD-L1+ macrophage, the number of PD-1+ Tregs within r um in the same FOV, summed over PD-L1+ macrophages (protein_lr_proximity.py: W = receptor+ Tregs within r of each macrophage, obs = W[ligand+].sum()); many-to-many, so a macrophage near two Tregs counts twice and a Treg near two macrophages counts twice; cells in different FOVs are never paired.')))
    R = pd.concat([HM.rename(columns={'estimate': 'median_FOV_diff_z'}), pd.DataFrame(PA), MR, pd.DataFrame(LAYOUT), BR, pd.DataFrame(PDROWS)],
                  ignore_index=True)
    cols = ['record', 'section', 'group', 'reference', 'marker', 'fov', 'n_cells', 'n_cells_ref', 'n_FOV', 'n_FOV_higher',
            'mean_z', 'mean_z_ref', 'median_FOV_diff_z', 'P', 'q_BH', 'sig', 'bh_family',
            'n_gate', 'n_gate_v7', 'n_protein', 'gate',
            'width_um', 'height_um', 'um_per_inch', 'map_width_in', 'map_height_in', 'scale_bar_um', 'rotation_deg',
            'radius_um', 'median_fov_obs_exp', 'median_ci95_lo', 'median_ci95_hi', 'n_FOV_obs_exp_0', 'n_boot', 'boot_seed',
            'obs_exp', 'obs_exp_pooled', 'pooled_ci95_lo', 'pooled_ci95_hi', 'perm_p', 'fovs_obs_exp_gt1', 'fovs_tested',
            'obs_pairs', 'exp_pairs', 'n_pdl1_mac', 'n_pd1_treg', 'n_pdl1_mac_paired', 'n_pd1_treg_paired', 'pair_distance_um',
            'note']
    R = R[cols]
    for c in R.columns:
        if (c.startswith('n_') or c in ('fov', 'radius_um', 'fovs_obs_exp_gt1', 'fovs_tested', 'obs_pairs', 'scale_bar_um',
                                        'rotation_deg', 'boot_seed')):
            R[c] = pd.to_numeric(R[c]).round().astype('Int64')
        elif R[c].dtype.kind == 'f':
            R[c] = R[c].map(lambda x: x if pd.isna(x) else float(f'{x:.3g}'))
    R['note'] = R['note'].fillna('')
    R.loc[0, 'note'] = ('UC rectum 385R (biopsy, all cells); CD ileum 375i (resection, 17-FOV subset, rna_qc_pass); protein '
                        "'_geo' z over all cells of the section; unit = FOV; two-sided Wilcoxon over FOV differences; BH over "
                        'the 12 pairs per section; star only if q < 0.05 and median > 0; no UC vs CD test')
    with open(f'{STEM}.csv', 'w') as fh:
        R.to_csv(fh, index=False)
    pd.set_option('display.width', 250, 'display.max_rows', 300)
    print(MR.drop(columns='gate').to_string())

    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.08)
    print('wrote', f'{STEM}.csv', f'{STEM}.png', 'size (in)', FW, FH, 'um per inch', round(UM_PER_IN, 1))

    # ── Panel files (same drawing functions, same sizes, same um-per-inch) ──
    os.makedirs(PANELS, exist_ok=True)


    def save_panel(name, w, h, draw, pad_bottom=0.0):
        f = plt.figure(figsize=(w + 3, h + 3 + pad_bottom))
        a = f.add_axes([1.5 / (w + 3), (1.5 + pad_bottom) / (h + 3 + pad_bottom), w / (w + 3), h / (h + 3 + pad_bottom)])
        draw(a)
        for ext in ('png',):
            f.savefig(f'{PANELS}/{name}.{ext}', dpi=300, bbox_inches='tight', pad_inches=0.08)
        plt.close(f)
        return name


    PF = []
    PF.append(save_panel('uc_rectum_rna_lineage_protein_dotplot', DOT_W, ROW_H, lambda a: draw_dot(a, 'UC', True)))
    PF.append(save_panel('cd_ileum_rna_lineage_protein_dotplot', DOT_W, ROW_H, lambda a: draw_dot(a, 'CD', True)))
    PF.append(save_panel('uc_rectum_protein_annotation_map', MAP_W['UC'], ROW_H, lambda a: draw_map(a, 'UC', True)))
    PF.append(save_panel('cd_ileum_protein_annotation_map', MAP_W['CD'], ROW_H, lambda a: draw_map(a, 'CD', True)))
    PW, PH = COL_W + 1, PD_H + 5.3                               # stat lines, key and colourbars stacked under the map
    fp = plt.figure(figsize=(PW, PH))
    draw_pd(fp, PW, PH, (0.5, 0.5), (0.5, 0.5 + PD_H + 0.6), (0.5, 0.5 + PD_H + 1.3), (0.5, 0.5 + PD_H + 2.5))
    for ext in ('png',):
        fp.savefig(f'{PANELS}/cd_ileum_pdl1_mac_pd1_treg_spatial_map.{ext}', dpi=300, bbox_inches='tight', pad_inches=0.08)
    plt.close(fp)
    PF.append('cd_ileum_pdl1_mac_pd1_treg_spatial_map')
    fk = plt.figure(figsize=(KEY_W + 1, KEY_H + 2))
    draw_keybar(fk, 0.5, 0.5, KEY_W + 1, KEY_H + 2)
    for ext in ('png',):
        fk.savefig(f'{PANELS}/protein_lineage_key_scale_bar.{ext}', dpi=300, bbox_inches='tight', pad_inches=0.08)
    plt.close(fk)
    PF.append('protein_lineage_key_scale_bar')
    print('panels:', PF)
