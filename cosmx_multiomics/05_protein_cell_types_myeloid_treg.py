#!/usr/bin/env python
"""CosMx spatial proteomics, UC vs CD: protein-gated cell types, myeloid activation and Treg receptors (Fig 7A-E).

Gates cell types on protein (percentile thresholds within each section), maps them, profiles composition by distance
from the epithelium, compares inflammatory and antigen-presentation proteins in macrophages and DCs (Cohen's d, CD - UC),
and compares ICOS, GITR and CD27 on Tregs and conventional CD4+ T cells. The gates and marker lists form Supp Table 12.

Inputs   UC object of 03_uc_rna_protein_processing.ipynb, CD object of 04_cd_rna_protein_processing.ipynb
Outputs  fig5_combined_v1.png (Fig 7A-E), fig7_cell_types_depth_profiles.csv (Fig 7A-B), fig7_myeloid_treg_markers.csv
         (Fig 7C-E)
Run      OPENBLAS_CORETYPE=SkylakeX python 05_protein_cell_types_myeloid_treg.py [png_path]
"""
import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
OUT = '/share/fsmresfiles/UC/revision_jenny/Fig5_spatial_protein'
PNG = sys.argv[1] if len(sys.argv) > 1 else f'{OUT}/fig5_combined_v1.png'
CSV_CT = f'{OUT}/fig7_cell_types_depth_profiles.csv'
CSV_MK = f'{OUT}/fig7_myeloid_treg_markers.csv'
plt.rcParams.update({'font.size': 15, 'axes.titlesize': 18, 'axes.labelsize': 17,
                     'xtick.labelsize': 16, 'ytick.labelsize': 16, 'legend.fontsize': 15})

# ---- load: UC rectum 385R (15 FOVs), CD ileum 375i (17-FOV subset) ----
CD_PATH = '/share/fsmresfiles/UC/AtoMx/CD_multiOmics/Cos4_315i_375i'
UC_PATH = '/share/fsmresfiles/UC/scOmics_model/data/Spatial_Protein_CosMx'
adata_rna_qc = sc.read_h5ad(f'{CD_PATH}/rna_norm_pr_geonorm_upd_anno.h5ad')
adata_rna_qc_375 = adata_rna_qc[(adata_rna_qc.obs['fov'] <= 29)].copy()
fovs = [6, 7, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29]
adata_rna_qc_375_sub = adata_rna_qc_375[adata_rna_qc_375.obs['fov'].isin(fovs)].copy()
adata_rna_qc_385 = sc.read_h5ad(f'{UC_PATH}/RNA1_Cos4_385R_381R_S2_norm_wpr_anno_jd.h5ad')
DATASETS = [('UC', adata_rna_qc_385), ('CD', adata_rna_qc_375_sub)]
px_per_um = np.sqrt(3266 / 47.25)
MIN_DC = 50   # strict DC gate needs >= 50 cells in each section, otherwise the relaxed gate is used
GATE_MARKERS = ['CD45', 'CD3', 'CD20', 'CD19', 'CD68', 'CD163', 'CD11c', 'CD14', 'CD15', 'CD16', 'CD123']


# ---- shared helpers: protein gates and statistics ----
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


def cell_masks(adata):
    """All protein-gated classes of one section (thresholds within section)."""
    lab, pos = gate_myeloid(adata)
    g = lambda p: adata.obs[f'{p}_geo'].values.astype(float)
    above = lambda p, q: g(p) > np.percentile(g(p), q)
    notabove = lambda p, q: g(p) <= np.percentile(g(p), q)
    mac = lab == 'Macrophage'
    treg = pos['CD3'] & above('CD4', 75) & above('FOXP3', 75) & above('CTLA4', 75)
    return {'Epithelium': above('EpCAM', 75),
            'Stroma': (above('Fibronectin', 60) & above('Vimentin', 50) & notabove('SMA', 85)
                       & notabove('CD45', 75) & notabove('EpCAM', 75) & notabove('CD31', 75)),
            'Macrophage': mac, 'DC strict': lab == 'DC',
            'DC relaxed': (above('CD11c', 75) & above('HLA-DR', 75) & notabove('CD68', 75) & notabove('CD14', 75)
                           & ~pos['CD3'] & ~pos['CD20'] & ~pos['CD19'] & ~mac),
            'B cell': (pos['CD20'] | pos['CD19']) & ~pos['CD3'],
            'Treg': treg, 'Other T': pos['CD3'] & ~treg, 'Conv CD4 T': pos['CD3'] & above('CD4', 75) & ~treg,
            '_hla': above('HLA-DR', 75), '_inos': above('iNOS', 75)}


def sig_of(q):
    return 'ns' if not np.isfinite(q) else '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'


def se(x):
    x = np.asarray(x, float)
    return np.std(x, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    return mannwhitneyu(a, b, alternative='two-sided')[1] if min(len(a), len(b)) >= 3 else np.nan


def bh(df):
    """BH within each 'family'."""
    df['q_BH'] = np.nan
    for fam, ix in df.groupby('family').groups.items():
        p = df.loc[ix, 'P']
        ok = p.notna()
        if ok.any():
            df.loc[p[ok].index, 'q_BH'] = multipletests(p[ok], method='fdr_bh')[1]
    df['sig'] = df['q_BH'].map(sig_of)
    return df


def cohen_d(c, u):
    n1, n2 = len(c), len(u)
    sp = np.sqrt(((n1 - 1) * c.var(ddof=1) + (n2 - 1) * u.var(ddof=1)) / (n1 + n2 - 2))
    return (c.mean() - u.mean()) / sp if sp > 0 else np.nan


def round_table(out, int_cols=()):
    for c in out.columns:
        if c.startswith('n_') or c in int_cols:
            out[c] = out[c].astype('Int64')
        elif out[c].dtype.kind == 'f':
            out[c] = out[c].map(lambda x: x if pd.isna(x) else float(f'{x:.3g}'))
    return out


MASKS = {cond: cell_masks(adata) for cond, adata in DATASETS}
N_STRICT = {c: int(MASKS[c]['DC strict'].sum()) for c in MASKS}
DC_GATE = 'DC strict' if min(N_STRICT.values()) >= MIN_DC else 'DC relaxed'
XY = {cond: np.asarray(adata.obsm['spatial_fov'])[:, :2].astype(float) / px_per_um for cond, adata in DATASETS}


# ---- computation: cell-type maps and depth profiles (Fig 7A-B) ----
def compute_cell_type_profiles():
    MIN_BIN, MIN_NFOV = 20, 3
    BINS = [0, 25, 50, 100, 200, np.inf]
    BIN_LAB = ['0–25', '25–50', '50–100', '100–200', '>200']
    FEATS = ['Stroma', 'APC', 'DC', 'Macrophage', 'Treg', 'Other T']
    PLOT_FEATS = ['Stroma', 'APC', 'Treg', 'Other T']
    FNAME = {'Stroma': 'Fibroblast-rich stroma', 'APC': 'All APCs (DC + macrophage + B)',
             'DC': 'DCs', 'Macrophage': 'Macrophages', 'Treg': 'Tregs (FOXP3⁺CTLA4⁺)',
             'Other T': 'Other T cells'}
    THR_MARKERS = {'EpCAM': [75], 'Fibronectin': [60], 'Vimentin': [50], 'SMA': [85], 'CD31': [75],
                   'CD45': [75], 'CD3': [75], 'CD20': [75], 'CD19': [75], 'CD68': [75], 'CD163': [75],
                   'CD11c': [75, 90], 'CD14': [75], 'CD15': [75], 'CD16': [75], 'CD123': [90],
                   'HLA-DR': [75], 'CD4': [75], 'FOXP3': [75], 'CTLA4': [75]}
    GATE_TXT = {
        'Epithelium': 'EpCAM > p75',
        'Stroma': 'Fibronectin > p60 & Vimentin > p50 & SMA <= p85 & CD45 <= p75 & EpCAM <= p75 & CD31 <= p75',
        'Macrophage': ('myeloid parent [CD45 > p75 & CD3,CD20,CD19 <= p75 & (CD68|CD163|CD11c|CD14) > p75, '
                       'minus neutrophil CD15&CD16 > p75 and pDC CD123 > p90 & CD11c <= p75] & '
                       '(CD68 > p75 | CD163 > p75) & CD11c <= p90'),
        'DC strict': 'myeloid parent & CD11c > p90 & CD68 <= p75 & CD14 <= p75 & CD163 <= p75',
        'DC relaxed': ('CD11c > p75 & HLA-DR > p75 & CD68 <= p75 & CD14 <= p75 & CD3,CD20,CD19 <= p75 '
                       '& not macrophage'),
        'B cell': '(CD20 > p75 | CD19 > p75) & CD3 <= p75',
        'Treg': 'CD3 > p75 & CD4 > p75 & FOXP3 > p75 & CTLA4 > p75',
        'Other T': 'CD3 > p75 & not Treg',
        'APC': 'DC (gate used) | macrophage | B cell'}
    KEYS = ['Epithelium', 'Stroma', 'Macrophage', 'DC strict', 'DC relaxed', 'B cell', 'Treg', 'Other T']

    def dist_to_epi(xy, epi):
        d, _ = cKDTree(xy[epi]).query(xy, k=2)
        out = d[:, 0].copy()
        out[epi] = d[epi, 1]
        return out / px_per_um

    D, rows = {}, []
    for cond, adata in DATASETS:
        xy = XY[cond]
        M = {k: MASKS[cond][k] for k in KEYS}
        g = lambda p: adata.obs[f'{p}_geo'].values.astype(float)
        thr = {f'{m} p{q}': np.percentile(g(m), q) for m, qs in THR_MARKERS.items() for q in qs}
        d = dist_to_epi(xy * px_per_um, M['Epithelium'])
        b = np.asarray(pd.cut(d, BINS, right=False, labels=BIN_LAB).astype(str))
        fov = adata.obs['fov'].values
        D[cond] = dict(xy=xy, d=d, bin=b, M=M, fov=fov, n=adata.n_obs)
        rows.append(dict(record='section', disease=cond, n_cells=adata.n_obs,
                         n_FOV=len(np.unique(fov)), dist_median_um=np.median(d),
                         dist_p99_um=np.percentile(d, 99), dist_max_um=d.max(),
                         pct_cells_gt100um=100 * np.mean(d > 100), pct_cells_gt200um=100 * np.mean(d > 200)))
        for k, val in thr.items():
            rows.append(dict(record='threshold (_geo, within section)', disease=cond, feature=k, value=val))
    for c in D:
        M = D[c]['M']
        M['DC'] = M[DC_GATE]
        M['APC'] = M['DC'] | M['Macrophage'] | M['B cell']
        for k in KEYS + ['APC']:
            gt = GATE_TXT[k]
            if k == 'APC':
                gt = gt.replace('gate used', DC_GATE)
            rows.append(dict(record='n_cells_class', disease=c, feature=k, gate=gt,
                             n_cells=int(M[k].sum()), pct_of_section=100 * M[k].mean(),
                             used_in_figure=(k != ('DC relaxed' if DC_GATE == 'DC strict' else 'DC strict'))))
        rows.append(dict(record='n_cells_overlap', disease=c, feature='DC strict & DC relaxed',
                         n_cells=int((M['DC strict'] & M['DC relaxed']).sum())))
        rows.append(dict(record='n_cells_overlap', disease=c, feature='DC used & Epithelium',
                         n_cells=int((M['DC'] & M['Epithelium']).sum())))
    # % of cells per distance bin per FOV
    prof = []
    for cond in ['UC', 'CD']:
        X = D[cond]
        for f in np.unique(X['fov']):
            for bl in BIN_LAB:
                ib = (X['fov'] == f) & (X['bin'] == bl)
                n = int(ib.sum())
                r = dict(disease=cond, fov=f, bin=bl, n=n)
                for k in FEATS:
                    r[k] = 100 * X['M'][k][ib].mean() if n >= MIN_BIN else np.nan
                prof.append(r)
    prof = pd.DataFrame(prof)
    prof_rows = []
    for cond in ['UC', 'CD']:
        for k in FEATS:
            for bl in BIN_LAB:
                p = prof[(prof.disease == cond) & (prof.bin == bl)]
                v = p[k].dropna().values
                prof_rows.append(dict(record='depth_profile (% of cells in bin; FOV mean +/- SE)',
                                      disease=cond, feature=FNAME[k], depth_bin_um=bl, n_FOV=len(v),
                                      n_cells=int(p.loc[p[k].notna(), 'n'].sum()),
                                      fov_mean=(v.mean() if len(v) else np.nan),
                                      fov_SE=(se(v) if len(v) else np.nan),
                                      fov_median=(np.median(v) if len(v) else np.nan),
                                      shown=len(v) >= MIN_NFOV))
    PR = pd.DataFrame(prof_rows)
    # FOV-paired Wilcoxon zone tests (BH over 12 tests; CSV only)
    ZONES = {'CD': [('0–50 µm', 0, 50), ('> 100 µm', 100, np.inf)],
             'UC': [('0–50 µm', 0, 50), ('50–100 µm', 50, 100)]}
    tst = []
    for cond in ['CD', 'UC']:
        X = D[cond]
        (za, a0, a1), (zb, b0, b1) = ZONES[cond]
        for k in FEATS:
            A, B, F = [], [], []
            for f in np.unique(X['fov']):
                inf = X['fov'] == f
                sa = inf & (X['d'] >= a0) & (X['d'] < a1)
                sb = inf & (X['d'] >= b0) & (X['d'] < b1)
                if min(sa.sum(), sb.sum()) >= MIN_BIN:
                    A.append(100 * X['M'][k][sa].mean()); B.append(100 * X['M'][k][sb].mean()); F.append(f)
            A, B = np.array(A), np.array(B)
            gt = GATE_TXT[DC_GATE if k == 'DC' else k]
            if k == 'APC':
                gt = gt.replace('gate used', DC_GATE)
            tst.append(dict(record='zone_test (FOV-paired Wilcoxon, BH over 12 tests; CSV only)', disease=cond,
                            feature=FNAME[k], gate=gt, zone_superficial=za, zone_deeper=zb,
                            n_FOV=len(A), fovs=';'.join(str(int(x)) for x in F),
                            fov_mean_superficial=A.mean(), fov_mean_deeper=B.mean(),
                            fov_median_superficial=np.median(A), fov_median_deeper=np.median(B),
                            mean_diff_deeper_minus_superficial=(B - A).mean(),
                            n_FOV_higher_deeper=int((B > A).sum()),
                            P=(wilcoxon(B, A)[1] if len(A) >= 3 and np.any(B != A) else np.nan)))
    tst = pd.DataFrame(tst)
    ok = tst['P'].notna()
    tst['q_BH'] = np.nan
    tst.loc[ok, 'q_BH'] = multipletests(tst.loc[ok, 'P'], method='fdr_bh')[1]
    tst['sig'] = tst['q_BH'].map(lambda q: 'n/a' if pd.isna(q) else sig_of(q))
    PR['in_figure'] = PR['feature'].isin([FNAME[k] for k in PLOT_FEATS])
    out = pd.concat([tst, PR, pd.DataFrame(rows)], ignore_index=True)
    lead = ['record', 'disease', 'feature', 'depth_bin_um', 'zone_superficial', 'zone_deeper']
    out = round_table(out[[c for c in lead if c in out] + [c for c in out if c not in lead]])
    out['notes'] = ''
    out.loc[0, 'notes'] = ('protein only; thresholds within section over all cells; distance = nearest '
                           'epithelial cell (EpCAM > p75; other than itself), um; unit = FOV; zone/bin '
                           f'needs >= {MIN_BIN} cells per FOV; profile point shown if >= {MIN_NFOV} FOVs; '
                           f'DC gate used = {DC_GATE} (strict needs >= {MIN_DC} cells in each section; '
                           f'strict n UC {N_STRICT["UC"]}, CD {N_STRICT["CD"]}); '
                           'UC = rectum 385R (15 FOVs), CD = ileum 375i (17-FOV subset); one mucosal biopsy '
                           'section per disease, so disease and location are confounded; descriptive; '
                           'zone tests are reported here only (not plotted); depth_profile rows with '
                           'in_figure = True are plotted in panel B')
    out.to_csv(CSV_CT, index=False)
    return dict(D=D, PR=PR, FNAME=FNAME, PLOT_FEATS=PLOT_FEATS, BIN_LAB=BIN_LAB)


# ---- computation: myeloid markers and Treg receptors (Fig 7C-E) ----
def compute_marker_panels():
    MIN_CELLS, MIN_TREG, MIN_CONV = 3, 3, 3
    TREG_MARKERS = ['ICOS', 'GITR', 'CD27']
    FOREST = {'Inflammatory': ['iNOS', 'IL-1b', 'IL-18'],
              'Antigen presentation / co-stimulation': ['HLA-DR', 'B7-H3', 'IDO1']}
    FOREST_SUB = ['Macrophage', 'DC']
    DAGGER = {'IL-1b'}
    N_BOOT = 5000
    GATE_TXT = {
        'Epithelium': 'EpCAM > p75',
        'Stroma': 'Fibronectin > p60 & Vimentin > p50 & SMA <= p85 & CD45, EpCAM, CD31 each <= p75',
        'Macrophage': 'myeloid parent & (CD68 > p75 | CD163 > p75) & CD11c <= p90',
        'DC strict': 'myeloid parent & CD11c > p90 & CD68 <= p75 & CD14 <= p75 & CD163 <= p75',
        'DC relaxed': ('CD11c > p75 & HLA-DR > p75 & CD68 <= p75 & CD14 <= p75 & CD3,CD20,CD19 <= p75 '
                       '& not macrophage'),
        'Myeloid': 'macrophage | DC',
        'HLA-DR-high myeloid': 'myeloid & HLA-DR > p75',
        'iNOS-high myeloid': 'myeloid & iNOS > p75',
        'Treg': 'CD3, CD4, FOXP3, CTLA4 each > p75',
        'Conv CD4 T': 'CD3 > p75 & CD4 > p75 & not Treg'}
    GATE_TXT['DC'] = GATE_TXT[DC_GATE] + ' (gate used, as in the cell-type maps)'
    KEYS = ['Epithelium', 'Stroma', 'Macrophage', 'DC strict', 'DC relaxed', 'Treg', 'Conv CD4 T']

    D, ncell = {}, []
    ALLP = sorted(set(TREG_MARKERS) | set(sum(FOREST.values(), [])))
    for cond, adata in DATASETS:
        M = {k: MASKS[cond][k] for k in KEYS}
        M['DC'] = M[DC_GATE]
        M['Myeloid'] = M['Macrophage'] | M['DC']
        M['HLA-DR-high myeloid'] = M['Myeloid'] & MASKS[cond]['_hla']
        M['iNOS-high myeloid'] = M['Myeloid'] & MASKS[cond]['_inos']
        fov = adata.obs['fov'].values.astype(int)
        V = {p: adata.obs[f'{p}_geo'].values.astype(float) for p in ALLP}
        D[cond] = dict(xy=XY[cond], M=M, fov=fov, V=V)
        for k, m in M.items():
            nf = sum((m & (fov == f)).sum() >= MIN_CELLS for f in np.unique(fov))
            ncell.append(dict(record='n_cells', disease=cond, subset=k, gate=GATE_TXT.get(k, ''), n_cells=int(m.sum()),
                              pct_of_section=100 * m.mean(), n_FOV_ge3=nf))
        ncell.append(dict(record='n_cells', disease=cond, subset='all cells', n_cells=adata.n_obs,
                          n_FOV=len(np.unique(fov))))
        ncell.append(dict(record='n_cells_overlap', disease=cond, subset='DC & Epithelium',
                          n_cells=int((M['DC'] & M['Epithelium']).sum())))
    # forest: per-FOV means, Cohen's d (CD - UC) with bootstrap CI, Mann-Whitney
    rng_b = np.random.default_rng(20261109)
    fr, forest_fov, perfov = [], [], []
    for s in FOREST_SUB:
        for fam, mk in FOREST.items():
            for p in mk:
                vals = {}
                for cond in ['UC', 'CD']:
                    X = D[cond]
                    qf = [f for f in np.unique(X['fov']) if (X['M'][s] & (X['fov'] == f)).sum() >= MIN_CELLS]
                    vals[cond] = np.array([X['V'][p][X['M'][s] & (X['fov'] == f)].mean() for f in qf])
                    for f, v in zip(qf, vals[cond]):
                        forest_fov.append(dict(record='forest per FOV mean', subset=s, protein=p, disease=cond,
                                               location='tissue', fov=f, value=v,
                                               n_cells=int((X['M'][s] & (X['fov'] == f)).sum())))
                u, c = vals['UC'], vals['CD']
                gate_mk = (s == 'DC' and p == 'HLA-DR' and DC_GATE == 'DC relaxed')
                d = cohen_d(c, u)
                bs = np.array([cohen_d(rng_b.choice(c, len(c)), rng_b.choice(u, len(u))) for _ in range(N_BOOT)])
                bs = bs[np.isfinite(bs)]
                fr.append(dict(record='forest tissue-wide UC vs CD (Cohen d, Mann-Whitney)', family=f'Forest | {s}',
                               subset=s, protein=p, function_group=fam, location='tissue', n_FOV_UC=len(u),
                               n_FOV_CD=len(c), mean_UC=u.mean(), mean_CD=c.mean(), median_UC=np.median(u),
                               median_CD=np.median(c), cohen_d_CD_minus_UC=d,
                               d_CI95_low=np.percentile(bs, 2.5), d_CI95_high=np.percentile(bs, 97.5),
                               P=(np.nan if gate_mk else mwu(c, u)),
                               notes=('gate marker, not tested (HLA-DR > p75 is part of the DC gate)' if gate_mk else
                                      ('IL-1b signal close to isotype background; interpret with caution'
                                       if p in DAGGER else ''))))
    FR = bh(pd.DataFrame(fr))
    FR.loc[FR['P'].isna(), 'sig'] = 'not tested'
    # Tregs vs conventional CD4+ T cells of the same FOV (paired Wilcoxon); Tregs UC vs CD (Mann-Whitney)
    tc, tpair, treg_n = [], {}, {}
    for cond in ['UC', 'CD']:
        X = D[cond]
        t, cv = X['M']['Treg'], X['M']['Conv CD4 T']
        qual = [f for f in np.unique(X['fov'])
                if (t & (X['fov'] == f)).sum() >= MIN_TREG and (cv & (X['fov'] == f)).sum() >= MIN_CONV]
        inq = np.isin(X['fov'], qual)
        treg_n[cond] = dict(n_FOV=len(qual), fovs=';'.join(map(str, qual)), n_treg=int((t & inq).sum()),
                            n_conv=int((cv & inq).sum()))
        for p in TREG_MARKERS:
            a = np.array([X['V'][p][t & (X['fov'] == f)].mean() for f in qual])
            b = np.array([X['V'][p][cv & (X['fov'] == f)].mean() for f in qual])
            tpair[(p, cond)] = (a, b)
            for f, va, vb in zip(qual, a, b):
                for cls, v, mm in [('Treg', va, t), ('Conv CD4 T', vb, cv)]:
                    perfov.append(dict(record='Treg panel per FOV mean (paired FOVs)', subset=cls, protein=p,
                                       disease=cond, location='tissue', fov=f, value=v,
                                       n_cells=int((mm & (X['fov'] == f)).sum())))
    for cond in ['UC', 'CD']:
        n = treg_n[cond]['n_FOV']
        for p in TREG_MARKERS:
            a, b = tpair[(p, cond)]
            tc.append(dict(record='Treg vs conv CD4 T (paired Wilcoxon signed-rank over FOVs)',
                           family=f'Treg vs conv CD4 T | {cond}', subset='Treg vs Conv CD4 T', protein=p, disease=cond,
                           location='tissue', n_FOV=n, fovs=treg_n[cond]['fovs'], n_cells_Treg=treg_n[cond]['n_treg'],
                           n_cells_conv=treg_n[cond]['n_conv'], mean_Treg=a.mean(), mean_conv=b.mean(),
                           median_Treg=np.median(a), median_conv=np.median(b), mean_diff=(a - b).mean(),
                           n_FOV_Treg_higher=int((a > b).sum()),
                           P=wilcoxon(a, b)[1] if np.any(a != b) else np.nan,
                           notes=(f'{n} paired FOVs: smallest attainable two-sided Wilcoxon P = {2 / 2 ** n:.4g}'
                                  if 2 / 2 ** n > 0.05 / len(TREG_MARKERS) else '')))
    for p in TREG_MARKERS:
        u, c = tpair[(p, 'UC')][0], tpair[(p, 'CD')][0]
        tc.append(dict(record='Treg UC vs CD (Mann-Whitney over FOVs)', family='Treg UC vs CD', subset='Treg', protein=p,
                       location='tissue', n_FOV_UC=len(u), n_FOV_CD=len(c), n_cells_UC=treg_n['UC']['n_treg'],
                       n_cells_CD=treg_n['CD']['n_treg'], mean_UC=u.mean(), mean_CD=c.mean(), median_UC=np.median(u),
                       median_CD=np.median(c), P=mwu(c, u), notes=''))
    TC = bh(pd.DataFrame(tc))
    TC['notes'] = TC['notes'].fillna('')
    TC['notes'] = [n + (f'; q = {q:.4g} before rounding (< 0.05)' if q < 0.05 <= float(f'{q:.3g}') else '')
                   for n, q in zip(TC['notes'], TC['q_BH'])]
    out = pd.concat([FR, TC, pd.DataFrame(perfov), pd.DataFrame(forest_fov), pd.DataFrame(ncell)],
                    ignore_index=True)
    lead = ['record', 'family', 'subset', 'function_group', 'protein', 'location', 'disease', 'fov']
    out = round_table(out[[c for c in lead if c in out] + [c for c in out if c not in lead]], int_cols=('fov',))
    out['notes'] = out['notes'].fillna('')
    out.loc[0, 'notes'] = ('protein only; values = geometric-mean-normalised protein intensity (_geo, a.u.); thresholds '
                           'within section over all cells (+ > p75); unit = FOV; '
                           f'>= {MIN_CELLS} subset cells per FOV (forest); Treg panel: FOVs with >= {MIN_TREG} Tregs and >= {MIN_CONV} '
                           'conv CD4 T, paired Wilcoxon signed-rank (Treg vs conv) and Mann-Whitney (Treg UC vs CD), BH over '
                           f'the {len(TREG_MARKERS)} Treg markers ({", ".join(TREG_MARKERS)}); '
                           'Mann-Whitney two-sided; BH within family; UC = rectum 385R (15 FOVs), CD = ileum 375i '
                           '(17-FOV subset); one biopsy section per disease, so disease and site are confounded; '
                           f'DC gate = {DC_GATE}, as in the cell-type maps (strict needs >= {MIN_DC} cells '
                           f'in each section; strict n UC {N_STRICT["UC"]}, CD {N_STRICT["CD"]}); '
                           'DC gate contains HLA-DR, so DC HLA-DR is not tested')
    out.to_csv(CSV_MK, index=False)
    return dict(D=D, FR=FR, TC=TC, tpair=tpair, FOREST=FOREST, FOREST_SUB=FOREST_SUB, TREG_MARKERS=TREG_MARKERS,
                DAGGER=DAGGER)


V_CT, V_MK = compute_cell_type_profiles(), compute_marker_panels()

# ---- figure: maps on the left; depth profiles, forest plot and Treg box plots on the right ----
COL = {'UC': '#457B9D', 'CD': '#E63946'}


def qtxt(q):
    if not np.isfinite(q):
        return 'q n/a'
    if q < 0.001:
        e = int(np.floor(np.log10(q)))
        return f'q = {q / 10 ** e:.1f}×10$^{{{e}}}$'
    t = f'{q:#.3g}'
    if q < 0.05 <= float(t):
        t = f'{q:#.4g}'
    return f'q = {t}'


D_CT, D_MK = V_CT['D'], V_MK['D']
ext = {c: (D_CT[c]['xy'].min(0), D_CT[c]['xy'].max(0)) for c in D_CT}
span = {c: ext[c][1] - ext[c][0] for c in D_CT}
H_UM = max(s[1] for s in span.values()) + 60
W = {c: span[c][0] + 40 for c in D_CT}
fig = plt.figure(figsize=(27, 16.5), layout='constrained')
fig_l, fig_r = fig.subfigures(1, 2, width_ratios=[0.8, 1])
sf_top, sf_bot = fig_l.subfigures(2, 1)
sf_lp, sf_for, sf_treg = fig_r.subfigures(3, 1, height_ratios=[1.8, 1, 1])


def map_row(sf, D, layers, handles, title, titles=False, **leg_kw):
    """UC and CD sections side by side at the same scale; layers = (mask, colour, size, zorder, edge, lw)."""
    for ax, cond in zip(sf.subplots(1, 2, width_ratios=[W['UC'], W['CD']]), ['UC', 'CD']):
        xy = D[cond]['xy'] - ext[cond][0]
        for m, c, s, z, ec, lw in layers(D[cond]['M']):
            ax.scatter(xy[m, 0], xy[m, 1], s=s, c=c, lw=lw, edgecolors=(ec if ec else 'none'), zorder=z,
                       rasterized=True)
        ax.set_xlim(-20, W[cond] - 20)
        ax.set_ylim(H_UM + 110, -40)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.plot([0, 500], [H_UM + 80] * 2, color='black', lw=4, solid_capstyle='butt')
        ax.text(250, H_UM + 50, '500 µm', ha='center', va='bottom', fontsize=16)
        if titles:
            ax.set_title(cond, fontsize=18)
    sf.legend(handles=handles, loc='outside right center', frameon=False, fontsize=15, handletextpad=0.3,
              title=title, title_fontsize=15, alignment='left', **leg_kw)


def dot(c, ms, label, ec=None, mew=0):
    return Line2D([], [], marker='o', ls='', color=c, mec=(ec if ec else c), mew=mew, ms=ms, label=label)


# top-left: cell-type maps (Fig 7A)
C_OTH, C_EPI = '#E6E6E6', '#B9B9B9'
C_B, C_T, C_STR, C_MAC, C_DC, C_TREG = '#7B4FA6', '#A6DBA0', '#3A86A8', '#F28E2B', '#C2185B', '#1B7837'


def layers_ct(M):
    shown = (M['Epithelium'] | M['Stroma'] | M['B cell'] | M['Other T'] | M['Treg']
             | M['Macrophage'] | M['DC'])
    return [(~shown, C_OTH, 1.0, 1, None, 0), (M['Epithelium'], C_EPI, 2.5, 2, None, 0),
            (M['B cell'], C_B, 3.5, 3, None, 0), (M['Other T'], C_T, 3.0, 4, None, 0),
            (M['Stroma'], C_STR, 11, 5, None, 0),
            (M['Macrophage'], C_MAC, 22, 6, 'black', 0.35),
            (M['DC'], C_DC, 26, 7, 'black', 0.4), (M['Treg'], C_TREG, 24, 8, None, 0)]


DC_LAB = ('DCs (CD11c⁺HLA-DR⁺\nCD68⁻CD14⁻Lin⁻)' if DC_GATE == 'DC relaxed'
          else 'DCs (CD45⁺Lin⁻CD11c-High\nCD68⁻CD14⁻CD163⁻)')
h7 = [dot(C_TREG, 11, 'Tregs (CD3⁺CD4⁺\nFOXP3⁺CTLA4⁺)'),
      dot(C_DC, 11, DC_LAB, 'black', 0.6),
      dot(C_MAC, 11, 'Macrophages (CD45⁺Lin⁻,\nCD68⁺ or CD163⁺,\nCD11c Not High)', 'black', 0.6),
      dot(C_STR, 10, 'Fibroblast-Rich Stroma\n(Fibronectin⁺Vimentin⁺\nSMA⁻CD45⁻EpCAM⁻CD31⁻)'),
      dot(C_T, 8, 'Other T Cells\n(CD3⁺, Not Treg)'),
      dot(C_B, 8, 'B Cells (CD20⁺ or CD19⁺,\nCD3⁻)'),
      dot(C_EPI, 8, 'Epithelium (EpCAM⁺)'),
      dot('#CFCFCF', 5, 'Other Cells')]
map_row(sf_top, D_CT, layers_ct, h7, 'Lin⁻ = CD3⁻CD19⁻CD20⁻', titles=True, labelspacing=0.75)

# bottom-left: HLA-DR / iNOS myeloid maps (Fig 7C)
C_EPI6, C_HLA, C_INOS, C_MY, C_OTH6 = '#C9C9C9', '#B8860B', '#D7191C', '#F3E3C3', '#ECECEC'
C_MY_EDGE = '#9E9E9E'


def layers_my(M):
    hla, ino, my = M['HLA-DR-high myeloid'], M['iNOS-high myeloid'], M['Myeloid']
    shown = M['Stroma'] | M['Epithelium'] | my
    return [(~shown, C_OTH6, 0.6, 1, None, 0), (M['Epithelium'] & ~M['Stroma'] & ~my, C_EPI6, 3, 2, None, 0),
            (M['Stroma'] & ~my, C_STR, 9, 3, None, 0), (my & ~hla & ~ino, C_MY, 24, 5, C_MY_EDGE, 0.6),
            (hla & ~ino, C_HLA, 30, 6, None, 0), (ino & ~hla, C_INOS, 34, 7, None, 0),
            (ino & hla, C_INOS, 34, 8, C_HLA, 1.6)]


h6 = [dot(C_STR, 9, 'Fibroblast-Rich Stroma'), dot(C_EPI6, 9, 'Epithelium (EpCAM)'),
      dot(C_HLA, 11, 'HLA-DR-High Myeloid'), dot(C_INOS, 11, 'iNOS-High Myeloid'),
      dot(C_INOS, 11, 'HLA-DR-High and iNOS-High', C_HLA, 2),
      dot(C_MY, 10, 'Other Myeloid', C_MY_EDGE, 0.8), dot('#BDBDBD', 4, 'Other Cells')]
map_row(sf_bot, D_MK, layers_my, h6, 'Myeloid = Macrophage or DC\n(High = Above p75 of Section)')

# top-right: depth profiles (Fig 7B)
PR, FNAME, BIN_LAB, PLOT_FEATS = V_CT['PR'], V_CT['FNAME'], V_CT['BIN_LAB'], V_CT['PLOT_FEATS']
FCOL = {'Stroma': C_STR, 'APC': '#555555', 'Treg': C_TREG, 'Other T': '#5AAE61'}
PTITLE = {'Stroma': 'Fibroblast-Rich Stroma', 'APC': 'All APCs (DC + Macrophage + B)',
          'Treg': 'Tregs (FOXP3⁺CTLA4⁺)', 'Other T': 'Other T Cells (Control)'}
xs = np.arange(len(BIN_LAB))
LS, MK = {'UC': '-', 'CD': '--'}, {'UC': 'o', 'CD': 's'}
hatch_who = set()
for i, (kf, ax) in enumerate(zip(PLOT_FEATS, sf_lp.subplots(2, 2).flat)):
    hatched = set()
    for cond in ['UC', 'CD']:
        p = PR[(PR.disease == cond) & (PR.feature == FNAME[kf])].set_index('depth_bin_um').loc[BIN_LAB]
        y, e = p.fov_mean.values.astype(float), p.fov_SE.values.astype(float)
        keep = p.shown.values.astype(bool)
        y[~keep], e[~keep] = np.nan, np.nan
        ax.errorbar(xs, y, yerr=e, color=FCOL[kf], ls=LS[cond], marker=MK[cond], ms=9, lw=2.6,
                    capsize=4, mfc=('white' if cond == 'CD' else FCOL[kf]), mew=1.9)
        hatched |= set(np.where(~keep)[0])
        if (~keep).any():
            hatch_who.add(cond)
    for j in hatched:
        ax.axvspan(j - 0.5, j + 0.5, facecolor='none', edgecolor='#B0B0B0', hatch='///', lw=0, zorder=0)
    ax.set_ylim(bottom=0)
    ax.set_xlim(-0.5, 4.5)
    ax.set_title(PTITLE[kf], fontsize=18, color=FCOL[kf], pad=5)
    ax.set_xticks(xs)
    ax.set_xticklabels(BIN_LAB if i >= 2 else [], rotation=35, ha='right')
    ax.tick_params(labelsize=16)
    ax.spines[['top', 'right']].set_visible(False)
who = ' and '.join(sorted(hatch_who)) if hatch_who else 'None'
lh = [Line2D([], [], color='black', ls='-', marker='o', ms=8, lw=2.2, label='UC'),
      Line2D([], [], color='black', ls='--', marker='s', mfc='white', ms=8, lw=2.2, label='CD'),
      Patch(facecolor='none', edgecolor='#B0B0B0', hatch='///', label=f'{who} < 3 FOVs (Not Drawn)')]
sf_lp.legend(handles=lh, loc='outside upper center', ncol=3, frameon=False, fontsize=16, columnspacing=1.6)
sf_lp.supylabel('% of Cells in Distance Bin\n(Mean ± SE Across FOVs)', fontsize=17, multialignment='center')
sf_lp.supxlabel('Distance to Nearest Epithelial Cell (µm)', fontsize=17)

# bottom-right: forest plot (Fig 7D)
FR, TC, tpair = V_MK['FR'], V_MK['TC'], V_MK['tpair']
FOREST, FOREST_SUB, DAGGER, TREG_MARKERS = V_MK['FOREST'], V_MK['FOREST_SUB'], V_MK['DAGGER'], V_MK['TREG_MARKERS']
LEAD = {'iNOS', 'HLA-DR'}
FAMDISP = {'Inflammatory': 'Inflammatory',
           'Antigen presentation / co-stimulation': 'Antigen Presentation /\nCo-Stimulation'}
DISP = {'NF-kB p65': 'NF-κB p65', 'IL-1b': 'IL-1β'}
TREG_ROLE = {'ICOS': 'Co-Stimulatory Receptor', 'GITR': 'Co-Stimulatory Receptor', 'CD27': 'Stability / Memory'}
rows_y, ylab, y = [], [], 0
for fam, mk in FOREST.items():
    ylab.append((y, fam, 'header'))
    y += 1
    for p in mk:
        rows_y.append((y, p))
        ylab.append((y, p, 'marker'))
        y += 1
    y += 0.4
ymax = y - 0.4
lim = np.nanmax(np.abs(FR.loc[FR['P'].notna(), ['d_CI95_low', 'd_CI95_high']].values))
lim = np.ceil(lim + 0.2)
STITLE = {'Macrophage': 'Macrophages', 'DC': 'DCs'}
for i, (s, ax) in enumerate(zip(FOREST_SUB, sf_for.subplots(1, 2))):
    sub = FR[FR.subset == s].set_index('protein')
    ax.axvline(0, color='#888888', lw=1.2, zorder=0)
    ykw = dict(xycoords=ax.get_yaxis_transform(), textcoords='offset points', va='center', fontsize=16)
    qkw = dict(xycoords=ax.get_yaxis_transform(), xytext=(8, 0), textcoords='offset points', ha='left',
               va='center', fontsize=15)
    for yy, p in rows_y:
        r_ = sub.loc[p]
        if not np.isfinite(r_['P']):
            ax.text(0, yy, 'Gate Marker, Not Tested', ha='center', va='center', fontsize=15, color='#777777',
                    style='italic', zorder=3, bbox=dict(facecolor='white', edgecolor='none', pad=1))
            ax.annotate('–', (1, yy), color='#777777', **qkw)
            continue
        sig = r_['q_BH'] < 0.05
        col = COL['CD'] if r_['cohen_d_CD_minus_UC'] > 0 else COL['UC']
        ax.plot([r_['d_CI95_low'], r_['d_CI95_high']], [yy, yy], color=col, lw=2.4, zorder=2)
        ax.plot(r_['cohen_d_CD_minus_UC'], yy, marker='o', ms=12, mfc=(col if sig else 'white'), mec=col,
                mew=2.2, zorder=3)
        ax.annotate(f"{qtxt(r_['q_BH'])} {r_['sig'] if sig else ''}".strip(), (1, yy), **qkw)
    for yy, lab, kind in ylab:
        if kind == 'header':
            ax.axhspan(yy - 0.5, yy + 0.5, color='#EFEFEF', lw=0, zorder=0)
            if i == 0:
                ax.annotate(FAMDISP[lab], (0, yy), xytext=(-7, 3 if '\n' in FAMDISP[lab] else 0), ha='right',
                            fontweight='bold', linespacing=0.95, **ykw)
        elif i == 0:
            ax.annotate(DISP.get(lab, lab), (0, yy), xytext=(-16, 0), ha='right',
                        fontweight=('bold' if lab in LEAD else 'normal'), **ykw)
            if lab in DAGGER:
                ax.annotate('†', (0, yy), xytext=(-5, 0), ha='right', color='#888888', **ykw)
    ax.set_ylim(ymax - 0.5, -0.5)
    ax.set_yticks([])
    ax.set_xlim(-lim, lim)
    ax.set_xticks(np.arange(-lim + 1, lim, 2))
    ax.set_xlabel("Cohen's d, CD − UC (95% CI)")
    ax.set_title(STITLE[s], fontsize=18)
    for txt, x0, ha, c in [('← Higher in UC', 0.02, 'left', 'UC'), ('Higher in CD →', 0.98, 'right', 'CD')]:
        ax.annotate(txt, (x0, 0), xycoords='axes fraction', xytext=(0, -54), textcoords='offset points',
                    ha=ha, va='top', fontsize=15, color=COL[c])
    ax.spines[['top', 'right', 'left']].set_visible(False)

# bottom-right: Treg box plots (Fig 7E)
POS = {('UC', 'conv'): 0.0, ('UC', 'Treg'): 1.0, ('CD', 'conv'): 2.5, ('CD', 'Treg'): 3.5}
TQ = TC.set_index(['family', 'protein'])
for p, ax in zip(TREG_MARKERS, sf_treg.subplots(1, 3)):
    allv = np.concatenate([np.concatenate(tpair[(p, c)]) for c in ['UC', 'CD']])
    for cond in ['UC', 'CD']:
        a, b = tpair[(p, cond)]
        x0, x1 = POS[(cond, 'conv')], POS[(cond, 'Treg')]
        bp = ax.boxplot([b, a], positions=[x0, x1], widths=0.62, patch_artist=True, showfliers=False,
                        medianprops=dict(color='black', linewidth=1.5))
        bp['boxes'][0].set_facecolor('white')
        bp['boxes'][0].set_edgecolor(COL[cond])
        bp['boxes'][0].set_linewidth(2)
        bp['boxes'][1].set_facecolor(COL[cond])
        bp['boxes'][1].set_alpha(0.85)
        for ya, yb in zip(a, b):
            ax.plot([x0, x1], [yb, ya], color='#9E9E9E', lw=0.9, zorder=2)
        for xx, v in [(x0, b), (x1, a)]:
            ax.scatter(np.full(len(v), xx), v, s=50, color='white', edgecolor='black', lw=0.9, zorder=3)
    lo, hi = allv.min(), allv.max()
    rg = hi - lo
    y1 = hi + 0.06 * rg
    for cond in ['UC', 'CD']:
        r_ = TQ.loc[(f'Treg vs conv CD4 T | {cond}', p)]
        x0, x1 = POS[(cond, 'conv')], POS[(cond, 'Treg')]
        ax.plot([x0, x0, x1, x1], [y1 - 0.02 * rg, y1, y1, y1 - 0.02 * rg], color='black', lw=1.2)
        lab = r_['sig'] if (cond == 'CD' or r_['q_BH'] < 0.05) else qtxt(r_['q_BH'])
        ax.text((x0 + x1) / 2, y1 + 0.01 * rg, lab, ha='center', va='bottom',
                fontsize=(16 if lab == r_['sig'] else 14))
    y2 = hi + 0.27 * rg
    r_ = TQ.loc[('Treg UC vs CD', p)]
    ax.plot([1, 1, 3.5, 3.5], [y2 - 0.02 * rg, y2, y2, y2 - 0.02 * rg], color='black', lw=1.2)
    ax.text(2.25, y2 + 0.01 * rg, f"Treg UC vs CD {r_['sig']}", ha='center', va='bottom', fontsize=15)
    ax.set_ylim(lo - 0.05 * rg, hi + 0.45 * rg)
    ax.set_title(f'Treg {p}\n{TREG_ROLE[p]}', fontsize=18)
    ax.set_xticks(list(POS.values()))
    ax.set_xticklabels(['Conv\nCD4', 'Treg', 'Conv\nCD4', 'Treg'], fontsize=16)
    ax.set_xlim(-0.6, 4.1)
    for cond, xm in [('UC', 0.5), ('CD', 3.0)]:
        ax.annotate(cond, (xm, 0), xycoords=ax.get_xaxis_transform(), xytext=(0, -50),
                    textcoords='offset points', ha='center', va='top', fontsize=17, color=COL[cond],
                    fontweight='bold')
    ax.set_ylabel(f'{p}\n(Mean per FOV, a.u.)')
    ax.spines[['top', 'right']].set_visible(False)

fig.savefig(PNG, dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.close(fig)
print('wrote', PNG, CSV_CT, CSV_MK)
