#!/usr/bin/env python
"""CosMx spatial proteomics, UC vs CD: protein-gated cell types, myeloid activation and Treg receptors.

One script draws the whole figure. It gates cell types on protein (percentile thresholds within each section), maps them,
profiles composition by distance from the epithelium, compares inflammatory and antigen-presentation proteins in macrophages and
DCs, and compares ICOS, GITR and CD27 on Tregs and conventional CD4+ T cells.

Panels made
  Fig 7A  spatial maps of protein-defined cell types
  Fig 7B  cell-type composition by distance from the closest epithelial cell (per FOV, mean +- SE)
  Fig 7C  macrophages and DCs coloured by HLA-DR and iNOS
  Fig 7D  iNOS, IL-1b, IL-18, HLA-DR, B7-H3, IDO1 in macrophages and DCs, Cohen's d (CD - UC) of per-FOV means
  Fig 7E  ICOS, GITR and CD27 on Tregs vs conventional CD4+ T cells from the same FOV
  Supp Table 12 (gates and marker lists are defined in this script)

Inputs
  /share/fsmresfiles/UC/scOmics_model/data/Spatial_Protein_CosMx/RNA1_Cos4_385R_381R_S2_norm_wpr_anno_jd.h5ad (UC section; 03_uc_rna_protein_processing.ipynb)
  /share/fsmresfiles/UC/AtoMx/CD_multiOmics/Cos4_315i_375i/rna_norm_pr_geonorm_upd_anno.h5ad (CD section; 04_cd_rna_protein_processing.ipynb)

Outputs
  /share/fsmresfiles/UC/revision_jenny/Fig5_spatial_protein/fig5_combined_v1.png (Fig 7A-E)
  /share/fsmresfiles/UC/revision_jenny/Fig5_spatial_protein/fig7_cell_types_depth_profiles.csv (gates, thresholds, cell counts,
    depth profiles and zone tests behind Fig 7A-B)
  /share/fsmresfiles/UC/revision_jenny/Fig5_spatial_protein/fig7_myeloid_treg_markers.csv (gates, cell counts, forest-plot and
    Treg statistics behind Fig 7C-E)

Run order of sections
  F7_protein_uc_cd

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist. Command-line arguments and environment variables used by a section are listed in its own
  docstring.
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
# Section F7_protein_uc_cd: protein-gated cell types, depth profiles, myeloid markers and Treg receptors, UC vs CD
# Panels: Fig 7A-E, Supp Table 12
# Note: optional argument: output PNG path (default: fig5_combined_v1.png in the output folder)
# ======================================================================================================================
if _run('F7_protein_uc_cd'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig5_spatial_protein/fig5_combined_v1.py'   # working folder of this step (relative paths resolve from here)
    """
Fig 7A-E as one matplotlib figure; every panel is drawn on one canvas, so all fonts are set in points.

Computation (two blocks, both run on the UC section 385R and the 17-FOV CD subset of section 375i):
  compute_cell_type_profiles  protein gates (percentile thresholds within each section, `_geo` values) for epithelium,
                              fibroblast-rich stroma, macrophages, DCs (strict gate if >= 50 cells in each section,
                              otherwise relaxed), B cells, Tregs and other T cells; distance of every cell to the nearest
                              epithelial cell; % of cells per distance bin per FOV (bins with >= 20 cells; a point is
                              shown when >= 3 FOVs contribute); FOV-paired Wilcoxon zone tests (BH; CSV only).
                              -> fig7_cell_types_depth_profiles.csv (Fig 7A-B)
  compute_marker_panels       the same gates plus HLA-DR-high / iNOS-high myeloid cells; per-FOV mean protein in
                              macrophages and DCs, Cohen's d (CD - UC) with 5,000-draw bootstrap 95% CI and Mann-Whitney
                              (BH within family); ICOS, GITR and CD27 on Tregs vs conventional CD4+ T cells of the same FOV
                              (paired Wilcoxon) and Tregs UC vs CD (Mann-Whitney), BH over the 3 markers.
                              -> fig7_myeloid_treg_markers.csv (Fig 7C-E)
Every plotted value comes from these two tables. The gates and marker lists form Supp Table 12.

Layout: 2 x 2.
  top-left      cell-type maps (titles 'UC' / 'CD'), legend to the right (legend title 'Lin- = CD3-CD19-CD20-')
  bottom-left   marker maps (no titles), legend to the right; same um per inch and same axis x positions as the
                top maps; both legends start at the same x
  top-right     4 depth profiles (2 x 2) with the shared y / x labels and the UC / CD key
  bottom-right  forest plot (Macrophages | DCs) above the 3 Treg box plots; both rows share left / right edges
Every axis is placed in inches: overhangs are measured on the rendered figure and the right column is sized so its top
and bottom match the left column; forest and Treg axes keep a constant height ratio with minimum heights 2.8 / 2.45 in
(the block boundary is printed at run time). Layout asserts at the end check that no blocks overlap.

Run:  OPENBLAS_CORETYPE=SkylakeX /data/ydn4687/geneformer_env/bin/python <this file> [png_path]
"""
    import io
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
    from matplotlib.transforms import Bbox
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

    # ── Load ──
    CD_PATH = '/share/fsmresfiles/UC/AtoMx/CD_multiOmics/Cos4_315i_375i'
    UC_PATH = '/share/fsmresfiles/UC/scOmics_model/data/Spatial_Protein_CosMx'
    adata_rna_qc = sc.read_h5ad(f'{CD_PATH}/rna_norm_pr_geonorm_upd_anno.h5ad')
    adata_rna_qc_375 = adata_rna_qc[(adata_rna_qc.obs['fov'] <= 29)].copy()
    fovs = [6, 7, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29]
    adata_rna_qc_375_sub = adata_rna_qc_375[adata_rna_qc_375.obs['fov'].isin(fovs)].copy()
    adata_rna_qc_385 = sc.read_h5ad(f'{UC_PATH}/RNA1_Cos4_385R_381R_S2_norm_wpr_anno_jd.h5ad')
    print('CD cells', adata_rna_qc_375_sub.n_obs, '| UC cells', adata_rna_qc_385.n_obs)
    DATASETS = [('UC', adata_rna_qc_385), ('CD', adata_rna_qc_375_sub)]
    px_per_um = np.sqrt(3266 / 47.25)


    # ════════════════════════ computation: cell-type maps and depth profiles ════════════════════════
    def compute_cell_type_profiles():
        MIN_BIN = 20
        MIN_NFOV = 3
        MIN_DC = 50
        BINS = [0, 25, 50, 100, 200, np.inf]
        BIN_LAB = ['0–25', '25–50', '50–100', '100–200', '>200']
        FEATS = ['Stroma', 'APC', 'DC', 'Macrophage', 'Treg', 'Other T']
        PLOT_FEATS = ['Stroma', 'APC', 'Treg', 'Other T']
        FNAME = {'Stroma': 'Fibroblast-rich stroma', 'APC': 'All APCs (DC + macrophage + B)',
                 'DC': 'DCs', 'Macrophage': 'Macrophages', 'Treg': 'Tregs (FOXP3⁺CTLA4⁺)',
                 'Other T': 'Other T cells'}
        GATE_MARKERS = ['CD45', 'CD3', 'CD20', 'CD19', 'CD68', 'CD163', 'CD11c', 'CD14', 'CD15', 'CD16',
                        'CD123']
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
            thr = {f'{m} p{q}': np.percentile(g(m), q) for m, qs in THR_MARKERS.items() for q in qs}
            return dict(Epithelium=epi, Stroma=stroma, Macrophage=mac, **{'DC strict': dc_strict,
                        'DC relaxed': dc_relax, 'B cell': bcell, 'Treg': treg, 'Other T': oth_t}), thr

        def dist_to_epi(xy, epi):
            d, _ = cKDTree(xy[epi]).query(xy, k=2)
            out = d[:, 0].copy()
            out[epi] = d[epi, 1]
            return out / px_per_um

        def sig_of(q):
            return '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'

        def se(x):
            x = np.asarray(x, float)
            return np.std(x, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan

        D, rows = {}, []
        for cond, adata in DATASETS:
            xy = np.asarray(adata.obsm['spatial_fov'])[:, :2].astype(float) / px_per_um
            M, thr = masks_of(adata)
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
        n_strict = {c: int(D[c]['M']['DC strict'].sum()) for c in D}
        DC_GATE = 'DC strict' if min(n_strict.values()) >= MIN_DC else 'DC relaxed'
        for c in D:
            M = D[c]['M']
            M['DC'] = M[DC_GATE]
            M['APC'] = M['DC'] | M['Macrophage'] | M['B cell']
            for k in ['Epithelium', 'Stroma', 'Macrophage', 'DC strict', 'DC relaxed', 'B cell', 'Treg',
                      'Other T', 'APC']:
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
        ZONES = {'CD': [('0–50 µm', 0, 50), ('> 100 µm', 100, np.inf)],
                 'UC': [('0–50 µm', 0, 50), ('50–100 µm', 50, 100)]}
        PAIR, tst = {}, []
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
                PAIR[(cond, k)] = (A, B)
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
        out = out[[c for c in lead if c in out] + [c for c in out if c not in lead]]
        for c in out.columns:
            if c.startswith('n_'):
                out[c] = out[c].astype('Int64')
            elif out[c].dtype.kind == 'f':
                out[c] = out[c].map(lambda x: x if pd.isna(x) else float(f'{x:.3g}'))
        out['notes'] = ''
        out.loc[0, 'notes'] = ('protein only; thresholds within section over all cells; distance = nearest '
                           'epithelial cell (EpCAM > p75; other than itself), um; unit = FOV; zone/bin '
                           f'needs >= {MIN_BIN} cells per FOV; profile point shown if >= {MIN_NFOV} FOVs; '
                           f'DC gate used = {DC_GATE} (strict needs >= {MIN_DC} cells in each section; '
                           f'strict n UC {n_strict["UC"]}, CD {n_strict["CD"]}); '
                           'UC = rectum 385R (15 FOVs), CD = ileum 375i (17-FOV subset); one mucosal biopsy '
                           'section per disease, so disease and location are confounded; descriptive; '
                           'zone tests are reported here only (not plotted); depth_profile rows with '
                           'in_figure = True are plotted in panel B')
        buf = io.StringIO()
        out.to_csv(buf, index=False)
        return dict(D=D, rows=rows, PR=PR, DC_GATE=DC_GATE, FNAME=FNAME, PLOT_FEATS=PLOT_FEATS, BIN_LAB=BIN_LAB,
                    csv=buf.getvalue())


    # ════════════════════════ computation: myeloid markers and Treg receptors ════════════════════════
    def compute_marker_panels():
        MIN_CELLS, MIN_TREG, MIN_FOV_TREG, MIN_CONV, MIN_DC = 3, 3, 4, 3, 50
        TREG_MARKERS = ['ICOS', 'GITR', 'CD27']
        FOREST = {'Inflammatory': ['iNOS', 'IL-1b', 'IL-18'],
                  'Antigen presentation / co-stimulation': ['HLA-DR', 'B7-H3', 'IDO1']}
        FOREST_SUB = ['Macrophage', 'DC']
        DAGGER = {'IL-1b'}
        N_BOOT = 5000
        GATE_MARKERS = ['CD45', 'CD3', 'CD20', 'CD19', 'CD68', 'CD163', 'CD11c', 'CD14', 'CD15', 'CD16', 'CD123']
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

        def gate(adata):
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
            lab, pos = gate(adata)
            g = lambda p: adata.obs[f'{p}_geo'].values.astype(float)
            above = lambda p, q: g(p) > np.percentile(g(p), q)
            notabove = lambda p, q: g(p) <= np.percentile(g(p), q)
            epi = above('EpCAM', 75)
            stroma = (above('Fibronectin', 60) & above('Vimentin', 50) & notabove('SMA', 85)
                      & notabove('CD45', 75) & notabove('EpCAM', 75) & notabove('CD31', 75))
            mac = lab == 'Macrophage'
            dc_strict = lab == 'DC'
            dc_relax = (above('CD11c', 75) & above('HLA-DR', 75) & notabove('CD68', 75) & notabove('CD14', 75)
                        & ~pos['CD3'] & ~pos['CD20'] & ~pos['CD19'] & ~mac)
            treg = above('CD3', 75) & above('CD4', 75) & above('FOXP3', 75) & above('CTLA4', 75)
            conv = pos['CD3'] & above('CD4', 75) & ~treg
            return {'Epithelium': epi, 'Stroma': stroma, 'Macrophage': mac, 'DC strict': dc_strict,
                    'DC relaxed': dc_relax, 'Treg': treg, 'Conv CD4 T': conv,
                    '_hla': above('HLA-DR', 75), '_inos': above('iNOS', 75)}

        def add_derived_masks(M, dc_gate):
            M['DC'] = M[dc_gate]
            M['Myeloid'] = M['Macrophage'] | M['DC']
            M['HLA-DR-high myeloid'] = M['Myeloid'] & M.pop('_hla')
            M['iNOS-high myeloid'] = M['Myeloid'] & M.pop('_inos')
            return M

        def sig_of(q):
            return 'ns' if not np.isfinite(q) else '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'

        def mwu(a, b):
            a, b = np.asarray(a, float), np.asarray(b, float)
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            return mannwhitneyu(a, b, alternative='two-sided')[1] if min(len(a), len(b)) >= 3 else np.nan

        def bh(df):
            df['q_BH'] = np.nan
            for fam, ix in df.groupby('family').groups.items():
                p = df.loc[ix, 'P']
                ok = p.notna()
                if ok.any():
                    df.loc[p[ok].index, 'q_BH'] = multipletests(p[ok], method='fdr_bh')[1]
            df['sig'] = df['q_BH'].map(sig_of)
            return df

        D, ncell = {}, []
        ALLP = sorted(set(TREG_MARKERS) | set(sum(FOREST.values(), [])))
        MM = {cond: masks_of(adata) for cond, adata in DATASETS}
        n_strict = {c: int(MM[c]['DC strict'].sum()) for c in MM}
        DC_GATE = 'DC strict' if min(n_strict.values()) >= MIN_DC else 'DC relaxed'
        GATE_TXT['DC'] = GATE_TXT[DC_GATE] + ' (gate used, as in the cell-type maps)'
        for cond, adata in DATASETS:
            xy = np.asarray(adata.obsm['spatial_fov'])[:, :2].astype(float) / px_per_um
            M = add_derived_masks(MM[cond], DC_GATE)
            fov = adata.obs['fov'].values.astype(int)
            V = {p: adata.obs[f'{p}_geo'].values.astype(float) for p in ALLP}
            D[cond] = dict(xy=xy, M=M, fov=fov, V=V)
            for k, m in M.items():
                nf = sum((m & (fov == f)).sum() >= MIN_CELLS for f in np.unique(fov))
                ncell.append(dict(record='n_cells', disease=cond, subset=k, gate=GATE_TXT.get(k, ''), n_cells=int(m.sum()),
                                  pct_of_section=100 * m.mean(), n_FOV_ge3=nf))
            ncell.append(dict(record='n_cells', disease=cond, subset='all cells', n_cells=adata.n_obs,
                              n_FOV=len(np.unique(fov))))
            ncell.append(dict(record='n_cells_overlap', disease=cond, subset='DC & Epithelium',
                              n_cells=int((M['DC'] & M['Epithelium']).sum())))
        perfov = []

        def cohen_d(c, u):
            n1, n2 = len(c), len(u)
            sp = np.sqrt(((n1 - 1) * c.var(ddof=1) + (n2 - 1) * u.var(ddof=1)) / (n1 + n2 - 2))
            return (c.mean() - u.mean()) / sp if sp > 0 else np.nan

        rng_b = np.random.default_rng(20261109)
        fr, forest_fov = [], []
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
        tc, tpair, treg_n = [], {}, {}
        for cond in ['UC', 'CD']:
            X = D[cond]
            t, cv = X['M']['Treg'], X['M']['Conv CD4 T']
            qual = [f for f in np.unique(X['fov'])
                    if (t & (X['fov'] == f)).sum() >= MIN_TREG and (cv & (X['fov'] == f)).sum() >= MIN_CONV]
            assert len(qual) >= MIN_FOV_TREG, (cond, 'paired FOVs', qual)
            inq = np.isin(X['fov'], qual)
            treg_n[cond] = dict(n_FOV=len(qual), fovs=';'.join(map(str, qual)), n_treg=int((t & inq).sum()),
                                n_conv=int((cv & inq).sum()), n_treg_total=int(t.sum()), n_conv_total=int(cv.sum()))
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
        out = out[[c for c in lead if c in out] + [c for c in out if c not in lead]]
        for c in out.columns:
            if c.startswith('n_') or c == 'fov':
                out[c] = out[c].astype('Int64')
            elif out[c].dtype.kind == 'f':
                out[c] = out[c].map(lambda x: x if pd.isna(x) else float(f'{x:.3g}'))
        if 'notes' not in out:
            out['notes'] = ''
        out['notes'] = out['notes'].fillna('')
        out.loc[0, 'notes'] = ('protein only; values = geometric-mean-normalised protein intensity (_geo, a.u.); thresholds '
                           'within section over all cells (+ > p75); unit = FOV; '
                           f'>= {MIN_CELLS} subset cells per FOV (forest); Treg panel: FOVs with >= {MIN_TREG} Tregs and >= {MIN_CONV} '
                           'conv CD4 T, paired Wilcoxon signed-rank (Treg vs conv) and Mann-Whitney (Treg UC vs CD), BH over '
                           f'the {len(TREG_MARKERS)} Treg markers ({", ".join(TREG_MARKERS)}); '
                           'Mann-Whitney two-sided; BH within family; UC = rectum 385R (15 FOVs), CD = ileum 375i '
                           '(17-FOV subset); one biopsy section per disease, so disease and site are confounded; '
                           f'DC gate = {DC_GATE}, as in the cell-type maps (strict needs >= {MIN_DC} cells '
                           f'in each section; strict n UC {n_strict["UC"]}, CD {n_strict["CD"]}); '
                           'DC gate contains HLA-DR, so DC HLA-DR is not tested')
        buf = io.StringIO()
        out.to_csv(buf, index=False)
        return dict(D=D, FR=FR, TC=TC, tpair=tpair, treg_n=treg_n, DC_GATE=DC_GATE, FOREST=FOREST,
                    FOREST_SUB=FOREST_SUB, TREG_MARKERS=TREG_MARKERS, DAGGER=DAGGER, csv=buf.getvalue())


    V_CT, V_MK = compute_cell_type_profiles(), compute_marker_panels()
    for V, path in [(V_CT, CSV_CT), (V_MK, CSV_MK)]:
        with open(path, 'w') as fh:
            fh.write(V['csv'])
        print('wrote', path)
    # both blocks must use the same DC gate, coordinates and cell-type masks
    assert V_CT['DC_GATE'] == V_MK['DC_GATE']
    for c in ['UC', 'CD']:
        assert np.array_equal(V_CT['D'][c]['xy'], V_MK['D'][c]['xy'])
        for k in ['Epithelium', 'Stroma', 'Macrophage', 'DC', 'Treg']:
            assert np.array_equal(V_CT['D'][c]['M'][k], V_MK['D'][c]['M'][k]), (c, k)

    # FOV and cell counts for the figure legend
    for c in ['UC', 'CD']:
        sec = [r for r in V_CT['rows'] if r['record'] == 'section' and r['disease'] == c][0]
        M = V_CT['D'][c]['M']
        print(f"COUNTS {c}: {sec['n_FOV']} FOVs, {sec['n_cells']:,} cells, DC n = {int(M['DC'].sum())}, "
          f"Treg n = {int(M['Treg'].sum())}")


    # ════════════════════════════════════════════ figure ════════════════════════════════════════════
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


    # sizes (inches). One map height is used for both map rows.
    ROW_IN = 0.77 * 22.3 / 3.40
    MAP_H_IN = 1.55 * ROW_IN
    MAP_GAP_IN, LEG_GAP_IN = 0.2, 0.15          # UC-CD map gap; CD map to legend
    MAP_ROW_GAP_IN = 0.30                       # top map row to bottom map row
    COL_GAP_IN = 0.45                           # left column (legends) to right column
    FOREST_W = 3.7                              # forest axis width
    FOREST_GAP_IN, TREG_GAP_IN, ROW_GAP_IN = 0.45, 0.45, 0.3
    Q_OFF_PT = 8
    LP_COL_GAP_IN = 0.45                        # depth profiles: col-1 axis right to col-2 tick labels
    LP_ROW_GAP_IN = 0.12                        # depth profiles: row-1 axis bottom to row-2 title
    RIGHT_GAP_IN = 0.35                         # top-right block to bottom-right block
    FOREST_FRAC = 0.74 / (0.74 + 0.66)          # forest : Treg axis height ratio
    FS_PTITLE, FS_TICK, FS_LAB, FS_KEY, FS_MLEG = 18, 16, 17, 16, 15

    D_CT, D_MK = V_CT['D'], V_MK['D']
    ext = {c: (D_CT[c]['xy'].min(0), D_CT[c]['xy'].max(0)) for c in D_CT}
    span = {c: ext[c][1] - ext[c][0] for c in D_CT}
    H_UM = max(s[1] for s in span.values()) + 60
    W = {c: span[c][0] + 40 for c in D_CT}
    UM_PER_IN = (H_UM + 150) / MAP_H_IN
    MAP_W_IN = {c: W[c] / UM_PER_IN for c in D_CT}

    FW, FH = 60.0, 40.0                          # scratch canvas; figure size set after measuring
    fig = plt.figure(figsize=(FW, FH))
    BOX = {}


    def put(ax, x, ytop, w, h):
        """Place an axis in inches: x from the left, ytop from the top of the figure."""
        BOX[ax] = (x, ytop, w, h)
        ax.set_position([x / FW, 1 - (ytop + h) / FH, w / FW, h / FH])


    def over(ax):
        """Overhang (inches) of the axis' tight bbox beyond its box: left, top, right, bottom."""
        bb = ax.get_tightbbox(fig.canvas.get_renderer())
        x, ytop, w, h = BOX[ax]
        dpi = fig.dpi
        return (x - bb.x0 / dpi, bb.y1 / dpi - (FH - ytop), bb.x1 / dpi - (x + w), (FH - ytop - h) - bb.y0 / dpi)


    def scale_bar(ax):
        ax.plot([0, 500], [H_UM + 80] * 2, color='black', lw=4, solid_capstyle='butt')
        ax.text(250, H_UM + 50, '500 µm', ha='center', va='bottom', fontsize=16)


    def map_frame(ax, cond):
        ax.set_xlim(-20, W[cond] - 20)
        ax.set_ylim(H_UM + 110, -40)
        ax.set_aspect('equal')
        ax.axis('off')


    # ── top-left: cell-type maps ──
    C_OTH, C_EPI = '#E6E6E6', '#B9B9B9'
    C_B, C_T, C_STR, C_MAC, C_DC, C_TREG = '#7B4FA6', '#A6DBA0', '#3A86A8', '#F28E2B', '#C2185B', '#1B7837'
    map_top, map_bot = {}, {}
    for cond in ['UC', 'CD']:
        ax = fig.add_axes([0, 0, 0.1, 0.1])
        put(ax, 1, 1, MAP_W_IN[cond], MAP_H_IN)
        X, M = D_CT[cond], D_CT[cond]['M']
        xy = X['xy'] - ext[cond][0]
        shown = (M['Epithelium'] | M['Stroma'] | M['B cell'] | M['Other T'] | M['Treg']
                 | M['Macrophage'] | M['DC'])
        lay = [(~shown, C_OTH, 1.0, 1, None, 0), (M['Epithelium'], C_EPI, 2.5, 2, None, 0),
               (M['B cell'], C_B, 3.5, 3, None, 0), (M['Other T'], C_T, 3.0, 4, None, 0),
               (M['Stroma'], C_STR, 11, 5, None, 0),
               (M['Macrophage'], C_MAC, 22, 6, 'black', 0.35),
               (M['DC'], C_DC, 26, 7, 'black', 0.4), (M['Treg'], C_TREG, 24, 8, None, 0)]
        for m, c, s, z, ec, lw in lay:
            ax.scatter(xy[m, 0], xy[m, 1], s=s, c=c, lw=lw, edgecolors=(ec if ec else 'none'), zorder=z,
                       rasterized=True)
        map_frame(ax, cond)
        scale_bar(ax)
        ax.set_title(cond, fontsize=18)
        map_top[cond] = ax
    DC_LAB = ('DCs (CD11c⁺HLA-DR⁺\nCD68⁻CD14⁻Lin⁻)' if V_CT['DC_GATE'] == 'DC relaxed'
              else 'DCs (CD45⁺Lin⁻CD11c-High\nCD68⁻CD14⁻CD163⁻)')
    LEG7 = [(C_TREG, None, 11, 'Tregs (CD3⁺CD4⁺\nFOXP3⁺CTLA4⁺)'),
            (C_DC, 'black', 11, DC_LAB),
            (C_MAC, 'black', 11, 'Macrophages (CD45⁺Lin⁻,\nCD68⁺ or CD163⁺,\nCD11c Not High)'),
            (C_STR, None, 10, 'Fibroblast-Rich Stroma\n(Fibronectin⁺Vimentin⁺\nSMA⁻CD45⁻EpCAM⁻CD31⁻)'),
            (C_T, None, 8, 'Other T Cells\n(CD3⁺, Not Treg)'),
            (C_B, None, 8, 'B Cells (CD20⁺ or CD19⁺,\nCD3⁻)'),
            (C_EPI, None, 8, 'Epithelium (EpCAM⁺)'),
            ('#CFCFCF', None, 5, 'Other Cells')]
    h7 = [Line2D([], [], marker='o', ls='', color=c, mec=(ec if ec else c), mew=(0.6 if ec else 0), ms=ms,
                 label=lab) for c, ec, ms, lab in LEG7]
    lax_t = fig.add_axes([0, 0, 0.1, 0.1])
    put(lax_t, 1, 1, 1, MAP_H_IN)
    lax_t.axis('off')
    lax_t.legend(handles=h7, loc='center left', bbox_to_anchor=(0, 0.5), frameon=False, fontsize=FS_MLEG,
                 handletextpad=0.3, labelspacing=0.75, borderaxespad=0, borderpad=0,
                 title='Lin⁻ = CD3⁻CD19⁻CD20⁻', title_fontsize=FS_MLEG, alignment='left')

    # ── bottom-left: marker maps, no titles ──
    C_EPI6, C_HLA, C_INOS, C_MY, C_OTH6 = '#C9C9C9', '#B8860B', '#D7191C', '#F3E3C3', '#ECECEC'
    C_MY_EDGE = '#9E9E9E'
    for cond in ['UC', 'CD']:
        ax = fig.add_axes([0, 0, 0.1, 0.1])
        put(ax, 1, 1, MAP_W_IN[cond], MAP_H_IN)
        X, M = D_MK[cond], D_MK[cond]['M']
        xy = X['xy'] - ext[cond][0]
        hla, ino, my = M['HLA-DR-high myeloid'], M['iNOS-high myeloid'], M['Myeloid']
        shown = M['Stroma'] | M['Epithelium'] | my
        lay = [(~shown, C_OTH6, 0.6, 1, None), (M['Epithelium'] & ~M['Stroma'] & ~my, C_EPI6, 3, 2, None),
               (M['Stroma'] & ~my, C_STR, 9, 3, None), (my & ~hla & ~ino, C_MY, 24, 5, C_MY_EDGE),
               (hla & ~ino, C_HLA, 30, 6, None), (ino & ~hla, C_INOS, 34, 7, None),
               (ino & hla, C_INOS, 34, 8, C_HLA)]
        for m, c, s, z, ec in lay:
            ax.scatter(xy[m, 0], xy[m, 1], s=s, c=c, lw=(0 if not ec else 0.6 if ec == C_MY_EDGE else 1.6),
                       edgecolors=(ec if ec else 'none'), zorder=z, rasterized=True)
        map_frame(ax, cond)
        scale_bar(ax)
        map_bot[cond] = ax
    h6 = [Line2D([], [], marker='o', ls='', color=C_STR, ms=9, label='Fibroblast-Rich Stroma'),
          Line2D([], [], marker='o', ls='', color=C_EPI6, ms=9, label='Epithelium (EpCAM)'),
          Line2D([], [], marker='o', ls='', color=C_HLA, ms=11, label='HLA-DR-High Myeloid'),
          Line2D([], [], marker='o', ls='', color=C_INOS, ms=11, label='iNOS-High Myeloid'),
          Line2D([], [], marker='o', ls='', color=C_INOS, mec=C_HLA, mew=2, ms=11, label='HLA-DR-High and iNOS-High'),
          Line2D([], [], marker='o', ls='', color=C_MY, mec=C_MY_EDGE, mew=0.8, ms=10, label='Other Myeloid'),
          Line2D([], [], marker='o', ls='', color='#BDBDBD', ms=4, label='Other Cells')]
    lax_b = fig.add_axes([0, 0, 0.1, 0.1])
    put(lax_b, 1, 1, 1, MAP_H_IN)
    lax_b.axis('off')
    lax_b.legend(handles=h6, loc='center left', bbox_to_anchor=(0, 0.5), frameon=False, fontsize=FS_MLEG,
                 handletextpad=0.3, borderaxespad=0, borderpad=0,
                 title='Myeloid = Macrophage or DC\n(High = Above p75 of Section)', title_fontsize=FS_MLEG,
                 alignment='left')

    # ── top-right: depth profiles ──
    PR, FNAME, BIN_LAB, PLOT_FEATS = V_CT['PR'], V_CT['FNAME'], V_CT['BIN_LAB'], V_CT['PLOT_FEATS']
    FCOL = {'Stroma': C_STR, 'APC': '#555555', 'Treg': C_TREG, 'Other T': '#5AAE61'}
    PTITLE = {'Stroma': 'Fibroblast-Rich Stroma', 'APC': 'All APCs (DC + Macrophage + B)',
              'Treg': 'Tregs (FOXP3⁺CTLA4⁺)', 'Other T': 'Other T Cells (Control)'}
    xs = np.arange(len(BIN_LAB))
    LS, MK = {'UC': '-', 'CD': '--'}, {'UC': 'o', 'CD': 's'}
    hatch_who, PAX = set(), []
    for i, kf in enumerate(PLOT_FEATS):
        ax = fig.add_axes([0, 0, 0.1, 0.1])
        put(ax, 20, 1, 5, 3)
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
        ax.set_title(PTITLE[kf], fontsize=FS_PTITLE, color=FCOL[kf], pad=5)
        ax.set_xticks(xs)
        ax.set_xticklabels(BIN_LAB if i >= 2 else [], rotation=35, ha='right')
        ax.tick_params(labelsize=FS_TICK)
        ax.spines[['top', 'right']].set_visible(False)
        PAX.append(ax)

    # ── bottom-right: forest plot ──
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
    ARROW_PT = 54            # 'Higher in UC / CD' arrows, 54 pt below the axis
    forest_axes = []
    for i, s in enumerate(FOREST_SUB):
        ax = fig.add_axes([0, 0, 0.1, 0.1])
        put(ax, 30, 5, FOREST_W, 3)
        forest_axes.append(ax)
        sub = FR[FR.subset == s].set_index('protein')
        ax.axvline(0, color='#888888', lw=1.2, zorder=0)
        qkw = dict(xycoords=ax.get_yaxis_transform(), xytext=(Q_OFF_PT, 0), textcoords='offset points', ha='left',
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
                    # 2-line header nudged 3 pt up (into the spacer above it) so it clears the lead marker below
                    ax.annotate(FAMDISP[lab], (0, yy), xycoords=ax.get_yaxis_transform(),
                                xytext=(-7, 3 if '\n' in FAMDISP[lab] else 0),
                                textcoords='offset points', ha='right', va='center', fontsize=16,
                                fontweight='bold', linespacing=0.95)
            elif i == 0:
                ax.annotate(DISP.get(lab, lab), (0, yy), xycoords=ax.get_yaxis_transform(), xytext=(-16, 0),
                            textcoords='offset points', ha='right', va='center', fontsize=16,
                            fontweight=('bold' if lab in LEAD else 'normal'))
                if lab in DAGGER:
                    ax.annotate('†', (0, yy), xycoords=ax.get_yaxis_transform(), xytext=(-5, 0),
                                textcoords='offset points', ha='right', va='center', fontsize=16, color='#888888')
        ax.set_ylim(ymax - 0.5, -0.5)
        ax.set_yticks([])
        ax.set_xlim(-lim, lim)
        ax.set_xticks(np.arange(-lim + 1, lim, 2))
        ax.set_xlabel("Cohen's d, CD − UC (95% CI)")
        ax.set_title(STITLE[s], fontsize=18)
        ax.annotate('← Higher in UC', (0.02, 0), xycoords='axes fraction', xytext=(0, -ARROW_PT),
                    textcoords='offset points', ha='left', va='top', fontsize=15, color=COL['UC'])
        ax.annotate('Higher in CD →', (0.98, 0), xycoords='axes fraction', xytext=(0, -ARROW_PT),
                    textcoords='offset points', ha='right', va='top', fontsize=15, color=COL['CD'])
        ax.spines[['top', 'right', 'left']].set_visible(False)

    # ── bottom-right: Treg box plots ──
    POS = {('UC', 'conv'): 0.0, ('UC', 'Treg'): 1.0, ('CD', 'conv'): 2.5, ('CD', 'Treg'): 3.5}
    TQ = TC.set_index(['family', 'protein'])
    GROUP_PT = 50            # 'UC' / 'CD' group labels, 50 pt below the axis
    treg_axes = []
    for i, p in enumerate(TREG_MARKERS):
        ax = fig.add_axes([0, 0, 0.1, 0.1])
        put(ax, 30, 15, 4, 3)
        treg_axes.append(ax)
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
            ax.annotate(cond, (xm, 0), xycoords=ax.get_xaxis_transform(), xytext=(0, -GROUP_PT),
                        textcoords='offset points', ha='center', va='top', fontsize=17, color=COL[cond],
                        fontweight='bold')
        ax.set_ylabel(f'{p}\n(Mean per FOV, a.u.)')
        ax.spines[['top', 'right']].set_visible(False)

    # ════════════════════════ packing: measure overhangs, then place in inches ════════════════════════
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    OV = {ax: over(ax) for ax in BOX}
    leg_w = {k: a.get_legend().get_window_extent(rend).width / fig.dpi for k, a in [('t', lax_t), ('b', lax_b)]}
    leg_h = {k: a.get_legend().get_window_extent(rend).height / fig.dpi for k, a in [('t', lax_t), ('b', lax_b)]}
    # text sizes of the shared depth-profile labels and key (measured on a scratch copy)
    _yl = fig.text(0.5, 0.5, '% of Cells in Distance Bin\n(Mean ± SE Across FOVs)', rotation=90, fontsize=FS_LAB)
    _xl = fig.text(0.5, 0.5, 'Distance to Nearest Epithelial Cell (µm)', fontsize=FS_LAB)
    who = ' and '.join(sorted(hatch_who)) if hatch_who else 'None'
    lh = [Line2D([], [], color='black', ls='-', marker='o', ms=8, lw=2.2, label='UC'),
          Line2D([], [], color='black', ls='--', marker='s', mfc='white', ms=8, lw=2.2, label='CD'),
          Patch(facecolor='none', edgecolor='#B0B0B0', hatch='///', label=f'{who} < 3 FOVs (Not Drawn)')]
    _key = fig.legend(handles=lh, loc='center', ncol=3, frameon=False, fontsize=FS_KEY, borderaxespad=0, borderpad=0,
                      columnspacing=1.6)
    fig.canvas.draw()
    YL_W = _yl.get_window_extent(rend).width / fig.dpi
    XL_H = _xl.get_window_extent(rend).height / fig.dpi
    KEY_H = _key.get_window_extent(rend).height / fig.dpi
    _yl.remove(); _xl.remove(); _key.remove()

    # ---- left column ----
    x_cd = MAP_W_IN['UC'] + MAP_GAP_IN
    x_leg = x_cd + MAP_W_IN['CD'] + LEG_GAP_IN
    LEFT_W = x_leg + max(leg_w.values())
    Tm = max(OV[map_top[c]][1] for c in map_top)                 # 'UC' / 'CD' titles
    y_map1 = Tm
    y_map2 = y_map1 + MAP_H_IN + MAP_ROW_GAP_IN
    Bm = max(0, max(OV[map_bot[c]][3] for c in map_bot))
    H_L = y_map2 + MAP_H_IN + Bm                                 # left column height (figure height)

    # ---- right column widths: forest row sets the width ----
    XR = LEFT_W + COL_GAP_IN
    Lf, Rm, Rd = OV[forest_axes[0]][0], OV[forest_axes[0]][2], OV[forest_axes[1]][2]
    W_R = Lf + FOREST_W + Rm + FOREST_GAP_IN + FOREST_W + Rd
    Lt = [OV[a][0] for a in treg_axes]
    Rt = [max(OV[a][2], 0) for a in treg_axes]
    tw = (W_R - sum(Lt) - sum(Rt) - 2 * TREG_GAP_IN) / 3
    Llp = [OV[a][0] for a in PAX]                                # y tick labels
    Rlp = [max(OV[a][2], 0) for a in PAX]
    lp_x0 = XR + YL_W + 0.1 + max(Llp[0], Llp[2])
    lp_w = (XR + W_R - max(Rlp[1], Rlp[3]) - lp_x0 - LP_COL_GAP_IN - max(Llp[1], Llp[3])) / 2
    lp_x1 = lp_x0 + lp_w + LP_COL_GAP_IN + max(Llp[1], Llp[3])

    # ---- right column heights: top = left top, bottom = left bottom ----
    Tlp = max(OV[a][1] for a in PAX)                             # profile titles
    Blp = max(OV[a][3] for a in PAX[2:])                         # rotated x tick labels
    B1lp = max(max(OV[a][3] for a in PAX[:2]), 0)               # row 1: tick marks / half of the lowest y label
    Tf, Bf = max(OV[a][1] for a in forest_axes), max(OV[a][3] for a in forest_axes)
    Tt, Bt = max(OV[a][1] for a in treg_axes), max(OV[a][3] for a in treg_axes)
    TOP_STATIC = KEY_H + 0.12 + Tlp + B1lp + LP_ROW_GAP_IN + Tlp + Blp + 0.06 + XL_H
    BOT_STATIC = Tf + Bf + ROW_GAP_IN + Tt + Bt
    # boundary: bottom-right block (forest title top) aligned with the top of the bottom map row
    y_bound = y_map2
    lp_h = (y_bound - RIGHT_GAP_IN - TOP_STATIC) / 2
    ft_h = H_L - y_bound - BOT_STATIC
    F_H, T_H = FOREST_FRAC * ft_h, (1 - FOREST_FRAC) * ft_h
    print(f'natural fit: profile axis h {lp_h:.2f} in, forest h {F_H:.2f} in, Treg h {T_H:.2f} in')
    # minimum readable heights; if the bottom-right block is short, move the boundary up (profiles shrink)
    F_MIN, T_MIN, LP_MIN = 2.8, 2.45, 2.35
    if F_H < F_MIN or T_H < T_MIN:
        need = max(F_MIN - F_H, 0) / FOREST_FRAC if F_H / FOREST_FRAC < T_H / (1 - FOREST_FRAC) \
            else max(T_MIN - T_H, 0) / (1 - FOREST_FRAC)
        need = max(need, (F_MIN - F_H) / FOREST_FRAC, (T_MIN - T_H) / (1 - FOREST_FRAC))
        shift = min(need, 2 * max(lp_h - LP_MIN, 0))
        y_bound -= shift
        lp_h -= shift / 2
        ft_h += shift
        F_H, T_H = FOREST_FRAC * ft_h, (1 - FOREST_FRAC) * ft_h
    print(f'used: boundary {y_map2 - y_bound:.2f} in above the bottom map row; profile axis h {lp_h:.2f}, '
      f'forest h {F_H:.2f}, Treg h {T_H:.2f} in')

    FW, FH = XR + W_R, H_L
    fig.set_size_inches(FW, FH)
    # left column
    for maps, y0 in [(map_top, y_map1), (map_bot, y_map2)]:
        put(maps['UC'], 0, y0, MAP_W_IN['UC'], MAP_H_IN)
        put(maps['CD'], x_cd, y0, MAP_W_IN['CD'], MAP_H_IN)
    put(lax_t, x_leg, y_map1, leg_w['t'], MAP_H_IN)
    put(lax_b, x_leg, y_map2, leg_w['b'], MAP_H_IN)
    # top-right
    y_lp1 = KEY_H + 0.12 + Tlp
    y_lp2 = y_lp1 + lp_h + B1lp + LP_ROW_GAP_IN + Tlp
    for i, ax in enumerate(PAX):
        put(ax, (lp_x0, lp_x1)[i % 2], (y_lp1, y_lp2)[i // 2], lp_w, lp_h)
    # bottom-right
    y_for = y_bound + Tf
    y_treg = y_for + F_H + Bf + ROW_GAP_IN + Tt
    put(forest_axes[0], XR + Lf, y_for, FOREST_W, F_H)
    put(forest_axes[1], XR + Lf + FOREST_W + Rm + FOREST_GAP_IN, y_for, FOREST_W, F_H)
    x_ = XR
    for a, l_, r_ in zip(treg_axes, Lt, Rt):
        put(a, x_ + l_, y_treg, tw, T_H)
        x_ += l_ + tw + r_ + TREG_GAP_IN
    # shared depth-profile labels and UC / CD key (figure fractions from inches)
    lp_mid_x = (lp_x0 + lp_x1 + lp_w) / 2
    YLAB = fig.text(XR / FW, 1 - (y_lp1 + y_lp2 + lp_h) / 2 / FH,
                    '% of Cells in Distance Bin\n(Mean ± SE Across FOVs)', rotation=90, ha='left', va='center',
                    fontsize=FS_LAB, multialignment='center')
    XLAB = fig.text(lp_mid_x / FW, 1 - (y_lp2 + lp_h + Blp + 0.06) / FH, 'Distance to Nearest Epithelial Cell (µm)',
                    ha='center', va='top', fontsize=FS_LAB)
    KEY = fig.legend(handles=lh, loc='upper center', bbox_to_anchor=(lp_mid_x / FW, 1.0), ncol=3, frameon=False,
                     fontsize=FS_KEY, borderaxespad=0, borderpad=0, columnspacing=1.6)

    # ════════════════════════ layout checks ════════════════════════
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()


    def bbin(a):
        """tight bbox in inches, y measured from the top."""
        bb = a.get_tightbbox(rend) if hasattr(a, 'get_tightbbox') and not isinstance(a, plt.Text) \
            else a.get_window_extent(rend)
        return (bb.x0 / fig.dpi, FH - bb.y1 / fig.dpi, bb.x1 / fig.dpi, FH - bb.y0 / fig.dpi)


    def union(items):
        b = [bbin(a) for a in items]
        return (min(x[0] for x in b), min(x[1] for x in b), max(x[2] for x in b), max(x[3] for x in b))


    BL = {'TL maps': union(list(map_top.values())), 'TL legend': bbin(lax_t.get_legend()),
          'BL maps': union(list(map_bot.values())), 'BL legend': bbin(lax_b.get_legend()),
          'TR': union(PAX + [YLAB, XLAB, KEY]), 'BR forest': union(forest_axes), 'BR treg': union(treg_axes)}
    for k, v in BL.items():
        print(f'{k:10s} x {v[0]:6.2f}-{v[2]:6.2f}  y(top-down) {v[1]:6.2f}-{v[3]:6.2f} in')
    P = {k: a.get_position() for k, a in [('UCt', map_top['UC']), ('UCb', map_bot['UC']),
                                          ('CDt', map_top['CD']), ('CDb', map_bot['CD'])]}
    assert abs(P['UCt'].x0 - P['UCb'].x0) < 1e-9 and abs(P['UCt'].x1 - P['UCb'].x1) < 1e-9, 'UC maps not aligned'
    assert abs(P['CDt'].x0 - P['CDb'].x0) < 1e-9 and abs(P['CDt'].x1 - P['CDb'].x1) < 1e-9, 'CD maps not aligned'
    for a in list(map_top.values()) + list(map_bot.values()):           # same µm per inch (equal aspect, same box)
        pos = a.get_position()
        assert abs(pos.height * FH - MAP_H_IN) < 1e-6
        assert abs((a.get_ylim()[0] - a.get_ylim()[1]) / (pos.height * FH) - UM_PER_IN) < 1e-6
        assert abs((a.get_xlim()[1] - a.get_xlim()[0]) / (pos.width * FW) - UM_PER_IN) < 1e-6
    assert abs(BL['TL legend'][0] - BL['BL legend'][0]) < 0.02, 'legends do not start at the same x'
    assert BL['TL legend'][0] > BL['TL maps'][2] and BL['BL legend'][0] > BL['BL maps'][2], 'legend over map'
    assert max(BL['TL legend'][2], BL['BL legend'][2]) < BL['TR'][0] and \
        max(BL['TL legend'][2], BL['BL legend'][2]) < BL['BR forest'][0], 'left column runs into right column'
    assert BL['TL maps'][3] < BL['BL maps'][1] + 0.05, 'map rows overlap'
    assert BL['TR'][3] < BL['BR forest'][1] and BL['BR forest'][3] < BL['BR treg'][1], 'right rows overlap'
    fe, te = BL['BR forest'], BL['BR treg']
    assert abs(fe[0] - te[0]) < 0.02 and abs(fe[2] - te[2]) < 0.02, ('forest / Treg rows not equal width', fe, te)
    _ptb = [bbin(a) for a in PAX]
    assert _ptb[0][3] < _ptb[2][1] and _ptb[1][3] < _ptb[3][1], 'profile rows overlap'
    assert _ptb[0][2] < _ptb[1][0] and _ptb[2][2] < _ptb[3][0], 'profile columns overlap'
    assert bbin(YLAB)[2] < min(_ptb[0][0], _ptb[2][0]), 'y label overlaps tick labels'
    assert bbin(XLAB)[1] > max(_ptb[2][3], _ptb[3][3]), 'x label overlaps tick labels'
    assert bbin(KEY)[3] < min(_ptb[0][1], _ptb[1][1]), 'key overlaps profile titles'
    print(f"top edges: left {BL['TL maps'][1]:.2f}, right {BL['TR'][1]:.2f}; bottom edges: left "
      f"{BL['BL maps'][3]:.2f}, right {BL['BR treg'][3]:.2f}; boundary: BL maps top {P['UCb'].y1:.3f} "
      f"(fig frac) vs forest top {1 - BL['BR forest'][1] / FH:.3f}")
    print(f'fig {FW:.2f} x {FH:.2f} in; map h {MAP_H_IN:.3f} in ({UM_PER_IN:.1f} µm/in); maps UC {MAP_W_IN["UC"]:.2f} '
      f'+ CD {MAP_W_IN["CD"]:.2f} in; legends {leg_w["t"]:.2f} / {leg_w["b"]:.2f} in (h {leg_h["t"]:.2f} / '
      f'{leg_h["b"]:.2f}); right col {W_R:.2f} in; Treg axis w {tw:.2f}; profile axis w {lp_w:.2f}')
    with open(PNG, 'wb') as fh:
        tb = fig.get_tightbbox(rend)
        fig.savefig(fh, format='png', dpi=200, pad_inches=0,
                    bbox_inches=Bbox.from_extents(tb.x0 - 0.08, tb.y0 - 0.08, tb.x1 + 0.12, tb.y1 + 0.08))
    plt.close(fig)
    print('Done ->', PNG)
