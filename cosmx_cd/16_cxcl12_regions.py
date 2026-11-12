#!/usr/bin/env python
"""CD CosMx 6K resections: CXCL12-High Regions, APC-rich tissue and the Treg Suppression Score.

Panels made
  Fig 5E   section 43563_A4: OGN+RSPO3+ fibroblast density, APCs and Treg CXCR4 over the CXCL12 field (CXCL12-High Region outline),
           and the Treg Suppression Score (20-Treg means)
  Fig 5F   one tissue piece of section 43563_A8: cellular neighbourhoods and Treg Suppression Score
  Supp 9E  CXCL12 by cell type (% expressing, level in expressing cells)
  Supp 9F  Treg Suppression Score by quintile of the local APC fraction in CXCL12-Low vs CXCL12-High Regions, and its rise
  also     APC abundance per cellular neighbourhood (table read by 17_neighbourhood_apc_mhcii.py, Supp 9C; S9E_cxcl12_by_celltype with RUN_FIG2=1)

Inputs
  /share/fsmresfiles/UC/AtoMx/CD_surgical_resection/ (combined/all_samples.parquet, combined/raw_counts_strict_tregs.csv, per-section Processed/norm_log1p.h5ad)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/
    D_shared_inputs/treg_definition_comparison_per_cell.csv (10_treg_suppression_score.py)
    B_resection_analyses/cellular_neighborhoods_v2/ and ligand_receptor_state/ (11_neighbourhoods_local_tables.py, 12_neighbourhood_treg_score.py)
    ligand_receptor_state/c_top_apc_quintile_per_section_exp.csv (per-section top-APC-quintile table; not written by a script of this repository)
  /share/fsmresfiles/UC/revision_jenny/updated_figures/cd_spatial_results_all.csv (15_apc_local_models.py, section spatial_results_table)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/A_fig4_panels/ (per-Treg CXCL12 region table)
  /share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/ and /share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/ (panels, CSVs)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/archive/supp_apc_abundance_by_neighbourhood.csv (APC table, RUN_FIG2=1)

Run order of sections
  F5E_cxcl12_region_tregs -> F5EF_section_maps -> S9F_apc_cxcl12_quintiles -> S9E_cxcl12_by_celltype
  F5EF_section_maps reads the per-Treg table of F5E_cxcl12_region_tregs (score colour limits). It runs once per figure:
  Fig 5E with T94_SEC=43563_A4 T94_MAKE=row, Fig 5F with T94_SEC=43563_A8 T94_MAKE=piece. Its shared map code (prep(), figure helpers,
  colour scales, panel drawing) is held in SEG_* string constants inside the section and executed there.

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its upstream
  outputs exist. Command-line arguments and environment variables used by a section are listed in its own docstring.
  Paths are the lab's absolute paths (to be replaced by placeholders before public release).
  Environment: see the README section "Environment" (OPENBLAS_CORETYPE=SkylakeX on Intel Sapphire Rapids nodes;
  F5EF_section_maps asserts OPENBLAS_CORETYPE=SkylakeX and OPENBLAS_NUM_THREADS=32).
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
# Section F5E_cxcl12_region_tregs: per-Treg CXCL12 / CXCR4 / Treg Suppression Score table of section 54033_B4 and a companion CXCL12 region figure
# Panels: Fig 5E-F (Treg Suppression Score colour limits, read by F5EF_section_maps)
# ======================================================================================================================
if _run('F5E_cxcl12_region_tregs'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/67e_cxcl12_region_treg_response_exp.py'   # working folder of this step (relative paths resolve from here)
    """Per-Treg table and companion CXCL12 region figure of CD resection section 54033_B4 (CosMx 6K).
Local CXCL12 of a Treg (lig_CXCL12_exp, 11_neighbourhoods_local_tables.py) = mean within-section CXCL12 z over the stromal cell types expected to
express CXCL12 (OGN+RSPO3+ Fib, ADAMDEC1+ Fib, Activated Fib, Crypt Top Fib, PLA2G2A+ ECM Fib, Myofibroblast, Pericyte, Endothelial; no SELENOP+
cells, immune cells or epithelium) within 50 um, >= 3 required. act_z is the Treg Suppression Score (-D6 signature score, z within section).
Table: every strict Treg of the section with rec_CXCR4, act_z, lig_CXCL12, lig_CXCL12_exp, apc_frac and the means of its 20 nearest Tregs
  ('CXCR4 smoothed', 'suppression score smoothed'); display coordinates. F5EF_section_maps sets its score colour limits from this table.
Figure (companion, not a paper panel): tissue pieces moved closer together for display (translations only; distances within a piece unchanged).
  Row 1: Gaussian-smoothed (sigma 50 um, 10 um bins) % of cells expressing CXCL12, OGN+RSPO3+ fibroblast density and APC density (DC, SELENOP+ Mac,
  CD9+ M2 Mac, Mo-Mac, B cell), each with the CXCL12-High outline (upper 30% of the smoothed CXCL12 field, closed and hole-filled).
  Row 2: Treg CXCR4 and Treg Suppression Score (20-Treg means), then the mean score in CXCL12-High vs CXCL12-Low Regions (upper vs lower
  within-section tertile of lig_CXCL12_exp) per within-section quintile of the local APC fraction, all 11 sections (tables of
  11_neighbourhoods_local_tables.py); in-panel P = region x APC quintile interaction; stars = BH q over the five quintiles.
Inputs: combined/all_samples.parquet, the section's norm_log1p.h5ad, D_shared_inputs/treg_definition_comparison_per_cell.csv,
  ligand_receptor_state/ (treg_local_ligand_receptor.csv.gz, treg_local_cxcl12_producer_exp.csv.gz, c_*_exp.csv, a_cxcr4_cxcl12_exp.csv,
  sensitivity_no_producer_exp.csv), cellular_neighborhoods_v2/cn_per_cell.parquet.
Outputs (A_fig4_panels/): fig4_cxcl12_region_tregs_exp.csv; fig4_cxcl12_region_treg_response_exp.png / .pdf + _legend.txt.
PANEL2=ogn (default and only accepted value): the second map shows OGN+RSPO3+ fibroblast density."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from anndata.io import read_elem
    from scipy import sparse
    from scipy.spatial import cKDTree
    from scipy.ndimage import gaussian_filter
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/A_fig4_panels'
    LR = f'{P4}/B_resection_analyses/ligand_receptor_state'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028
    SEC, FOLDER = '54033_B4', 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133'
    APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']
    BIN, SIG, KNN = 10, 5.0, 20
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa = sa[sa.section == SEC]; sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
    with h5py.File(f'{R}/{FOLDER}/Processed/norm_log1p.h5ad', 'r') as h:
        ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
    nm = pd.Index(SEC + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    keep = nm.isin(sa.index); gl = [g for g in ['CXCL12', 'CXCR4'] if g in v]
    E = pd.DataFrame(X[np.where(keep)[0]][:, [v.get_loc(g) for g in gl]].toarray(), index=nm[keep], columns=gl); del X; gc.collect()
    cells = sa.loc[E.index].copy(); print('cells', len(cells), '| genes', gl, flush=True)
    # tissue pieces squeezed and stitched as for the CN map of 12_neighbourhood_treg_score.py (display only; distances within a piece are unchanged)
    GAP, SBIN = 60, 50
    def squeeze(v):
        b = np.floor((v - v.min()) / SBIN).astype(int); occ = np.zeros(b.max() + 2, bool); occ[np.unique(b)] = True
        shift = np.zeros(len(occ)); run = 0; cut = 0
        for i in range(len(occ)):
            if not occ[i]:
                run += 1
                if run * SBIN > GAP: cut += SBIN
            else: run = 0
            shift[i] = cut
        return lambda x: x - shift[np.clip(np.floor((np.asarray(x) - v.min()) / SBIN).astype(int), 0, len(shift) - 1)]
    fx, fy = squeeze(cells.x.values), squeeze(cells.y.values)
    PCn = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0)
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0); T = T[T.section == SEC].copy()
    T = T.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left')   # lig_CXCL12_exp, n_exp_producers_50um (11_neighbourhoods_local_tables.py)
    T = T.join(sa[['x', 'y']], how='inner', rsuffix='_s')
    cells['x'], cells['y'] = fx(cells.x.values), fy(cells.y.values); T['x'], T['y'] = fx(T.x.values), fy(T.y.values)
    from scipy import ndimage
    gx = np.floor((cells.x - cells.x.min()) / SBIN).astype(int).values; gy = np.floor((cells.y - cells.y.min()) / SBIN).astype(int).values
    occ = np.zeros((gy.max() + 1, gx.max() + 1), bool); occ[gy, gx] = True
    lab, npc = ndimage.label(ndimage.binary_dilation(occ, iterations=1), structure=np.ones((3, 3)))
    cells['piece'] = lab[gy, gx]; order = list(cells.piece.value_counts().index)
    SHIFT = {order[0]: np.zeros(2)}; placed = cells[cells.piece == order[0]][['x', 'y']].values
    for pc in order[1:]:
        P = cells[cells.piece == pc][['x', 'y']].values; tree = cKDTree(placed[::5]); c0 = placed.mean(0)
        d = c0 - P.mean(0); Lh = np.linalg.norm(d); u = d / Lh; best = np.zeros(2); t = 0.0
        while t < Lh:
            t2 = t + 25; sh = u * t2
            if tree.query(P[::5] + sh, k=1)[0].min() < GAP: break
            best, t = sh, t2
        SHIFT[pc] = best; placed = np.vstack([placed, P + best])
    print('tissue pieces', npc, '| shifts um', {int(k): np.round(v).tolist() for k, v in SHIFT.items()}, flush=True)
    sh = np.vstack([SHIFT[p] for p in cells.piece.values]); cells['x'] += sh[:, 0]; cells['y'] += sh[:, 1]
    tp = cells.piece.values[cKDTree(cells[['x', 'y']].values - sh).query(T[['x', 'y']].values)[1]]
    shT = np.vstack([SHIFT[p] for p in tp]); T['x'] += shT[:, 0]; T['y'] += shT[:, 1]
    x0, y0 = cells.x.min() - 3 * SIG * BIN, cells.y.min() - 3 * SIG * BIN
    nx, ny = int((cells.x.max() - x0) / BIN) + 4 * int(SIG) + 4, int((cells.y.max() - y0) / BIN) + 4 * int(SIG) + 4
    ix, iy = ((cells.x - x0) / BIN).astype(int).values, ((cells.y - y0) / BIN).astype(int).values
    allc = np.zeros((ny, nx)); np.add.at(allc, (iy, ix), 1); Asm = gaussian_filter(allc, SIG); mask = Asm > 0.05 * Asm.max()
    mask = ndimage.binary_fill_holes(mask)   # small cell-free gaps inside the tissue (lumen, tears)
    # would leave the field undefined (white specks, each with its own contour loop); filled, the white outline marks only the region boundary.
    XC, YC = x0 + (np.arange(nx) + 0.5) * BIN, y0 + (np.arange(ny) + 0.5) * BIN
    EXT = (x0, x0 + nx * BIN, y0 + ny * BIN, y0)          # imshow origin='upper': row 0 sits at y0, so the axis is already y-down
    def field(vals):
        num = np.zeros((ny, nx)); np.add.at(num, (iy, ix), vals.astype(float))
        return np.where(mask, gaussian_filter(num, SIG) / np.where(Asm > 0, Asm, 1) * 100, np.nan)
    CNp = pq.read_table(f'{P4}/B_resection_analyses/cellular_neighborhoods_v2/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas()
    cells['CN'] = CNp.CN_k10.reindex(cells.index).values
    print('CN10 cells in this section: {:,} / {:,}'.format(int((cells.CN == 10).sum()), len(cells)), flush=True)
    P2 = os.environ.get('PANEL2', 'ogn'); assert P2 == 'ogn', 'run with PANEL2=ogn'   # OGN+RSPO3+ density is shown directly: no clustering step between
    # the panel and the claim. The CN10 enrichment inside the CXCL12-High outline is printed below for reference.
    F = {'CXCL12': field(E.CXCL12.values > 0), 'APC': field(cells.ct.isin(APC).values),
         'P2': field((cells.CN == 10).values if P2 == 'cn10' else (cells.ct == 'OGN+RSPO3+ Fib').values)}
    thr = np.nanpercentile(F['CXCL12'], 70)
    CXhi = ndimage.binary_closing(np.where(np.isnan(F['CXCL12']), 0.0, F['CXCL12'] >= thr).astype(bool), np.ones((3, 3)), iterations=2)
    CXhi = ndimage.binary_fill_holes(CXhi).astype(float)
    inside = CXhi[iy, ix] > 0
    OGNin = (cells.ct[inside] == 'OGN+RSPO3+ Fib').mean() * 100; OGNout = (cells.ct[~inside] == 'OGN+RSPO3+ Fib').mean() * 100
    CNin = (cells.CN[inside] == 10).mean() * 100; CNout = (cells.CN[~inside] == 10).mean() * 100
    APCin = cells.ct[inside].isin(APC).mean() * 100; APCout = cells.ct[~inside].isin(APC).mean() * 100
    print('CXCL12-high contour at {:.1f}% expressing (70th pct)'.format(thr), flush=True)
    print('  inside vs outside: CN10 {:.1f}% vs {:.1f}% ({:.1f}x) | OGN+RSPO3+ Fib {:.2f}% vs {:.2f}% ({:.1f}x) | APC {:.1f}% vs {:.1f}% ({:.2f}x) | cells {:,} vs {:,}'.format(
        CNin, CNout, CNin / CNout, OGNin, OGNout, OGNin / OGNout, APCin, APCout, APCin / APCout, int(inside.sum()), int((~inside).sum())), flush=True)
    tree = cKDTree(T[['x', 'y']].values); _, nb = tree.query(T[['x', 'y']].values, k=min(KNN, len(T)))
    for src, dst in [('rec_CXCR4', 'CXCR4 smoothed'), ('act_z', 'suppression score smoothed')]:
        vv = T[src].values; T[dst] = np.nanmean(np.where(np.isnan(vv[nb]), np.nan, vv[nb]), axis=1)
    with open(f'{OUT}/fig4_cxcl12_region_tregs_exp.csv', 'w') as fh:
        T[['x', 'y', 'rec_CXCR4', 'act_z', 'lig_CXCL12', 'lig_CXCL12_exp', 'n_exp_producers_50um', 'apc_frac', 'CXCR4 smoothed', 'suppression score smoothed']].rename(
            columns={'act_z': 'treg_suppression_score_z'}).to_csv(fh, float_format='%.4g')
    print('Tregs plotted:', len(T), flush=True)
    C = pd.read_csv(f'{LR}/c_apc_matched_cxcl12_high_vs_low_exp.csv').rename(columns={'suppression_cxcl12_high': 'act_cxcl12_high', 'suppression_cxcl12_low': 'act_cxcl12_low'})
    PL = pd.read_csv(f'{LR}/c_pooled_apc_stratified_exp.csv')
    IX = pd.read_csv(f'{LR}/c_interaction_apc_by_cxcl12_countadj_exp.csv').set_index('term')
    AP = pd.read_csv(f'{LR}/a_cxcr4_cxcl12_exp.csv'); AP = AP[AP.cxcl12 == 'lig_CXCL12_exp'].set_index('model')
    SEN = pd.read_csv(f'{LR}/sensitivity_no_producer_exp.csv')
    S2 = SEN[SEN.analysis.str.startswith('S2 interaction')].set_index('group')
    NP = pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0); N_NOPROD, N_ALL = int(NP.lig_CXCL12_exp.isna().sum()), len(NP)
    def pfmt(p, tex=True):
        if p >= 0.001: return f'{p:.2g}'
        ex = int(np.floor(np.log10(p))); m = p / 10 ** ex
        if m >= 9.95: m, ex = 1.0, ex + 1
        return '{:.1f} \u00d7 10$^{{{}}}$'.format(m, ex) if tex else '{:.1f}\u00d710{}'.format(m, str(ex).translate(str.maketrans('-0123456789', '\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079')))
    HI_SIG = IX.p['hi'] < 0.05 and IX.coef['hi'] < 0; IX_SIG = IX.p['hi:apcz'] < 0.05 and IX.coef['hi:apcz'] < 0   # the statistic this panel actually claims
    def stars(q): return 'n.s.' if q >= 0.05 else ('*' if q >= 0.01 else ('**' if q >= 0.001 else '***'))
    plt.rcParams.update({'font.size': 15})
    # 2 x 3 with generous gaps (six maps in a row make the subtitles collide; the two rows can be cropped apart)
    # explicit layout in inches: five maps of identical size (column width = row height x tissue aspect) with equal gaps;
    # row 2 adds a spacer for the bar chart's y-axis label/ticks, then the bar chart; row 1 is left-aligned with a filler
    FH = 12.2; TOP, BOT, HSP = 0.88, 0.11, 0.30
    RH = FH * (TOP - BOT) / (2 + HSP)                                  # row height (in)
    MW = RH * abs(EXT[1] - EXT[0]) / abs(EXT[3] - EXT[2])              # map width that exactly fits the tissue at aspect='equal'
    GAP, SPC, BW = 0.67, 0.95, 4.2                                     # map-map gap, map-to-bar spacer, bar chart width (in)
    ROW2 = 2 * MW + GAP + SPC + BW; FILL = ROW2 - (3 * MW + 2 * GAP)
    LM, RM = 0.55, 0.15; FW = LM + ROW2 + RM
    fig = plt.figure(figsize=(FW, FH))
    outer = fig.add_gridspec(2, 1, left=LM / FW, right=1 - RM / FW, top=TOP, bottom=BOT, hspace=HSP)
    g1 = outer[0].subgridspec(1, 6, wspace=0, width_ratios=[MW, GAP, MW, GAP, MW, FILL])
    g2 = outer[1].subgridspec(1, 5, wspace=0, width_ratios=[MW, GAP, MW, SPC, BW])
    axes = np.array([[fig.add_subplot(g1[0, c]) for c in (0, 2, 4)],
                     [fig.add_subplot(g2[0, c]) for c in (0, 2, 4)]])
    # the stitched tissue is ~0.5 as wide as tall, so the map axes are sized to it; the gap is then real spacing, not empty axes
    def cbar(ax, im, label, loc=(-0.12, 0.265), w=0.34):   # inside the panel, in white margin beside the tissue, clear of the row below
        cax = ax.inset_axes([loc[0], loc[1], w, 0.018])   # row 1 in the concave notch at the left; row 2 to the right of the tissue
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        cax.set_title(label, fontsize=11.5, pad=3)          # label ABOVE the bar, inside the panel, so it cannot reach the row below
        cb.ax.tick_params(labelsize=11, length=2.5, pad=1.5); cb.outline.set_linewidth(0.6)
        cb.ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        cb.ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '\u2212')))   # '0' not '0.0': keeps ticks apart
    def frame(ax):
        ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3]); ax.set_aspect('equal'); ax.axis('off')
    MAPS = [('CXCL12', 'CXCL12 Ligand', 'viridis', '% cells expressing'),
            ('P2', 'Neighbourhood CN10\n(OGN\u207aRSPO3\u207a/PLA2G2A\u207a Fib, Endothelial)' if P2 == 'cn10' else 'OGN\u207aRSPO3\u207a Fibroblasts', 'viridis', '% of cells'),
            ('APC', 'Antigen-Presenting Cells\n(DC, Macrophage, B Cell)', 'magma', '% of cells')]
    for j, (k, ttl, cmap, cl) in enumerate(MAPS):
        ax = axes[0, j]
        ax.imshow(np.where(mask, 0.0, np.nan), extent=EXT, cmap='Greys', vmin=0, vmax=6, interpolation='nearest')
        im = ax.imshow(F[k], extent=EXT, cmap=cmap, vmin=0, vmax=np.nanpercentile(F[k], 99), interpolation='bilinear')
        ax.contour(XC, YC, CXhi, levels=[0.5], colors='white', linewidths=1.3)
        ax.set_title(ttl, fontsize=15, pad=10); frame(ax); cbar(ax, im, cl)
    for j, (k, ttl, cmap) in enumerate([('CXCR4 smoothed', 'Treg CXCR4\n(The CXCL12 Receptor)', 'viridis'),
                                        ('suppression score smoothed', 'Treg Suppression Score', 'RdBu_r')]):
        ax = axes[1, j]; vv = T[k].values
        ax.scatter(cells.x, cells.y, s=0.25, c='#ebebeb', rasterized=True, lw=0)
        if cmap == 'RdBu_r':
            lim = np.nanpercentile(np.abs(vv - np.nanmean(vv)), 97); lo, hi = np.nanmean(vv) - lim, np.nanmean(vv) + lim
        else: lo, hi = np.nanpercentile(vv, 2), np.nanpercentile(vv, 98)
        sc = ax.scatter(T.x, T.y, c=vv, cmap=cmap, s=24, lw=0.2, edgecolor='0.35', vmin=lo, vmax=hi, rasterized=True)
        ax.contour(XC, YC, CXhi, levels=[0.5], colors='k', linewidths=1.3)
        ax.set_title(ttl, fontsize=15, pad=10); frame(ax); cbar(ax, sc, 'z (20-Treg mean)', loc=(-0.12, 0.265), w=0.46)
    ax = axes[1, 2]; w = 0.38; xs = np.arange(len(C))
    ax.bar(xs - w / 2, C.act_cxcl12_low, w, color='#2166ac', edgecolor='k', lw=0.7, label='CXCL12-Low Region')
    ax.bar(xs + w / 2, C.act_cxcl12_high, w, color='#b2182b', edgecolor='k', lw=0.7, label='CXCL12-High Region')
    ymax = max(C.act_cxcl12_low.max(), C.act_cxcl12_high.max()); ymin = min(C.act_cxcl12_low.min(), C.act_cxcl12_high.min())
    for xi, r in zip(xs, C.itertuples()):
        top = max(r.act_cxcl12_low, r.act_cxcl12_high, 0.0)   # label above the zero line too, so it never sits on a negative bar
        if r.q < 0.05:                                        # only significant quintiles are marked
            ax.text(xi, top + 0.02 * (ymax - ymin) + 0.012, stars(r.q), ha='center', va='bottom', fontsize=16)
    ax.axhline(0, color='0.3', lw=1); ax.set_xticks(xs); ax.set_xticklabels([f'{v*100:.0f}%' for v in C.mean_apc_frac], fontsize=13.5)
    ax.set_ylim(ymin - 0.08 * (ymax - ymin) - 0.03, ymax + 0.26 * (ymax - ymin))
    ax.set_xlabel('Antigen-Presenting Cells within 50 \u03bcm\n(Quintile of Local APC Fraction)', fontsize=14.5)
    ax.set_ylabel('Treg Suppression Score (z)', fontsize=15, labelpad=1)
    ax.set_title('Equal Sender Exposure\nCXCL12-High vs CXCL12-Low Regions', fontsize=16)
    ax.legend(frameon=False, fontsize=12.5, loc='upper left', handlelength=1.1, borderpad=0.1)
    pv = IX.p['hi:apcz']
    ax.text(0.03, 0.70, ('Weaker APC Response\nin CXCL12-High Regions\n' if IX_SIG else 'APC Response, CXCL12-High\nvs CXCL12-Low Regions\n') + '$P$ = ' + pfmt(pv) + ('' if pv < 0.05 else ' (n.s.)'),
            transform=ax.transAxes, ha='left', va='top', fontsize=11)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    TAG = '' if P2 == 'ogn' else '_CN10'
    STEM = f'{OUT}/fig4_cxcl12_region_treg_response_exp{TAG}'
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.06)
    fig.savefig(f'{STEM}.pdf', bbox_inches='tight', pad_inches=0.06); plt.close(fig)
    assert P2 == 'ogn'
    c5 = C.set_index('apc_quintile').loc[C.apc_quintile.max()]
    a2 = AP.loc['CXCR4 ~ CXCL12 + Treg Suppression Score + library size']
    s2_ok = S2.p['hi:apcz'] < 0.05 and S2.coef['hi:apcz'] < 0
    leg = ('Tregs in CXCL12-High Regions show a weaker response to local antigen-presenting cells (APCs). '
       'Maps, CD resection section {sec} (CosMx 6K, {ncell:,} cells; tissue pieces repositioned for display only): Gaussian-smoothed fraction of cells '
       'expressing CXCL12, OGN⁺RSPO3⁺ fibroblast density and APC density (DC, macrophage, B cell), and the {ntreg:,} strict Tregs coloured by CXCR4 '
       'and by the Treg Suppression Score (negated D6 pseudotime signature, z-scored within section; each Treg shown as the mean of its 20 nearest Tregs), '
       'with the upper 30% of the smoothed CXCL12 field outlined as a display aid. '
       'For all statistics (11 resection sections), local CXCL12 of each Treg was the mean within-section CXCL12 z-score of stromal cells expected to '
       'express CXCL12 (fibroblasts, myofibroblasts, pericytes and endothelial cells) within 50 μm, computed for the {nok:,} Tregs with at least three such '
       'cells ({nnp:,} Tregs with fewer were excluded{s2txt}); Treg CXCR4 rose with this local CXCL12 in a linear mixed model with a section random intercept '
       '(adjusted for the Treg Suppression Score and library size, P = {a2p}). '
       'Bar chart, mean Treg Suppression Score in CXCL12-High versus CXCL12-Low Regions (upper versus lower within-section tertile of local CXCL12; {nbar:,} Tregs) '
       'within quintiles of the local APC fraction, taken within section (x axis, mean fraction per quintile); in the highest APC quintile the score was '
       '{lo5:+.2f} SD in CXCL12-Low and {hi5:+.2f} SD in CXCL12-High Regions. '
       'In this top quintile the section mean was lower in CXCL12-High Regions in 10 of 11 sections. '  # 10 of 11: ligand_receptor_state/c_top_apc_quintile_per_section_exp.csv (difference_high_minus_low < 0)
       'The in-panel P is the CXCL12 region-by-APC quintile interaction from a library-size-adjusted linear mixed model with a section random intercept; '
       'stars, Benjamini-Hochberg q across the five quintiles (* < 0.05, ** < 0.01, *** < 0.001), shown for significant quintiles only.').format(
           sec=SEC, ncell=len(cells), ntreg=len(T), nok=N_ALL - N_NOPROD, nnp=N_NOPROD,
           s2txt=(', and counting them as CXCL12-Low did not change the interaction, P = ' + pfmt(S2.p['hi:apcz'], False)) if s2_ok else
                 (', and counting them as CXCL12-Low gave interaction P = ' + pfmt(S2.p['hi:apcz'], False)),
           a2p=pfmt(a2.p, False), nbar=int(PL.n[0]), lo5=c5.act_cxcl12_low, hi5=c5.act_cxcl12_high).replace('-0.', '\u22120.')
    for w in ['ctivation', 'erritor', 'lunting']: assert w not in leg, w
    with open(f'{STEM}_legend.txt', 'w') as fh: fh.write(leg + chr(10))
    print(chr(10) + leg)


# ======================================================================================================================
# Section F5EF_section_maps: maps of CD resection sections 43563_A4 (Fig 5E, four-panel row) and 43563_A8 (Fig 5F, one tissue piece)
# Panels: Fig 5E, 5F
# Note: Fig 5E: T94_SEC=43563_A4 T94_MAKE=row; Fig 5F: T94_SEC=43563_A8 T94_MAKE=piece (environment variables)
# ======================================================================================================================
if _run('F5EF_section_maps'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/94_plain_maps_A4_A8.py'   # working folder of this step (relative paths resolve from here)
    """Maps of CD resection sections (CosMx 6K; 11 sections, 4 patients). One run per section (T94_SEC, default 43563_A8) and output set (T94_MAKE,
comma-separated; defaults: 43563_A8 cn,ogn,row; 43563_A4 row; 54033_B4 / 43563_A11 fov).
  Fig 5E   T94_SEC=43563_A4 T94_MAKE=row -> fig4_A4_four_panel_plain.png: four maps of the same size, each with its colourbar and key on its right:
           OGN+RSPO3+ fibroblast density (% of cells; white CXCL12-High Region outline) | APCs (DC, macrophage, B cell) over the CXCL12 field
           (% of cells expressing; CXCL12-High Region border) | Treg CXCR4 (20-Treg mean) over the CXCL12 field | Treg Suppression Score (20-Treg mean).
  Fig 5F   T94_SEC=43563_A8 T94_MAKE=piece -> fig4_A8_cn_treg_two_panel_piece.png: one display unit of 43563_A8 alone (the largest unit whose centroid
           lies left of and below the section's cell median: the large lower-left piece with a CN8 aggregate), cellular neighbourhoods (CN_k10; full CN
           names from cn_k10_celltype_log2_enrichment.csv, key in two columns below) | Treg Suppression Score (90 um Treg dots over the non-Treg cells
           drawn as small light-grey points with the CN-panel marker size).
  Other T94_MAKE values write companion maps (not paper panels): cn (CN map of the 43563_A8 display crop), ogn (OGN+RSPO3+ density of the crop),
           row for 43563_A8, two / two_full / two_cropped (CN | score of the full section or the crop), piececx (APCs over the CXCL12 field for the
           Fig 5F piece); for 54033_B4 (FOVs 160-221) and 43563_A11 (FOVs 1-68, plus 193-229 only if in the same tissue fragment): fov (CN | score of
           one tissue piece in tissue coordinates) and, for 54033_B4, fovcx (APCs over the CXCL12 field of that piece).
Display: tissue pieces are moved closer together (translations only; distances within a piece unchanged). The 43563_A8 crop (T67U_V8C=1, set
  for 43563_A8) leaves out three pieces (the two bottom-right ones and the right-most middle one) and re-packs the other five.
Fields: per tissue piece, Gaussian-smoothed (sigma 50 um, 10 um grid) % of cells; CXCL12-High Region = upper 30% of the smoothed CXCL12 field,
  closed, hole-filled, fragments < 0.25 mm2 dropped (display only; statistics use the per-Treg definitions). Treg values = mean of the 20 nearest
  Tregs (itself included) in tissue coordinates. Score colour limits = mean +- the 97th percentile of the absolute deviation of the 20-Treg-mean
  score in fig4_cxcl12_region_tregs_exp.csv (section F5E_cxcl12_region_tregs). CXCL12 colour scale per section (median / 75th / 98th percentile).
Shared map code: prep() and its constants, the 43563_A8 re-pack, figure helpers, colour scales, panels 1-3 and a render check are held in the SEG_*
  string constants below and executed in two namespaces (prep_ns, helper_ns); prep() also recomputes the per-section top-APC-quintile numbers from
  the Treg tables and asserts them against c_top_apc_quintile_per_section_exp.csv.
Inputs: combined/all_samples.parquet, combined/raw_counts_strict_tregs.csv, each section's norm_log1p.h5ad; ligand_receptor_state/
  (treg_local_ligand_receptor.csv.gz, treg_local_cxcl12_producer_exp.csv.gz, c_*_exp.csv incl. c_top_apc_quintile_per_section_exp.csv);
  cellular_neighborhoods_v2/ (cn_per_cell.parquet, cn_treg_state_enrichment_all.csv, cn_k10_celltype_log2_enrichment.csv);
  A_fig4_panels/fig4_cxcl12_region_tregs_exp.csv.
Outputs: /share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/ (T94_OUTDIR=<folder> writes elsewhere); existing files are not overwritten.
Other environment variables: T94_EPI=ct|cn (epithelial cells of the two-tone background of the 43563_A11 piece: cell types or CN1 / CN6).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 T94_SEC=<section> T94_MAKE=<outputs> RUN_SECTIONS=F5EF_section_maps python <this file>
"""
    import os, tempfile, types, textwrap
    HERE = os.path.dirname(os.path.abspath(__file__))
    SEC = os.environ.get('T94_SEC', '43563_A8'); assert SEC in ('43563_A8', '43563_A4', '54033_B4', '43563_A11'), SEC; A8 = SEC == '43563_A8'
    FOVSEL = {'54033_B4': ([(160, 221)], 'B4'), '43563_A11': ([(1, 68), (193, 229)], 'A11')}   # pieces by FOV; later FOV blocks only if in the same tissue fragment
    MAKE = set(os.environ.get('T94_MAKE', 'cn,ogn,row' if A8 else ('fov' if SEC in FOVSEL else 'row')).split(','))
    assert MAKE <= {'cn', 'ogn', 'row', 'two', 'two_full', 'two_cropped', 'piece', 'piececx', 'fov', 'fovcx'} and (A8 or MAKE <= ({'fov', 'fovcx'} if SEC in FOVSEL else {'row', 'two'})), MAKE
    for k_ in ('T67U_V8_FULL', 'T67U_V8B', 'T67U_V8C', 'T67U_V8_PAIR', 'T67U_V7_CN4', 'T67U_CELLS', 'T67U_PCT', 'T67U_CN2', 'T67U_EPS'): os.environ.pop(k_, None)
    os.environ.update(T67U_VER='v8', T67U_SEC=SEC, T67U_V7_TMP=tempfile.gettempdir(), PANEL2='ogn', **({'T67U_V8C': '1'} if A8 else {}))   # settings of the shared map code (43563_A8: display crop)
    # ---- shared map code (string constants), executed below in the namespaces prep_ns and helper_ns
    SRC_U, SRC_T = 'F5EF_section_maps:prep', 'F5EF_section_maps:helpers'   # compile() labels (tracebacks)
    # constants, data, squeeze(), stitch(), prep()
    SEG_PREP = r'''"""Shared map code: constants, data, squeeze(), stitch() and prep() (per-piece fields, tissue masks, CXCL12-High Region, display units)."""
import os, gc, warnings; warnings.filterwarnings('ignore')
assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, PowerNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from PIL import Image
from anndata.io import read_elem
from scipy import sparse, stats, ndimage
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
OUT = f'{P4}/A_fig4_panels'
LR = f'{P4}/B_resection_analyses/ligand_receptor_state'
R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028
FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
          '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
          '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
          '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
          '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
          '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']
BIN, SIG = 10, 5.0                                                 # 10 um bins, sigma 5 bins = 50 um for the cell-composition fields
NW_SIG, MIN_W, CX_PCT, APC_PCT, VLIM = 200, 0.5, 70, 80, 0.6        # Treg-field sigma (um; chosen for readability), min kernel-weighted Treg count, percentiles, score limits
MIN_AREA = 2 * np.pi * NW_SIG ** 2 / BIN ** 2                     # CXCL12-High outline patches smaller than the Treg kernel's effective area (2 pi sigma^2) are dropped
APC_MIN_AREA = 2 * np.pi * 75 ** 2 / BIN ** 2                     # APC-Rich mask minimum area (0.035 mm2), fixed: does not follow the Treg kernel
MIN_TREG_SEC, NSHOW, PREF_CONTEXT = 400, 3, '54033_B4'
C_LOW, C_HIGH, C_MASK_A, C_OUT_B, C_TIS = '#2166ac', '#b2182b', '#f2f2f2', '#d9d9d9', '#bdbdbd'
C_CXFILL, C_TIS2, C_APCLINE, LW_APC = '#d6c4ec', '#f0f0f0', '#404040', 0.8   # example maps: CXCL12-High fill, tissue, APC-Rich line (pt)
# fill #d6c4ec: the hue of #e6ddf0 one step darker, so it separates from the #f0f0f0 tissue
DOT_S, DOT_CLIP = 34, 1.0                                         # example maps: Treg dot size (pt^2, uniform), score colour clip (z)
APC_GAMMA, APC_VMAX_PCT = 0.6, 98                                 # APC panel colour scale
INS_W, INS_H = 1.45, 1.35                                         # footprint (in) of the per-section inset incl. its labels
EXPECTED_PRODUCERS = ['OGN+RSPO3+ Fib', 'ADAMDEC1+ Fib', 'Activated Fib', 'Crypt Top Fib', 'PLA2G2A+ ECM Fib', 'Myofibroblast', 'Pericyte', 'Endothelial']   # as 11_neighbourhoods_local_tables.py
RAD, MIN_PROD = 50, 3                                             # as 11_neighbourhoods_local_tables.py
C_Q5HI, C_Q5LO, Q5_ALPHA = '#9e7cc1', '#5aae61', 0.55             # (not drawn here)
C_EDGE_HI, C_EDGE_LO, EDGE_LW, DOT_S_K = '#54278f', '#1b7837', 2.2, 80   # example dots: edge = own CXCL12 class, uniform size (pt^2)
C_CXFILL_67H, C_APCDOT, C_BGCELL = '#e6ddf0', '#08306b', '#e4e4e4'    # APC panel: smoothed CXCL12-High fill, APC dots, other cells
APC_S_67H, MAP_W_67H = 0.2, 2.11                                 # reference APC dot size (pt^2) on a 2.11 in wide map
FINE = 2.5   # um; grid of the test-defined regions in the shown sections: at a 10 um grid the test class can differ between a Treg and its bin
             # centre (<= 7 um away), because local CXCL12 averages few producer cells; a 2.5 um grid keeps every Treg within 1.8 um of its bin centre.
             # Other sections: 10 um.
LW_OUT, LW_TIS = 1.0, 0.6                                          # the one black line (pt); tissue edge
P2 = os.environ.get('PANEL2', 'ogn'); assert P2 == 'ogn', 'run with PANEL2=ogn'
sig3 = lambda v: float(f'{v:.3g}') if pd.notna(v) else v
def r3(df):                                                        # 3 significant figures, counts exact
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == 'f' and not (df[c].dropna() % 1 == 0).all(): df[c] = df[c].map(sig3)
    return df
def pfmt(p, tex=True):
    if p >= 0.001: return f'{p:.2g}'
    ex = int(np.floor(np.log10(p))); m = p / 10 ** ex
    if m >= 9.95: m, ex = 1.0, ex + 1
    return '{:.1f} × 10$^{{{}}}$'.format(m, ex) if tex else '{:.1f}×10{}'.format(m, str(ex).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')))
sgn = lambda v: f'{v:+.2f}'.replace('-', '−')
# ---------------- per-section top-APC-quintile numbers (recomputed from the Treg tables and asserted) ----------------
C = pd.read_csv(f'{LR}/c_apc_matched_cxcl12_high_vs_low_exp.csv').rename(columns={'suppression_cxcl12_high': 'act_cxcl12_high', 'suppression_cxcl12_low': 'act_cxcl12_low'})
PL = pd.read_csv(f'{LR}/c_pooled_apc_stratified_exp.csv')
IX = pd.read_csv(f'{LR}/c_interaction_apc_by_cxcl12_countadj_exp.csv').set_index('term')
PS = pd.read_csv(f'{LR}/c_top_apc_quintile_per_section_exp.csv')
SM = PS[PS.part.str.contains('per-section means')].set_index('section')[['n_low', 'n_high', 'mean_low', 'mean_high', 'difference_high_minus_low']]
GMr = PS[PS.part.str.contains('mean of section means')].set_index('group')
GM = {k: tuple(GMr.loc[f'CXCL12-{k.title()} Region', ['mean', 'ci_low', 'ci_high']]) for k in ['low', 'high']}
SGNr = PS[PS.group.eq('sections with lower CXCL12-High mean')].iloc[0]; WXr = PS[PS.group.eq('Wilcoxon signed-rank, section means')].iloc[0]
N_SEC, N_SEC_LOWER = len(SM), int(SGNr.n_sections_lower); assert N_SEC == 11 == int(SGNr.n_sections) == int(WXr.n_sections)
TA = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)[['section', 'act_z', 'apc_frac']]
TA = TA.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left')
cidx = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0, usecols=[0]).index.intersection(TA.index); assert len(cidx) == len(TA)
TA = TA.loc[cidx]                                                # Treg order of raw_counts_strict_tregs.csv (rank(method='first') ties depend on it)
TA['apc_q'] = TA.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
TA['t_exp'] = TA.groupby('section').lig_CXCL12_exp.transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))
SMr = TA[(TA.apc_q == 4) & TA.t_exp.isin(['low', 'high'])].groupby(['section', 't_exp'], observed=True).act_z.agg(['mean', 'size']).unstack()
SMr.columns = [f'{a}_{b}' for a, b in SMr.columns]; SMr = SMr.loc[SM.index]
close = lambda a, b: np.allclose(np.asarray(a, float), np.asarray(b, float), rtol=5e-3, atol=1e-4)   # the csv carries 3 significant figures
assert (SMr.size_low.values == SM.n_low.values).all() and (SMr.size_high.values == SM.n_high.values).all()
assert close(SMr.mean_low, SM.mean_low) and close(SMr.mean_high, SM.mean_high)
assert N_SEC_LOWER == int((SMr.mean_high < SMr.mean_low).sum()) and close(stats.binomtest(N_SEC_LOWER, N_SEC, 0.5).pvalue, SGNr.p)
assert close(stats.wilcoxon(SMr.mean_high, SMr.mean_low).pvalue, WXr.p)
tq = stats.t.ppf(0.975, N_SEC - 1)
for k in ['low', 'high']:
    m_, h_ = SMr[f'mean_{k}'].mean(), tq * SMr[f'mean_{k}'].std() / np.sqrt(N_SEC); assert close([m_, m_ - h_, m_ + h_], GM[k]), k
c5 = C.set_index('apc_quintile').loc[C.apc_quintile.max()]
assert int(c5.n_low) == int(SM.n_low.sum()) and int(c5.n_high) == int(SM.n_high.sum())
INT_P, TOP_Q = float(IX.p['hi:apcz']), float(c5.q); assert int(IX.n['hi:apcz']) == int(PL.n[0])
# ---------------- section selection (pre-stated rule; reads only effects and Treg counts) ----------------
SEL = pd.DataFrame({'n_strict_tregs': TA.section.value_counts().reindex(SMr.index), 'n_low': SMr.size_low, 'n_high': SMr.size_high,
                    'mean_low': SMr.mean_low, 'mean_high': SMr.mean_high, 'effect_high_minus_low': SMr.mean_high - SMr.mean_low})
assert close(SEL.effect_high_minus_low, SM.difference_high_minus_low)
MED = float(SEL.effect_high_minus_low.median())
SEL['abs_distance_to_median'] = (SEL.effect_high_minus_low - MED).abs(); SEL['candidate'] = SEL.n_strict_tregs >= MIN_TREG_SEC
cand = SEL[SEL.candidate].sort_values(['abs_distance_to_median', 'n_strict_tregs'], ascending=[True, False])
SEL['candidate_rank'] = pd.Series(np.arange(1, len(cand) + 1), index=cand.index).reindex(SEL.index)
CHOSEN = list(cand.index[:NSHOW]); CONTEXT = PREF_CONTEXT if PREF_CONTEXT in CHOSEN else CHOSEN[0]
SHOW = [CONTEXT] + [s for s in CHOSEN if s != CONTEXT]
SEL['chosen'] = SEL.index.isin(CHOSEN); SEL['panel_order'] = pd.Series(np.arange(1, NSHOW + 1), index=SHOW).reindex(SEL.index); SEL['context_row'] = SEL.index == CONTEXT
SEL = SEL.sort_values(['candidate', 'abs_distance_to_median'], ascending=[False, True])
print('median top-APC-quintile effect, 11 sections: {:+.4f}'.format(MED)); print(SEL.round(4).to_string(), flush=True)
print('chosen (closest to the median among >= {} strict Tregs): {} | context row: {} | row-2 order: {}'.format(MIN_TREG_SEC, CHOSEN, CONTEXT, SHOW), flush=True)
# ---------------- map settings and per-section preparation ----------------
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import distance_transform_edt
EXR = SEL[(SEL.n_strict_tregs >= 400) & (SEL.n_low >= 15) & (SEL.n_high >= 15)].sort_values('effect_high_minus_low')
assert list(EXR.index[:3]) == ['43563_A11', '43563_A4', '54033_B5'], list(EXR.index[:3])
SEC = os.environ.get('T67U_SEC', '43563_A4'); VER = os.environ.get('T67U_VER', 'v1'); assert VER in ('v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8'), VER
V8 = VER == 'v8'; VER = 'v7' if V8 else VER                                       # T67U_VER=v8 uses the code path of VER v7 with the V8 flags
V7_DERIVE = ['43563_A4', '11808_A4_Cos3', '54033_B4', '54033_B5', '43563_A11']   # the five sections the area cutoffs are derived from
V7_SECS = V7_DERIVE + ['43563_A8', '43563_A1']                                     # plus 43563_A8 and 43563_A1
V7_STACK = ['54033_B4', '43563_A4', '43563_A8', '43563_A1']                         # score stack, top to bottom
V7V = os.environ.get('T67U_V7_CN4', ''); assert V7V in ('', 'a', 'b') and not (V8 and V7V), V7V   # gate-2 variants: '' CN8 and CN4 >= 10%; 'a' CN4 >= 5%; 'b' CN8-rich (>= C8) exempt
V7TAG = 'v8' if V8 else 'v7' + V7V
V8FULL = V8 and os.environ.get('T67U_V8_FULL') == '1'                              # two-panel stack, rows A4, B4
if V8FULL: V7TAG = 'v8_full'
V8C = V8 and os.environ.get('T67U_V8C') == '1'                                   # 43563_A8 display crop + re-pack
V8B = V8 and (os.environ.get('T67U_V8B') == '1' or V8C); assert not V8B or (SEC == '43563_A8' and not V8FULL), 'T67U_V8B / T67U_V8C are for 43563_A8 only'
# T67U_V8C: areas, labels and cutoffs of the WHOLE section; display only: three tissue pieces are left out (the bottom-right two and the right-most middle one,
# no large APC aggregates) and the remaining pieces are re-packed with the same stitch() (pure per-piece translations, rounded to the 10 um grid);
# on-plot numbers use the Tregs of the displayed pieces; whole-section numbers are printed.
# T67U_V8B (43563_A8 only): minimum piece area 0.01 mm2 and no crypt / epithelial rule (the small round APC blobs are DC-rich pericryptal lamina
# propria, ~50% crypt share); otherwise the same density rule, CN8 + CN4 rule and merge
V8PAIR = os.environ.get('T67U_V8_PAIR', ''); assert V8PAIR in ('', 'A8') and (V8FULL or not V8PAIR), V8PAIR   # two-panel stack A4 + 43563_A8
if V8: V7_STACK = ['54033_B4', '43563_A4'] if not V8FULL else (['43563_A4', '43563_A8'] if os.environ.get('T67U_V8_PAIR') == 'A8' else ['43563_A4', '54033_B4'])   # stacks: score-only B4 + A4; two-panel A4 + B4
V7_TMP = os.environ.get('T67U_V7_TMP'); assert VER != 'v7' or SEC == 'all' or (V7_TMP and os.path.isdir(V7_TMP)), 'set T67U_V7_TMP to a scratch folder'
assert SEC in ('43563_A4', '11808_A4_Cos3', '54033_B4') or (VER == 'v7' and SEC in V7_SECS + ['all', 'stack']), SEC
LEAD = SEC; SECS = [SEC] if SEC not in ('all', 'stack') else list(V7_DERIVE); SHOW = SECS   # one section per run (T67U_SEC=all: the five derivation sections)
if VER == 'v7' and SEC == 'stack':   # T67U_SEC=stack: the per-section Treg Suppression Score panels, top to bottom, left-aligned on white, 60 px apart
    IMS = [Image.open(f'{V7_TMP}/{s_}_{V7TAG}_score.png').convert('RGB') for s_ in V7_STACK]
    W_, G_ = max(im.width for im in IMS), 60; CAN = Image.new('RGB', (W_, sum(im.height for im in IMS) + G_ * (len(IMS) - 1)), 'white'); y_ = 0
    for im in IMS: CAN.paste(im, (0, y_)); y_ += im.height + G_
    PNGS = (f'{OUT}/fig4_cxcl12_apc_clusters_A4_A8_v8.png' if V8PAIR else f'{OUT}/fig4_cxcl12_apc_clusters_2samples_v8.png') if V8FULL else f'{OUT}/fig4_treg_score_apc_areas_{V7TAG}.png'
    with open(PNGS, 'wb') as fh: CAN.save(fh, format='PNG', dpi=(300, 300))
    print('wrote', PNGS, [im.size for im in IMS], CAN.size, flush=True); raise SystemExit(0)
if VER == 'v1': CL_SIG_UM, CL_PCT, CL_MIN_APC, CL_MIN_MM2, TREG_NEAR_UM, CL_DEN_X, CL_CLOSE_UM = 20, float(os.environ.get('T67U_PCT', 90)), 20, 0.01, 20, 2.0, 10
else:   # band-scale APC areas: sigma 40 um, 85th percentile (fallback 80 only), 40 um closing + hole filling, >= 0.05 mm2, >= 50 APCs, no density floor;
        # Tregs counted 'in area' when inside it; smoothed outlines; outlines on the score panel only; smaller Treg dots
    CL_SIG_UM, CL_PCT, CL_MIN_APC, CL_MIN_MM2, TREG_NEAR_UM, CL_DEN_X, CL_CLOSE_UM = 40, float(os.environ.get('T67U_PCT', 85)), 50, 0.05, 0, 0.0, 40
if VER == 'v3':   # intermediate setting: sigma 30 um, 90th pct, 30 um closing + holes filled, >= 0.03 mm2, >= 40 APCs
    CL_SIG_UM, CL_PCT, CL_MIN_APC, CL_MIN_MM2, TREG_NEAR_UM, CL_DEN_X, CL_CLOSE_UM = 30, float(os.environ.get('T67U_PCT', 90)), 40, 0.03, 0, 0.0, 30
if VER == 'v4':   # dense APC aggregates by DBSCAN on APC coordinates (tissue um): eps 15 um (pre-stated fallback 20), min_samples 15, >= 150 APCs
    CL_SIG_UM, CL_PCT, CL_MIN_APC, CL_MIN_MM2, TREG_NEAR_UM, CL_DEN_X, CL_CLOSE_UM = 20, 90.0, 300, 0.0, 0, 0.0, 10   # >= 300 APCs per blob
DB_TOP = 10                                                         # at most the top 10 blobs by APC count
DB_EPS, DB_MIN, BLOB_SIG_BINS, BLOB_FRAC = float(os.environ.get('T67U_EPS', 15)), 15, 1.5, 0.2
if VER in ('v5', 'v6'):   # APC areas from the cellular neighborhoods: CN8 (lymphoid aggregates) + CN4 (DC niche); optional + CN2 (T67U_CN2=1)
    CL_SIG_UM, CL_PCT, CL_MIN_APC, CL_MIN_MM2, TREG_NEAR_UM, CL_DEN_X, CL_CLOSE_UM = 20, 90.0, 0, 0.05, 0, 0.0, 10
CN_SET = [8, 4] + ([2] if os.environ.get('T67U_CN2') == '1' else []); CN_SIG_UM, CN_FRAC = 20, 0.5
if VER == 'v6': CN_SET = [8, 4, 2]                               # CN8 + CN4 + CN2; keep an area only if CN8 + CN4 >= 25% of its CN8 + CN4 + CN2 cells
CN8_SHARE_MIN = 0.25; DROPPED = []; NCAND = 0
V7_SIG_UM, V7_Q, V7_CLOSE_UM, V7_MIN_MM2, V7_MERGE_UM, V7_CRYPT_Q = 30, 0.25, 40, 0.05, 50, 0.90   # APC density sigma, density quantile, closing, min area, merge gap, crypt quantile
V7_CN8_SIG_UM, V7_CN8_LEVEL, V7_OVAL_SOLID = 20, 0.5, 0.7                                          # CN8 ovals = large compact CN8 aggregates
V7_T, V7_C8, V7_CR = 4.28, 0.345, 0.080                                    # area cutoffs derived with T67U_SEC=all (five sections pooled; 3 significant figures)
V8_T, V8_C8, V8_CR, V8_G = 3.13, 0.303, 0.108, 5.67                     # density-rule cutoffs derived with T67U_VER=v8 T67U_SEC=all
V7_APC_PCT, V7_MIN_CN, V7_MIN_SMALL = 90, 0.10, 0.02                                           # gates: APC fraction >= section 90th percentile; CN8 and CN4 each >= 10%
if VER == 'v7': CL_SIG_UM, CL_PCT, CL_MIN_APC, CL_MIN_MM2, TREG_NEAR_UM, CL_DEN_X, CL_CLOSE_UM = 20, 90.0, 0, V7_MIN_MM2, 0, 0.0, 10   # (same values as the CN-based settings)
V2P = VER in ('v2', 'v3', 'v4', 'v5', 'v6', 'v7')                     # band-scale display and per-area rules
CELLS_ = os.environ.get('T67U_CELLS', 'apc'); assert CELLS_ in ('apc', 'dc') and (CELLS_ == 'apc' or VER == 'v2'), CELLS_   # APC areas (the DC-only setting is not used)
FOCUS = APC if CELLS_ == 'apc' else ['DC']; FNAME = 'APC' if CELLS_ == 'apc' else 'DC'   # cell types of the areas
DOT_AREA_X = 1.0 if VER == 'v1' else 0.47                         # Treg dot area factor
C_CLH, C_CLL, CL_LW = '#c51b7d', '#1a9850', 1.8
CLEAN_R_UM, CLEAN_MIN_MM2, KNN, DOT_UM = 50, 0.05, 20, 250
C_TISSUE, C_APCDOT, C_GREEN, LW_REG = '#ebebeb', '#08306b', '#1b7837', 1.0
APC_S_REF, IN_PER_MM_REF = 0.48, 3.28 / 7.61                       # reference APC dot size and map scale (density-matched APC dots)
CXCMAP = LinearSegmentedColormap.from_list('cxcl12_purple', ['#faf6fd', '#8e5fc0'])   # near-white -> medium purple (capped so dark-blue dots stay visible)
VAR = 'v6'  # map colour style used ('' and 'v1'-'v5': other background styles)
assert VAR in ('', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6'), VAR
V5_CMAP = LinearSegmentedColormap.from_list('v5_pink', ['#ffffff', '#dd3497', '#ae017e'])   # white -> strong pink -> magenta top
V5_LO_PCT = float(os.environ.get('T67S_V5_LO', 50)); V5_GAMMA = 1.8                        # vmin = section median (raise to 60-70 if the band does not stand out)
import matplotlib.patheffects as pe
V4_CMAP = LinearSegmentedColormap.from_list('v4_pink', ['#fff0f6', '#e7298a'])   # one continuous pink ramp over all tissue, never darker than mid magenta
C_BELOW = '#ececec'                                                            # tissue below the negative-probe background
C_V3_HI, C_V3_TIS = '#f6b3d9', '#efefef'                          # two-tone style: CXCL12-High Region light magenta-pink, other tissue light gray
B_HI = os.environ.get('T67S_HI', '#c51b7d'); B_LO, B_LW = '#1a9850', 2.2   # panel-4 borders: APC-rich & CXCL12-high / -low
CXCR4_CMAP_V2 = LinearSegmentedColormap.from_list('blue_cyan_white', ['#08306b', '#2171b5', '#00e5ff', '#ffffff'])
SA = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
SA['section'] = SA.patient_sec.astype(str); SA.index = [p + '_' + str(i) for p, i in zip(SA.section, SA.index)]
SA = SA[SA.section.isin(SECS)]; SA['x'], SA['y'] = SA.CenterX_global_px * UM, SA.CenterY_global_px * UM; SA['ct'] = SA.cell_type_general.astype(str)
TL = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
GAP, SBIN = 60, 50
def squeeze(v):                                                    # gap squeeze, as in F5E_cxcl12_region_tregs and the CN map of 12_neighbourhood_treg_score.py
    b = np.floor((v - v.min()) / SBIN).astype(int); occ = np.zeros(b.max() + 2, bool); occ[np.unique(b)] = True
    shift = np.zeros(len(occ)); run = 0; cut = 0
    for i in range(len(occ)):
        if not occ[i]:
            run += 1
            if run * SBIN > GAP: cut += SBIN
        else: run = 0
        shift[i] = cut
    return lambda x: x - shift[np.clip(np.floor((np.asarray(x) - v.min()) / SBIN).astype(int), 0, len(shift) - 1)]
def stitch(xy, SEC):                                               # piece stitching, with the 43563_A4 narrow-piece step of the CN map of 12_neighbourhood_treg_score.py
    xy = pd.DataFrame(xy, columns=['x', 'y'])
    gx = np.floor((xy.x - xy.x.min()) / SBIN).astype(int).values; gy = np.floor((xy.y - xy.y.min()) / SBIN).astype(int).values
    occ = np.zeros((gy.max() + 1, gx.max() + 1), bool); occ[gy, gx] = True
    lab, npc = ndimage.label(ndimage.binary_dilation(occ, iterations=1), structure=np.ones((3, 3)))
    piece = lab[gy, gx]; order = list(pd.Series(piece).value_counts().index)
    SHIFT = {order[0]: np.zeros(2)}; placed = xy.values[piece == order[0]]
    for pc in order[1:]:
        P = xy.values[piece == pc]; tree = cKDTree(placed[::5]); c0 = placed.mean(0)
        d = c0 - P.mean(0); Lh = np.linalg.norm(d); u = d / Lh; best = np.zeros(2); t = 0.0
        if V8 and tree.query(P[::5], k=1)[0].min() < GAP:              # a piece that overlaps placed pieces at its start position
            cand_ = []                                                 # (43563_A8: piece 4 under piece 8) is moved to the nearest clear position
            for ang in np.arange(0, 2 * np.pi, np.pi / 8):
                v_ = np.array([np.cos(ang), np.sin(ang)])
                for r_ in np.arange(25, 20001, 25):
                    if tree.query(P[::5] + v_ * r_, k=1)[0].min() >= GAP: cand_.append((r_, v_ * r_)); break
            best = min(cand_, key=lambda z: z[0])[1]; Lh = 0.0
            print('stitch: piece {} overlapped placed pieces at its start position; moved by ({:.0f}, {:.0f}) um'.format(pc, *best), flush=True)
        while t < Lh:
            t2 = t + 25; sh = u * t2
            if tree.query(P[::5] + sh, k=1)[0].min() < GAP: break
            best, t = sh, t2
        SHIFT[pc] = best; placed = np.vstack([placed, P + best])
    if SEC == '43563_A4':
        narrow = min(order[1:], key=lambda p: np.ptp(xy.x.values[piece == p]))
        others = np.vstack([xy.values[piece == p] + SHIFT[p] for p in order if p != narrow])[::3]; tree = cKDTree(others)
        Pn = xy.values[piece == narrow][::3]; big = order[0]; Bp = xy.values[piece == big] + SHIFT[big]
        base = np.array([Bp[:, 0].min() - GAP - Pn[:, 0].max(), Bp[:, 1].max() - Pn[:, 1].max()])
        best = None
        for dyu in np.arange(0, 4001, 50):
            for dxr in np.arange(0, 4001, 50):
                sh_ = base + np.array([dxr, -dyu])
                if tree.query(Pn + sh_, k=1)[0].min() < GAP: continue
                sc = dyu + dxr
                if best is None or sc > best[0]: best = (sc, sh_)
        SHIFT[narrow] = best[1]
    return piece, SHIFT, npc
def prep(LEAD):
    # ---- cells, CXCL12, display transform (one per section) ----
    sa = SA[SA.section == LEAD]
    with h5py.File(f'{R}/{FOLDER[LEAD]}/Processed/norm_log1p.h5ad', 'r') as h:
        ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
    nm = pd.Index(LEAD + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    keep = nm.isin(sa.index)
    E = pd.DataFrame(X[np.where(keep)[0]][:, [v.get_loc('CXCL12')]].toarray(), index=nm[keep], columns=['CXCL12']); del X; gc.collect()
    cells = sa.loc[E.index].copy(); print(LEAD, 'cells', len(cells), flush=True)
    # negative-probe background: mean over the 20 standard negative probes of the % of cells with >= 1 count of that probe, in this section's map cells.
    # nFeature_negprobes (obs) = number of the 20 probes detected per cell (checked against the exprMat flat file: equal for all 357,188 cells), so
    # mean(nFeature_negprobes) / 20 x 100 = the mean per-probe detection rate (2.53% here; flat file per probe 1.4-4.8%).
    OBN = pd.Series(ob['nFeature_negprobes'].values, index=nm).reindex(cells.index); assert OBN.notna().all()
    NEG_BG = float(OBN.mean() / 20 * 100); print('negative-probe background: {:.2f}% of cells per probe'.format(NEG_BG), flush=True)
    ox, oy = cells.x.values.copy(), cells.y.values.copy()
    fx, fy = squeeze(ox), squeeze(oy); sx, sy = fx(ox), fy(oy)
    piece, SHIFT, npc = stitch(np.c_[sx, sy], LEAD); cells['piece'] = piece
    shv = np.vstack([SHIFT[p] for p in piece]); cells['x'], cells['y'] = sx + shv[:, 0], sy + shv[:, 1]
    DISP = pd.DataFrame({'dx': cells.x.values - ox, 'dy': cells.y.values - oy, 'piece': piece}, index=cells.index)
    # the gap squeeze can move parts of one stitched piece by different amounts (a cell-free band crossing the whole section inside a piece, e.g.
    # 43563_A11); every (piece, translation) pair is therefore a display unit, and all per-piece steps below run per unit (pure translations, asserted)
    unit = pd.factorize(pd.Series(list(zip(piece, np.round(DISP.dx.values, 3), np.round(DISP.dy.values, 3)))))[0] + 1
    NUNIT = int(unit.max()); print('display units (piece x translation):', NUNIT, 'from', npc, 'pieces', flush=True)
    piece = unit; cells['piece'] = unit; DISP['piece'] = unit
    PIECE_T = DISP.groupby('piece')[['dx', 'dy']].agg(['mean', 'std', 'min', 'max'])
    for p_, r in PIECE_T.iterrows():                                   # each piece is moved by ONE translation (so per-piece smoothing in display coordinates = in tissue coordinates)
        assert r[('dx', 'max')] - r[('dx', 'min')] < 1e-6 and r[('dy', 'max')] - r[('dy', 'min')] < 1e-6, f'piece {p_} is not a pure translation'
    TR = {int(p_): (float(r[('dx', 'mean')]), float(r[('dy', 'mean')])) for p_, r in PIECE_T.iterrows()}
    print('pieces', npc, '| translation per piece (um):', {k: (round(a), round(b)) for k, (a, b) in TR.items()}, flush=True)
    cn = pq.read_table(f'{P4}/B_resection_analyses/cellular_neighborhoods_v2/cn_per_cell.parquet', columns=['section', 'x_um', 'y_um', '__index_level_0__']).to_pandas()   # guard: same layout as the CN map of 12_neighbourhood_treg_score.py
    cn = cn[cn.section == LEAD]; f5x, f5y = squeeze(cn.x_um.values), squeeze(cn.y_um.values); s5x, s5y = f5x(cn.x_um.values), f5y(cn.y_um.values)
    pc5, SH5, _ = stitch(np.c_[s5x, s5y], LEAD); sh5 = np.vstack([SH5[p] for p in pc5])
    D45 = pd.DataFrame({'dx': s5x + sh5[:, 0] - cn.x_um.values, 'dy': s5y + sh5[:, 1] - cn.y_um.values}, index=cn.index)
    com = DISP.index.intersection(D45.index); dd = np.hypot(*(DISP.loc[com, ['dx', 'dy']].values - D45.loc[com].values).T)
    print('CN-map layout comparison: {:,} of {:,} cells in cn_per_cell; 99th pct offset difference {:.1f} um'.format(len(com), len(cells), np.percentile(dd, 99) if len(com) else np.nan), flush=True)
    if LEAD == '43563_A4': assert len(com) > 0.95 * len(cells) and np.percentile(dd, 99) <= 25   # that layout is built for 43563_A4; asserted there only
    # Tregs: exact positions; display = tissue position + their piece's translation; 20-Treg means in TISSUE coordinates (true tissue neighbors)
    T = TL[TL.section == LEAD].copy(); T = T.join(sa[['x', 'y']], how='inner', rsuffix='_s'); T['ox'], T['oy'] = T.x_s.values, T.y_s.values
    T['piece'] = cells.piece.reindex(T.index).values; assert T.piece.notna().all()
    T['x'] = T.ox + T.piece.map(lambda p: TR[int(p)][0]); T['y'] = T.oy + T.piece.map(lambda p: TR[int(p)][1])
    _, nb = cKDTree(T[['ox', 'oy']].values).query(T[['ox', 'oy']].values, k=min(KNN, len(T)))
    T['score_knn'] = np.nanmean(T.act_z.values[nb], axis=1); T['cxcr4_knn'] = np.nanmean(T.rec_CXCR4.values[nb], axis=1)
    # ---- per-piece fields on the display grid (translations are exact, so each piece is smoothed alone: no spill across pieces or gaps) ----
    x0, y0 = cells.x.min() - 3 * SIG * BIN, cells.y.min() - 3 * SIG * BIN
    nx, ny = int((cells.x.max() - x0) / BIN) + 4 * int(SIG) + 4, int((cells.y.max() - y0) / BIN) + 4 * int(SIG) + 4
    ix, iy = ((cells.x - x0) / BIN).astype(int).values, ((cells.y - y0) / BIN).astype(int).values
    XC, YC = x0 + (np.arange(nx) + 0.5) * BIN, y0 + (np.arange(ny) + 0.5) * BIN; EXT = (x0, x0 + nx * BIN, y0 + ny * BIN, y0)
    PIDS = sorted(TR); ASM, OCC = {}, {}
    for p_ in PIDS:
        a = np.zeros((ny, nx)); sel = piece == p_; np.add.at(a, (iy[sel], ix[sel]), 1); OCC[p_] = a > 0; ASM[p_] = gaussian_filter(a, SIG, mode='constant')
    AMAX = max(a.max() for a in ASM.values())
    MASK0 = {p_: ndimage.binary_fill_holes(ASM[p_] > 0.05 * AMAX) for p_ in PIDS}     # mask rule of F5E_cxcl12_region_tregs, per piece (holes filled within the piece only)
    cover = sum(m.astype(int) for m in MASK0.values()); N_CONTESTED = int((cover > 1).sum())
    DIST = np.stack([distance_transform_edt(~OCC[p_]) for p_ in PIDS]); DIST[np.stack([~MASK0[p_] for p_ in PIDS])] = np.inf
    OWN = np.argmin(DIST, axis=0); ANY = cover > 0
    MASK = {p_: ANY & (OWN == k) for k, p_ in enumerate(PIDS)}         # contested pixels go to the piece with the nearest own cell
    TISSUE = np.zeros((ny, nx), bool)
    for m in MASK.values(): TISSUE |= m
    assert sum(m.astype(int) for m in MASK.values()).max() <= 1, 'piece masks overlap'
    LOST = {p_: (int((OCC[p_] & MASK0[p_] & ~MASK[p_]).sum()), int((OCC[p_] & ~MASK0[p_]).sum()), int(OCC[p_].sum())) for p_ in PIDS}
    print('occupied bins per piece: (lost to another piece, below the density threshold, total):', LOST, flush=True)
    OCCN = sum(o.astype(int) for o in OCC.values()); SHARED = {p_: int((OCC[p_] & MASK0[p_] & ~MASK[p_] & (OCCN > 1)).sum()) for p_ in PIDS}   # a 10 um bin holding cells of two units
    print('occupied bins shared with another display unit (tie, given to one of them):', {k: v for k, v in SHARED.items() if v}, flush=True)
    for p_ in PIDS: assert LOST[p_][0] == SHARED[p_], f'piece {p_}: own occupied bins given to another piece'
    print('per-piece masks: {} pixels claimed by two pieces before assignment; after: disjoint'.format(N_CONTESTED), flush=True)
    def field(vals):                                                   # % per piece, masked to that piece's own tissue; NaN (white) elsewhere
        out = np.full((ny, nx), np.nan)
        for p_ in PIDS:
            num = np.zeros((ny, nx)); sel = piece == p_; np.add.at(num, (iy[sel], ix[sel]), vals[sel].astype(float))
            f_ = gaussian_filter(num, SIG, mode='constant') / np.where(ASM[p_] > 0, ASM[p_], 1) * 100; out[MASK[p_]] = f_[MASK[p_]]
        return out
    FCX = field(E.CXCL12.values > 0); FAPC = field(cells.ct.isin(APC).values); FOGN = field((cells.ct == 'OGN+RSPO3+ Fib').values)
    for F_ in (FCX, FAPC, FOGN): assert (np.isnan(F_) == ~TISSUE).all(), 'a field has values outside tissue or holes inside'
    def per_piece(m, fn):
        out = np.zeros((ny, nx), bool)
        for p_ in PIDS: out |= fn(m & MASK[p_]) & MASK[p_]
        return out
    def region(F_, pct, min_area):                                    # region(), morphology per piece
        thr = np.nanpercentile(F_[TISSUE], pct)
        def one(m):
            m = ndimage.binary_fill_holes(ndimage.binary_closing(m, np.ones((3, 3)), iterations=2)); lab_, n_ = ndimage.label(m, structure=np.ones((3, 3)))
            return np.r_[False, ndimage.sum(m, lab_, np.arange(1, n_ + 1)) >= min_area][lab_]
        return thr, per_piece(np.nan_to_num(F_) >= thr, one)
    rad = int(round(CLEAN_R_UM / BIN)); yy, xx = np.mgrid[-rad:rad + 1, -rad:rad + 1]; DISK = (xx ** 2 + yy ** 2) <= rad ** 2; MINB = int(round(CLEAN_MIN_MM2 * 1e6 / BIN ** 2))
    def clean_one(m):                                                  # display cleaning
        m = ndimage.binary_opening(ndimage.binary_closing(m, structure=DISK), structure=DISK); lab_, n_ = ndimage.label(m, structure=np.ones((3, 3)))
        return np.r_[False, ndimage.sum(m, lab_, np.arange(1, n_ + 1)) >= MINB][lab_]
    clean = lambda m: per_piece(m, clean_one)
    thr_cx, CX = region(FCX, CX_PCT, MIN_AREA); thr_apc, AP = region(FAPC, APC_PCT, APC_MIN_AREA)
    CXd, APd = clean(CX), clean(AP); ZONE_HI = clean(CXd & APd); ZONE_LO = clean(APd & ~CXd)
    assert not (ZONE_HI & ZONE_LO).any()
    tix, tiy = ((T.x - x0) / BIN).astype(int).values, ((T.y - y0) / BIN).astype(int).values
    T['in_cxcl12_high'] = CXd[tiy, tix]; T['in_apc_rich_cxcl12_high'] = ZONE_HI[tiy, tix]; T['in_apc_rich_cxcl12_low'] = ZONE_LO[tiy, tix]
    t_ = T[['act_z']].join(TA[['apc_q', 't_exp']], how='left'); q_ = t_[t_.apc_q == 4]
    hi_, lo_ = q_[q_.t_exp == 'high'].act_z, q_[q_.t_exp == 'low'].act_z
    assert len(hi_) == SM.n_high[LEAD] and len(lo_) == SM.n_low[LEAD] and close(hi_.mean(), SM.mean_high[LEAD]) and close(lo_.mean(), SM.mean_low[LEAD])
    zh, zl = T[T.in_apc_rich_cxcl12_high], T[T.in_apc_rich_cxcl12_low]
    mm2 = BIN ** 2 / 1e6
    DESC = dict(n_tregs=len(T), n_apc_rich_cxcl12_high=len(zh), mean_score_apc_rich_cxcl12_high=zh.act_z.mean(), mean_knn_apc_rich_cxcl12_high=zh.score_knn.mean(),
                n_apc_rich_cxcl12_low=len(zl), mean_score_apc_rich_cxcl12_low=zl.act_z.mean(), mean_knn_apc_rich_cxcl12_low=zl.score_knn.mean(),
                mean_score_other=T[~T.in_apc_rich_cxcl12_high & ~T.in_apc_rich_cxcl12_low].act_z.mean(),
                cxcr4_knn_in_cxcl12_high=T[T.in_cxcl12_high].cxcr4_knn.mean(), cxcr4_knn_outside=T[~T.in_cxcl12_high].cxcr4_knn.mean(), n_tregs_in_cxcl12_high=int(T.in_cxcl12_high.sum()),
                tissue_mm2=TISSUE.sum() * mm2, cxcl12_high_mm2=CXd.sum() * mm2, apc_rich_mm2=APd.sum() * mm2, apc_rich_cxcl12_high_mm2=ZONE_HI.sum() * mm2,
                apc_rich_cxcl12_low_mm2=ZONE_LO.sum() * mm2, cxcl12_field_threshold_pct=thr_cx, apc_field_threshold_pct=thr_apc, contested_pixels_before_assignment=N_CONTESTED)
    print(pd.Series(DESC).round(3).to_string(), flush=True)
    return dict(cells=cells, T=T, EXT=EXT, XC=XC, YC=YC, TISSUE=TISSUE, FCX=FCX, FOGN=FOGN, CXd=CXd, APd=APd, ZONE_HI=ZONE_HI, ZONE_LO=ZONE_LO, DESC=DESC, TR=TR, NUNIT=NUNIT,
                MASK=MASK, PIDS=PIDS)
'''
    # 43563_A8 display crop: three pieces left out, the others re-packed (pure translations on the 10 um grid)
    SEG_A8_REPACK = r'''unit_ = cells.piece.values; OXY = SA.loc[cells.index, ['x', 'y']].values                 # tissue um
OC = {int(p_): OXY[unit_ == p_].mean(0) for p_ in PIDS}
REMOVE = []
for ref_ in [(10734, 4605), (13479, 10820), (12143, 8081)]:                              # right-most middle piece; the two bottom-right pieces (tissue-coordinate centroids)
    p_ = min(OC, key=lambda q: np.hypot(*(OC[q] - np.array(ref_)))); assert np.hypot(*(OC[p_] - np.array(ref_))) < 300, ref_; REMOVE.append(p_)
KEEP = [p_ for p_ in PIDS if p_ not in REMOVE]; kc_ = np.isin(unit_, KEEP)
print('A8 crop: display units removed {} (cells {}), kept {}'.format(REMOVE, [int((unit_ == p_).sum()) for p_ in REMOVE], KEEP), flush=True)
fx_, fy_ = squeeze(OXY[kc_, 0]), squeeze(OXY[kc_, 1]); sx_, sy_ = fx_(OXY[kc_, 0]), fy_(OXY[kc_, 1])
pc_, SH_, _ = stitch(np.c_[sx_, sy_], SEC); shv_ = np.vstack([SH_[q] for q in pc_])
NX_, NY_ = sx_ + shv_[:, 0], sy_ + shv_[:, 1]
DLT = {}
for p_ in KEEP:                                                                      # per kept unit: one translation (asserted), rounded to the 10 um grid
    m_ = unit_[kc_] == p_; ddx, ddy = NX_[m_] - cells.x.values[kc_][m_], NY_[m_] - cells.y.values[kc_][m_]
    assert np.ptp(ddx) < 1e-6 and np.ptp(ddy) < 1e-6, p_
    DLT[p_] = (int(round(ddx[0] / BIN)), int(round(ddy[0] / BIN)))
print('A8 re-pack shifts (bins of 10 um):', DLT, flush=True)
rr_all, cc_all = [], []
for p_ in KEEP:
    r_, c_ = np.nonzero(MASK[p_]); rr_all.append(r_ + DLT[p_][1]); cc_all.append(c_ + DLT[p_][0])
PAD = 3 * int(SIG); rmin_, cmin_ = min(a.min() for a in rr_all) - PAD, min(a.min() for a in cc_all) - PAD
NYc, NXc = max(a.max() for a in rr_all) - rmin_ + PAD + 1, max(a.max() for a in cc_all) - cmin_ + PAD + 1
x0c, y0c = x0_ + cmin_ * BIN, y0_ + rmin_ * BIN; EXTc = (x0c, x0c + NXc * BIN, y0c + NYc * BIN, y0c)
XCc, YCc = x0c + (np.arange(NXc) + 0.5) * BIN, y0c + (np.arange(NYc) + 0.5) * BIN
# tissue masks reach ~100 um beyond the cells; after re-packing, pixels claimed by two pieces go to the piece with the nearest own cell (as prep())
NM, OCn = {}, {}
for p_ in KEEP:
    r_, c_ = np.nonzero(MASK[p_]); a_ = np.zeros((NYc, NXc), bool); a_[r_ + DLT[p_][1] - rmin_, c_ + DLT[p_][0] - cmin_] = True; NM[p_] = a_
    m_ = unit_ == p_; o_ = np.zeros((NYc, NXc), bool)
    o_[((cells.y.values[m_] - y0_) / BIN).astype(int) + DLT[p_][1] - rmin_, ((cells.x.values[m_] - x0_) / BIN).astype(int) + DLT[p_][0] - cmin_] = True; OCn[p_] = o_
DSTn = np.stack([np.where(NM[p_], distance_transform_edt(~OCn[p_]), np.inf) for p_ in KEEP]); OWNn = np.array(KEEP)[np.argmin(DSTn, axis=0)]
ANYn = np.isfinite(DSTn).any(axis=0); NCONT = int((sum(NM[p_].astype(int) for p_ in KEEP) > 1).sum())
for p_ in KEEP: assert not (OCn[p_] & ANYn & (OWNn != p_) & ~np.any([OCn[q] for q in KEEP if q != p_], axis=0)).any(), f'piece {p_}: own occupied bins lost'
print('A8 re-pack: {} pixels claimed by two pieces, given to the piece with the nearest own cell'.format(NCONT), flush=True)
def remap(A, fill):
    out = np.full((NYc, NXc), fill, dtype=A.dtype)
    for p_ in KEEP:
        r_, c_ = np.nonzero(MASK[p_]); rn, cn_ = r_ + DLT[p_][1] - rmin_, c_ + DLT[p_][0] - cmin_; ok_ = OWNn[rn, cn_] == p_
        out[rn[ok_], cn_[ok_]] = A[r_[ok_], c_[ok_]]
    return out
'''
    # 43563_A8 display crop: shifted cell and Treg positions
    SEG_A8_SHIFT = r'''def shift_um(u): return np.array([DLT[int(q)][0] * BIN if int(q) in DLT else np.nan for q in u]), np.array([DLT[int(q)][1] * BIN if int(q) in DLT else np.nan for q in u])
cdx, cdy = shift_um(unit_); CXn, CYn = cells.x.values + cdx, cells.y.values + cdy
from scipy.spatial import cKDTree as _KD
for p_ in KEEP:                                                                      # pieces stay >= 50 um apart after rounding
    m_ = kc_ & (unit_ == p_); o_ = kc_ & (unit_ != p_)
    assert _KD(np.c_[CXn[o_], CYn[o_]][::3]).query(np.c_[CXn[m_], CYn[m_]][::3], k=1)[0].min() >= 50, p_
Tk = T.piece.astype(int).isin(KEEP).values; tdx, tdy = shift_um(T.piece.values); TXn, TYn = T.x.values + tdx, T.y.values + tdy
'''
    # figure helpers (frame, outline, flat, dot_s, scalebar, vbar, vkey, guards)
    SEG_HELPERS = r'''# ---------------- figure helpers ----------------
plt.rcParams.update({'font.size': 14})
asp = lambda e: abs(e[1] - e[0]) / abs(e[3] - e[2])
def ax_in(l, b, w, h): return fig.add_axes([l / FW, b / FH, w / FW, h / FH])
def frame(ax, d): E_ = d['EXT']; ax.set_xlim(E_[0], E_[1]); ax.set_ylim(E_[2], E_[3]); ax.set_aspect('equal'); ax.axis('off'); ax.set_facecolor('white')
def outline(ax, d, m, col, lw=LW_REG): ax.contour(d['XC'], d['YC'], m.astype(float), levels=[0.5], colors=col, linewidths=lw, zorder=6)
def flat(ax, d, m, col, z=0): ax.imshow(np.where(m, 1.0, np.nan), extent=d['EXT'], cmap=ListedColormap([col]), interpolation='nearest', zorder=z)
def dot_s(d, w_in): return (DOT_UM / 1000 * (w_in / abs(d['EXT'][1] - d['EXT'][0]) * 1000) * 72) ** 2
def scalebar(ax, d):
    m = ndimage.binary_dilation(d['TISSUE'], iterations=3); ny_, nx_ = m.shape; kw, kh = int(np.ceil(1300 / BIN)), int(np.ceil(450 / BIN))
    C_ = [n_ for n_, (r0, c0) in [('upper right', (0, nx_ - kw)), ('upper left', (0, 0)), ('lower right', (ny_ - kh, nx_ - kw)), ('lower left', (ny_ - kh, 0))]
          if not m[r0:r0 + kh, c0:c0 + kw].any()]
    c_ = C_[0] if C_ else 'lower left'; anc = {'upper left': (0, 1), 'upper right': (1, 1), 'lower left': (0, 0), 'lower right': (1, 0)}[c_]
    sb = AnchoredSizeBar(ax.transData, 1000, '1 mm', loc=c_, bbox_to_anchor=anc, bbox_transform=ax.transAxes, pad=0, borderpad=0.1, sep=2, frameon=False,
                         size_vertical=40, color='k', fontproperties=fm.FontProperties(size=11.5)); ax.add_artist(sb); return sb, c_
def score_panel(ax, d, w_in):
    T_ = d['T']; flat(ax, d, d['TISSUE'], C_TISSUE, z=0)
    im = ax.scatter(T_.x, T_.y, c=T_.score_knn, cmap='RdBu_r', vmin=VLO, vmax=VHI, s=dot_s(d, w_in), lw=0.2, edgecolor='0.35', alpha=0.85, rasterized=True, zorder=2)
    outline(ax, d, d['ZONE_HI'], B_HI, B_LW); outline(ax, d, d['ZONE_LO'], B_LO, B_LW); frame(ax, d); return im
def vbar(x, ytop, h, im, label, ext_=None, ticks=None):
    cax = ax_in(x + 0.05, ytop - h, 0.16, h); cb = fig.colorbar(im, cax=cax, orientation='vertical', **({'extend': ext_} if ext_ else {}))
    cb.ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    if ticks is not None: cb.set_ticks(ticks)
    else: cb.ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '\u2212'))); cb.ax.tick_params(labelsize=11.5, length=2.5, pad=1.5)
    cb.outline.set_linewidth(0.6); cax.set_title(label, fontsize=12, pad=5, loc='left', linespacing=1.15); return cax
def vkey(x, ytop, handles):
    return fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(x / FW, ytop / FH), frameon=False, fontsize=11.5, handlelength=1.1, handletextpad=0.35,
                      labelspacing=0.5, borderaxespad=0, borderpad=0)
SCORE_LABEL = 'Treg\nSuppression\nScore (z,\n20-Treg mean)'
ZONE_KEY = lambda: [Line2D([], [], color=B_HI, lw=B_LW, label='APC-Rich,\nCXCL12-High'), Line2D([], [], color=B_LO, lw=B_LW, label='APC-Rich,\nCXCL12-Low')]
def minitable(x_in, ytop_in, s_, fs=11):
    """2 x 2 mean Treg Suppression Score (n) + the per-Treg test line; text only, columns placed at fixed offsets (in)."""
    q = QT[s_]; f = lambda k: '{} ({})'.format(sgn(q[k][0]) if np.isfinite(q[k][0]) else 'n/a', q[k][1]); lh = fs * 1.35 / 72; TX = []
    rows = [('Mean score (n)', 'APC-Rich', 'APC-Poor'), ('CXCL12-High', f(('high', 'rich')), f(('high', 'poor'))), ('CXCL12-Low', f(('low', 'rich')), f(('low', 'poor')))]
    for i, r in enumerate(rows):
        for j, (t_, cx_) in enumerate(zip(r, [0.0, 1.30, 2.40])):
            TX.append(fig.text((x_in + cx_) / FW, (ytop_in - i * lh) / FH, t_, ha='left', va='top', fontsize=fs, color='0.15' if i == 0 or j == 0 else 'k'))
    TX.append(fig.text(x_in / FW, (ytop_in - 3.25 * lh) / FH, 'Test, most-APC Tregs, local CXCL12:\nHigh {} (n = {}) vs Low {} (n = {})'.format(
        sgn(SM.mean_high[s_]), int(SM.n_high[s_]), sgn(SM.mean_low[s_]), int(SM.n_low[s_])), ha='left', va='top', fontsize=fs, linespacing=1.3))
    return TX
def table_ok(TX):                                                  # cells of one mini-table must not touch each other
    rd = fig.canvas.get_renderer(); bb = [t.get_window_extent(rd) for t in TX]
    for i in range(len(bb)):
        for j in range(i + 1, len(bb)): assert not bb[i].overlaps(bb[j]), ('mini-table cells overlap', TX[i].get_text(), TX[j].get_text())
def guards(ITEMS, KEYS, OWNER, extra_pairs=()):
    fig.canvas.draw(); rd = fig.canvas.get_renderer()
    for i, k in enumerate(KEYS):
        ITEMS.append((f'k{i}', k.get_window_extent(rd)))
        if not isinstance(k, matplotlib.legend.Legend): ITEMS += [(f'k{i}', k.title.get_window_extent(rd))] + [(f'k{i}', t.get_window_extent(rd)) for t in k.get_yticklabels() if t.get_text()]
    for a_, (o1, b1) in enumerate(ITEMS):
        for o2, b2 in ITEMS[a_ + 1:]:
            if o1 == o2 or OWNER.get(o1) == o2 or OWNER.get(o2) == o1 or {o1, o2} in extra_pairs: continue
            if b1.overlaps(b2): raise AssertionError(f'overlap {o1} {o2}')
    for k_ in KEYS:
        if isinstance(k_, matplotlib.legend.Legend): continue
        lo_, hi_ = sorted(k_.get_ylim()); tb_ = sorted([t.get_window_extent(rd) for t, v_ in zip(k_.get_yticklabels(), k_.get_yticks()) if lo_ - 1e-9 <= v_ <= hi_ + 1e-9 and t.get_text()], key=lambda b: b.y0)
        assert all(a.y1 < b.y0 for a, b in zip(tb_, tb_[1:])), 'colorbar tick labels collide'
    for o, b in ITEMS: assert b.x0 >= fig.bbox.x0 - 1 and b.x1 <= fig.bbox.x1 + 1 and b.y0 >= fig.bbox.y0 - 1, f'{o} runs off the figure'
'''
    # colour scales of the CXCR4 and CXCL12 maps
    SEG_SCALES = r'''C4LO, C4HI = np.percentile(T.cxcr4_knn, 2), np.percentile(T.cxcr4_knn, 98)
CXLO, CX75, CXHI = np.nanpercentile(FCX[TISSUE], 50), np.nanpercentile(FCX[TISSUE], 75), np.nanpercentile(FCX[TISSUE], 98)   # CXCL12 scale: section median / 75th / 98th percentile
V6_CMAP = LinearSegmentedColormap.from_list('v6_pink', [(0.0, '#ffffff'), (CXLO / CXHI, '#fde0ef'), (CX75 / CXHI, '#f768a1'), (1.0, '#ae017e')])
CXKW = dict(cmap=V6_CMAP, norm=matplotlib.colors.Normalize(vmin=0, vmax=CXHI))
'''
    # panels 1-3 of the four-panel row
    SEG_PANELS123 = r'''ax = AX[0]; imO = ax.imshow(FOGN, extent=A['EXT'], cmap='viridis', vmin=0, vmax=np.nanpercentile(FOGN[TISSUE], 99), interpolation='nearest', zorder=1)
outline(ax, A, CXd, 'white'); ax.set_title('OGN\u207aRSPO3\u207a Fibroblasts', fontsize=15, pad=6); frame(ax, A)
def cx_background(ax):
    im_ = ax.imshow(FCX, extent=A['EXT'], interpolation='nearest', zorder=1, **CXKW)
    ax.contour(A['XC'], A['YC'], TISSUE.astype(float), levels=[0.5], colors='0.72', linewidths=0.6, zorder=2); return im_
ipm = CW / abs(A['EXT'][1] - A['EXT'][0]) * 1000; APC_S = APC_S_REF * (ipm / IN_PER_MM_REF) ** 2
ax = AX[1]; imC = cx_background(ax); cc = A['cells']; ap_ = cc.ct.isin(APC).values
ax.scatter(cc.x.values[ap_], cc.y.values[ap_], s=APC_S, c=C_APCDOT, lw=0, rasterized=True, zorder=3)
outline(ax, A, CXd, '#333333'); ax.set_title('Antigen-Presenting Cells\n(DC, Macrophage, B Cell)', fontsize=15, pad=6); frame(ax, A)
ax = AX[2]; cx_background(ax); o4 = np.argsort(T.cxcr4_knn.values)
im4 = ax.scatter(T.x.values[o4], T.y.values[o4], c=T.cxcr4_knn.values[o4], cmap='YlOrBr', vmin=C4LO, vmax=C4HI, s=dot_s(A, CW), edgecolors='0.15', linewidths=0.4, rasterized=True, zorder=3)
outline(ax, A, CXd, '#333333'); ax.set_title('Treg CXCR4\n(The CXCL12 Receptor)', fontsize=15, pad=6); frame(ax, A)
'''
    # render check: outside tissue is pure white
    SEG_RENDER_CHECK = r'''for k, (F_, kw_) in enumerate([(FOGN, dict(cmap='viridis', vmin=0, vmax=np.nanpercentile(FOGN[TISSUE], 99))), (FCX, CXKW)]):   # render check: outside tissue white
    f2 = plt.figure(figsize=(CW, RH), dpi=100); a2 = f2.add_axes([0, 0, 1, 1]); a2.imshow(F_, extent=A['EXT'], interpolation='nearest', **kw_)
    a2.set_xlim(A['EXT'][0], A['EXT'][1]); a2.set_ylim(A['EXT'][2], A['EXT'][3]); a2.axis('off'); f2.canvas.draw(); img = np.asarray(f2.canvas.buffer_rgba())[..., :3]; plt.close(f2)
    tis = ndimage.binary_dilation(np.asarray(Image.fromarray(TISSUE.astype(np.uint8) * 255).resize((img.shape[1], img.shape[0]), Image.NEAREST)) > 0, iterations=2)
    assert (img[~tis] == 255).all(), f'non-white pixels outside tissue in field {k}'
print('render check ok: outside tissue is pure white', flush=True)
'''
    U = types.ModuleType('prep_ns'); U.__file__ = SRC_U
    exec(compile(SEG_PREP, SRC_U, 'exec'), U.__dict__)                       # constants, data, squeeze(), stitch(), prep()
    assert U.V8 and U.SEC == SEC and U.V8C == A8 and U.DOT_AREA_X == 0.47 and not U.V8FULL
    TAG = 'A8' if A8 else 'A4'
    OUTD = os.environ.get('T94_OUTDIR', '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches')
    OUTS = dict(cn=f'{U.OUT}/fig4_cn_map_A8_plain.png', ogn=f'{U.OUT}/fig4_ogn_density_A8.png', row=f'{U.OUT}/fig4_{TAG}_four_panel_plain.png',
                two_full=f'{U.OUT}/fig4_A8_cn_treg_two_panel_full.png', two_cropped=f'{U.OUT}/fig4_A8_cn_treg_two_panel_cropped.png',
                piece=f'{U.OUT}/fig4_A8_cn_treg_two_panel_piece.png', piececx=f'{U.OUT}/fig4_A8_piece_cxcl12_apc.png')
    if not A8: OUTS = dict(row=OUTS['row'], two_full=f'{U.OUT}/fig4_A4_cn_treg_two_panel.png')   # 43563_A4: full layout only
    if SEC in FOVSEL: OUTS = dict(fov=f'{U.OUT}/fig4_{FOVSEL[SEC][1]}_piece_cn_treg_two_panel.png', **({'fovcx': f'{U.OUT}/fig4_B4_piece_cxcl12_apc.png'} if SEC == '54033_B4' else {}))
    OUTS = {k: v.replace(U.OUT, OUTD) for k, v in OUTS.items() if k.split('_')[0] in MAKE or k in MAKE}
    for f_ in OUTS.values(): assert not os.path.exists(f_), f'{f_} exists: not overwritten'
    import numpy as np, pandas as pd, matplotlib, seaborn as sns, pyarrow.parquet as pq
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scipy import ndimage
    from scipy.spatial import cKDTree
    d = U.prep(SEC)
    E67 = pd.read_csv(f'{U.OUT}/fig4_cxcl12_region_tregs_exp.csv', index_col=0)        # score colour limits from the Treg table of F5E_cxcl12_region_tregs
    vb = E67['suppression score smoothed'].values; LIM = np.nanpercentile(np.abs(vb - np.nanmean(vb)), 97); VLO, VHI = np.nanmean(vb) - LIM, np.nanmean(vb) + LIM
    assert abs(VLO + 0.685) < 2e-3 and abs(VHI - 0.723) < 2e-3
    if A8:   # ---- 43563_A8 display: three pieces left out, the other five re-packed by stitch() (pure translations on the 10 um grid); fields moved, not re-smoothed
        U.__dict__.update(cells=d['cells'], T=d['T'], MASK=d['MASK'], PIDS=d['PIDS'], TISSUE=d['TISSUE'], x0_=d['EXT'][0], y0_=d['EXT'][3])
        exec(compile(SEG_A8_REPACK, SRC_U, 'exec'), U.__dict__)
        exec(compile(SEG_A8_SHIFT, SRC_U, 'exec'), U.__dict__)
        rm, kc, Tk = U.remap, U.kc_, U.Tk
        A = dict(EXT=U.EXTc, XC=U.XCc, YC=U.YCc, TISSUE=rm(d['TISSUE'], False), CXd=rm(d['CXd'], False), FOGN=rm(d['FOGN'], np.nan), FCX=rm(d['FCX'], np.nan),
                 cells=d['cells'][kc].assign(x=U.CXn[kc], y=U.CYn[kc]), T=d['T'][Tk].assign(x=U.TXn[Tk], y=U.TYn[Tk]))
        assert A['TISSUE'].sum() == int(U.ANYn.sum()) and all((np.isnan(A[k]) == ~A['TISSUE']).all() for k in ('FOGN', 'FCX'))
        print('A8 display: pieces kept {} | cells {} of {} | Tregs {} of {}'.format(U.KEEP, len(A['cells']), len(d['cells']), len(A['T']), len(d['T'])), flush=True)
    else: A = {k: d[k] for k in ('EXT', 'XC', 'YC', 'TISSUE', 'CXd', 'FOGN', 'FCX', 'cells', 'T')}
    # ---- figure helpers and this section's colour scales (on the map shown) ----
    H = types.ModuleType('helper_ns'); H.__file__ = SRC_T; H.__dict__.update({k: v for k, v in U.__dict__.items() if not k.startswith('__')})
    exec(compile(SEG_HELPERS, SRC_T, 'exec'), H.__dict__)   # frame, outline, flat, dot_s, scalebar, vbar, vkey, guards
    H.__dict__.update(A=A, LEAD=SEC, T=A['T'], FOGN=A['FOGN'], FCX=A['FCX'], TISSUE=A['TISSUE'], CXd=A['CXd'], VLO=VLO, VHI=VHI)
    exec(compile(SEG_SCALES, SRC_T, 'exec'), H.__dict__)
    OGN_VMAX = float(np.nanpercentile(A['FOGN'][A['TISSUE']], 99))                    # OGN vmax rule of panel 1
    print('colour scales (map shown): OGN+RSPO3+ 0-{:.2f}% | CXCL12 0-{:.2f}% (median {:.2f}, p75 {:.2f}) | Treg CXCR4 {:.3f} to {:.3f} | whole-section OGN p99 {:.2f}%, '
      'CXCL12 p98 {:.2f}%'.format(OGN_VMAX, H.CXHI, H.CXLO, H.CX75, H.C4LO, H.C4HI, np.nanpercentile(d['FOGN'][d['TISSUE']], 99), np.nanpercentile(d['FCX'][d['TISSUE']], 98)), flush=True)
    C_BG, C_BG_EPI, C_BG_OTH, TREG_UM2 = '0.82', '0.62', '0.86', 90.0                    # piece backgrounds + 90 um Treg dots
    FOVSTYLE = {'54033_B4': dict(s=0.375, two_tone=False, treg_um=75.0), '43563_A11': dict(s=0.25, two_tone=True, treg_um=90.0)}   # per-piece marker sizes
    EPI_DEF = os.environ.get('T94_EPI', 'ct'); assert EPI_DEF in ('ct', 'cn'), EPI_DEF
    EPI_CT = ['Enterocyte', 'Secretory Epi', 'Microfold', 'Absorptive Epi', 'Absoprtive Epi', 'Stem']   # epithelial types ('Absoprtive' = the label's spelling)
    def bg_single(ax, cc_, tidx, s_):
        """non-Treg cells as tiny light-grey points (no edges, rasterized, under the Treg dots)"""
        nt = ~cc_.index.isin(tidx) & (cc_.ct.values != 'Treg')
        ax.scatter(cc_.x.values[nt], cc_.y.values[nt], s=s_, c=C_BG, rasterized=True, lw=0, alpha=1.0, zorder=1); return int(nt.sum())
    def bg_two_tone(ax, cc_, tidx, cn_, s_):
        """non-Treg cells as tiny points: epithelium medium grey on top of all other cells in light grey (no edges, rasterized, under the Treg dots)"""
        nt = ~cc_.index.isin(tidx) & (cc_.ct.values != 'Treg'); ep = cc_.ct.isin(EPI_CT).values if EPI_DEF == 'ct' else np.isin(cn_, [1, 6])
        ax.scatter(cc_.x.values[nt & ~ep], cc_.y.values[nt & ~ep], s=s_, c=C_BG_OTH, rasterized=True, lw=0, alpha=1.0, zorder=1)
        ax.scatter(cc_.x.values[nt & ep], cc_.y.values[nt & ep], s=s_, c=C_BG_EPI, rasterized=True, lw=0, alpha=1.0, zorder=1.5)
        return int((nt & ep).sum()), int((nt & ~ep).sum())
    def trim_ticks(caxs):   # keep only the ticks inside the colour range (the ones drawn), so the guards do not test hidden off-range tick labels against the keys
        for c_ in caxs: lo_, hi_ = sorted(c_.get_ylim()); c_.set_yticks([t for t in c_.get_yticks() if lo_ - 1e-9 <= t <= hi_ + 1e-9])
    def save(fig, path, pad):
        assert not os.path.exists(path), f'{path} exists: not overwritten'
        with open(path, 'wb') as fh: fig.savefig(fh, format='png', dpi=300, bbox_inches='tight', pad_inches=pad)
        plt.close(fig); print('wrote', path, flush=True)
    if 'row' in MAKE:
        # ---------------- four-panel row (FW 22.5 in, four maps CW x RH, colourbar column RC on the right of each) ----------------
        plt.rcParams.update({'font.size': 14})
        FW = 22.5; LM, RM, G, RC = 0.55, 0.10, 0.15, 1.45; T1, BU = 0.80, 0.25
        CW = (FW - LM - RM - 4 * RC - 4 * G) / 4; RH = CW / H.asp(A['EXT']); FH = T1 + RH + BU; y1 = BU
        fig = plt.figure(figsize=(FW, FH)); H.__dict__.update(fig=fig, FW=FW, FH=FH, CW=CW, RH=RH)
        X1 = [LM + j * (CW + RC + G) for j in range(4)]; AX = [H.ax_in(X1[j], y1, CW, RH) for j in range(4)]; H.AX = AX
        exec(compile(SEG_PANELS123, SRC_T, 'exec'), H.__dict__)   # panels 1-3
        ax = AX[3]; H.flat(ax, A, A['TISSUE'], U.C_TISSUE, z=0); Tq = A['T']                 # panel 4: plain score panel
        ipm = CW / abs(A['EXT'][1] - A['EXT'][0]) * 1000; DOT_S = 0.55 * U.DOT_AREA_X * (U.DOT_UM / 1000 * ipm * 72) ** 2
        imS = ax.scatter(Tq.x, Tq.y, c=Tq.score_knn, cmap='RdBu_r', vmin=VLO, vmax=VHI, s=DOT_S, lw=0, alpha=0.85, rasterized=True, zorder=2); H.frame(ax, A)
        ax.set_title('Treg Suppression Score', fontsize=15, pad=6)
        SB1 = H.scalebar(AX[0], A)
        LBL = fig.text(0.22 / FW, (y1 + RH / 2) / FH, SEC, rotation=90, ha='center', va='center', fontsize=15)
        XR = [X1[j] + CW + 0.08 for j in range(4)]; CBH = 0.42 * RH; CBT = y1 + RH - 0.75; KY = CBT - CBH - 0.30
        KEYB = [H.vbar(XR[0], CBT, CBH, H.imO, 'OGN⁺RSPO3⁺\n(% of cells)'),
                H.vbar(XR[1], CBT, CBH, H.imC, 'CXCL12\n(% of cells\nexpressing)', 'max', [t for t in np.arange(0, H.CXHI + 1e-9, 20)]),
                H.vbar(XR[2], CBT, CBH, H.im4, 'Treg CXCR4\n(z, 20-Treg\nmean)', 'both'),
                H.vbar(XR[3], CBT, CBH, imS, H.SCORE_LABEL, 'both'),
                H.vkey(XR[0], KY, [Line2D([], [], color='0.6', lw=1.8, label='CXCL12-High\nRegion (white\noutline)')]),
                H.vkey(XR[1], KY, [Line2D([], [], ls='', marker='o', ms=4.5, mfc=U.C_APCDOT, mec='none', label='APC'), Line2D([], [], color='#333333', lw=U.LW_REG, label='CXCL12-High\nRegion')]),
                H.vkey(XR[2], KY, [Line2D([], [], color='#333333', lw=U.LW_REG, label='CXCL12-High\nRegion')])]
        trim_ticks(KEYB[:4]); fig.canvas.draw(); rd = fig.canvas.get_renderer()
        IT = [(f'm{i}', a.get_window_extent(rd)) for i, a in enumerate(AX)] + [(f't{i}', a.title.get_window_extent(rd)) for i, a in enumerate(AX)] + \
             [('sb1', SB1[0].get_window_extent(rd)), ('lbl', LBL.get_window_extent(rd))]
        H.guards(IT, KEYB, {'sb1': 'm0'})                                                   # no overlaps, no colourbar tick collisions, nothing off the figure
        mb = [a.get_window_extent(rd) for a in AX]; assert all(mb[i].x1 < mb[i + 1].x0 for i in range(3)), 'panels overlap'
        assert max(abs(b.width - mb[0].width) for b in mb) < 1e-6 and max(abs(b.height - mb[0].height) for b in mb) < 1e-6, 'panels differ in size'
        print('four-panel row {:.2f} x {:.2f} in | each map {:.3f} x {:.3f} in | scale bar {}'.format(FW, FH, CW, RH, SB1[1]), flush=True)
        save(fig, OUTS['row'], 0.06)
        exec(compile(SEG_RENDER_CHECK, SRC_T, 'exec'), H.__dict__)   # outside tissue pure white
    if 'ogn' in MAKE:
        # ---------------- OGN+RSPO3+ density alone (panel 1 without the outline) ----------------
        plt.rcParams.update({'font.size': 14})
        Hc = 6.6; MWc = H.asp(A['EXT']) * Hc; LMc, KC2, RMc = 0.55, 1.75, 0.10; FWp, FHp = LMc + MWc + KC2 + RMc, 0.85 + Hc + 0.15
        fig = plt.figure(figsize=(FWp, FHp)); H.__dict__.update(fig=fig, FW=FWp, FH=FHp)
        ax = H.ax_in(LMc, 0.15, MWc, Hc)
        imO = ax.imshow(A['FOGN'], extent=A['EXT'], cmap='viridis', vmin=0, vmax=OGN_VMAX, interpolation='nearest', zorder=1)
        ax.set_title('OGN⁺RSPO3⁺ Fibroblasts', fontsize=15, pad=6); H.frame(ax, A)
        SBo = H.scalebar(ax, A); LBo = fig.text(0.22 / FWp, (0.15 + Hc / 2) / FHp, SEC, rotation=90, ha='center', va='center', fontsize=15)
        KO = [H.vbar(LMc + MWc + 0.08, 0.15 + Hc - 0.75, 0.42 * Hc, imO, 'OGN⁺RSPO3⁺\n(% of cells)')]
        trim_ticks(KO); fig.canvas.draw(); rd = fig.canvas.get_renderer()
        H.guards([('m', ax.get_window_extent(rd)), ('t', ax.title.get_window_extent(rd)), ('sb', SBo[0].get_window_extent(rd)), ('lbl', LBo.get_window_extent(rd))], KO, {'sb': 'm'})
        print('OGN map {:.2f} x {:.2f} in | scale bar {}'.format(MWc, Hc, SBo[1]), flush=True)
        save(fig, OUTS['ogn'], 0.06)
    if MAKE & {'cn', 'two', 'two_full', 'two_cropped', 'piece', 'piececx', 'fov', 'fovcx'}:   # ---- CN names (short, as the Fig 5B CN map), tab10 palette; CN_k10 per cell ----
        CND = f'{U.P4}/B_resection_analyses/cellular_neighborhoods_v2'
        CNS = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas().CN_k10
        NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
        NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
        TRIM = [('OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
                ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
                ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
                ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
        def short_cn(x):
            for a, b in TRIM: x = x.replace(a, b)
            return x
        NMS = {k: short_cn(v) for k, v in NM.items()}; assert sorted(NMS) == list(range(1, 11)) and all(NMS[k].startswith(f'CN{k}: ') for k in NMS)
        PAL = dict(zip(range(1, 11), sns.color_palette('tab10', 10).as_hex()))
        ENR = pd.read_csv(f'{CND}/cn_k10_celltype_log2_enrichment.csv', index_col=0); NAMES = {i + 1: n for i, n in enumerate(ENR.index)}   # full CN names
        assert all(NAMES[c].startswith(f'CN{c}: ') for c in NAMES) and len(NAMES) == 10
    if 'cn' in MAKE:
        # ---------------- CN map (style of the Fig 5B CN map of 14_treg_score_maps.py, section F5B_cn_section_map, without boxes) ----------------
        CNv = CNS.reindex(A['cells'].index)
        print('cells shown without a CN (not drawn: only CN cells are drawn): {} of {}'.format(int(CNv.isna().sum()), len(CNv)), flush=True)
        plt.rcParams.update({'font.size': 17})
        fig, ax = plt.subplots(figsize=(11, 11)); cx_, cy_ = A['cells'].x.values, A['cells'].y.values; cn_ = CNv.values
        for k in range(1, 11): m_ = cn_ == k; ax.scatter(cx_[m_], cy_[m_], s=0.25, c=PAL[k], rasterized=True, lw=0, alpha=1.0)
        ax.set_aspect('equal'); ax.invert_yaxis(); ax.axis('off')
        ok_ = ~np.isnan(cn_); xm, ym = cx_[ok_].min(), cy_[ok_].min()
        ax.plot([xm, xm + 1000], [ym + 150] * 2, 'k-', lw=4, solid_capstyle='butt'); TB = ax.text(xm + 500, ym + 60, '1 mm', ha='center', va='bottom', fontsize=15)
        assert not ((cx_ <= xm + 1000) & (cy_ <= ym + 250)).any(), 'scale bar would lie on tissue'
        ax.set_title(f'CD Resection {SEC}', fontsize=17)
        fig.canvas.draw(); XRT = ax.transData.transform([[cx_.max(), cy_.min()]])[0][0] / fig.bbox.width   # right edge of the tissue (figure fraction)
        # the Fig 5B CN map anchors its key at 0.88 of the figure width (wide 43563_A4 map); this map is narrower, so the key starts 0.5 in right of the tissue
        KC = fig.legend(handles=[Line2D([0], [0], marker='o', ls='', mfc=PAL[k], mec=PAL[k], ms=13, label=NMS[k]) for k in range(1, 11)], loc='center left',
                        bbox_to_anchor=(XRT + 0.5 / 11, 0.5), ncol=1, frameon=False, fontsize=15, handletextpad=0.5, labelspacing=0.6)
        fig.canvas.draw(); rd = fig.canvas.get_renderer(); tx0, tx1 = ax.transData.transform([[cx_.min(), cy_.min()], [cx_.max(), cy_.max()]])
        TIS_B = matplotlib.transforms.Bbox.from_extents(min(tx0[0], tx1[0]), min(tx0[1], tx1[1]), max(tx0[0], tx1[0]), max(tx0[1], tx1[1]))
        IT = [('tissue', TIS_B), ('title', ax.title.get_window_extent(rd)), ('key', KC.get_window_extent(rd)), ('bar', TB.get_window_extent(rd))]
        for i_, (o1, b1) in enumerate(IT):
            for o2, b2 in IT[i_ + 1:]:
                if {o1, o2} != {'tissue', 'bar'} and b1.overlaps(b2): raise AssertionError(f'overlap {o1} {o2}')   # the bar sits in the empty top-left corner (asserted above)
        save(fig, OUTS['cn'], 0.03)

    if MAKE & {'two', 'two_full', 'two_cropped'}:
        # ---------------- two-panel figures: CN map | Treg Suppression Score (plain), the same map size; full section and the 43563_A8 crop ----------------
        FULL = dict(EXT=d['EXT'], TISSUE=d['TISSUE'], cells=d['cells'], T=d['T'])          # prep(): all display units
        for key_, L_, Hp in [t_ for t_ in [('two_full', FULL, 5.2)] + ([('two_cropped', A, 6.6)] if A8 else []) if t_[0] in OUTS]:              # map heights: full section 5.2 in, crop 6.6 in
            plt.rcParams.update({'font.size': 14})
            MWp = H.asp(L_['EXT']) * Hp; LMp, KCN, Gp, KCB, RMp, T1p, BUp = 0.55, 3.60, 0.30, 1.75, 0.10, 0.85, 0.15
            FWp = LMp + MWp + KCN + Gp + MWp + KCB + RMp; FHp = T1p + Hp + BUp; X2 = LMp + MWp + KCN + Gp
            fig = plt.figure(figsize=(FWp, FHp)); H.__dict__.update(fig=fig, FW=FWp, FH=FHp)
            A1, A2 = H.ax_in(LMp, BUp, MWp, Hp), H.ax_in(X2, BUp, MWp, Hp)
            cc_ = L_['cells']; cn_ = CNS.reindex(cc_.index).values                             # CN map dots (s 0.25, tab10, CN1 first); cells without a CN not drawn
            for k in range(1, 11): m_ = cn_ == k; A1.scatter(cc_.x.values[m_], cc_.y.values[m_], s=0.25, c=PAL[k], rasterized=True, lw=0, alpha=1.0)
            H.frame(A1, L_); A1.set_title('Cellular Neighborhoods\n(k = 10)', fontsize=15, pad=6)
            Tq = L_['T']; bg_ = ~cc_.index.isin(Tq.index) & (cc_.ct.values != 'Treg')         # score: plain panel
            H.flat(A2, L_, L_['TISSUE'], U.C_TISSUE, z=0)
            ipm = MWp / abs(L_['EXT'][1] - L_['EXT'][0]) * 1000
            DOT_S = 0.55 * U.DOT_AREA_X * (U.DOT_UM / 1000 * ipm * 72) ** 2                   # plain score dot
            imS = A2.scatter(Tq.x, Tq.y, c=Tq.score_knn, cmap='RdBu_r', vmin=VLO, vmax=VHI, s=DOT_S, lw=0, alpha=0.85, rasterized=True, zorder=2); H.frame(A2, L_)
            A2.set_title('Treg Suppression Score', fontsize=15, pad=6)
            SBS = [H.scalebar(a_, L_) for a_ in (A1, A2)]
            LBL = fig.text(0.22 / FWp, (BUp + Hp / 2) / FHp, SEC, rotation=90, ha='center', va='center', fontsize=15)
            KCNL = fig.legend(handles=[Line2D([0], [0], marker='o', ls='', mfc=PAL[k], mec=PAL[k], ms=9, label=NMS[k]) for k in range(1, 11)], loc='upper left',
                              bbox_to_anchor=((LMp + MWp + 0.30) / FWp, (BUp + Hp - 0.05) / FHp), frameon=False, fontsize=12, handlelength=1.0, handletextpad=0.4,
                              labelspacing=0.55, borderaxespad=0, borderpad=0)
            CB = [H.vbar(X2 + MWp + 0.08, BUp + Hp - 0.75, 0.40 * Hp, imS, H.SCORE_LABEL, 'both')]
            trim_ticks(CB); fig.canvas.draw(); rd = fig.canvas.get_renderer()
            IT = [('m1', A1.get_window_extent(rd)), ('m2', A2.get_window_extent(rd)), ('t1', A1.title.get_window_extent(rd)), ('t2', A2.title.get_window_extent(rd)),
                  ('sb1', SBS[0][0].get_window_extent(rd)), ('sb2', SBS[1][0].get_window_extent(rd)), ('lbl', LBL.get_window_extent(rd))]
            H.guards(IT, CB + [KCNL], {'sb1': 'm1', 'sb2': 'm2'})                               # no overlaps, no tick collisions, nothing off the figure
            b1_, b2_ = A1.get_window_extent(rd), A2.get_window_extent(rd); assert abs(b1_.width - b2_.width) < 1e-6 and abs(b1_.height - b2_.height) < 1e-6
            print('non-Treg cells drawn {} | cells typed Treg but not in the strict Treg set (not drawn) {}'.format(int(bg_.sum()), int(((cc_.ct.values == 'Treg') & ~cc_.index.isin(Tq.index)).sum())), flush=True)
            print('{}: maps {:.2f} x {:.2f} in | cells {} ({} without a CN) | Tregs {} | Treg dot {:.0f} um = {:.2f} pt across | scale bar {} | CN key {:.2f} in wide'.format(
                  key_, MWp, Hp, len(cc_), int(np.isnan(cn_).sum()), len(Tq), np.sqrt(DOT_S) / 72 / ipm * 1000, np.sqrt(DOT_S), SBS[0][1], KCNL.get_window_extent(rd).width / fig.dpi), flush=True)
            save(fig, OUTS[key_], 0.06)
    if MAKE & {'piece', 'piececx'}:
        # ---------------- one display unit alone: the large lower-left piece of the full layout (CN8 aggregate, dark-red Treg cluster) ----------------
        cf, Tf = d['cells'], d['T']; uu = cf.piece.values.astype(int)
        UT = pd.DataFrame({'unit': uu, 'x': cf.x.values, 'y': cf.y.values, 'cn8': CNS.reindex(cf.index).values == 8}).groupby('unit').agg(
            n_cells=('x', 'size'), x_um=('x', 'mean'), y_um=('y', 'mean'), n_cn8=('cn8', 'sum'))
        UT = UT.join(Tf.groupby(Tf.piece.astype(int)).agg(n_tregs=('act_z', 'size'), mean_raw_score=('act_z', 'mean'), max_score_20treg=('score_knn', 'max')))
        UT['in_v8c'] = UT.index.isin(U.KEEP)
        print(UT.round(3).to_string(), flush=True)
        LL = UT[(UT.x_um < np.median(cf.x.values)) & (UT.y_um > np.median(cf.y.values))]; PU = int(LL.n_cells.idxmax())   # y grows downward on the maps
        tp_ = Tf[Tf.piece.astype(int) == PU]
        print('piece: display unit {} (of {}; kept in the 43563_A8 crop: {}) | {} cells, {} CN8 cells | n = {} Tregs, mean raw Treg Suppression Score {:+.3f} (section: n = {}, mean {:+.3f})'.format(
              PU, len(UT), PU in U.KEEP, int(UT.n_cells[PU]), int(UT.n_cn8[PU]), len(tp_), tp_.act_z.mean(), len(Tf), Tf.act_z.mean()), flush=True)
        M_ = d['MASK'][PU]; r_, c_ = np.nonzero(M_); PADp = 3 * int(U.SIG); ny_, nx_ = M_.shape
        r0, r1, c0, c1 = max(r_.min() - PADp, 0), min(r_.max() + PADp + 1, ny_), max(c_.min() - PADp, 0), min(c_.max() + PADp + 1, nx_)
        x0f, y0f = d['EXT'][0], d['EXT'][3]
        LP = dict(EXT=(x0f + c0 * U.BIN, x0f + c1 * U.BIN, y0f + r1 * U.BIN, y0f + r0 * U.BIN), TISSUE=M_[r0:r1, c0:c1], cells=cf[uu == PU], T=tp_)
        assert (LP['cells'].x.between(LP['EXT'][0], LP['EXT'][1]) & LP['cells'].y.between(LP['EXT'][3], LP['EXT'][2])).all()
        plt.rcParams.update({'font.size': 14})
        Hp = 6.0; MWp = H.asp(LP['EXT']) * Hp; LMp, Gp, KCB, RMp, T1p, KEYH = 0.55, 0.40, 1.75, 0.10, 0.85, 1.45
        FWp = LMp + MWp + Gp + MWp + KCB + RMp; FHp = T1p + Hp + KEYH; BUp = KEYH; X2 = LMp + MWp + Gp
        fig = plt.figure(figsize=(FWp, FHp)); H.__dict__.update(fig=fig, FW=FWp, FH=FHp)
        A1, A2 = H.ax_in(LMp, BUp, MWp, Hp), H.ax_in(X2, BUp, MWp, Hp)
        ipm = MWp / abs(LP['EXT'][1] - LP['EXT'][0]) * 1000
        cc_ = LP['cells']; cn_ = CNS.reindex(cc_.index).values; s_cn = (9 / 1000 * ipm * 72) ** 2   # CN dot of the Fig 5B map (0.25 pt2 at 0.77 in/mm) at the same physical size, 9 um
        for k in range(1, 11): m_ = cn_ == k; A1.scatter(cc_.x.values[m_], cc_.y.values[m_], s=s_cn, c=PAL[k], rasterized=True, lw=0, alpha=1.0)
        H.frame(A1, LP); A1.set_title('Cellular Neighborhoods\n(k = 10)', fontsize=15, pad=6)
        bg_ = ~cc_.index.isin(tp_.index) & (cc_.ct.values != 'Treg')
        print('background cells (single grey):', bg_single(A2, cc_, tp_.index, s_cn), flush=True)   # non-Treg cells, CN-panel marker
        DOT_S = (TREG_UM2 / 1000 * ipm * 72) ** 2                                            # 90 um Treg dots
        imS = A2.scatter(tp_.x, tp_.y, c=tp_.score_knn, cmap='RdBu_r', vmin=VLO, vmax=VHI, s=DOT_S, lw=0, alpha=0.85, rasterized=True, zorder=2); H.frame(A2, LP)
        A2.set_title('Treg Suppression Score', fontsize=15, pad=6)
        SBS = [H.scalebar(a_, LP) for a_ in (A1, A2)]
        LBL = fig.text(0.22 / FWp, (BUp + Hp / 2) / FHp, SEC, rotation=90, ha='center', va='center', fontsize=15)
        dot = lambda c, lab_: Line2D([], [], ls='', marker='o', ms=8, mfc=c, mec='none', label=lab_)   # key style
        KCNL = fig.legend(handles=[dot(PAL[k], NAMES[k]) for k in range(1, 11)], ncol=2, loc='upper left', bbox_to_anchor=(LMp / FWp, (BUp - 0.20) / FHp), frameon=False,
                          fontsize=11.5, handletextpad=0.2, labelspacing=0.32, columnspacing=1.6, borderaxespad=0, borderpad=0)
        CB = [H.vbar(X2 + MWp + 0.08, BUp + Hp - 0.75, 0.40 * Hp, imS, H.SCORE_LABEL, 'both')]
        trim_ticks(CB); fig.canvas.draw(); rd = fig.canvas.get_renderer()
        IT = [('m1', A1.get_window_extent(rd)), ('m2', A2.get_window_extent(rd)), ('t1', A1.title.get_window_extent(rd)), ('t2', A2.title.get_window_extent(rd)),
              ('sb1', SBS[0][0].get_window_extent(rd)), ('sb2', SBS[1][0].get_window_extent(rd)), ('lbl', LBL.get_window_extent(rd))]
        H.guards(IT, CB + [KCNL], {'sb1': 'm1', 'sb2': 'm2'})
        b1_, b2_ = A1.get_window_extent(rd), A2.get_window_extent(rd); assert abs(b1_.width - b2_.width) < 1e-6 and abs(b1_.height - b2_.height) < 1e-6
        print('non-Treg cells drawn {} | cells typed Treg but not in the strict Treg set (not drawn) {}'.format(int(bg_.sum()), int(((cc_.ct.values == 'Treg') & ~cc_.index.isin(tp_.index)).sum())), flush=True)
        print('piece figure {:.2f} x {:.2f} in | maps {:.2f} x {:.2f} in | Treg dot {:.0f} um = {:.2f} pt across | scale bar {} | key {:.2f} x {:.2f} in'.format(FWp, FHp, MWp, Hp,
              np.sqrt(DOT_S) / 72 / ipm * 1000, np.sqrt(DOT_S), SBS[0][1],
              KCNL.get_window_extent(rd).width / fig.dpi, KCNL.get_window_extent(rd).height / fig.dpi), flush=True)
        if 'piece' in OUTS: save(fig, OUTS['piece'], 0.06)
        else: plt.close(fig)
    if 'piececx' in OUTS:
        # ---------------- A8 piece: APCs over the CXCL12 field (panel 2 style), as fig4_B4_piece_cxcl12_apc.png; fields of the whole section, unit PU's own tissue ----------------
        G7 = dict(H.__dict__); G7.update(T=d['T'], FCX=d['FCX'], TISSUE=d['TISSUE'])         # colour scales of the WHOLE section (not the crop)
        exec(compile(SEG_SCALES, SRC_T, 'exec'), G7)
        plt.rcParams.update({'font.size': 14}); BUc = 0.15
        FWc = LMp + MWp + KCB + RMp; FHc = T1p + Hp + BUc                                    # the piece figure's map size; colourbar and key on the right
        fig = plt.figure(figsize=(FWc, FHc)); H.__dict__.update(fig=fig, FW=FWc, FH=FHc); ax = H.ax_in(LMp, BUc, MWp, Hp)
        ok_ = d['MASK'][PU]                                                                  # display coordinates of the full layout, as the piece figure
        imC = ax.imshow(np.where(ok_, d['FCX'], np.nan), extent=d['EXT'], interpolation='nearest', zorder=1, **G7['CXKW'])   # as cx_background()
        ax.contour(d['XC'], d['YC'], ok_.astype(float), levels=[0.5], colors='0.72', linewidths=0.6, zorder=2)
        ax.contour(d['XC'], d['YC'], (d['CXd'] & ok_).astype(float), levels=[0.5], colors='#333333', linewidths=U.LW_REG, zorder=6)   # CXCL12-High Region
        ap_ = cc_.ct.isin(U.APC).values; APC_S = U.APC_S_REF * (ipm / U.IN_PER_MM_REF) ** 2   # density-matched APC dots
        ax.scatter(cc_.x.values[ap_], cc_.y.values[ap_], s=APC_S, c=U.C_APCDOT, lw=0, rasterized=True, zorder=3); H.frame(ax, LP)
        ax.set_title('Antigen-Presenting Cells\n(DC, Macrophage, B Cell)', fontsize=15, pad=6)
        SBc = H.scalebar(ax, LP); LBc = fig.text(0.22 / FWc, (BUc + Hp / 2) / FHc, SEC, rotation=90, ha='center', va='center', fontsize=15)
        CBHc, CBTc = 0.40 * Hp, BUc + Hp - 0.75
        KX = [H.vbar(LMp + MWp + 0.08, CBTc, CBHc, imC, 'CXCL12\n(% of cells\nexpressing)', 'max', [t for t in np.arange(0, G7['CXHI'] + 1e-9, 20)]),
              H.vkey(LMp + MWp + 0.08, CBTc - CBHc - 0.30, [Line2D([], [], ls='', marker='o', ms=4.5, mfc=U.C_APCDOT, mec='none', label='APC'),
                                                            Line2D([], [], color='#333333', lw=U.LW_REG, label='CXCL12-High\nRegion')])]
        trim_ticks(KX[:1]); fig.canvas.draw(); rd = fig.canvas.get_renderer()
        H.guards([('m', ax.get_window_extent(rd)), ('t', ax.title.get_window_extent(rd)), ('sb', SBc[0].get_window_extent(rd)), ('lbl', LBc.get_window_extent(rd))], KX, {'sb': 'm'})
        print('A8 piece CXCL12 / APC map {:.2f} x {:.2f} in | APCs {} | CXCL12 scale 0-{:.1f}% (whole section; median {:.1f}, p75 {:.1f}) | scale bar {}'.format(MWp, Hp,
              int(ap_.sum()), G7['CXHI'], G7['CXLO'], G7['CX75'], SBc[1]), flush=True)
        save(fig, OUTS['piececx'], 0.06)
    if MAKE & {'fov', 'fovcx'}:
        # ---------------- one tissue piece chosen by FOV, tissue coordinates, the cropped two-panel style ----------------
        BLK, TAGF = FOVSEL[SEC]
        sa_ = U.SA[U.SA.section == SEC].copy(); sa_['fov'] = sa_.index.str.rsplit('_', n=1).str[-1].astype(int)   # cell id = section_cell_fov
        GB = 50; gx = ((sa_.x.values - sa_.x.min()) // GB).astype(int); gy = ((sa_.y.values - sa_.y.min()) // GB).astype(int)   # piece rule of the CN map (50 um grid, 1-bin dilation, 8-neighbour components)
        occ = np.zeros((gy.max() + 1, gx.max() + 1), bool); occ[gy, gx] = True
        labc = ndimage.label(ndimage.binary_dilation(occ, iterations=1), structure=np.ones((3, 3)))[0][gy, gx]
        m0 = sa_.fov.between(*BLK[0]).values; keep_ = m0.copy(); comp0 = set(np.unique(labc[m0]))
        print('FOV {}-{}: {} cells, fragment(s) {} ({} cells in that fragment)'.format(*BLK[0], int(m0.sum()), sorted(comp0), int(np.isin(labc, list(comp0)).sum())), flush=True)
        for a_, b_ in BLK[1:]:
            mb = sa_.fov.between(a_, b_).values; same = set(np.unique(labc[mb])) <= comp0
            dmin = cKDTree(np.c_[sa_.x.values[m0], sa_.y.values[m0]]).query(np.c_[sa_.x.values[mb], sa_.y.values[mb]], k=1)[0].min()
            print('FOV {}-{}: {} cells, fragment(s) {} | nearest cell to FOVs {}-{}: {:.0f} um | same fragment: {} -> {}'.format(a_, b_, int(mb.sum()), sorted(set(np.unique(labc[mb]))),
                  *BLK[0], dmin, same, 'shown together' if same else 'not shown'), flush=True)
            if same: keep_ |= mb
        cc_ = sa_[keep_]; Tq = d['T'][d['T'].index.isin(cc_.index)].copy(); Tq['x'], Tq['y'] = Tq.ox.values, Tq.oy.values   # tissue coordinates
        assert np.allclose(Tq.x.values, cc_.loc[Tq.index, 'x'].values) and np.allclose(Tq.y.values, cc_.loc[Tq.index, 'y'].values)
        print('piece {}: {} cells, n = {} Tregs, mean raw Treg Suppression Score {:+.3f}'.format(TAGF, len(cc_), len(Tq), Tq.act_z.mean()), flush=True)
        PADf, BIN_ = 150, U.BIN; x0, y0 = cc_.x.min() - PADf, cc_.y.min() - PADf
        nxf, nyf = int(np.ceil((cc_.x.max() + PADf - x0) / BIN_)), int(np.ceil((cc_.y.max() + PADf - y0) / BIN_))
        oc = np.zeros((nyf, nxf), bool); oc[((cc_.y.values - y0) // BIN_).astype(int), ((cc_.x.values - x0) // BIN_).astype(int)] = True
        LF = dict(EXT=(x0, x0 + nxf * BIN_, y0 + nyf * BIN_, y0), TISSUE=ndimage.binary_dilation(oc, iterations=5))   # tissue mask: scale-bar placement only
        plt.rcParams.update({'font.size': 14})
        asp_ = H.asp(LF['EXT']); Hp = min(6.6, 7.5 / asp_); MWp = asp_ * Hp                    # the cropped figure's map height (6.6 in), maps at most 7.5 in wide
        LMp, KCN, Gp, KCB, RMp, T1p, BUp = 0.55, 4.90, 0.30, 1.75, 0.10, 0.85, 0.15          # KCN: room for the full CN names
        FWp = LMp + MWp + KCN + Gp + MWp + KCB + RMp; FHp = T1p + Hp + BUp; X2 = LMp + MWp + KCN + Gp
        fig = plt.figure(figsize=(FWp, FHp)); H.__dict__.update(fig=fig, FW=FWp, FH=FHp)
        A1, A2 = H.ax_in(LMp, BUp, MWp, Hp), H.ax_in(X2, BUp, MWp, Hp)
        cn_ = CNS.reindex(cc_.index).values
        FS_ = FOVSTYLE[SEC]
        for k in range(1, 11): m_ = cn_ == k; A1.scatter(cc_.x.values[m_], cc_.y.values[m_], s=FS_['s'], c=PAL[k], rasterized=True, lw=0, alpha=1.0)
        H.frame(A1, LF); A1.set_title('Cellular Neighborhoods\n(k = 10)', fontsize=15, pad=6)
        print('background cells:', bg_two_tone(A2, cc_, Tq.index, cn_, FS_['s']) if FS_['two_tone'] else bg_single(A2, cc_, Tq.index, FS_['s']), flush=True)   # CN-panel marker
        ipm = MWp / abs(LF['EXT'][1] - LF['EXT'][0]) * 1000; DOT_S = (FS_['treg_um'] / 1000 * ipm * 72) ** 2   # Treg dots
        imS = A2.scatter(Tq.x, Tq.y, c=Tq.score_knn, cmap='RdBu_r', vmin=VLO, vmax=VHI, s=DOT_S, lw=0, alpha=0.85, rasterized=True, zorder=2); H.frame(A2, LF)
        A2.set_title('Treg Suppression Score', fontsize=15, pad=6)
        SBS = [H.scalebar(a_, LF) for a_ in (A1, A2)]
        LBL = fig.text(0.22 / FWp, (BUp + Hp / 2) / FHp, SEC, rotation=90, ha='center', va='center', fontsize=15)
        KCNL = fig.legend(handles=[Line2D([0], [0], marker='o', ls='', mfc=PAL[k], mec=PAL[k], ms=9, label=NAMES[k]) for k in range(1, 11)], loc='upper left',
                          bbox_to_anchor=((LMp + MWp + 0.30) / FWp, (BUp + Hp - 0.05) / FHp), frameon=False, fontsize=12, handlelength=1.0, handletextpad=0.4,
                          labelspacing=0.55, borderaxespad=0, borderpad=0)
        CB = [H.vbar(X2 + MWp + 0.08, BUp + Hp - 0.75, 0.40 * Hp, imS, H.SCORE_LABEL, 'both')]
        trim_ticks(CB); fig.canvas.draw(); rd = fig.canvas.get_renderer()
        IT = [('m1', A1.get_window_extent(rd)), ('m2', A2.get_window_extent(rd)), ('t1', A1.title.get_window_extent(rd)), ('t2', A2.title.get_window_extent(rd)),
              ('sb1', SBS[0][0].get_window_extent(rd)), ('sb2', SBS[1][0].get_window_extent(rd)), ('lbl', LBL.get_window_extent(rd))]
        H.guards(IT, CB + [KCNL], {'sb1': 'm1', 'sb2': 'm2'})
        b1_, b2_ = A1.get_window_extent(rd), A2.get_window_extent(rd); assert abs(b1_.width - b2_.width) < 1e-6 and abs(b1_.height - b2_.height) < 1e-6
        print('fov figure {:.2f} x {:.2f} in | maps {:.2f} x {:.2f} in | Treg dot {:.0f} um | scale bar {} | key {:.2f} x {:.2f} in'.format(FWp, FHp, MWp, Hp,
              np.sqrt(DOT_S) / 72 / ipm * 1000, SBS[0][1], KCNL.get_window_extent(rd).width / fig.dpi, KCNL.get_window_extent(rd).height / fig.dpi), flush=True)
        if 'fov' in OUTS: save(fig, OUTS['fov'], 0.06)
        else: plt.close(fig)
    if 'fovcx' in OUTS:
        # ---------------- B4 piece: APCs over the CXCL12 field (panel 2 style); whole-section fields per display unit, shown in tissue coordinates ----------------
        cd_ = d['cells']; inpc = cd_.index.isin(cc_.index); units = sorted(set(cd_.piece.values[inpc].astype(int)))
        x0d, y0d = d['EXT'][0], d['EXT'][3]; nyd, nxd = d['TISSUE'].shape
        plt.rcParams.update({'font.size': 14})
        FWc = LMp + MWp + KCB + RMp; FHc = T1p + Hp + BUp                                    # the two-panel figure's map size; colourbar and key on the right
        fig = plt.figure(figsize=(FWc, FHc)); H.__dict__.update(fig=fig, FW=FWc, FH=FHc); ax = H.ax_in(LMp, BUp, MWp, Hp)
        for u in units:                                                                      # each display unit = one exact translation of the tissue
            dx, dy = d['TR'][u]; su = cd_.piece.values.astype(int) == u
            iy_, ix_ = ((cd_.y.values[su] - y0d) // BIN_).astype(int), ((cd_.x.values[su] - x0d) // BIN_).astype(int)
            occ_u, occ_p = np.zeros((nyd, nxd), bool), np.zeros((nyd, nxd), bool); occ_u[iy_, ix_] = True; occ_p[iy_[inpc[su]], ix_[inpc[su]]] = True
            near = ndimage.distance_transform_edt(~occ_u, return_distances=False, return_indices=True)
            ok_ = d['MASK'][u] & occ_p[near[0], near[1]]                                     # the unit's tissue pixels whose nearest cell is a piece cell
            print('unit {}: {} cells, {} of them in the piece; translation ({:.0f}, {:.0f}) um; tissue pixels shown {} of {}'.format(u, int(su.sum()), int(inpc[su].sum()), dx, dy,
                  int(ok_.sum()), int(d['MASK'][u].sum())), flush=True)
            Eu = (d['EXT'][0] - dx, d['EXT'][1] - dx, d['EXT'][2] - dy, d['EXT'][3] - dy); XCu, YCu = d['XC'] - dx, d['YC'] - dy   # display -> tissue coordinates
            imC = ax.imshow(np.where(ok_, d['FCX'], np.nan), extent=Eu, interpolation='nearest', zorder=1, **H.CXKW)            # as cx_background()
            ax.contour(XCu, YCu, ok_.astype(float), levels=[0.5], colors='0.72', linewidths=0.6, zorder=2)
            ax.contour(XCu, YCu, (d['CXd'] & ok_).astype(float), levels=[0.5], colors='#333333', linewidths=U.LW_REG, zorder=6)   # CXCL12-High Region
        ap_ = cc_.ct.isin(U.APC).values; APC_S = U.APC_S_REF * (ipm / U.IN_PER_MM_REF) ** 2   # density-matched APC dots
        ax.scatter(cc_.x.values[ap_], cc_.y.values[ap_], s=APC_S, c=U.C_APCDOT, lw=0, rasterized=True, zorder=3); H.frame(ax, LF)
        ax.set_title('Antigen-Presenting Cells\n(DC, Macrophage, B Cell)', fontsize=15, pad=6)
        SBc = H.scalebar(ax, LF); LBc = fig.text(0.22 / FWc, (BUp + Hp / 2) / FHc, SEC, rotation=90, ha='center', va='center', fontsize=15)
        CBHc, CBTc = 0.40 * Hp, BUp + Hp - 0.75
        KX = [H.vbar(LMp + MWp + 0.08, CBTc, CBHc, imC, 'CXCL12\n(% of cells\nexpressing)', 'max', [t for t in np.arange(0, H.CXHI + 1e-9, 20)]),
              H.vkey(LMp + MWp + 0.08, CBTc - CBHc - 0.30, [Line2D([], [], ls='', marker='o', ms=4.5, mfc=U.C_APCDOT, mec='none', label='APC'),
                                                            Line2D([], [], color='#333333', lw=U.LW_REG, label='CXCL12-High\nRegion')])]
        trim_ticks(KX[:1]); fig.canvas.draw(); rd = fig.canvas.get_renderer()
        H.guards([('m', ax.get_window_extent(rd)), ('t', ax.title.get_window_extent(rd)), ('sb', SBc[0].get_window_extent(rd)), ('lbl', LBc.get_window_extent(rd))], KX, {'sb': 'm'})
        print('CXCL12 / APC map {:.2f} x {:.2f} in | APCs {} | CXCL12 scale 0-{:.1f}% (whole section) | scale bar {}'.format(MWp, Hp, int(ap_.sum()), H.CXHI, SBc[1]), flush=True)
        save(fig, OUTS['fovcx'], 0.06)


# ======================================================================================================================
# Section S9F_apc_cxcl12_quintiles: Treg Suppression Score by local APC quintile in CXCL12-Low vs CXCL12-High Regions, and its rise
# Panels: Supp 9F
# ======================================================================================================================
if _run('S9F_apc_cxcl12_quintiles'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/96b_quintile_cxcl12_bars_with_rise.py'   # working folder of this step (relative paths resolve from here)
    """Treg Suppression Score by local APC fraction in CXCL12-Low vs CXCL12-High Regions (CD resections, CosMx 6K; 11 sections, 6,554 strict Tregs).
Left: mean Treg Suppression Score (act_z: -D6 signature score, z within section) +- SE in CXCL12-Low vs CXCL12-High Regions per within-section quintile
  of the local APC fraction (apc_frac: APCs = DC, SELENOP+ Mac, CD9+ M2 Mac, Mo-Mac, B cell among all cells within 50 um); stars = BH q of per-quintile
  Mann-Whitney tests over the five quintiles (significant quintiles only). Regions = lower / upper within-section tertile of local CXCL12
  (lig_CXCL12_exp of 11_neighbourhoods_local_tables.py: stromal cells expected to express CXCL12 within 50 um, >= 3 required).
Right ('Rise'): mean score in the highest minus the lowest quintile, per region; error bar = +- SE of the rise = sqrt(SE_Q1^2 + SE_Q5^2).
  Difference P (CLUSTERED_P = True): OLS act_z ~ hi + top + hi:top + cz on Q1 + Q5 Tregs (top = highest quintile, cz = within-section z of log1p
  library size) with section-clustered SEs (t, G - 1 df), plus the number of sections whose per-section estimate has the pooled sign. With
  CLUSTERED_P = False the panel shows the LMM P of act_z ~ hi * top + cz + (1 | section) and the in-panel region x quintile interaction P.
Checks: the APC quintile values and interaction P reproduce c_apc_matched_cxcl12_high_vs_low_exp.csv and c_interaction_apc_by_cxcl12_countadj_exp.csv
  (11_neighbourhoods_local_tables.py); the clustered P and section agreement equal the rise-difference row of cd_spatial_results_all.csv
  (15_apc_local_models.py, section spatial_results_table). The same statistics are computed for the local DC fraction; FIGS = ['apc'] draws the APC figure.
Inputs: ligand_receptor_state/ (treg_local_ligand_receptor.csv.gz, treg_local_cxcl12_producer_exp.csv.gz, c_*_exp.csv), combined/all_samples.parquet,
  combined/raw_counts_strict_tregs.csv, each section's norm_log1p.h5ad (cell index), updated_figures/cd_spatial_results_all.csv.
Output: SuppFig9A-F_niche_cxcl12_apc/fig4_apc_quintile_cxcl12_bars_with_rise.png + _legend.txt (S96_OUT=<folder> writes elsewhere).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 RUN_SECTIONS=S9F_apc_cxcl12_quintiles python <this file>"""
    import os, sys, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from anndata.io import read_elem
    from scipy import stats
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    LR = f'{P4}/B_resection_analyses/ligand_receptor_state'
    CLUSTERED_P = True                                               # section-clustered difference P on the rise panel; no in-panel interaction P
    OUT = os.environ.get('S96_OUT', '/share/fsmresfiles/UC/revision_jenny/' + ('SuppFig9A-F_niche_cxcl12_apc' if CLUSTERED_P else 'Fig4A-F_treg_suppression_niches'))          # S96_OUT: optional output folder
    STEMS = {'dc': f'{OUT}/fig4_dc_quintile_cxcl12_bars_with_rise', 'apc': f'{OUT}/fig4_apc_quintile_cxcl12_bars_with_rise'}
    FIGS = ['apc']                                                   # figures drawn (APC version)
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50
    APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']; DC = ['DC']
    COL_HIGH, COL_LOW = '#F7B325', '#1a9850'                          # CXCL12-High orange (as the CXCL12-High ellipses of Fig 5E), CXCL12-Low green
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    def close(a, b, rtol=5e-3, atol=1e-4): return bool(np.isclose(a, b, rtol=rtol, atol=atol))
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
    # ---------------- Tregs and local DC count ----------------
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    T = T.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left'); assert T.n_exp_producers_50um.notna().all()
    N_TREG, N_SEC, N_NOPROD = len(T), T.section.nunique(), int(T.lig_CXCL12_exp.isna().sum())
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
    for c in ['dc_n', 'apc_n2', 'nb_n2']: T[c] = np.nan
    for s_, f in FOLDER.items():
        sub = sa[sa.section == s_]
        with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h: ob = read_elem(h['obs'])
        nm = pd.Index(s_ + '_' + ob.index.astype(str)); cs = sub.loc[sub.index.intersection(nm)]
        tg = T.index[T.section == s_]
        nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg, ['x', 'y']].values, r=RAD)
        ctv = cs.ct.values; dm, am = np.isin(ctv, DC), np.isin(ctv, APC)
        T.loc[tg, 'nb_n2'] = [len(i) for i in nb]; T.loc[tg, 'dc_n'] = [dm[i].sum() for i in nb]; T.loc[tg, 'apc_n2'] = [am[i].sum() for i in nb]
    del sa; gc.collect()
    assert (T.nb_n == T.nb_n2).all() and (T.apc_n == T.apc_n2).all(), 'neighbourhood or APC count differs from treg_local_ligand_receptor.csv.gz'
    T['dc_frac'] = T.dc_n / T.nb_n.replace(0, np.nan)
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    assert len(cnt) == len(T)
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()   # Treg order of raw_counts_strict_tregs.csv (rank(method='first') ties depend on it)
    zs = lambda v: (v - v.mean()) / v.std()
    T['cz'] = T.groupby('section').cz.transform(zs)
    T['apc_q'] = T.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    T['dc_q'] = T.groupby('section').dc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    assert (T.groupby('section').apc_q.nunique() == 5).all() and (T.groupby('section').dc_q.nunique() == 5).all() and T.dc_q.notna().all()
    T['t_exp'] = T.groupby('section').lig_CXCL12_exp.transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))
    def lmm(d, f, terms):
        m = smf.mixedlm(f, d, groups=d.section).fit(reml=True, method='lbfgs')
        return {t: dict(coef=m.params[t], se=m.bse[t], p=m.pvalues[t], n=len(d)) for t in terms}
    def per_quintile(qcol, fcol):
        sub = T[T.t_exp.isin(['low', 'high']) & T[qcol].notna()]
        out = []
        for q, g in sub.groupby(qcol):
            hi, lo = g[g.t_exp == 'high'].act_z.dropna(), g[g.t_exp == 'low'].act_z.dropna()
            if len(hi) < 20 or len(lo) < 20: continue
            u = stats.mannwhitneyu(hi, lo); sd = np.sqrt((hi.var() * (len(hi) - 1) + lo.var() * (len(lo) - 1)) / (len(hi) + len(lo) - 2))
            out.append(dict(quintile=int(q) + 1, mean_frac=g[fcol].mean(), n_high=len(hi), n_low=len(lo),
                            suppression_cxcl12_high=hi.mean(), suppression_cxcl12_low=lo.mean(),
                            sem_high=hi.std() / np.sqrt(len(hi)), sem_low=lo.std() / np.sqrt(len(lo)), delta_sd=(hi.mean() - lo.mean()) / sd, p=u.pvalue))
        out = pd.DataFrame(out); out['q'] = multipletests(out.p, method='fdr_bh')[1]; return out
    def interaction(qcol):
        d = T[T.t_exp.isin(['low', 'high']) & T[qcol].notna()].copy(); d['hi'] = (d.t_exp == 'high').astype(float)
        d['qz'] = (d[qcol] - d[qcol].mean()) / d[qcol].std(); d = d[['section', 'act_z', 'hi', 'qz', 'cz']].dropna().rename(columns={'act_z': 'y'})
        return lmm(d, 'y ~ hi * qz + cz', ['hi', 'qz', 'hi:qz', 'cz'])
    def did(qcol):   # difference-in-differences, lowest (Q1) and highest (Q5) quintile Tregs only
        d = T[T.t_exp.isin(['low', 'high']) & T[qcol].isin([0, 4])].copy()
        d['hi'] = (d.t_exp == 'high').astype(float); d['top'] = (d[qcol] == 4).astype(float)
        d = d[['section', 'act_z', 'hi', 'top', 'cz']].dropna().rename(columns={'act_z': 'y'})
        return lmm(d, 'y ~ hi * top + cz', ['hi', 'top', 'hi:top', 'cz'])
    def rises(Q):
        out = {}
        for k in ['low', 'high']:
            out[k] = (Q[f'suppression_cxcl12_{k}'][5] - Q[f'suppression_cxcl12_{k}'][1], np.sqrt(Q[f'sem_{k}'][1] ** 2 + Q[f'sem_{k}'][5] ** 2))
        return out
    # ---------------- APC quintiles (asserted against the tables of 11_neighbourhoods_local_tables.py) ----------------
    CQ = pd.read_csv(f'{LR}/c_apc_matched_cxcl12_high_vs_low_exp.csv').set_index('apc_quintile')
    IX = pd.read_csv(f'{LR}/c_interaction_apc_by_cxcl12_countadj_exp.csv').set_index('term')
    QA = per_quintile('apc_q', 'apc_frac').set_index('quintile'); assert list(QA.index) == list(CQ.index) == [1, 2, 3, 4, 5]
    for qq in QA.index:
        a_, b_ = QA.loc[qq], CQ.loc[qq]
        assert a_.n_high == b_.n_high and a_.n_low == b_.n_low, qq
        for c in ['suppression_cxcl12_high', 'suppression_cxcl12_low', 'delta_sd', 'p', 'q']: assert close(a_[c], b_[c]), (qq, c)
        assert close(a_.mean_frac, b_.mean_apc_frac), qq
    IA = interaction('apc_q')
    assert close(IA['hi:qz']['coef'], IX.coef['hi:apcz']) and close(IA['hi:qz']['p'], IX.p['hi:apcz']) and IA['hi:qz']['n'] == int(IX.n['hi:apcz'])
    print('check passed: APC bar values and interaction P = {:.3g} reproduced'.format(IA['hi:qz']['p']), flush=True)
    # ---------------- DC quintiles ----------------
    QD = per_quintile('dc_q', 'dc_frac').set_index('quintile'); assert list(QD.index) == [1, 2, 3, 4, 5]
    ID = interaction('dc_q')
    # ---------------- rise and difference-in-differences ----------------
    RES = {}
    for key, Q, I, qcol in [('dc', QD, ID, 'dc_q'), ('apc', QA, IA, 'apc_q')]:
        rz, dd = rises(Q), did(qcol)
        RES[key] = dict(Q=Q, I=I, rise=rz, did=dd)
        print('{}: rise CXCL12-Low {:+.3f} +- {:.3f} | CXCL12-High {:+.3f} +- {:.3f} | difference {:+.3f} | DiD hi:top coef {:+.3f} P = {:.3g} (n {:,})'.format(
            key.upper(), rz['low'][0], rz['low'][1], rz['high'][0], rz['high'][1], rz['high'][0] - rz['low'][0], dd['hi:top']['coef'], dd['hi:top']['p'],
            dd['hi:top']['n']), flush=True)
    # CLUSTERED_P: section-clustered P and section sign agreement of the APC difference-in-differences (as in 15_apc_local_models.py)
    if CLUSTERED_P:
        d = T[T.t_exp.isin(['low', 'high']) & T.apc_q.isin([0, 4])].copy()
        d['hi'] = (d.t_exp == 'high').astype(float); d['top'] = (d.apc_q == 4).astype(float)
        d = d[['section', 'act_z', 'hi', 'top', 'cz']].dropna(); d['hixtop'] = d.hi * d.top
        FC = 'act_z ~ hi + top + hixtop + cz'
        o = smf.ols(FC, d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]), use_t=True)
        CL_P = 2 * stats.t.sf(abs(o.params['hixtop'] / o.bse['hixtop']), d.section.nunique() - 1)
        sg = np.sign(RES['apc']['did']['hi:top']['coef']); AG, NS = 0, 0
        for _, g in d.groupby('section'):
            if g.hixtop.std() == 0 or len(g) < 5: continue
            NS += 1; AG += int(np.sign(smf.ols(FC, g).fit().params['hixtop']) == sg)
        CS = pd.read_csv('/share/fsmresfiles/UC/revision_jenny/updated_figures/cd_spatial_results_all.csv')
        cr = CS[(CS.block == 'APC x CXCL12') & CS.feature.str.startswith('Rise')]; assert len(cr) == 1; cr = cr.iloc[0]
        assert int(cr.n_tregs) == len(d) == RES['apc']['did']['hi:top']['n'] and close(cr.p, RES['apc']['did']['hi:top']['p'])
        assert np.isclose(CL_P, cr.clustered_p, rtol=5e-3, atol=0) and cr.sections_same_sign == f'{AG}/{NS}', (CL_P, cr.clustered_p, AG, NS)
        assert f'{CL_P:.2g}' == '0.025' and (AG, NS) == (9, 11)
        print(f'check passed: APC rise difference clustered P = {CL_P:.3g} (CSV {cr.clustered_p}), {AG}/{NS} sections (CSV {cr.sections_same_sign})', flush=True)
    # ---------------- figures ----------------
    def pfmt(p):
        if p >= 0.001: return f'{p:.2g}'
        ex = int(np.floor(np.log10(p))); m = p / 10 ** ex
        if m >= 9.95: m, ex = 1.0, ex + 1
        return '{:.1f} × 10$^{{{}}}$'.format(m, ex)
    def stars(q): return 'n.s.' if q >= 0.05 else ('*' if q >= 0.01 else ('**' if q >= 0.001 else '***'))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 15, 'axes.linewidth': 1.2})
    EKW = dict(ecolor='k', elinewidth=1.2, capsize=3, capthick=1.2)
    LAB = {'dc': dict(ab='DC', title='Treg Suppression Score\nby Local DC Fraction',
                      xlab='Dendritic Cells within 50 μm\n(Quintile of Local DC Fraction)'),
           'apc': dict(ab='APC', title='Treg Suppression Score\nby Local APC Fraction',
                       xlab='Antigen-Presenting Cells within 50 μm\n(Quintile of Local APC Fraction)')}
    for key in FIGS:
        Q, INT_P, rz, dd = RES[key]['Q'], RES[key]['I']['hi:qz']['p'], RES[key]['rise'], RES[key]['did']
        L = LAB[key]
        fig = plt.figure(figsize=(8.7, 5.4))
        gs = fig.add_gridspec(1, 2, width_ratios=[4.6, 1.75], wspace=0.48)
        ax, axr = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
        # ---- left: quintile bars ----
        C = Q.reset_index(); xs = np.arange(len(C)); w = 0.38
        elo, ehi = C.sem_low, C.sem_high                              # mean +- 1 SE
        ax.bar(xs - w / 2, C.suppression_cxcl12_low, w, yerr=elo, color=COL_LOW, edgecolor='k', lw=0.7, label='CXCL12-Low Region', error_kw=EKW, zorder=2)
        ax.bar(xs + w / 2, C.suppression_cxcl12_high, w, yerr=ehi, color=COL_HIGH, edgecolor='k', lw=0.7, label='CXCL12-High Region', error_kw=EKW, zorder=2)
        ymax = max((C.suppression_cxcl12_low + elo).max(), (C.suppression_cxcl12_high + ehi).max(), 0.0)
        ymin = min((C.suppression_cxcl12_low - elo).min(), (C.suppression_cxcl12_high - ehi).min(), 0.0)
        for xi, r in zip(xs, C.itertuples()):
            top = max(r.suppression_cxcl12_low + r.sem_low, r.suppression_cxcl12_high + r.sem_high, 0.0)
            if r.q < 0.05: ax.text(xi, top + 0.02 * (ymax - ymin), stars(r.q), ha='center', va='bottom', fontsize=16)
        ax.axhline(0, color='0.3', lw=1, zorder=1); ax.set_xticks(xs)
        TL = [f'{v*100:.0f}%' for v in C.mean_frac]; assert len(set(TL)) == len(TL), TL
        ax.set_xticklabels(TL, fontsize=13.5)
        ax.set_ylim(ymin - 0.06 * (ymax - ymin), ymax + 0.30 * (ymax - ymin))
        ax.set_xlabel(L['xlab'], fontsize=14.5)
        ax.set_ylabel('Treg Suppression Score (z)', fontsize=15, labelpad=1)
        ax.set_title(L['title'], fontsize=16, fontweight='normal')
        ax.legend(frameon=False, fontsize=12.5, loc='upper left', handlelength=1.1, borderpad=0.1)
        if not CLUSTERED_P: ax.text(0.03, 0.80, 'CXCL12 Region × {} Quintile\nInteraction $P$ = '.format(L['ab']) + pfmt(INT_P) + ('' if INT_P < 0.05 else ' (n.s.)'),
                transform=ax.transAxes, ha='left', va='top', fontsize=12.5)
        # ---- right: rise, lowest -> highest quintile ----
        rv = [rz['low'][0], rz['high'][0]]; re_ = [rz['low'][1], rz['high'][1]]
        axr.bar([0, 1], rv, 0.62, yerr=re_, color=[COL_LOW, COL_HIGH], edgecolor='k', lw=0.7, error_kw=EKW, zorder=2)
        rtop = max(v + e for v, e in zip(rv, re_)); rbot = min(0.0, min(v - e for v, e in zip(rv, re_)))
        rng = rtop - rbot; yb = rtop + 0.07 * rng
        axr.plot([0, 0, 1, 1], [yb - 0.03 * rng, yb, yb, yb - 0.03 * rng], color='k', lw=1.1)
        pd_ = CL_P if CLUSTERED_P else dd['hi:top']['p']               # section-clustered P + sections agreeing
        PTXT = 'Difference\n$P$ = ' + pfmt(pd_) + ('' if pd_ < 0.05 else ' (n.s.)') + (f';\n{AG}/{NS} sections' if CLUSTERED_P else '')
        PT = axr.text(0.5, yb + 0.03 * rng, PTXT, ha='center', va='bottom', fontsize=12.5)
        axr.set_ylim(rbot - 0.04 * rng if rbot < 0 else 0, rtop + 0.42 * rng)
        axr.axhline(0, color='0.3', lw=1, zorder=1)
        axr.set_xlim(-0.6, 1.6); axr.set_xticks([0, 1]); axr.set_xticklabels(['Low', 'High'], fontsize=13.5)
        axr.set_xlabel('CXCL12\nRegion', fontsize=14.5)
        axr.set_ylabel('Rise in Treg Suppression Score,\nLowest → Highest Quintile (SD)', fontsize=13.5, labelpad=2)
        axr.set_title('Rise', fontsize=16, fontweight='normal')
        for a in (ax, axr):
            for sp in ['top', 'right']: a.spines[sp].set_visible(False)
        fig.canvas.draw(); RD = fig.canvas.get_renderer()
        # overlap guards: P text vs the rise panel's y tick labels and its own y axis; rise y label vs the left panel
        pb = PT.get_window_extent(RD)
        assert pb.x0 > axr.get_window_extent(RD).x0 and not any(pb.overlaps(t.get_window_extent(RD)) for t in axr.get_yticklabels() if t.get_text()), 'P text overlaps'
        assert axr.yaxis.label.get_window_extent(RD).x0 > ax.get_tightbbox(RD).x1, 'rise y label overlaps the left panel'
        assert CLUSTERED_P or not ax.texts[-1].get_window_extent(RD).overlaps(ax.get_legend().get_window_extent(RD)), 'interaction text overlaps legend'
        for a in (ax, axr):
            assert all(t.get_fontweight() in ('normal', 400) for t in [a.title, a.xaxis.label, a.yaxis.label] + a.texts)
        with open(f'{STEMS[key]}.png', 'wb') as fh: fig.savefig(fh, format='png', dpi=300, bbox_inches='tight', pad_inches=0.06)
        plt.close(fig); print('wrote', STEMS[key] + '.png (PNG only)', flush=True)
        # ---- legend ----
        mn = lambda v: f'{v:+.2f}'.replace('-', '−')
        frac = ('local dendritic cell (DC) fraction, DCs among all cells within 50 μm' if key == 'dc' else
                'local antigen-presenting cell (APC) fraction, APCs (DCs, macrophages and B cells) among all cells within 50 μm')
        NSIG = int((Q.q < 0.05).sum())
        startxt = ('Stars, Benjamini-Hochberg q of per-quintile Mann-Whitney tests across the five quintiles (* < 0.05, ** < 0.01, *** < 0.001), shown '
               'for significant quintiles only.' if NSIG else
                   'Per-quintile Mann-Whitney tests, Benjamini-Hochberg-adjusted across the five quintiles, reached q < 0.05 in no quintile (highest '
               'quintile, P = {:.2g}, q = {:.2g}), so no stars are shown.'.format(Q.p[5], Q.q[5]))
        leg = ('Treg Suppression Score by Local {ab} Fraction. Left, Treg Suppression Score (negated D6 pseudotime signature, z-scored within section; '
           'mean ± SE) of strict Tregs in Crohn\'s disease resections (CosMx 6K, {ns} sections) in CXCL12-Low (green, n = {nl:,}) and CXCL12-High '
           'Regions (orange, n = {nh:,}), across within-section quintiles of the {frac} (x axis, mean fraction per quintile). Regions are the lower '
           'and upper within-section tertiles of local CXCL12, the mean within-section CXCL12 z-score of stromal cells expected to express CXCL12 '
           '(fibroblasts, myofibroblasts, pericytes and endothelial cells) within 50 μm, requiring at least three such cells ({nnp:,} Tregs with '
           'fewer were excluded). The in-panel P is the CXCL12 region-by-{ab} quintile interaction from a library-size-adjusted linear mixed model '
           'with a section random intercept. {startxt} Right, rise in the score from the lowest to the highest {ab} quintile (mean in the highest '
           'minus mean in the lowest quintile; ± SE, the square root of the summed squared SEs): {rl} ± {sl} SD in CXCL12-Low and {rh} ± {sh} SD in '
           'CXCL12-High Regions. Difference P, CXCL12 region-by-quintile term of a difference-in-differences linear mixed model on Tregs of the lowest '
           'and highest quintiles only ({nd:,} Tregs), adjusted for library size, with a section random intercept.').format(
               ab=L['ab'], ns=N_SEC, nl=int(Q.n_low.sum()), nh=int(Q.n_high.sum()), frac=frac, nnp=N_NOPROD, startxt=startxt,
               rl=mn(rz['low'][0]), sl=f"{rz['low'][1]:.2f}", rh=mn(rz['high'][0]), sh=f"{rz['high'][1]:.2f}", nd=dd['hi:top']['n'])
        if CLUSTERED_P:   # no in-panel interaction P; Difference P is the section-clustered P
            o1 = ('The in-panel P is the CXCL12 region-by-{ab} quintile interaction from a library-size-adjusted linear mixed model with a section '
              'random intercept. ').format(ab=L['ab'])
            o2 = 'linear mixed model on Tregs of the lowest and highest quintiles only ({:,} Tregs), adjusted for library size, with a section random intercept.'.format(dd['hi:top']['n'])
            assert leg.count(o1) == 1 and leg.count(o2) == 1
            leg = leg.replace(o1, '').replace(o2, ('model on Tregs of the lowest and highest quintiles only ({:,} Tregs), adjusted for library size, with standard errors '
              'clustered by section (t distribution, {} df); {}/{} sections, sections whose per-section estimate has the sign of the pooled '
              'estimate.').format(dd['hi:top']['n'], N_SEC - 1, AG, NS))
        for w_ in ['ctivation', 'erritor', 'lunting', 'olocaliz', 'olocalis', 'uiescent', 'serial', '(blue', '(red', '(yellow']: assert w_ not in leg, w_
        for w_ in ['neighbour', 'colour', 'signalling', 'standardis', 'characteris', 'behaviour', 'centre', 'favour', 'analys', 'tumour', 'modell']: assert w_ not in leg.lower(), w_
        with open(f'{STEMS[key]}_legend.txt', 'w') as fh: fh.write(leg + '\n')
        print(chr(10) + leg + chr(10), flush=True)


# ======================================================================================================================
# Section S9E_cxcl12_by_celltype: CXCL12 expression by cell type (Supp 9E); with RUN_FIG2=1 also the APC abundance per CN table
# Panels: Supp 9E; APC abundance per CN table (input of 17_neighbourhood_apc_mhcii.py, Supp 9C)
# Note: ends with sys.exit(0) unless RUN_FIG2=1, so it is the last section of this file
# ======================================================================================================================
if _run('S9E_cxcl12_by_celltype'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/80_supp_cxcl12_and_apc_abundance.py'   # working folder of this step (relative paths resolve from here)
    """CXCL12 expression by cell type and APC abundance per cellular neighbourhood (CD resections, CosMx 6K; 11 sections, 4 patients).
FIGURE 1 (Supp 9E) supp_cxcl12_expression_by_cell_type: per section, cells of norm_log1p.h5ad present in combined/all_samples.parquet, cell_type_general
  from all_samples.parquet, expressing = CXCL12 > 0. Per section and cell type: % expressing and mean CXCL12 in expressing cells; mean +- 95% CI (t)
  across sections with >= MIN_CELLS cells of the type; per-section rank; OGN+RSPO3+ Fib vs each other type: paired t-test on the per-section values
  over shared sections, BH over all comparisons. SELENOP+ Fib is excluded from the figure and the statistics (SELENOP+ Mac is shown and tested); the
  All Cells reference keeps every cell. 'Absoprtive Epi' (spelling of the annotation; BEST4+ cells in the resections) is shown as 'BEST4+ Epi'.
FIGURE 2 / TABLE (RUN_FIG2=1) supp_apc_abundance_by_neighbourhood: APC abundance per cellular neighbourhood (CN_k10 of
  cellular_neighborhoods_v2/cn_per_cell.parquet, 1-based), all cells (not Treg-centred). APCs = cell_type_general DC, B cell, CD9+ M2 Mac, Mo-Mac,
  SELENOP+ Mac (SELENOP+ Fib is not an APC), asserted to match the cn_per_cell cell_type vocabulary (DC, B, CD9 Mac, Mo-Mac + Inf Mo-Mac, Mac S +
  Mac S_M_P + Mac S_M_S). Per section x CN (>= MIN_CN_CELLS cells): APC % = 100 x APCs / all cells, also per APC group (DC, B cell,
  Macrophage = SELENOP+ Mac + CD9+ M2 Mac + Mo-Mac). Linear mixed model APC % ~ C(CN, reference CN10) + (1 | section) (REML), BH over the 9
  contrasts; number of sections in which CN10 is lower than the CN. CN labels as the Fig 5B CN map. The CSV is an input of
  17_neighbourhood_apc_mhcii.py (sections cn_apc_mhcii_table and cn10_apc_matched_score).
Outputs: SuppFig9A-F_niche_cxcl12_apc/supp_cxcl12_expression_by_cell_type.{png,csv} + _legend.txt; with RUN_FIG2=1 also
  panels/fig4_and_related/archive/supp_apc_abundance_by_neighbourhood.{png,csv} + _legend.txt (on-disk folder name).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 [RUN_FIG2=1] RUN_SECTIONS=S9E_cxcl12_by_celltype python <this file>
     (CACHE80=<path.pkl> stores / re-uses the per-cell CXCL12 + CN table)."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
    from anndata.io import read_elem
    from scipy import sparse, stats as sst
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    from statsmodels.tools.sm_exceptions import ConvergenceWarning; warnings.simplefilter('ignore', ConvergenceWarning)
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUTA = f'{P4}/A_fig4_panels'
    LR = f'{P4}/B_resection_analyses/ligand_receptor_state'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    STEM1 = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cxcl12_expression_by_cell_type'; STEM2 = f'{P4}/archive/supp_apc_abundance_by_neighbourhood'   # STEM2: APC abundance per CN (input of 17_neighbourhood_apc_mhcii.py)
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    EXCL = ['SELENOP+ Fib']                                      # not shown, not tested (SELENOP+ Mac is shown)
    OGN = 'OGN+RSPO3+ Fib'; HI = '#b2182b'; HIL = '#e7a3a3'; GR = '#7f7f7f'; GRL = '#c4c4c4'
    DISP = {'OGN+RSPO3+ Fib': 'OGN⁺RSPO3⁺ Fib', 'ADAMDEC1+ Fib': 'ADAMDEC1⁺ Fib', 'PLA2G2A+ ECM Fib': 'PLA2G2A⁺ ECM Fib', 'CD9+ M2 Mac': 'CD9⁺ M2 Mac',
            'Absoprtive Epi': 'BEST4⁺ Epi', 'Gd T': 'γδ T', 'SELENOP+ Mac': 'SELENOP⁺ Mac'}
    APC = ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac']   # SELENOP+ Mac is an APC
    APC_SHARED = {'DC': 'DC', 'B': 'B cell', 'CD9 Mac': 'CD9+ M2 Mac', 'Mo-Mac': 'Mo-Mac', 'Inf Mo-Mac': 'Mo-Mac',
                  'Mac S': 'SELENOP+ Mac', 'Mac S_M_P': 'SELENOP+ Mac', 'Mac S_M_S': 'SELENOP+ Mac'}
    APC_GRP = {'DC': 'DC', 'B cell': 'B cell', 'CD9+ M2 Mac': 'Macrophage', 'Mo-Mac': 'Macrophage', 'SELENOP+ Mac': 'Macrophage'}
    APCG = ['DC', 'B cell', 'Macrophage']                       # composition groups (all macrophages in one category)
    APC_COL = {'DC': '#08519c', 'B cell': '#6baed6', 'Macrophage': '#31a354'}
    MIN_CELLS = 20; MIN_CN_CELLS = 100
    sig3 = lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    def ci95(v):
        v = np.asarray(v, float); v = v[np.isfinite(v)]; n = len(v)
        if n < 2: return (v.mean() if n else np.nan), np.nan, np.nan, n
        h = sst.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n); return v.mean(), v.mean() - h, v.mean() + h, n
    def wtxt(path, s):
        with open(path, 'w') as fh: fh.write(s)
    # ======================================================= per-cell table (CXCL12 per cell, plus CN)
    CACHE = os.environ.get('CACHE80')
    if CACHE and os.path.exists(CACHE):
        C = pd.read_pickle(CACHE)
    else:
        sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
        sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
        sa['ct'] = sa.cell_type_general.astype(str); sa = sa[['section', 'ct']]
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['section', 'cell_type', 'CN_k10', '__index_level_0__']).to_pandas()
        parts = []
        for s_, f in FOLDER.items():
            sub = sa[sa.section == s_]
            with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
                ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
            nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsc() if sparse.issparse(X) else sparse.csc_matrix(X)
            keep = nm.isin(sub.index)
            e = pd.Series(np.asarray(X[:, v.get_loc('CXCL12')].toarray()).ravel()[keep], index=nm[keep]); del X, ob; gc.collect()
            e = e.loc[e.index.intersection(sub.index)]                    # cells present in all_samples.parquet
            cs = sub.loc[e.index].copy(); cs['expr'] = e.values
            parts.append(cs); print(s_, 'cells', len(cs), flush=True)
        C = pd.concat(parts); del parts, sa; gc.collect()
        ok = C.index.isin(cn.index); print('cells with a CN:', ok.sum(), '/', len(C), '| cn_per_cell rows:', len(cn), flush=True)
        cn = cn.reindex(C.index)
        assert (cn.section.dropna().astype(str).values == C.section[cn.section.notna()].values).all()
        C['CN'] = cn.CN_k10.values; C['cn_cell_type'] = cn.cell_type.values
        if CACHE: pd.to_pickle(C, CACHE)
    C['pos'] = C.expr > 0
    # ======================================================= FIGURE 1: CXCL12 by cell type
    pool = C.groupby('ct').agg(n_cells=('pos', 'size'), n_expressing=('pos', 'sum'))
    pool['pct_expressing'] = 100 * pool.n_expressing / pool.n_cells
    pool['mean_expr_in_expressing'] = C[C.pos].groupby('ct').expr.mean()
    tot_pct = 100 * C.pos.mean()
    g = C.groupby(['section', 'ct'])
    PS = pd.DataFrame({'n_cells': g.size(), 'n_expressing': g.pos.sum()})
    PS['pct_expressing'] = 100 * PS.n_expressing / PS.n_cells
    PS['mean_expr_in_expressing'] = C[C.pos].groupby(['section', 'ct']).expr.mean().reindex(PS.index)
    PS = PS.reset_index(); PS = PS[~PS.ct.isin(EXCL)]
    PSu = PS[PS.n_cells >= MIN_CELLS].copy()
    sec_all = (100 * C.groupby('section').pos.mean()); all_m, all_lo, all_hi, _ = ci95(sec_all.values)
    pv = PSu.pivot(index='section', columns='ct', values='pct_expressing')
    rank = pv.rank(axis=1, ascending=False, method='min')
    CTS = [c for c in pool.index if c not in EXCL]
    rows = []
    for ct in CTS:
        m, lo, hi, n = ci95(pv[ct]) if ct in pv else (np.nan,) * 3 + (0,)
        me = PSu[PSu.ct == ct].mean_expr_in_expressing
        m2, lo2, hi2, n2 = ci95(me)
        rows.append(dict(cell_type_general=ct, label=DISP.get(ct, ct), n_cells=int(pool.loc[ct, 'n_cells']), n_expressing=int(pool.loc[ct, 'n_expressing']),
                         pooled_pct_expressing=pool.loc[ct, 'pct_expressing'], n_sections=n, mean_pct_expressing=m, ci_low=lo, ci_high=hi,
                         fold_over_all_cells=m / all_m, pooled_mean_expr_in_expressing=pool.loc[ct, 'mean_expr_in_expressing'],
                         mean_expr_in_expressing=m2, expr_ci_low=lo2, expr_ci_high=hi2,
                         n_sections_ranked_top=int((rank[ct] == 1).sum()) if ct in rank else 0,
                         median_rank_across_sections=rank[ct].median() if ct in rank else np.nan))
    T1 = pd.DataFrame(rows).sort_values('mean_pct_expressing', ascending=False).reset_index(drop=True)
    # OGN vs each other type, paired over shared sections
    tst = []
    for ct in T1.cell_type_general:
        if ct == OGN: tst.append((np.nan, np.nan, np.nan)); continue
        d = pv[[OGN, ct]].dropna(); t, p = sst.ttest_rel(d[OGN], d[ct])
        tst.append(((d[OGN] - d[ct]).mean(), p, f'{int((d[OGN] > d[ct]).sum())}/{len(d)}'))
    T1['ogn_minus_type_pct'] = [x[0] for x in tst]; T1['p_paired_t_vs_ogn'] = [x[1] for x in tst]
    T1['n_sections_ogn_higher'] = [x[2] for x in tst]
    m_ = T1.p_paired_t_vs_ogn.notna(); T1['q_vs_ogn'] = np.nan
    T1.loc[m_, 'q_vs_ogn'] = multipletests(T1.loc[m_, 'p_paired_t_vs_ogn'], method='fdr_bh')[1]
    T1['rank'] = np.arange(1, len(T1) + 1)
    lv = PSu.pivot(index='section', columns='ct', values='mean_expr_in_expressing'); pl = []
    for ct in T1.cell_type_general:
        if ct == OGN: pl.append(np.nan); continue
        d = lv[[OGN, ct]].dropna(); pl.append(sst.ttest_rel(d[OGN], d[ct])[1])
    T1['p_paired_t_level_vs_ogn'] = pl; T1['q_level_vs_ogn'] = np.nan
    T1.loc[m_, 'q_level_vs_ogn'] = multipletests(T1.loc[m_, 'p_paired_t_level_vs_ogn'], method='fdr_bh')[1]
    lv_hi = T1[(T1.q_level_vs_ogn < 0.05) & (T1.mean_expr_in_expressing < T1.mean_expr_in_expressing.iloc[0])]
    lv_ns = T1[(T1.q_level_vs_ogn >= 0.05)]
    print('level: OGN higher (q<0.05) than', len(lv_hi), 'types; ns vs', list(lv_ns.cell_type_general), flush=True)
    assert T1.cell_type_general.iloc[0] == OGN
    n_top = int(T1.loc[T1.cell_type_general == OGN, 'n_sections_ranked_top'].iloc[0]); n_sec = int(T1.loc[T1.cell_type_general == OGN, 'n_sections'].iloc[0])
    second = T1.iloc[1]
    print('== FIGURE 1 table'); print(T1.drop(columns=['label']).round(3).to_string()); print('OGN top in', n_top, '/', n_sec, 'sections', flush=True)
    ALLR = pd.DataFrame([dict(cell_type_general='ALL CELLS (reference; includes SELENOP+ Fib)', label='All Cells', n_cells=len(C), n_expressing=int(C.pos.sum()),
                              pooled_pct_expressing=tot_pct, n_sections=len(sec_all), mean_pct_expressing=all_m, ci_low=all_lo, ci_high=all_hi, fold_over_all_cells=1.0,
                              pooled_mean_expr_in_expressing=C.expr[C.pos].mean())])
    PSo = PS.rename(columns={'ct': 'cell_type_general'}).assign(row_type='section', used_in_mean=lambda d: d.n_cells >= MIN_CELLS)
    OUT1 = pd.concat([T1.assign(row_type='summary'), ALLR.assign(row_type='summary'), PSo], ignore_index=True)
    lead = ['row_type', 'rank', 'cell_type_general', 'label', 'section']
    OUT1 = OUT1[lead + [c for c in OUT1 if c not in lead]].apply(sig3)
    OUT1.to_csv(f'{STEM1}.csv', index=False); print('wrote', f'{STEM1}.csv', flush=True)
    # ---- figure 1
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    n = len(T1); Y = {c: i for i, c in enumerate(T1.cell_type_general)}
    FH = 0.315 * n + 1.9
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.9, FH), sharey=True, gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.13))
    rng = np.random.default_rng(1)
    for ct, i in Y.items():
        hi_ = ct == OGN; col, colL = (HI, HIL) if hi_ else (GR, GRL)
        r = T1.iloc[i]
        v = PSu[PSu.ct == ct]
        for a_, val, m, lo, hi in [(ax, v.pct_expressing, r.mean_pct_expressing, r.ci_low, r.ci_high),
                                   (bx, v.mean_expr_in_expressing, r.mean_expr_in_expressing, r.expr_ci_low, r.expr_ci_high)]:
            a_.barh(i, m, height=0.72, color=colL, edgecolor='none', zorder=1)
            a_.scatter(val, i + rng.uniform(-0.2, 0.2, len(val)), s=10, color=col, alpha=0.75, linewidth=0, zorder=2)
            a_.plot([lo, hi], [i, i], color='black' if not hi_ else HI, lw=1.8, zorder=3, solid_capstyle='butt')
    ax.axvline(all_m, color='#252525', ls='--', lw=1.4, zorder=0)
    ax.text(all_m + 1, n - 0.6, 'All Cells', ha='left', va='bottom', fontsize=14, color='#252525')
    ax.set_yticks(range(n)); ax.set_yticklabels(T1.label, fontsize=13.5); ax.set_ylim(n - 0.4, -0.6)
    for tl, ct in zip(ax.get_yticklabels(), T1.cell_type_general):
        if ct == OGN: tl.set_color(HI); tl.set_fontsize(15)   # colour only, no bold
    ax.set_xlabel('Cells Expressing CXCL12 (%)', fontsize=15.5); bx.set_xlabel('Mean CXCL12\n(log-normalised)', fontsize=15.5)
    xm = np.nanmax(PSu.pct_expressing); ax.set_xlim(0, xm * 1.04); ax.set_xticks(np.arange(0, xm * 1.04, 20))
    bx.set_xlim(0, np.nanmax(PSu.mean_expr_in_expressing) * 1.04); bx.set_xticks([0, 0.5, 1.0]); bx.set_xticklabels(['0', '0.5', '1.0'])
    for a_ in (ax, bx): a_.spines[['top', 'right']].set_visible(False); a_.tick_params(axis='x', labelsize=13.5); a_.margins(y=0)
    bx.tick_params(axis='y', length=0)
    ax.set_title('% Expressing CXCL12', loc='left', fontsize=16, pad=5)
    bx.set_title('Level in CXCL12-\nExpressing Cells', loc='left', fontsize=16, pad=5)
    # significance (inside the plot): OGN vs every other type
    qmax = T1.q_vs_ogn.max(); nsig = int((T1.q_vs_ogn < 0.05).sum()); ncomp = int(T1.q_vs_ogn.notna().sum())
    hi_pct = T1[(T1.q_vs_ogn < 0.05) & (T1.ogn_minus_type_pct > 0)]; n_hi_pct = len(hi_pct); q_hi_max = hi_pct.q_vs_ogn.max()
    ns_pct = T1[T1.q_vs_ogn >= 0.05]
    print('% expressing: OGN higher (q<0.05) than', n_hi_pct, '/', ncomp, 'types; max q', q_hi_max, '; ns vs', list(ns_pct.cell_type_general), flush=True)
    o = T1.iloc[0]
    TXA = (f'OGN⁺RSPO3⁺ Fib: {o.mean_pct_expressing:.0f}%\n({o.mean_pct_expressing / all_m:.1f}× all cells)\n'
       f'Top in {n_top}/{n_sec} sections\nHigher than {n_hi_pct}/{ncomp} types,\nq {"< 0.001" if q_hi_max < 0.001 else f"≤ {q_hi_max:.2g}"}')
    TXB = f'Level: higher than\n{len(lv_hi)}/{int(m_.sum())} types, q < 0.05'
    h = [Line2D([], [], marker='o', ls='', color=GR, ms=4.5, label='One Section'), Line2D([], [], color='black', lw=1.8, label='Mean ± 95% CI'),
         Line2D([], [], color='#252525', ls='--', lw=1.4, label='All Cells (Mean of Sections)')]
    TITLE1 = 'CXCL12 Expression by Cell Type'
    assert n_hi_pct + len(ns_pct) == ncomp and n_top >= 0.8 * n_sec     # the in-plot statements must hold in the data
    fig.suptitle(TITLE1, fontsize=18, x=0.01, ha='left', y=0.995)
    fig.legend(handles=h, loc='upper left', bbox_to_anchor=(0.01, 0.975), ncol=3, frameon=False, fontsize=13.5, handletextpad=0.4, columnspacing=1.4)
    fig.subplots_adjust(top=1 - 1.3 / FH, bottom=0.9 / FH, left=0.21, right=0.99)
    # the two statistics boxes, stacked with equal width in panel a, placed as far left as the data allow (data coords)
    from matplotlib.patches import FancyBboxPatch
    rend = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    PADX, PADY, GAP, MARG = 1.0, 0.3, 0.45, 1.5
    tt = [ax.text(0, 0, t_, ha='left', va='top', fontsize=13.5, color=HI, zorder=6, linespacing=1.3) for t_ in (TXA, TXB)]
    ext = [inv.transform(t_.get_window_extent(rend)) for t_ in tt]                 # [[x0, y0], [x1, y1]] in data units (y inverted axis)
    W = max(e[1, 0] - e[0, 0] for e in ext); HS = [abs(e[1, 1] - e[0, 1]) for e in ext]
    SPAN = sum(h_ + 2 * PADY for h_ in HS) + GAP
    # right edge of everything drawn in each row: dots, CI, bar, the All Cells dashed line, and the 'All Cells' label (last row)
    cov = PSu.merge(T1[['cell_type_general', 'ci_high']], left_on='ct', right_on='cell_type_general')
    RE = pd.Series({c: max(cov[cov.ct == c].pct_expressing.max(), T1.set_index('cell_type_general').loc[c, ['ci_high', 'mean_pct_expressing']].max(), all_m) for c in Y})
    lab = [t_ for t_ in ax.texts if t_.get_text() == 'All Cells'][0]
    RE.iloc[-1] = max(RE.iloc[-1], inv.transform(lab.get_window_extent(rend))[1, 0])
    best = None
    for top in np.arange(1.5, n - 0.5 - SPAN + 1e-9, 0.25):                        # possible top edges (row units); never the OGN row
        rows_in = [c for c, i in Y.items() if top - 0.5 < i < top + SPAN + 0.5]
        xl = RE[rows_in].max() + MARG
        if xl + W + 2 * PADX <= ax.get_xlim()[1] and (best is None or xl < best[0] - 1e-9): best = (xl, top, rows_in)
    assert best is not None, 'no free space for the statistics boxes'
    x0, Y_TOP, rows_in = best; y = Y_TOP
    for t_, h_ in zip(tt, HS):
        t_.set_position((x0 + PADX, y + PADY))
        ax.add_patch(FancyBboxPatch((x0, y), W + 2 * PADX, h_ + 2 * PADY, boxstyle='round,pad=0,rounding_size=0.6', mutation_aspect=0.25,
                                    facecolor='white', edgecolor=HI, lw=1, zorder=5))
        y += h_ + 2 * PADY + GAP
    mx = RE[rows_in].max()
    print('box width {:.1f} units'.format(W + 2 * PADX)); print('stat boxes: x from {:.1f}% to {:.1f}% ; rows {} ; max data value in those rows {:.1f}%'.format(x0, x0 + W + 2 * PADX, rows_in, mx), flush=True)
    assert mx < x0 - 0.5, 'statistics boxes overlap data'
    fig.savefig(f'{STEM1}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig); print('wrote', f'{STEM1}.png', flush=True)
    LEVB = ('The same paired test on these levels (BH over {} comparisons) finds them higher in OGN⁺RSPO3⁺ fibroblasts than in {} of the other types '
        '(q < 0.05; not different from {}). ').format(int(m_.sum()), len(lv_hi), ' and '.join(DISP.get(c, c) for c in lv_ns.cell_type_general) if len(lv_ns) else 'none')
    qtx = '< 0.001' if q_hi_max < 0.001 else f'≤ {q_hi_max:.2g}'
    PCTB = (f'all q {qtx}' if len(ns_pct) == 0 else
            f'higher than {n_hi_pct} of the {ncomp} types (q {qtx}; not different from ' + ' and '.join(DISP.get(c, c) for c in ns_pct.cell_type_general) + ')')
    leg1 = ('CXCL12 expression by cell type in CD resections (CosMx 6K; 11 sections, 4 patients; {:,} cells). '
        'Expression is the log-normalised CXCL12 value of each cell (norm_log1p), and a cell is counted as expressing if CXCL12 > 0. '
        'Cell types are the cell_type_general annotation; SELENOP+ fibroblasts are not shown. BEST4⁺ Epi = BEST4⁺ absorptive epithelial cells. '
        'Left: percentage of cells of each type that express CXCL12. Dots: one section each (cell types with at least {} cells in the section); '
        'bars and black lines: mean and 95% CI (t distribution) across sections. The dashed line is the all-cell percentage (all cell types, including SELENOP+ fibroblasts), '
        'averaged over the sections ({:.1f}%; pooled over all cells {:.1f}%). '
        'OGN⁺RSPO3⁺ fibroblasts express CXCL12 in {:.1f}% of cells (95% CI {:.1f}–{:.1f}%; {:.1f}-fold the all-cell percentage) and rank first in {} of {} sections; '
        'the next cell type, {}, expresses it in {:.1f}% of cells. OGN⁺RSPO3⁺ fibroblasts were compared with each other cell type by a paired t-test on the per-section percentages '
        '(sections where both types have at least {} cells), with Benjamini–Hochberg correction over the {} comparisons: {}. '
        'Right: mean log-normalised CXCL12 in the expressing cells of each type, per section, same order. {}'
        'Per-section values, ranks and statistics are in the accompanying CSV.').format(
            len(C), MIN_CELLS, all_m, tot_pct, o.mean_pct_expressing, o.ci_low, o.ci_high, o.mean_pct_expressing / all_m, n_top, n_sec,
            second.label.replace('⁺', '+'), second.mean_pct_expressing, MIN_CELLS, ncomp, PCTB, LEVB)
    wtxt(f'{STEM1}_legend.txt', leg1 + '\n')
    # ======================================================= FIGURE 2: APC abundance per neighbourhood
    if os.environ.get('RUN_FIG2') != '1':   # figure 2 and its table are written only with RUN_FIG2=1
        print('figure 2 skipped (set RUN_FIG2=1 to write it)', flush=True); sys.exit(0)
    CC = C[C.CN.notna()].copy(); CC['CN'] = CC.CN.astype(int)
    assert sorted(CC.CN.unique()) == list(range(1, 11))
    mp = CC.cn_cell_type.map(APC_SHARED)
    apc_ct = CC.ct.where(CC.ct.isin(APC))
    agree = (mp.fillna('-') == apc_ct.fillna('-')).mean(); print('APC call agreement cell_type_general vs cn_per_cell vocabulary: {:.5f}'.format(agree), flush=True)
    assert agree > 0.999
    CC['apc'] = apc_ct; CC['apcg'] = apc_ct.map(APC_GRP)
    cnt = CC.groupby(['section', 'CN']).size().rename('n_cells').to_frame()
    for a in APCG: cnt[f'n_{a}'] = CC[CC.apcg == a].groupby(['section', 'CN']).size().reindex(cnt.index).fillna(0).astype(int)
    for a in APC: cnt[f'n_celltype_{a}'] = CC[CC.apc == a].groupby(['section', 'CN']).size().reindex(cnt.index).fillna(0).astype(int)   # record only
    cnt['n_APC'] = cnt[[f'n_{a}' for a in APCG]].sum(1); assert (cnt.n_APC == cnt[[f'n_celltype_{a}' for a in APC]].sum(1)).all()
    cnt['apc_pct'] = 100 * cnt.n_APC / cnt.n_cells
    for a in APCG: cnt[f'pct_{a}'] = 100 * cnt[f'n_{a}'] / cnt.n_cells
    cnt = cnt.reset_index(); cnt['used'] = cnt.n_cells >= MIN_CN_CELLS
    U = cnt[cnt.used].copy()
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [('OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]   # as the Fig 5B CN map
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {k: short_cn(v) for k, v in NM.items()}
    assert LAB[10] == 'CN10: OGN⁺RSPO3⁺/PLA2G2A⁺ Fib' and LAB[4] == 'CN4: DC / Mac / CD9 Mac'
    # mixed model, reference CN10
    MOD = smf.mixedlm('apc_pct ~ C(CN, Treatment(reference=10))', U, groups=U.section)
    try: mm = MOD.fit(reml=True, method='lbfgs'); FIT = 'lbfgs'
    except np.linalg.LinAlgError: mm = MOD.fit(reml=True, method='powell'); FIT = 'powell'   # lbfgs hits a singular matrix on these data
    mm_nm = MOD.fit(reml=True, method='nm'); assert np.allclose(mm.params, mm_nm.params, atol=0.02), 'mixed model optimisers disagree'
    print('mixed model fitted with', FIT, '| section variance', float(mm.cov_re.iloc[0, 0]), flush=True)
    ci = mm.conf_int()
    st = []
    for c in range(1, 10):
        t = f'C(CN, Treatment(reference=10))[T.{c}]'
        w = U.pivot(index='section', columns='CN', values='apc_pct')[[c, 10]].dropna()
        st.append(dict(CN=c, estimate_vs_CN10=mm.params[t], estimate_ci_low=ci.loc[t, 0], estimate_ci_high=ci.loc[t, 1], p=mm.pvalues[t],
                       n_sections_paired=len(w), n_sections_CN10_lower=int((w[10] < w[c]).sum())))
    ST = pd.DataFrame(st); ST['q'] = multipletests(ST.p, method='fdr_bh')[1]
    rows = []
    for c in range(1, 11):
        v = U[U.CN == c]; m, lo, hi, nn = ci95(v.apc_pct)
        r = dict(CN=c, CN_label=LAB[c], n_sections=nn, n_cells=int(cnt[cnt.CN == c].n_cells.sum()), n_APC=int(cnt[cnt.CN == c].n_APC.sum()),
                 mean_apc_pct=m, ci_low=lo, ci_high=hi, pooled_apc_pct=100 * cnt[cnt.CN == c].n_APC.sum() / cnt[cnt.CN == c].n_cells.sum())
        for a in APCG: r[f'mean_pct_{a}'] = v[f'pct_{a}'].mean()
        for a in APCG: r[f'share_of_APCs_{a}'] = 100 * v[f'n_{a}'].sum() / max(v.n_APC.sum(), 1)
        rows.append(r)
    T2 = pd.DataFrame(rows).merge(ST, on='CN', how='left')
    T2['rank_low_to_high'] = T2.mean_apc_pct.rank(method='min').astype(int)
    print('== FIGURE 2 table'); print(T2.round(3).to_string(), flush=True)
    r10 = T2.set_index('CN').loc[10]; oth = T2[T2.CN != 10]
    n_lower_sig = int(((oth.q < 0.05) & (oth.estimate_vs_CN10 > 0)).sum()); n_higher_sig = int(((oth.q < 0.05) & (oth.estimate_vs_CN10 < 0)).sum())
    SPARSE = r10.rank_low_to_high <= 3 and n_lower_sig >= 6
    TITLE2 = 'APC Abundance by Neighbourhood'   # neutral title (SPARSE is printed only)
    print('CN10 rank (1 = fewest APCs):', r10.rank_low_to_high, '| sparse statement supported:', SPARSE, '| CNs with significantly more APCs than CN10:', n_lower_sig,
          '| fewer:', n_higher_sig, '| title:', TITLE2, flush=True)
    OUT2 = pd.concat([T2.assign(row_type='summary', section='all'),
                      cnt.assign(row_type='section', CN_label=cnt.CN.map(LAB)).rename(columns={'apc_pct': 'mean_apc_pct', 'used': 'used_in_mean'})], ignore_index=True)
    lead = ['row_type', 'CN', 'CN_label', 'section']
    OUT2 = OUT2[lead + [c for c in OUT2 if c not in lead]].apply(sig3)
    OUT2.to_csv(f'{STEM2}.csv', index=False); print('wrote', f'{STEM2}.csv', flush=True)
    # ---- figure 2 (horizontal: rows = CN1..CN10, a = APC % with sections, b = stacked APC-type composition)
    CNO = list(range(1, 11)); Yc = {c: i for i, c in enumerate(CNO)}
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.5, 6.6), sharey=True, gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.05))
    rng = np.random.default_rng(2)
    for c in CNO:
        i = Yc[c]; col, colL = (HI, HIL) if c == 10 else ('#08519c', '#9ecae1')
        v = U[U.CN == c].apc_pct.values; r = T2.set_index('CN').loc[c]
        ax.scatter(v, i + rng.uniform(-0.18, 0.18, len(v)), s=30, color=colL, edgecolor='white', linewidth=0.4, zorder=2)
        ax.plot([r.ci_low, r.ci_high], [i, i], color=col, lw=2.6, zorder=3, solid_capstyle='butt')
        ax.scatter([r.mean_apc_pct], [i], s=90, color=col, zorder=4, edgecolor='white', linewidth=0.8)
    xmax = U.apc_pct.max()
    for c in CNO:
        if c == 10: continue
        r = T2.set_index('CN').loc[c]
        if r.q < 0.05:
            x = max(U[U.CN == c].apc_pct.max(), r.ci_high)
            ax.text(x + 0.02 * xmax, Yc[c], stars(r.q), ha='left', va='center', fontsize=15, color='#252525', fontweight='bold')
    ax.text(max(U[U.CN == 10].apc_pct.max(), r10.ci_high) + 0.02 * xmax, Yc[10], 'Ref.', ha='left', va='center', fontsize=14, color=HI, fontweight='bold')
    ax.set_xlim(0, xmax * 1.12)
    ax.set_yticks(range(len(CNO))); ax.set_yticklabels([LAB[c] for c in CNO], fontsize=14); ax.set_ylim(len(CNO) - 0.45, -0.55)
    for tl, c in zip(ax.get_yticklabels(), CNO):
        if c == 10: tl.set_color(HI); tl.set_fontweight('bold')
    ax.set_xlabel('APCs (% of All Cells in the Neighbourhood)', fontsize=15)
    ax.set_title('APCs as % of Cells', loc='left', fontsize=16, fontweight='bold', pad=5)
    left = np.zeros(len(CNO))
    for a in APCG:
        w = T2.set_index('CN').loc[CNO, f'mean_pct_{a}'].values
        bx.barh(range(len(CNO)), w, left=left, height=0.7, color=APC_COL[a], edgecolor='white', linewidth=0.6, label=DISP.get(a, a))
        left += w
    bx.set_xlim(0, left.max() * 1.05)
    bx.set_xlabel('Mean APC Composition (% of All Cells)', fontsize=15)
    bx.set_title('APC Composition', loc='left', fontsize=16, fontweight='bold', pad=5)
    bx.tick_params(axis='y', length=0)
    bx.legend(loc='center right', bbox_to_anchor=(1.0, 0.42), frameon=False, fontsize=13.5, handlelength=1.1, handletextpad=0.4, borderaxespad=0.2)
    for a_ in (ax, bx): a_.spines[['top', 'right']].set_visible(False); a_.tick_params(axis='x', labelsize=13.5)
    h = [Line2D([], [], marker='o', ls='', color='#9ecae1', ms=6.5, label='One Section'), Line2D([], [], marker='o', color='#08519c', ms=9, lw=2.6, label='Mean ± 95% CI'),
         Line2D([], [], ls='', label='* q < 0.05, *** q < 0.001 vs CN10 (Mixed Model, BH)')]
    fig.suptitle(TITLE2, fontsize=18, fontweight='bold', x=0.01, ha='left', y=1.0)
    fig.legend(handles=h, loc='upper left', bbox_to_anchor=(0.01, 0.955), ncol=3, frameon=False, fontsize=13.5, handletextpad=0.4, columnspacing=1.4)
    fig.subplots_adjust(top=0.84, left=0.285, right=0.99, bottom=0.1)
    fig.savefig(f'{STEM2}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig); print('wrote', f'{STEM2}.png', flush=True)
    s10 = ', '.join(f'{DISP.get(a, a)} {r10[f"mean_pct_{a}"]:.1f}%' for a in APCG)
    leg2 = ('Antigen-presenting cell (APC) abundance in each spatial cellular neighbourhood, normalised to the size of the neighbourhood (all cells, not centred on Tregs). '
        'CD resections (CosMx 6K; 11 sections, 4 patients; {:,} cells with a neighbourhood). Neighbourhoods are the k = 10 cellular neighbourhoods of Fig 4 (CN1–CN10). '
        'APCs are DCs, B cells, CD9⁺ M2 macrophages, Mo-Macs (including inflammatory Mo-Macs) and SELENOP⁺ macrophages; SELENOP⁺ fibroblasts are not APCs. '
        'Left: APCs as a percentage of all cells of the neighbourhood, per section (dots; section × neighbourhood combinations with at least {} cells), '
        'with the mean and 95% CI (t distribution) across sections. Statistics: linear mixed model APC % ~ neighbourhood (reference CN10) + (1 | section); '
        'stars give the Benjamini–Hochberg q over the 9 contrasts with CN10 (* q < 0.05, ** q < 0.01, *** q < 0.001; non-significant contrasts are unmarked). '
        'Right: mean percentage of all cells in the neighbourhood made up by each APC group (DC, B cell, and Macrophage = SELENOP⁺, CD9⁺ M2 and Mo-Mac macrophages; stacked; the total equals the mean on the left). '
        'CN10 (OGN⁺RSPO3⁺/PLA2G2A⁺ fibroblast neighbourhood): APCs are {:.1f}% of cells (95% CI {:.1f}–{:.1f}%), rank {} of 10 from lowest; {}. '
        '{} of the 9 other neighbourhoods have significantly more APCs than CN10 and {} significantly fewer. '
        'Per-section counts and the model estimates are in the accompanying CSV.').format(
            len(CC), MIN_CN_CELLS, r10.mean_apc_pct, r10.ci_low, r10.ci_high, int(r10.rank_low_to_high), s10, n_lower_sig, n_higher_sig)
    wtxt(f'{STEM2}_legend.txt', leg2 + '\n')
    print('done', flush=True)
