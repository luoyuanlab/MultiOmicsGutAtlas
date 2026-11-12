#!/usr/bin/env python
"""CD CosMx: tissue maps of the Treg Suppression Score.

Panels made
  Fig 5A  WTX biopsy maps (inflamed and non-inflamed CD) with enlarged high and low Treg suppression regions
  Fig 5B  resection section 43563_A4 coloured by cellular neighbourhood, with rule-selected high- and low-score windows

Inputs
  /share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/3 Raw data/R6221_<id>/ (WTX cell polygons)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/
    D_shared_inputs/treg_definition_comparison_per_cell.csv (from 10_treg_suppression_score.py)
    B_resection_analyses/cellular_neighborhoods_v2/ (cn_per_cell.parquet, discovery_biopsy/cn_per_cell_biopsy.csv.gz, cn_treg_state_enrichment_all.csv)
    and B_resection_analyses/ligand_receptor_state/treg_local_ligand_receptor.csv.gz (from 11_neighbourhoods_local_tables.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/ (panels, window and Treg CSVs, legend text)

Run order of sections
  F5A_wtx_biopsy_maps -> F5B_cn_section_map

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
# Section F5A_wtx_biopsy_maps: WTX biopsy maps coloured by the Treg Suppression Score, with enlarged high and low suppression regions
# Panels: Fig 5A
# ======================================================================================================================
if _run('F5A_wtx_biopsy_maps'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/42_wtx_biopsy_treg_score_popouts.py'   # working folder of this step (relative paths resolve from here)
    """WTX biopsy maps of the Treg Suppression Score with enlarged high- and low-suppression windows.
Pieces: Inflamed CD 2 (i0476_slide1) and Non-inflamed CD (i0476_slide3). Tregs: core biopsy Tregs coloured by the continuous Treg Suppression Score
(-D6 signature score, z within piece; red = activated, blue = low-suppressive).
Highlighted cell types: DC, B cell, OGN+RSPO3+ fibroblast and PLA2G2A+ fibroblast (PLA2G2A+ Fib share the OGN+RSPO3+ neighbourhood in the
resections); macrophages are not highlighted; epithelium is mid grey (blue and red belong to the Treg score scale).
Windows: with BOXES=placed (default) the boxes in PLACED (x0, y0, w, h in um) are drawn, each placed on a cluster of high- or low-scoring Tregs.
With any other BOXES value the windows come from a rule inside each region of REG: rectangular windows on a 25-um grid with >= 200 cells and
>= 3 Tregs; among windows whose mean Treg Suppression Score z has the intended sign (> 0 high, < 0 low suppression), the one with the largest
directional sum of z (RULE=sum, default; rewards both more Tregs and a stronger state) or with the most Tregs (RULE=count).
Outputs (Fig4A-F_treg_suppression_niches/): wtx_biopsy_treg_score_popouts.png/.pdf, wtx_biopsy_treg_score_popouts_windows.csv (one row per window:
position, size, cell, Treg and key cell-type counts, mean score) and wtx_biopsy_treg_score_popouts_tregs.csv (every core Treg with its score and window)."""
    import numpy as np, pandas as pd, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches'; UM = 0.12028
    RAW = '/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/3 Raw data'; STEP = 25
    import os; RULE = os.environ.get('RULE', 'sum')
    # region (x0, x1, y0, y1 in um, y from the top) and window (w, h) per zoom
    REG = {'i0476_slide1': ('Inflamed CD 2', [('upper left', (300, 1500, 11600, 12600), (500, 400), 'act'), ('lower left', (200, 1300, 12600, 13700), (500, 400), 'low')]),
           'i0476_slide3': ('Non-inflamed CD', [('top', (1250, 2418, 300, 1300), (600, 450), 'act'), ('lower middle', (1250, 2418, 1500, 2700), (500, 450), 'low')])}
    # boxes drawn with BOXES=placed (x0, y0, w, h in um), each placed on a cluster of high- or low-scoring Tregs inside the tissue:
    PLACED = {'i0476_slide1': [('upper left', 'act', 650, 11780, 675, 440),   # x0 650 = position with the most DC + B cells among those scanned (162); 14 Tregs, mean z +0.51
                               ('bottom left', 'low', 270, 13380, 480, 480)],   # square; right edge aligned with the I0453 low box
             'i0476_slide3': [('top', 'act', 1630, 540, 520, 470),   # left edge set to leave out the low-suppressive Treg at the tissue edge; holds all 18 Tregs of the top aggregate
                               ('bottom', 'low', 1740, 2560, 420, 420)]}
    CROP_TOP = {'i0476_slide3': 520}   # overview cropped above y = 520 um (sparse detached cells)
    PID = {'i0476_slide1': 'I0476\n(inflamed CD)', 'i0476_slide3': 'I0453\n(non-inflamed CD)'}   # internal sample identifiers (patient column)
    ZT = {'act': 'High Treg Suppression Region', 'low': 'Low Treg Suppression Region'}
    TT = ['CD4 Trm', 'Gd T', 'CD8 GZMK- Trm', 'CD8 GZMK+ Trm', 'CD4 Effector', 'NK', 'Resting T', 'Exhausted CD4 Tfh', 'Tfh']
    STR = ['Activated Fib', 'Myofibroblast', 'Pericyte', 'ADAMDEC1 Fib', 'Endothelial', 'Crypt Top Fib', 'Neuron']
    GRP = {'DC': '#1a9850', 'B cell': '#7b3294', 'OGN\u207aRSPO3\u207a Fib': '#f46d43', 'PLA2G2A\u207a Fib': '#d4a017',
           'Epithelium': '#a6a6a6', 'Other immune': '#ececec', 'Other stromal': '#f5efe0', 'Other': '#f5f5f5'}
    KEY = ['DC', 'B cell', 'OGN\u207aRSPO3\u207a Fib', 'PLA2G2A\u207a Fib']
    CM = plt.get_cmap('RdBu_r'); NORM = Normalize(-2, 2); ZC = {'act': '#b2182b', 'low': '#2166ac'}
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); PC = PC[(PC.cohort == 'biopsy') & (PC.is_core_treg.astype(str) == 'True')]
    b = pd.read_csv(f'{P4}/B_resection_analyses/cellular_neighborhoods_v2/discovery_biopsy/cn_per_cell_biopsy.csv.gz', index_col=0); b = b[b.section.isin(REG)].copy()
    b['z'] = -PC.signature_score_D6.reindex(b.index); b['z'] = b.groupby('section').z.transform(lambda v: (v - v.mean()) / v.std()); ct = b.ct.astype(str)
    b['grp'] = np.select([ct == 'DC', ct == 'B', ct == 'OGN_RSPO3 Fib', ct == 'PLA2G2A Fib', ct == 'Epithelial', ct.isin(TT + ['Plasma', 'Mac S', 'CD9 Mac', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac',
                          'Inf Mo-Mac', 'Granulocyte', 'Neutrophil', 'Cycling']), ct.isin(STR)], ['DC', 'B cell', 'OGN\u207aRSPO3\u207a Fib', 'PLA2G2A\u207a Fib', 'Epithelium', 'Other immune', 'Other stromal'], 'Other')
    plt.rcParams.update({'font.size': 17})
    fig = plt.figure(figsize=(23, 17)); gs = fig.add_gridspec(2, 3, width_ratios=[0.62, 1.2, 1.2], wspace=0.05, hspace=0.1, bottom=0.075, top=0.97, left=0.03, right=0.98); wins = []
    for r_, (p, (nm, regions)) in enumerate(REG.items()):
        s = b[b.section == p]; T = s[s.z.notna()]; X, Y = s.x.values, s.y.values; isT = s.z.notna().values; zz = s.z.fillna(0).values
        xm, xM, ym, yM = X.min(), X.max(), Y.min(), Y.max(); ZW = []
        for lab, (rx0, rx1, ry0, ry1), (ww, wh), kind in regions:
            best = None
            for x0 in np.arange(rx0, rx1 - ww + 1, STEP):
                mx = (X >= x0) & (X < x0 + ww)
                for y0 in np.arange(ry0, ry1 - wh + 1, STEP):
                    m = mx & (Y >= y0) & (Y < y0 + wh); mt = m & isT
                    if m.sum() < 200 or mt.sum() < 3: continue
                    mz = zz[mt].mean() * (1 if kind == 'act' else -1)
                    if mz <= 0: continue
                    key = (int(mt.sum()), mz) if RULE == 'count' else (mz * mt.sum(), int(mt.sum()))   # 'sum' = directional sum of z (more Tregs AND stronger state)
                    if best is None or key > best[0]: best = (key, x0, y0)
            ZW.append((lab, kind, best[1], best[2], ww, wh))
        ZW = PLACED[p] if os.environ.get('BOXES', 'placed') == 'placed' else ZW
        ax0 = fig.add_subplot(gs[r_, 0])
        bg = s[~s.grp.isin(KEY) & ~isT]
        for g_, sz in [('Other', 1), ('Other immune', 1.5), ('Other stromal', 1.5), ('Epithelium', 2.5)]:
            q = bg[bg.grp == g_]; ax0.scatter(q.x, q.y, s=sz, c=GRP[g_], rasterized=True, lw=0)
        for k in KEY:
            q = s[s.grp == k]; ax0.scatter(q.x, q.y, s=6, c=GRP[k], rasterized=True, lw=0)
        q = T.reindex(T.z.abs().sort_values().index); sc = ax0.scatter(q.x, q.y, s=60, c=q.z, cmap=CM, norm=NORM, edgecolor='k', lw=0.6, zorder=5)
        for i, (lab, kind, x0, y0, ww, wh) in enumerate(ZW):
            ax0.add_patch(Rectangle((x0, y0), ww, wh, fill=False, ec=ZC[kind], lw=3, zorder=6))
        ax0.set_aspect('equal'); ax0.set_ylim(yM + 220, CROP_TOP.get(p, ym) - 20); ax0.axis('off'); ax0.set_anchor('NW'); ax0.set_title(f'Biopsy {PID[p]}', fontsize=18, loc='center')   # overviews left-aligned, title centred over each tissue
        ax0.plot([xm, xm + 500], [yM + 60] * 2, 'k-', lw=3.5); ax0.text(xm + 250, yM + 190, '500 \u00b5m', ha='center', fontsize=16)
        allw = pd.Index([])
        for (_, _, x0, y0, ww, wh) in ZW: allw = allw.union(s.index[(s.x >= x0) & (s.x < x0 + ww) & (s.y >= y0) & (s.y < y0 + wh)])
        pat = p.split('_')[0].upper()   # polygon files are per run (i0476 run holds the I0453 pieces)
        poly = pd.concat([c[c.cell.isin(allw)] for c in pd.read_csv(f'{RAW}/R6221_{pat}/wtx_{pat}-polygons.csv.gz', chunksize=2_000_000)])
        for c_, (lab, kind, x0, y0, ww, wh) in enumerate(ZW):
            x1, y1 = x0 + ww, y0 + wh; w = s[(s.x >= x0) & (s.x < x1) & (s.y >= y0) & (s.y < y1)]; wt = w[w.z.notna()]; pw = poly[poly.cell.isin(w.index)]
            cen = pw.groupby('cell')[['x_global_px', 'y_global_px']].mean() * UM; cm = cen.index.intersection(w.index)
            dd = np.hypot(cen.loc[cm].x_global_px - w.loc[cm].x, cen.loc[cm].y_global_px - w.loc[cm].y); assert dd.median() < 3 and len(cm) > 0.9 * len(w)
            ax = fig.add_subplot(gs[r_, c_ + 1]); verts, faces, edges, lws = [], [], [], []
            for nm_, g in pw.groupby('cell', sort=False):
                tr = nm_ in wt.index; fc = CM(NORM(wt.z[nm_])) if tr else GRP[s.grp.get(nm_, 'Other')]
                verts.append(np.c_[g.x_global_px.values, g.y_global_px.values] * UM); faces.append(fc); edges.append('k' if tr else 'white'); lws.append(1.6 if tr else 0.2)
            ax.add_collection(PolyCollection(verts, facecolors=faces, edgecolors=edges, linewidths=lws))
            ax.scatter(wt.x, wt.y, s=330, facecolors='none', edgecolors='k', linewidths=5, zorder=6)   # black under-ring keeps near-zero (white) scores visible
            ax.scatter(wt.x, wt.y, s=330, facecolors='none', edgecolors=[CM(NORM(v)) for v in wt.z], linewidths=3, zorder=7)
            ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_aspect('equal'); ax.set_anchor('NW'); ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_edgecolor(ZC[kind]); sp.set_linewidth(3.5)
            ax.plot([x0 + 15, x0 + 115], [y1 - 15] * 2, 'k-', lw=3.5, zorder=8); ax.text(x0 + 65, y1 - 34, '100 \u00b5m', ha='center', va='bottom', fontsize=16, zorder=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))
            hi, lo = int((wt.z > 0).sum()), int((wt.z <= 0).sum())
            ax.set_title(ZT[kind], fontsize=19, color=ZC[kind])
            wins.append(dict(piece=p, piece_label=nm, region=lab, kind=kind, x0_um=x0, y0_um=y0, width_um=ww, height_um=wh, n_cells=len(w), n_tregs=len(wt), mean_treg_suppression_score_z=wt.z.mean(),
                             n_above_mean=hi, n_below_mean=lo, n_DC=int((w.ct == 'DC').sum()), n_B=int((w.ct == 'B').sum()), n_OGN_fib=int((w.ct == 'OGN_RSPO3 Fib').sum()),
                             n_PLA2G2A_fib=int((w.ct == 'PLA2G2A Fib').sum()), n_epithelial=int((w.ct == 'Epithelial').sum())))
    H = [Line2D([0], [0], marker='s', ls='', mfc=GRP[k], mec=GRP[k], ms=15, label=k) for k in GRP if k != 'Other']
    fig.legend(handles=H, loc='lower left', ncol=4, frameon=False, fontsize=17, bbox_to_anchor=(0.03, -0.005), columnspacing=1.2, handletextpad=0.3, borderaxespad=0, borderpad=0)
    cax = fig.add_axes([0.62, 0.045, 0.28, 0.014]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal')
    cb.set_label('Treg Suppression Score (z within biopsy)\n\u2190 low suppression          high suppression \u2192', fontsize=17); cb.ax.tick_params(labelsize=15)
    fig.savefig(f'{OUT}/wtx_biopsy_treg_score_popouts.png', dpi=300, bbox_inches='tight', pad_inches=0.03); fig.savefig(f'{OUT}/wtx_biopsy_treg_score_popouts.pdf', bbox_inches='tight', pad_inches=0.03)
    TT_ = b[b.z.notna()][['section', 'x', 'y', 'z']].rename(columns={'section': 'piece', 'x': 'x_um', 'y': 'y_um', 'z': 'treg_suppression_score_z'})
    TT_['biopsy'] = TT_.piece.map({k: v.split('\n')[0] for k, v in PID.items()}); TT_['region'] = ''
    for w_ in wins:
        m_ = (TT_.piece == w_['piece']) & TT_.x_um.between(w_['x0_um'], w_['x0_um'] + w_['width_um'], inclusive='left') & TT_.y_um.between(w_['y0_um'], w_['y0_um'] + w_['height_um'], inclusive='left')
        TT_.loc[m_, 'region'] = ZT[w_['kind']]
    TT_.to_csv(f'{OUT}/wtx_biopsy_treg_score_popouts_tregs.csv', float_format='%.4g')
    pd.DataFrame(wins).to_csv(f'{OUT}/wtx_biopsy_treg_score_popouts_windows.csv', index=False, float_format='%.3g'); print(pd.DataFrame(wins).round(2).to_string(index=False))


# ======================================================================================================================
# Section F5B_cn_section_map: resection section 43563_A4 coloured by cellular neighbourhood, with rule-selected high and low Treg Suppression Score windows
# Panels: Fig 5B
# ======================================================================================================================
if _run('F5B_cn_section_map'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/73_cn8_tls_zoom.py'   # working folder of this step (relative paths resolve from here)
    """Cellular-neighbourhood map of CD resection section 43563_A4 with rule-selected windows of high and low Treg Suppression Score.
Inputs: cellular_neighborhoods_v2/cn_per_cell.parquet (CN_k10, 1-based; x_um / y_um = all_samples px x 0.12028) and cn_treg_state_enrichment_all.csv
(CN names), from 11_neighbourhoods_local_tables.py; ligand_receptor_state/treg_local_ligand_receptor.csv.gz (act_z = Treg Suppression Score,
-signature_score_D6 z-scored within section over all Tregs).
Window = the extent of field of view 29 of the section (W x H um). Windows are slid over the whole section on a 25-um grid; windows with >= 15 Tregs are
ranked by the mean act_z of their Tregs (t-based 95% CI over Tregs; spatial correlation is not modelled and the windows are selected extremes).
Greedy selection of the top 3 (H1-H3) and bottom 3 (L1-L3) windows, each new window at least one window width (x) or height (y) from every chosen
window -> fig4_cn_section_map_multibox_boxes.csv, with the CN composition of each window and the fraction of its cells in the central CN8 aggregate
(largest connected CN8 component on a 50-um grid; definitions 'dil' = occupied bins with 1-bin dilation, 'any' = occupied bins, 'core' = bins with
>= 3 cells and >= 50% CN8 cells).
Boxes shown: High 1 = H1, Low 1 = L2, Low 2 = L3, and High 2 = the first further high window of the same greedy ranking (continued beyond 3; H3 skipped)
with >= 15 Tregs, a 95% CI excluding 0 and a centre at least two window widths from the other three boxes (H4 on this section, asserted).
Display: tissue pieces are moved closer together (translation only; distances within a piece unchanged); a window is drawn at its translated
centre with true size; the piece of a window = the majority piece of its cells.
Outputs (Fig4A-F_treg_suppression_niches/): fig4_cn_section_map_boxes.png/.pdf + _legend.txt, fig4_cn_section_map_boxes.csv (the four boxes shown),
fig4_cn_section_map_boxes_high2_search.csv (windows examined for High 2), fig4_cn_section_map_multibox_boxes.csv (top-3 / bottom-3 windows)."""
    import os
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX', 'run with OPENBLAS_CORETYPE=SkylakeX (wrong matmul on gpu11 otherwise)'
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    from scipy import ndimage, stats
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; BR = f'{P4}/B_resection_analyses'
    SEC, HIGH_FOV, BIN, STEP = '43563_A4', '29', 50, 25
    c = pq.read_table(f'{BR}/cellular_neighborhoods_v2/cn_per_cell.parquet', columns=['section', 'fov', 'x_um', 'y_um', 'cell_type', 'CN_k10', '__index_level_0__']).to_pandas()
    T = pd.read_csv(f'{BR}/ligand_receptor_state/treg_local_ligand_receptor.csv.gz', index_col=0)
    assert T.index.isin(c.index).all(); c['act_z'] = T.act_z.reindex(c.index); c['is_treg'] = c.act_z.notna()
    assert (c.loc[c.is_treg, 'cell_type'] == 'Treg').all()
    print('all CN8 Tregs (11 sections): mean act_z', round(c[c.is_treg & (c.CN_k10 == 8)].act_z.mean(), 3), 'n', int((c.is_treg & (c.CN_k10 == 8)).sum()))

    def components(s, mode):
        gx = np.floor((s.x_um - s.x_um.min()) / BIN).astype(int).values; gy = np.floor((s.y_um - s.y_um.min()) / BIN).astype(int).values
        m = (s.CN_k10 == 8).values; tot = np.zeros((gy.max() + 1, gx.max() + 1)); n8 = tot.copy()
        np.add.at(tot, (gy, gx), 1); np.add.at(n8, (gy[m], gx[m]), 1)
        if mode == 'dil': lab = ndimage.label(ndimage.binary_dilation(n8 > 0, iterations=1), structure=np.ones((3, 3)))[0]
        elif mode == 'any': lab = ndimage.label(n8 > 0, structure=np.ones((3, 3)))[0]
        else: lab = ndimage.label((n8 / np.maximum(tot, 1) >= 0.5) & (tot >= 3), structure=np.ones((3, 3)))[0]
        comp = np.where(m, lab[gy, gx], 0)   # component id of CN8 cells, 0 otherwise
        return comp, lab, gx, gy

    def mci(v):
        v = np.asarray(v, float); n = len(v); mu = v.mean(); h = stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n)
        return mu, mu - h, mu + h, stats.ttest_1samp(v, 0).pvalue

    s0 = c[c.section == SEC]
    fh = s0[s0.fov.astype(str) == HIGH_FOV]; W, H = np.ptp(fh.x_um.values), np.ptp(fh.y_um.values)
    # ---- window search over the whole section; window = W x H (the extent of FOV 29)
    # Windows slid on a 25-um grid with >= 15 Tregs, ranked by mean act_z; greedy top-3 high and bottom-3 low, each new window at least one
    # window width along x or one window height along y from every chosen window -> fig4_cn_section_map_multibox_boxes.csv
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from scipy.spatial import cKDTree
    OUT = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches'; CND = f'{BR}/cellular_neighborhoods_v2'
    MIN_T_MB, NBOX = 15, 3
    s = s0.copy(); X, Y, TR, Z, CN = s.x_um.values, s.y_um.values, s.is_treg.values, s.act_z.fillna(0).values, s.CN_k10.values
    tx, ty, tz = X[TR], Y[TR], Z[TR]; o = np.argsort(ty); tx, ty, tz = tx[o], ty[o], tz[o]
    xs = np.arange(X.min() - W / 2, X.max() - W / 2 + STEP, STEP); ys = np.arange(Y.min() - H / 2, Y.max() - H / 2 + STEP, STEP)
    win = []
    for x0 in xs:
        k = (tx >= x0) & (tx <= x0 + W)
        if k.sum() < MIN_T_MB: continue
        yy, zz = ty[k], tz[k]; cz = np.r_[0, np.cumsum(zz)]; cz2 = np.r_[0, np.cumsum(zz ** 2)]
        a = np.searchsorted(yy, ys, 'left'); b = np.searchsorted(yy, ys + H, 'right'); n = b - a; g = n >= MIN_T_MB
        for y0, a_, b_, n_ in zip(ys[g], a[g], b[g], n[g]): win.append((x0, y0, n_, (cz[b_] - cz[a_]) / n_))
    win = pd.DataFrame(win, columns=['x0', 'y0', 'n', 'mean']); print('windows with >=', MIN_T_MB, 'Tregs:', len(win))
    def pick(df, chosen):
        out = []
        for _, r in df.iterrows():
            if all(abs(r.x0 - c[0]) >= 2 * W or abs(r.y0 - c[1]) >= 2 * H for c in chosen + out): out.append((r.x0, r.y0))
            if len(out) == NBOX: break
        return out
    hi_ = pick(win.sort_values('mean', ascending=False), []); lo_ = pick(win.sort_values('mean'), hi_)
    # fraction of each window's cells in the central (largest) CN8 aggregate under the three definitions of components()
    AGG = {}
    for mode in ['dil', 'any', 'core']:
        comp, *_ = components(s0, mode); big = pd.Series(comp[comp > 0]).value_counts().index[0]; AGG[mode] = comp == big
    recs = []
    for lab_, boxes in [('H', hi_), ('L', lo_)]:
        for i, (x0, y0) in enumerate(boxes, 1):
            w = (X >= x0) & (X <= x0 + W) & (Y >= y0) & (Y <= y0 + H); zt = s.act_z.values[w & TR]; mu, lo, hi, p = mci(zt)
            comp_ = pd.Series(CN[w]).value_counts(normalize=True).sort_index()
            recs.append(dict(box=f'{lab_}{i}', region='High Treg Suppression' if lab_ == 'H' else 'Low Treg Suppression', centre_x_um=x0 + W / 2, centre_y_um=y0 + H / 2,
                             width_um=W, height_um=H, n_cells=int(w.sum()), n_tregs=len(zt), mean_act_z=mu, ci95_lo=lo, ci95_hi=hi, p_vs0=p,
                             dominant_CN=int(comp_.idxmax()), dominant_CN_frac=comp_.max(), **{f'frac_CN{k}': comp_.get(k, 0.0) for k in range(1, 11)},
                             **{f'frac_cells_in_central_CN8_aggregate_{m}': AGG[m][w].mean() for m in AGG}))
    MB = pd.DataFrame(recs)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [('OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    NMS = {k: short_cn(v) for k, v in NM.items()}
    MB['dominant_CN_name'] = MB.dominant_CN.map(NMS)
    MB.to_csv(f'{OUT}/fig4_cn_section_map_multibox_boxes.csv', index=False, float_format='%.3g'); print(MB.round(3).to_string(index=False))
    H3, L3 = MB[MB.box.str.startswith('H')], MB[MB.box.str.startswith('L')]
    print(f'High boxes: dominant CN8 {int((H3.dominant_CN == 8).sum())}/3, CN4 {int((H3.dominant_CN == 4).sum())}/3 | Low boxes: CN10 {int((L3.dominant_CN == 10).sum())}/3, CN5 {int((L3.dominant_CN == 5).sum())}/3')
    print('High boxes with >=50% cells in the central CN8 aggregate:', {m: H3.box[H3[f'frac_cells_in_central_CN8_aggregate_{m}'] >= 0.5].tolist() for m in AGG})

    # ---- display transform: tissue pieces moved closer together (translation only), as in the CN section maps of 12_neighbourhood_treg_score.py
    GAP = 60
    def squeeze(v):
        b = np.floor((v - v.min()) / BIN).astype(int); occ = np.zeros(b.max() + 2, bool); occ[np.unique(b)] = True
        shift = np.zeros(len(occ)); run = 0; cut = 0
        for i in range(len(occ)):
            if not occ[i]:
                run += 1
                if run * BIN > GAP: cut += BIN
            else: run = 0
            shift[i] = cut
        return lambda x: x - shift[np.clip(np.floor((np.asarray(x) - v.min()) / BIN).astype(int), 0, len(shift) - 1)]
    m_ = s0[['x_um', 'y_um', 'CN_k10']].copy(); fx, fy = squeeze(m_.x_um.values), squeeze(m_.y_um.values)
    m_ = m_.assign(x_um=fx(m_.x_um.values), y_um=fy(m_.y_um.values))
    gx = np.floor((m_.x_um - m_.x_um.min()) / BIN).astype(int).values; gy = np.floor((m_.y_um - m_.y_um.min()) / BIN).astype(int).values
    occ = np.zeros((gy.max() + 1, gx.max() + 1), bool); occ[gy, gx] = True
    labp, npc = ndimage.label(ndimage.binary_dilation(occ, iterations=1), structure=np.ones((3, 3)))
    m_ = m_.assign(piece=labp[gy, gx]); order = list(m_.piece.value_counts().index)
    SHIFT = {order[0]: np.zeros(2)}; placed = m_[m_.piece == order[0]][['x_um', 'y_um']].values
    for pc in order[1:]:
        Pp = m_[m_.piece == pc][['x_um', 'y_um']].values; tree = cKDTree(placed[::5]); c0 = placed.mean(0); d = c0 - Pp.mean(0); L = np.linalg.norm(d); u = d / L
        best = np.zeros(2); t = 0.0
        while t < L:
            t2 = t + 25; sh = u * t2
            if tree.query(Pp[::5] + sh, k=1)[0].min() < GAP: break
            best, t = sh, t2
        SHIFT[pc] = best; placed = np.vstack([placed, Pp + best])
    narrow = min(order[1:], key=lambda p: np.ptp(m_[m_.piece == p].x_um.values))
    others = np.vstack([m_[m_.piece == p][['x_um', 'y_um']].values + SHIFT[p] for p in order if p != narrow])[::3]; tree = cKDTree(others)
    Pn = m_[m_.piece == narrow][['x_um', 'y_um']].values[::3]; big = order[0]; Bp = m_[m_.piece == big][['x_um', 'y_um']].values + SHIFT[big]
    base = np.array([Bp[:, 0].min() - GAP - Pn[:, 0].max(), Bp[:, 1].max() - Pn[:, 1].max()]); best = None
    for dyu in np.arange(0, 4001, 50):
        for dxr in np.arange(0, 4001, 50):
            sh_ = base + np.array([dxr, -dyu])
            if tree.query(Pn + sh_, k=1)[0].min() < GAP: continue
            sc = dyu + dxr
            if best is None or sc > best[0]: best = (sc, sh_)
    SHIFT[narrow] = best[1]
    shv = np.vstack([SHIFT[p] for p in m_.piece.values]); m_ = m_.assign(dx=m_.x_um.values + shv[:, 0], dy=m_.y_um.values + shv[:, 1])
    print('pieces', npc, '| shifts um', {int(k): np.round(v).tolist() for k, v in SHIFT.items()})
    PIECE = m_.piece.values
    def to_disp(cx, cy, w):   # raw-frame window centre -> display centre (majority piece of the window's cells)
        pc = pd.Series(PIECE[w]).mode().iloc[0]; return float(fx(cx)) + SHIFT[pc][0], float(fy(cy)) + SHIFT[pc][1]
    PAL = dict(zip(range(1, 11), sns.color_palette('tab10', 10).as_hex())); HC, LC = '#b2182b', '#2166ac'
    plt.rcParams.update({'font.size': 17})
    def cn_map(ax, alpha=1.0, grey=False):
        for k in range(1, 11):
            q = m_[m_.CN_k10 == k]; ax.scatter(q.dx, q.dy, s=0.25, c='#cfcfcf' if grey else PAL[k], rasterized=True, lw=0, alpha=alpha)
        ax.set_aspect('equal'); ax.invert_yaxis(); ax.axis('off')
        xm, ym = m_.dx.min(), m_.dy.min(); ax.plot([xm, xm + 1000], [ym + 150] * 2, 'k-', lw=4, solid_capstyle='butt'); ax.text(xm + 500, ym + 60, '1 mm', ha='center', va='bottom', fontsize=15)
    CNH = [Line2D([0], [0], marker='o', ls='', mfc=PAL[k], mec=PAL[k], ms=13, label=NMS[k]) for k in range(1, 11)]

    # ---- boxes shown: High 1 = H1, Low 1 = L2, Low 2 = L3 (display labels by rank). High 2 = the first further-ranked high window of the same
    # greedy rule (continued beyond 3; H3 skipped, its CI includes 0) with >= 15 Tregs, 95% CI excluding 0 and centre >= 2 window widths from
    # High 1 (H1), Low 1 (L2) and Low 2 (L3); its distance to FOV 29 is recorded in the High 2 search table.
    def pick_n(df, chosen, n):
        out = []
        for _, r in df.iterrows():
            if all(abs(r.x0 - c[0]) >= 2 * W or abs(r.y0 - c[1]) >= 2 * H for c in chosen + out): out.append((r.x0, r.y0))
            if len(out) == n: break
        return out
    hi_ext = pick_n(win.sort_values('mean', ascending=False), [], 40); assert [tuple(map(float, b)) for b in hi_ext[:3]] == [tuple(map(float, b)) for b in hi_]
    def box_rec(lbl, x0, y0):
        w = (X >= x0) & (X <= x0 + W) & (Y >= y0) & (Y <= y0 + H); zt = s.act_z.values[w & TR]; mu, lo, hi, p = mci(zt)
        comp_ = pd.Series(CN[w]).value_counts(normalize=True).sort_index()
        return dict(box=lbl, region='High Treg Suppression' if lbl[0] == 'H' else 'Low Treg Suppression', centre_x_um=x0 + W / 2, centre_y_um=y0 + H / 2,
                    width_um=W, height_um=H, n_cells=int(w.sum()), n_tregs=len(zt), mean_act_z=mu, ci95_lo=lo, ci95_hi=hi, p_vs0=p,
                    dominant_CN=int(comp_.idxmax()), dominant_CN_name=NMS[int(comp_.idxmax())], dominant_CN_frac=comp_.max(), **{f'frac_CN{k}': comp_.get(k, 0.0) for k in range(1, 11)})
    SEL3 = MB[MB.box.isin(['H1', 'L2', 'L3'])]; fc = SEL3[['centre_x_um', 'centre_y_um']].values
    fov29 = s0[s0.fov.astype(str) == HIGH_FOV]; f29c = np.array([fov29.x_um.min() + W / 2, fov29.y_um.min() + H / 2])
    cand = []; chosen = None
    for k, (x0, y0) in enumerate(hi_ext, 1):
        if k in (1, 3): continue
        r = box_rec(f'H{k}', x0, y0); d = np.hypot(fc[:, 0] - r['centre_x_um'], fc[:, 1] - r['centre_y_um'])
        r.update({f'dist_to_{b}_um': v for b, v in zip(SEL3.box, d)}); r['dist_to_FOV29_um'] = float(np.hypot(*(f29c - [r['centre_x_um'], r['centre_y_um']])))
        r['ci_excludes_0'] = bool(r['ci95_lo'] > 0); r['far_enough'] = bool(d.min() >= 2 * W); cand.append(r)
        if r['ci_excludes_0'] and r['far_enough'] and r['n_tregs'] >= MIN_T_MB: chosen = r; break
    print(pd.DataFrame(cand)[['box', 'n_tregs', 'mean_act_z', 'ci95_lo', 'ci95_hi', 'dominant_CN', 'dist_to_H1_um', 'dist_to_L2_um', 'dist_to_L3_um', 'dist_to_FOV29_um', 'ci_excludes_0', 'far_enough']].round(3).to_string(index=False))
    assert chosen is not None, 'no High window satisfies the rule'
    pd.DataFrame(cand).to_csv(f'{OUT}/fig4_cn_section_map_boxes_high2_search.csv', index=False, float_format='%.3g')
    H2L = chosen['box']; KEEP = {'H1': 'High 1', H2L: 'High 2', 'L2': 'Low 1', 'L3': 'Low 2'}
    FB = pd.concat([MB[MB.box.isin(['H1', 'L2', 'L3'])], pd.DataFrame([{k: v for k, v in chosen.items() if k in MB.columns}])]).set_index('box').loc[list(KEEP)].reset_index()
    FB.insert(1, 'display_label', FB.box.map(KEEP)); FB = FB.rename(columns={'box': 'rank_label'})
    cxy = FB[['centre_x_um', 'centre_y_um']].values
    for j, lb in enumerate(FB.display_label): FB[f'dist_to_{lb.replace(" ", "")}_um'] = np.hypot(cxy[:, 0] - cxy[j, 0], cxy[:, 1] - cxy[j, 1])
    print('chosen High 2:', H2L)
    assert H2L == 'H4' and KEEP == {'H1': 'High 1', 'H4': 'High 2', 'L2': 'Low 1', 'L3': 'Low 2'}, KEEP   # the rule selects H4 as High 2 on this section
    FB[['rank_label', 'display_label', 'region', 'centre_x_um', 'centre_y_um', 'width_um', 'height_um', 'n_cells', 'n_tregs', 'mean_act_z', 'ci95_lo', 'ci95_hi', 'p_vs0',
        'dominant_CN', 'dominant_CN_name', 'dominant_CN_frac'] + [f'frac_CN{k}' for k in range(1, 11)] + [c for c in FB.columns if c.startswith('dist_to_')]].rename(columns={'mean_act_z': 'mean_treg_suppression_score_z'}).to_csv(f'{OUT}/fig4_cn_section_map_boxes.csv', index=False, float_format='%.3g')
    LABEL_BELOW = {'Low 1'}
    fig, ax = plt.subplots(figsize=(11, 11)); cn_map(ax)
    for _, r in FB.iterrows():
        w = (X >= r.centre_x_um - W / 2) & (X <= r.centre_x_um + W / 2) & (Y >= r.centre_y_um - H / 2) & (Y <= r.centre_y_um + H / 2)
        cx, cy = to_disp(r.centre_x_um, r.centre_y_um, w); hi_b = r.rank_label.startswith('H'); col = HC if hi_b else LC
        ax.add_patch(Rectangle((cx - W / 2, cy - H / 2), W, H, fill=False, ec=col, lw=3.5, ls='-' if hi_b else (0, (3, 1.5)), zorder=5))
        # box label + mean Treg Suppression Score (z) of the box's Tregs (true minus, 2 decimals) + n Tregs (small), same colour as the box.
        # Labels sit above the box, except Low 1 (below: Low 2 sits directly above it).
        val = f'{r.mean_act_z:+.2f}'.replace('-', '\u2212'); bb = dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1)
        if r.display_label in LABEL_BELOW:
            ax.annotate(f'{r.display_label}\n{val}', (cx, cy + H / 2), xytext=(0, -8), textcoords='offset points', ha='center', va='top', fontsize=17, fontweight='bold', color=col, zorder=6, bbox=bb, linespacing=1.05)
            ax.annotate(f'n = {r.n_tregs} Tregs', (cx, cy + H / 2), xytext=(0, -45), textcoords='offset points', ha='center', va='top', fontsize=12, color=col, zorder=6, bbox=bb)
        else:
            ax.annotate(f'n = {r.n_tregs} Tregs', (cx, cy - H / 2), xytext=(0, 4), textcoords='offset points', ha='center', va='bottom', fontsize=12, color=col, zorder=6, bbox=bb)
            ax.annotate(f'{r.display_label}\n{val}', (cx, cy - H / 2), xytext=(0, 20), textcoords='offset points', ha='center', va='bottom', fontsize=17, fontweight='bold', color=col, zorder=6, bbox=bb, linespacing=1.05)
    ax.set_title(f'CD Resection {SEC}', fontsize=17)
    Hh = CNH + [Line2D([0], [0], marker='s', ls='', mfc='none', mec=HC, mew=3, ms=14, label='High Treg Suppression Region'),
                Line2D([0], [0], color=LC, lw=3, ls=(0, (3, 1.5)), label='Low Treg Suppression Region'),
                Line2D([0], [0], ls='', marker='', label='Box value: mean Treg Suppression\nScore (z) of the Tregs in the box')]
    fig.legend(handles=Hh, loc='center left', bbox_to_anchor=(0.88, 0.5), ncol=1, frameon=False, fontsize=15, handletextpad=0.5, labelspacing=0.6)
    fig.savefig(f'{OUT}/fig4_cn_section_map_boxes.png', dpi=300, bbox_inches='tight', pad_inches=0.03); fig.savefig(f'{OUT}/fig4_cn_section_map_boxes.pdf', bbox_inches='tight', pad_inches=0.03); plt.close(fig)
    mm = lambda v, fmt='.2f': format(v, fmt).replace('-', '\u2212')   # true minus in numbers only
    desc = '; '.join(f'{r.display_label}: n = {r.n_tregs} Tregs, mean {mm(r.mean_act_z, "+.2f")} (95% CI {mm(r.ci95_lo)} to {mm(r.ci95_hi)}), dominant {r.dominant_CN_name}' for _, r in FB.iterrows())
    open(f'{OUT}/fig4_cn_section_map_boxes_legend.txt', 'w').write(
        f'CD resection section {SEC} (CosMx 6K), every cell coloured by its spatial cellular neighbourhood (k = 10). Tissue pieces are moved closer together '
    'for display (translation only; distances within each piece unchanged; scale bar, 1 mm). Boxes are {:.0f} x {:.0f} um windows (the size of one field of view). '
    'Box selection: a window of this size was slid over the whole section on a 25 um grid; windows with at least 15 Tregs ({} windows) were ranked by the mean '
    'Treg Suppression Score of the Tregs inside (Treg Suppression Score = the scRNA-seq-derived late-versus-early Treg pseudotime signature score, sign reversed, z-scored over all Tregs within the section; section mean = 0; higher = more suppressive, activated state); '
    'the top 3 (high) and bottom 3 (low) windows were picked greedily, each at least one window width from every other picked window. The four boxes shown are a '
    'display subset of rule-ranked windows (solid red, high Treg suppression; dashed blue, low Treg suppression; numbered by rank): High 1 = the top-ranked '
    'window; Low 1 and Low 2 = the 2nd and 3rd lowest-ranked windows; High 2 = the first further-ranked high window (greedy selection continued beyond the top 3 under '
    'the same separation rule) whose 95% CI excludes 0 and whose centre lies at least two window widths from the other three boxes ({}; the 2nd and 3rd ranked high '
    'windows were not shown because their CIs included 0). The value printed under each box label is the mean Treg Suppression Score (z) of the Tregs inside that box, and the small text gives the number of Tregs in the box. {}. 95% CIs are t-based '
    'over Tregs and do not account for spatial correlation; windows are selected extremes, so their means overstate typical regional differences. '
    'The top-3 / bottom-3 windows are in fig4_cn_section_map_multibox_boxes.csv, the High 2 windows examined in fig4_cn_section_map_boxes_high2_search.csv, '
    'and the four shown in fig4_cn_section_map_boxes.csv.\n'.format(W, H, len(win), 'ranked ' + H2L, desc))
    print(FB[['rank_label', 'display_label', 'n_tregs', 'mean_act_z', 'ci95_lo', 'ci95_hi', 'dominant_CN_name', 'dominant_CN_frac']].round(3).to_string(index=False))
