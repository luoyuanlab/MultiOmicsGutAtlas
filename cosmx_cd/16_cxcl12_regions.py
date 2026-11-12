#!/usr/bin/env python
"""CD CosMx 6K resections: CXCL12-High Regions, antigen-presenting cells (APCs) and the Treg Suppression Score.

Panels: Fig 5E (section 43563_A4 maps), Fig 5F (one tissue piece of 43563_A8), Supp 9E (CXCL12 by cell type),
        Supp 9F (score by local APC quintile in CXCL12-Low vs CXCL12-High Regions).
Inputs: resection combined/ tables and per-section norm_log1p.h5ad; ligand_receptor_state/ and cellular_neighborhoods_v2/ (11, 12).
Outputs: A_fig4_panels/fig4_cxcl12_region_tregs_exp.csv; Fig 5E-F and Supp 9E-F PNGs (Supp 9E also its CSV);
         archive/supp_apc_abundance_by_neighbourhood.csv (APC share per CN, read by 17).
Run: python 16_cxcl12_regions.py, or RUN_SECTIONS=<id>[,<id>...] for single sections (file order).
"""
import os, gc, warnings
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from anndata.io import read_elem
from scipy import sparse, ndimage, stats
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
warnings.filterwarnings('ignore')

_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


RJ = '/share/fsmresfiles/UC/revision_jenny'
P4 = f'{RJ}/cd_treg_rebuild_rowfix/panels/fig4_and_related'; OUTA = f'{P4}/A_fig4_panels'
LR = f'{P4}/B_resection_analyses/ligand_receptor_state'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028
FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
          '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
          '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
          '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
          '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
          '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
APC = ['DC', 'SELENOP+ Mac', 'CD9+ M2 Mac', 'Mo-Mac', 'B cell']
GAP, SBIN = 60, 50                     # display layout: gaps > 60 um squeezed on a 50 um grid; tissue pieces moved together (translations only)


def load_cells(sections):
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa = sa[sa.section.isin(sections)]; sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
    return sa


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


def stitch(xy, sec, v8):
    """Tissue pieces (50 um grid components) moved towards the largest piece; v8: overlapping pieces moved to the nearest clear position."""
    xy = pd.DataFrame(xy, columns=['x', 'y'])
    gx = np.floor((xy.x - xy.x.min()) / SBIN).astype(int).values; gy = np.floor((xy.y - xy.y.min()) / SBIN).astype(int).values
    occ = np.zeros((gy.max() + 1, gx.max() + 1), bool); occ[gy, gx] = True
    lab, npc = ndimage.label(ndimage.binary_dilation(occ, iterations=1), structure=np.ones((3, 3)))
    piece = lab[gy, gx]; order = list(pd.Series(piece).value_counts().index)
    SHIFT = {order[0]: np.zeros(2)}; placed = xy.values[piece == order[0]]
    for pc in order[1:]:
        P = xy.values[piece == pc]; tree = cKDTree(placed[::5]); c0 = placed.mean(0)
        d = c0 - P.mean(0); Lh = np.linalg.norm(d); u = d / Lh; best = np.zeros(2); t = 0.0
        if v8 and tree.query(P[::5], k=1)[0].min() < GAP:
            cand_ = []
            for ang in np.arange(0, 2 * np.pi, np.pi / 8):
                v_ = np.array([np.cos(ang), np.sin(ang)])
                for r_ in np.arange(25, 20001, 25):
                    if tree.query(P[::5] + v_ * r_, k=1)[0].min() >= GAP: cand_.append((r_, v_ * r_)); break
            best = min(cand_, key=lambda z: z[0])[1]; Lh = 0.0
        while t < Lh:
            t2 = t + 25; sh = u * t2
            if tree.query(P[::5] + sh, k=1)[0].min() < GAP: break
            best, t = sh, t2
        SHIFT[pc] = best; placed = np.vstack([placed, P + best])
    if sec == '43563_A4':                                                   # narrow piece placed left of the largest one, as the CN map of 14
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


# ======================================================================================================================
# Section F5E_cxcl12_region_tregs: per-Treg table of section 54033_B4 (20-Treg mean score sets the colour limits of Fig 5E-F)
# ======================================================================================================================
if _run('F5E_cxcl12_region_tregs'):
    SEC, KNN = '54033_B4', 20
    sa = load_cells([SEC])
    with h5py.File(f'{R}/{FOLDER[SEC]}/Processed/norm_log1p.h5ad', 'r') as h: ob = read_elem(h['obs'])
    nm = pd.Index(SEC + '_' + ob.index.astype(str)); cells = sa.loc[nm[nm.isin(sa.index)]].copy()
    fx, fy = squeeze(cells.x.values), squeeze(cells.y.values)
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0); T = T[T.section == SEC].copy()
    T = T.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left')
    T = T.join(sa[['x', 'y']], how='inner', rsuffix='_s')
    cells['x'], cells['y'] = fx(cells.x.values), fy(cells.y.values); T['x'], T['y'] = fx(T.x.values), fy(T.y.values)
    piece, SHIFT, npc = stitch(cells[['x', 'y']].values, SEC, v8=False); cells['piece'] = piece
    sh = np.vstack([SHIFT[p] for p in cells.piece.values]); cells['x'] += sh[:, 0]; cells['y'] += sh[:, 1]
    tp = cells.piece.values[cKDTree(cells[['x', 'y']].values - sh).query(T[['x', 'y']].values)[1]]
    shT = np.vstack([SHIFT[p] for p in tp]); T['x'] += shT[:, 0]; T['y'] += shT[:, 1]
    _, nb = cKDTree(T[['x', 'y']].values).query(T[['x', 'y']].values, k=min(KNN, len(T)))
    for src, dst in [('rec_CXCR4', 'CXCR4 smoothed'), ('act_z', 'suppression score smoothed')]:
        vv = T[src].values; T[dst] = np.nanmean(np.where(np.isnan(vv[nb]), np.nan, vv[nb]), axis=1)
    T[['x', 'y', 'rec_CXCR4', 'act_z', 'lig_CXCL12', 'lig_CXCL12_exp', 'n_exp_producers_50um', 'apc_frac', 'CXCR4 smoothed', 'suppression score smoothed']].rename(
        columns={'act_z': 'treg_suppression_score_z'}).to_csv(f'{OUTA}/fig4_cxcl12_region_tregs_exp.csv', float_format='%.4g')
    print('tissue pieces', npc, '| Tregs', len(T), flush=True)


# ======================================================================================================================
# Section F5EF_section_maps: Fig 5E (four maps of 43563_A4) and Fig 5F (CNs and Treg Suppression Score of one tissue piece of 43563_A8)
# ======================================================================================================================
if _run('F5EF_section_maps'):
    import seaborn as sns
    import matplotlib.font_manager as fm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from scipy.ndimage import distance_transform_edt
    OUTD = f'{RJ}/Fig4A-F_treg_suppression_niches'
    BIN, SIG = 10, 5.0                                                   # 10 um grid, sigma 50 um for the cell-composition fields
    NW_SIG, CX_PCT = 200, 70                                             # CXCL12-High Region = upper 30% of the smoothed CXCL12 field
    MIN_AREA = 2 * np.pi * NW_SIG ** 2 / BIN ** 2                        # outline patches smaller than 2 pi sigma^2 are dropped
    CLEAN_R_UM, CLEAN_MIN_MM2, KNN, DOT_UM, DOT_AREA_X = 50, 0.05, 20, 250, 0.47
    C_TISSUE, C_APCDOT, LW_REG = '#ebebeb', '#08306b', 1.0
    APC_S_REF, IN_PER_MM_REF = 0.48, 3.28 / 7.61                         # density-matched APC dots
    C_BG, TREG_UM2 = '0.82', 90.0                                        # Fig 5F: background cells, 90 um Treg dots
    SCORE_LABEL = 'Treg\nSuppression\nScore (z,\n20-Treg mean)'
    SA = load_cells(['43563_A4', '43563_A8'])
    TL = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    vb = pd.read_csv(f'{OUTA}/fig4_cxcl12_region_tregs_exp.csv', index_col=0)['suppression score smoothed'].values
    LIM = np.nanpercentile(np.abs(vb - np.nanmean(vb)), 97); VLO, VHI = np.nanmean(vb) - LIM, np.nanmean(vb) + LIM

    def prep(LEAD):
        """Display layout of one section, per-piece smoothed fields (% of cells), CXCL12-High Region and 20-Treg means (tissue coordinates)."""
        sa = SA[SA.section == LEAD]
        with h5py.File(f'{R}/{FOLDER[LEAD]}/Processed/norm_log1p.h5ad', 'r') as h:
            ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
        nm = pd.Index(LEAD + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
        keep = nm.isin(sa.index)
        E = pd.DataFrame(X[np.where(keep)[0]][:, [v.get_loc('CXCL12')]].toarray(), index=nm[keep], columns=['CXCL12']); del X; gc.collect()
        cells = sa.loc[E.index].copy()
        ox, oy = cells.x.values.copy(), cells.y.values.copy()
        fx, fy = squeeze(ox), squeeze(oy); sx, sy = fx(ox), fy(oy)
        piece, SHIFT, _ = stitch(np.c_[sx, sy], LEAD, v8=True)
        shv = np.vstack([SHIFT[p] for p in piece]); cells['x'], cells['y'] = sx + shv[:, 0], sy + shv[:, 1]
        DISP = pd.DataFrame({'dx': cells.x.values - ox, 'dy': cells.y.values - oy}, index=cells.index)
        # display unit = (piece, translation): the gap squeeze can move parts of one piece by different amounts
        piece = pd.factorize(pd.Series(list(zip(piece, np.round(DISP.dx.values, 3), np.round(DISP.dy.values, 3)))))[0] + 1
        cells['piece'] = piece; DISP['piece'] = piece
        TR = {int(p_): (float(r.dx), float(r.dy)) for p_, r in DISP.groupby('piece')[['dx', 'dy']].mean().iterrows()}
        T = TL[TL.section == LEAD].copy(); T = T.join(sa[['x', 'y']], how='inner', rsuffix='_s'); T['ox'], T['oy'] = T.x_s.values, T.y_s.values
        T['piece'] = cells.piece.reindex(T.index).values
        T['x'] = T.ox + T.piece.map(lambda p: TR[int(p)][0]); T['y'] = T.oy + T.piece.map(lambda p: TR[int(p)][1])
        _, nb = cKDTree(T[['ox', 'oy']].values).query(T[['ox', 'oy']].values, k=min(KNN, len(T)))
        T['score_knn'] = np.nanmean(T.act_z.values[nb], axis=1); T['cxcr4_knn'] = np.nanmean(T.rec_CXCR4.values[nb], axis=1)
        # per-piece fields on the display grid; each piece is smoothed alone, contested pixels go to the piece with the nearest own cell
        x0, y0 = cells.x.min() - 3 * SIG * BIN, cells.y.min() - 3 * SIG * BIN
        nx, ny = int((cells.x.max() - x0) / BIN) + 4 * int(SIG) + 4, int((cells.y.max() - y0) / BIN) + 4 * int(SIG) + 4
        ix, iy = ((cells.x - x0) / BIN).astype(int).values, ((cells.y - y0) / BIN).astype(int).values
        XC, YC = x0 + (np.arange(nx) + 0.5) * BIN, y0 + (np.arange(ny) + 0.5) * BIN; EXT = (x0, x0 + nx * BIN, y0 + ny * BIN, y0)
        PIDS = sorted(TR); ASM, OCC = {}, {}
        for p_ in PIDS:
            a = np.zeros((ny, nx)); sel = piece == p_; np.add.at(a, (iy[sel], ix[sel]), 1); OCC[p_] = a > 0; ASM[p_] = gaussian_filter(a, SIG, mode='constant')
        AMAX = max(a.max() for a in ASM.values())
        MASK0 = {p_: ndimage.binary_fill_holes(ASM[p_] > 0.05 * AMAX) for p_ in PIDS}
        cover = sum(m.astype(int) for m in MASK0.values())
        DIST = np.stack([distance_transform_edt(~OCC[p_]) for p_ in PIDS]); DIST[np.stack([~MASK0[p_] for p_ in PIDS])] = np.inf
        OWN = np.argmin(DIST, axis=0); ANY = cover > 0
        MASK = {p_: ANY & (OWN == k) for k, p_ in enumerate(PIDS)}
        TISSUE = np.zeros((ny, nx), bool)
        for m in MASK.values(): TISSUE |= m

        def field(vals):
            out = np.full((ny, nx), np.nan)
            for p_ in PIDS:
                num = np.zeros((ny, nx)); sel = piece == p_; np.add.at(num, (iy[sel], ix[sel]), vals[sel].astype(float))
                f_ = gaussian_filter(num, SIG, mode='constant') / np.where(ASM[p_] > 0, ASM[p_], 1) * 100; out[MASK[p_]] = f_[MASK[p_]]
            return out

        def per_piece(m, fn):
            out = np.zeros((ny, nx), bool)
            for p_ in PIDS: out |= fn(m & MASK[p_]) & MASK[p_]
            return out

        def one(m):
            m = ndimage.binary_fill_holes(ndimage.binary_closing(m, np.ones((3, 3)), iterations=2)); lab_, n_ = ndimage.label(m, structure=np.ones((3, 3)))
            return np.r_[False, ndimage.sum(m, lab_, np.arange(1, n_ + 1)) >= MIN_AREA][lab_]
        rad = int(round(CLEAN_R_UM / BIN)); yy, xx = np.mgrid[-rad:rad + 1, -rad:rad + 1]; DISK = (xx ** 2 + yy ** 2) <= rad ** 2; MINB = int(round(CLEAN_MIN_MM2 * 1e6 / BIN ** 2))

        def clean_one(m):
            m = ndimage.binary_opening(ndimage.binary_closing(m, structure=DISK), structure=DISK); lab_, n_ = ndimage.label(m, structure=np.ones((3, 3)))
            return np.r_[False, ndimage.sum(m, lab_, np.arange(1, n_ + 1)) >= MINB][lab_]
        FCX = field(E.CXCL12.values > 0); FOGN = field((cells.ct == 'OGN+RSPO3+ Fib').values)
        CXd = per_piece(per_piece(np.nan_to_num(FCX) >= np.nanpercentile(FCX[TISSUE], CX_PCT), one), clean_one)
        return dict(cells=cells, T=T, EXT=EXT, XC=XC, YC=YC, TISSUE=TISSUE, FCX=FCX, FOGN=FOGN, CXd=CXd, MASK=MASK)

    # ---- figure helpers: one row of axes with fixed widths (inches), map axes and key axes alternate
    asp = lambda e: abs(e[1] - e[0]) / abs(e[3] - e[2])
    def grid(widths, h):
        fig, axs = plt.subplots(1, len(widths), figsize=(sum(widths), h), gridspec_kw=dict(width_ratios=widths, wspace=0, left=0, right=1, bottom=0, top=1))
        for k in axs[1::2]: k.axis('off')
        return fig, axs
    def frame(ax, d, title):
        E_ = d['EXT']; ax.set_xlim(E_[0], E_[1]); ax.set_ylim(E_[2], E_[3]); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title, fontsize=15, pad=6)
    def side_label(ax, t): ax.annotate(t, (0, 0.5), xycoords='axes fraction', xytext=(-24, 0), textcoords='offset points', rotation=90, ha='center', va='center', fontsize=15)
    def outline(ax, d, m, col, lw=LW_REG): ax.contour(d['XC'], d['YC'], m.astype(float), levels=[0.5], colors=col, linewidths=lw, zorder=6)
    def flat(ax, d, m, col, z=0): ax.imshow(np.where(m, 1.0, np.nan), extent=d['EXT'], cmap=ListedColormap([col]), interpolation='nearest', zorder=z)

    def scalebar(ax, d):                                                 # 1 mm bar in the first tissue-free corner
        m = ndimage.binary_dilation(d['TISSUE'], iterations=3); ny_, nx_ = m.shape; kw, kh = int(np.ceil(1300 / BIN)), int(np.ceil(450 / BIN))
        C_ = [n_ for n_, (r0, c0) in [('upper right', (0, nx_ - kw)), ('upper left', (0, 0)), ('lower right', (ny_ - kh, nx_ - kw)), ('lower left', (ny_ - kh, 0))]
              if not m[r0:r0 + kh, c0:c0 + kw].any()]
        c_ = C_[0] if C_ else 'lower left'; anc = {'upper left': (0, 1), 'upper right': (1, 1), 'lower left': (0, 0), 'lower right': (1, 0)}[c_]
        ax.add_artist(AnchoredSizeBar(ax.transData, 1000, '1 mm', loc=c_, bbox_to_anchor=anc, bbox_transform=ax.transAxes, pad=0, borderpad=0.1, sep=2, frameon=False,
                                      size_vertical=40, color='k', fontproperties=fm.FontProperties(size=11.5)))

    def vbar(kax, top, h, im, label, ext_=None, ticks=None):            # colourbar in the key axes (top, h: fractions of the key axes)
        cax = kax.inset_axes([0.08, top - h, 0.1, h]); cb = kax.figure.colorbar(im, cax=cax, orientation='vertical', **({'extend': ext_} if ext_ else {}))
        cb.ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        if ticks is not None: cb.set_ticks(ticks)
        else: cb.ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−'))); cb.ax.tick_params(labelsize=11.5, length=2.5, pad=1.5)
        cb.outline.set_linewidth(0.6); cax.set_title(label, fontsize=12, pad=5, loc='left', linespacing=1.15)
        lo_, hi_ = sorted(cax.get_ylim()); cax.set_yticks([t for t in cax.get_yticks() if lo_ - 1e-9 <= t <= hi_ + 1e-9])   # only ticks inside the range

    def vkey(kax, y, handles):
        kax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.05, y), frameon=False, fontsize=11.5, handlelength=1.1, handletextpad=0.35,
                   labelspacing=0.5, borderaxespad=0, borderpad=0)

    # ---------------- Fig 5E: four maps of 43563_A4, each with its colourbar and key on the right ----------------
    SEC = '43563_A4'; A = prep(SEC); T = A['T']; FOGN, FCX, TISSUE, CXd = A['FOGN'], A['FCX'], A['TISSUE'], A['CXd']
    C4LO, C4HI = np.percentile(T.cxcr4_knn, 2), np.percentile(T.cxcr4_knn, 98)
    CXLO, CX75, CXHI = np.nanpercentile(FCX[TISSUE], 50), np.nanpercentile(FCX[TISSUE], 75), np.nanpercentile(FCX[TISSUE], 98)
    V6_CMAP = LinearSegmentedColormap.from_list('v6_pink', [(0.0, '#ffffff'), (CXLO / CXHI, '#fde0ef'), (CX75 / CXHI, '#f768a1'), (1.0, '#ae017e')])
    CXKW = dict(cmap=V6_CMAP, norm=matplotlib.colors.Normalize(vmin=0, vmax=CXHI))
    plt.rcParams.update({'font.size': 14})
    CW, RH = 3.86, 3.86 / asp(A['EXT'])                                  # map width and height (inches)
    fig, axs = grid([CW, 1.6] * 4, RH); AX, KX = axs[0::2], axs[1::2]
    ipm = CW / abs(A['EXT'][1] - A['EXT'][0]) * 1000; APC_S = APC_S_REF * (ipm / IN_PER_MM_REF) ** 2   # inches per mm on the maps

    def cx_background(ax):
        im_ = ax.imshow(FCX, extent=A['EXT'], interpolation='nearest', zorder=1, **CXKW)
        ax.contour(A['XC'], A['YC'], TISSUE.astype(float), levels=[0.5], colors='0.72', linewidths=0.6, zorder=2); return im_
    ax = AX[0]; imO = ax.imshow(FOGN, extent=A['EXT'], cmap='viridis', vmin=0, vmax=np.nanpercentile(FOGN[TISSUE], 99), interpolation='nearest', zorder=1)
    outline(ax, A, CXd, 'white'); frame(ax, A, 'OGN⁺RSPO3⁺ Fibroblasts')
    ax = AX[1]; imC = cx_background(ax); cc = A['cells']; ap_ = cc.ct.isin(APC).values
    ax.scatter(cc.x.values[ap_], cc.y.values[ap_], s=APC_S, c=C_APCDOT, lw=0, rasterized=True, zorder=3)
    outline(ax, A, CXd, '#333333'); frame(ax, A, 'Antigen-Presenting Cells\n(DC, Macrophage, B Cell)')
    ax = AX[2]; cx_background(ax); o4 = np.argsort(T.cxcr4_knn.values)
    im4 = ax.scatter(T.x.values[o4], T.y.values[o4], c=T.cxcr4_knn.values[o4], cmap='YlOrBr', vmin=C4LO, vmax=C4HI, s=(DOT_UM / 1000 * ipm * 72) ** 2, edgecolors='0.15', linewidths=0.4, rasterized=True, zorder=3)
    outline(ax, A, CXd, '#333333'); frame(ax, A, 'Treg CXCR4\n(The CXCL12 Receptor)')
    ax = AX[3]; flat(ax, A, TISSUE, C_TISSUE, z=0)
    DOT_S = 0.55 * DOT_AREA_X * (DOT_UM / 1000 * ipm * 72) ** 2
    imS = ax.scatter(T.x, T.y, c=T.score_knn, cmap='RdBu_r', vmin=VLO, vmax=VHI, s=DOT_S, lw=0, alpha=0.85, rasterized=True, zorder=2); frame(ax, A, 'Treg Suppression Score')
    scalebar(AX[0], A); side_label(AX[0], SEC)
    CBT = 1 - 0.75 / RH; KY = CBT - 0.42 - 0.30 / RH                     # colourbars 0.75 in below the map top, keys 0.30 in below them
    vbar(KX[0], CBT, 0.42, imO, 'OGN⁺RSPO3⁺\n(% of cells)')
    vbar(KX[1], CBT, 0.42, imC, 'CXCL12\n(% of cells\nexpressing)', 'max', [t for t in np.arange(0, CXHI + 1e-9, 20)])
    vbar(KX[2], CBT, 0.42, im4, 'Treg CXCR4\n(z, 20-Treg\nmean)', 'both')
    vbar(KX[3], CBT, 0.42, imS, SCORE_LABEL, 'both')
    vkey(KX[0], KY, [Line2D([], [], color='0.6', lw=1.8, label='CXCL12-High\nRegion (white\noutline)')])
    vkey(KX[1], KY, [Line2D([], [], ls='', marker='o', ms=4.5, mfc=C_APCDOT, mec='none', label='APC'), Line2D([], [], color='#333333', lw=LW_REG, label='CXCL12-High\nRegion')])
    vkey(KX[2], KY, [Line2D([], [], color='#333333', lw=LW_REG, label='CXCL12-High\nRegion')])
    fig.savefig(f'{OUTD}/fig4_A4_four_panel_plain.png', dpi=300, bbox_inches='tight', pad_inches=0.06); plt.close(fig)

    # ---------------- Fig 5F: the large lower-left display unit of 43563_A8 alone: CNs | Treg Suppression Score ----------------
    SEC = '43563_A8'; d = prep(SEC); cf, Tf = d['cells'], d['T']; uu = cf.piece.values.astype(int)
    CNS = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas().CN_k10
    PAL = dict(zip(range(1, 11), sns.color_palette('tab10', 10).as_hex()))
    NAMES = {i + 1: n for i, n in enumerate(pd.read_csv(f'{CND}/cn_k10_celltype_log2_enrichment.csv', index_col=0).index)}   # full CN names
    UT = pd.DataFrame({'unit': uu, 'x': cf.x.values, 'y': cf.y.values}).groupby('unit').agg(n_cells=('x', 'size'), x_um=('x', 'mean'), y_um=('y', 'mean'))
    PU = int(UT[(UT.x_um < np.median(cf.x.values)) & (UT.y_um > np.median(cf.y.values))].n_cells.idxmax())   # y grows downward on the maps
    tp_ = Tf[Tf.piece.astype(int) == PU]
    M_ = d['MASK'][PU]; r_, c_ = np.nonzero(M_); PADp = 3 * int(SIG); ny_, nx_ = M_.shape
    r0, r1, c0, c1 = max(r_.min() - PADp, 0), min(r_.max() + PADp + 1, ny_), max(c_.min() - PADp, 0), min(c_.max() + PADp + 1, nx_)
    x0f, y0f = d['EXT'][0], d['EXT'][3]
    LP = dict(EXT=(x0f + c0 * BIN, x0f + c1 * BIN, y0f + r1 * BIN, y0f + r0 * BIN), TISSUE=M_[r0:r1, c0:c1], cells=cf[uu == PU], T=tp_)
    plt.rcParams.update({'font.size': 14})
    Hp = 6.0; MWp = asp(LP['EXT']) * Hp
    fig, (A1, _g, A2, K2) = grid([MWp, 0.40, MWp, 1.75], Hp)
    ipm = MWp / abs(LP['EXT'][1] - LP['EXT'][0]) * 1000
    cc_ = LP['cells']; cn_ = CNS.reindex(cc_.index).values; s_cn = (9 / 1000 * ipm * 72) ** 2   # 9 um CN dots
    for k in range(1, 11): m_ = cn_ == k; A1.scatter(cc_.x.values[m_], cc_.y.values[m_], s=s_cn, c=PAL[k], rasterized=True, lw=0, alpha=1.0)
    frame(A1, LP, 'Cellular Neighborhoods\n(k = 10)')
    nt = ~cc_.index.isin(tp_.index) & (cc_.ct.values != 'Treg')                                 # non-Treg cells, light grey
    A2.scatter(cc_.x.values[nt], cc_.y.values[nt], s=s_cn, c=C_BG, rasterized=True, lw=0, alpha=1.0, zorder=1)
    imS = A2.scatter(tp_.x, tp_.y, c=tp_.score_knn, cmap='RdBu_r', vmin=VLO, vmax=VHI, s=(TREG_UM2 / 1000 * ipm * 72) ** 2, lw=0, alpha=0.85, rasterized=True, zorder=2)
    frame(A2, LP, 'Treg Suppression Score')
    for a_ in (A1, A2): scalebar(a_, LP)
    side_label(A1, SEC)
    fig.legend(handles=[Line2D([], [], ls='', marker='o', ms=8, mfc=PAL[k], mec='none', label=NAMES[k]) for k in range(1, 11)], ncol=2, loc='upper left',
               bbox_to_anchor=(0, -0.20 / Hp), frameon=False, fontsize=11.5, handletextpad=0.2, labelspacing=0.32, columnspacing=1.6, borderaxespad=0, borderpad=0)
    vbar(K2, 1 - 0.75 / Hp, 0.40, imS, SCORE_LABEL, 'both')
    fig.savefig(f'{OUTD}/fig4_A8_cn_treg_two_panel_piece.png', dpi=300, bbox_inches='tight', pad_inches=0.06); plt.close(fig)
    print('Fig 5F: display unit', PU, '|', len(cc_), 'cells,', len(tp_), 'Tregs', flush=True)


# ======================================================================================================================
# Section S9F_apc_cxcl12_quintiles: score by local APC quintile in CXCL12-Low vs CXCL12-High Regions, and its rise (Supp 9F)
# ======================================================================================================================
if _run('S9F_apc_cxcl12_quintiles'):
    COL_HIGH, COL_LOW = '#F7B325', '#1a9850'
    T = pd.read_csv(f'{LR}/treg_local_ligand_receptor.csv.gz', index_col=0)
    T = T.join(pd.read_csv(f'{LR}/treg_local_cxcl12_producer_exp.csv.gz', index_col=0), how='left')
    cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); cnt = cnt.loc[cnt.index.intersection(T.index)]
    T = T.loc[cnt.index]; T['cz'] = np.log1p(cnt.sum(1).values); del cnt; gc.collect()   # Treg order of the count table (rank ties depend on it)
    T['cz'] = T.groupby('section').cz.transform(lambda v: (v - v.mean()) / v.std())
    T['apc_q'] = T.groupby('section').apc_frac.transform(lambda v: pd.qcut(v, 5, labels=False, duplicates='drop'))
    T['t_exp'] = T.groupby('section').lig_CXCL12_exp.transform(lambda v: pd.qcut(v.rank(method='first'), 3, labels=['low', 'mid', 'high']))
    # per quintile: CXCL12-High vs -Low Tregs, Mann-Whitney, BH over quintiles
    sub = T[T.t_exp.isin(['low', 'high']) & T.apc_q.notna()]; out = []
    for q, g in sub.groupby('apc_q'):
        hi, lo = g[g.t_exp == 'high'].act_z.dropna(), g[g.t_exp == 'low'].act_z.dropna()
        if len(hi) < 20 or len(lo) < 20: continue
        out.append(dict(quintile=int(q) + 1, mean_frac=g.apc_frac.mean(), suppression_cxcl12_high=hi.mean(), suppression_cxcl12_low=lo.mean(),
                        sem_high=hi.std() / np.sqrt(len(hi)), sem_low=lo.std() / np.sqrt(len(lo)), p=stats.mannwhitneyu(hi, lo).pvalue))
    Q = pd.DataFrame(out); Q['q'] = multipletests(Q.p, method='fdr_bh')[1]; Q = Q.set_index('quintile')
    rz = {k: (Q[f'suppression_cxcl12_{k}'][5] - Q[f'suppression_cxcl12_{k}'][1], np.sqrt(Q[f'sem_{k}'][1] ** 2 + Q[f'sem_{k}'][5] ** 2)) for k in ['low', 'high']}
    # difference in rise: Q1 + Q5 Tregs, act_z ~ hi * top + cz; LMM with section random intercept, and OLS with section-clustered SEs (t, G - 1 df)
    d = T[T.t_exp.isin(['low', 'high']) & T.apc_q.isin([0, 4])].copy()
    d['hi'] = (d.t_exp == 'high').astype(float); d['top'] = (d.apc_q == 4).astype(float); d = d[['section', 'act_z', 'hi', 'top', 'cz']].dropna()
    mm = smf.mixedlm('y ~ hi * top + cz', d.rename(columns={'act_z': 'y'}), groups=d.section).fit(reml=True, method='lbfgs')
    d['hixtop'] = d.hi * d.top; FC = 'act_z ~ hi + top + hixtop + cz'
    o = smf.ols(FC, d).fit(cov_type='cluster', cov_kwds=dict(groups=pd.factorize(d.section)[0]), use_t=True)
    CL_P = 2 * stats.t.sf(abs(o.params['hixtop'] / o.bse['hixtop']), d.section.nunique() - 1)
    sg = np.sign(mm.params['hi:top']); AG, NS = 0, 0
    for _, g in d.groupby('section'):
        if g.hixtop.std() == 0 or len(g) < 5: continue
        NS += 1; AG += int(np.sign(smf.ols(FC, g).fit().params['hixtop']) == sg)
    print('rise CXCL12-Low {:+.3f} +- {:.3f} | High {:+.3f} +- {:.3f} | LMM P = {:.3g} (n {:,}) | clustered P = {:.3g}, {}/{} sections'.format(
        rz['low'][0], rz['low'][1], rz['high'][0], rz['high'][1], mm.pvalues['hi:top'], len(d), CL_P, AG, NS), flush=True)

    def pfmt(p):
        if p >= 0.001: return f'{p:.2g}'
        ex = int(np.floor(np.log10(p))); m = p / 10 ** ex
        if m >= 9.95: m, ex = 1.0, ex + 1
        return '{:.1f} × 10$^{{{}}}$'.format(m, ex)
    stars = lambda q: 'n.s.' if q >= 0.05 else ('*' if q >= 0.01 else ('**' if q >= 0.001 else '***'))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 15, 'axes.linewidth': 1.2, 'axes.spines.top': False, 'axes.spines.right': False})
    EKW = dict(ecolor='k', elinewidth=1.2, capsize=3, capthick=1.2)
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(8.7, 5.4), gridspec_kw=dict(width_ratios=[4.6, 1.75], wspace=0.48))
    C = Q.reset_index(); xs = np.arange(len(C)); w = 0.38
    for k, dx, col, lab in [('low', -w / 2, COL_LOW, 'CXCL12-Low Region'), ('high', w / 2, COL_HIGH, 'CXCL12-High Region')]:
        ax.bar(xs + dx, C[f'suppression_cxcl12_{k}'], w, yerr=C[f'sem_{k}'], color=col, edgecolor='k', lw=0.7, label=lab, error_kw=EKW, zorder=2)
    up, dn = np.maximum(C.suppression_cxcl12_low + C.sem_low, C.suppression_cxcl12_high + C.sem_high), np.minimum(C.suppression_cxcl12_low - C.sem_low, C.suppression_cxcl12_high - C.sem_high)
    ymax, ymin = max(up.max(), 0.0), min(dn.min(), 0.0)
    for xi, top, q in zip(xs, np.maximum(up, 0.0), C.q):
        if q < 0.05: ax.text(xi, top + 0.02 * (ymax - ymin), stars(q), ha='center', va='bottom', fontsize=16)
    ax.axhline(0, color='0.3', lw=1, zorder=1); ax.set_xticks(xs); ax.set_xticklabels([f'{v*100:.0f}%' for v in C.mean_frac], fontsize=13.5)
    ax.set_ylim(ymin - 0.06 * (ymax - ymin), ymax + 0.30 * (ymax - ymin))
    ax.set_xlabel('Antigen-Presenting Cells within 50 μm\n(Quintile of Local APC Fraction)', fontsize=14.5)
    ax.set_ylabel('Treg Suppression Score (z)', fontsize=15, labelpad=1)
    ax.set_title('Treg Suppression Score\nby Local APC Fraction', fontsize=16)
    ax.legend(frameon=False, fontsize=12.5, loc='upper left', handlelength=1.1, borderpad=0.1)
    rv, re_ = np.array([rz['low'][0], rz['high'][0]]), np.array([rz['low'][1], rz['high'][1]])
    axr.bar([0, 1], rv, 0.62, yerr=re_, color=[COL_LOW, COL_HIGH], edgecolor='k', lw=0.7, error_kw=EKW, zorder=2)
    rtop, rbot = (rv + re_).max(), min(0.0, (rv - re_).min()); rng = rtop - rbot; yb = rtop + 0.07 * rng
    axr.plot([0, 0, 1, 1], [yb - 0.03 * rng, yb, yb, yb - 0.03 * rng], color='k', lw=1.1)
    axr.text(0.5, yb + 0.03 * rng, 'Difference\n$P$ = ' + pfmt(CL_P) + ('' if CL_P < 0.05 else ' (n.s.)') + f';\n{AG}/{NS} sections', ha='center', va='bottom', fontsize=12.5)
    axr.set_ylim(rbot - 0.04 * rng if rbot < 0 else 0, rtop + 0.42 * rng); axr.axhline(0, color='0.3', lw=1, zorder=1)
    axr.set_xlim(-0.6, 1.6); axr.set_xticks([0, 1]); axr.set_xticklabels(['Low', 'High'], fontsize=13.5)
    axr.set_xlabel('CXCL12\nRegion', fontsize=14.5)
    axr.set_ylabel('Rise in Treg Suppression Score,\nLowest → Highest Quintile (SD)', fontsize=13.5, labelpad=2)
    axr.set_title('Rise', fontsize=16)
    fig.savefig(f'{RJ}/SuppFig9A-F_niche_cxcl12_apc/fig4_apc_quintile_cxcl12_bars_with_rise.png', dpi=300, bbox_inches='tight', pad_inches=0.06); plt.close(fig)


# ======================================================================================================================
# Section S9E_cxcl12_by_celltype: CXCL12 by cell type (Supp 9E, figure + CSV) and APC share of cells per CN (table read by 17)
# ======================================================================================================================
if _run('S9E_cxcl12_by_celltype'):
    from statsmodels.tools.sm_exceptions import ConvergenceWarning; warnings.simplefilter('ignore', ConvergenceWarning)
    STEM1 = f'{RJ}/SuppFig9A-F_niche_cxcl12_apc/supp_cxcl12_expression_by_cell_type'; STEM2 = f'{P4}/archive/supp_apc_abundance_by_neighbourhood'
    EXCL = ['SELENOP+ Fib']                                      # not shown, not tested (SELENOP+ Mac is shown)
    OGN = 'OGN+RSPO3+ Fib'; HI = '#b2182b'; HIL = '#e7a3a3'; GR = '#7f7f7f'; GRL = '#c4c4c4'
    DISP = {'OGN+RSPO3+ Fib': 'OGN⁺RSPO3⁺ Fib', 'ADAMDEC1+ Fib': 'ADAMDEC1⁺ Fib', 'PLA2G2A+ ECM Fib': 'PLA2G2A⁺ ECM Fib', 'CD9+ M2 Mac': 'CD9⁺ M2 Mac',
            'Absoprtive Epi': 'BEST4⁺ Epi', 'Gd T': 'γδ T', 'SELENOP+ Mac': 'SELENOP⁺ Mac'}
    APC5 = ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac']
    APC_GRP = {'DC': 'DC', 'B cell': 'B cell', 'CD9+ M2 Mac': 'Macrophage', 'Mo-Mac': 'Macrophage', 'SELENOP+ Mac': 'Macrophage'}
    APCG = ['DC', 'B cell', 'Macrophage']
    MIN_CELLS = 20; MIN_CN_CELLS = 100
    sig3 = lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)

    def ci95(v):
        v = np.asarray(v, float); v = v[np.isfinite(v)]; n = len(v)
        if n < 2: return (v.mean() if n else np.nan), np.nan, np.nan, n
        h = stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n); return v.mean(), v.mean() - h, v.mean() + h, n
    # ---- per-cell CXCL12 and CN (CACHE80=<path.pkl> stores / re-uses this table)
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
            e = e.loc[e.index.intersection(sub.index)]
            cs = sub.loc[e.index].copy(); cs['expr'] = e.values; parts.append(cs)
        C = pd.concat(parts); del parts, sa; gc.collect()
        cn = cn.reindex(C.index); C['CN'] = cn.CN_k10.values; C['cn_cell_type'] = cn.cell_type.values
        if CACHE: pd.to_pickle(C, CACHE)
    C['pos'] = C.expr > 0
    # ---- CXCL12 by cell type: % expressing and level in expressing cells, per section; OGN+RSPO3+ Fib vs each type, paired t-test, BH
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
    rows = []
    for ct in [c for c in pool.index if c not in EXCL]:
        m, lo, hi, n = ci95(pv[ct]) if ct in pv else (np.nan,) * 3 + (0,)
        m2, lo2, hi2, n2 = ci95(PSu[PSu.ct == ct].mean_expr_in_expressing)
        rows.append(dict(cell_type_general=ct, label=DISP.get(ct, ct), n_cells=int(pool.loc[ct, 'n_cells']), n_expressing=int(pool.loc[ct, 'n_expressing']),
                         pooled_pct_expressing=pool.loc[ct, 'pct_expressing'], n_sections=n, mean_pct_expressing=m, ci_low=lo, ci_high=hi,
                         fold_over_all_cells=m / all_m, pooled_mean_expr_in_expressing=pool.loc[ct, 'mean_expr_in_expressing'],
                         mean_expr_in_expressing=m2, expr_ci_low=lo2, expr_ci_high=hi2,
                         n_sections_ranked_top=int((rank[ct] == 1).sum()) if ct in rank else 0,
                         median_rank_across_sections=rank[ct].median() if ct in rank else np.nan))
    T1 = pd.DataFrame(rows).sort_values('mean_pct_expressing', ascending=False).reset_index(drop=True)
    tst = []
    for ct in T1.cell_type_general:
        if ct == OGN: tst.append((np.nan, np.nan, np.nan)); continue
        d = pv[[OGN, ct]].dropna(); t, p = stats.ttest_rel(d[OGN], d[ct])
        tst.append(((d[OGN] - d[ct]).mean(), p, f'{int((d[OGN] > d[ct]).sum())}/{len(d)}'))
    T1['ogn_minus_type_pct'] = [x[0] for x in tst]; T1['p_paired_t_vs_ogn'] = [x[1] for x in tst]
    T1['n_sections_ogn_higher'] = [x[2] for x in tst]
    m_ = T1.p_paired_t_vs_ogn.notna(); T1['q_vs_ogn'] = np.nan
    T1.loc[m_, 'q_vs_ogn'] = multipletests(T1.loc[m_, 'p_paired_t_vs_ogn'], method='fdr_bh')[1]
    T1['rank'] = np.arange(1, len(T1) + 1)
    lv = PSu.pivot(index='section', columns='ct', values='mean_expr_in_expressing'); pl = []
    for ct in T1.cell_type_general:
        if ct == OGN: pl.append(np.nan); continue
        d = lv[[OGN, ct]].dropna(); pl.append(stats.ttest_rel(d[OGN], d[ct])[1])
    T1['p_paired_t_level_vs_ogn'] = pl; T1['q_level_vs_ogn'] = np.nan
    T1.loc[m_, 'q_level_vs_ogn'] = multipletests(T1.loc[m_, 'p_paired_t_level_vs_ogn'], method='fdr_bh')[1]
    lv_hi = T1[(T1.q_level_vs_ogn < 0.05) & (T1.mean_expr_in_expressing < T1.mean_expr_in_expressing.iloc[0])]
    n_top = int(T1.loc[T1.cell_type_general == OGN, 'n_sections_ranked_top'].iloc[0]); n_sec = int(T1.loc[T1.cell_type_general == OGN, 'n_sections'].iloc[0])
    ALLR = pd.DataFrame([dict(cell_type_general='ALL CELLS (reference; includes SELENOP+ Fib)', label='All Cells', n_cells=len(C), n_expressing=int(C.pos.sum()),
                              pooled_pct_expressing=tot_pct, n_sections=len(sec_all), mean_pct_expressing=all_m, ci_low=all_lo, ci_high=all_hi, fold_over_all_cells=1.0,
                              pooled_mean_expr_in_expressing=C.expr[C.pos].mean())])
    PSo = PS.rename(columns={'ct': 'cell_type_general'}).assign(row_type='section', used_in_mean=lambda d: d.n_cells >= MIN_CELLS)
    OUT1 = pd.concat([T1.assign(row_type='summary'), ALLR.assign(row_type='summary'), PSo], ignore_index=True)
    lead = ['row_type', 'rank', 'cell_type_general', 'label', 'section']
    OUT1[lead + [c for c in OUT1 if c not in lead]].apply(sig3).to_csv(f'{STEM1}.csv', index=False)
    # ---- figure (Supp 9E)
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    n = len(T1); Y = {c: i for i, c in enumerate(T1.cell_type_general)}
    FH = 0.315 * n + 1.9
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.9, FH), sharey=True, gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.13))
    rng = np.random.default_rng(1)
    for ct, i in Y.items():
        hi_ = ct == OGN; col, colL = (HI, HIL) if hi_ else (GR, GRL)
        r = T1.iloc[i]; v = PSu[PSu.ct == ct]
        for a_, val, m, lo, hi in [(ax, v.pct_expressing, r.mean_pct_expressing, r.ci_low, r.ci_high),
                                   (bx, v.mean_expr_in_expressing, r.mean_expr_in_expressing, r.expr_ci_low, r.expr_ci_high)]:
            a_.barh(i, m, height=0.72, color=colL, edgecolor='none', zorder=1)
            a_.scatter(val, i + rng.uniform(-0.2, 0.2, len(val)), s=10, color=col, alpha=0.75, linewidth=0, zorder=2)
            a_.plot([lo, hi], [i, i], color='black' if not hi_ else HI, lw=1.8, zorder=3, solid_capstyle='butt')
    ax.axvline(all_m, color='#252525', ls='--', lw=1.4, zorder=0)
    ax.text(all_m + 1, n - 0.6, 'All Cells', ha='left', va='bottom', fontsize=14, color='#252525')
    ax.set_yticks(range(n)); ax.set_yticklabels(T1.label, fontsize=13.5); ax.set_ylim(n - 0.4, -0.6)
    for tl, ct in zip(ax.get_yticklabels(), T1.cell_type_general):
        if ct == OGN: tl.set_color(HI); tl.set_fontsize(15)
    ax.set_xlabel('Cells Expressing CXCL12 (%)', fontsize=15.5); bx.set_xlabel('Mean CXCL12\n(log-normalised)', fontsize=15.5)
    xm = np.nanmax(PSu.pct_expressing); ax.set_xlim(0, xm * 1.04); ax.set_xticks(np.arange(0, xm * 1.04, 20))
    bx.set_xlim(0, np.nanmax(PSu.mean_expr_in_expressing) * 1.04); bx.set_xticks([0, 0.5, 1.0]); bx.set_xticklabels(['0', '0.5', '1.0'])
    for a_ in (ax, bx): a_.spines[['top', 'right']].set_visible(False); a_.tick_params(axis='x', labelsize=13.5); a_.margins(y=0)
    bx.tick_params(axis='y', length=0)
    ax.set_title('% Expressing CXCL12', loc='left', fontsize=16, pad=5)
    bx.set_title('Level in CXCL12-\nExpressing Cells', loc='left', fontsize=16, pad=5)
    ncomp = int(T1.q_vs_ogn.notna().sum())
    hi_pct = T1[(T1.q_vs_ogn < 0.05) & (T1.ogn_minus_type_pct > 0)]; n_hi_pct = len(hi_pct); q_hi_max = hi_pct.q_vs_ogn.max()
    o = T1.iloc[0]
    TXA = (f'OGN⁺RSPO3⁺ Fib: {o.mean_pct_expressing:.0f}%\n({o.mean_pct_expressing / all_m:.1f}× all cells)\n'
           f'Top in {n_top}/{n_sec} sections\nHigher than {n_hi_pct}/{ncomp} types,\nq {"< 0.001" if q_hi_max < 0.001 else f"≤ {q_hi_max:.2g}"}')
    TXB = f'Level: higher than\n{len(lv_hi)}/{int(m_.sum())} types, q < 0.05'
    h = [Line2D([], [], marker='o', ls='', color=GR, ms=4.5, label='One Section'), Line2D([], [], color='black', lw=1.8, label='Mean ± 95% CI'),
         Line2D([], [], color='#252525', ls='--', lw=1.4, label='All Cells (Mean of Sections)')]
    fig.suptitle('CXCL12 Expression by Cell Type', fontsize=18, x=0.01, ha='left', y=0.995)
    fig.legend(handles=h, loc='upper left', bbox_to_anchor=(0.01, 0.975), ncol=3, frameon=False, fontsize=13.5, handletextpad=0.4, columnspacing=1.4)
    fig.subplots_adjust(top=1 - 1.3 / FH, bottom=0.9 / FH, left=0.21, right=0.99)
    # the two statistics boxes, stacked in the lower part of panel a, right of the All Cells line
    for t_, yb in ((TXA, 0.135), (TXB, 0.045)):
        ax.text(all_m + 2.5, yb, t_, transform=ax.get_xaxis_transform(), ha='left', va='bottom', fontsize=13.5, color=HI, zorder=6, linespacing=1.3,
                bbox=dict(boxstyle='round,pad=0.45', facecolor='white', edgecolor=HI, lw=1))
    fig.savefig(f'{STEM1}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    # ---- APC share of all cells per CN; mixed model APC % ~ C(CN, ref CN10) + (1 | section), BH over 9 contrasts
    CC = C[C.CN.notna()].copy(); CC['CN'] = CC.CN.astype(int)
    apc_ct = CC.ct.where(CC.ct.isin(APC5))
    CC['apc'] = apc_ct; CC['apcg'] = apc_ct.map(APC_GRP)
    cnt = CC.groupby(['section', 'CN']).size().rename('n_cells').to_frame()
    for a in APCG: cnt[f'n_{a}'] = CC[CC.apcg == a].groupby(['section', 'CN']).size().reindex(cnt.index).fillna(0).astype(int)
    for a in APC5: cnt[f'n_celltype_{a}'] = CC[CC.apc == a].groupby(['section', 'CN']).size().reindex(cnt.index).fillna(0).astype(int)
    cnt['n_APC'] = cnt[[f'n_{a}' for a in APCG]].sum(1)
    cnt['apc_pct'] = 100 * cnt.n_APC / cnt.n_cells
    for a in APCG: cnt[f'pct_{a}'] = 100 * cnt[f'n_{a}'] / cnt.n_cells
    cnt = cnt.reset_index(); cnt['used'] = cnt.n_cells >= MIN_CN_CELLS
    U = cnt[cnt.used].copy()
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [('OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]   # CN labels of the Fig 5B map

    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {k: short_cn(v) for k, v in NM.items()}
    MOD = smf.mixedlm('apc_pct ~ C(CN, Treatment(reference=10))', U, groups=U.section)
    try: mm = MOD.fit(reml=True, method='lbfgs')
    except np.linalg.LinAlgError: mm = MOD.fit(reml=True, method='powell')   # lbfgs hits a singular matrix on these data
    ci = mm.conf_int(); st = []
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
    OUT2 = pd.concat([T2.assign(row_type='summary', section='all'),
                      cnt.assign(row_type='section', CN_label=cnt.CN.map(LAB)).rename(columns={'apc_pct': 'mean_apc_pct', 'used': 'used_in_mean'})], ignore_index=True)
    lead = ['row_type', 'CN', 'CN_label', 'section']
    OUT2[lead + [c for c in OUT2 if c not in lead]].apply(sig3).to_csv(f'{STEM2}.csv', index=False)
