#!/usr/bin/env python
"""CD CosMx 6K resections: figures of the ligand+ sender -> receptor+ strict Treg proximity (no statistic recomputed).

Panels: Fig 6C (per-section O/E dot plot of the 8 co-stimulatory/suppressive receptor module pairs, log axis, median dot, BH stars),
        Fig 6D (whole section 43563_A4: CD80+ Mac -> CTLA4+ and CD80+ Mac -> CD28+ strict Tregs within 30 um, CN tint),
        Supp 11B (per-pair bars: sections with O/E > 1, observed pairs; CN of the Treg in each pair vs all strict Tregs),
        Supp 11C (whole section 54033_B4, same two pairs).
Adapted from revision scripts Fig6C_spatial_lr_dotplot_maps.py and SuppFig10B_spatial_lr_support_niche_maps.py.
Rows: perm_q_bh_treg ascending, then sections with O/E > 1 descending, then sections with zero pairs ascending, then observed pairs
  descending. Maps: tissue pieces (50 um occupancy grid, 1-bin dilation, 8-neighbour components) packed by translation only; pairs
  found in the original coordinates. Supp 11C section = highest mean O/E (CTLA4, CD28) among sections with >= 20 pairs for both.
Inputs: SuppFig10B_spatial_lr_pairspecific_stats.csv and cache/spatial_lr/lr_genes_lognorm_allcells.parquet (18); resection
  combined/all_samples.parquet; treg_definition_comparison_per_cell.csv (10); cellular_neighborhoods_v2/ (11, 12).
Outputs: Fig6_treg_antitnf_taurus/Fig6C_spatial_lr_dotplot_maps.png, .csv (Fig 6C-D);
  SuppFig11_treg_ligand_senders_spatial_lr/SuppFig11B_spatial_lr_support_niche_maps.png, .csv and SuppFig10B_pair_cn_counts.csv (Supp 11B-C).
Sections (after 18): F6CD_dotplot_maps, S11BC_support_niche_maps.
Run: python 19_spatial_lr_figures.py, or RUN_SECTIONS=<id>[,<id>...] for single sections (file order).
"""
import os
import numpy as np, pandas as pd, pyarrow.parquet as pq, seaborn as sns
from scipy import ndimage
from scipy.spatial import cKDTree
from scipy.signal import fftconvolve
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple

_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


S10 = '/share/fsmresfiles/UC/revision_jenny/SuppFig11_treg_ligand_senders_spatial_lr'
FIG6 = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus'
CACHE = f'{S10}/cache/spatial_lr'
R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
P4 = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related'
CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
UM, RAD, LO, HI, BIN, GAPB = 0.12028, 30, 0.4, 12.0, 50, 2
PADUM = 80                                                 # map margin (um)

# ---------------- shared: the 8 primary pairs at 30 um (strict Treg receivers), row order, medians ----------------
d = pd.read_csv(f'{S10}/SuppFig10B_spatial_lr_pairspecific_stats.csv')
px = d[(d.table == 'proximity') & (d.radius_um == 30)]
POOL = px[(px.row == 'pooled') & (px.role == 'primary')].set_index('pair'); SEC = px[(px.row == 'section') & (px.role == 'primary')]
PAIRS8 = ['CD86->CD28', 'CD86->CTLA4', 'CD80->CD28', 'CD80->CTLA4', 'IL10->IL10RA', 'CD274->PDCD1', 'TNFSF18->TNFRSF18', 'PDCD1LG2->PDCD1']
SENDER = {'Mac': 'Mac', 'Stromal': 'Str'}
ALIAS = {'CD274': 'PD-L1', 'PDCD1LG2': 'PD-L2', 'TNFSF18': 'GITRL', 'TNFRSF18': 'GITR', 'PDCD1': 'PD-1'}
al = lambda g: ALIAS.get(g, g)
_g = SEC[SEC.pair.isin(PAIRS8)].groupby('pair')[['obs_total_treg', 'exp_total_treg']].sum()
NOBS = _g.obs_total_treg.round().astype(int)               # observed pairs within 30 um, summed over the 11 sections
MED = {}                                                   # median of the 11 per-section O/E (sections with no observed pair = 0)
for pr in PAIRS8:
    s = SEC[SEC.pair == pr].sort_values('section'); x = s.obs_exp_treg.values.astype(float)
    MED[pr] = float(np.median(x))
NZERO = SEC[SEC.pair.isin(PAIRS8)].groupby('pair').obs_total_treg.apply(lambda x: int((x == 0).sum()))
SIG = [p for p in PAIRS8 if POOL.loc[p].perm_q_bh_treg < 0.05]
ROWS = (POOL.loc[PAIRS8].assign(ngt=lambda t: t.sections_obs_exp_gt1_treg, nzero=lambda t: NZERO.loc[t.index], nobs=lambda t: NOBS.loc[t.index])
        .sort_values(['perm_q_bh_treg', 'ngt', 'nzero', 'nobs'], ascending=[True, False, True, False]).index.tolist())
Q8 = pd.Series(multipletests(POOL.loc[ROWS].perm_p_treg, method='fdr_bh')[1], index=ROWS)
stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
lab = lambda pr: f'{al(pr.split("->")[0])} ({SENDER[POOL.loc[pr].sender]}) → {al(pr.split("->")[1])}'
PAL = dict(zip(range(1, 11), sns.color_palette('tab10', 10).as_hex()))   # CN_k10 colours (as the Fig 5 niche map)


def cn_names():
    """Short CN_k10 names (CN1-CN10)."""
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [('OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]

    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    return {k: short_cn(v) for k, v in NM.items()}


def load_cells(genes):
    """Resection cells (um coordinates), the cached log-norm genes, strict Tregs, Mac senders and CN_k10."""
    sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'patient_sec', 'cell_type_general',
                                                                      '__index_level_0__']).to_pandas()
    sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
    sa['x'], sa['y'] = sa.CenterX_global_px.values * UM, sa.CenterY_global_px.values * UM
    E = pd.read_parquet(f'{CACHE}/lr_genes_lognorm_allcells.parquet', columns=genes).loc[sa.index]
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0, usecols=[0, 1]); TREG = PC.index[PC.cohort == 'resection']
    sa['treg'] = sa.index.isin(TREG); sa['mac'] = sa.cell_type_general.astype(str).isin(['SELENOP+ Mac', 'CD9+ M2 Mac']) & ~sa.treg
    CNc = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas().CN_k10
    sa['CN'] = CNc.reindex(sa.index).values
    return sa, E


def pairs(cs, lig, rec):
    """ligand+ Mac / receptor+ strict Treg masks and all (sender, receiver) positional index pairs within RAD um (original coordinates)."""
    XY = cs[['x', 'y']].values; Es = E.loc[cs.index]
    l = cs.mac.values & (Es[lig].values > 0); r = cs.treg.values & (Es[rec].values > 0); ir = np.flatnonzero(r)
    nb = cKDTree(XY[r]).query_ball_point(XY[l], r=RAD)
    return l, r, np.array([(i, ir[j]) for i, js in zip(np.flatnonzero(l), nb) for j in js], dtype=int).reshape(-1, 2)


def pack(cs):
    """Pieces = 50 um occupancy grid, 1-bin dilation, 8-neighbour components. Packing (translation only, in whole 50 um bins): the largest
    piece stays; each other piece (largest first) goes to the free position (>= GAPB bins from every placed piece's dilated footprint, found by FFT
    correlation of occupancy masks) that gives the smallest square-ish bounding box (min of max(W, H), then of W x H)."""
    x, y = cs.x.values, cs.y.values
    gx = np.floor((x - x.min()) / BIN).astype(int); gy = np.floor((y - y.min()) / BIN).astype(int)
    occ = np.zeros((gy.max() + 1, gx.max() + 1), bool); occ[gy, gx] = True
    lab_, npc = ndimage.label(ndimage.binary_dilation(occ, iterations=1), structure=np.ones((3, 3)))
    piece = lab_[gy, gx]; order = list(pd.Series(piece).value_counts().index)
    masks = {}
    for pc in order:
        yy, xx = np.nonzero(lab_ == pc); masks[pc] = (yy.min(), xx.min(), (lab_ == pc)[yy.min():yy.max() + 1, xx.min():xx.max() + 1])
    N = sum(max(m[2].shape) for m in masks.values()) + 20; canvas = np.zeros((2 * N, 2 * N), bool)
    pos = {}; y0, x0, m0 = masks[order[0]]; oy, ox = N - m0.shape[0] // 2, N - m0.shape[1] // 2
    canvas[oy:oy + m0.shape[0], ox:ox + m0.shape[1]] |= m0; pos[order[0]] = (oy, ox)
    bb = [oy, ox, oy + m0.shape[0], ox + m0.shape[1]]
    for pc in order[1:]:
        _, _, m = masks[pc]; h, w = m.shape
        blocked = ndimage.binary_dilation(canvas, iterations=GAPB)
        ov = fftconvolve(blocked.astype(float), m[::-1, ::-1].astype(float), mode='valid')   # ov[i, j] = overlap with m placed at (i, j)
        free = ov < 0.5; I, J = np.nonzero(free)
        W = np.maximum(bb[3], J + w) - np.minimum(bb[1], J); H = np.maximum(bb[2], I + h) - np.minimum(bb[0], I)
        k = np.lexsort((W * H, np.maximum(W, H)))[0]; oy, ox = I[k], J[k]
        canvas[oy:oy + h, ox:ox + w] |= m; pos[pc] = (oy, ox)
        bb = [min(bb[0], oy), min(bb[1], ox), max(bb[2], oy + h), max(bb[3], ox + w)]
    SH = {pc: np.array([(pos[pc][1] - masks[pc][1]) * BIN, (pos[pc][0] - masks[pc][0]) * BIN], float) for pc in order}
    sh = np.vstack([SH[p] for p in piece])
    return np.c_[x, y] + sh, piece, npc, SH


# ======================================================================================================================
# Section F6CD_dotplot_maps: Fig 6C O/E dot plot (left) and Fig 6D whole-section maps of 43563_A4 (right), one PNG + CSV
# ======================================================================================================================
if _run('F6CD_dotplot_maps'):
    STEM = f'{FIG6}/Fig6C_spatial_lr_dotplot_maps'
    MSEC = '43563_A4'; MPAIRS = ['CD80->CTLA4', 'CD80->CD28']; REC = {p: p.split('->')[1] for p in MPAIRS}
    sa, E = load_cells(['CD80', 'CD28', 'CTLA4'])
    sa = sa[sa.section == MSEC].copy(); E = E.loc[sa.index]
    print(f'{MSEC} cells with a CN: {sa.CN.notna().mean():.4f}', flush=True)

    D, piece, NPC, _ = pack(sa); EXT = (np.ptp(D[:, 0]), np.ptp(D[:, 1])); EXT0 = (np.ptp(sa.x), np.ptp(sa.y))
    MAP = {}
    for p in MPAIRS:
        l, r, P = pairs(sa, 'CD80', REC[p])
        rr = SEC[(SEC.pair == p) & (SEC.section == MSEC)].iloc[0]
        MAP[p] = dict(l=l, r=r, P=P, pl=np.unique(P[:, 0]), pr=np.unique(P[:, 1]), oe=float(rr.obs_exp_treg))
    print('maps:', {p: (len(m['P']), m['oe'], int(m['l'].sum()), len(m['pl']), int(m['r'].sum()), len(m['pr'])) for p, m in MAP.items()}, NPC, EXT, flush=True)

    # ---- figure ----
    FT = 12; plt.rcParams.update({'font.size': FT, 'font.family': 'DejaVu Sans', 'axes.linewidth': 1.1})
    C_MAC, C_REC = '#F28E2B', '#1B7837'                       # macrophage orange / receptor+ Treg green
    C_MED = {True: C_MAC, False: '#222222'}; C_SEC = {True: '#F9C08C', False: '#9A9A9A'}
    TINT = 0.58                                               # CN_k10 tab10 hue mixed 58% toward white
    WHITE = lambda c, f=TINT: tuple(np.array(matplotlib.colors.to_rgb(c)) * (1 - f) + f)
    PALE = {k_: WHITE(PAL[k_]) for k_ in PAL}; C_NOCN = '#E0E0E0'
    S_PMAC, S_PREC, S_UMAC, S_UREC, A_UNP, LW_REC, LW_MAC, LW_UREC = 14, 40, 3, 10, 0.3, 1.2, 0.5, 0.6
    TOPB, BOTB = 0.36, 0.74                                   # top band (dot-plot key), bottom band (x tick labels + 2-line x title / map key + bar)
    QH = 2.8                                                  # dot-plot axes height
    MAPH = 3.05                                               # map height
    FH = BOTB + MAPH + 0.22                                   # map labels above the maps; dot-plot key in the same top band
    YT, LPAD = 0.62, 1.95                                     # y-title band, row-label band
    ZW, ZG, LW = 0.55, 0.22, 2.4                              # '0' strip width, break gap, log-axis width
    STW, MG = 0.42, 0.18                                      # star column, gap before the maps
    QX = YT + LPAD; LX = QX + ZW + ZG
    SC = MAPH / (EXT[1] + 2 * PADUM)                          # inches per um (one scale for both maps)
    MW = (EXT[0] + 2 * PADUM) * SC; MGAP = 0.12
    MX0 = LX + LW + STW + MG; FW = MX0 + 2 * MW + MGAP + 0.04
    QB = BOTB
    fig = plt.figure(figsize=(FW, FH))
    fx, fy = lambda v: v / FW, lambda v: v / FH
    ax0 = fig.add_axes([fx(QX), fy(QB), fx(ZW), fy(QH)]); ax = fig.add_axes([fx(LX), fy(QB), fx(LW), fy(QH)], sharey=ax0)
    rng = np.random.default_rng(1); PLOT = []
    for i, pr in enumerate(ROWS):
        sg = pr in SIG; y = -i; s = SEC[SEC.pair == pr]; xs = s.obs_exp_treg.values; jy = y + rng.uniform(-0.15, 0.15, len(xs)); z = xs == 0
        ax.scatter(xs[~z], jy[~z], s=40, color=C_SEC[sg], alpha=0.95, lw=0, zorder=2)
        ax0.scatter(np.zeros(z.sum()), jy[z], s=40, color=C_SEC[sg], alpha=0.95, lw=0, zorder=2)
        PLOT += [dict(pair=pr, section=sec_, obs_exp_treg=v_, y_jittered=yj, plotted_in='0 (No Pairs) strip' if zz else 'log axis')
                 for sec_, v_, yj, zz in zip(s.section, xs, jy, z)]
        ax.plot(MED[pr], y, 'o', ms=12, mfc=C_MED[sg], mec='black' if sg else C_MED[sg], mew=1.3 if sg else 0, zorder=4)
        st = stars(POOL.loc[pr].perm_q_bh_treg)
        if st: ax.annotate(st, (1.0, y), xycoords=('axes fraction', 'data'), xytext=(6, 0), textcoords='offset points', ha='left', va='center_baseline',
                           fontsize=FT + 4, fontweight='bold', color='#222222', annotation_clip=False)
    ax.axvline(1, color='#777777', ls='--', lw=1.2, zorder=1)
    ax.set_xscale('log', base=2); ax.set_xlim(LO, HI)
    TK = [0.5, 1, 2, 4, 8]; ax.set_xticks(TK); ax.set_xticklabels([f'{t:g}' for t in TK]); ax.minorticks_off()
    ax0.set_facecolor('#EFEFEF'); ax0.set_xlim(-1, 1); ax0.set_xticks([0]); ax0.set_xticklabels(['0\n(No Pairs)'], fontsize=FT - 1, linespacing=0.95)
    ax0.set_yticks(-np.arange(len(ROWS))); ax0.set_yticklabels([lab(p) for p in ROWS], fontsize=FT)
    ax0.set_ylim(-len(ROWS) + 0.5, 0.5); ax0.spines[['top', 'right', 'left']].set_visible(False); ax0.tick_params(axis='y', length=0, pad=6)
    ax.spines[['top', 'right', 'left']].set_visible(False); ax.tick_params(axis='y', length=0, labelleft=False)
    for a_, xe in ((ax0, 1), (ax, 0)):
        a_.plot([xe], [0], transform=a_.transAxes, marker=[(-0.5, -1), (0.5, 1)], ms=14, mew=1.8, color='black', ls='', clip_on=False, zorder=10)
    fig.text(fx(LX + LW / 2), fy(0.02), 'Ligand⁺ Sender Near Receptor⁺ Treg\n(Observed / Expected, Log Scale)', ha='center', va='bottom', fontsize=FT - 1,
             linespacing=1.0)
    fig.text(fx(0.02), fy(QB + QH / 2), 'Co-Stimulatory/Suppressive\nLigand (Sender) → Treg Receptor', rotation=90,
             ha='left', va='center', fontsize=FT - 1, linespacing=1.05)
    mk = lambda c, ec, ms, mew: Line2D([], [], marker='o', ls='', mfc=c, mec=ec, ms=ms, mew=mew)
    h = [(mk(C_MED[True], 'black', 9, 1.1), mk(C_MED[False], C_MED[False], 9, 0)), (mk(C_SEC[True], C_SEC[True], 6, 0), mk(C_SEC[False], C_SEC[False], 6, 0))]
    leg = fig.legend(handles=h, labels=['Median Across Sections', 'One Section'], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
                     loc='center right', ncol=2, frameon=False, fontsize=FT - 1, bbox_to_anchor=(fx(LX + LW + STW), fy(QB + QH + 0.2)), borderaxespad=0,
                     columnspacing=1.0, handletextpad=0.3, handlelength=1.9)
    fig.canvas.draw(); _lb = leg.get_window_extent(fig.canvas.get_renderer())
    note = fig.text(_lb.x0 / fig.dpi / FW - fx(0.15), fy(QB + QH + 0.2), 'Orange: q < 0.05', ha='right', va='center', fontsize=FT - 1, color='#444444')

    # maps
    MAX = {}
    for k, p in enumerate(MPAIRS):
        m = MAP[p]; cx = MX0 + k * (MW + MGAP)
        a = fig.add_axes([fx(cx), fy(QB), fx(MW), fy(MAPH)]); MAX[p] = a
        S_BG = 0.35 * (SC / (4.8 / 10310)) ** 2
        bgc = [PALE[int(c_)] if pd.notna(c_) else matplotlib.colors.to_rgb(C_NOCN) for c_ in sa.CN.values]
        a.scatter(D[:, 0], D[:, 1], s=S_BG, c=bgc, lw=0, rasterized=True, zorder=1)
        il, ir = np.flatnonzero(m['l']), np.flatnonzero(m['r']); ul, ur = np.setdiff1d(il, m['pl']), np.setdiff1d(ir, m['pr'])
        a.scatter(D[ul, 0], D[ul, 1], s=S_UMAC, color=C_MAC, alpha=A_UNP, lw=0, zorder=2, rasterized=True)
        a.scatter(D[ur, 0], D[ur, 1], s=S_UREC, facecolors='none', edgecolors=C_REC, alpha=A_UNP, lw=LW_UREC, zorder=2, rasterized=True)
        a.scatter(D[m['pl'], 0], D[m['pl'], 1], s=S_PMAC, color=C_MAC, ec='black', lw=LW_MAC, zorder=3, rasterized=True)
        a.scatter(D[m['pr'], 0], D[m['pr'], 1], s=S_PREC, facecolors='none', edgecolors=C_REC, lw=LW_REC, zorder=4, rasterized=True)
        a.set_xlim(D[:, 0].min() - PADUM, D[:, 0].min() - PADUM + EXT[0] + 2 * PADUM); a.set_ylim(D[:, 1].max() + PADUM, D[:, 1].min() - PADUM)
        a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); [sp.set_visible(False) for sp in a.spines.values()]
        a.text(0.0, 1.0, f'CD80 → {REC[p]} Pairs', transform=a.transAxes, ha='left', va='bottom', fontsize=FT - 1)
    # one 1 mm bar + one marker key under the maps
    by0 = BOTB - 0.22; bx0 = LX + LW + 0.1
    bar = fig.add_artist(Line2D([fx(bx0), fx(bx0 + 1000 * SC)], [fy(by0), fy(by0)], color='black', lw=4, solid_capstyle='butt'))
    bart = fig.text(fx(bx0 + 1000 * SC + 0.08), fy(by0), '1 mm', ha='left', va='center', fontsize=FT - 1)
    hk = [Line2D([], [], marker='o', ls='', mfc=C_MAC, mec='black', mew=LW_MAC, ms=np.sqrt(S_PMAC)),
          Line2D([], [], marker='o', ls='', mfc='none', mec=C_REC, mew=LW_REC, ms=np.sqrt(S_PREC)),
          (Line2D([], [], marker='o', ls='', color=C_MAC, alpha=0.45, mew=0, ms=np.sqrt(S_UMAC)),
           Line2D([], [], marker='o', ls='', mfc='none', mec=C_REC, alpha=0.45, mew=LW_UREC, ms=np.sqrt(S_UREC)))]
    fig.canvas.draw(); _bt = bart.get_window_extent(fig.canvas.get_renderer())
    mleg = fig.legend(handles=hk, labels=['CD80⁺ Macrophage', 'Receptor⁺ Treg', 'Faded = Not Paired'], ncol=3,
                      handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)}, loc='center left', frameon=False, fontsize=FT - 1,
                      bbox_to_anchor=(_bt.x1 / fig.dpi / FW + fx(0.12), fy(by0)), handletextpad=0.2, columnspacing=0.6, borderaxespad=0, handlelength=1.3)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04)

    # ---- CSV ----
    sig3 = lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and pd.notna(v) else v
    A = pd.DataFrame(dict(part='dot plot, pair', plot_order=np.arange(1, 9), row_label=[lab(p) for p in ROWS], pair=ROWS,
                          sender=[POOL.loc[p].sender for p in ROWS], receivers='receptor+ strict Tregs', radius_um=RAD,
                          obs_exp_median_treg=[MED[p] for p in ROWS], perm_p_treg=POOL.loc[ROWS].perm_p_treg.values,
                          perm_q_bh_treg=POOL.loc[ROWS].perm_q_bh_treg.values, stars=[stars(POOL.loc[p].perm_q_bh_treg) for p in ROWS],
                          plotted_colour=['macrophage orange (#F28E2B; q < 0.05)' if p in SIG else 'grey' for p in ROWS],
                          sections_obs_exp_gt1_treg=POOL.loc[ROWS].sections_obs_exp_gt1_treg.astype(int).values, sections_tested_treg=11,
                          sections_zero_pairs_treg=NZERO.loc[ROWS].values, obs_total_treg=NOBS.loc[ROWS].values,
                          note='Large dot = median of the 11 per-section O/E (sections with no observed pair = 0, drawn in the 0 (No Pairs) strip); stars = '
                               'perm_q_bh_treg (BH across the 14 tested rows of spatial_lr_pairspecific.csv; ** < 0.01); row order and values as Supp 10B v13 '
                               '(spatial_lr_supp10b_simple.csv). Ligands = co-stimulatory/suppressive receptor module of the Fig 6A Treg activation programs.'))
    B = pd.DataFrame(PLOT).assign(part='dot plot, section')
    B['plot_order'] = B.pair.map({p: i + 1 for i, p in enumerate(ROWS)}); B['row_label'] = B.pair.map(lab); B = B.sort_values(['plot_order', 'section'])
    M = pd.DataFrame([dict(part='map', plot_order=k + 1, row_label=f'CD80 → {REC[p]} Pairs', pair=p, sender='Mac', receivers='receptor+ strict Tregs',
                           radius_um=RAD, section=MSEC, obs_exp_treg=MAP[p]['oe'], obs_total_treg=len(MAP[p]['P']),
                           map_n_cells_section=len(sa), map_n_cd80_pos_mac_section=int(MAP[p]['l'].sum()), map_n_cd80_pos_mac_paired=len(MAP[p]['pl']),
                           map_n_receptor_pos_treg_section=int(MAP[p]['r'].sum()), map_n_receptor_pos_treg_paired=len(MAP[p]['pr']),
                           map_n_tissue_pieces=NPC, map_extent_original_um=f'{EXT0[0]:.0f} x {EXT0[1]:.0f}', map_extent_packed_um=f'{EXT[0]:.0f} x {EXT[1]:.0f}',
                           map_scale_1mm_in=1000 * SC, map_cropped=False, map_background=f'all cells coloured by CN_k10 (tab10 mixed {TINT:.0%} toward white; cells without a CN light grey)',
                           note=(f'Whole section {MSEC} (author-selected example, Supp 10B top row), no crop; pieces packed by translation only (pack() of '
                                 f'spatial_lr_supp10b_simple.py). Paired (full opacity) = CD80+ Mac with >= 1 {REC[p]}+ strict Treg within 30 um, or the Treg with '
                                 f'>= 1 CD80+ Mac within 30 um (original coordinates). CD80+ Mac = orange dot ({S_PMAC} pt^2, black edge {LW_MAC} pt); {REC[p]}+ '
                                 f'Treg = green ring ({S_PREC} pt^2, {LW_REC} pt); unpaired cells {S_UMAC} / {S_UREC} pt^2 at {A_UNP:.0%} alpha. One scale and one 1 mm '
                                 'bar for both maps. Pair counts asserted equal to spatial_lr_pairspecific.csv and spatial_lr_supp10b_simple.csv.'))
                      for k, p in enumerate(MPAIRS)])
    C = pd.concat([A, B, M], ignore_index=True)
    lead = ['part', 'plot_order', 'row_label', 'pair', 'sender', 'receivers', 'radius_um', 'section']
    C = C[[c for c in lead if c in C] + [c for c in C if c not in lead and c != 'note'] + ['note']]
    for c in [c for c in C if c.startswith(('map_n', 'sections_')) or c in ('radius_um', 'obs_total_treg', 'plot_order')]:
        C[c] = pd.to_numeric(C[c]).round().astype('Int64')
    for c in C.select_dtypes('float').columns: C[c] = C[c].map(sig3)
    with open(f'{STEM}.csv', 'w') as fh: C.to_csv(fh, index=False)
    print(C[C.part != 'dot plot, section'].drop(columns=['note']).to_string(), flush=True)
    plt.close(fig)


# ======================================================================================================================
# Section S11BC_support_niche_maps: Supp 11B per-pair support bars and Treg CN bar; Supp 11C whole-section maps of 54033_B4
# ======================================================================================================================
if _run('S11BC_support_niche_maps'):
    OUT = S10
    STEM = f'{OUT}/SuppFig11B_spatial_lr_support_niche_maps'
    print('row order:', [(lab(p), POOL.loc[p].perm_q_bh_treg, POOL.loc[p].obs_exp_treg) for p in ROWS], flush=True)

    # ---------------- cells, CN, pair counts for all 11 sections ----------------
    sa, E = load_cells(['CD80', 'CD28', 'CTLA4', 'CD274', 'PDCD1'])
    print(f'cells with a CN: {sa.CN.notna().mean():.4f}; strict Tregs with a CN: {sa.CN[sa.treg].notna().mean():.4f}', flush=True)
    NMS = cn_names()
    SECS = sorted(sa.section.unique())
    CS, PR = {}, {}
    for s_ in SECS:
        cs = sa[sa.section == s_].copy(); CS[s_] = cs
        for pr, (lig, rec) in {'CD80->CD28': ('CD80', 'CD28'), 'CD80->CTLA4': ('CD80', 'CTLA4'), 'CD274->PDCD1': ('CD274', 'PDCD1')}.items():
            PR[(s_, pr)] = pairs(cs, lig, rec)

    # ---------------- CN counts CSV (counts only) ----------------
    cnrows = []
    CNLIST = list(range(1, 11)) + ['none']
    cnkey = lambda v: 'none' if pd.isna(v) else int(v)
    NPAIRS3 = ['CD80->CTLA4', 'CD80->CD28', 'CD274->PDCD1']
    for pr in NPAIRS3:
        acc_p, acc_t = {c: 0 for c in CNLIST}, {c: 0 for c in CNLIST}
        for s_ in SECS + ['pooled']:
            if s_ != 'pooled':
                cs = CS[s_]; P = PR[(s_, pr)][2]; cn = cs.CN.values
                cp = pd.Series([cnkey(cn[j]) for j in P[:, 1]], dtype=object).value_counts().to_dict() if len(P) else {}
                ct = pd.Series([cnkey(v) for v in cn[cs.treg.values]], dtype=object).value_counts().to_dict()
                for c in CNLIST: acc_p[c] += cp.get(c, 0); acc_t[c] += ct.get(c, 0)
            else:
                cp, ct = acc_p, acc_t
            npair, ntr = sum(cp.values()), sum(ct.values())
            for c in CNLIST:
                cnrows.append(dict(pair=pr, row_label=lab(pr) + ('⁺ Tregs'), section=s_, CN=c, CN_label=NMS.get(c, 'No CN assigned'),
                                   n_pairs_within_30um_by_treg_CN=cp.get(c, 0), pct_of_pairs=100 * cp.get(c, 0) / npair if npair else np.nan,
                                   n_strict_tregs=ct.get(c, 0), pct_of_strict_tregs=100 * ct.get(c, 0) / ntr if ntr else np.nan,
                                   n_pairs_section_total=npair, n_strict_tregs_section_total=ntr))
    CN_T = pd.DataFrame(cnrows)
    for c in ['pct_of_pairs', 'pct_of_strict_tregs']: CN_T[c] = CN_T[c].map(lambda v: float(f'{v:.3g}') if pd.notna(v) else v)
    CN_T['note'] = ('Counts only, no test. Each pair (ligand+ Mac, receptor+ strict Treg within 30 um; same pairs as obs_total_treg in '
                    'spatial_lr_pairspecific.csv) is assigned the CN_k10 of its Treg (a Treg near two ligand+ Macs counts twice, once per pair); reference = CN of all '
                    'strict Tregs in the same section(s). v9: CD80->CTLA4 added; the pooled rows are the niche bar of Supp 10B (All Tregs = pct_of_strict_tregs, '
                    'pair rows = pct_of_pairs; every CN its own segment; v10: segments >= 15% labelled, v9: >= 10%). Counts unchanged from v9 (v14: row_label uses Fig 6 sender names).')
    with open(f'{OUT}/SuppFig10B_pair_cn_counts.csv', 'w') as fh: CN_T.to_csv(fh, index=False)
    pp = CN_T[CN_T.section == 'pooled']
    for pr in NPAIRS3:
        t = pp[pp.pair == pr].sort_values('pct_of_pairs', ascending=False).head(4)
        print(pr, 'top CNs of pairs (pooled):', [(r.CN_label, r.n_pairs_within_30um_by_treg_CN, r.pct_of_pairs, r.pct_of_strict_tregs) for r in t.itertuples()], flush=True)

    # ---------------- map section ----------------
    MPAIRS = ['CD80->CTLA4', 'CD80->CD28']                  # map columns (receptor)
    REC = {p: p.split('->')[1] for p in MPAIRS}
    R1MIN = 20
    W_ = pd.concat({REC[p]: SEC[SEC.pair == p].set_index('section')[['obs_total_treg', 'obs_exp_treg']] for p in MPAIRS}, axis=1)
    W_.columns = [f'{b_}_{a_}' for a_, b_ in W_.columns]       # obs_total_treg_CTLA4, obs_exp_treg_CTLA4, ...
    W_['min_pairs'] = W_[[f'obs_total_treg_{REC[p]}' for p in MPAIRS]].min(1); W_['tot_pairs'] = W_[[f'obs_total_treg_{REC[p]}' for p in MPAIRS]].sum(1)
    W_['min_oe'] = W_[[f'obs_exp_treg_{REC[p]}' for p in MPAIRS]].min(1); W_['mean_oe'] = W_[[f'obs_exp_treg_{REC[p]}' for p in MPAIRS]].mean(1)
    print(W_.sort_values('tot_pairs', ascending=False).round(3).to_string(), flush=True)
    # highest mean O/E among sections with >= 20 pairs for both receptors (54033_B4)
    E1 = W_[W_.min_pairs >= R1MIN].sort_values(['mean_oe', 'tot_pairs'], ascending=False)
    ROW1 = E1.index[0]
    cand = lambda t: '; '.join(f"{s_} CTLA4 {r.obs_total_treg_CTLA4:.0f} pairs O/E {r.obs_exp_treg_CTLA4:.3g}, CD28 {r.obs_total_treg_CD28:.0f} pairs O/E "
                               f"{r.obs_exp_treg_CD28:.3g}, total {r.tot_pairs:.0f}" for s_, r in t.iterrows())
    ROW1_RULE = (f'54033_B4 (v9-v13 rule, kept): highest mean O/E (CTLA4, CD28) among sections with >= {R1MIN} observed pairs for both receptors. Eligible in rank order: ' + cand(E1)
                 + f'. Shown: {ROW1}.')
    print(ROW1_RULE, flush=True)
    cs = CS[ROW1]; D, piece, npc, SH = pack(cs); ext = (np.ptp(D[:, 0]), np.ptp(D[:, 1]))
    print(f'{ROW1}: {npc} pieces; extent {np.ptp(cs.x):.0f} x {np.ptp(cs.y):.0f} -> {ext[0]:.0f} x {ext[1]:.0f} um', flush=True)
    MAP = {}
    for p in MPAIRS:                                     # paired = has a partner within 30 um (original coordinates)
        l, r, P = PR[(ROW1, p)]
        MAP[p] = dict(D=D, l=l, r=r, P=P, npc=npc, ext=ext, ext0=(np.ptp(cs.x), np.ptp(cs.y)), pl=np.unique(P[:, 0]), pr=np.unique(P[:, 1]),
                      oe=W_.loc[ROW1, f'obs_exp_treg_{REC[p]}'])
    print('mapped panels:', [(p, ROW1, m['oe'], len(m['P']), len(m['pl']), len(m['pr'])) for p, m in MAP.items()], flush=True)

    # ---------------- figure ----------------
    FT = 14; plt.rcParams.update({'font.size': FT, 'font.family': 'DejaVu Sans', 'axes.linewidth': 1.2})
    C_MAC, C_REC, C_BAR = '#F28E2B', '#1B7837', '#4D4D4D'
    TINT = 0.58                                              # CN tint 58% toward white
    WHITE = lambda c, f=TINT: tuple(np.array(matplotlib.colors.to_rgb(c)) * (1 - f) + f)
    PALE = {k_: WHITE(PAL[k_]) for k_ in PAL}; C_NOCN = '#E0E0E0'
    S_PMAC, S_PREC, S_UMAC, S_UREC, A_UNP, LW_REC, LW_MAC, LW_UREC = 14, 40, 5, 18, 0.25, 1.2, 0.4, 0.9
    # layout (inches)
    LPAD, BW, BG = 2.75, 2.25, 0.45                          # row-label band, bar column width, gap between bar columns
    HDR, RWH, SKH, GAP, NRH, NXL = 0.78, 0.47, 0.42, 0.22, 0.55, 0.78   # bar headers, bar row height, significance key, gap, niche row height, niche x band
    NBW, NNW = 5.15, 1.05                                    # niche bar width, room for 'n = 6,554'
    LCW = LPAD + NBW + NNW                                   # left column width
    FIG_H = HDR + 8 * RWH + SKH + GAP + 4 * NRH + NXL
    MX, LBL, KEYH, MGAP = LCW + 0.3, 0.36, 1.72, 0.15       # maps x, corner-label band, key band under the maps, gap between maps
    MPH = FIG_H - LBL - KEYH                                 # map axes height
    m0 = MAP[MPAIRS[0]]
    SC = MPH / (m0['ext'][1] + 2 * PADUM)                    # inches per um (both maps, same section)
    MPW = (m0['ext'][0] + 2 * PADUM) * SC
    FIG_W = MX + 2 * MPW + MGAP
    S_BG = 0.35 * (SC / (4.8 / 10310)) ** 2                  # background dot area scaled with the map scale
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    fx, fy = lambda v: v / FIG_W, lambda v: v / FIG_H

    # 1. per-pair support bars (sections_obs_exp_gt1_treg / 11 and the pooled observed pairs)
    QH = 8 * RWH; QB = FIG_H - HDR - QH
    NGT1 = {p: int(POOL.loc[p].sections_obs_exp_gt1_treg) for p in ROWS}; NT = {p: int(POOL.loc[p].sections_tested_treg) for p in ROWS}
    bc = lambda p: C_MAC if p in SIG else C_BAR
    ab1 = fig.add_axes([fx(LPAD), fy(QB), fx(BW), fy(QH)]); ab2 = fig.add_axes([fx(LPAD + BW + BG), fy(QB), fx(BW), fy(QH)], sharey=ab1)
    BH = 0.62
    for i, p in enumerate(ROWS):
        ab1.barh(-i, NGT1[p], height=BH, color=bc(p), zorder=2)
        ab1.text(NGT1[p] + 0.3, -i, f'{NGT1[p]}/{NT[p]}', va='center', ha='left', fontsize=FT, zorder=3)
        ab2.barh(-i, NOBS[p], height=BH, color=bc(p), zorder=2)
        ab2.text(NOBS[p] + 8, -i, f'{NOBS[p]}', va='center', ha='left', fontsize=FT, zorder=3)
    ab1.set_xlim(0, 13.2); ab2.set_xlim(0, 375)
    ab1.set_ylim(-len(ROWS) + 0.5, 0.5)
    ab1.set_yticks(-np.arange(len(ROWS))); ab1.set_yticklabels([lab(p) for p in ROWS], fontsize=FT + 1)
    for a_, t_ in ((ab1, 'Sections With\nO/E > 1'), (ab2, 'Observed\nPairs (n)')):
        a_.set_title(t_, fontsize=FT + 1, pad=8, linespacing=1.05)
        a_.tick_params(axis='y', length=0, pad=8); a_.set_xticks([]); a_.spines[['top', 'right', 'bottom']].set_visible(False)
        a_.spines['left'].set_color('#555555')
    ab2.tick_params(axis='y', labelleft=False)
    sk = fig.legend(handles=[Patch(color=C_MAC), Patch(color=C_BAR)], labels=['Significant (q < 0.05)', 'Not Significant'], ncol=2, loc='center left',
                    frameon=False, fontsize=FT - 1, bbox_to_anchor=(fx(LPAD - 0.05), fy(QB - SKH / 2 - 0.02)), borderaxespad=0, handlelength=1.2,
                    handletextpad=0.4, columnspacing=1.2)

    # 2. niche bar (pooled CN counts of SuppFig10B_pair_cn_counts.csv; every CN its own segment; labels on segments >= 15%)
    NB_ROWS = [('All Tregs', None)] + [(f'{al(p.split("->")[0])} → {al(p.split("->")[1])} Pairs', p) for p in NPAIRS3]
    NBLAB = 15
    NBB = NXL; NBHt = 4 * NRH
    nb = fig.add_axes([fx(LPAD), fy(NBB), fx(NBW), fy(NBHt)])
    NBT, NBLBL = [], []
    for i, (lb, p) in enumerate(NB_ROWS):
        t = pp[pp.pair == (p or NPAIRS3[0])].set_index('CN')
        cnt = t.n_strict_tregs if p is None else t.n_pairs_within_30um_by_treg_CN
        cnt = cnt.loc[list(range(1, 11))].astype(float); tot = cnt.sum(); pc = 100 * cnt / tot; left = 0.0
        for k_ in range(1, 11):
            w = pc.loc[k_]
            nb.barh(-i, w, left=left, height=0.78, color=PAL[k_], ec='white', lw=0.6, zorder=2)
            if w >= NBLAB:
                lum = np.dot(matplotlib.colors.to_rgb(PAL[k_]), [0.299, 0.587, 0.114])
                nb.text(left + w / 2, -i, f'CN{k_}\n{w:.0f}%', ha='center', va='center', fontsize=FT - 1, linespacing=0.95,
                        color='white' if lum < 0.55 else 'black', zorder=3)
                NBLBL.append((lb, k_, f'{w:.0f}%'))
            left += w
        nb.text(101.5, -i, f'n = {int(tot):,}', ha='left', va='center', fontsize=FT - 1, color='#333333')
        NBT.append(dict(row=lb, pair=p or 'all strict Tregs', n=int(tot), **{f'CN{k_}_pct': pc.loc[k_] for k_ in range(1, 11)}))
    nb.set_xlim(0, 100); nb.set_ylim(-len(NB_ROWS) + 0.5, 0.5)
    nb.set_yticks(-np.arange(len(NB_ROWS))); nb.set_yticklabels([r[0] for r in NB_ROWS], fontsize=FT + 1)
    nb.set_xticks([0, 25, 50, 75, 100]); nb.set_xticklabels(['0', '25', '50', '75', '100'], fontsize=FT - 1)
    nb.set_xlabel('Tregs in Each Cellular Neighborhood (%)', fontsize=FT, labelpad=3)
    nb.tick_params(axis='y', length=0, pad=8); nb.spines[['top', 'right', 'left']].set_visible(False)
    NBT = pd.DataFrame(NBT); print(NBT.round(1).to_string(), flush=True)

    # 3. maps: 54033_B4 whole section, CD80+ Mac -> CTLA4+ Treg and -> CD28+ Treg; one scale, one 1 mm bar
    MB = KEYH
    NOCN = {}
    for c, p in enumerate(MPAIRS):
        m = MAP[p]; D = m['D']; ax_l = MX + c * (MPW + MGAP)
        a = fig.add_axes([fx(ax_l), fy(MB), fx(MPW), fy(MPH)])
        cn = CS[ROW1].CN.values; NOCN[p] = int(pd.isna(cn).sum())
        bgc = [PALE[int(c_)] if pd.notna(c_) else matplotlib.colors.to_rgb(C_NOCN) for c_ in cn]     # cells without a CN light grey
        a.scatter(D[:, 0], D[:, 1], s=S_BG, c=bgc, lw=0, rasterized=True, zorder=1)
        il, ir = np.flatnonzero(m['l']), np.flatnonzero(m['r']); ul, ur = np.setdiff1d(il, m['pl']), np.setdiff1d(ir, m['pr'])
        a.scatter(D[ul, 0], D[ul, 1], s=S_UMAC, color=C_MAC, alpha=A_UNP, lw=0, zorder=2, rasterized=True)
        a.scatter(D[ur, 0], D[ur, 1], s=S_UREC, facecolors='none', edgecolors=C_REC, alpha=A_UNP, lw=LW_UREC, zorder=2, rasterized=True)
        a.scatter(D[m['pl'], 0], D[m['pl'], 1], s=S_PMAC, color=C_MAC, ec='black', lw=LW_MAC, zorder=3, rasterized=True)
        a.scatter(D[m['pr'], 0], D[m['pr'], 1], s=S_PREC, facecolors='none', edgecolors=C_REC, lw=LW_REC, zorder=4, rasterized=True)
        a.set_xlim(D[:, 0].min() - PADUM, D[:, 0].min() - PADUM + m0['ext'][0] + 2 * PADUM); a.set_ylim(D[:, 1].max() + PADUM, D[:, 1].min() - PADUM)
        a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); [sp.set_visible(False) for sp in a.spines.values()]
        a.text(0.0, 1.0, f'CD80 → {REC[p]} Pairs', transform=a.transAxes, ha='left', va='bottom', fontsize=FT + 1)   # corner label above the map
    print('cells without a CN in', ROW1, NOCN, 'of', len(CS[ROW1]), flush=True)
    # keys under the maps: CN key (2 columns; full + tinted swatch per CN), then 1 mm bar + marker key (1 column)
    hc = [(Line2D([], [], marker='s', ls='', mfc=PAL[k_], mec='none', ms=13), Line2D([], [], marker='s', ls='', mfc=PALE[k_], mec='#888888', mew=0.4, ms=13))
          for k_ in range(1, 11)]                           # each entry: full tab10 swatch (niche bar) + 58% tint swatch (maps)
    cleg = fig.legend(handles=hc, labels=[NMS[k_] for k_ in range(1, 11)], loc='upper left', ncol=2, frameon=False, fontsize=FT - 1, handletextpad=0.3,
                      handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)}, labelspacing=0.28, borderaxespad=0, handlelength=1.9, columnspacing=1.0, bbox_to_anchor=(fx(MX), fy(KEYH - 0.12)),
                      title='Cellular Neighborhood', title_fontsize=FT, alignment='left')
    fig.canvas.draw(); _cb = cleg.get_window_extent(fig.canvas.get_renderer())
    kx = _cb.x1 / fig.dpi + 0.35                             # marker key column right of the CN key
    hk = [Line2D([], [], marker='o', ls='', mfc=C_MAC, mec='black', mew=LW_MAC, ms=np.sqrt(S_PMAC)),
          Line2D([], [], marker='o', ls='', mfc='none', mec=C_REC, mew=LW_REC, ms=np.sqrt(S_PREC)),
          (Line2D([], [], marker='o', ls='', color=C_MAC, alpha=0.45, mew=0, ms=np.sqrt(S_UMAC)),
           Line2D([], [], marker='o', ls='', mfc='none', mec=C_REC, alpha=0.45, mew=LW_UREC, ms=np.sqrt(S_UREC)))]
    mleg = fig.legend(handles=hk, labels=['CD80⁺ Macrophage', 'Receptor⁺ Treg', 'Faded = Not Paired'], ncol=1,
                      handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)}, loc='upper left', frameon=False, fontsize=FT - 1,
                      bbox_to_anchor=(fx(kx), fy(KEYH - 0.5)), handletextpad=0.4, labelspacing=0.35, borderaxespad=0, handlelength=1.5)
    by0 = KEYH - 0.26
    fig.add_artist(Line2D([fx(kx + 0.05), fx(kx + 0.05 + 1000 * SC)], [fy(by0), fy(by0)], color='black', lw=4, solid_capstyle='butt'))
    fig.text(fx(kx + 0.05 + 1000 * SC + 0.1), fy(by0), '1 mm', ha='left', va='center', fontsize=FT - 1)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.06)
    print(f'figure {FIG_W:.2f} x {FIG_H:.2f} in (before tight crop); map scale {1000 * SC:.3f} in per mm; map {MPW:.2f} x {MPH:.2f} in', flush=True)
    plt.close(fig)

    # ---------------- CSV ----------------
    ORDER_RULE = ('Rows top to bottom (same order as Fig 6): perm_q_bh_treg ascending, then sections_obs_exp_gt1_treg descending, then sections_zero_pairs_treg '
                  'ascending, then obs_total_treg descending. The three BH-significant pairs share q = 0.00932, so sections with O/E > 1 sets their order '
                  '(CD80 -> CTLA4 9/11, CD80 -> CD28 8/11, PD-L1 -> PD-1 7/11); the q = 0.874 tie (CD86 -> CTLA4 6/11, IL10 -> IL10RA 3/11) is also broken by '
                  'sections with O/E > 1; the zero-pair and observed-pair tie-breakers are not needed for any row.')
    T = POOL.loc[ROWS].reset_index()
    T = pd.DataFrame(dict(part='per-pair support bars, pooled (11 sections, 4 patients)', plot_order=np.arange(1, 9), row_label=[lab(p) for p in ROWS], pair=ROWS,
                          sender=T.sender, receivers='receptor+ strict Tregs', radius_um=RAD,
                          bar_sections_obs_exp_gt1=[f'{NGT1[p]}/{NT[p]}' for p in ROWS], bar_observed_pairs_n=NOBS.loc[ROWS].values,
                          bar_colour=['macrophage orange (#F28E2B; perm_q_bh_treg < 0.05)' if p in SIG else 'dark grey (#4D4D4D; perm_q_bh_treg >= 0.05)' for p in ROWS],
                          perm_p_treg=T.perm_p_treg, perm_q_bh_treg=T.perm_q_bh_treg, perm_q_bh_treg_8pairs=Q8.loc[ROWS].values,
                          sections_obs_exp_gt1_treg=T.sections_obs_exp_gt1_treg, sections_tested_treg=T.sections_tested_treg,
                          sections_zero_pairs_treg=NZERO.loc[ROWS].values, obs_total_treg=NOBS.loc[ROWS].values, exp_total_treg=_g.exp_total_treg.loc[ROWS].values,
                          obs_exp_median_treg=[MED[p] for p in ROWS], obs_exp_pooled_treg=T.obs_exp_treg, sign_test_p_treg=T.sign_test_p_treg,
                          row_order_rule=ORDER_RULE,
                          note='v14: the O/E dot plot (median per-section O/E, obs_exp_median_treg) is in Fig 6 (fig6_spatial_lr_combined); here only the two bar columns '
                               'are plotted. bar_sections_obs_exp_gt1 = sections with O/E > 1 (of 11); bar_observed_pairs_n = obs_total_treg = ligand+ sender - '
                               'receptor+ strict Treg pairs within 30 um summed over the 11 sections. perm_q_bh_treg = BH across the 14 tested rows of '
                               'spatial_lr_pairspecific.csv (null: ligand+ status permuted among the section sender cells, 1,000 permutations, seed 0); '
                               'orange bars = perm_q_bh_treg < 0.05. Values unchanged from v13 (nothing recomputed).'))
    S_ = SEC[SEC.pair.isin(ROWS)].copy(); S_['o'] = S_.pair.map({p: i for i, p in enumerate(ROWS)}); S_ = S_.sort_values(['o', 'section'])
    S2 = pd.DataFrame(dict(part='per section (source of the bars)', plot_order=S_.o.values + 1, row_label=[lab(p) for p in S_.pair], pair=S_.pair.values,
                           sender=S_.sender.values, receivers='receptor+ strict Tregs', radius_um=RAD, section=S_.section.values, obs_exp_treg=S_.obs_exp_treg.values,
                           obs_total_treg=S_.obs_total_treg.values, exp_total_treg=S_.exp_total_treg.values, n_receptor_pos_tregs=S_.n_receptor_pos_tregs.values,
                           counted_in_sections_obs_exp_gt1=S_.obs_exp_treg.values > 1))
    NB_ = NBT.rename(columns={'row': 'row_label'}).assign(part='niche bar, pooled (11 sections, 4 patients)', plot_order=np.arange(1, 5), radius_um=RAD,
                                                          receivers='strict Tregs (All Tregs) or the receptor+ strict Treg of each pair',
                                                          n_niche_bar=lambda t: t.n,
                                                          niche_labelled_segments=['; '.join(f'CN{k_} {v}' for r_, k_, v in NBLBL if r_ == lb) for lb, _ in NB_ROWS],
                                                          note='Percent of Tregs in each CN_k10 (pooled counts of spatial_lr_supp10b_pair_cn_counts.csv: All Tregs = '
                                                               'pct_of_strict_tregs, pair rows = pct_of_pairs; a Treg near two ligand+ Macs counts once per pair). '
                                                               'Every CN its own segment; segments >= 15% labelled. Counts, no test.').drop(columns='n')

    def mapnote(p):
        return (f'v14 map: {ROW1} whole section (no crop), the second example section (43563_A4 is in Fig 6). {ROW1_RULE} '
                f'Markers (Supp 10B v10 sizes): CD80+ Mac = filled orange dot (area {S_PMAC} pt^2, black edge {LW_MAC} pt); {REC[p]}+ strict Treg = hollow green ring '
                f'(area {S_PREC} pt^2, line {LW_REC} pt) drawn on top. Paired (full opacity) = CD80+ Mac with >= 1 {REC[p]}+ strict Treg within 30 um, or the '
                f'Treg with >= 1 CD80+ Mac within 30 um (original coordinates); unpaired cells smaller ({S_UMAC} / {S_UREC} pt^2) at {A_UNP:.0%} alpha. Pieces packed '
                f'by translation only (50 um grid components, nested by FFT mask correlation into the smallest square-ish box, >= {GAPB} bins apart; as v13). '
                f'One scale for both maps ({1000 * SC:.3g} in per mm on the figure), one 1 mm bar. Background colour = CN_k10, tab10 mixed {TINT:.0%} toward white '
                f'(as Fig 6); cells without a CN light grey ({NOCN[p]} cells).')
    M = pd.DataFrame([dict(part='map', plot_order=k + 1, row_label=f'CD80 → {REC[p]} Pairs', pair=p, sender='Mac', receivers='receptor+ strict Tregs', radius_um=RAD,
                           section=ROW1, obs_exp_treg=MAP[p]['oe'], obs_total_treg=len(MAP[p]['P']), n_receptor_pos_tregs=int(MAP[p]['r'].sum()),
                           map_corner_label=f'CD80 → {REC[p]} Pairs', map_n_cells_section=len(CS[ROW1]),
                           map_n_cd80_pos_mac_section=int(MAP[p]['l'].sum()), map_n_cd80_pos_mac_paired=len(MAP[p]['pl']),
                           map_n_receptor_pos_treg_paired=len(MAP[p]['pr']), map_selection_rule=ROW1_RULE,
                           map_cropped=False, map_scale_1mm_in=1000 * SC, map_n_tissue_pieces=MAP[p]['npc'],
                           map_extent_original_um=f'{MAP[p]["ext0"][0]:.0f} x {MAP[p]["ext0"][1]:.0f}',
                           map_extent_packed_um=f'{MAP[p]["ext"][0]:.0f} x {MAP[p]["ext"][1]:.0f}', note=mapnote(p))
                      for k, p in enumerate(MPAIRS)])
    C = pd.concat([T, S2, NB_, M], ignore_index=True)
    lead = ['part', 'plot_order', 'row_label', 'pair', 'sender', 'receivers', 'radius_um', 'section']
    C = C[lead + [c for c in C if c not in lead and c != 'note'] + ['note']]
    for c in [c for c in C if c.startswith(('map_n', 'sections_', 'n_')) or c in ('radius_um', 'obs_total_treg', 'plot_order', 'bar_observed_pairs_n')]:
        C[c] = pd.to_numeric(C[c]).round().astype('Int64')
    for c in C.select_dtypes('float').columns: C[c] = C[c].map(lambda v: float(f'{v:.3g}') if pd.notna(v) else v)
    with open(f'{STEM}.csv', 'w') as fh: C.to_csv(fh, index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    print(C[C.part != 'per section (source of the bars)'].drop(columns=['note', 'row_order_rule']).to_string(), flush=True)
