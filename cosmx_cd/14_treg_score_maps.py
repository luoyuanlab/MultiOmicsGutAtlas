#!/usr/bin/env python
"""CD CosMx: tissue maps of the Treg Suppression Score.

Panels   Fig 5A WTX biopsy maps with enlarged high and low Treg suppression regions;
         Fig 5B resection section 43563_A4 coloured by cellular neighbourhood with rule-selected high- and low-score windows
Inputs   WTX cell polygons (AtoMx raw data); treg_definition_comparison_per_cell.csv (10); cellular_neighborhoods_v2/ and
         ligand_receptor_state/treg_local_ligand_receptor.csv.gz (11)
Outputs  /share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/ (panels, window and Treg CSVs)
Run      python <file>, or RUN_SECTIONS=F5A_wtx_biopsy_maps|F5B_cn_section_map python <file>
"""
import os
import numpy as np, pandas as pd, pyarrow.parquet as pq, seaborn as sns
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from scipy import ndimage, stats
from scipy.spatial import cKDTree
_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; BR = f'{P4}/B_resection_analyses'
OUT = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches'
HC, LC = '#b2182b', '#2166ac'   # high / low Treg suppression


# ======================================================================================================================
# Section F5A_wtx_biopsy_maps: WTX biopsy maps coloured by the Treg Suppression Score, with enlarged high and low suppression regions
# ======================================================================================================================
if _run('F5A_wtx_biopsy_maps'):
    UM = 0.12028; RAW = '/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/3 Raw data'
    PIECES = {'i0476_slide1': 'Inflamed CD 2', 'i0476_slide3': 'Non-inflamed CD'}
    # windows (region, kind, x0, y0, w, h in um), each placed on a cluster of high- or low-scoring Tregs inside the tissue
    PLACED = {'i0476_slide1': [('upper left', 'act', 650, 11780, 675, 440), ('bottom left', 'low', 270, 13380, 480, 480)],
              'i0476_slide3': [('top', 'act', 1630, 540, 520, 470), ('bottom', 'low', 1740, 2560, 420, 420)]}
    CROP_TOP = {'i0476_slide3': 520}   # overview cropped above y = 520 um (sparse detached cells)
    PID = {'i0476_slide1': 'I0476\n(inflamed CD)', 'i0476_slide3': 'I0453\n(non-inflamed CD)'}
    ZT = {'act': 'High Treg Suppression Region', 'low': 'Low Treg Suppression Region'}
    TT = ['CD4 Trm', 'Gd T', 'CD8 GZMK- Trm', 'CD8 GZMK+ Trm', 'CD4 Effector', 'NK', 'Resting T', 'Exhausted CD4 Tfh', 'Tfh']
    STR = ['Activated Fib', 'Myofibroblast', 'Pericyte', 'ADAMDEC1 Fib', 'Endothelial', 'Crypt Top Fib', 'Neuron']
    GRP = {'DC': '#1a9850', 'B cell': '#7b3294', 'OGN⁺RSPO3⁺ Fib': '#f46d43', 'PLA2G2A⁺ Fib': '#d4a017',
           'Epithelium': '#a6a6a6', 'Other immune': '#ececec', 'Other stromal': '#f5efe0', 'Other': '#f5f5f5'}
    KEY = ['DC', 'B cell', 'OGN⁺RSPO3⁺ Fib', 'PLA2G2A⁺ Fib']
    CM = plt.get_cmap('RdBu_r'); NORM = Normalize(-2, 2); ZC = {'act': HC, 'low': LC}
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); PC = PC[(PC.cohort == 'biopsy') & (PC.is_core_treg.astype(str) == 'True')]
    b = pd.read_csv(f'{BR}/cellular_neighborhoods_v2/discovery_biopsy/cn_per_cell_biopsy.csv.gz', index_col=0); b = b[b.section.isin(PIECES)].copy()
    b['z'] = -PC.signature_score_D6.reindex(b.index); b['z'] = b.groupby('section').z.transform(lambda v: (v - v.mean()) / v.std()); ct = b.ct.astype(str)
    b['grp'] = np.select([ct == 'DC', ct == 'B', ct == 'OGN_RSPO3 Fib', ct == 'PLA2G2A Fib', ct == 'Epithelial', ct.isin(TT + ['Plasma', 'Mac S', 'CD9 Mac', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac',
                          'Inf Mo-Mac', 'Granulocyte', 'Neutrophil', 'Cycling']), ct.isin(STR)], ['DC', 'B cell', 'OGN⁺RSPO3⁺ Fib', 'PLA2G2A⁺ Fib', 'Epithelium', 'Other immune', 'Other stromal'], 'Other')
    plt.rcParams.update({'font.size': 17})
    fig = plt.figure(figsize=(23, 17)); gs = fig.add_gridspec(2, 3, width_ratios=[0.62, 1.2, 1.2], wspace=0.05, hspace=0.1, bottom=0.075, top=0.97, left=0.03, right=0.98); wins = []
    for r_, (p, nm) in enumerate(PIECES.items()):
        s = b[b.section == p]; T = s[s.z.notna()]; isT = s.z.notna().values
        xm, ym, yM = s.x.min(), s.y.min(), s.y.max(); ZW = PLACED[p]
        ax0 = fig.add_subplot(gs[r_, 0])
        bg = s[~s.grp.isin(KEY) & ~isT]
        for g_, sz in [('Other', 1), ('Other immune', 1.5), ('Other stromal', 1.5), ('Epithelium', 2.5)]:
            q = bg[bg.grp == g_]; ax0.scatter(q.x, q.y, s=sz, c=GRP[g_], rasterized=True, lw=0)
        for k in KEY:
            q = s[s.grp == k]; ax0.scatter(q.x, q.y, s=6, c=GRP[k], rasterized=True, lw=0)
        q = T.reindex(T.z.abs().sort_values().index); sc = ax0.scatter(q.x, q.y, s=60, c=q.z, cmap=CM, norm=NORM, edgecolor='k', lw=0.6, zorder=5)
        for lab, kind, x0, y0, ww, wh in ZW:
            ax0.add_patch(Rectangle((x0, y0), ww, wh, fill=False, ec=ZC[kind], lw=3, zorder=6))
        ax0.set_aspect('equal'); ax0.set_ylim(yM + 220, CROP_TOP.get(p, ym) - 20); ax0.axis('off'); ax0.set_anchor('NW'); ax0.set_title(f'Biopsy {PID[p]}', fontsize=18, loc='center')
        ax0.plot([xm, xm + 500], [yM + 60] * 2, 'k-', lw=3.5); ax0.text(xm + 250, yM + 190, '500 µm', ha='center', fontsize=16)
        allw = pd.Index([])
        for (_, _, x0, y0, ww, wh) in ZW: allw = allw.union(s.index[(s.x >= x0) & (s.x < x0 + ww) & (s.y >= y0) & (s.y < y0 + wh)])
        pat = p.split('_')[0].upper()   # polygon files are per run (the i0476 run holds the I0453 pieces)
        poly = pd.concat([c[c.cell.isin(allw)] for c in pd.read_csv(f'{RAW}/R6221_{pat}/wtx_{pat}-polygons.csv.gz', chunksize=2_000_000)])
        for c_, (lab, kind, x0, y0, ww, wh) in enumerate(ZW):
            x1, y1 = x0 + ww, y0 + wh; w = s[(s.x >= x0) & (s.x < x1) & (s.y >= y0) & (s.y < y1)]; wt = w[w.z.notna()]; pw = poly[poly.cell.isin(w.index)]
            ax = fig.add_subplot(gs[r_, c_ + 1]); verts, faces, edges, lws = [], [], [], []
            for nm_, g in pw.groupby('cell', sort=False):
                tr = nm_ in wt.index; fc = CM(NORM(wt.z[nm_])) if tr else GRP[s.grp.get(nm_, 'Other')]
                verts.append(np.c_[g.x_global_px.values, g.y_global_px.values] * UM); faces.append(fc); edges.append('k' if tr else 'white'); lws.append(1.6 if tr else 0.2)
            ax.add_collection(PolyCollection(verts, facecolors=faces, edgecolors=edges, linewidths=lws))
            ax.scatter(wt.x, wt.y, s=330, facecolors='none', edgecolors='k', linewidths=5, zorder=6)   # black under-ring keeps near-zero (white) scores visible
            ax.scatter(wt.x, wt.y, s=330, facecolors='none', edgecolors=[CM(NORM(v)) for v in wt.z], linewidths=3, zorder=7)
            ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_aspect('equal'); ax.set_anchor('NW'); ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_edgecolor(ZC[kind]); sp.set_linewidth(3.5)
            ax.plot([x0 + 15, x0 + 115], [y1 - 15] * 2, 'k-', lw=3.5, zorder=8); ax.text(x0 + 65, y1 - 34, '100 µm', ha='center', va='bottom', fontsize=16, zorder=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))
            ax.set_title(ZT[kind], fontsize=19, color=ZC[kind])
            wins.append(dict(piece=p, piece_label=nm, region=lab, kind=kind, x0_um=x0, y0_um=y0, width_um=ww, height_um=wh, n_cells=len(w), n_tregs=len(wt), mean_treg_suppression_score_z=wt.z.mean(),
                             n_above_mean=int((wt.z > 0).sum()), n_below_mean=int((wt.z <= 0).sum()), n_DC=int((w.ct == 'DC').sum()), n_B=int((w.ct == 'B').sum()), n_OGN_fib=int((w.ct == 'OGN_RSPO3 Fib').sum()),
                             n_PLA2G2A_fib=int((w.ct == 'PLA2G2A Fib').sum()), n_epithelial=int((w.ct == 'Epithelial').sum())))
    H = [Line2D([0], [0], marker='s', ls='', mfc=GRP[k], mec=GRP[k], ms=15, label=k) for k in GRP if k != 'Other']
    fig.legend(handles=H, loc='lower left', ncol=4, frameon=False, fontsize=17, bbox_to_anchor=(0.03, -0.005), columnspacing=1.2, handletextpad=0.3, borderaxespad=0, borderpad=0)
    cax = fig.add_axes([0.62, 0.045, 0.28, 0.014]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal')
    cb.set_label('Treg Suppression Score (z within biopsy)\n← low suppression          high suppression →', fontsize=17); cb.ax.tick_params(labelsize=15)
    fig.savefig(f'{OUT}/wtx_biopsy_treg_score_popouts.png', dpi=300, bbox_inches='tight', pad_inches=0.03); fig.savefig(f'{OUT}/wtx_biopsy_treg_score_popouts.pdf', bbox_inches='tight', pad_inches=0.03)
    TT_ = b[b.z.notna()][['section', 'x', 'y', 'z']].rename(columns={'section': 'piece', 'x': 'x_um', 'y': 'y_um', 'z': 'treg_suppression_score_z'})
    TT_['biopsy'] = TT_.piece.map({k: v.split('\n')[0] for k, v in PID.items()}); TT_['region'] = ''
    for w_ in wins:
        m_ = (TT_.piece == w_['piece']) & TT_.x_um.between(w_['x0_um'], w_['x0_um'] + w_['width_um'], inclusive='left') & TT_.y_um.between(w_['y0_um'], w_['y0_um'] + w_['height_um'], inclusive='left')
        TT_.loc[m_, 'region'] = ZT[w_['kind']]
    TT_.to_csv(f'{OUT}/wtx_biopsy_treg_score_popouts_tregs.csv', float_format='%.4g')
    pd.DataFrame(wins).to_csv(f'{OUT}/wtx_biopsy_treg_score_popouts_windows.csv', index=False, float_format='%.3g')


# ======================================================================================================================
# Section F5B_cn_section_map: section 43563_A4 coloured by cellular neighbourhood, with rule-selected high and low Treg Suppression Score windows
# Window = extent of FOV 29, slid on a 25-um grid; windows with >= 15 Tregs ranked by mean score; greedy top 3 high and bottom 3 low, each
# >= 2 window widths (x) or heights (y) from every chosen window. Shown: High 1 = H1, Low 1 = L2, Low 2 = L3, High 2 = next high window of the
# same ranking (H2, H3 skipped) with >= 15 Tregs, 95% CI excluding 0 and centre >= 2 window widths from the other three boxes.
# Display: tissue pieces moved closer together (translation only); each window drawn at its translated centre with true size.
# ======================================================================================================================
if _run('F5B_cn_section_map'):
    SEC, HIGH_FOV, BIN, STEP, GAP, MIN_T = '43563_A4', '29', 50, 25, 60, 15
    CND = f'{BR}/cellular_neighborhoods_v2'
    c = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['section', 'fov', 'x_um', 'y_um', 'cell_type', 'CN_k10', '__index_level_0__']).to_pandas()
    T = pd.read_csv(f'{BR}/ligand_receptor_state/treg_local_ligand_receptor.csv.gz', index_col=0)
    c['act_z'] = T.act_z.reindex(c.index); c['is_treg'] = c.act_z.notna()

    def mci(v):   # mean, t-based 95% CI over Tregs, one-sample t-test P
        v = np.asarray(v, float); n = len(v); mu = v.mean(); h = stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n)
        return mu, mu - h, mu + h, stats.ttest_1samp(v, 0).pvalue

    s0 = c[c.section == SEC]; s = s0.copy()
    fh = s0[s0.fov.astype(str) == HIGH_FOV]; W, H = np.ptp(fh.x_um.values), np.ptp(fh.y_um.values)
    X, Y, TR, Z, CN = s.x_um.values, s.y_um.values, s.is_treg.values, s.act_z.fillna(0).values, s.CN_k10.values
    tx, ty, tz = X[TR], Y[TR], Z[TR]; o = np.argsort(ty); tx, ty, tz = tx[o], ty[o], tz[o]
    xs = np.arange(X.min() - W / 2, X.max() - W / 2 + STEP, STEP); ys = np.arange(Y.min() - H / 2, Y.max() - H / 2 + STEP, STEP)
    win = []
    for x0 in xs:
        k = (tx >= x0) & (tx <= x0 + W)
        if k.sum() < MIN_T: continue
        yy, zz = ty[k], tz[k]; cz = np.r_[0, np.cumsum(zz)]
        a = np.searchsorted(yy, ys, 'left'); b = np.searchsorted(yy, ys + H, 'right'); n = b - a; g = n >= MIN_T
        for y0, a_, b_, n_ in zip(ys[g], a[g], b[g], n[g]): win.append((x0, y0, n_, (cz[b_] - cz[a_]) / n_))
    win = pd.DataFrame(win, columns=['x0', 'y0', 'n', 'mean'])

    def pick_n(df, chosen, n):
        out = []
        for _, r in df.iterrows():
            if all(abs(r.x0 - c[0]) >= 2 * W or abs(r.y0 - c[1]) >= 2 * H for c in chosen + out): out.append((r.x0, r.y0))
            if len(out) == n: break
        return out

    hi_ext = pick_n(win.sort_values('mean', ascending=False), [], 40); lo_ = pick_n(win.sort_values('mean'), hi_ext[:3], 3)
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

    def box_rec(lbl, x0, y0):
        w = (X >= x0) & (X <= x0 + W) & (Y >= y0) & (Y <= y0 + H); zt = s.act_z.values[w & TR]; mu, lo, hi, p = mci(zt)
        comp_ = pd.Series(CN[w]).value_counts(normalize=True).sort_index()
        return dict(box=lbl, region='High Treg Suppression' if lbl[0] == 'H' else 'Low Treg Suppression', centre_x_um=x0 + W / 2, centre_y_um=y0 + H / 2,
                    width_um=W, height_um=H, n_cells=int(w.sum()), n_tregs=len(zt), mean_act_z=mu, ci95_lo=lo, ci95_hi=hi, p_vs0=p,
                    dominant_CN=int(comp_.idxmax()), dominant_CN_name=NMS[int(comp_.idxmax())], dominant_CN_frac=comp_.max(), **{f'frac_CN{k}': comp_.get(k, 0.0) for k in range(1, 11)})

    SEL3 = pd.DataFrame([box_rec('H1', *hi_ext[0]), box_rec('L2', *lo_[1]), box_rec('L3', *lo_[2])]); fc = SEL3[['centre_x_um', 'centre_y_um']].values
    fov29 = s0[s0.fov.astype(str) == HIGH_FOV]; f29c = np.array([fov29.x_um.min() + W / 2, fov29.y_um.min() + H / 2])
    cand = []; chosen = None
    for k, (x0, y0) in enumerate(hi_ext, 1):
        if k in (1, 3): continue
        r = box_rec(f'H{k}', x0, y0); d = np.hypot(fc[:, 0] - r['centre_x_um'], fc[:, 1] - r['centre_y_um'])
        r.update({f'dist_to_{b}_um': v for b, v in zip(SEL3.box, d)}); r['dist_to_FOV29_um'] = float(np.hypot(*(f29c - [r['centre_x_um'], r['centre_y_um']])))
        r['ci_excludes_0'] = bool(r['ci95_lo'] > 0); r['far_enough'] = bool(d.min() >= 2 * W); cand.append(r)
        if r['ci_excludes_0'] and r['far_enough'] and r['n_tregs'] >= MIN_T: chosen = r; break
    pd.DataFrame(cand).to_csv(f'{OUT}/fig4_cn_section_map_boxes_high2_search.csv', index=False, float_format='%.3g')
    KEEP = {'H1': 'High 1', chosen['box']: 'High 2', 'L2': 'Low 1', 'L3': 'Low 2'}
    FB = pd.concat([SEL3, pd.DataFrame([{k: v for k, v in chosen.items() if k in SEL3.columns}])]).set_index('box').loc[list(KEEP)].reset_index()
    FB.insert(1, 'display_label', FB.box.map(KEEP)); FB = FB.rename(columns={'box': 'rank_label'})
    cxy = FB[['centre_x_um', 'centre_y_um']].values
    for j, lb in enumerate(FB.display_label): FB[f'dist_to_{lb.replace(" ", "")}_um'] = np.hypot(cxy[:, 0] - cxy[j, 0], cxy[:, 1] - cxy[j, 1])
    FB[['rank_label', 'display_label', 'region', 'centre_x_um', 'centre_y_um', 'width_um', 'height_um', 'n_cells', 'n_tregs', 'mean_act_z', 'ci95_lo', 'ci95_hi', 'p_vs0',
        'dominant_CN', 'dominant_CN_name', 'dominant_CN_frac'] + [f'frac_CN{k}' for k in range(1, 11)] + [c for c in FB.columns if c.startswith('dist_to_')]].rename(
        columns={'mean_act_z': 'mean_treg_suppression_score_z'}).to_csv(f'{OUT}/fig4_cn_section_map_boxes.csv', index=False, float_format='%.3g')

    # display transform: empty bands wider than GAP removed along x and y, then pieces translated towards the largest piece
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
    PIECE = m_.piece.values

    def to_disp(cx, cy, w):   # raw-frame window centre -> display centre (majority piece of the window's cells)
        pc = pd.Series(PIECE[w]).mode().iloc[0]; return float(fx(cx)) + SHIFT[pc][0], float(fy(cy)) + SHIFT[pc][1]

    PAL = dict(zip(range(1, 11), sns.color_palette('tab10', 10).as_hex()))
    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(figsize=(11, 11))
    for k in range(1, 11):
        q = m_[m_.CN_k10 == k]; ax.scatter(q.dx, q.dy, s=0.25, c=PAL[k], rasterized=True, lw=0, alpha=1.0)
    ax.set_aspect('equal'); ax.invert_yaxis(); ax.axis('off')
    xm, ym = m_.dx.min(), m_.dy.min(); ax.plot([xm, xm + 1000], [ym + 150] * 2, 'k-', lw=4, solid_capstyle='butt'); ax.text(xm + 500, ym + 60, '1 mm', ha='center', va='bottom', fontsize=15)
    for _, r in FB.iterrows():
        w = (X >= r.centre_x_um - W / 2) & (X <= r.centre_x_um + W / 2) & (Y >= r.centre_y_um - H / 2) & (Y <= r.centre_y_um + H / 2)
        cx, cy = to_disp(r.centre_x_um, r.centre_y_um, w); hi_b = r.rank_label.startswith('H'); col = HC if hi_b else LC
        ax.add_patch(Rectangle((cx - W / 2, cy - H / 2), W, H, fill=False, ec=col, lw=3.5, ls='-' if hi_b else (0, (3, 1.5)), zorder=5))
        # label + mean score (z) of the box's Tregs + n Tregs; above the box, except Low 1 (below: Low 2 sits directly above it)
        val = f'{r.mean_act_z:+.2f}'.replace('-', '−'); bb = dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1)
        if r.display_label == 'Low 1':
            ax.annotate(f'{r.display_label}\n{val}', (cx, cy + H / 2), xytext=(0, -8), textcoords='offset points', ha='center', va='top', fontsize=17, fontweight='bold', color=col, zorder=6, bbox=bb, linespacing=1.05)
            ax.annotate(f'n = {r.n_tregs} Tregs', (cx, cy + H / 2), xytext=(0, -45), textcoords='offset points', ha='center', va='top', fontsize=12, color=col, zorder=6, bbox=bb)
        else:
            ax.annotate(f'n = {r.n_tregs} Tregs', (cx, cy - H / 2), xytext=(0, 4), textcoords='offset points', ha='center', va='bottom', fontsize=12, color=col, zorder=6, bbox=bb)
            ax.annotate(f'{r.display_label}\n{val}', (cx, cy - H / 2), xytext=(0, 20), textcoords='offset points', ha='center', va='bottom', fontsize=17, fontweight='bold', color=col, zorder=6, bbox=bb, linespacing=1.05)
    ax.set_title(f'CD Resection {SEC}', fontsize=17)
    Hh = [Line2D([0], [0], marker='o', ls='', mfc=PAL[k], mec=PAL[k], ms=13, label=NMS[k]) for k in range(1, 11)]
    Hh += [Line2D([0], [0], marker='s', ls='', mfc='none', mec=HC, mew=3, ms=14, label='High Treg Suppression Region'),
           Line2D([0], [0], color=LC, lw=3, ls=(0, (3, 1.5)), label='Low Treg Suppression Region'),
           Line2D([0], [0], ls='', marker='', label='Box value: mean Treg Suppression\nScore (z) of the Tregs in the box')]
    fig.legend(handles=Hh, loc='center left', bbox_to_anchor=(0.88, 0.5), ncol=1, frameon=False, fontsize=15, handletextpad=0.5, labelspacing=0.6)
    fig.savefig(f'{OUT}/fig4_cn_section_map_boxes.png', dpi=300, bbox_inches='tight', pad_inches=0.03); fig.savefig(f'{OUT}/fig4_cn_section_map_boxes.pdf', bbox_inches='tight', pad_inches=0.03); plt.close(fig)
