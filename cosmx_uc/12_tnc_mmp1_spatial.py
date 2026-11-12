#!/usr/bin/env python
"""UC CosMx (WTX, 6K): TNC in OGN+RSPO3+ fibroblasts near MMP1+ macrophages, and TNC vs its upstream regulators.

Panels made
  Fig 3F   TNC in OGN+RSPO3+ fibroblasts and MMP1 in Mac S+M+S+ (WTX maps; inset, logistic regression with sample fixed effects)
  Supp 6E  TNC vs its top NicheNet upstream regulators within OGN+RSPO3+ fibroblasts, against 300 random genes

Inputs
  /share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad (09_mac_fib_spatial_input.py)
  /share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/whole_trans/Processed_merged/v7/h5ad_obj/qc_wtx_v7_merged_anno.h5ad
  /share/fsmresfiles/UC/revision_jenny/Fig2E_I_J_ligand_regulator_cooccurrence_neighbourhood/intermediate/fig2c_per_piece.csv (11_ligand_regulator_cooccurrence.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/Fig2F_tnc_mmp1_spatial/ (panel, pooled and per-piece test CSVs)
  /share/fsmresfiles/UC/revision_jenny/SuppFig6E_tnc_regulators/ (panel and CSVs)

Run order of sections
  F3F_tnc_near_mmp1_test -> F3F_tnc_mmp1_map -> S6E_tnc_regulator_correlations -> S6E_tnc_regulator_plot

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
# Section F3F_tnc_near_mmp1_test: TNC detection in OGN+RSPO3+ fibroblasts near vs far from MMP1+ SELENOP+MMP9+SPP1+ macrophages
# Panels: Fig 3F (inset statistics)
# ======================================================================================================================
if _run('F3F_tnc_near_mmp1_test'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig2F_tnc_mmp1_spatial/01_tnc_near_mmp1_test.py'   # working folder of this step (relative paths resolve from here)
    """Compares TNC detection in OGN+RSPO3+ fibroblasts near vs far from MMP1+ SELENOP+MMP9+SPP1+ macrophages.
Inflamed UC pieces from the Fig 3E cohort (CosMx WTX slide x patient; CosMx 6K, one section per block). Per piece with >=1 MMP1+
macrophage: OGN+RSPO3+ fibroblasts split by distance to the nearest MMP1+ macrophage (<= 50 um 'near' vs > 50 um 'far');
TNC detection rate and mean TNC per group. Pieces with >= 3 near fibroblasts are compared pairwise (near vs far within piece);
pooled: logistic regression of TNC detection on near/far with piece fixed effects.
Outputs: tnc_near_mmp1_by_piece.csv, tnc_near_mmp1_pooled.csv (inset of the Fig 3F map)"""
    import numpy as np, pandas as pd, anndata as ad, statsmodels.formula.api as smf
    from scipy.spatial import cKDTree
    from scipy.stats import wilcoxon
    O = '/share/fsmresfiles/UC/revision_jenny/Fig2F_tnc_mmp1_spatial'; PX = 0.12028; MAX_UM = 50; MIN_NEAR = 3
    def two(a):
        g = a[:, [list(a.var_names).index(x) for x in ('TNC', 'MMP1')]].to_memory(); return g, (g.X.toarray() if hasattr(g.X, 'toarray') else np.asarray(g.X))
    cells = []
    def add(cohort, piece, patient, xy, ct, X):
        fib = np.where(ct == 'OGN+RSPO3+ Fib')[0]; mac = np.where((ct == 'Mac S+M+S+') & (X[:, 1] > 0))[0]
        if len(mac) == 0 or len(fib) == 0: return
        d = cKDTree(xy[mac]).query(xy[fib])[0]
        cells.append(pd.DataFrame(dict(cohort=cohort, piece=piece, patient=patient, near=(d <= MAX_UM).astype(int), dist_um=d,
                                       tnc=X[fib, 0], tnc_pos=(X[fib, 0] > 0).astype(int))))
    W = ad.read_h5ad('/share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/whole_trans/Processed_merged/v7/h5ad_obj/qc_wtx_v7_merged_anno.h5ad', backed='r')
    g, X = two(W); o = g.obs; xy = np.asarray(g.obsm['spatial_fov']) * PX
    for (s, p), _ in o.groupby([o.slide.astype(str), o.patient.astype(str)], observed=True):
        if p == 'I0346': continue                                                       # healthy control
        i = np.where((o.slide.astype(str) == s).values & (o.patient.astype(str) == p).values)[0]
        add('WTX', f'{s}|{p}', p, xy[i], o.aucell_cell_type_short.astype(str).values[i], X[i])
    S = ad.read_h5ad('/share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad', backed='r')
    g, X = two(S); o = g.obs
    keep = o.cohort.astype(str).str.contains('Serial').values & (o.condition.astype(str) == 'UC').values
    for b in o[keep].block.astype(str).unique():
        i = np.where(keep & (o.block.astype(str) == b).values)[0]
        add('6K', b, o.patient.astype(str).values[i][0], o[['x_um', 'y_um']].values[i], o.cell_type_short.astype(str).values[i], X[i])
    D = pd.concat(cells, ignore_index=True)
    P = D.groupby(['cohort', 'piece', 'patient', 'near']).agg(n=('tnc', 'size'), tnc_pos_pct=('tnc_pos', lambda v: 100 * v.mean()), tnc_mean=('tnc', 'mean')).unstack('near')
    P.columns = [f'{a}_{"near" if b else "far"}' for a, b in P.columns]; P = P.reset_index()
    ok = P[P.n_near >= MIN_NEAR].dropna(subset=['tnc_pos_pct_far'])
    print(P.round(2).to_string()); print(f'\npieces with >= {MIN_NEAR} near fibroblasts: {len(ok)} from {ok.patient.nunique()} patients')
    print('TNC+ % near > far in', int((ok.tnc_pos_pct_near > ok.tnc_pos_pct_far).sum()), 'of', len(ok),
          '| Wilcoxon P =', round(wilcoxon(ok.tnc_pos_pct_near, ok.tnc_pos_pct_far).pvalue, 4) if len(ok) >= 5 else 'n<5')
    Dk = D[D.piece.isin(P.piece)]
    m = smf.logit('tnc_pos ~ near + C(piece)', Dk).fit(disp=0)
    print(f'pooled (piece fixed effects): OR near vs far = {np.exp(m.params["near"]):.2f} (95% CI {np.exp(m.conf_int().loc["near", 0]):.2f}-{np.exp(m.conf_int().loc["near", 1]):.2f}), P = {m.pvalues["near"]:.3g}; n fibroblasts = {len(Dk)}, near = {Dk.near.sum()}')
    P.round(3).to_csv(f'{O}/tnc_near_mmp1_by_piece.csv', index=False)
    pd.DataFrame([dict(OR=np.exp(m.params['near']), CI_lo=np.exp(m.conf_int().loc['near', 0]), CI_hi=np.exp(m.conf_int().loc['near', 1]), P=m.pvalues['near'], n_fib=len(Dk), n_near=int(Dk.near.sum()), n_pieces=len(ok), n_patients=ok.patient.nunique(), n_pieces_up=int((ok.tnc_pos_pct_near > ok.tnc_pos_pct_far).sum()))]).to_csv(f'{O}/tnc_near_mmp1_pooled.csv', index=False)


# ======================================================================================================================
# Section F3F_tnc_mmp1_map: TNC / MMP1 co-localization maps of both CosMx WTX slides with a shared colour scale and boxed near pairs
# Panels: Fig 3F
# ======================================================================================================================
if _run('F3F_tnc_mmp1_map'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig2F_tnc_mmp1_spatial/02_tnc_mmp1_spatial_combined.py'   # working folder of this step (relative paths resolve from here)
    """TNC / MMP1 co-localization panel: both CosMx WTX slides (i0261_i0331_slide4 = I0261, i0261_i0331_slide1 = I0331) side by
side with ONE shared colour scale per gene (5th-99th percentile of the target cells pooled across both slides) and data-defined
boxes. Measured expression; expressing cells enlarged and drawn highest on top; no smoothing.
Boxes: every TNC+ OGN+RSPO3+ fibroblast whose nearest MMP1+ SELENOP+MMP9+SPP1+ macrophage lies within MAX_UM (global coordinates,
0.12028 um/px), boxed around the pair with PAD_UM padding; overlapping boxes merged. Scale bar instead of pixel axes.
Inset: pooled OR from tnc_near_mmp1_pooled.csv (section F3F_tnc_near_mmp1_test).
Output: fig2_tnc_mmp1_combined.png/.pdf, tnc_mmp1_boxed_pairs.csv"""
    import numpy as np, pandas as pd, anndata as ad, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    import matplotlib.cm as cm, matplotlib.colors as mcolors
    from matplotlib.patches import Rectangle, Patch
    from scipy.spatial import cKDTree
    O = '/share/fsmresfiles/UC/revision_jenny/Fig2F_tnc_mmp1_spatial'
    PX, MAX_UM, PAD_UM, BAR_UM = 0.12028, 50, 60, 500
    S_BG, S_ZERO, S_EXPR, A_EXPR = 4, 4, 26, 1.0     # point sizes (background, non-expressing target, expressing) and alpha
    SLIDES = ['i0261_i0331_slide4', 'i0261_i0331_slide1']; BOX = '#2E8B57'   # patient I0261, patient I0331 (one patient per slide)
    trunc = lambda c, lo=0.12, hi=1.0: mcolors.LinearSegmentedColormap.from_list('t', c(np.linspace(lo, hi, 256)))
    CM = {'TNC': trunc(cm.Reds, 0.15), 'MMP1': trunc(cm.Blues, 0.15)}   # full light-to-dark range; thin outline keeps pale dots visible on grey
    TGT = {'TNC': 'OGN+RSPO3+ Fib', 'MMP1': 'Mac S+M+S+'}
    A = ad.read_h5ad('/share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/whole_trans/Processed_merged/v7/h5ad_obj/qc_wtx_v7_merged_anno.h5ad', backed='r')
    G = A[:, [list(A.var_names).index(g) for g in ('TNC', 'MMP1')]].to_memory(); G = G[G.obs.slide.isin(SLIDES).values]
    X = G.X.toarray() if hasattr(G.X, 'toarray') else np.asarray(G.X); E = {'TNC': X[:, 0], 'MMP1': X[:, 1]}
    ct = G.obs.aucell_cell_type_short.astype(str).values; sl = G.obs.slide.astype(str).values; C = np.asarray(G.obsm['spatial_fov'])
    MASK = {'TNC': ct == TGT['TNC'], 'MMP1': ct == TGT['MMP1']}
    LIM = {g: (np.quantile(E[g][MASK[g] & (E[g] > 0)], 0.05), np.quantile(E[g][MASK[g]], 0.99)) for g in E}   # shared across slides; colour scale over detected cells
    def boxes(idx):
        f = idx[MASK['TNC'][idx] & (E['TNC'][idx] > 0)]; m = idx[MASK['MMP1'][idx] & (E['MMP1'][idx] > 0)]
        d, j = cKDTree(C[m]).query(C[f]); keep = d * PX <= MAX_UM; pad = PAD_UM / PX; bx = []; rows = []
        for fi, mi, dd in zip(f[keep], m[j[keep]], d[keep]):
            p = np.vstack([C[fi], C[mi]]); bx.append([*(p.min(0) - pad), *(p.max(0) + pad)])
            rows.append(dict(slide=sl[fi], fib_TNC=E['TNC'][fi], mac_MMP1=E['MMP1'][mi], distance_um=dd * PX))
        merged = True
        while merged and len(bx) > 1:                                   # merge overlapping boxes
            merged = False
            for a in range(len(bx)):
                for b in range(a + 1, len(bx)):
                    A_, B_ = bx[a], bx[b]
                    if not (A_[2] < B_[0] or B_[2] < A_[0] or A_[3] < B_[1] or B_[3] < A_[1]):
                        bx[a] = [min(A_[0], B_[0]), min(A_[1], B_[1]), max(A_[2], B_[2]), max(A_[3], B_[3])]; bx.pop(b); merged = True; break
                if merged: break
        return bx, rows
    plt.rcParams.update({'font.size': 13, 'figure.dpi': 300})
    fig = plt.figure(figsize=(14, 13)); gs = fig.add_gridspec(4, 4, width_ratios=[12, 12, 0.5, 3.6], height_ratios=[1, 1, 0.8, 0.2], wspace=0.45, hspace=0.3)

    allrows = []; sc = {}
    for k, s in enumerate(SLIDES):
        ax = fig.add_subplot(gs[:, k]); idx = np.where(sl == s)[0]
        ax.scatter(C[idx, 0], C[idx, 1], c='#e3e3e3', s=S_BG, linewidths=0, rasterized=True)
        for g in ('TNC', 'MMP1'):
            t = idx[MASK[g][idx]]; z = t[E[g][t] <= 0]; e = t[E[g][t] > 0]; e = e[np.argsort(E[g][e])]
            ax.scatter(C[z, 0], C[z, 1], c='#e3e3e3', s=S_ZERO, linewidths=0, rasterized=True)   # not detected: drawn as background
            sc[g] = ax.scatter(C[e, 0], C[e, 1], c=E[g][e], s=S_EXPR, alpha=A_EXPR, linewidths=0.4, edgecolors='#555555', cmap=CM[g], vmin=LIM[g][0], vmax=LIM[g][1], zorder=3)
        bx, rows = boxes(idx); allrows += rows
        for x0, y0, x1, y1 in bx: ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor=BOX, lw=2.2, zorder=4))
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); [v.set_visible(False) for v in ax.spines.values()]
        x0, y0 = C[idx, 0].min(), C[idx, 1].min(); L = BAR_UM / PX
        ax.plot([x0, x0 + L], [y0 - 500, y0 - 500], color='black', lw=3); ax.text(x0 + L / 2, y0 - 750, f'{BAR_UM} µm', ha='center', va='top', fontsize=12)
        ax.set_title(f"{ {'i0261_i0331_slide4': 'I0261', 'i0261_i0331_slide1': 'I0331'}[s] }", fontsize=14)
    cb1 = fig.colorbar(sc['TNC'], cax=fig.add_subplot(gs[0, 2])); cb1.set_label('TNC\n(OGN⁺RSPO3⁺ Fib)', fontsize=12, labelpad=8)
    cb2 = fig.colorbar(sc['MMP1'], cax=fig.add_subplot(gs[1, 2])); cb2.set_label('MMP1\n(SELENOP⁺MMP9⁺SPP1⁺ Mac)', fontsize=12, labelpad=8)
    lax = fig.add_subplot(gs[2, 2:]); lax.axis('off')
    # pooled within-sample test (section F3F_tnc_near_mmp1_test): logistic regression with sample fixed effects
    ST = pd.read_csv(f'{O}/tnc_near_mmp1_pooled.csv').iloc[0]
    SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    pfmt = lambda p: f'{p:.3f}' if p >= 0.001 else (lambda m, e: f'{float(m):.1f} × 10{str(int(e)).translate(SUP)}')(*f'{p:.1e}'.split('e'))
    fig.text(0.47, 0.84, f'TNC⁺ odds near (≤ {MAX_UM} µm)\nvs far from MMP1⁺ Mac:\nOR {ST.OR:.1f} (95% CI {ST.CI_lo:.2f}–{ST.CI_hi:.2f})\nP = {pfmt(ST.P)}', ha='center', va='top', fontsize=12)
    lax.legend(handles=[Patch(facecolor='none', edgecolor=BOX, lw=2.2, label=f'TNC⁺ OGN⁺RSPO3⁺ Fib\nwith an MMP1⁺\nSELENOP⁺MMP9⁺SPP1⁺ Mac\nwithin {MAX_UM} µm')],
               loc='upper left', frameon=False, fontsize=11.5, handlelength=1.4)
    fig.savefig(f'{O}/fig2_tnc_mmp1_combined.png', bbox_inches='tight', pad_inches=0.5); fig.savefig(f'{O}/fig2_tnc_mmp1_combined.pdf', bbox_inches='tight', pad_inches=0.5)
    R = pd.DataFrame(allrows).round(3); R.to_csv(f'{O}/tnc_mmp1_boxed_pairs.csv', index=False); print(R.to_string()); print('limits', LIM)


# ======================================================================================================================
# Section S6E_tnc_regulator_correlations: within-fibroblast correlation of TNC with its NicheNet regulators and MMP1-route genes, against random genes
# Panels: Supp Fig 6E (correlations)
# ======================================================================================================================
if _run('S6E_tnc_regulator_correlations'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig6E_tnc_regulators/01_pathway_spatial.py'   # working folder of this step (relative paths resolve from here)
    """Correlations along the MMP1 -> receptor (F2R/PAR1, CD44) -> signalling -> TNC route and between TNC and its top-7 NicheNet
regulators in OGN+RSPO3+ fibroblasts. Same cohort, pairing and pooling as Fig 3E (11_ligand_regulator_cooccurrence.py):
inflamed UC CosMx WTX (I0261, I0331; per slide) + CosMx 6K (4 patients; per piece, middle section). Pieces with >= 10
OGN+RSPO3+ fibroblasts and >= 5 Mac S+M+S+.
  (a) MMP1 in the nearest Mac S+M+S+ vs each route gene in the fibroblast;
  (b) within the fibroblast, each route gene / regulator vs TNC.
Per-piece Spearman r pooled across pieces by random-effects meta-analysis (Fisher z, var 1.06/(n-3)). Within-fibroblast
correlations are ranked against TNC's correlation with 300 random genes detected in >= 5% of OGN+RSPO3+ fibroblasts
(background_percentile).
Outputs: pathway_spatial_mmp1_tnc.csv, pathway_spatial_tnc_background.csv."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc
    from scipy import stats
    from scipy.spatial import cKDTree
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig6E_tnc_regulators'; PX=0.12028
    def lognorm(a):
        a=a.copy(); sc.pp.normalize_total(a); sc.pp.log1p(a); return a
    S=sc.read_h5ad('/share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad')
    S=lognorm(S[(S.obs['cohort'].astype(str).str.contains('Serial'))&(S.obs['condition'].astype(str)=='UC')].copy())
    S.obs['ct']=S.obs['cell_type_short'].astype(str); S.obs['piece']=S.obs['block'].astype(str); S.obsm['xy']=S.obs[['x_um','y_um']].values; S.obs['cohort2']='Serial 6K'
    W=sc.read_h5ad('/share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/whole_trans/Processed_merged/v7/h5ad_obj/qc_wtx_v7_merged_anno.h5ad')
    W=W[W.obs['patient'].astype(str).isin(['I0261','I0331'])].copy(); W.obs['ct']=W.obs['aucell_cell_type_short'].astype(str)
    W.obs['piece']=W.obs['slide'].astype(str); W.obsm['xy']=W.obsm['spatial_fov']*PX; W.obs['cohort2']='WTX'
    def col(a,mask,g):
        x=a[mask,g].X; return x.toarray().ravel() if hasattr(x,'toarray') else np.asarray(x).ravel()
    rows=[]; BG=[]
    REGS=['SMAD3','JUN','STAT3','SMAD4','FLI1','SMAD2','RELA']   # TNC's top-7 regulators in the NicheNet gene-regulatory layer (by weight)
    NODES=sorted(set(['F2R','CD44','MAPK3','MAPK1','JUN','FOS','ETS1','SMAD3','RELA','TNC']+REGS))
    RNG=np.random.default_rng(0)
    for a in (W,S):
        for pc in a.obs['piece'].unique():
            pm=(a.obs['piece']==pc).values; om=pm&(a.obs['ct']=='OGN+RSPO3+ Fib').values; mm=pm&(a.obs['ct']=='Mac S+M+S+').values
            if om.sum()<10 or mm.sum()<5: continue
            _,nn=cKDTree(a.obsm['xy'][mm]).query(a.obsm['xy'][om],k=1)
            d=pd.DataFrame({'piece':f"{a.obs['cohort2'].iat[0]}|{pc}",'MMP1':col(a,mm,'MMP1')[nn]})
            for g in NODES: d[g]=col(a,om,g) if g in a.var_names else np.nan
            rows.append(d)
            fx=a[om].X; fx=fx.toarray() if hasattr(fx,'toarray') else np.asarray(fx)
            BG.append((f"{a.obs['cohort2'].iat[0]}|{pc}",pd.DataFrame(fx,columns=a.var_names)))
    D=pd.concat(rows,ignore_index=True)
    def meta(rs,ns):
        z=np.arctanh(np.clip(rs,-0.999,0.999)); v=1.06/(np.array(ns)-3); w=1/v; zf=(w*z).sum()/w.sum(); Q=(w*(z-zf)**2).sum(); k=len(z)
        t2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())) if k>1 else 0; ws=1/(v+t2); zr=(ws*z).sum()/ws.sum(); se=np.sqrt(1/ws.sum())
        return np.tanh(zr),np.tanh(zr-1.96*se),np.tanh(zr+1.96*se),2*stats.norm.sf(abs(zr/se)),k
    def pooled(x,y):
        rs,ns=[],[]
        for pc,g in D.groupby('piece'):
            g=g[[x,y]].dropna()
            if len(g)>=10 and g[x].nunique()>1 and g[y].nunique()>1: rs.append(stats.spearmanr(g[x],g[y])[0]); ns.append(len(g))
        return meta(np.array(rs),ns) if rs else (np.nan,)*4+(0,)
    out=[]
    for g in NODES:   # (a) MMP1 in the nearest macrophage vs each pathway node in the fibroblast
        r,lo,hi,p,k=pooled('MMP1',g); out.append({'test':'MMP1 (nearest Mac S+M+S+) vs fibroblast gene','gene':g,'r':r,'r_lo':lo,'r_hi':hi,'P':p,'pieces':k})
    for g in sorted(set(['F2R','CD44','MAPK3','MAPK1','JUN','FOS','ETS1','SMAD3','RELA']+REGS)):   # (b) within the fibroblast, node vs TNC
        r,lo,hi,p,k=pooled(g,'TNC'); out.append({'test':'within fibroblast: gene vs TNC','gene':g,'r':r,'r_lo':lo,'r_hi':hi,'P':p,'pieces':k})
    # specificity: TNC vs 300 random genes present in both panels and detected in >= 5% of OGN+RSPO3+ fibroblasts (pooled the same way)
    common=set.intersection(*[set(df.columns) for _,df in BG]); det=pd.concat([df[list(common)]>0 for _,df in BG]).mean()
    pool=[g for g in det[det>=0.05].index if g not in set(NODES)]; bg_genes=list(RNG.choice(pool,size=min(300,len(pool)),replace=False))
    bg=[]
    for g in bg_genes:
        rs,ns=[],[]
        for pc,df in BG:
            if len(df)>=10 and df[g].nunique()>1 and df['TNC'].nunique()>1: rs.append(stats.spearmanr(df[g],df['TNC'])[0]); ns.append(len(df))
        if rs: bg.append({'gene':g,'r':meta(np.array(rs),ns)[0]})
    BGD=pd.DataFrame(bg); BGD.to_csv(f'{O}/pathway_spatial_tnc_background.csv',index=False)
    for o in out:
        if o['test'].startswith('within'): o['background_percentile']=100*(BGD.r<o['r']).mean()
    print('background genes',len(BGD),'median r',round(BGD.r.median(),4),'95th pct',round(BGD.r.quantile(0.95),4))
    T=pd.DataFrame(out); T.to_csv(f'{O}/pathway_spatial_mmp1_tnc.csv',index=False); print(T.round(4).to_string(index=False)); print('SAVED')


# ======================================================================================================================
# Section S6E_tnc_regulator_plot: bar plot of TNC vs its upstream regulators within OGN+RSPO3+ fibroblasts (inflamed UC, CosMx WTX + 6K)
# Panels: Supp Fig 6E
# ======================================================================================================================
if _run('S6E_tnc_regulator_plot'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig6E_tnc_regulators/02_fig_pathway.py'   # working folder of this step (relative paths resolve from here)
    """TNC vs its upstream regulators within OGN+RSPO3+ fibroblasts (inflamed UC, CosMx WTX + 6K).
Regulators = TNC's top-7 NicheNet gene-regulatory regulators. Numbers above bars = percentile of each regulator's r among TNC's r with 300
random detected genes (set-level P = 0.0001). Bars: pooled Spearman r (random effects across samples), 95% CI;
colors P < 0.05 / < 0.1 / >= 0.1. MMP1 vs fibroblast CD44 / F2R / TNC is reported in pathway_spatial_mmp1_tnc.csv only.
Inputs from section S6E_tnc_regulator_correlations. Output: fig_mmp1_tnc_pathway.png."""
    import numpy as np, pandas as pd, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig6E_tnc_regulators'
    T=pd.read_csv(f'{O}/pathway_spatial_mmp1_tnc.csv'); BG=pd.read_csv(f'{O}/pathway_spatial_tnc_background.csv')
    PP=pd.read_csv('/share/fsmresfiles/UC/revision_jenny/Fig2E_I_J_ligand_regulator_cooccurrence_neighbourhood/intermediate/fig2c_per_piece.csv').drop_duplicates(['cohort','slide'])
    npat=len({p for c,sl in zip(PP.cohort,PP.slide) for p in (['I0261','I0331'] if c=='WTX' else [sl.split('_')[0]])})   # WTX slides hold I0261 + I0331
    col=lambda p: '#D62828' if p<0.05 else '#F77F00' if p<0.1 else '#FCBF49'
    plt.rcParams.update({'font.size':11,'figure.dpi':300,'axes.spines.top':False,'axes.spines.right':False})
    B=T[T.test.str.startswith('within')&T.gene.isin(['SMAD3','JUN','STAT3','SMAD4','FLI1','SMAD2','RELA'])].sort_values('r',ascending=False).reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(6.4,4.2)); x=np.arange(len(B)); lo,hi=BG.r.quantile([0.025,0.975])
    ax.bar(x,B.r,color=[col(p) for p in B.P],edgecolor='black',linewidth=1,zorder=2)
    ax.errorbar(x,B.r,yerr=np.vstack([B.r-B.r_lo,B.r_hi-B.r]),fmt='none',ecolor='black',capsize=3,zorder=3)
    for i,(h,pc) in enumerate(zip(B.r_hi,B.background_percentile)): ax.text(i,max(h,0)+0.006,f'{pc:.0f}th',ha='center',va='bottom',fontsize=9)
    ax.axhline(0,color='black',lw=0.8); ax.set_xticks(x); ax.set_xticklabels(B.gene,fontsize=10.5); ax.set_ylim(top=float(B.r_hi.max())*1.2)
    ax.set_xlabel('TNC upstream regulator (NicheNet)'); ax.set_ylabel('Spearman r with TNC\nin OGN⁺RSPO3⁺ fibroblasts')
    ax.set_title('TNC vs Its Upstream Regulators in OGN⁺RSPO3⁺ Fibroblasts\nInflamed UC Patients',fontsize=11.5)
    ax.text(0.99,0.98,'Numbers: percentile vs TNC correlation\nwith 300 random genes',transform=ax.transAxes,ha='right',va='top',fontsize=8.5,color='#777777')
    ax.legend(handles=[Patch(facecolor='#D62828',edgecolor='black',label='P < 0.05'),Patch(facecolor='#F77F00',edgecolor='black',label='P < 0.1'),
                       Patch(facecolor='#FCBF49',edgecolor='black',label='P ≥ 0.1')],
              loc='upper center',bbox_to_anchor=(0.5,-0.2),ncol=3,frameon=False,fontsize=9.5)
    plt.tight_layout(); fig.savefig(f'{O}/fig_mmp1_tnc_pathway.png',bbox_inches='tight'); plt.close(fig); print('SAVED')
