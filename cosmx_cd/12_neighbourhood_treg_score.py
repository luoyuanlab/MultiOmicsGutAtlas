#!/usr/bin/env python
"""CD CosMx 6K resections: cellular-neighbourhood composition and the Treg Suppression Score inside vs outside each neighbourhood.

Panels made
  Fig 5C left  Treg Suppression Score inside vs outside each CN (LMM, BH over the tested CNs)
  Supp 9A, 9C  inputs (CN composition table; per-CN score table, read by 13, 15 and 17)

Inputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/ (D_shared_inputs/ from 10_treg_suppression_score.py;
    B_resection_analyses/cellular_neighborhoods_v2/ from 11_neighbourhoods_local_tables.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/A_fig4_panels/ (fig4_cn_composition_heatmap.csv/.png,
    fig4_cn_treg_state_forest.png/.csv/_legend.txt)

Run order of sections
  cn_composition_heatmap -> F5C_cn_forest

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist. Command-line arguments and environment variables used by a section are listed in its
  own docstring.
  Paths are the lab's absolute paths; edit them to match the local data location.
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
# Section cn_composition_heatmap: cell-type composition of the cellular neighbourhoods (log2 enrichment vs all cells)
# Panels: Supp 9A input (CN composition heatmap CSV)
# ======================================================================================================================
if _run('cn_composition_heatmap'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/cn_composition_heatmap.py'   # working folder of this step (relative paths resolve from here)
    """Cell-type composition of the k = 10 cellular neighbourhoods (11_neighbourhoods_local_tables.py, section cellular_neighbourhoods):
log2 enrichment of each cell type in each neighbourhood vs all cells (log2 of the type's fraction in the neighbourhood / its fraction
in all cells, pseudocount 1e-4), from cn_per_cell.parquet (CN_k10, 1-based; names from cn_treg_state_enrichment_all.csv, 0-based),
with Mo-Mac + Inf Mo-Mac, CD8 GZMK+ / GZMK- Trm and Tfh + Exhausted CD4 Tfh merged; columns grouped by lineage, all cell types shown,
compact neighbourhood labels.
Outputs panels/fig4_and_related/A_fig4_panels/fig4_cn_composition_heatmap.csv (unclipped values; read by 13 and 17) and .png/.pdf
(colour scale clipped at +/-3)."""
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.patches import Rectangle
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/A_fig4_panels'
    CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN\u207aRSPO3\u207a Fib').replace('PLA2G2A Fib', 'PLA2G2A\u207a Fib') for k, v in NM.items()}
    # compact neighbourhood labels (same abbreviations as the forest plot and module heatmap)
    TRIM = [('OGN\u207aRSPO3\u207a Fib / PLA2G2A\u207a Fib / Endothelial', 'OGN\u207aRSPO3\u207a/PLA2G2A\u207a Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    NMS = {k: short_cn(v) for k, v in NM.items()}
    # merge Mo-Mac + Inf Mo-Mac, CD8 GZMK+/- Trm, Tfh + Exhausted CD4 Tfh; enrichment from cell counts
    # (log2 of the type's fraction in the neighbourhood / its fraction in all cells)
    MERGE = {'Inf Mo-Mac': 'Mo-Mac', 'CD8 GZMK+ Trm': 'CD8 Trm', 'CD8 GZMK- Trm': 'CD8 Trm', 'Exhausted CD4 Tfh': 'Tfh'}
    cc = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['cell_type', 'CN_k10']).to_pandas()
    def enrich(lab):
        t = pd.crosstab(cc.CN_k10, lab); fr = t.div(t.sum(1), axis=0); return np.log2((fr + 1e-4) / (t.sum(0) / t.values.sum() + 1e-4))   # pseudocount as in cellular_neighbourhoods
    E = enrich(cc.cell_type.replace(MERGE)); E.index = [NMS[i] for i in E.index]   # compact labels, as in the forest and module plots
    LIN = {'Epithelial': ['Epithelial'],
           'Stromal': ['OGN_RSPO3 Fib', 'PLA2G2A Fib', 'ADAMDEC1 Fib', 'Activated Fib', 'Crypt Top Fib', 'Myofibroblast', 'Pericyte', 'Endothelial', 'Neuron'],
           'Myeloid': ['DC', 'CD9 Mac', 'Mac S', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac', 'Neutrophil', 'Granulocyte'],
           'Lymphoid': ['Treg', 'CD4 Trm', 'CD4 Effector', 'Tfh', 'Resting T', 'CD8 Trm', 'Gd T', 'NK', 'B', 'Plasma'],
           'Other': ['Cycling']}
    cols = [c for v in LIN.values() for c in v if c in E.columns]; assert set(cols) == set(E.columns), set(E.columns) - set(cols)
    H = E[cols].clip(-3, 3); H.columns = [c.replace('OGN_RSPO3 Fib', 'OGN\u207aRSPO3\u207a Fib').replace('PLA2G2A Fib', 'PLA2G2A\u207a Fib') for c in cols]
    E[cols].to_csv(f'{OUT}/fig4_cn_composition_heatmap.csv', float_format='%.3g')
    PAL = dict(zip(range(1, 11), sns.color_palette('tab10', 10).as_hex()))
    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(figsize=(12.6, 4.6))
    sns.heatmap(H, cmap='RdBu_r', center=0, vmin=-3, vmax=3, ax=ax, linewidths=0.4, linecolor='white',
                cbar_kws={'label': 'log2 enrichment\nvs all cells', 'shrink': 0.9, 'pad': 0.012})
    ax.set_xticklabels(H.columns, rotation=60, ha='right', fontsize=16); ax.set_yticklabels(H.index, rotation=0, fontsize=16)
    for i, t in enumerate(ax.get_yticklabels()): t.set_color('k')
    for i in range(10): ax.add_patch(Rectangle((-0.6, i + 0.08), 0.45, 0.84, color=PAL[i + 1], clip_on=False, transform=ax.transData))
    ax.tick_params(axis='y', pad=22)
    x = 0
    for lin, cs in LIN.items():
        n = len([c for c in cs if c in E.columns]); ax.axvline(x, color='k', lw=1.5) if x else None
        ax.text(x + n / 2, -0.25, lin, ha='center', va='bottom', fontsize=17, fontweight='bold'); x += n
    ax.set_xlabel(''); ax.set_ylabel('Spatial Cellular\nNeighbourhood', fontsize=17)
    fig.savefig(f'{OUT}/fig4_cn_composition_heatmap.png', dpi=300, bbox_inches='tight', pad_inches=0.03); fig.savefig(f'{OUT}/fig4_cn_composition_heatmap.pdf', bbox_inches='tight', pad_inches=0.03); plt.close(fig)


# ======================================================================================================================
# Section F5C_cn_forest: Treg Suppression Score inside vs outside each cellular neighbourhood (forest plot)
# Panels: Fig 5C left, Supp 9C input (Treg Suppression Score inside vs outside each CN)
# ======================================================================================================================
if _run('F5C_cn_forest'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/F5C_cn_forest.py'   # working folder of this step (relative paths resolve from here)
    """Neighbourhood forest plot, CD resections (all Tregs, statistics on the continuous Treg Suppression Score).
Neighbourhoods: cellular neighbourhoods, k = 10 (11_neighbourhoods_local_tables.py, section cellular_neighbourhoods;
cellular_neighborhoods_v2/cn_per_cell.parquet CN_k10, 1-based; names from cn_treg_state_enrichment_all.csv, whose CN column is 0-based).
Treg Suppression Score (code: activation_z) = -signature_score_D6 z-scored within section (D_shared_inputs/treg_definition_comparison_per_cell.csv,
all resection Tregs).
Per neighbourhood: linear mixed model activation_z ~ in_CN + (1 | section); x = coefficient (SD, Tregs inside minus outside the neighbourhood) with
95% CI; BH across the tested neighbourhoods; sections agreeing with the pooled direction / sections with >= 3 Tregs inside the neighbourhood
(CSV and legend).
Testing rule (pre-specified; same few-Treg exclusion as the other Treg-level Fig 5 panels): a neighbourhood is tested only if it has
>= 100 Tregs in total AND >= 10 Tregs in at least 5 sections; other neighbourhoods are dropped from the plot and the BH adjustment, but keep their
numbers in the CSV (tested = False, mixed_q empty).
Sensitivity (CSV only): CMH OR activated vs low-suppressive for the median split and the thirds split.
Figure ~7.3 in wide (x range fitted to the CIs and stars; CN labels 14 pt, axis labels 16 pt); layout asserts keep the stars clear of the
CIs and the zero line and the text free of overlaps.
Output panels/fig4_and_related/A_fig4_panels/fig4_cn_treg_state_forest.png/.csv/_legend.txt.
Run: OPENBLAS_CORETYPE=SkylakeX /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.ticker import MultipleLocator
    from statsmodels.stats.contingency_tables import StratifiedTable
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/A_fig4_panels'; CNV = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    MIN_TREG_TOTAL = 100; MIN_TREG_SEC = 10; MIN_SEC = 5   # pre-stated testing rule
    def stars(q): return '' if not np.isfinite(q) or q >= 0.05 else ('*' if q >= 0.01 else ('**' if q >= 0.001 else '***'))
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); pr = PC[PC.cohort == 'resection']
    cn = pq.read_table(f'{CNV}/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas()
    NM = pd.read_csv(f'{CNV}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name
    T = pr[['section', 'D6_median', 'D6_tertile']].assign(CN=cn.loc[pr.index, 'CN_k10'].values, act=-pr.signature_score_D6.values)
    T['activation_z'] = T.groupby('section').act.transform(lambda v: (v - v.mean()) / v.std())
    def cmh(flag, col):
        tabs = []
        for s, g in T.groupby('section'):
            xx = flag[g.index]; ge, gl = g[col] == 'Early', g[col] == 'Late'
            tabs.append(np.array([[xx[ge].sum(), (~xx[ge]).sum()], [xx[gl].sum(), (~xx[gl]).sum()]], float) + 0.5)
        return StratifiedTable(tabs).oddsratio_pooled
    rows = []
    for c in sorted(T.CN.unique()):
        d = T.assign(inCN=(T.CN == c).astype(float)); m = smf.mixedlm('activation_z ~ inCN', d, groups=d.section).fit(reml=True, method='lbfgs'); b = m.params['inCN']; k = n = 0
        for _, g in d.groupby('section'):
            if (g.inCN == 1).sum() >= 3: n += 1; k += int(np.sign(g[g.inCN == 1].activation_z.mean() - g[g.inCN == 0].activation_z.mean()) == np.sign(b))
        rows.append(dict(CN=int(c), CN_name=NM[c - 1], n_tregs_in_CN=int((T.CN == c).sum()), activation_z_in_minus_out=b, ci_lo=b - 1.96 * m.bse['inCN'], ci_hi=b + 1.96 * m.bse['inCN'],
                         mixed_p=m.pvalues['inCN'], sections_agreeing=k, sections=n, n_sections_ge10_tregs=int((d[d.inCN == 1].groupby('section').size() >= MIN_TREG_SEC).sum()),
                         pct_activated_in_CN=100 * (T[T.D6_median == 'Early'].CN == c).mean(), pct_low_suppressive_in_CN=100 * (T[T.D6_median == 'Late'].CN == c).mean(),
                         sens_cmh_or_activated_vs_low_median=cmh(T.CN == c, 'D6_median'), sens_cmh_or_activated_vs_low_thirds=cmh(T.CN == c, 'D6_tertile')))
    r = pd.DataFrame(rows); r['tested'] = (r.n_tregs_in_CN >= MIN_TREG_TOTAL) & (r.n_sections_ge10_tregs >= MIN_SEC)
    r['mixed_q'] = np.nan; r.loc[r.tested, 'mixed_q'] = multipletests(r.loc[r.tested, 'mixed_p'], method='fdr_bh')[1]   # BH over tested CNs only
    FAIL = r[~r.tested]; print('Not tested (few Tregs):', FAIL[['CN', 'CN_name', 'n_tregs_in_CN', 'n_sections_ge10_tregs']].to_string(index=False))
    csv_txt = r.rename(columns={'activation_z_in_minus_out': 'treg_suppression_score_z_in_minus_out'}).to_csv(index=False, float_format='%.3g')   # reader-facing column name
    open(f'{OUT}/fig4_cn_treg_state_forest.csv', 'w').write(csv_txt)
    # compact neighbourhood labels (same abbreviations as the module heatmap); full names stay in the CSV
    TRIM = [('OGN_RSPO3 Fib / PLA2G2A Fib / Endothelial', 'OGN\u207aRSPO3\u207a/PLA2G2A\u207a Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    r['name'] = [short_cn(x) for x in r.CN_name]
    rp = r[r.tested].reset_index(drop=True)   # plotted rows: tested neighbourhoods only
    # ---- plot (~7.3 in wide; x range fitted to the CIs and stars, which stay beside the CIs)
    FS_CN, FS_TICK, FS_LAB, FS_STAR = 14, 14, 16, 18
    plt.rcParams.update({'font.size': FS_LAB})
    FW, FH = 7.3, 4.1 - 0.29 * (len(r) - len(rp)); y = np.arange(len(rp)) * 0.8   # constant per-row height
    AXW = 3.25
    fig = plt.figure(figsize=(FW, FH)); ax = fig.add_axes([(FW - AXW - 0.35) / FW, 0.11, AXW / FW, 0.77])
    col = np.where(rp.mixed_q >= 0.05, '#9e9e9e', np.where(rp.activation_z_in_minus_out > 0, '#c26a77', '#5b7fa6'))
    ax.hlines(y, rp.ci_lo, rp.ci_hi, color='0.35', lw=2.0, zorder=2); ax.scatter(rp.activation_z_in_minus_out, y, c=col, s=150, zorder=3, edgecolor='k', lw=0.7)
    ST = [ax.text(hi + 0.02 if o >= 0 else lo - 0.02, yi, stars(q), ha='left' if o >= 0 else 'right', va='center', fontsize=FS_STAR)
          for yi, lo, hi, o, q in zip(y, rp.ci_lo, rp.ci_hi, rp.activation_z_in_minus_out, rp.mixed_q)]
    ax.axvline(0, color='0.4', lw=1, ls='--'); ax.set_yticks(y); ax.set_yticklabels(rp.name, fontsize=FS_CN); ax.tick_params(axis='y', pad=2)
    ax.set_ylabel('Spatial Neighborhood', fontsize=FS_LAB, labelpad=4); ax.invert_yaxis(); ax.set_ylim(y[-1] + 0.45, y[0] - 0.45)
    ax.set_xlabel('Treg Suppression Score (z)\nIn vs Outside Neighborhood\n\u2190 Low-suppressive   Activated \u2192', fontsize=FS_LAB)
    ax.tick_params(axis='x', labelsize=FS_TICK); ax.xaxis.set_major_locator(MultipleLocator(0.2))
    sns.despine(ax=ax)
    x0, x1 = rp.ci_lo.min() - 0.02, rp.ci_hi.max() + 0.02                 # x range: the CIs plus the stars beside them, 3 pt margin
    for it in range(50):
        ax.set_xlim(x0, x1); fig.canvas.draw(); R_ = fig.canvas.get_renderer(); inv = ax.transData.inverted(); m = 3 * fig.dpi / 72
        bx = [t.get_window_extent(R_) for t in ST if t.get_text()]
        n0 = min([rp.ci_lo.min() - 0.02] + [inv.transform((b.x0 - m, 0))[0] for b in bx]); n1 = max([rp.ci_hi.max() + 0.02] + [inv.transform((b.x1 + m, 0))[0] for b in bx])
        if abs(n0 - x0) < 1e-6 and abs(n1 - x1) < 1e-6: break
        x0, x1 = n0, n1
    assert it < 49; ax.set_xlim(x0, x1); fig.canvas.draw(); R_ = fig.canvas.get_renderer()
    # ---- layout checks: stars inside the axes and clear of dots / CIs / zero line; no text overlaps; nothing bold; saved width 7-7.5 in
    AB = ax.get_window_extent(R_); px = fig.dpi / 72
    ov = lambda a, b: a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1
    for t, yi, lo, hi, o in zip(ST, y, rp.ci_lo, rp.ci_hi, rp.activation_z_in_minus_out):
        if not t.get_text(): continue
        b = t.get_window_extent(R_); assert b.x0 >= AB.x0 and b.x1 <= AB.x1, ('star outside axes', t.get_text())
        (cl, _), (ch, _) = ax.transData.transform((lo, yi)), ax.transData.transform((hi, yi)); z0 = ax.transData.transform((0, yi))[0]
        assert (b.x0 > ch + 1.0 * px) if o >= 0 else (b.x1 < cl - 1.0 * px), 'star on its CI'
        assert not (b.x0 <= z0 <= b.x1), 'star on the zero line'
    TX = [ax.xaxis.label, ax.yaxis.label] + list(ax.get_yticklabels()) + [t for t in ax.get_xticklabels() if t.get_text()] + [t for t in ST if t.get_text()]
    BX = [t.get_window_extent(R_) for t in TX]
    for i in range(len(BX)):
        for j in range(i + 1, len(BX)): assert not ov(BX[i], BX[j]), ('text overlap', TX[i].get_text(), TX[j].get_text())
    assert all(t.get_fontweight() in ('normal', 400) for t in TX)
    EXTRA = [ax.xaxis.label]   # matplotlib >= 3.6 collapses a long x label in the axes' tight bbox; include it explicitly so it is never cropped
    tb = fig.get_tightbbox(R_, bbox_extra_artists=EXTRA); W_IN, H_IN = tb.width + 2 * 0.03, tb.height + 2 * 0.03
    xb = ax.xaxis.label.get_window_extent(R_); assert tb.x0 <= xb.x0 / fig.dpi + 1e-6 and tb.x1 >= xb.x1 / fig.dpi - 1e-6, 'x label outside the saved area'
    assert 7.0 <= W_IN <= 7.5, W_IN
    print(f'x range {x0:.3f} to {x1:.3f}; figure {W_IN:.2f} x {H_IN:.2f} in', flush=True)
    fig.savefig(f'{OUT}/fig4_cn_treg_state_forest.png', dpi=300, bbox_inches='tight', pad_inches=0.03, bbox_extra_artists=EXTRA); plt.close(fig)   # PNG only
    # ---- legend: section-level consistency of the significant neighbourhoods read from the CSV (sections with >= 3 Tregs in the neighbourhood)
    RC = pd.read_csv(f'{OUT}/fig4_cn_treg_state_forest.csv')
    SG = RC[RC.tested & (RC.mixed_q < 0.05)].sort_values('CN')
    assert len(SG) and (SG.sections_agreeing <= SG.sections).all()
    parts, last = [], None
    for x in SG.itertuples():
        d_ = 'higher' if x.treg_suppression_score_z_in_minus_out > 0 else 'lower'
        parts.append((f'{d_} in ' if d_ != last else '') + f'{int(x.sections_agreeing)}/{int(x.sections)}' + (' sections' if not parts else '') + f' for CN{int(x.CN)}'); last = d_
    cons = parts[0] if len(parts) == 1 else ', '.join(parts[:-1]) + ' and ' + parts[-1]
    FL = RC[~RC.tested].sort_values('CN')
    excl = '; '.join(f'CN{int(x.CN)}, {int(x.n_tregs_in_CN)} Tregs, {int(x.n_sections_ge10_tregs)} section{"s" if x.n_sections_ge10_tregs != 1 else ""} with \u226510 Tregs' for x in FL.itertuples())
    open(f'{OUT}/fig4_cn_treg_state_forest_legend.txt', 'w').write(
        f'Treg Suppression Score across spatial cellular neighbourhoods in CD resections (CosMx 6K; 11 sections from 4 patients; all {len(T):,} Tregs). '
    'Neighbourhoods: k-means (k = 10) on 20-nearest-neighbour cell-type composition, named by their top enriched cell types. Dot, difference in the '
    'continuous Treg Suppression Score (scRNA-seq-derived late-versus-early pseudotime signature, sign reversed, z-scored within section) between Tregs '
    'inside and outside the neighbourhood, from a linear mixed model with section random intercept; bar, 95% CI; colour, direction when significant '
    '(grey otherwise); stars, Benjamini-Hochberg-adjusted across the tested neighbourhoods (* q < 0.05, ** q < 0.01, *** q < 0.001). '
    f'Section-level consistency of the significant neighbourhoods: the score inside the neighbourhood was {cons} (sections with \u22653 Tregs in the neighbourhood). '
    f'Neighbourhoods with fewer than 100 Tregs or fewer than 5 sections with \u226510 Tregs ({excl}) were not tested. Neighbourhood names are abbreviated; full names and the '
    'split-based sensitivity odds ratios (median and thirds) are in the source data.\n')
    print(r[['name', 'tested', 'n_tregs_in_CN', 'n_sections_ge10_tregs', 'activation_z_in_minus_out', 'mixed_q', 'sections_agreeing', 'sections', 'sens_cmh_or_activated_vs_low_median', 'sens_cmh_or_activated_vs_low_thirds']].round(4).to_string(index=False))
