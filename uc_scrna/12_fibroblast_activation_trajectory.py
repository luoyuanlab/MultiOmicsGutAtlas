#!/usr/bin/env python
"""UC scRNA-seq atlas: fibroblast activation trajectory.

Panels: Fig 2H (DPT pseudotime UMAP), Fig 2I (PROGENy TGF-b / NF-kB / TNF-a along the trajectory),
        Fig 2J (Palantir pseudotime vs disease severity by root); supp tables for 2I, 2J and the ADAMDEC1+ fraction.
Inputs: UC atlas h5ad and subcluster4/ fibroblast UMAP + Harmony embedding (04_fibroblast_subclustering.R);
        adamdec1_sample_proportions.csv (per-sample ADAMDEC1+ fraction of connective cells).
Outputs: Fig1H-J_fibroblast_trajectory/ (pseudotime CSV, intermediate/ tables, panels, supp tables).
Run: python <file>, or RUN_SECTIONS=<id>[,<id>...] python <file>.
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

_RUN_SECTIONS = [s.strip() for s in os.environ.get('RUN_SECTIONS', '').split(',') if s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


M = '/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip'
ATLAS = f'{M}/pdated_annoed_all_cell_types.h5ad'
O = '/share/fsmresfiles/UC/revision_jenny/Fig1H-J_fibroblast_trajectory'
BLACK = '#1b1b1f'
PERICYTES = ['CD36+ Peri', 'RERGL+ Contr Peri']
PATHS = ['TGFb', 'NFkB', 'TNFa']; LAB = {'TGFb': 'TGF-β', 'NFkB': 'NF-κB', 'TNFa': 'TNF-α'}
BASE = 'ADAMDEC1+ Fib'; STATES = ['OGN+RSPO3+ Fib', 'T-interact Fib', 'Activ Fib', 'CXCL5+ Activ Fib']; ORDER = [BASE] + STATES
MIN_CELLS = 10
SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
REDS = LinearSegmentedColormap.from_list('reds_trunc', plt.cm.Reds(np.linspace(0.12, 1.0, 256)))   # light red -> dark red


def sig3(x):
    return f'{x:#.3g}'.rstrip('.')


def pfmt(p):
    if p >= 0.001:
        return sig3(p)
    m, e = f'{p:.2e}'.split('e'); return f'{m} × 10{str(int(e)).translate(SUP)}'


def ocol(f, c):
    """One obs column of an open h5ad (h5py), as strings."""
    import h5py
    g = f['obs'][c]
    if isinstance(g, h5py.Group) and 'categories' in g:
        cats = [x.decode() if isinstance(x, bytes) else x for x in g['categories'][:]]
        return pd.Categorical.from_codes(g['codes'][:], cats).astype(str)
    return np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]])


def patient_means(df, cols):
    """Per patient x cell type mean (units with >= MIN_CELLS cells)."""
    g = df.groupby(['patient', 'ct']); pm = g[cols].mean()
    return pm[g.size() >= MIN_CELLS].reset_index()


def paired_tests(pm):
    """Each trajectory state vs ADAMDEC1+ Fib in the same patients: paired Wilcoxon, BH over pathways x states."""
    rows = []
    for st in STATES:
        a = pm[pm.ct == BASE].set_index('patient'); b = pm[pm.ct == st].set_index('patient'); common = a.index.intersection(b.index)
        for p in PATHS:
            x = b.loc[common, p]; y = a.loc[common, p]
            rows.append([p, st, len(common), y.mean(), x.mean(), int((x > y).sum()), stats.wilcoxon(x, y).pvalue])
    R = pd.DataFrame(rows, columns=['pathway', 'state', 'n', 'base_mean', 'state_mean', 'n_higher', 'p'])
    R['q'] = multipletests(R['p'], method='fdr_bh')[1]
    return R


# Section F2HJ_pseudotime: DPT (ADAMDEC1+ root) and Palantir from four candidate roots, pericytes excluded (Fig 2H-J)
if _run('F2HJ_pseudotime'):
    import scanpy as sc, palantir
    scrna = sc.read_h5ad(ATLAS)
    fib = scrna[scrna.obs['cell_category'] == 'Connective'].copy()
    umap = pd.read_csv(f'{M}/subcluster4/umap_iter4_anno.csv', index_col=0)
    harm = pd.read_csv(f'{M}/subcluster4/con_harmony_iter4.csv', index_col=0)
    fib.obsm['X_umap'] = umap[['umapharmonyconiter4_1', 'umapharmonyconiter4_2']].to_numpy()
    fib.obsm['X_harmony'] = harm.to_numpy()
    fib = fib[~fib.obs['cell_type_short'].isin(PERICYTES)].copy()

    palantir.utils.run_diffusion_maps(fib, n_components=5, pca_key='X_harmony')
    palantir.utils.determine_multiscale_space(fib)
    palantir.utils.run_magic_imputation(fib)

    def centroid_cell(ad, ct):
        cells = ad.obs_names[ad.obs['cell_type_short'] == ct]; idx = ad.obs_names.get_indexer(cells)
        coords = ad.obsm['X_umap'][idx]; c = coords.mean(0)
        return cells[np.argmin(((coords - c) ** 2).sum(1))]

    pt = {}
    for name in ['ADAMDEC1+ Fib', 'Cycl Fib', 'OGN+RSPO3+ Fib', 'Activ Fib']:
        # use_early_cell_as_start=True: otherwise Palantir starts every root from the same automatically chosen cell
        pr = palantir.core.run_palantir(fib, centroid_cell(fib, name), num_waypoints=500, use_early_cell_as_start=True)
        pt[f'Palantir\n{name}'] = pd.Series(pr.pseudotime).reindex(fib.obs_names).values
    sc.pp.neighbors(fib, use_rep='X_harmony', n_neighbors=20); sc.tl.diffmap(fib)
    fib.uns['iroot'] = int(np.flatnonzero(fib.obs['cell_type_short'] == 'ADAMDEC1+ Fib')[0]); sc.tl.dpt(fib)
    pt['DPT\nADAMDEC1+ Fib'] = fib.obs['dpt_pseudotime'].values

    ptdf = pd.DataFrame(pt, index=fib.obs_names); ptdf['cell_type_short'] = fib.obs['cell_type_short'].values
    ptdf.to_csv(f'{O}/fib_pseudotime_all_methods.csv')


# Section adamdec1_fraction_mixed_model: ADAMDEC1+ fraction by condition, y ~ condition (+ C(study)) + (1|patient), REML
if _run('adamdec1_fraction_mixed_model'):
    import h5py
    import statsmodels.formula.api as smf
    CONDS = ['Healthy', 'Non-inflamed', 'Inflamed']
    with h5py.File(ATLAS, 'r') as f:
        s2p = (pd.DataFrame({'sample': ocol(f, 'sample'), 'patient': ocol(f, 'Patient'), 'study': ocol(f, 'orig.ident')})
               .drop_duplicates('sample').set_index('sample'))
    prop = pd.read_csv(f'{O}/adamdec1_sample_proportions.csv', index_col=0)
    prop = prop[prop['n_conn'] >= 20].dropna(subset=['prop_conn', 'condition']).copy(); prop['sample'] = prop.index.astype(str)
    d = prop.reset_index(drop=True).rename(columns={'prop_conn': 'y'})

    def fit(formula, d):
        last = None
        for meth in ['lbfgs', 'powell', 'nm', 'bfgs']:
            try:
                m = smf.mixedlm(formula, d, groups=d['patient']).fit(reml=True, method=meth)
                if np.all(np.isfinite(m.bse.values)): return m, meth
            except Exception as e: last = e
        raise RuntimeError(f'all optimizers failed: {last}')
    rows = []
    d = d[d['condition'].isin(CONDS)].join(s2p, on='sample')
    d['condition'] = pd.Categorical(d['condition'], CONDS)
    for model, extra in [('patient RE', ''), ('patient RE + study fixed', ' + C(study)')]:
        for ref, levs in [('Healthy', ['Non-inflamed', 'Inflamed']), ('Non-inflamed', ['Inflamed'])]:
            try: m, meth = fit(f'y ~ C(condition, Treatment("{ref}")){extra}', d)
            except RuntimeError as e: print('  FAILED', model, ref, e); continue
            for lev in levs:
                t = f'C(condition, Treatment("{ref}"))[T.{lev}]'
                rows.append(['adamdec1_fraction', model + f' [{meth}]', f'{lev} vs {ref}', m.params[t], m.bse[t], m.pvalues[t],
                             d.patient.nunique(), d.study.nunique(), len(d)])
    R = pd.DataFrame(rows, columns=['measure', 'model', 'contrast', 'estimate', 'se', 'p', 'n_patients', 'n_studies', 'n_samples'])
    R.to_csv(f'{O}/intermediate/fig1_patient_level_mixedmodel_v2.csv', index=False)


# Section F2J_disease_alignment_by_root: Spearman(Palantir pseudotime, severity) per root, patient bootstrap CI, BH (Fig 2J)
if _run('F2J_disease_alignment_by_root'):
    import h5py
    plt.rcParams.update({'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 12, 'axes.titleweight': 'normal',
        'axes.spines.top': False, 'axes.spines.right': False, 'font.family': 'sans-serif', 'figure.dpi': 300,
        'text.color': BLACK, 'axes.labelcolor': BLACK, 'xtick.color': BLACK, 'ytick.color': BLACK, 'axes.edgecolor': BLACK})
    pt = pd.read_csv(f'{O}/fib_pseudotime_all_methods.csv', index_col=0); pt.index = pt.index.astype(str)
    with h5py.File(ATLAS, 'r') as f:
        meta = pd.DataFrame({'condition': ocol(f, 'condition'), 'grp': ocol(f, 'Patient')},
                            index=ocol(f, '_index').astype(str)).reindex(pt.index)
    sev = meta['condition'].map({'Healthy': 0, 'Non-inflamed': 1, 'Inflamed': 2})
    roots = [c for c in pt.columns if c.startswith('Palantir')]
    rng = np.random.default_rng(0); NB = 2000; rows = []
    for c in roots:
        d = pd.DataFrame({'x': pt[c].values, 's': sev.values, 'g': meta['grp'].values}).dropna(); d = d[np.isfinite(d['x'])]
        est = stats.spearmanr(d['x'], d['s']).correlation
        groups = d['g'].unique(); idx = {g: np.flatnonzero(d['g'].values == g) for g in groups}
        X = d['x'].values; S = d['s'].values; bs = []
        for b in range(NB):
            j = np.concatenate([idx[g] for g in rng.choice(groups, size=len(groups), replace=True)])
            if len(np.unique(S[j])) < 2: continue
            bs.append(stats.spearmanr(X[j], S[j]).correlation)
        lo, hi = np.percentile(bs, [2.5, 97.5]); bs = np.array(bs)
        pb = max(2 * min((bs <= 0).mean(), (bs >= 0).mean()), 1 / len(bs))
        rows.append([c.split('\n')[1].replace('+ Fib', '+').replace(' Fib', '').strip(), est, lo, hi, pb, len(bs), len(groups)])
    m = pd.DataFrame(rows, columns=['root', 'rho', 'lo', 'hi', 'p_boot', 'n_boot', 'n_patients']).sort_values('rho', ascending=False).reset_index(drop=True)
    m['q'] = multipletests(m['p_boot'], method='fdr_bh')[1]
    fl = lambda v, n: f'< {1/n:.1g}' if v <= 1 / n else f'{v:#.3g}'.rstrip('.')
    m['p_txt'] = [fl(v, n) for v, n in zip(m.p_boot, m.n_boot)]; m['q_txt'] = [fl(v, n) if v > 1 / n else f'{v:#.3g}' for v, n in zip(m.q, m.n_boot)]
    m.to_csv(f'{O}/intermediate/disease_alignment_bootstrap_ci_v6.csv', index=False)
    y = np.arange(len(m))[::-1]
    col = ['#2E86AB' if 'ADAMDEC1' in r else ('#7FB069' if v >= 0 else '#C1436D') for r, v in zip(m['root'], m['rho'])]
    fig, ax = plt.subplots(figsize=(5.8, 3.2)); ax.axvline(0, color='0.6', lw=.9, zorder=0)
    xmax_txt = max(m['hi'].max(), 0) + 0.06
    for k, (yi, (_, r)) in enumerate(zip(y, m.iterrows())):
        ax.plot([r['lo'], r['hi']], [yi, yi], color=col[k], lw=2.4, solid_capstyle='round', zorder=2)
        ax.plot([0, r['rho']], [yi, yi], color=col[k], lw=1.0, alpha=.45, zorder=1)
        ax.scatter([r['rho']], [yi], color=col[k], s=60, zorder=3, edgecolor='white', linewidth=.8)
        ax.text(r['rho'], yi + 0.26, f"{r['rho']:.2f}", va='bottom', ha='center', fontsize=9.5, color=BLACK)
        qs = ('q ≤ 0.001' if r['p_boot'] <= 1 / r['n_boot'] else f"q = {r['q']:#.3g}".rstrip('.'))
        ax.text(xmax_txt, yi, qs, va='center', ha='left', fontsize=8.5, color=BLACK)
    ax.set_yticks(y); ax.set_yticklabels(m['root'], fontsize=11)
    ax.set_xlabel('Palantir Pseudotime vs Disease Severity\n(Spearman rho, 95% CI across patients)')
    ax.set_title('Disease Alignment by Root')
    pad = 0.05; ax.set_xlim(min(m['lo'].min(), 0) - pad, xmax_txt + 0.17); ax.margins(y=0.22)
    plt.tight_layout(); fig.savefig(f'{O}/fig1_disease_alignment_by_root.png', bbox_inches='tight'); plt.close(fig)


# Section F2J_root_supp_tables: supp tables of the ADAMDEC1+ fraction model (BH over 3 contrasts) and disease alignment by root
if _run('F2J_root_supp_tables'):
    r = pd.read_csv(f'{O}/intermediate/fig1_patient_level_mixedmodel_v2.csv')
    r = r[r.measure == 'adamdec1_fraction'].copy(); r['model'] = r['model'].str.replace(r' \[.*\]', '', regex=True)
    main = r[r.model == 'patient RE'].copy(); main['q'] = multipletests(main['p'], method='fdr_bh')[1]; main = main.set_index('contrast')
    adj = r[r.model == 'patient RE + study fixed'].set_index('contrast')
    CON = {'Non-inflamed vs Healthy': 'Non-inflamed vs healthy', 'Inflamed vs Healthy': 'Inflamed vs healthy',
           'Inflamed vs Non-inflamed': 'Inflamed vs non-inflamed'}
    t1 = pd.DataFrame({'Comparison': list(CON.values()),
      'Difference in ADAMDEC1⁺ fraction': [sig3(main.loc[c, 'estimate']) for c in CON],
      'SE': [sig3(main.loc[c, 'se']) for c in CON],
      'P value': [pfmt(main.loc[c, 'p']) for c in CON],
      'FDR q value': [pfmt(main.loc[c, 'q']) for c in CON],
      'P value (adjusted for study)': [pfmt(adj.loc[c, 'p']) for c in CON],
      'Patients': [int(main.loc[c, 'n_patients']) for c in CON],
      'Samples': [int(main.loc[c, 'n_samples']) for c in CON]})
    t1.to_csv(f'{O}/supp_table_adamdec1_fraction.csv', index=False)
    b = pd.read_csv(f'{O}/intermediate/disease_alignment_bootstrap_ci_v6.csv', dtype={'p_txt': str, 'q_txt': str})
    b.loc[b['p_txt'].astype(str).str.startswith('<'), 'q_txt'] = '≤ 0.001'  # p at bootstrap floor (1/2000)
    t2 = pd.DataFrame({'Pseudotime root': b['root'], 'Spearman rho': b['rho'].map(sig3),
      '95% CI lower': b['lo'].map(sig3), '95% CI upper': b['hi'].map(sig3),
      'Bootstrap P (rho ≠ 0)': b['p_txt'], 'FDR q value': b['q_txt'], 'Patients': b['n_patients']})
    t2.to_csv(f'{O}/supp_table_disease_alignment_by_root.csv', index=False)


# Section F2I_progeny_signalling_stats: PROGENy (decoupler ULM) per fibroblast; states vs ADAMDEC1+ (paired Wilcoxon, BH over 12)
if _run('F2I_progeny_signalling_stats'):
    CACHE = f'{O}/intermediate/progeny_scores_per_cell.csv.gz'   # per-cell scores, also read by F2I_signalling_heatmap
    if os.path.exists(CACHE):
        df = pd.read_csv(CACHE)
    else:
        import scanpy as sc, decoupler as dc
        scrna = sc.read_h5ad(ATLAS)
        fib = scrna[scrna.obs['cell_category'] == 'Connective'].copy()
        fib = fib[~fib.obs['cell_type_short'].isin(PERICYTES)].copy()
        raw = fib.raw.to_adata()[fib.obs_names]; raw.var_names = fib.raw.var['features'].astype(str).values
        raw = raw[:, ~raw.var_names.duplicated()].copy(); sc.pp.normalize_total(raw, target_sum=1e4); sc.pp.log1p(raw)
        dc.mt.ulm(data=raw, net=dc.op.progeny(organism='human'))
        A = dc.pp.get_obsm(raw, key='score_ulm'); df = pd.DataFrame(A.X, columns=list(A.var_names))[PATHS]
        df['ct'] = fib.obs['cell_type_short'].values; df['patient'] = fib.obs['Patient'].values
        df.to_csv(CACHE, index=False)
    R = paired_tests(patient_means(df, PATHS))
    nm = lambda c: c.replace(' Fib', '').replace('+', '⁺') + ' fibroblasts' if c != 'T-interact Fib' else 'T cell-interacting fibroblasts'
    T = pd.DataFrame({'Pathway': R['pathway'].map(LAB), 'Comparison': R['state'].map(nm) + ' vs ADAMDEC1⁺ fibroblasts',
      'Patients (paired)': R['n'], 'ADAMDEC1⁺ activity (mean)': R['base_mean'].map(sig3), 'Comparison activity (mean)': R['state_mean'].map(sig3),
      'Patients higher than ADAMDEC1⁺': R['n_higher'].astype(str) + ' of ' + R['n'].astype(str),
      'P value (paired Wilcoxon, two-sided)': R['p'].map(pfmt), 'FDR q value': R['q'].map(pfmt)})
    T.to_csv(f'{O}/supp_table_signaling.csv', index=False)


# Section F2H_pseudotime_umap: DPT pseudotime (ADAMDEC1+ root) on the iteration-4 stromal UMAP (Fig 2H)
if _run('F2H_pseudotime_umap'):
    import anndata as ad
    pt = pd.read_csv(f'{O}/fib_pseudotime_all_methods.csv', index_col=0); pt.index = pt.index.astype(str)
    # UMAP rows follow the atlas connective cells by position (as in F2HJ_pseudotime)
    o = ad.read_h5ad(ATLAS, backed='r').obs
    o = o[o['cell_category'].astype(str) == 'Connective']
    um = pd.read_csv(f'{M}/subcluster4/umap_iter4_anno.csv', index_col=0)
    um = pd.DataFrame(um[['umapharmonyconiter4_1', 'umapharmonyconiter4_2']].to_numpy(), index=o.index.astype(str), columns=['x', 'y'])
    d = pt[['DPT\nADAMDEC1+ Fib', 'cell_type_short']].join(um, how='inner')
    d.columns = ['dpt', 'ct', 'x', 'y']; d = d.sort_values('dpt')                     # high pseudotime drawn on top
    vmax = np.nanpercentile(d.dpt, 99)
    plt.rcParams.update({'font.size': 12, 'figure.dpi': 300})
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    sca = ax.scatter(d.x, d.y, c=d.dpt, cmap=REDS, vmin=0, vmax=vmax, s=1.5, lw=0, rasterized=True)
    UMAP_LAB = {'ADAMDEC1+ Fib': 'ADAMDEC1⁺ Fib', 'OGN+RSPO3+ Fib': 'OGN⁺RSPO3⁺ Fib', 'CXCL5+ Activ Fib': 'CXCL5⁺ Fib',
                'HHIP+ Myofib': 'Myofib', 'VSTM2A+ Crypt Top Fib': 'Crypt Top Fib'}
    for k, v in UMAP_LAB.items():
        s = d[d.ct == k]
        if len(s): ax.text(s.x.median(), s.y.median(), v, ha='center', va='center', fontsize=16, color=BLACK,
                           bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.6))
    ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(False) for sp in ax.spines.values()]
    cb = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02); cb.set_label('Fibroblast Pseudotime', fontsize=17); cb.ax.tick_params(labelsize=14)
    fig.tight_layout(); fig.savefig(f'{O}/fig1_fibroblast_pseudotime_umap.png', bbox_inches='tight')


# Section F2I_signalling_heatmap: pathway activity heatmap (z per pathway, stars = BH q) with a mean-pseudotime strip (Fig 2I)
if _run('F2I_signalling_heatmap'):
    import seaborn as sns
    df = pd.read_csv(f'{O}/intermediate/progeny_scores_per_cell.csv.gz')
    pt = pd.read_csv(f'{O}/fib_pseudotime_all_methods.csv', index_col=0)
    df['dpt'] = pt['DPT\nADAMDEC1+ Fib'].values                                       # same cells, same order
    pm = patient_means(df, PATHS + ['dpt'])
    R = paired_tests(pm)
    cm = pm.groupby('ct')[PATHS + ['dpt']].mean().loc[ORDER]
    z = (cm[PATHS] - cm[PATHS].mean()) / cm[PATHS].std()
    star = lambda q: '***' if q < 0.01 else '**' if q < 0.05 else '*' if q < 0.1 else ''
    ann = pd.DataFrame('', index=PATHS, columns=ORDER)
    for p in PATHS:
        ann.loc[p, BASE] = f'{cm.loc[BASE, p]:.2f}\nref'
        for st in STATES: ann.loc[p, st] = f'{cm.loc[st, p]:.2f}{star(R[(R.pathway == p) & (R.state == st)].q.iloc[0])}'
    cols = [c.replace('+ Fib', '+').replace(' Fib', '').replace('+', '⁺') for c in ORDER]
    plt.rcParams.update({'font.size': 15, 'axes.titlesize': 17, 'axes.labelsize': 16, 'figure.dpi': 300,
                         'text.color': BLACK, 'axes.labelcolor': BLACK, 'xtick.color': BLACK, 'ytick.color': BLACK})
    fig = plt.figure(figsize=(8.6, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 3], width_ratios=[30, 1], hspace=0.08, wspace=0.05)
    axp, axh, cbp, cbh = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    tc = lambda v, lo, hi: 'white' if (v - lo) / (hi - lo) > 0.6 else BLACK                               # readable text on dark cells
    P = pd.DataFrame([cm['dpt'].values], index=['Pseudotime'], columns=cols)
    sns.heatmap(P, ax=axp, cmap=REDS, vmin=0, vmax=float(P.values.max()), cbar_ax=cbp, linewidths=.5, annot=False, xticklabels=False)
    for j, v in enumerate(P.values[0]): axp.text(j + .5, .5, f'{v:.3f}', ha='center', va='center', fontsize=14, color=tc(v, 0, P.values.max()))
    axp.set_yticklabels(['Pseudotime'], rotation=0); axp.set_ylabel(''); cbp.set_ylabel(''); cbp.tick_params(labelsize=11)
    axp.set_title('Signaling Activity Along Fibroblast Activation')
    zz = z.T; zz.index = [LAB[i] for i in zz.index]; zz.columns = cols
    sns.heatmap(zz, ax=axh, cmap='Reds', cbar_ax=cbh, linewidths=.5, annot=False)
    lo, hi = float(zz.values.min()), float(zz.values.max())
    for i in range(zz.shape[0]):
        for j in range(zz.shape[1]): axh.text(j + .5, i + .5, ann.values[i, j], ha='center', va='center', fontsize=14, color=tc(zz.values[i, j], lo, hi))
    cbh.set_ylabel('Activity (z per Pathway)', fontsize=14); cbh.tick_params(labelsize=12)
    axh.set_xlabel('Fibroblast Activation Trajectory'); axh.set_ylabel('')
    axh.set_xticklabels(cols, rotation=25, ha='right', rotation_mode='anchor', fontsize=14); axh.set_yticklabels(axh.get_yticklabels(), rotation=0, fontsize=16)
    axh.text(1.0, -0.47, 'Paired Wilcoxon vs ADAMDEC1⁺ (same patients), BH-FDR: * q<0.1, ** q<0.05, *** q<0.01', transform=axh.transAxes, ha='right', va='top', fontsize=10.5, color='#777777')
    fig.savefig(f'{O}/fig1_functional_pathways_pseudotime.png', bbox_inches='tight')
