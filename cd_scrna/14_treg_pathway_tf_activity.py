#!/usr/bin/env python
"""CD scRNA-seq Tregs: pathway and transcription-factor activity along pseudotime.

Panels made
  Fig 4G   Treg activation and suppressive programs, MYC/metabolism and hypoxia/Wnt/inflammation programs along pseudotime
           (last vs first quintile, patient-paired Wilcoxon, BH)
  Supp 8F  CollecTRI TF activity along pseudotime (top falling and rising TFs, sign test, driver targets)

Inputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/treg_pseudotime.h5ad (12_treg_pseudotime.py)
  decoupler resources (MSigDB Hallmark, PROGENy, CollecTRI) and the local snapshots in
    /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/ and /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/archive/csv/
  utils/http_user_agent.py (helper module in this folder)

Outputs
  /share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/ (Fig 4G panel and CSVs)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/extended_data_fig8/ext8_tf_activity_along_pseudotime.csv (also read by treatment_validation/03)
  /share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness/ (Supp 8F panel)

Run order of sections
  F4G_pathways -> S8F_tf_correlation -> S8F_tf_heatmap

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist. Command-line arguments and environment variables used by a section are listed in its own
  docstring.
  Paths are the lab's absolute paths (to be replaced by placeholders before public release).
  Environment: see the README section "Environment" (OPENBLAS_CORETYPE=SkylakeX on Intel Sapphire Rapids nodes; the
  environment asserts of each section are kept as comments).
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
# Section F4G_pathways: Pathway activity per pseudotime quintile (patient-paired) and the Fig 4G panel.
# Panels: Fig 4G
# Note: imports cd_scrna/utils/http_user_agent.py (sends a browser User-Agent with the decoupler resource downloads)
# ======================================================================================================================
if _run('F4G_pathways'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/F4G_pathways.py'   # working folder of this step (relative paths resolve from here)
    """Pathway activity per pseudotime quintile, allowing non-monotonic trends (dc.mt.ulm on the 'log1p' layer of the
Treg pseudotime object). Per pathway:
  - per-patient mean activity in each pseudotime bin (patients with >= MIN_CELLS Tregs in that bin);
  - mixed model on the patient-bin means: score ~ linear + quadratic orthogonal-polynomial bin terms + (1 | patient)
    (statsmodels MixedLM, ML fit). Overall bin effect = likelihood-ratio test vs score ~ 1 + (1 | patient) (2 df);
    linear term = monotonic trend, quadratic term = rise-then-fall (negative) or fall-then-rise (positive);
  - each bin vs bin 0: patient-paired Wilcoxon signed-rank;
  - BH FDR within each resource and test family.
Resources: PROGENy, MSigDB hallmark (all pathways), a suppression panel specified before testing (SUPP_SETS below:
Treg-intrinsic suppression / Treg identity programs from MSigDB) and an inflammation / Treg-instability panel
(INFL_SETS) - all of them are reported. The suppression panel is also scored with AUCell (dc.mt.aucell, per-cell gene
ranks) and the same statistics are reported.
Fig 4G: three groups of pathways (Treg activation and suppressive function, MYC / metabolism, hypoxia / Wnt /
inflammation), mean activity in 50 equal-size windows of rank-scaled pseudotime; last vs first quintile, patient-paired
Wilcoxon, BH.
Outputs: Fig3G_pathway_stats_all_gene_sets.csv, Fig3G_pathway_patient_bin_means.csv, Fig3G_pathways_along_pseudotime.png
and its _stats.csv / _values.csv.
Run: OPENBLAS_CORETYPE=SkylakeX RUN_SECTIONS=F4G_pathways python <this file>
"""
    import os, sys, warnings
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX', 'run with OPENBLAS_CORETYPE=SkylakeX'
    warnings.filterwarnings('ignore')
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
    sys.path.insert(0, _os.path.join(_THIS_DIR, 'utils'))   # http_user_agent.py ships in cd_scrna/utils/
    import numpy as np, pandas as pd, anndata as ad, decoupler as dc
    import statsmodels.formula.api as smf
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    OUT = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals'
    MIN_CELLS = 3
    BINS = [0, 1, 2, 3, 4]
    plt.rcParams.update({'font.size': 11})

    SUPP_SETS = {
        'HALLMARK_IL2_STAT5_SIGNALING': 'IL-2/STAT5 signaling (Hallmark)',
        'HALLMARK_TGF_BETA_SIGNALING': 'TGF-β signaling (Hallmark)',
        'REACTOME_CO_INHIBITION_BY_CTLA4': 'Co-inhibition by CTLA4 (Reactome)',
        'REACTOME_RUNX1_AND_FOXP3_CONTROL_THE_DEVELOPMENT_OF_REGULATORY_T_LYMPHOCYTES_TREGS': 'RUNX1/FOXP3 Treg program (Reactome)',
        'GOBP_NEGATIVE_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION': 'Neg. regulation of αβ T cell activation (GO)',
        'GOBP_NEGATIVE_REGULATION_OF_T_CELL_PROLIFERATION': 'Neg. regulation of T cell proliferation (GO)',
        'GOBP_T_CELL_TOLERANCE_INDUCTION': 'T cell tolerance induction (GO)',
        'GSE15659_TREG_VS_TCONV_UP': 'Human Treg vs Tconv up (GSE15659)',
        'GSE15659_NONSUPPRESSIVE_TCELL_VS_ACTIVATED_TREG_DN': 'Suppressive activated Treg vs non-suppressive T (GSE15659)',
        'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY': 'TNFR2 non-canonical NF-κB (Reactome)',
    }

    # Inflammatory / Treg-instability programs (all reported)
    INFL_SETS = {
        'HALLMARK_INFLAMMATORY_RESPONSE': 'Inflammatory response (Hallmark)',
        'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNF-α signaling via NF-κB (Hallmark)',
        'HALLMARK_INTERFERON_GAMMA_RESPONSE': 'IFN-γ response (Hallmark)',
        'HALLMARK_IL6_JAK_STAT3_SIGNALING': 'IL-6/JAK/STAT3 signaling (Hallmark)',
        'GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE': 'Pos. regulation of inflammatory response (GO)',
        'GOBP_T_HELPER_17_TYPE_IMMUNE_RESPONSE': 'Th17 immune response (GO)',
        'GOBP_T_HELPER_1_TYPE_IMMUNE_RESPONSE': 'Th1 immune response (GO)',
        'GOBP_T_CELL_MEDIATED_CYTOTOXICITY': 'T cell-mediated cytotoxicity (GO)',
        'GSE15659_NONSUPPRESSIVE_TCELL_VS_ACTIVATED_TREG_UP': 'Non-suppressive T vs activated Treg up (GSE15659)',
        'GSE15659_TREG_VS_TCONV_DN': 'Tconv vs Treg up (GSE15659)',
        'TH17_LIKE_FIG4D': 'Th17-like module (Fig. 4D genes)',
    }
    TH17_GENES = ['RORC', 'CCR6', 'IL17A', 'IL17F']


    def get_net(name):
        try:
            import http_user_agent  # noqa: F401
            return getattr(dc.op, name)(organism='human')
        except Exception as e:
            print(f'{name}: download failed ({e}); using cached snapshot')
            return pd.read_csv(f'{D}/archive/csv/{name}_human_snapshot.csv.gz')


    def get_msigdb_sets(sets):
        f = f'{D}/archive/csv/msigdb_human_snapshot.csv.gz'
        if os.path.exists(f):
            m = pd.read_csv(f, usecols=['genesymbol', 'geneset'])
        else:
            import http_user_agent  # noqa: F401
            m = dc.op.resource('MSigDB', organism='human')
        m = m[m['geneset'].isin(sets)].drop_duplicates(['geneset', 'genesymbol'])
        return m.rename(columns={'geneset': 'source', 'genesymbol': 'target'})[['source', 'target']]


    def stars(q):
        return '' if not np.isfinite(q) or q >= 0.05 else ('*' if q >= 0.01 else ('**' if q >= 0.001 else '***'))


    clus0 = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    o = clus0.obs
    o['pt_bin'] = o['pt_bin'].astype(int)
    supp_net = get_msigdb_sets(list(SUPP_SETS))
    assert set(supp_net['source']) == set(SUPP_SETS), set(SUPP_SETS) - set(supp_net['source'])
    scores = {}
    for res, net in [('PROGENy', get_net('progeny')), ('MSigDB_hallmark', get_net('hallmark')), ('Suppression_panel', supp_net)]:
        dc.mt.ulm(data=clus0, net=net, layer='log1p')
        scores[res] = dc.pp.get_obsm(clus0, key='score_ulm')
    infl_net = pd.concat([get_msigdb_sets([k for k in INFL_SETS if k != 'TH17_LIKE_FIG4D']),
                          pd.DataFrame({'source': 'TH17_LIKE_FIG4D', 'target': TH17_GENES})], ignore_index=True)
    dc.mt.ulm(data=clus0, net=infl_net, layer='log1p', tmin=4)   # tmin=4 keeps the 4-gene Th17-like module
    scores['Inflammation_panel'] = dc.pp.get_obsm(clus0, key='score_ulm')
    dc.mt.aucell(data=clus0, net=supp_net, layer='log1p')
    scores['Suppression_panel_AUCell'] = dc.pp.get_obsm(clus0, key='score_aucell')
    for k, v in scores.items():
        print(k, v.shape)

    # orthogonal polynomial contrasts for 5 equally spaced bins
    P = np.polynomial.polynomial.polyvander(np.array(BINS, float), 2)
    Q, _ = np.linalg.qr(P)
    POLY = pd.DataFrame({'lin': Q[:, 1] * np.sign(Q[-1, 1]), 'quad': Q[:, 2] * np.sign(Q[0, 2])}, index=BINS)


    def binwise_stats(scad, resource):
        X = pd.DataFrame(np.asarray(scad.X), index=scad.obs_names, columns=scad.var_names)
        rows, means = [], []
        for v in X.columns:
            d = pd.DataFrame({'patient': o['patient'].astype(str).values, 'bin': o['pt_bin'].values, 'x': X[v].values})
            g = d.groupby(['patient', 'bin'])['x'].agg(['mean', 'size']).reset_index()
            g = g[g['size'] >= MIN_CELLS].rename(columns={'mean': 'y'})
            g = g.join(POLY, on='bin')
            full = smf.mixedlm('y ~ lin + quad', g, groups=g['patient']).fit(reml=False)
            null = smf.mixedlm('y ~ 1', g, groups=g['patient']).fit(reml=False)
            lr = 2 * (full.llf - null.llf)
            row = {'name': v, 'n_patients': g['patient'].nunique(), 'n_patient_bins': len(g),
                   'lrt_p_bin_effect': stats.chi2.sf(lr, 2),
                   'coef_linear': full.params['lin'], 'p_linear': full.pvalues['lin'],
                   'coef_quadratic': full.params['quad'], 'p_quadratic': full.pvalues['quad']}
            W = g.pivot(index='patient', columns='bin', values='y')
            pm = W.mean()
            row['peak_bin'] = int(pm.idxmax()); row['trough_bin'] = int(pm.idxmin())
            for b in BINS:
                row[f'patient_mean_bin{b}'] = pm.get(b, np.nan)
            for b in BINS[1:]:
                pr = W[[0, b]].dropna()
                dlt = pr[b] - pr[0]
                row[f'n_pairs_bin{b}_vs_bin0'] = len(pr)
                row[f'median_delta_bin{b}_vs_bin0'] = dlt.median()
                row[f'n_up_bin{b}_vs_bin0'] = int((dlt > 0).sum())
                row[f'wilcoxon_p_bin{b}_vs_bin0'] = stats.wilcoxon(pr[b], pr[0]).pvalue if len(pr) >= 3 else np.nan
            rows.append(row)
            pmn = W.mean(); sem = W.sem()
            means.append(pd.DataFrame({'resource': resource, 'name': v, 'pt_bin': BINS,
                                       'mean_of_patient_means': [pmn.get(b) for b in BINS],
                                       'sem_across_patients': [sem.get(b) for b in BINS],
                                       'n_patients': [W[b].notna().sum() if b in W else 0 for b in BINS]}))
        S = pd.DataFrame(rows)
        for c in ['lrt_p_bin_effect', 'p_linear', 'p_quadratic']:
            S[c.replace('p_', 'q_', 1) if c.startswith('p_') else c.replace('_p_', '_q_')] = multipletests(S[c], method='fdr_bh')[1]
        pc = [f'wilcoxon_p_bin{b}_vs_bin0' for b in BINS[1:]]
        qv = multipletests(S[pc].fillna(1).values.ravel(), method='fdr_bh')[1].reshape(len(S), len(pc))
        for i, c in enumerate(pc):
            S[c.replace('_p_', '_q_')] = qv[:, i]
        def shape(r):
            if r['lrt_q_bin_effect'] >= 0.05:
                return 'n.s.'
            lin = r['q_linear'] < 0.05; quad = r['q_quadratic'] < 0.05
            if quad and r['coef_quadratic'] < 0 and r['peak_bin'] in (1, 2, 3):
                return 'rise then fall'
            if quad and r['coef_quadratic'] > 0 and r['trough_bin'] in (1, 2, 3):
                return 'fall then rise'
            if lin:
                return 'increasing' if r['coef_linear'] > 0 else 'decreasing'
            return 'non-linear'
        S['trajectory_shape'] = S.apply(shape, axis=1)
        S.insert(0, 'resource', resource)
        return S, pd.concat(means, ignore_index=True)


    ALL, MEANS = [], []
    for res, sc in scores.items():
        S, M = binwise_stats(sc, res)
        ALL.append(S); MEANS.append(M)
    ALL = pd.concat(ALL, ignore_index=True); MEANS = pd.concat(MEANS, ignore_index=True)
    ALL.to_csv(f'{OUT}/Fig3G_pathway_stats_all_gene_sets.csv', index=False, float_format='%.3g')
    MEANS.to_csv(f'{OUT}/Fig3G_pathway_patient_bin_means.csv', index=False, float_format='%.3g')

    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
    show = ['resource', 'name', 'n_patients', 'lrt_q_bin_effect', 'coef_linear', 'q_linear', 'coef_quadratic', 'q_quadratic',
            'peak_bin', 'trajectory_shape'] + [f'wilcoxon_q_bin{b}_vs_bin0' for b in BINS[1:]]
    for res in scores:
        print(ALL[ALL.resource == res][show].sort_values('lrt_q_bin_effect').round(4).to_string(index=False))


    print('done')

    NWIN = 50   # ~59 Tregs per window (rank-scaled diffusion pseudotime, the Fig 4E display scale; equal cells per window)

    # ---- Fig 4G panel: three groups, rows sorted by p within group.
    GROUPS = [
        ('Treg Activation & Suppressive Function', [('Suppression_panel', 'HALLMARK_IL2_STAT5_SIGNALING', 'IL-2/STAT5 signaling (Hallmark)'),
                              ('Suppression_panel', 'REACTOME_RUNX1_AND_FOXP3_CONTROL_THE_DEVELOPMENT_OF_REGULATORY_T_LYMPHOCYTES_TREGS',
                               'RUNX1/FOXP3 Treg control (Reactome)'),
                              ('Suppression_panel', 'REACTOME_CO_INHIBITION_BY_CTLA4', 'Co-inhibition by CTLA4 (Reactome)'),
                              ('Suppression_panel', 'GOBP_NEGATIVE_REGULATION_OF_T_CELL_PROLIFERATION', 'Neg. regulation of T cell proliferation (GO)'),
                              ('Suppression_panel', 'GOBP_T_CELL_TOLERANCE_INDUCTION', 'T cell tolerance induction (GO)'),
                              ('Suppression_panel', 'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY', 'TNFR2 non-canonical NF-κB (Reactome)')]),
        ('MYC / metabolism', [('MSigDB_hallmark', 'MYC_TARGETS_V1', 'MYC targets V1 (Hallmark)'),
                              ('MSigDB_hallmark', 'OXIDATIVE_PHOSPHORYLATION', 'Oxidative phosphorylation (Hallmark)'),
                              ('MSigDB_hallmark', 'MTORC1_SIGNALING', 'mTORC1 signaling (Hallmark)')]),
        ('Hypoxia / Wnt / inflammation', [('PROGENy', 'Hypoxia', 'Hypoxia (PROGENy)'),
                              ('PROGENy', 'WNT', 'Wnt (PROGENy)'),
                              ('MSigDB_hallmark', 'TNFA_SIGNALING_VIA_NFKB', 'TNF-α signaling via NF-κB (Hallmark)'),
                              ('MSigDB_hallmark', 'INFLAMMATORY_RESPONSE', 'Inflammatory response (Hallmark)')]),
    ]
    r_ = o['dpt_pseudotime'].rank(pct=True).values
    win_ = np.minimum((r_ * NWIN).astype(int), NWIN - 1)
    blocks, stat_rows = [], []
    for gname, items in GROUPS:
        rows = []
        for res, n, disp in items:
            x = np.asarray(scores[res][:, [n]].X, float).ravel(); x = (x - x.mean()) / x.std()
            st = ALL[(ALL.resource == res) & (ALL.name == n)].iloc[0]
            d, p, q = st['median_delta_bin4_vs_bin0'], st['wilcoxon_p_bin4_vs_bin0'], st['wilcoxon_q_bin4_vs_bin0']
            rows.append(dict(disp=disp, vals=pd.Series(x).groupby(win_).mean().values, p=p,
                             sig=f"{'↑' if d > 0 else '↓'}{stars(q) if stars(q) else ' n.s.'}"))
            stat_rows.append({'group': gname, 'pathway': disp, 'resource': res, 'gene_set': n,
                              'n_patients_paired': int(st['n_pairs_bin4_vs_bin0']), 'n_patients_up': int(st['n_up_bin4_vs_bin0']),
                              'median_delta_last_vs_first_quintile': d, 'p_last_vs_first_quintile': p, 'q_last_vs_first_quintile_BH': q,
                              'q_bin_effect_mixed_model_BH': st['lrt_q_bin_effect'], 'trajectory_shape': st['trajectory_shape']})
        blocks.append((gname, sorted(rows, key=lambda r: r['p'])))   # strongest (smallest p) first within group
    nr = [len(b[1]) for b in blocks]
    fig = plt.figure(figsize=(7.4, 0.33 * sum(nr) + 0.25 * len(nr) + 0.7))
    gs = fig.add_gridspec(len(nr), 1, height_ratios=nr, hspace=0.22, top=0.93)
    for k, (gname, rows) in enumerate(blocks):
        ax = fig.add_subplot(gs[k])
        M = np.vstack([r['vals'] for r in rows])
        im = ax.imshow(M, cmap='RdBu_r', vmin=-0.8, vmax=0.8, aspect='auto', interpolation='nearest', extent=(0, 1, len(rows) - 0.5, -0.5))
        for x in (0.2, 0.4, 0.6, 0.8):
            ax.axvline(x, color='k', lw=0.6, ls=':', alpha=0.6)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([r['disp'] for r in rows], fontsize=10.5)
        for j, r in enumerate(rows):
            ax.text(1.015, j, r['sig'], transform=ax.get_yaxis_transform(), va='center', ha='left', fontsize=10.5)
        ax.set_title(gname, loc='left', fontsize=10.5, fontweight='bold', pad=3)
        if k == 0:
            ax.text(1.015, 1.0, 'BH FDR q', transform=ax.transAxes, va='bottom', ha='left', fontsize=9.5, fontweight='bold')
        if k < len(blocks) - 1:
            ax.set_xticks([])
        else:
            ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_xlabel('Diffusion Pseudotime (rank-scaled)', fontsize=12)
    fig.suptitle('Pathway Activity Along Treg Pseudotime', fontsize=13, y=1.0)
    cax = fig.add_axes([1.03, 0.3, 0.018, 0.4]); cb = fig.colorbar(im, cax=cax); cb.set_label('Activity (z-score)', fontsize=9)
    fig.savefig(f'{OUT}/Fig3G_pathways_along_pseudotime.png', dpi=300, bbox_inches='tight'); plt.close(fig)
    pd.DataFrame(stat_rows).to_csv(f'{OUT}/Fig3G_pathways_along_pseudotime_stats.csv', index=False, float_format='%.3g')
    vals = pd.DataFrame([r['vals'] for _, rows in blocks for r in rows], index=[r['disp'] for _, rows in blocks for r in rows],
                        columns=[f'window_{i + 1}' for i in range(NWIN)])
    vals.to_csv(f'{OUT}/Fig3G_pathways_along_pseudotime_values.csv', float_format='%.3g')
    print('combined done')


# ======================================================================================================================
# Section S8F_tf_correlation: CollecTRI TF activity vs Treg pseudotime (within-patient Spearman, sign test).
# Panels: Supp Fig 8F (values); also read by treatment_validation/03
# ======================================================================================================================
if _run('S8F_tf_correlation'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/figure3/S8F_tf_correlation.py'   # working folder of this step (relative paths resolve from here)
    """TFs most correlated with Treg pseudotime. CollecTRI ULM activity (layer log1p) for all TFs; within-patient Spearman
rho with diffusion pseudotime (patients with >= 30 Tregs), median across patients, n patients per sign, sign test, BH.
Eligible: >= 10 expressed targets and TF gene detected in >= 5% of Tregs (reported both with and without this filter).
Plot: top 10 positive and top 10 negative by median rho (eligible, q < 0.05), heatmap of mean activity in 50 equal-size
windows of rank-scaled pseudotime, z-scored per TF; right column = median within-patient rho and n patients in that
direction. Outputs panels/figure3/tf_pseudotime_correlation/ (tf_pseudotime_correlation.csv is read by S8F_tf_heatmap)."""
    import os, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    import numpy as np, pandas as pd, anndata as ad, decoupler as dc, matplotlib
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = f'{D}/panels/figure3/tf_pseudotime_correlation'; os.makedirs(OUT, exist_ok=True)
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad'); o = t.obs; pt = o['dpt_pseudotime'].values; pat = o['patient'].astype(str).values
    big = [p for p in np.unique(pat) if (pat == p).sum() >= 30]
    net = pd.read_csv(f'{D}/objects/collectri_human_snapshot.csv.gz'); dc.mt.ulm(data=t, net=net, layer='log1p'); A = dc.pp.get_obsm(t, key='score_ulm')
    X = pd.DataFrame(np.asarray(A.X), index=A.obs_names, columns=A.var_names)
    L = t.layers['log1p']; det = pd.Series(np.asarray((L > 0).mean(0)).ravel() * 100, index=t.var_names)
    ntarg = net[net.target.isin(t.var_names)].groupby('source').size()
    rows = []
    for tf in X.columns:
        rh = np.array([stats.spearmanr(X[tf].values[pat == p], pt[pat == p])[0] for p in big]); rh = rh[np.isfinite(rh)]
        npos = int((rh > 0).sum())
        rows.append(dict(tf=tf, n_targets=int(ntarg.get(tf, 0)), tf_gene_detected_pct=det.get(tf, np.nan), median_rho=np.median(rh), n_patients=len(rh),
                         n_positive=npos, n_negative=len(rh) - npos, p_sign=stats.binomtest(npos, len(rh)).pvalue, rho_cells=stats.spearmanr(X[tf].values, pt)[0]))
    S = pd.DataFrame(rows); S['q_sign'] = multipletests(S.p_sign, method='fdr_bh')[1]
    S['eligible'] = (S.n_targets >= 10) & (S.tf_gene_detected_pct >= 5)
    S = S.reindex(S.median_rho.abs().sort_values(ascending=False).index); S.to_csv(f'{OUT}/tf_pseudotime_correlation.csv', index=False, float_format='%.3g')
    pd.set_option('display.width', 220)
    c = ['tf', 'n_targets', 'tf_gene_detected_pct', 'median_rho', 'n_positive', 'n_negative', 'q_sign', 'rho_cells']
    print('== all TFs, top 15 by |median within-patient rho|'); print(S.head(15)[c].round(3).to_string(index=False))
    E = S[S.eligible & (S.q_sign < 0.05)]
    up = E[E.median_rho > 0].head(10); dn = E[E.median_rho < 0].head(10)
    print('== eligible rising'); print(up[c].round(3).to_string(index=False)); print('== eligible falling'); print(dn[c].round(3).to_string(index=False))
    print('== reference'); print(S[S.tf.isin(['RORC', 'AP1', 'FOXP3', 'STAT5A', 'STAT5B', 'FOXO3', 'RFX5', 'MYC'])][c].round(3).to_string(index=False))
    sel = pd.concat([dn.sort_values('median_rho'), up.sort_values('median_rho')])
    NW = 50; win = np.minimum((pd.Series(pt).rank(pct=True).values * NW).astype(int), NW - 1)
    M = np.vstack([pd.Series(((X[tf] - X[tf].mean()) / X[tf].std()).values).groupby(win).mean().values for tf in sel.tf])
    fig, ax = plt.subplots(figsize=(7.6, 0.3 * len(sel) + 1.4))
    im = ax.imshow(M, cmap='RdBu_r', vmin=-0.8, vmax=0.8, aspect='auto', interpolation='nearest', extent=(0, 1, len(sel) - 0.5, -0.5))
    ax.axhline(len(dn) - 0.5, color='white', lw=4)
    for x in (0.2, 0.4, 0.6, 0.8): ax.axvline(x, color='k', lw=0.6, ls=':', alpha=0.6)
    ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel.tf, fontsize=10)
    for j, r in enumerate(sel.itertuples()):
        ax.text(1.015, j, f"ρ {r.median_rho:+.2f}  ({r.n_negative if r.median_rho < 0 else r.n_positive}/{r.n_patients})", transform=ax.get_yaxis_transform(), va='center', fontsize=9)
    ax.text(1.015, 1.0, 'median within-patient ρ\n(patients in direction)', transform=ax.transAxes, va='bottom', fontsize=8.5, fontweight='bold')
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_xlabel('Diffusion Pseudotime (rank-scaled)', fontsize=12)
    ax.set_title('TF Activity Most Correlated With Treg Pseudotime', fontsize=13)
    cax = fig.add_axes([1.16, 0.3, 0.02, 0.4]); cb = fig.colorbar(im, cax=cax); cb.set_label('Activity (z-score)', fontsize=9)
    fig.savefig(f'{OUT}/tf_pseudotime_correlation_top.png', dpi=300, bbox_inches='tight'); plt.close(fig)
    sel.to_csv(f'{OUT}/tf_pseudotime_correlation_plotted.csv', index=False, float_format='%.3g')


# ======================================================================================================================
# Section S8F_tf_heatmap: Top falling and rising TFs along pseudotime with their driver targets.
# Panels: Supp Fig 8F
# ======================================================================================================================
if _run('S8F_tf_heatmap'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness/S8F_tf_heatmap.py'   # working folder of this step (relative paths resolve from here)
    """Supp 8F TF panel: TFs whose transcript is detected in >= 1% of Tregs, with >= 10 CollecTRI targets and sign-test
q < 0.05; top 10 each direction by median within-patient rho. Group marks: red = suppressive / activated program,
falling; blue = low-suppressive / Th17-like program, rising; top driver targets per TF.
q = Benjamini-Hochberg-adjusted p of the per-patient sign test across all 676 CollecTRI TFs (section S8F_tf_correlation).
Driver targets: per target gene (layer log1p), per-patient mean in bin 0 and bin 4 (>= 3 Tregs each) -> paired
difference, Wilcoxon, standardised difference = median patient difference / gene SD across cells; contribution to the
TF's ULM change = CollecTRI mode of regulation (+1 activation / -1 repression) x standardised difference. Drivers =
targets whose contribution has the TF's direction, BH q < 0.05 within the TF's targets, gene detected in >= 1% of Tregs;
top 4 by |contribution|; if none pass, top 3 with nominal p < 0.05 shown grey with a dagger.
'⊣' = repressed target (moves opposite to the TF).
Outputs: ext8_tf_activity_along_pseudotime.png, _drivers.csv, _legend.txt (SuppFig8A-I folder) and
panels/extended_data_fig8/ext8_tf_activity_along_pseudotime.csv (read by treatment_validation/03)."""
    import os, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    import numpy as np, pandas as pd, anndata as ad, decoupler as dc, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from scipy import stats, sparse
    from statsmodels.stats.multitest import multipletests
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = '/share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness'
    C = pd.read_csv(f'{D}/panels/figure3/tf_pseudotime_correlation/tf_pseudotime_correlation.csv').set_index('tf')
    CF = C[(C.tf_gene_detected_pct >= 1) & (C.n_targets >= 10) & (C.q_sign < 0.05)]
    GROUPS = [('Suppressive / activated TFs (falling)', list(CF.sort_values('median_rho').index[:10]), '#c0392b'),
              ('Low-suppressive / Th17-like TFs (rising)', list(CF.sort_values('median_rho', ascending=False).index[:10]), '#2e6da4')]
    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad'); o = t.obs; pt = o['dpt_pseudotime'].values; b = o['pt_bin'].astype(int).values; pat = o['patient'].astype(str).values
    net = pd.read_csv(f'{D}/objects/collectri_human_snapshot.csv.gz'); dc.mt.ulm(data=t, net=net, layer='log1p'); A = dc.pp.get_obsm(t, key='score_ulm')
    X = pd.DataFrame(np.asarray(A.X), index=A.obs_names, columns=A.var_names)
    L = t.layers['log1p']; L = L.tocsc() if sparse.issparse(L) else sparse.csc_matrix(L)
    cache = {}
    def gene_stats(gn):
        if gn in cache: return cache[gn]
        x = np.asarray(L[:, t.var_names.get_loc(gn)].todense()).ravel()
        d = pd.DataFrame({'p': pat, 'b': b, 'x': x}); g = d.groupby(['p', 'b']).x.agg(['mean', 'size']).reset_index(); g = g[g['size'] >= 3]
        W = g.pivot(index='p', columns='b', values='mean'); pr = W[[0, 4]].dropna() if {0, 4} <= set(W.columns) else pd.DataFrame(columns=[0, 4])
        dl = pr[4] - pr[0]; sd = x.std() if x.std() > 0 else np.nan
        out = dict(n_patients=len(pr), median_patient_delta=dl.median(), mean_patient_delta=dl.mean(), n_up=int((dl > 0).sum()),
                   p=stats.wilcoxon(pr[4], pr[0]).pvalue if len(pr) >= 5 and dl.abs().sum() > 0 else np.nan,
                   std_delta=dl.median() / sd if sd == sd else np.nan, gene_sd=sd, pct_detected=float((x > 0).mean() * 100))
        cache[gn] = out; return out
    drows = []
    for _, tfs, _ in GROUPS:
        for tf in tfs:
            sgn = np.sign(C.loc[tf, 'median_rho']); tg = net[(net.source == tf) & net.target.isin(t.var_names)]
            R = pd.DataFrame([dict(tf=tf, gene=r.target, mor=r.weight, **gene_stats(r.target)) for r in tg.itertuples()])
            R['contribution'] = R.mor * R.std_delta
            R['contribution_mean'] = R.mor * R.mean_patient_delta / R.gene_sd
            ok = R.p.notna(); R['q'] = np.nan
            if ok.any(): R.loc[ok, 'q'] = multipletests(R.loc[ok, 'p'], method='fdr_bh')[1]
            R['is_driver'] = (np.sign(R.contribution) == sgn) & (R.q < 0.05) & (R.pct_detected >= 1)
            drows.append(R)
    DR = pd.concat(drows); DR['abs_contrib'] = DR.contribution.abs()
    DR.sort_values(['tf', 'is_driver', 'abs_contrib'], ascending=[True, False, False]).drop(columns='abs_contrib').to_csv(f'{OUT}/ext8_tf_activity_along_pseudotime_drivers.csv', index=False, float_format='%.3g')
    TOP = {tf: d[d.is_driver].sort_values('abs_contrib', ascending=False).head(4) for tf, d in DR.groupby('tf')}
    NOM = {tf: d[(np.sign(d.contribution) == np.sign(C.loc[tf, 'median_rho'])) & (d.p < 0.05) & (d.pct_detected >= 1)].sort_values('abs_contrib', ascending=False).head(3)
           for tf, d in DR.groupby('tf')}   # fallback when no target passes BH q < 0.05 (shown grey, dagger)
    def gene_items(tf):
        """>= 2 targets: p < 0.05 targets moving in the TF's direction (black, up to 4, by contribution), padded with the next largest
    direction-consistent contributors (grey; ranked by median-based contribution, ties at 0 broken by mean-based contribution)."""
        sg = np.sign(C.loc[tf, 'median_rho']); d = DR[(DR.tf == tf) & (DR.pct_detected >= 1)].copy()
        d['key'] = np.where(d.contribution.fillna(0) != 0, d.contribution.fillna(0), 0) * sg
        d['key2'] = d.contribution_mean.fillna(0) * sg
        sig = d[(d.key > 0) & (d.p < 0.05)].sort_values('key', ascending=False).head(4)
        rest = d[~d.gene.isin(sig.gene) & ((d.key > 0) | ((d.key == 0) & (d.key2 > 0)))].sort_values(['key', 'key2'], ascending=False)
        pad = rest.head(max(0, 2 - len(sig)))
        return [(('⊣' if m < 0 else '') + g, 'k') for g, m in zip(sig.gene, sig.mor)] + [(('⊣' if m < 0 else '') + g, 'grey') for g, m in zip(pad.gene, pad.mor)]
    def lab_genes(tf):
        it = gene_items(tf); return ', '.join(g for g, _ in it), any(c == 'grey' for _, c in it)
    def draw_genes(ax, j, items):
        t = ax.text(1.215, j, items[0][0] + (',' if len(items) > 1 else ''), transform=ax.get_yaxis_transform(), va='center', fontsize=9, style='italic', color=items[0][1])
        for k, (g, c) in enumerate(items[1:], 1):
            t = ax.annotate(' ' + g + (',' if k < len(items) - 1 else ''), xy=(1, 0.5), xycoords=t, va='center', ha='left', fontsize=9, style='italic', color=c)
    NW = 50; win = np.minimum((pd.Series(pt).rank(pct=True).values * NW).astype(int), NW - 1)
    def stars(q): return '***' if q < 0.001 else '**' if q < 0.01 else '*'
    nr = [len(g) for _, g, _ in GROUPS]
    fig = plt.figure(figsize=(9.6, 0.36 * sum(nr) + 0.25 * len(nr) + 0.9)); gs = fig.add_gridspec(len(nr), 1, height_ratios=nr, hspace=0.12, top=0.92)
    rows = []
    for k, (gname, tfs, col) in enumerate(GROUPS):
        ax = fig.add_subplot(gs[k])
        M = np.vstack([pd.Series(((X[tf] - X[tf].mean()) / X[tf].std()).values).groupby(win).mean().values for tf in tfs])
        im = ax.imshow(M, cmap='RdBu_r', vmin=-0.8, vmax=0.8, aspect='auto', interpolation='nearest', extent=(0, 1, len(tfs) - 0.5, -0.5))
        for x in (0.2, 0.4, 0.6, 0.8): ax.axvline(x, color='k', lw=0.6, ls=':', alpha=0.6)
        ax.set_yticks(range(len(tfs))); ax.set_yticklabels(tfs, fontsize=10.5)
        for j, tf in enumerate(tfs):
            r = C.loc[tf]; nd = int(r.n_negative if r.median_rho < 0 else r.n_positive)
            ax.text(1.015, j, f"{'↓' if r.median_rho < 0 else '↑'}{stars(r.q_sign)}  ρ {r.median_rho:+.2f} ({nd}/{int(r.n_patients)})", transform=ax.get_yaxis_transform(), va='center', fontsize=9.5)
            gl, nom = lab_genes(tf); draw_genes(ax, j, gene_items(tf))
            rows.append(dict(group=gname, tf=tf, median_rho=r.median_rho, n_in_direction=nd, n_patients=int(r.n_patients), q_sign=r.q_sign, top_driver_targets=gl))
        ax.set_title(gname, loc='left', fontsize=10.5, fontweight='bold', pad=3, color=col)
        if k == 0:
            ax.text(1.015, 1.0, 'Median ρ (patients)', transform=ax.transAxes, va='bottom', fontsize=9, fontweight='bold')
            ax.text(1.215, 1.0, 'Top driver targets', transform=ax.transAxes, va='bottom', fontsize=9, fontweight='bold')
        if k < len(GROUPS) - 1: ax.set_xticks([])
        else: ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_xlabel('Diffusion Pseudotime (rank-scaled)', fontsize=12)
    fig.text(0.04, 0.5, 'Transcription Factor', rotation=90, va='center', ha='center', fontsize=12)   # title and footnote are in the legend text
    cax = fig.add_axes([1.39, 0.3, 0.014, 0.4]); cb = fig.colorbar(im, cax=cax); cb.set_label('TF activity (z-score)', fontsize=9)
    open(f'{OUT}/ext8_tf_activity_along_pseudotime_legend.txt', 'w').write("Rows, CollecTRI transcription factors (TFs) whose transcript is detected in >=1% of Tregs and that have >=10 targets; top 10 per direction by median within-patient Spearman rho between TF activity (decoupler univariate linear model on log-normalized expression) and diffusion pseudotime. Columns, 50 equal-size windows of rank-scaled pseudotime; colour, TF activity z-scored across Tregs. rho, median within-patient correlation (patients with that sign / patients tested, patients with >=30 Tregs); stars, Benjamini-Hochberg-adjusted two-sided sign test across patients over all 676 TFs (*q < 0.05, **q < 0.01, ***q < 0.001). Top driver targets, targets changing in the TF's direction between pseudotime bins 0 and 4 (per-patient means, >=3 Tregs per bin; Wilcoxon signed-rank p < 0.05; black), ranked by contribution to the activity change (mode of regulation x median per-patient change / gene SD); when fewer than two targets reach p < 0.05, the next largest contributors are shown in grey; the symbol before a gene marks a repressed target (moves opposite to the TF)." + '\n')
    fig.savefig(f'{OUT}/ext8_tf_activity_along_pseudotime.png', dpi=300, bbox_inches='tight'); plt.close(fig)
    pd.DataFrame(rows).to_csv(f'{D}/panels/extended_data_fig8/ext8_tf_activity_along_pseudotime.csv', index=False, float_format='%.3g'); print(pd.DataFrame(rows).round(3).to_string(index=False))
