#!/usr/bin/env python
"""CD CosMx 6K resections: cellular-neighbourhood (CN) composition and the Treg Suppression Score inside vs outside each CN.

Inputs  panels/fig4_and_related/D_shared_inputs/treg_definition_comparison_per_cell.csv (10);
        B_resection_analyses/cellular_neighborhoods_v2/ cn_per_cell.parquet, cn_treg_state_enrichment_all.csv (11)
Outputs panels/fig4_and_related/A_fig4_panels/fig4_cn_composition_heatmap.csv (Supp 9A input; read by 13, 17),
        fig4_cn_treg_state_forest.csv (Fig 5C left source data; read by 13, 15, 17)
Sections cn_composition_heatmap -> F5C_cn_forest (RUN_SECTIONS=<id>,... runs a subset)
"""
import os as _os
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


import numpy as np, pandas as pd, pyarrow.parquet as pq
D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/A_fig4_panels'
CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name   # 0-based CN names


# ---- Section cn_composition_heatmap: log2 enrichment of each cell type per CN (k = 10) vs all cells, pseudocount 1e-4;
#      Mo-Mac + Inf Mo-Mac, CD8 GZMK+/- Trm and Tfh + Exhausted CD4 Tfh merged; columns grouped by lineage; compact CN labels
if _run('cn_composition_heatmap'):
    NMS = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [('OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    for k in NMS:
        for a, b in TRIM: NMS[k] = NMS[k].replace(a, b)
    MERGE = {'Inf Mo-Mac': 'Mo-Mac', 'CD8 GZMK+ Trm': 'CD8 Trm', 'CD8 GZMK- Trm': 'CD8 Trm', 'Exhausted CD4 Tfh': 'Tfh'}
    cc = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['cell_type', 'CN_k10']).to_pandas()
    t = pd.crosstab(cc.CN_k10, cc.cell_type.replace(MERGE)); fr = t.div(t.sum(1), axis=0)
    E = np.log2((fr + 1e-4) / (t.sum(0) / t.values.sum() + 1e-4)); E.index = [NMS[i] for i in E.index]
    LIN = {'Epithelial': ['Epithelial'],
           'Stromal': ['OGN_RSPO3 Fib', 'PLA2G2A Fib', 'ADAMDEC1 Fib', 'Activated Fib', 'Crypt Top Fib', 'Myofibroblast', 'Pericyte', 'Endothelial', 'Neuron'],
           'Myeloid': ['DC', 'CD9 Mac', 'Mac S', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac', 'Neutrophil', 'Granulocyte'],
           'Lymphoid': ['Treg', 'CD4 Trm', 'CD4 Effector', 'Tfh', 'Resting T', 'CD8 Trm', 'Gd T', 'NK', 'B', 'Plasma'],
           'Other': ['Cycling']}
    cols = [c for v in LIN.values() for c in v if c in E.columns]
    E[cols].to_csv(f'{OUT}/fig4_cn_composition_heatmap.csv', float_format='%.3g')


# ---- Section F5C_cn_forest: per CN, LMM activation_z ~ in_CN + (1 | section) on all resection Tregs (activation_z = Treg
#      Suppression Score = -signature_score_D6 z within section); tested if >= 100 Tregs and >= 10 Tregs in >= 5 sections; BH over
#      tested CNs; sections agreeing in direction (>= 3 Tregs in the CN); CMH odds ratios (median and thirds splits) as sensitivity
if _run('F5C_cn_forest'):
    import warnings; warnings.filterwarnings('ignore')
    from statsmodels.stats.contingency_tables import StratifiedTable
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    MIN_TREG_TOTAL = 100; MIN_TREG_SEC = 10; MIN_SEC = 5   # pre-stated testing rule
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); pr = PC[PC.cohort == 'resection']
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas()
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
    r.rename(columns={'activation_z_in_minus_out': 'treg_suppression_score_z_in_minus_out'}).to_csv(f'{OUT}/fig4_cn_treg_state_forest.csv', index=False, float_format='%.3g')
    print(r[['CN', 'tested', 'n_tregs_in_CN', 'activation_z_in_minus_out', 'mixed_q', 'sections_agreeing', 'sections']].round(4).to_string(index=False))
