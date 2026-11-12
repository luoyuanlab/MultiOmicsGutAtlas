#!/usr/bin/env python
"""UC scRNA-seq atlas: NicheNet ligands from OGN+RSPO3+ fibroblasts to SELENOP+MMP9+SPP1+ (SMS) macrophages.

Panels: Supp 5C (NicheNet ligand activity, fibroblast enrichment, patient-level ligand-response test);
        supp tables of the NicheNet ligands and ligand rank stability.
Inputs: UC atlas h5ad; NicheNet ligand_target_matrix_nsga2r_final.rds and lr_network_human_21122021.rds.
Outputs: SuppFig5C_nichenet_fib_to_mac_ligands/ (DE table, intermediate/ caches, CSVs, panel).
Run: python <file>, or RUN_SECTIONS=<id>[,<id>...] python <file>.
"""
import os
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
warnings.filterwarnings('ignore')

_RUN_SECTIONS = [s.strip() for s in os.environ.get('RUN_SECTIONS', '').split(',') if s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


ATLAS = '/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad'
O = '/share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands'
NN = '/share/fsmresfiles/ydn4687/nichenet'
TARGET = 'Mac S+M+S+'; REF = ['Mac S+XS-', 'Mac S+SG+']       # SMS macs vs pooled resident SELENOP+ macs
PB_CACHE = f'{O}/intermediate/intermediate_patient_pseudobulk_counts.csv.gz'
FIB = ['ADAMDEC1+ Fib', 'VSTM2A+ Crypt Top Fib', 'OGN+RSPO3+ Fib', 'NRG1+ Crypt Top Fib', 'Activ Fib', 'CXCL5+ Activ Fib',
       'T-interact Fib', 'SELENOP+ Fib']
SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')


def sig3(v):
    return f'{v:#.3g}'.rstrip('.')


def pfmt(v):
    if v >= 0.001:
        return sig3(v)
    m, e = f'{v:.2e}'.split('e'); return f'{m} × 10{str(int(e)).translate(SUP)}'


def read_rds(name):
    import pyreadr
    return list(pyreadr.read_r(f'{NN}/{name}').values())[0]


def expressed(scrna, ct, thr=0.10):
    """Genes detected in > thr of the cells of one cell type."""
    sub = scrna[np.flatnonzero((scrna.obs['cell_type_short'] == ct).values)]
    X = sub.to_memory().X if scrna.isbacked else sub.X
    return set(scrna.var_names[np.asarray((X > 0).mean(0)).flatten() > thr])


def nichenet_inputs(scrna):
    """Fibroblast ligands with an SMS-expressed receptor, the SMS-expressed NicheNet background, ltm and the LR pairs."""
    fib_expr = expressed(scrna, 'OGN+RSPO3+ Fib'); mac_expr = expressed(scrna, TARGET)
    ltm = read_rds('ligand_target_matrix_nsga2r_final.rds')
    lr = read_rds('lr_network_human_21122021.rds'); lr.columns = [c.lower() for c in lr.columns]
    lr_ok = lr[lr['from'].isin(fib_expr) & lr['to'].isin(mac_expr)]; ligs = sorted(set(lr_ok['from']) & set(ltm.columns))
    bg = sorted(set(ltm.index) & mac_expr)
    return ltm, lr_ok, ligs, bg, mac_expr


def ligand_activity(ltm, bg, ligs, geneset):
    """NicheNet ligand activity: AUPR / AUROC of each ligand's regulatory potential for the gene set."""
    resp = np.array([g in geneset for g in bg], int); rows = []
    for L in ligs:
        tv = ltm.loc[bg, L].values.astype(float)
        if len(np.unique(tv)) < 2: continue
        rows.append([L, average_precision_score(resp, tv), roc_auc_score(resp, tv)])
    return pd.DataFrame(rows, columns=['ligand', 'aupr', 'auroc']).sort_values('aupr', ascending=False).reset_index(drop=True), resp


# Section S5C_pseudobulk_de: patient pseudobulk DESeq2, SMS vs pooled resident SELENOP+ macs (>= 10 cells, genes >= 10 counts)
if _run('S5C_pseudobulk_de'):
    import scanpy as sc, scipy.sparse as sp
    from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
    scrna = sc.read_h5ad(ATLAS)
    mye = scrna[scrna.obs['cell_category'] == 'Myeloid'].copy()
    sub = mye[mye.obs['cell_type_short'].isin([TARGET] + REF)].copy()
    raw = sub.raw.to_adata()[sub.obs_names]; genes = sub.raw.var['features'].astype(str).values
    X = raw.X; X = X.toarray() if sp.issparse(X) else np.asarray(X)
    cnt = pd.DataFrame(X, index=sub.obs_names, columns=genes); cnt = cnt.loc[:, ~cnt.columns.duplicated()]
    meta = sub.obs[['Patient', 'cell_type_short']].copy(); meta['grp'] = np.where(meta['cell_type_short'] == TARGET, 'target', 'ref')
    cnt['pb'] = meta['Patient'].astype(str) + '__' + meta['grp']; pb = cnt.groupby('pb').sum()
    pbmeta = pd.DataFrame({'grp': [i.split('__')[1] for i in pb.index], 'patient': [i.split('__')[0] for i in pb.index]}, index=pb.index)
    ncell = meta.groupby(meta['Patient'].astype(str) + '__' + meta['grp']).size()
    keep = ncell[ncell >= 10].index; pb = pb.loc[pb.index.isin(keep)]; pbmeta = pbmeta.loc[pb.index]
    pb = pb.loc[:, pb.sum(0) >= 10].round().astype(int)
    dds = DeseqDataSet(counts=pb, metadata=pbmeta, design_factors='grp', ref_level=['grp', 'ref'], quiet=True)
    dds.deseq2(); st = DeseqStats(dds, contrast=['grp', 'target', 'ref'], quiet=True); st.summary()
    st.results_df.sort_values('padj').to_csv(f'{O}/pseudobulk_DE_homeoSp_pt.csv')


# Section S5C_geneset_size_sensitivity: NicheNet ligand ranks for four definitions of the SMS target gene set
if _run('S5C_geneset_size_sensitivity'):
    import scanpy as sc
    de = pd.read_csv(f'{O}/pseudobulk_DE_homeoSp_pt.csv', index_col=0).dropna(subset=['padj'])
    scrna = sc.read_h5ad(ATLAS, backed='r')
    ltm, _, ligs, bg, _ = nichenet_inputs(scrna)
    up = de[de.log2FoldChange > 0]
    SETS = {'padj<0.1 & LFC>0.5 (default)': up[(up.padj < 0.1) & (up.log2FoldChange > 0.5)].index,
            'padj<0.05 & LFC>1': up[(up.padj < 0.05) & (up.log2FoldChange > 1)].index,
            'top 300 by LFC (padj<0.05)': up[up.padj < 0.05].sort_values('log2FoldChange', ascending=False).index[:300],
            'top 100 by LFC (padj<0.05)': up[up.padj < 0.05].sort_values('log2FoldChange', ascending=False).index[:100]}
    out = []
    for name, gs in SETS.items():
        gs = set(gs) & set(bg)
        a, resp = ligand_activity(ltm, bg, ligs, gs)
        a['rank'] = a.index + 1; a['set'] = name; a['n_genes'] = len(gs); a['aupr_baseline'] = resp.mean(); out.append(a)
    pd.concat(out).to_csv(f'{O}/intermediate/intermediate_ligand_rank_by_geneset_size.csv', index=False)


# Section S5C_nichenet_ranking: NicheNet ligand activity for the SMS program (padj<0.05 & log2FC>1); rank stability table
if _run('S5C_nichenet_ranking'):
    import scanpy as sc
    de = pd.read_csv(f'{O}/pseudobulk_DE_homeoSp_pt.csv', index_col=0).dropna(subset=['padj'])
    geneset = set(de[(de.padj < 0.05) & (de.log2FoldChange > 1)].index)
    scrna = sc.read_h5ad(ATLAS, backed='r')
    ltm, lr_ok, ligs, bg, mac_expr = nichenet_inputs(scrna)
    act, _ = ligand_activity(ltm, bg, ligs, geneset)
    fibs = scrna[np.flatnonzero(scrna.obs['cell_type_short'].isin(FIB).values)].to_memory()
    isogn = (fibs.obs['cell_type_short'] == 'OGN+RSPO3+ Fib').values
    mexpr = lambda mask, g: float(np.asarray(fibs[mask, g].X.mean()).flatten()[0])
    act['ogn_mean'] = [mexpr(isogn, g) for g in act.ligand]; act['otherfib_mean'] = [mexpr(~isogn, g) for g in act.ligand]
    act['fib_enriched'] = act['ogn_mean'] > act['otherfib_mean']
    act['mac_receptors'] = [','.join(sorted(set(lr_ok[lr_ok['from'] == g]['to']) & mac_expr)) for g in act.ligand]
    act['rank'] = act.index + 1
    act.to_csv(f'{O}/nichenet_ligand_activity_patient.csv', index=False)
    keep = act.head(15)['ligand'].tolist()
    st = pd.read_csv(f'{O}/intermediate/intermediate_ligand_rank_by_geneset_size.csv')
    piv = st[st.ligand.isin(keep)].pivot(index='ligand', columns='set', values='rank').reindex(keep)
    order = ['padj<0.05 & LFC>1', 'padj<0.1 & LFC>0.5 (default)', 'top 300 by LFC (padj<0.05)', 'top 100 by LFC (padj<0.05)']
    nm = {'padj<0.05 & LFC>1': 'Rank, padj<0.05 and LFC>1 (main)', 'padj<0.1 & LFC>0.5 (default)': 'Rank, padj<0.1 and LFC>0.5',
          'top 300 by LFC (padj<0.05)': 'Rank, top 300 genes', 'top 100 by LFC (padj<0.05)': 'Rank, top 100 genes'}
    piv = piv[order].rename(columns=nm).reset_index().rename(columns={'ligand': 'Ligand'})
    piv.to_csv(f'{O}/supp_table_ligand_rank_stability.csv', index=False)


# Section S5C_patient_bootstrap: patient-bootstrap stability of the ligand ranks (100 resamples of patients, DE redone each time)
if _run('S5C_patient_bootstrap'):
    import scanpy as sc, scipy.sparse as sp
    from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
    NBOOT = 100; rng = np.random.default_rng(0)
    act = pd.read_csv(f'{O}/nichenet_ligand_activity_patient.csv'); LIGS = act['ligand'].tolist()
    scrna = sc.read_h5ad(ATLAS)
    ob = scrna.obs
    if os.path.exists(PB_CACHE):
        PB = pd.read_csv(PB_CACHE, index_col=0)
    else:
        sub = scrna[ob['cell_type_short'].isin([TARGET] + REF).values].copy()
        raw = sub.raw.to_adata()[sub.obs_names]; genes = sub.raw.var['features'].astype(str).values
        grp = np.where(sub.obs['cell_type_short'] == TARGET, 'target', 'ref'); key = sub.obs['Patient'].astype(str).values + '__' + grp
        keys = pd.Series(key); n = keys.value_counts(); keep = n[n >= 10].index
        rows = {}
        X = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
        for k in keep:
            rows[k] = np.asarray(X[np.flatnonzero(keys.values == k)].sum(0)).flatten()
        PB = pd.DataFrame(rows, index=genes).T; PB = PB.loc[:, ~PB.columns.duplicated()]
        PB.to_csv(PB_CACHE)
    units = pd.Series(PB.index); pat = units.str.split('__').str[0]; grp = units.str.split('__').str[1]
    patients = sorted(pat.unique())
    mac_expr = expressed(scrna, TARGET)                       # NicheNet background as in S5C_nichenet_ranking
    ltm = read_rds('ligand_target_matrix_nsga2r_final.rds')
    bg = sorted(set(ltm.index) & mac_expr); Mlt = ltm.loc[bg, LIGS].values.astype(float)

    def rank_once(units_idx, tag):
        C = PB.iloc[units_idx].copy(); C.index = [f'{tag}{i}' for i in range(len(C))]
        meta = pd.DataFrame({'grp': grp.iloc[units_idx].values}, index=C.index)
        if (meta.grp == 'target').sum() < 3 or (meta.grp == 'ref').sum() < 3: return None
        C = C.loc[:, C.sum(0) >= 10].round().astype(int)
        dds = DeseqDataSet(counts=C, metadata=meta, design_factors='grp', ref_level=['grp', 'ref'], quiet=True); dds.deseq2()
        st = DeseqStats(dds, contrast=['grp', 'target', 'ref'], quiet=True); st.summary(); r = st.results_df.dropna(subset=['padj'])
        gs = set(r[(r.padj < 0.05) & (r.log2FoldChange > 1)].index); resp = np.array([g in gs for g in bg], int)
        if resp.sum() < 10: return None
        au = np.array([average_precision_score(resp, Mlt[:, j]) for j in range(len(LIGS))])
        return pd.Series(au, index=LIGS).rank(ascending=False, method='min')
    ranks = {}
    for b in range(NBOOT):
        pick = rng.choice(patients, size=len(patients), replace=True)
        idx = np.concatenate([np.flatnonzero(pat.values == p) for p in pick])
        r = rank_once(idx, f'b{b}_')
        if r is not None: ranks[b] = r
    R = pd.DataFrame(ranks)
    stab = pd.DataFrame({'top5': (R <= 5).mean(1) * 100, 'top10': (R <= 10).mean(1) * 100, 'median_rank': R.median(1)})
    act = act.set_index('ligand').join(stab).reset_index(); act.to_csv(f'{O}/nichenet_ligand_activity_patient.csv', index=False)


# Section S5C_ligand_response_test: ULM ligand-response score per patient pseudobulk (top 250 NicheNet targets),
# SMS vs resident: paired Wilcoxon and score ~ group + (1|patient), BH across ligands
if _run('S5C_ligand_response_test'):
    import anndata as ad, decoupler as dc
    import statsmodels.formula.api as smf
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    NTOP = 250
    PB = pd.read_csv(PB_CACHE, index_col=0)
    units = pd.Series(PB.index, index=PB.index); pat = units.str.split('__').str[0]; grp = units.str.split('__').str[1]
    cpm = PB.div(PB.sum(1), axis=0) * 1e6; L = np.log1p(cpm)
    L = L.loc[:, (PB > 0).mean(0) >= 0.2]
    act = pd.read_csv(f'{O}/nichenet_ligand_activity_patient.csv'); LIGS = act['ligand'].tolist()
    ltm = read_rds('ligand_target_matrix_nsga2r_final.rds')
    genes = [g for g in ltm.index if g in L.columns]; rows = []
    for lg in LIGS:
        w = ltm.loc[genes, lg].sort_values(ascending=False).head(NTOP)
        rows += [(lg, g, float(v)) for g, v in w.items() if v > 0]
    net = pd.DataFrame(rows, columns=['source', 'target', 'weight'])
    A = ad.AnnData(L.values.astype(float), obs=pd.DataFrame({'patient': pat.values, 'group': grp.values}, index=L.index), var=pd.DataFrame(index=L.columns))
    dc.mt.ulm(data=A, net=net, tmin=10)
    S = dc.pp.get_obsm(A, key='score_ulm'); S = pd.DataFrame(S.X, index=A.obs_names, columns=S.var_names)
    S['patient'] = pat.values; S['group'] = grp.values
    both = [p for p, g in S.groupby('patient') if set(g.group) == {'target', 'ref'}]
    out = []
    for lg in [c for c in S.columns if c not in ('patient', 'group')]:
        w = S[S.patient.isin(both)].pivot(index='patient', columns='group', values=lg); d = w['target'] - w['ref']
        pw = stats.wilcoxon(w['target'], w['ref']).pvalue
        est, pm = np.nan, np.nan; t = 'C(group, Treatment("ref"))[T.target]'
        for meth in ['lbfgs', 'powell', 'nm', 'bfgs']:
            try:
                m = smf.mixedlm(f'Q("{lg}") ~ C(group, Treatment("ref"))', S, groups=S['patient']).fit(reml=True, method=meth)
                if np.isfinite(m.bse[t]): est, pm = m.params[t], m.pvalues[t]; break
            except Exception: pass
        out.append([lg, len(both), d.mean(), int((d > 0).sum()), pw, est, pm])
    R = pd.DataFrame(out, columns=['ligand', 'patients_paired', 'mean_diff_sms_minus_resident', 'patients_higher_in_sms', 'p_paired_wilcoxon', 'mixed_estimate', 'p_mixed'])
    R['q_paired'] = multipletests(R['p_paired_wilcoxon'], method='fdr_bh')[1]
    R['q_mixed'] = multipletests(R['p_mixed'].fillna(1), method='fdr_bh')[1]
    R = act[['ligand', 'rank', 'aupr']].merge(R, on='ligand').sort_values('rank')
    R.to_csv(f'{O}/ligand_response_patient_test.csv', index=False)


# Section S5C_ligand_activity_plot: OGN+RSPO3+ specificity per ligand (per-patient log2FC, paired Wilcoxon, BH; enriched =
# median log2FC >= 0.5 and q < 0.1), NicheNet ligand supp table and the Supp 5C bar plot (top 25; stars = ligand-response q)
if _run('S5C_ligand_activity_plot'):
    import scanpy as sc, scipy.sparse as sp, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    LFC_MIN = 0.5; Q_MAX = 0.10; NTOP = 25; MIN = 10; PC = 0.1
    act = pd.read_csv(f'{O}/nichenet_ligand_activity_patient.csv').sort_values('rank')
    scrna = sc.read_h5ad(ATLAS, backed='r')
    idx = np.flatnonzero(scrna.obs['cell_type_short'].isin(FIB).values)
    ligs = [g for g in act['ligand'] if g in scrna.var_names]
    f = scrna[idx, ligs].to_memory()
    X = f.X.toarray() if sp.issparse(f.X) else np.asarray(f.X)
    E = np.expm1(X); isogn = (f.obs['cell_type_short'] == 'OGN+RSPO3+ Fib').values; pat = f.obs['Patient'].astype(str).values
    rows = []
    for p in np.unique(pat):
        m = pat == p; a = m & isogn; b = m & ~isogn
        if a.sum() >= MIN and b.sum() >= MIN: rows.append((p, E[a].mean(0), E[b].mean(0)))
    OG = np.vstack([r[1] for r in rows]); OT = np.vstack([r[2] for r in rows]); LFC = np.log2((OG + PC) / (OT + PC))
    res = []
    for j, g in enumerate(ligs):
        d = OG[:, j] - OT[:, j]
        p = stats.wilcoxon(OG[:, j], OT[:, j]).pvalue if np.any(d != 0) else 1.0
        res.append([g, np.median(LFC[:, j]), int((LFC[:, j] > 0).sum()), len(rows), p])
    EN = pd.DataFrame(res, columns=['ligand', 'median_log2FC_ogn_vs_other_fib', 'n_patients_higher_in_ogn', 'n_patients', 'p'])
    EN['q'] = multipletests(EN['p'], method='fdr_bh')[1]; EN['enriched'] = (EN.median_log2FC_ogn_vs_other_fib >= LFC_MIN) & (EN.q < Q_MAX)
    EN.to_csv(f'{O}/ligand_fib_enrichment.csv', index=False)
    act = act.drop(columns=[c for c in ['fib_enriched', 'ogn_mean', 'otherfib_mean'] if c in act]).merge(EN, on='ligand', how='left')
    LR = pd.read_csv(f'{O}/ligand_response_patient_test.csv')[['ligand', 'patients_paired', 'patients_higher_in_sms', 'q_paired']]
    act = act.merge(LR, on='ligand', how='left')
    top = act.head(NTOP).copy(); n = 100
    pd.DataFrame({'Rank': top['rank'], 'Ligand': top['ligand'], 'AUPR': top['aupr'].map(sig3), 'AUROC': top['auroc'].map(sig3),
      f'Top 10 in patient resamples (% of {n})': top['top10'].round(0).astype(int), f'Top 5 in patient resamples (% of {n})': top['top5'].round(0).astype(int),
      'Median rank in resamples': top['median_rank'].map(lambda v: f'{v:g}'),
      'Patients with higher ligand response in SMS mac': top['patients_higher_in_sms'].astype(int).astype(str) + ' of ' + top['patients_paired'].astype(int).astype(str),
      'Ligand response FDR q (paired Wilcoxon)': top['q_paired'].map(pfmt),
      'log2FC, OGN⁺RSPO3⁺ vs other fibroblasts (median over patients)': top['median_log2FC_ogn_vs_other_fib'].map(sig3),
      'Enrichment FDR q': top['q'].map(pfmt), 'Enriched in OGN⁺RSPO3⁺ fibroblasts': np.where(top['enriched'], 'Yes', 'No'),
      'Receptors on SELENOP⁺MMP9⁺SPP1⁺ macrophages': top['mac_receptors']}).to_csv(f'{O}/supp_table_nichenet_ligands.csv', index=False)
    # figure
    BLACK = '#1b1b1f'; TEAL = '#1C7C74'; GRAY = '#B9C2C7'
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'axes.labelsize': 15, 'axes.titleweight': 'normal',
        'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 300, 'text.color': BLACK,
        'axes.labelcolor': BLACK, 'xtick.color': BLACK, 'ytick.color': BLACK, 'axes.edgecolor': BLACK})
    t = top.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(5.6, 0.24 * NTOP + 2.0))
    ax.barh(range(len(t)), t['aupr'], color=[TEAL if e else GRAY for e in t['enriched']], height=0.5)
    ax.set_yticks(range(len(t))); ax.set_yticklabels(t['ligand'], fontsize=13)
    star = lambda q: '***' if q < 0.01 else '**' if q < 0.05 else '*' if q < 0.1 else ''
    xa = t['aupr'].max() * 1.04
    for i, (a, au, q) in enumerate(zip(t['aupr'], t['auroc'], t['q_paired'])):
        ax.text(a + 0.002, i, star(q), va='center', ha='left', fontsize=13); ax.text(xa + 0.012, i, f'{au:.2f}', va='center', ha='left', fontsize=12)
    ax.text(xa + 0.012, len(t) - 0.55, 'AUROC', va='bottom', ha='left', fontsize=12)
    ax.set_xlabel('Ligand Activity (AUPR)'); ax.set_ylabel('Ligand')
    ax.set_title('OGN+RSPO3+Fib Ligands Predicting SELENOP+MMP9+SPP1+ Mac Programs', fontsize=12)
    ax.set_xlim(0, t['aupr'].max() * 1.32); ax.set_xticks([0, 0.1, 0.2]); ax.set_ylim(-0.6, len(t) + 0.3)
    ax.legend(handles=[Patch(color='none', label='Fibroblast specificity:'), Patch(color=TEAL, label='OGN RSPO3 enriched'), Patch(color=GRAY, label='Not enriched')],
              loc='upper center', bbox_to_anchor=(0.42, -0.1), ncol=3, frameon=False, fontsize=11, handlelength=1.0, handletextpad=0.4, columnspacing=0.9)
    plt.tight_layout(); fig.savefig(f'{O}/fig_nichenet_ligand_activity_patient.png', bbox_inches='tight'); plt.close(fig)
