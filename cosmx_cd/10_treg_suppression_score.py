#!/usr/bin/env python
"""CD CosMx (WTX biopsies, 6K resections): antigen-presentation Treg LR module and Treg Suppression Score per Treg.

Treg Suppression Score = scRNA-seq late- vs early-pseudotime Treg signature, sign reversed, z-scored within section downstream.
Inputs  objects/treg_pseudotime.h5ad (cd_scrna/12); objects/nichenet_outputs/ (cd_scrna/16, immune_sender_ranking);
        objects/label_transfer/ (09); WTX norm_anno_tregbin.h5ad and lr_liana_i0476_1_nn_db.h5ad (02); raw_counts_strict_tregs.csv (05)
Outputs panels/fig4_and_related/D_shared_inputs/: fig4_ag_treg_module_lrs.csv, treg_D6_signature_genes.csv (Supp Table 10),
        treg_definition_comparison_per_cell.csv (read by 11-17 and treatment_validation/03)
Sections antigen_presentation_module -> treg_suppression_score (RUN_SECTIONS=<id>,... runs a subset)
"""
import os as _os
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# ---- Section antigen_presentation_module: top-25 NicheNet immune-sender -> late-Treg LR pairs that rank in the top 1,000 by
#      Moran's I in the LIANA object of CD_B4 slide 1, ordered by Moran's I rank
if _run('antigen_presentation_module'):
    import warnings
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, h5py
    from anndata.io import read_elem

    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = f'{D}/panels/fig4_and_related'
    WTX = '/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged'
    RANK = f'{D}/objects/nichenet_outputs/nichenet_destab_late_treg_ligand_strength_imm_logpr_best_sender.csv'

    with h5py.File(f'{WTX}/h5ad_obj/lr_liana_i0476_1_nn_db.h5ad', 'r') as h:
        bvar = read_elem(h['var'])
    mor = bvar.sort_values('morans', ascending=False); mor['morans_rank'] = np.arange(len(mor))

    def ag_module(rank_csv, n_top=25, max_rank=1000):
        r = pd.read_csv(rank_csv).sort_values('priority_rank_geom', ascending=False)
        top = list(dict.fromkeys(f'{a}^{b}' for a, b in zip(r['from'], r['to'])))[:n_top]
        keep = [p for p in top if p in mor.index and mor.loc[p, 'morans_rank'] < max_rank]
        return sorted(keep, key=lambda p: mor.loc[p, 'morans_rank'])

    ag_lrs = ag_module(RANK)
    pd.DataFrame([dict(lr_pair=p, morans_rank_CD_B4_slide1=int(mor.loc[p, 'morans_rank']), morans=mor.loc[p, 'morans'])
                  for p in ag_lrs]).to_csv(f'{OUT}/D_shared_inputs/fig4_ag_treg_module_lrs.csv', index=False, float_format='%.3g')


# ---- Section treg_suppression_score: signature = Wilcoxon pseudotime bin 4 vs bin 0 (q < 0.05, |log2FC| >= 0.5, >= 10% detected,
#      on panel, curated-module LR genes excluded; top 50 per side); CosMx score = score_genes(late) - score_genes(early) on
#      log1p(CP10k); thirds (Late/Mid/Early) and median split within section
if _run('treg_suppression_score'):
    import numpy as np, pandas as pd, anndata as ad, scanpy as sc
    from scipy import sparse
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = f'{D}/panels/fig4_and_related'; LT = f'{D}/objects/label_transfer'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'
    WTX = '/share/fsmresfiles/UC/AtoMx/CD_6k_wholetrans/whole_trans/Processed_merged'
    # LR pairs of the curated modules (pro-fibrotic/ECM, pro-inflammatory, immunomodulatory) and the antigen-presentation Treg module
    MODULE_LRS = (['TIMP2^MMP2', 'COL3A1^DDR2', 'COL17A1^ITGB1', 'COL9A1^ITGA2', 'COL27A1^ITGB1', 'TGM2^SDC4'] +
                  ['TNFSF13^TNFRSF1A', 'TNFSF13^TNFRSF14', 'IL6^IL6ST', 'IL1A^IL1RAP', 'CCL2^CCR3', 'CCL19^CCR7', 'CCL25^CCR9',
                   'CCL25^ACKR4', 'CCL11^ACKR4'] +
                  ['LGALS1^PTPRC', 'IL12A^CD28', 'CD58^CD2', 'PTPRC^CD22', 'CD22^PTPRC'] +
                  pd.read_csv(f'{OUT}/D_shared_inputs/fig4_ag_treg_module_lrs.csv')['lr_pair'].tolist())
    LR_GENES = sorted({g for p in MODULE_LRS for g in p.split('^')})

    # ---- signature: bin 4 vs bin 0 Wilcoxon test on the scRNA-seq Tregs (non-control cells, layer log1p)
    t_ = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    f_ = t_[t_.obs['condition_summary'].astype(str) != 'Control'].copy(); f_.X = f_.layers['log1p']
    f_.obs['pt_bin'] = f_.obs['pt_bin'].astype(int).astype(str).astype('category')
    sc.tl.rank_genes_groups(f_, groupby='pt_bin', groups=['4'], reference='0', method='wilcoxon')
    DE = sc.get.rank_genes_groups_df(f_, group='4').set_index('names')
    Xl = f_.X.tocsr() if sparse.issparse(f_.X) else sparse.csr_matrix(f_.X)
    for b_ in ('0', '4'):
        m_ = (f_.obs['pt_bin'] == b_).values
        DE[f'pct_bin{b_}'] = pd.Series(np.asarray((Xl[m_] > 0).mean(0)).ravel(), index=f_.var_names).reindex(DE.index).values
    del t_, f_, Xl
    SIG_ROWS = []

    def d6_sets(panel, cohort, n=50):
        ok = (DE.pvals_adj < 0.05) & DE.index.isin(list(panel)) & ~DE.index.isin(LR_GENES)
        up4 = DE[ok & (DE.logfoldchanges >= 0.5) & (DE.pct_bin4 >= 0.10)].sort_values('scores', ascending=False).head(n)
        up0 = DE[ok & (DE.logfoldchanges <= -0.5) & (DE.pct_bin0 >= 0.10)].sort_values('scores', ascending=True).head(n)
        for side, d_ in (('late (up in bin 4)', up4), ('early (up in bin 0)', up0)):
            for g_, r_ in d_.iterrows():
                SIG_ROWS.append(dict(cohort=cohort, set=side, gene=g_, wilcoxon_score=r_.scores, log2FC_bin4_vs_bin0=r_.logfoldchanges,
                                     q=r_.pvals_adj, pct_bin0=r_.pct_bin0, pct_bin4=r_.pct_bin4))
        return list(up0.index), list(up4.index)


    def d6_labels(Lfull, early_set, late_set, section, min_n=15):
        ad_ = ad.AnnData(Lfull.values.astype(np.float32), obs=pd.DataFrame(index=Lfull.index), var=pd.DataFrame(index=Lfull.columns))
        sc.tl.score_genes(ad_, late_set, score_name='late', random_state=0)
        sc.tl.score_genes(ad_, early_set, score_name='early', random_state=0)
        s_ = pd.Series(ad_.obs['late'].values - ad_.obs['early'].values, index=Lfull.index)
        n_sec = section.map(section.value_counts())
        grp_ = section.where(n_sec >= min_n, '__cohort_small_sections__')
        q1 = s_.groupby(grp_).transform(lambda x: x.quantile(1 / 3)); q2 = s_.groupby(grp_).transform(lambda x: x.quantile(2 / 3))
        tert = pd.Series(np.where(s_ > q2, 'Late', np.where(s_ <= q1, 'Early', 'Mid')), index=s_.index)
        med = s_.groupby(section).transform('median')
        return s_, tert, pd.Series(np.where(s_ > med, 'Late', 'Early'), index=s_.index)

    percell = []
    # ---- resections: strict Tregs (label-transfer query; raw counts written for the transfer)
    RC = f'{R}/combined/raw_counts_strict_tregs.csv'
    nx_ = pd.read_csv(f'{LT}/resection_strict_treg_pt_transfer_labels.csv', index_col=0)   # query metadata: nCount_RNA = library size
    Xall = pd.read_csv(RC, index_col=0)                               # all panel genes (score_genes control pool)
    idx = Xall.index.intersection(nx_.index); Xall = Xall.loc[idx]
    section = pd.Series(['_'.join(i.split('_')[:3 if i.startswith('11808') else 2]) for i in idx], index=idx)   # = patient_sec (11808_A4 run twice: Cos3/Cos4)
    Lall = np.log1p(Xall.div(nx_.loc[idx, 'nCount_RNA'], axis=0) * 1e4); del Xall
    e6, l6 = d6_sets(Lall.columns, 'resection')
    s6, tert, med = d6_labels(Lall, e6, l6, section); del Lall
    percell.append(pd.DataFrame({'cohort': 'resection', 'section': section, 'is_core_treg': True, 'signature_score_D6': s6,
                                 'D6_tertile': tert, 'D6_median': med}))

    # ---- biopsies: all AUCell Tregs of the CosMx WTX object (core Tregs flagged)
    a = ad.read_h5ad(f'{WTX}/h5ad_obj/norm_anno_tregbin.h5ad', backed='r')
    ti = np.where((a.obs['cell_type_general'] == 'Treg').values)[0]
    sub = a[ti].to_memory()                                           # X = raw counts (integers) in this object
    lib = np.asarray(sub.X.sum(1)).ravel()
    cells = sub.obs['cell'].astype(str).values
    nbb = pd.read_csv(f'{LT}/biopsy_alltreg_treg_pt_transfer_labels.csv', index_col=0)   # biopsy Tregs of the label-transfer query
    idxb = pd.Index(cells)[pd.Index(cells).isin(nbb.index)]
    ob_ = sub.obs.set_index(pd.Index(cells))
    secb = pd.Series(ob_.loc[idxb, 'slide'].astype(str).values, idxb)
    Mall = sub.X.toarray() if hasattr(sub.X, 'toarray') else np.asarray(sub.X)  # all WTX genes (raw counts in X)
    Lball = pd.DataFrame(np.log1p(Mall / lib[:, None] * 1e4), index=cells, columns=list(sub.var_names)).loc[idxb]; del Mall
    e6, l6 = d6_sets(Lball.columns, 'biopsy')
    s6b, tertb, medb = d6_labels(Lball, e6, l6, secb); del Lball
    core = (ob_.loc[idxb, 'cell_type_treg'] == 'Core Treg').values
    percell.append(pd.DataFrame({'cohort': 'biopsy', 'section': secb, 'is_core_treg': core, 'signature_score_D6': s6b,
                                 'D6_tertile': tertb, 'D6_median': medb}))

    SG = pd.DataFrame(SIG_ROWS); SG.to_csv(f'{OUT}/D_shared_inputs/treg_D6_signature_genes.csv', index=False, float_format='%.3g')
    pd.concat(percell).to_csv(f'{OUT}/D_shared_inputs/treg_definition_comparison_per_cell.csv', float_format='%.4g')
