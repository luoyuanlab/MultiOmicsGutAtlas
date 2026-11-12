#!/usr/bin/env python
"""CITE-seq / TCR: TRBV usage and repertoire of early- vs late-pseudotime Tregs.

Panels made
  Supp 8I       TRBV usage in early (bins 0-1) vs late (bins 2-4) pseudotime Tregs with a productive TRB chain (chi-square, Fisher)
  Supp Table 9  TCR repertoire by pseudotime group

Inputs
  /share/fsmresfiles/UC/Schipma/cite-seq_combined/ (annotation table and TCR contigs; Treg query from 04_treg_subset_export.ipynb)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer/cite_notstrict_treg_pt_transfer_labels_rowfix.csv
    (pseudotime-bin labels of the CITE-seq Tregs, cosmx_cd/09_treg_label_transfer.R)

Outputs
  /share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness/ext8I_trbv_usage.png, .csv (Supp 8I)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/D_shared_inputs/supp_table9_tcr_repertoire.csv (Supp Table 9)
    plus supp_table9_tcr_repertoire_variants.csv / supp_table9_tcr_tests_variants.csv (all metrics and tests) in the same folder
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer/cite_treg_tcr_per_cell_rowfix.csv (per-cell TCR summary)

Run order of sections
  S8I_treg_tcr_repertoire

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
# Section S8I_treg_tcr_repertoire: TRBV usage and TCR repertoire of early- vs late-pseudotime CITE-seq Tregs
# Panels: Supp Fig 8I, Supp Table 9
# ======================================================================================================================
if _run('S8I_treg_tcr_repertoire'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness/07_extfig8L_supp_table9_tcr.py'   # working folder of this step (relative paths resolve from here)
    """TRBV usage (Supp Fig 8I) and TCR repertoire table (Supp Table 9) of early- vs late-pseudotime Tregs, CITE-seq.

Cells: T/NK cells with RNA_annotation == Protein_annotation and Protein_annotation == 'tnk' that have TCR contigs; the
  AUCell 'Treg' subset (745 cells = the label-transfer query, 04_treg_subset_export.ipynb); high-confidence contigs.
Labels: pseudotime bins transferred from the CD scRNA-seq Treg reference (cosmx_cd/09_treg_label_transfer.R);
  Early = bins 0-1, Late = bins 2-4 (Methods).
Per cell: has_productive, has_TRA/TRB (productive), pairing = productive TRA and TRB, best contig per chain
  (productive, then UMIs, reads).
Clonotype identity: 10x raw_clonotype_id is numbered per lane (12 lanes, '-1-0'..'-1-11'), so clonotypes are keyed by
  lane + raw_clonotype_id (diversity, top clone); clonal sharing = lane-aware TRB CDR3 nucleotide match.
Metrics (cells with >= 1 productive contig): unique clonotypes; Shannon entropy of clonotypes, also with the larger
  group downsampled to the smaller one over 1,000 draws (seed 42); top-clone fraction and size (cells); % cells in
  expanded clonotypes (>= 2 cells); clonal sharing; mean TRB CDR3 length (cells with a TRB).
Tests: Fisher's exact test for paired fraction and expanded-cell fraction; Mann-Whitney U for TRB CDR3 length; global
  chi-square on all TRBV genes. TRBV panel: the 10 most frequent TRBV genes, per-gene Fisher's exact test, BH-FDR.
Outputs: ext8I_trbv_usage.png/.csv (Supp 8I); D_shared_inputs/supp_table9_tcr_repertoire.csv (Supp Table 9; full
  metrics and tests in supp_table9_tcr_repertoire_variants.csv and supp_table9_tcr_tests_variants.csv);
  objects/label_transfer/cite_treg_tcr_per_cell_rowfix.csv (per-cell TCR summary).
Run with OPENBLAS_CORETYPE=SkylakeX /data/ydn4687/geneformer_env/bin/python."""
    import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
    from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUT = f'{D}/panels/fig4_and_related'; OUT8 = '/share/fsmresfiles/UC/revision_jenny/SuppFig8A-I_treg_pseudotime_robustness'
    C = '/share/fsmresfiles/UC/Schipma/cite-seq_combined'; T = f'{C}/annotation_for_cosmx/third_round'

    # ---- quality filter, T/NK cells with TCR contigs
    anno = pd.read_csv(f'{T}/annotation_table.csv')
    good = set(anno.loc[anno.RNA_annotation == anno.Protein_annotation, 'Unnamed: 0'])
    tcr = pd.read_csv(f'{C}/TCR_merged.csv')
    tcr_good = tcr[tcr.barcode.isin(good)].merge(anno, left_on='barcode', right_on='Unnamed: 0', how='left')
    tcr_good = tcr_good[tcr_good.Protein_annotation == 'tnk']

    # ---- Treg set = label-transfer query (745 AUCell Tregs); Early = bins 0-1, Late = bins 2-4
    new = pd.read_csv(f'{D}/objects/label_transfer/cite_notstrict_treg_pt_transfer_labels_rowfix.csv', index_col=0)
    assert len(new) == 745
    LABELS = {'early_bins0-1_vs_late_bins2-4': np.where(new['predicted.id'] <= 1, 'Early', 'Late')}
    LABELS = {k: pd.Series(v, index=new.index) for k, v in LABELS.items()}
    tcr_treg0 = tcr_good[tcr_good.barcode.isin(new.index)].copy()
    tcr_treg0['lane'] = tcr_treg0.barcode.str.extract(r'(-\d+-\d+)$')[0]
    print('Treg contigs', len(tcr_treg0), '| Treg cells with contigs', tcr_treg0.barcode.nunique())


    def pick_best_contig(x):
        x = x.copy(); x['productive_bool'] = x['productive'] == True
        return x.sort_values(['productive_bool', 'umis', 'reads'], ascending=[False, False, False]).iloc[0]


    def summarize_cell(x):
        has_TRA = ((x['productive'] == True) & (x['chain'] == 'TRA')).any(); has_TRB = ((x['productive'] == True) & (x['chain'] == 'TRB')).any()
        tra = pick_best_contig(x[x.chain == 'TRA']) if (x.chain == 'TRA').any() else None
        trb = pick_best_contig(x[x.chain == 'TRB']) if (x.chain == 'TRB').any() else None
        fn = lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan
        return pd.Series({'has_productive': (x['productive'] == True).any(), 'has_TRA': has_TRA, 'has_TRB': has_TRB,
                          'has_pairing': has_TRA & has_TRB, 'raw_clonotype_id': fn(x['raw_clonotype_id']), 'lane': x['lane'].iloc[0],
                          'dataset': x['dataset'].iloc[0],
                          'TRB_v_gene': trb['v_gene'] if trb is not None else np.nan, 'TRB_j_gene': trb['j_gene'] if trb is not None else np.nan,
                          'TRB_cdr3': trb['cdr3'] if trb is not None else np.nan,
                          'TRB_cdr3_nt': trb['cdr3_nt'] if (trb is not None and trb['productive'] == True) else np.nan,
                          'TRA_v_gene': tra['v_gene'] if tra is not None else np.nan, 'TRA_cdr3': tra['cdr3'] if tra is not None else np.nan})


    df = tcr_treg0[tcr_treg0.high_confidence == True]
    cells = df.groupby('barcode', sort=False).apply(summarize_cell)
    cells['clonotype_lane_safe'] = cells.lane + '|' + cells.raw_clonotype_id.astype(str)
    cells['clone_trb_cdr3nt_lane'] = np.where(cells.TRB_cdr3_nt.notna(), cells.lane + '|' + cells.TRB_cdr3_nt.astype(str), np.nan)
    for k, v in LABELS.items(): cells[k] = v.reindex(cells.index).values
    cells.to_csv(f'{D}/objects/label_transfer/cite_treg_tcr_per_cell_rowfix.csv')


    def shannon(x):
        p = x.value_counts(); p = p / p.sum(); return float(-(p * np.log(p)).sum())


    def metrics(cells, lab):
        allc = cells.copy(); allc['g'] = allc[lab]
        prod = allc[allc.has_productive == True].copy()
        E, L = prod[prod.g == 'Early'], prod[prod.g == 'Late']
        out = {}
        for key, col in [('lane_safe', 'clonotype_lane_safe'), ('raw_id_collides', 'raw_clonotype_id')]:
            n_small = min(len(E), len(L)); rng = np.random.seed(42)
            big = E if len(E) >= len(L) else L
            boot = [shannon(big.sample(n_small, replace=False)[col]) for _ in range(1000)]
            for g, x in [('Early', E), ('Late', L)]:
                vc = x[col].value_counts()
                out[(g, f'unique_clonotypes_{key}')] = len(vc)
                out[(g, f'shannon_{key}')] = shannon(x[col])
                out[(g, f'shannon_downsampled_{key}')] = np.mean(boot) if x is big else shannon(x[col])
                out[(g, f'shannon_downsampled_sd_{key}')] = np.std(boot) if x is big else 0.0
                out[(g, f'top_clone_fraction_{key}')] = vc.max() / vc.sum()
                out[(g, f'top_clone_size_cells_{key}')] = int(vc.max())
                out[(g, f'pct_cells_in_expanded_clones_{key}')] = 100 * (x[col].map(vc) >= 2).mean()
            out[('Late', f'clonal_sharing_pct_late_cells_in_early_clonotypes_{key}')] = 100 * L[col].isin(set(E[col])).mean()
            out[('Late', f'n_late_cells_shared_{key}')] = int(L[col].isin(set(E[col])).sum())
        # lane-aware TRB CDR3nt sharing
        Lt, Et = L[L.clone_trb_cdr3nt_lane.notna()], E[E.clone_trb_cdr3nt_lane.notna()]
        sh = Lt.clone_trb_cdr3nt_lane.isin(set(Et.clone_trb_cdr3nt_lane))
        out[('Late', 'clonal_sharing_pct_late_cells_TRB_CDR3nt_lane_aware')] = 100 * sh.mean()
        out[('Late', 'n_late_cells_shared_TRB_CDR3nt_lane_aware')] = f'{int(sh.sum())}/{len(Lt)}'
        sE = Et.clone_trb_cdr3nt_lane.isin(set(Lt.clone_trb_cdr3nt_lane))
        out[('Early', 'clonal_sharing_pct_early_cells_TRB_CDR3nt_lane_aware')] = 100 * sE.mean()
        for g, x, a in [('Early', E, allc[allc.g == 'Early']), ('Late', L, allc[allc.g == 'Late'])]:
            out[(g, 'n_tregs_query')] = int((cells[lab] == g).sum())
            out[(g, 'n_tregs_with_tcr_contigs')] = len(a)
            out[(g, 'n_tregs_productive_tcr')] = len(x)
            out[(g, 'paired_TRA_TRB_fraction')] = a.has_pairing.astype(bool).mean()
            out[(g, 'mean_TRB_CDR3_aa_length')] = x.TRB_cdr3.dropna().str.len().mean()
            out[(g, 'n_datasets')] = x.dataset.nunique()
        r = pd.Series(out).unstack(0)
        # tests
        aE, aL = allc[allc.g == 'Early'], allc[allc.g == 'Late']
        pr = fisher_exact([[aE.has_pairing.astype(bool).sum(), (~aE.has_pairing.astype(bool)).sum()],
                           [aL.has_pairing.astype(bool).sum(), (~aL.has_pairing.astype(bool)).sum()]])[1]
        vcE, vcL = E.clonotype_lane_safe.map(E.clonotype_lane_safe.value_counts()) >= 2, L.clonotype_lane_safe.map(L.clonotype_lane_safe.value_counts()) >= 2
        pe = fisher_exact([[vcE.sum(), (~vcE).sum()], [vcL.sum(), (~vcL).sum()]])[1]
        pc = mannwhitneyu(E.TRB_cdr3.dropna().str.len(), L.TRB_cdr3.dropna().str.len()).pvalue
        tv = pd.crosstab(prod.TRB_v_gene, prod.g); chi = chi2_contingency(tv)
        tests = {'fisher_p_paired_fraction': pr, 'fisher_p_expanded_cells_lane_safe': pe, 'mannwhitney_p_cdr3_length': pc,
                 'chi2_TRBV': chi[0], 'chi2_TRBV_df': chi[2], 'chi2_TRBV_p': chi[1]}
        return r, tests, prod


    tables, tests_all = [], []
    for lab in LABELS:
        r, t, prod = metrics(cells, lab)
        r = r.reset_index().rename(columns={'index': 'metric'}); r.insert(0, 'labels', lab); tables.append(r)
        tests_all.append(dict(labels=lab, **t))
        if lab == 'early_bins0-1_vs_late_bins2-4': prod_main = prod
    var = pd.concat(tables); var.to_csv(f'{OUT}/D_shared_inputs/supp_table9_tcr_repertoire_variants.csv', index=False, float_format='%.3g')
    tests_all = pd.DataFrame(tests_all); tests_all.to_csv(f'{OUT}/D_shared_inputs/supp_table9_tcr_tests_variants.csv', index=False, float_format='%.3g')
    pd.set_option('display.width', 220); print(var.to_string(index=False)); print(tests_all.to_string(index=False))

    # ---- Supp Table 9
    m = var[var.labels == 'early_bins0-1_vs_late_bins2-4'].set_index('metric')
    t0 = tests_all.iloc[0]
    rows = [('# Tregs (label-transfer query)', 'n_tregs_query', ''), ('# Tregs with productive TCR', 'n_tregs_productive_tcr', ''),
            ('# Unique clonotypes (lane-safe)', 'unique_clonotypes_lane_safe', ''),
            ('Paired TRA-TRB fraction', 'paired_TRA_TRB_fraction', f"Fisher p = {t0.fisher_p_paired_fraction:.3g}"),
            ('Mean TRB CDR3 length (aa)', 'mean_TRB_CDR3_aa_length', f"Mann-Whitney p = {t0.mannwhitney_p_cdr3_length:.3g}"),
            ('Shannon diversity (observed)', 'shannon_lane_safe', ''),
            ('Shannon diversity (larger group downsampled, 1,000 draws)', 'shannon_downsampled_lane_safe', ''),
            ('Top clone fraction', 'top_clone_fraction_lane_safe', ''), ('Top clone size (cells)', 'top_clone_size_cells_lane_safe', ''),
            ('% cells in expanded clonotypes (>= 2 cells)', 'pct_cells_in_expanded_clones_lane_safe', f"Fisher p = {t0.fisher_p_expanded_cells_lane_safe:.3g}"),
            ('Clonal sharing: % late cells whose TRB CDR3 (lane-aware) is also in early', 'clonal_sharing_pct_late_cells_TRB_CDR3nt_lane_aware', ''),
            ('TRBV usage', None, f"chi-square = {t0.chi2_TRBV:.3g}, df = {int(t0.chi2_TRBV_df)}, p = {t0.chi2_TRBV_p:.3g}")]
    st9 = pd.DataFrame([dict(Metric=a, **{'Early-pseudotime Treg': m.loc[b, 'Early'] if b else '', 'Late-pseudotime Treg': m.loc[b, 'Late'] if b else ''}, Test=c)
                        for a, b, c in rows])
    st9.loc[st9.Metric.str.startswith('Shannon diversity (larger'), 'Early-pseudotime Treg'] = m.loc['shannon_downsampled_lane_safe', 'Early']
    shared = m.loc['n_late_cells_shared_TRB_CDR3nt_lane_aware', 'Late']
    st9.loc[st9.Metric.str.startswith('Clonal sharing'), 'Early-pseudotime Treg'] = ''
    st9.loc[st9.Metric.str.startswith('Clonal sharing'), 'Test'] = f'{shared} late cells'


    def fmt(v):
        try: f = float(v); return str(int(f)) if f.is_integer() else f'{f:.3g}'
        except (TypeError, ValueError): return v


    for c in ['Early-pseudotime Treg', 'Late-pseudotime Treg']: st9[c] = st9[c].map(fmt)
    st9.to_csv(f'{OUT}/D_shared_inputs/supp_table9_tcr_repertoire.csv', index=False)
    print(st9.to_string(index=False))

    # ---- Supp 8I TRBV usage
    prod = prod_main[prod_main.TRB_v_gene.notna()]   # cells with a TRB V gene (denominator of the percentages)
    n = prod.g.value_counts()
    freq = pd.crosstab(prod.TRB_v_gene, prod.g, normalize='columns') * 100; cnt = pd.crosstab(prod.TRB_v_gene, prod.g)
    top = cnt.sum(1).sort_values(ascending=False).head(10).index
    pv = [fisher_exact([[cnt.loc[g, 'Early'], n['Early'] - cnt.loc[g, 'Early']], [cnt.loc[g, 'Late'], n['Late'] - cnt.loc[g, 'Late']]])[1] for g in cnt.index]
    fdr = pd.Series(multipletests(pv, method='fdr_bh')[1], index=cnt.index)
    src = pd.DataFrame({'n_early': cnt.Early, 'n_late': cnt.Late, 'pct_early': freq.Early, 'pct_late': freq.Late,
                        'fisher_p': pd.Series(pv, index=cnt.index), 'fdr_bh': fdr}).sort_values('n_early', ascending=False)
    src['in_panel_top10'] = src.index.isin(top)
    src.to_csv(f'{OUT8}/ext8I_trbv_usage.csv', float_format='%.3g')
    top = sorted(top, key=lambda g: (int(''.join(c for c in g.split('-')[0] if c.isdigit()) or 0), g))
    long = freq.loc[top].reset_index().melt(id_vars='TRB_v_gene', var_name='pt_bin', value_name='pct')
    long['pt_bin'] = long.pt_bin.map({'Early': f"Early (n = {n['Early']})", 'Late': f"Late (n = {n['Late']})"})
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    sns.barplot(data=long, x='TRB_v_gene', y='pct', hue='pt_bin', palette=['#dfc6c7', '#885e7c'], ax=ax)
    for i, g in enumerate(top):
        if fdr[g] < 0.05: ax.text(i, freq.loc[g].max() + 0.3, '*', ha='center', fontsize=14)
    ax.set_xticks(range(len(top))); ax.set_xticklabels(top, rotation=90)
    ax.set_ylabel('Percentage of Cells', fontsize=11); ax.set_xlabel('TRBV Gene', fontsize=11)
    ax.set_title(f'TRBV Usage in Tregs\nAll TRBV Genes χ² p = {t0.chi2_TRBV_p:.2g}, No Gene FDR < 0.05' if (fdr >= 0.05).all()
                 else f'TRBV Usage in Tregs\nAll TRBV Genes χ² p = {t0.chi2_TRBV_p:.2g}, * FDR < 0.05', fontsize=11, loc='center')
    ax.legend(title='Pseudotime Group', fontsize=9, title_fontsize=9, frameon=True, loc='upper left', bbox_to_anchor=(1.01, 1))
    sns.despine(ax=ax); plt.tight_layout(); fig.savefig(f'{OUT8}/ext8I_trbv_usage.png', dpi=300, bbox_inches='tight')
    print(src.head(12).round(3).to_string())
