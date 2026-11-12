#!/usr/bin/env python
"""TAURUS CD scRNA-seq (anti-TNF): senders of the ligands of the four Treg activation programs.

Panels made
  Supp 10  per-ligand sender contributions to the Remission - Non-remission difference in change

Inputs
  /share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny/raw_norm_counts_anno.h5ad.gz
  /share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/cache/cn10_cytokine_cells.parquet and cn10_cycling_lineage_check.csv (04_taurus_cytokine_sources.py)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/cellchatdb_human_interactions_snapshot.csv

Outputs
  /share/fsmresfiles/UC/revision_jenny/SuppFig10_treg_activation_ligand_senders/ (CSV, cache/, panel and legend)

Run order of sections
  04_taurus_cytokine_sources.py (cytokine_celltype_screen, F6B_coarse_screen) -> S10_ligand_sources -> S10_ligand_senders

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist. Command-line arguments and environment variables used by a section are listed in its own
  docstring.
  Paths are the lab's absolute paths (to be replaced by placeholders before public release).
  Environment: see the README section "Environment" (OPENBLAS_CORETYPE=SkylakeX on Intel Sapphire Rapids nodes;
  OPENBLAS_NUM_THREADS=32 for identical model output).
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
# Section S10_ligand_sources: senders of the Treg activation program ligands, all-cell models and decomposition of the change
# Panels: Supp Fig 10 (models)
# ======================================================================================================================
if _run('S10_ligand_sources'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig10_treg_activation_ligand_senders/63_treg_activation_ligand_sources.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq (numbers only, no figure; Supp Fig 10 source data): which cell groups SEND the ligands of the four
Treg activation programs scored in Tregs in Fig 6A, and do those sources persist in non-remitters?
Programs -> ligand families: IL-2/STAT5 = IL2; TNFR2 non-canonical NF-kB = TNF, LTA; co-stimulatory/suppressive receptor module = CD80, CD86,
  ICOSLG, TNFSF4 (OX40L), TNFSF18 (GITRL), CD70, CD274 (PD-L1), PDCD1LG2 (PD-L2), CD40LG, IL10; TNFRSF9 / 4-1BB = TNFSF9 (4-1BBL).
As F6B_coarse_screen (04_taurus_cytokine_sources.py): cell loading (raw_norm_counts_anno.h5ad.gz, X = log1p(CP10k), layers['counts'] = raw
  UMIs), coarse-group mapping (COARSE; Cycling -> Epithelial by the cached data check; Neuron = no group), eligibility rule, R model
  y ~ post * remission + (1 | patient) + (1 | patient:timepoint) (lme4 + lmerTest REML, Satterthwaite df, y standardised, bobyqa retry), the
  four contrasts, and the abundance x expression decomposition (two-factor midpoint split, observed biopsy values).
Cache: the cell cache of 04_taurus_cytokine_sources.py (Fig6_treg_antitnf_taurus/supplementary/cache/cn10_cytokine_cells.parquet) is read
  only; the ligand columns (log-norm + raw counts) are built once with the same loader into cache/treg_activation_ligand_cells.parquet (this
  folder) and joined on the cell index (asserted equal).
Part 1 (sources): per ligand x coarse group, % of pre-treatment transcripts (raw UMI sum; 13 groups + Neuron = 100%), % of pre-treatment group
  cells expressing (> 0), mean log-norm; plus the top 3 fine types (aucell_cell_type_short) by transcript share.
  Main source group(s) = transcript sources (>= 20% of a ligand's pre-treatment transcripts, at least the top group) plus differential
  drivers (groups carrying >= 25% of the gross differential change in part 3: |group contribution| / sum over groups of |contribution|).
Part 2 (all-cell model): per biopsy, mean log-norm over all its cells; model above; BH within the ligand family per contrast.
  Also, for each main source group: within-group expression (biopsies with >= 10 cells of the group) and the group's abundance
  (% of all biopsy cells), same model; BH within the ligand family per contrast (expression) / across the tested groups (abundance).
Part 3 (decomposition, as F6B_coarse_screen): all-cell change per remission group = sum over groups of the change in mean(fraction x within-group mean);
  each group's change and the differential change (remission - non-remission) split into abundance and expression parts.
Output: treg_activation_ligand_sources.csv (this folder; 3 significant figures). No figure.
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time, subprocess, tempfile, shutil
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'; R = '/share/fsmresfiles/UC/revision_jenny'
    F6S = f'{R}/Fig6_treg_antitnf_taurus/supplementary'; HERE = f'{R}/SuppFig10_treg_activation_ligand_senders'   # cell cache of 04_taurus_cytokine_sources.py in F6S
    OUT = f'{HERE}/treg_activation_ligand_sources.csv'; CACHE = f'{F6S}/cache/cn10_cytokine_cells.parquet'
    CACHE2 = f'{HERE}/cache/treg_activation_ligand_cells.parquet'; CYC = f'{F6S}/cache/cn10_cycling_lineage_check.csv'
    KEYS = ['Patient', 'Remission_status', 'Treatment']; MIN_N = 10; MIN_PRE_PAT = 3; MAIN_SHARE = 20; DRV_PCT = 25
    FAM = {'IL-2/STAT5': ['IL2'], 'TNFR2 non-canonical NF-kB': ['TNF', 'LTA'],
           'Co-stimulatory/suppressive receptor module': ['CD80', 'CD86', 'ICOSLG', 'TNFSF4', 'TNFSF18', 'CD70', 'CD274', 'PDCD1LG2', 'CD40LG', 'IL10'],
           'TNFRSF9 / 4-1BB': ['TNFSF9']}
    ALIAS = {'TNFSF4': 'OX40L', 'TNFSF18': 'GITRL', 'CD274': 'PD-L1', 'PDCD1LG2': 'PD-L2', 'TNFSF9': '4-1BBL', 'CD40LG': 'CD40L'}
    GENES = [g for v in FAM.values() for g in v]; G2F = {g: f for f, v in FAM.items() for g in v}
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500); pd.set_option('display.max_columns', 50)
    assert os.path.exists(CACHE), f'{CACHE} missing: run 04_taurus_cytokine_sources.py (cytokine_celltype_screen) first'

    # ================================================================ cells (cytokine cell cache) + ligand columns (same loader, own cache)
    C = pd.read_parquet(CACHE)[['sample_id'] + KEYS + ['aucell_cell_type_short', 'aucell_cell_category_from_type']]; log('cells', C.shape)
    if not os.path.exists(CACHE2):
        import scanpy as sc
        ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); log('loaded', ad.shape)
        assert ad.raw is None and np.allclose(np.expm1(ad.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)        # X = log1p(CP10k)
        miss = [g for g in GENES if g not in ad.var_names]; assert not miss, f'ligands not in var_names: {miss}'
        assert ad.obs_names.is_unique and 'counts' in ad.layers
        gi = [ad.var_names.get_loc(g) for g in GENES]
        cnt = ad.layers['counts'][:, gi].toarray(); lg = ad.X[:, gi].toarray()
        assert np.allclose(cnt.astype(float), np.round(cnt)) and ((cnt > 0) == (lg > 0)).all()          # raw integer UMIs, same zeros
        L = pd.DataFrame(index=ad.obs_names.copy())
        for j, g in enumerate(GENES): L[g] = lg[:, j]; L[f'{g}_counts'] = cnt[:, j]
        del ad
        with open(CACHE2, 'wb') as fh: L.to_parquet(fh)
        log('cached', CACHE2, L.shape)
    L = pd.read_parquet(CACHE2)
    assert L.index.equals(C.index), 'ligand cache and cytokine cell cache cell order differ'
    C = pd.concat([C, L[[c for g in GENES for c in (g, f'{g}_counts')]].astype(np.float64)], axis=1)   # float64: exact decomposition sums
    C['biopsy'] = C.sample_id
    assert C.groupby('biopsy').Patient.nunique().max() == 1 and C.groupby('biopsy').Treatment.nunique().max() == 1
    assert C.groupby('Patient').Remission_status.nunique().max() == 1

    # ================================================================ coarse mapping (as F6B_coarse_screen)
    FIBS = ['ADAMDEC1+ Activ Fib', 'ADAMDEC1+ Fib', 'Activ Fib', 'GREM2+ Myofib', 'HHIP+ Myofib', 'OGN+RSPO3+ Fib', 'PLA2G2A+ ECM Fib',
            'SELENOP+ Fib', 'T-interact Fib', 'VSTM2A+ Crypt Top Fib']
    EPIS = ['Absorp Progenitor (Col)', 'Absorp Progenitor (SI)', 'BEST4+ cell', 'Enterocyte', 'Enteroendocrine', 'Goblet', 'Imm Enterocyte',
            'M cell', 'Paneth', 'Stem', 'Tuft']
    COARSE = {'DC': ['cDC1', 'cDC2', 'CCR7+ DC', 'Trans cDC2/Mac'],
              'Macrophage': ['Mac S+XS-', 'Mac S+M+S+', 'Mac S+SG+', 'Mac S+M+P+', 'CD9+ Mac'],
              'Monocyte-derived Mac': ['Mo-Mac', 'Inf Mo-Mac'],
              'Granulocyte/Neutrophil': ['Granulocyte', 'Neutrophil'],
              'B cell': ['B cell'], 'Plasma': ['Plasma'],
              'CD4 T': ['CD4 Effector', 'CD4 Trm(Prolif)', 'CD4 Trm+D42', 'Exhaus CD4 Tfh', 'Tfh', 'Treg'],
              'CD8 T': ['CD8 GZMK+ Trm', 'CD8 GZMK- Trm', 'CD8 GZMK- Trm(Prolif)'],
              'Other T/NK': ['Gd T', 'Gd T (Vd2g9+)', 'NK', 'Resting T'],
              'Fibroblast': FIBS, 'Pericyte': ['CD36+ Peri', 'RERGL+ Contr Peri'], 'Endothelial': ['Endothelial'],
              'Epithelial': EPIS + ['Cycling']}
    NOGRP = 'Neuron (no coarse group)'; ALL = 'All cells'
    CY = pd.read_csv(CYC).set_index('item').value
    assert max(CY['% of Cycling cells expressing CD3D'], CY['% of Cycling cells expressing CD3E']) < 5 and CY['% of Cycling cells expressing EPCAM'] > 50
    T2G = {t: g for g, ts in COARSE.items() for t in ts}
    types = set(C.aucell_cell_type_short)
    assert set(T2G) <= types and types - set(T2G) == {'Neuron'}
    C['grp'] = C.aucell_cell_type_short.map(T2G).fillna(NOGRP)
    GROUPS = list(COARSE)
    MEMBERS = {g: '; '.join(v) for g, v in COARSE.items()}; MEMBERS[ALL] = 'every cell'; MEMBERS[NOGRP] = 'Neuron'
    GRP = [('rem_pre', 'Remission', 'Pre'), ('rem_post', 'Remission', 'Post'), ('nonrem_pre', 'Non_Remission', 'Pre'), ('nonrem_post', 'Non_Remission', 'Post')]
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}
    BM = C.groupby('biopsy')[KEYS].first().assign(n_all=C.groupby('biopsy').size())
    assert BM.n_all.min() >= MIN_N
    log('biopsies', len(BM), 'patients', BM.Patient.nunique())

    # ================================================================ part 1: sources
    P1 = []
    for g in GENES:
        for name in GROUPS + [NOGRP]:
            m = pre & (C.grp == name).values
            P1.append(dict(part='1 source', program=G2F[g], feature=g, cell_group=name, share_of_transcripts_pre=100 * C.loc[m, f'{g}_counts'].sum() / TOT_PRE[g],
                           pct_cells_expressing_pre=100 * (C.loc[m, g] > 0).mean(), mean_expression_pre=C.loc[m, g].mean(),
                           pct_of_all_cells_pre=100 * m.sum() / pre.sum(), n_cells_pre=int(m.sum()), total_umi_pre=TOT_PRE[g]))
        P1.append(dict(part='1 source', program=G2F[g], feature=g, cell_group=ALL, share_of_transcripts_pre=100.0, pct_cells_expressing_pre=100 * (C.loc[pre, g] > 0).mean(),
                       mean_expression_pre=C.loc[pre, g].mean(), pct_of_all_cells_pre=100.0, n_cells_pre=int(pre.sum()), total_umi_pre=TOT_PRE[g]))
    P1 = pd.DataFrame(P1)
    assert np.allclose(P1[P1.cell_group != ALL].groupby('feature').share_of_transcripts_pre.sum(), 100)
    P1['source_rank'] = P1[P1.cell_group.isin(GROUPS)].groupby('feature').share_of_transcripts_pre.rank(ascending=False, method='first')
    P1['main_source'] = P1.cell_group.isin(GROUPS) & ((P1.share_of_transcripts_pre >= MAIN_SHARE) | (P1.source_rank == 1))
    MAIN = {g: list(P1[(P1.feature == g) & P1.main_source].sort_values('source_rank').cell_group) for g in GENES}
    # top fine types
    FT = []
    Cp = C[pre]
    for g in GENES:
        s = Cp.groupby('aucell_cell_type_short')[f'{g}_counts'].sum() / TOT_PRE[g] * 100; s = s.sort_values(ascending=False).head(3)
        pe = Cp.groupby('aucell_cell_type_short')[g].apply(lambda x: 100 * (x > 0).mean())
        for i, (t, v) in enumerate(s.items()):
            FT.append(dict(part='1b top fine types', program=G2F[g], feature=g, cell_group=T2G.get(t, NOGRP), fine_type=t, source_rank=i + 1,
                           share_of_transcripts_pre=v, pct_cells_expressing_pre=pe[t], n_cells_pre=int((Cp.aucell_cell_type_short == t).sum())))
    FT = pd.DataFrame(FT)
    for g in GENES:
        s = P1[(P1.feature == g) & P1.cell_group.isin(GROUPS)].sort_values('share_of_transcripts_pre', ascending=False).head(4)
        f = FT[FT.feature == g]
        print(f'{g:9s} UMI {TOT_PRE[g]:8.0f} expr {P1[(P1.feature == g) & (P1.cell_group == ALL)].pct_cells_expressing_pre.iloc[0]:.2f}% | '
              + '; '.join(f'{r.cell_group} {r.share_of_transcripts_pre:.0f}% ({r.pct_cells_expressing_pre:.1f}% expr)' for r in s.itertuples())
              + ' | fine: ' + '; '.join(f'{r.fine_type} {r.share_of_transcripts_pre:.0f}%' for r in f.itertuples()), flush=True)

    # ================================================================ part 3: decomposition (as F6B_coarse_screen)
    DEC = []
    cnt = C.groupby(['biopsy', 'grp']).size().unstack(fill_value=0).reindex(BM.index)
    frac = cnt.div(BM.n_all, axis=0)
    for f in GENES:
        mg = C.groupby(['biopsy', 'grp'])[f].mean().unstack().reindex(BM.index)
        contrib = (frac * mg).fillna(0)
        allm = C.groupby('biopsy')[f].mean().reindex(BM.index)
        assert np.allclose(contrib.sum(1), allm, atol=1e-10)
        for rs, rlab in [('Remission', 'remission'), ('Non_Remission', 'non-remission')]:
            b0 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Pre')]; b1 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Post')]
            tot0, tot1 = allm[b0].mean(), allm[b1].mean(); dtot = tot1 - tot0
            for g in GROUPS + [NOGRP]:
                C0, C1 = contrib.loc[b0, g].mean(), contrib.loc[b1, g].mean(); F0, F1 = frac.loc[b0, g].mean(), frac.loc[b1, g].mean()
                W0 = C0 / F0 if F0 > 0 else 0.0; W1 = C1 / F1 if F1 > 0 else 0.0
                dab = (F1 - F0) * (W0 + W1) / 2; dex = (W1 - W0) * (F0 + F1) / 2
                assert abs(dab + dex - (C1 - C0)) < 1e-12
                DEC.append(dict(feature=f, program=G2F[f], cell_group=g, remission_group=rlab, contribution_pre=C0, contribution_post=C1, contribution_to_change=C1 - C0,
                                pct_of_all_cell_change=100 * (C1 - C0) / dtot if dtot != 0 else np.nan, abundance_part=dab, expression_part=dex,
                                mean_fraction_pre=100 * F0, mean_fraction_post=100 * F1, weighted_group_mean_pre=W0, weighted_group_mean_post=W1,
                                all_cell_mean_pre=tot0, all_cell_mean_post=tot1, all_cell_change=dtot, n_biopsies_pre=len(b0), n_biopsies_post=len(b1)))
    DEC = pd.DataFrame(DEC)
    chk = DEC.groupby(['feature', 'remission_group']).agg(s=('contribution_to_change', 'sum'), t=('all_cell_change', 'first'))
    assert np.allclose(chk.s, chk.t, atol=1e-12), chk
    W = DEC.pivot_table(index=['feature', 'program', 'cell_group'], columns='remission_group',
                        values=['contribution_to_change', 'abundance_part', 'expression_part', 'all_cell_change'])
    DD = pd.DataFrame({k: W[k]['remission'] - W[k]['non-remission'] for k in ['contribution_to_change', 'abundance_part', 'expression_part', 'all_cell_change']}).reset_index()
    DD['pct_of_all_cell_change'] = 100 * DD.contribution_to_change / DD.all_cell_change
    DD['remission_group'] = 'differential (remission − non-remission)'
    DEC = pd.concat([DEC, DD], ignore_index=True)
    DEC['part'] = '3 decomposition'
    # differential drivers: groups carrying >= DRV_PCT% of the gross (sum of |group contributions|) differential change of a ligand
    DIFF = DEC[DEC.remission_group.str.startswith('differential') & DEC.cell_group.isin(GROUPS)].copy()
    DIFF['pct_of_gross_differential'] = 100 * DIFF.contribution_to_change.abs() / DIFF.groupby('feature').contribution_to_change.transform(lambda x: x.abs().sum())
    DRIVER = {g: list(DIFF[(DIFF.feature == g) & (DIFF.pct_of_gross_differential >= DRV_PCT)].sort_values('pct_of_gross_differential', ascending=False).cell_group) for g in GENES}
    TSRC = {g: list(MAIN[g]) for g in GENES}
    MAIN = {g: TSRC[g] + [x for x in DRIVER[g] if x not in TSRC[g]] for g in GENES}
    ROLE = lambda g, cg: '; '.join(r for r, ok in [('transcript source', cg in TSRC[g]), ('differential driver', cg in DRIVER[g])] if ok)
    for g in GENES: print(f'{g:9s} transcript sources {TSRC[g]}; differential drivers {DRIVER[g]}')
    # ================================================================ part 2: model input (all cells; main-source expression + abundance)
    def pseudobulk(D, col):
        return D.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'),
                                       Treatment=('Treatment', 'first'), n_cells=(col, 'size'), y=(col, 'mean')).reset_index()
    GINFO = {}
    for name in GROUPS + [ALL]:
        m = (C.grp == name).values if name != ALL else np.ones(len(C), bool)
        nb = C[m].groupby('biopsy').size().reindex(BM.index, fill_value=0); keep = BM[nb >= MIN_N]
        np_ = {nm: keep[(keep.Remission_status == r_) & (keep.Treatment == t_)].Patient.nunique() for nm, r_, t_ in GRP}
        nbx = {nm: int(((keep.Remission_status == r_) & (keep.Treatment == t_)).sum()) for nm, r_, t_ in GRP}
        GINFO[name] = dict(nb=nb, keep=keep, eligible=np_['rem_pre'] >= MIN_PRE_PAT and np_['nonrem_pre'] >= MIN_PRE_PAT and nbx['rem_post'] >= 1 and nbx['nonrem_post'] >= 1)
    LONG, DESC = [], []
    def add(key, PB, meta):
        d = dict(meta)
        for nm, r_, t_ in GRP:
            s = PB[(PB.Remission_status == r_) & (PB.Treatment == t_)].y
            d[f'mean_of_biopsy_means_{nm}'] = s.mean(); d[f'sem_of_biopsy_means_{nm}'] = s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else np.nan
            d[f'n_biopsies_model_{nm}'] = len(s); d[f'n_patients_model_{nm}'] = PB[(PB.Remission_status == r_) & (PB.Treatment == t_)].Patient.nunique()
        d['estimable'] = bool(PB.y.std() > 0)
        d['feature_key'] = key; DESC.append(d)
        if d['estimable']:
            LONG.append(PB.assign(feature_key=key, pt=PB.Patient + '_' + PB.Treatment, post=(PB.Treatment == 'Post').astype(float),
                                  rem=(PB.Remission_status == 'Remission').astype(float))[['feature_key', 'Patient', 'pt', 'post', 'rem', 'y']])
    for g in GENES:
        add(f'{ALL}||{g}||expr', pseudobulk(C, g), dict(part='2 all-cell model', program=G2F[g], feature=g, cell_group=ALL, measure='all-cell biopsy mean expression',
            unit='log-normalised expression, log1p(CP10k)', biopsies_used='all biopsies'))
        for name in MAIN[g]:
            gi_ = GINFO[name]
            if not gi_['eligible']: continue
            D = C[(C.grp == name).values & C.biopsy.isin(gi_['keep'].index).values]
            add(f'{name}||{g}||expr', pseudobulk(D, g), dict(part='2b main-source within-group model', program=G2F[g], feature=g, cell_group=name,
                measure='within-group biopsy mean expression', unit='log-normalised expression, log1p(CP10k)', biopsies_used=f'biopsies with >= {MIN_N} cells of the group'))
    ABG = sorted({n for g in GENES for n in MAIN[g] if GINFO[n]['eligible']}, key=GROUPS.index)
    for name in ABG:
        PB = BM.assign(y=100 * GINFO[name]['nb'] / BM.n_all).reset_index()[['biopsy'] + KEYS + ['y']]
        add(f'{name}||Abundance||abund', PB, dict(part='2c main-source abundance model', program='(all)', feature='Abundance', cell_group=name,
            measure='% of all biopsy cells in the group', unit='% of all cells in the biopsy (estimate in percentage points)', biopsies_used='all biopsies'))
    DESC = pd.DataFrame(DESC); LONG = pd.concat(LONG, ignore_index=True)
    log('models to fit', LONG.feature_key.nunique(), '; not estimable:', DESC[~DESC.estimable].feature_key.tolist())

    # ================================================================ lme4 / lmerTest (as F6B_coarse_screen)
    RCODE = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); D <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); out <- list()
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1),
           baseline_difference = c(rem = 1))
for (cc in unique(D$feature_key)) {
  d <- D[D$feature_key == cc, ]; d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
  y <- d$y; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  opt <- "nloptwrap (default)"
  m <- tryCatch(lmer(ys ~ post * rem + (1 | Patient) + (1 | pt), data = d, REML = TRUE), error = function(e) NULL)
  if (is.null(m)) { opt <- "bobyqa (default optimizer failed)"
    m <- tryCatch(lmer(ys ~ post * rem + (1 | Patient) + (1 | pt), data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
                  error = function(e) NULL) }
  if (is.null(m)) { for (k in names(Ls)) out[[length(out) + 1]] <- data.frame(feature_key = cc, contrast_key = k, optimizer = "fit failed"); next }
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m)); stopifnot(length(fe) == 4)
  for (k in names(Ls)) {
    L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
    ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
    out[[length(out) + 1]] <- data.frame(feature_key = cc, contrast_key = k, estimate = ct$Estimate * sy, estimate_sd = ct$Estimate,
      se = ct$`Std. Error` * sy, df = ct$df, statistic = ct$`t value`, p = ct$`Pr(>|t|)`, var_patient = vc$vcov[vc$grp == "Patient"] * sy^2,
      var_patient_timepoint = vc$vcov[vc$grp == "pt"] * sy^2, var_residual = vc$vcov[vc$grp == "Residual"] * sy^2, singular = isSingular(m),
      converged = !any(grepl("converge|gradient|Hessian|eigenvalue", c(m@optinfo$conv$lme4$messages, m@optinfo$warnings))), optimizer = opt)
  }
}
nm <- unique(unlist(lapply(out, names)))
out <- lapply(out, function(x) { for (n in setdiff(nm, names(x))) x[[n]] <- NA; x[nm] })
write.csv(do.call(rbind, out), a[2], row.names = FALSE)
"""
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        LONG.to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        log(f'fitting {LONG.feature_key.nunique()} models in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/r_out.csv'], check=True)
        RR = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(RR) == 4 * LONG.feature_key.nunique() and (RR.p.notna() | (RR.optimizer == 'fit failed')).all()
    log('optimizer use:', RR.drop_duplicates('feature_key').optimizer.value_counts().to_dict())
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)', 'baseline_difference': 'Baseline difference (remission − non-remission, pre)'}
    RR['contrast'] = RR.contrast_key.map(CON)
    tq = stats.t.ppf(0.975, RR.df); RR['ci_low'] = RR.estimate - tq * RR.se; RR['ci_high'] = RR.estimate + tq * RR.se
    RR = DESC.merge(RR, on='feature_key', how='left')
    RR['q_bh'] = np.nan; RR['bh_family'] = ''
    for (pt_, pr, k), i in RR[RR.p.notna()].groupby(['part', 'program', 'contrast_key']).groups.items():
        RR.loc[i, 'q_bh'] = multipletests(RR.loc[i, 'p'], method='fdr_bh')[1]
        RR.loc[i, 'bh_family'] = (f'{pr}: ligands in the family, this contrast' if pt_ == '2 all-cell model' else
                                  f'{pr}: ligand x main-source tests in the family, this contrast' if pt_.startswith('2b') else 'main-source groups, this contrast')
    RR['n_tests_in_bh_family'] = RR.groupby(['part', 'program', 'contrast_key']).q_bh.transform('count')
    for pt_ in RR.part.unique():
        print(f'\n== {pt_}')
        S = RR[RR.part == pt_].pivot_table(index=['program', 'cell_group', 'feature'], columns='contrast_key', values=['estimate', 'p', 'q_bh'], sort=False)
        for (pr, cg, f), r in S.iterrows():
            print(f'{f:9s} {cg:22s} ' + ' | '.join(f"{k[:12]} {r[('estimate', k)]:+.3g} P {r[('p', k)]:.3g} q {r[('q_bh', k)]:.3g}" for k in CON))

    for f in GENES:   # observed all-cell means reproduce the model input
        a_ = RR[(RR.part == '2 all-cell model') & (RR.feature == f)].iloc[0]
        for rl, nm in [('remission', 'rem'), ('non-remission', 'nonrem')]:
            d_ = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].iloc[0]
            assert np.isclose(d_.all_cell_mean_pre, a_[f'mean_of_biopsy_means_{nm}_pre']) and np.isclose(d_.all_cell_mean_post, a_[f'mean_of_biopsy_means_{nm}_post'])
    DEC['main_source'] = [cg in MAIN[f] for f, cg in zip(DEC.feature, DEC.cell_group)]
    DEC['source_role'] = [ROLE(f, cg) if cg in GROUPS else '' for f, cg in zip(DEC.feature, DEC.cell_group)]
    DEC = DEC.merge(DIFF[['feature', 'cell_group', 'pct_of_gross_differential']], on=['feature', 'cell_group'], how='left')
    DEC['pct_abundance_of_group_change'] = 100 * DEC.abundance_part / DEC.contribution_to_change
    DEC['unit'] = 'log-normalised expression, log1p(CP10k), all-cell biopsy mean scale; mean_fraction in % of all cells'
    DEC['method'] = ('observed means over biopsies (no model): contribution = mean over biopsies of fraction × within-group mean; '
                 'abundance / expression parts = two-factor midpoint split, add up to the contribution')
    print('\n== decomposition, main sources (transcript sources + differential drivers)')
    for f in GENES:
        for cg in MAIN[f]:
            s = DEC[(DEC.feature == f) & (DEC.cell_group == cg)].set_index('remission_group')
            print(f'{f:9s} {cg:22s} ' + ' | '.join(f"{rl[:12]} dC {s.loc[rl, 'contribution_to_change']:+.3g} ({s.loc[rl, 'pct_of_all_cell_change']:.0f}% of all) ab {s.loc[rl, 'abundance_part']:+.3g} ex {s.loc[rl, 'expression_part']:+.3g}"
                                                for rl in s.index) + f" | frac rem {s.loc['remission', 'mean_fraction_pre']:.1f}->{s.loc['remission', 'mean_fraction_post']:.1f}, nonrem {s.loc['non-remission', 'mean_fraction_pre']:.1f}->{s.loc['non-remission', 'mean_fraction_post']:.1f}")

    # ================================================================ CSV
    P1['main_source'] = P1.main_source.astype('boolean')
    RR['source_role'] = [ROLE(f, cg) if cg in GROUPS and f in GENES else '' for f, cg in zip(RR.feature, RR.cell_group)]
    P1['source_role'] = [ROLE(f, cg) if cg in GROUPS else '' for f, cg in zip(P1.feature, P1.cell_group)]
    P1['main_source'] = [cg in MAIN[f] for f, cg in zip(P1.feature, P1.cell_group)]
    ALLR = pd.concat([P1, FT, RR.drop(columns=['feature_key']), DEC], ignore_index=True)
    ALLR['ligand_alias'] = ALLR.feature.map(ALIAS).fillna('')
    ALLR['mapping'] = ALLR.cell_group.map(MEMBERS)
    ALLR.loc[ALLR.part.str.startswith('2'), 'model_formula'] = 'biopsy value ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df; no covariates'
    ALLR['main_source_rule'] = (f'main source = transcript source (coarse group with >= {MAIN_SHARE}% of the ligand pre-treatment transcripts, raw UMI, or the top group) '
                            f'or differential driver (group carrying >= {DRV_PCT}% of the gross differential change: |group contribution| / sum of |group contributions|)')
    COLS = ['part', 'program', 'feature', 'ligand_alias', 'cell_group', 'fine_type', 'mapping', 'main_source', 'source_role', 'source_rank', 'share_of_transcripts_pre',
            'pct_cells_expressing_pre', 'mean_expression_pre', 'pct_of_all_cells_pre', 'n_cells_pre', 'total_umi_pre', 'measure', 'contrast', 'estimate', 'ci_low',
            'ci_high', 'se', 'estimate_sd', 'df', 'statistic', 'p', 'q_bh', 'bh_family', 'n_tests_in_bh_family', 'estimable'] \
           + [f'{a}_{nm}' for a in ['mean_of_biopsy_means', 'sem_of_biopsy_means', 'n_biopsies_model', 'n_patients_model'] for nm, _, _ in GRP] \
           + ['remission_group', 'contribution_pre', 'contribution_post', 'contribution_to_change', 'pct_of_all_cell_change', 'abundance_part', 'expression_part',
              'pct_abundance_of_group_change', 'pct_of_gross_differential', 'mean_fraction_pre', 'mean_fraction_post', 'weighted_group_mean_pre', 'weighted_group_mean_post', 'all_cell_mean_pre',
              'all_cell_mean_post', 'all_cell_change', 'n_biopsies_pre', 'n_biopsies_post', 'biopsies_used', 'var_patient', 'var_patient_timepoint', 'var_residual',
              'singular', 'converged', 'optimizer', 'unit', 'model_formula', 'main_source_rule', 'method']
    assert set(ALLR.columns) <= set(COLS) | {'contrast_key'}, set(ALLR.columns) - set(COLS)
    RES = ALLR.reindex(columns=COLS)
    with open(OUT, 'w') as fh:
        RES.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, RES.shape)


# ======================================================================================================================
# Section S10_ligand_senders: dot plot of the Treg activation program ligand senders (figure and legend; plot only)
# Panels: Supp Fig 10
# Note: refuses to overwrite an existing PNG / legend; pass --check <dir> to write them to another folder
# ======================================================================================================================
if _run('S10_ligand_senders'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig10_treg_activation_ligand_senders/64_SuppFig10_ligand_senders.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq, Supplemental Figure 10: cell groups that send the cognate ligands of the four Fig 6A Treg activation
programs (IL-2/STAT5, TNFR2 non-canonical NF-kB, co-stimulatory receptor module, TNFRSF4 / OX40), and which group carries the
remission - non-remission difference in change. ~11 x 6 in; y: ligand -> Treg receptor [Fig 6A program]; x: sender cell group (45-degree
headers); P strip: P (biopsy-level mixed model); colour bar and dot-size legend stacked on the right.
Inputs (read; no model is re-fitted):
  treg_activation_ligand_sources.csv (S10_ligand_sources, this folder): part '1 source' (share of pre-treatment transcripts per coarse group),
    part '3 decomposition' (per-group change in the all-cell mean, per remission group), part '2 all-cell model' (change / difference P),
    part '2c main-source abundance model' (only for the biopsy / patient counts in the legend). TNFSF9 (4-1BBL; 4-1BB is not a Fig 6A
    program) is in the CSV but is dropped here. TNFSF4 is tested under the co-stimulatory module; it is tagged here as the OX40
    program ligand (Fig 6A). No BH q value is shown (only per-ligand all-cell P), so no q is recomputed.
  cd_treg_rebuild_rowfix/objects/cellchatdb_human_interactions_snapshot.csv (the CellChatDB snapshot used to build the co-stimulatory
    receptor module): every receptor written next to a ligand is asserted to be a CellChatDB receptor of that ligand.
Dot plot. Rows = the 13 ligands of the four Fig 6A programs. Columns = the coarse groups of F6B_coarse_screen, myeloid first: DC, Macrophage, Mo-Mac, B cell, Plasma,
  CD4 T, CD8 T, Stromal (= Endothelial + Fibroblast (all subtypes) + Pericyte; summed; the decomposition is additive across groups), and Other
  (= Granulocyte/Neutrophil, Other T/NK, Epithelial, Neuron; summed) only if any ligand has > 5% of its pre-treatment transcripts there.
  Dot size = % of the ligand's pre-treatment raw UMI from the group.
  Dot colour = group contribution to the difference in change = 100 * (dC_rem / pre_rem - dC_nonrem / pre_nonrem), where dC = change
  (post - pre) in the group's mean contribution (fraction x within-group mean) to the all-cell biopsy mean and pre = that remission group's
  pre-treatment all-cell mean (the % scale of Fig 6B).
  Rows ordered by the plotted sender column with the largest |dot colour| (column order), then by difference P; groups of rows separated by
  thin light-grey lines (no labels). After each ligand: receptor(s) on Tregs and the Fig 6A program tag (grey).
  Right strip: all-cell change P in remission, in non-remission and the difference P (part 2 all-cell model), 2 significant figures,
  bold if P < 0.05.
Layout checks (asserted): ligand / receptor labels do not overlap each other or the dots; column and strip headers do not overlap; P values do
  not overlap; dots never touch a neighbouring row / column; legend items do not overlap each other or the P strip.
Outputs: SuppFig10_treg_activation_ligand_senders.png, SuppFig10_treg_activation_ligand_senders_legend.txt (this folder; never overwritten).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file> [--check <outdir>]"""
    import os, sys, warnings
    warnings.filterwarnings('ignore')
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import re
    import numpy as np, pandas as pd, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize
    from decimal import Decimal, ROUND_HALF_UP

    R = '/share/fsmresfiles/UC/revision_jenny'; F6S = f'{R}/Fig6_treg_antitnf_taurus/supplementary'; HERE = f'{R}/SuppFig10_treg_activation_ligand_senders'
    CSV = f'{HERE}/treg_activation_ligand_sources.csv'
    CCDB = f'{R}/cd_treg_rebuild_rowfix/objects/cellchatdb_human_interactions_snapshot.csv'
    STEM = f'{HERE}/SuppFig10_treg_activation_ligand_senders'
    if '--check' in sys.argv: STEM = f"{sys.argv[sys.argv.index('--check') + 1]}/SuppFig10_treg_activation_ligand_senders"; os.makedirs(os.path.dirname(STEM), exist_ok=True)
    for f in [f'{STEM}.png', f'{STEM}_legend.txt']: assert not os.path.exists(f), f'refusing to overwrite {f}'

    # Fig 6A programs: (tag, legend name, ligands, program name of the ligands in the S10_ligand_sources CSV)
    COSTIM_CSV = 'Co-stimulatory/suppressive receptor module'
    PROG = [('IL-2/STAT5', 'IL-2/STAT5', ['IL2'], 'IL-2/STAT5'),
            ('TNFR2', 'TNFR2 non-canonical NF-κB', ['TNF', 'LTA'], 'TNFR2 non-canonical NF-kB'),
            ('co-stim', 'co-stimulatory receptor module', ['CD80', 'CD86', 'CD274', 'PDCD1LG2', 'IL10', 'TNFSF18', 'CD40LG', 'CD70', 'ICOSLG'], COSTIM_CSV),
            ('OX40', 'TNFRSF4 (OX40)', ['TNFSF4'], COSTIM_CSV)]
    DROPPED = ['TNFSF9']                                                                # 4-1BBL: TNFRSF9 / 4-1BB is not a Fig 6A program
    TAG = {g: t for t, _, gs, _ in PROG for g in gs}
    ALIAS = {'TNFSF4': 'OX40L', 'TNFSF18': 'GITRL', 'CD274': 'PD-L1', 'PDCD1LG2': 'PD-L2', 'CD40LG': 'CD40L'}
    # receptor(s) on Tregs (program receptors; the co-stimulatory pathways are the CellChatDB pathways of the receptor module)
    REC = {'IL2': ['IL2RA', 'IL2RB', 'IL2RG'], 'TNF': ['TNFRSF1B'], 'LTA': ['TNFRSF1B'], 'CD80': ['CD28', 'CTLA4'], 'CD86': ['CD28', 'CTLA4'],
           'CD274': ['PDCD1'], 'PDCD1LG2': ['PDCD1'], 'IL10': ['IL10RA', 'IL10RB'], 'TNFSF18': ['TNFRSF18'], 'CD40LG': ['CD40'],
           'TNFSF4': ['TNFRSF4'], 'CD70': ['CD27'], 'ICOSLG': ['ICOS']}
    GENES = [g for _, _, gs, _ in PROG for g in gs]
    assert set(REC) == set(GENES) and len(GENES) == 13 and not set(DROPPED) & set(GENES)
    def rlab(rs):   # IL2RA/IL2RB/IL2RG -> IL2RA/B/G when the receptors differ only by the last letter
        cp = os.path.commonprefix(rs)
        return rs[0] + '/' + '/'.join(r[len(cp):] for r in rs[1:]) if len(rs) > 1 and cp and all(len(r) - len(cp) == 1 for r in rs) else '/'.join(rs)
    assert rlab(REC['IL2']) == 'IL2RA/B/G' and rlab(REC['CD80']) == 'CD28/CTLA4' and rlab(REC['IL10']) == 'IL10RA/B'
    COLS = ['DC', 'Macrophage', 'Monocyte-derived Mac', 'B cell', 'Plasma', 'CD4 T', 'CD8 T', 'Stromal', 'Other']
    CLAB = {'Macrophage': 'Mac', 'Monocyte-derived Mac': 'Mo-Mac'}                      # short column headers
    STROMAL = ['Endothelial', 'Fibroblast', 'Pericyte']
    OTHER = ['Granulocyte/Neutrophil', 'Other T/NK', 'Epithelial', 'Neuron (no coarse group)']
    MYELOID = 3
    OTHER_MIN = 5.0                                                                     # % of a ligand's transcripts needed to plot the Other column

    # ================================================================ CSV (S10_ligand_sources)
    d = pd.read_csv(CSV)
    assert set(d.feature.dropna()) - {'Abundance'} == set(GENES) | set(DROPPED), 'S10_ligand_sources tested a different ligand set'
    for tag, _, gs, prog in PROG: assert set(d[d.feature.isin(gs)].program) == {prog}, tag
    d = d[~d.feature.isin(DROPPED)].copy()                                               # TNFSF9 out of every table below
    def col(s): return s.where(~s.isin(STROMAL), 'Stromal').where(~s.isin(OTHER), 'Other')
    S = d[(d.part == '1 source') & (d.cell_group != 'All cells')].copy(); S['col'] = col(S.cell_group)
    assert set(S.cell_group) == set(COLS[:-2]) | set(STROMAL) | set(OTHER)
    SHARE = S.groupby(['feature', 'col']).share_of_transcripts_pre.sum().unstack().reindex(index=GENES, columns=COLS)
    assert np.allclose(SHARE.sum(1), 100, atol=0.5), SHARE.sum(1)
    OTH_SHOW = bool((SHARE['Other'] > OTHER_MIN).any())
    OTH_MAX = SHARE['Other'].max(); OTH_MAX_G = SHARE['Other'].idxmax()
    OTH_GRP_MAX = S[S.cell_group.isin(OTHER)].share_of_transcripts_pre.max()
    print(f"Other (granulocyte/neutrophil, other T/NK, epithelial, neuron): max {OTH_MAX:.2f}% ({OTH_MAX_G}); plotted: {OTH_SHOW}", flush=True)
    PCOLS = COLS if OTH_SHOW else COLS[:-1]
    DC_ = d[d.part == '3 decomposition'].copy(); DC_['col'] = col(DC_.cell_group)
    def contrib(rg):
        x = DC_[DC_.remission_group == rg]
        pre = x.groupby('feature').all_cell_mean_pre.first(); chg = x.groupby('feature').all_cell_change.first()
        c = x.groupby(['feature', 'col']).contribution_to_change.sum().unstack().reindex(index=GENES, columns=COLS)
        assert np.allclose(c.sum(1), chg.reindex(GENES), rtol=2e-2, atol=1e-6), rg
        return c, pre.reindex(GENES), chg.reindex(GENES)
    cR, preR, chR = contrib('remission'); cN, preN, chN = contrib('non-remission')
    VAL = 100 * (cR.div(preR, axis=0) - cN.div(preN, axis=0))
    NET = 100 * (chR / preR - chN / preN)
    assert np.allclose(VAL.sum(1), NET, atol=1.0), (VAL.sum(1) - NET)
    A = d[d.part == '2 all-cell model'].pivot_table(index='feature', columns='contrast', values='p').reindex(GENES)
    PC = {'rem': 'Change in remission (post − pre)', 'non': 'Change in non-remission (post − pre)', 'diff': 'Differential change (remission − non-remission)'}
    assert A[list(PC.values())].notna().all().all()
    AB = d[d.part == '2c main-source abundance model'].set_index(['cell_group', 'contrast'])
    NB = {nm: int(AB.iloc[0][f'n_biopsies_model_{nm}']) for nm in ['rem_pre', 'rem_post', 'nonrem_pre', 'nonrem_post']}
    NP = {nm: int(AB.iloc[0][f'n_patients_model_{nm}']) for nm in ['rem_pre', 'rem_post', 'nonrem_pre', 'nonrem_post']}
    for k_ in AB.index: assert all(int(AB.loc[k_][f'n_biopsies_model_{nm}']) == NB[nm] for nm in NB), k_
    print('biopsies', NB, 'patients', NP, flush=True)
    TOT113 = d[(d.part == '1 source') & (d.feature == 'TNFSF18') & (d.cell_group == 'All cells')].total_umi_pre.iloc[0]

    # ================================================================ row order: sender column with the largest |contribution to the difference|
    TOPC = VAL[PCOLS].abs().idxmax(axis=1)
    if not OTH_SHOW:   # the omitted remainder must not carry a larger contribution than the top plotted column
        gap = VAL['Other'].abs() - VAL[PCOLS].abs().max(1)
        print('omitted Other |contribution| minus top plotted |contribution| (max):', round(gap.max(), 2), flush=True)
        assert (gap < 0).all(), gap[gap >= 0]
    ORD = pd.DataFrame({'top': TOPC, 'ci': TOPC.map(PCOLS.index), 'p': A[PC['diff']]}).rename_axis('g').reset_index()
    ORD = ORD.sort_values(['ci', 'p', 'g'], kind='mergesort')
    GENES = ORD.g.tolist(); GRPS = [c for c in PCOLS if c in set(TOPC)]
    print('\nrow order (top-contributing sender column, then difference P)')
    print(ORD.assign(top_val=[VAL.loc[g, t] for g, t in zip(ORD.g, ORD.top)], tag=ORD.g.map(TAG)).to_string(index=False))

    # ================================================================ receptors: each must be a CellChatDB receptor of its ligand
    cc = pd.read_csv(CCDB)
    spl = lambda s: [t.strip() for t in re.split(r'[,;]', str(s)) if t.strip() and t.strip() != 'nan']
    for g, rs in REC.items():
        x = cc[cc['ligand.symbol'].map(lambda s: g in spl(s))]
        db = {t for r in x['receptor.symbol'] for t in spl(r)}
        assert set(rs) <= db, (g, rs, db)
        extra = sorted(db - set(rs))
        print(f'{g}: shown {"/".join(rs)}' + (f'; other CellChatDB receptors not shown: {", ".join(extra)}' if extra else ''))

    pd.set_option('display.width', 250)
    print('\nshare of pre-treatment transcripts (%)\n', SHARE.reindex(GENES).round(1).to_string())
    print('\ncontribution to the difference in change (% of pre-treatment mean)\n', VAL.reindex(GENES).round(1).to_string())
    print('\nall-cell P\n', A.reindex(GENES)[list(PC.values())].to_string(), flush=True)

    # ================================================================ figure
    def sig2(p):   # 2 significant figures, half-up on the CSV's 3-significant-figure value, trailing zero kept (0.998 -> 1.0)
        d_ = Decimal(repr(float(p))); r_ = d_.quantize(Decimal(1).scaleb(d_.adjusted() - 1), rounding=ROUND_HALF_UP)
        if r_.adjusted() > d_.adjusted(): r_ = r_.quantize(Decimal(1).scaleb(r_.adjusted() - 1))
        return str(r_)
    assert sig2(0.998) == '1.0' and sig2(0.0498) == '0.050' and sig2(0.00521) == '0.0052' and sig2(0.0109) == '0.011'
    fmtv = lambda p: '< 0.001' if p < 0.001 else sig2(p)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.linewidth': 1.0})
    BLUE, RED, K, GR = '#2166ac', '#b2182b', '#111111', '#8c8c8c'
    VMAX = 75.0; assert VAL[PCOLS].abs().max().max() <= VMAX
    CMAP = plt.get_cmap('RdBu'); NORM = Normalize(-VMAX, VMAX)
    # ================================================================ layout: compact ~11 x 6 in
    FS_HD, FS_LIG, FS_REC, FS_P, FS_LEG = 13, 13, 10, 12.5, 11.5                    # column/strip headers, ligand, receptor+tag, P values, legends
    SMAX = 500.0; size = lambda s: SMAX * np.asarray(s) / 100.0                         # marker area (pt^2) at 100%; diameter = sqrt(area)
    GAP = 0.5                                                                           # extra space (rows) between groups of rows (separator line)
    ypos, seps, y = {}, [], 0.0
    for i, g in enumerate(GENES):
        if i and TOPC[g] != TOPC[GENES[i - 1]]: seps.append(y - 0.5 + GAP / 2); y += GAP
        ypos[g] = y; y += 1.0
    YMAX = y - 1.0
    UNIT, COLW, PGAP = 0.32, 0.35, 0.13                                                 # inches per row, per dot column; gap between P columns
    WIN, HIN, DPI = 16.0, 9.0, 300                                                      # canvas; saved with a tight bounding box
    NCOL = len(PCOLS); XLO, XHI = -0.5, NCOL - 0.45
    AX0, BOTA = 4.4, 1.0; AXW = COLW * (XHI - XLO); AXH = UNIT * (YMAX + 1.1); TOPA = BOTA + AXH
    fig = plt.figure(figsize=(WIN, HIN), dpi=DPI)
    axin = lambda x, y_, w, h: fig.add_axes([x / WIN, y_ / HIN, w / WIN, h / HIN])
    axA = axin(AX0, BOTA, AXW, AXH)
    SX0 = AX0 + AXW + 0.05
    axS = axin(SX0, BOTA, 2.0, AXH); axS.sharey(axA)                                   # width reset below from the widest P value per column
    for g in GENES:
        for j, c in enumerate(PCOLS):
            s = SHARE.loc[g, c]
            if s > 0: axA.scatter(j, ypos[g], s=size(s), c=[CMAP(NORM(VAL.loc[g, c]))], edgecolors='#444444', linewidths=0.6, zorder=3)
    axA.set_xlim(XLO, XHI); axA.set_ylim(YMAX + 0.55, -0.55); axA.set_xticks([]); axA.set_yticks([])
    for sp in axA.spines.values(): sp.set_visible(False)
    for j in range(NCOL): axA.axvline(j, color='#eeeeee', lw=0.8, zorder=0)
    for g in GENES: axA.axhline(ypos[g], color='#eeeeee', lw=0.8, zorder=0)
    axA.axvline(MYELOID - 0.5, color=GR, lw=1.0, ls=(0, (3, 3)), zorder=1)
    HOFF = (-3, 4)                                                                      # header anchor offset (pt) from the column top
    HD = [axA.annotate(CLAB.get(c, c), xy=(j, 1), xycoords=('data', 'axes fraction'), xytext=HOFF, textcoords='offset points', rotation=45,
                       ha='left', va='bottom', rotation_mode='anchor', fontsize=FS_HD, color=K, annotation_clip=False) for j, c in enumerate(PCOLS)]
    # ligand (black) and receptor + program tag (grey, smaller, same line) label columns, left of the dots
    yt_ = axA.get_yaxis_transform()
    LIG = {g: axA.text(0, ypos[g], f'{g} ({ALIAS[g]})' if g in ALIAS else g, transform=yt_, ha='left', va='center', fontsize=FS_LIG, color=K) for g in GENES}
    RCP = {g: axA.text(0, ypos[g], f'→ {rlab(REC[g])} [{TAG[g]}]', transform=yt_, ha='left', va='center', fontsize=FS_REC, color='#7a7a7a') for g in GENES}
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); axw = axA.get_window_extent(rend).width
    wl = max(t.get_window_extent(rend).width for t in LIG.values()); wr = max(t.get_window_extent(rend).width for t in RCP.values())
    PT = DPI / 72 / axw                                                                  # 1 pt in axes fraction
    XR = -(wr / axw + 6 * PT); XL = XR - 5 * PT                                         # receptor column (left-aligned), ligand right edge
    for g in GENES: LIG[g].set_x(XL); LIG[g].set_ha('right'); RCP[g].set_x(XR)
    axS.axis('off')
    STRIP = [('rem', 'Rem P', BLUE), ('non', 'Non-rem P', RED), ('diff', 'Diff P', K)]
    PV, PVC = [], []
    for j, (k, hd, cc_) in enumerate(STRIP):
        for g in GENES:
            p = A.loc[g, PC[k]]
            PV.append(axS.text(0, ypos[g], fmtv(p), ha='center', va='center', fontsize=FS_P, color=cc_, fontweight='bold' if p < 0.05 else 'normal')); PVC.append(j)
    fig.canvas.draw()
    PW = [max(PV[i].get_window_extent(rend).width for i in range(len(PV)) if PVC[i] == j) / DPI for j in range(3)]   # widest value per column (in)
    CX = [0.06 + PW[0] / 2]; [CX.append(CX[-1] + PW[j - 1] / 2 + PGAP + PW[j] / 2) for j in (1, 2)]
    SW = CX[2] + PW[2] / 2 + 0.04
    axS.set_position([SX0 / WIN, BOTA / HIN, SW / WIN, AXH / HIN]); axS.set_xlim(0, SW)
    for i, t in enumerate(PV): t.set_x(CX[PVC[i]])
    SH = [axS.annotate(hd, xy=(CX[j], 1), xycoords=('data', 'axes fraction'), xytext=HOFF, textcoords='offset points', rotation=45, ha='left',
                       va='bottom', rotation_mode='anchor', fontsize=FS_HD, color=cc_, annotation_clip=False) for j, (k, hd, cc_) in enumerate(STRIP)]
    for g in GENES: axS.axhline(ypos[g], color='#f3f3f3', lw=0.8, zorder=0)
    XS0 = XL - wl / axw - 4 * PT
    fig.canvas.draw()
    XS1 = axA.transAxes.inverted().transform(axS.transAxes.transform((1, 0)))[0]
    for ys in seps: axA.add_artist(Line2D([XS0, XS1], [ys, ys], transform=yt_, color='#c8c8c8', lw=0.9, clip_on=False, zorder=0))
    # legends: right of the P strip, stacked (colour bar on top, dot-size scale below), starting at the top of the dot rows
    LX = SX0 + SW + 0.2
    CBT = fig.text(LX / WIN, TOPA / HIN, 'Contribution to\ndifference in\nchange (% of\npre-treatment\nmean;\nremission −\nnon-remission)',
                   ha='left', va='top', fontsize=FS_LEG, linespacing=1.05)
    fig.canvas.draw(); cbt_y0 = CBT.get_window_extent(rend).y0 / DPI                     # inches
    CBH = 1.2; cax = axin(LX, cbt_y0 - 0.12 - CBH, 0.16, CBH)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation='vertical', ticks=[-75, -50, -25, 0, 25, 50, 75])
    cb.ax.tick_params(labelsize=FS_LEG, length=3); cb.outline.set_linewidth(0.6)
    fig.canvas.draw(); cax_y0 = cax.get_tightbbox(rend).y0 / DPI
    SZT = fig.text(LX / WIN, (cax_y0 - 0.22) / HIN, '% of pre-\ntreatment\ntranscripts', ha='left', va='top', fontsize=FS_LEG, linespacing=1.05)
    fig.canvas.draw(); szt_y0 = SZT.get_window_extent(rend).y0 / DPI
    SZS = [10, 25, 50, 75]; SSP = 0.28                                                  # inches between size-legend dots
    sax = axin(LX, szt_y0 - 0.08 - SSP * len(SZS), 1.0, SSP * len(SZS)); sax.axis('off'); sax.set_xlim(0, 1.0); sax.set_ylim(SSP * len(SZS), 0)
    SZL = []
    for i, s in enumerate(SZS):
        sax.scatter(0.16, SSP * (i + 0.5), s=size(s), facecolors='white', edgecolors='#444444', linewidths=0.8, clip_on=False)
        SZL.append(sax.text(0.40, SSP * (i + 0.5), f'{s}%', ha='left', va='center', fontsize=FS_LEG))
    # axis titles: y (rotated, left of the ligand labels), x over the sender columns and a P-strip title, both above the 45-degree headers
    FS_AT = 13.5
    YT = axA.text(XS0 - 5 * PT, 0.5, 'Treg Activation Program (Ligand → Treg Receptor)', transform=axA.transAxes, rotation=90, ha='right', va='center', fontsize=FS_AT)
    fig.canvas.draw()
    ytop = max(h.get_window_extent(rend).y1 for h in HD + SH)
    yT = fig.transFigure.inverted().transform((0, ytop + 6 * DPI / 72))[1]
    xT = fig.transFigure.inverted().transform((axA.transData.transform(((NCOL - 1) / 2, 0))[0], 0))[0]
    xP = fig.transFigure.inverted().transform((axS.transData.transform((SW / 2, 0))[0], 0))[0]
    XT = fig.text(xT, yT, 'Sender cell group', ha='center', va='bottom', fontsize=FS_AT)
    PT_ = fig.text(xP, yT, 'P (biopsy-level mixed model)', ha='center', va='bottom', fontsize=FS_AT)

    # layout checks
    fig.canvas.draw()
    bb = lambda a: a.get_window_extent(rend)
    def ov(a, b, pad=1.0): return not (a.x1 + pad <= b.x0 or b.x1 + pad <= a.x0 or a.y1 + pad <= b.y0 or b.y1 + pad <= a.y0)
    LH = 1.2 * FS_HD                                                                     # line height (pt) of the 45-degree headers
    # parallel 45-degree headers: perpendicular distance between neighbouring baselines must exceed one line height
    assert COLW * 72 * np.sin(np.radians(45)) > LH, 'column headers overlap'
    hx = [axA.transData.transform((j, 0))[0] for j in range(NCOL)] + [axS.transData.transform((CX[j], 0))[0] for j in range(3)]
    assert min(np.diff(hx)) / DPI * 72 * np.sin(np.radians(45)) > LH, 'column / strip headers overlap'
    LAB = list(LIG.values()) + list(RCP.values())
    for i, t in enumerate(LAB):
        for t2 in LAB[i + 1:]: assert not ov(bb(t), bb(t2), pad=0.0), (t.get_text(), t2.get_text())
    ax0px = axA.transData.transform((0, 0))[0]
    for t in LAB: assert bb(t).x1 < ax0px - 8, ('label reaches the first dot column', t.get_text())
    for g in GENES: assert bb(LIG[g]).x1 + 4 <= bb(RCP[g]).x0, g
    for i, t in enumerate(PV):                                                           # P values: >= 8 pt clear space between columns
        for t2 in PV[i + 1:]: assert not ov(bb(t), bb(t2), pad=8 * DPI / 72), (t.get_text(), t2.get_text())
    d_pt = lambda s: np.sqrt(size(s)) + 0.6                                              # dot diameter (pt) incl. edge
    rowpt = UNIT * 72; colpt = COLW * 72; DMAX = d_pt(SHARE[PCOLS].max().max())
    for gi, g in enumerate(GENES):                                                       # dots never touch a neighbouring row / column
        for j, c in enumerate(PCOLS):
            if j + 1 < NCOL: assert (d_pt(SHARE.loc[g, c]) + d_pt(SHARE.loc[g, PCOLS[j + 1]])) / 2 < colpt, (g, c)
            if gi + 1 < len(GENES):
                g2 = GENES[gi + 1]; assert (d_pt(SHARE.loc[g, c]) + d_pt(SHARE.loc[g2, c])) / 2 < rowpt * (ypos[g2] - ypos[g]), (g, g2, c)
    assert d_pt(25) - d_pt(10) > 3.5, 'size legend 10% vs 25% not distinguishable'
    xmaxP = max(bb(t).x1 for t in PV); LEGB = [bb(CBT), cax.get_tightbbox(rend), bb(SZT), bb(sax)] + [bb(t) for t in SZL]
    for b1 in LEGB: assert b1.x0 > xmaxP + 4, 'legend overlaps the P strip'
    for a_, b_ in [(0, 1), (1, 2), (2, 3), (0, 2)]: assert not ov(LEGB[a_], LEGB[b_]), ('legend items overlap', a_, b_)
    for h in SH: assert not ov(bb(h), bb(CBT)), 'strip header overlaps the colour-bar title'
    assert not ov(bb(XT), bb(PT_), pad=10), 'x-axis title overlaps the P-strip title'
    for t in [XT, PT_]:
        assert not ov(bb(t), bb(CBT)), 'axis title overlaps the colour-bar title'
        for h in HD + SH: assert not ov(bb(t), bb(h)), ('axis title overlaps a column header', t.get_text(), h.get_text())
    assert bb(YT).x1 < min(bb(t).x0 for t in LIG.values()) - 4, 'y-axis title overlaps the ligand labels'
    assert bb(YT).y0 >= axA.get_window_extent(rend).y0 - 1 and bb(YT).y1 <= axA.get_window_extent(rend).y1 + 1, 'y-axis title longer than the rows'
    assert bb(sax).y0 >= axA.get_window_extent(rend).y0 - 1, 'legend column extends below the dot rows'
    fig.savefig(f'{STEM}.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08, facecolor='white')
    print('wrote', f'{STEM}.png', f'dot max diameter {DMAX:.1f} pt, row {rowpt:.1f} pt, column {colpt:.1f} pt', flush=True)

    # ================================================================ legend
    LCOL = {'DC': 'DC', 'Macrophage': 'macrophage', 'Monocyte-derived Mac': 'Mo-Mac', 'B cell': 'B cell', 'Plasma': 'plasma',
            'CD4 T': 'CD4 T', 'CD8 T': 'CD8 T', 'Stromal': 'stromal', 'Other': 'other'}
    T4 = S[S.feature == 'TNFSF4'].set_index('cell_group').share_of_transcripts_pre
    TSH4 = T4[['CD4 T', 'CD8 T', 'Other T/NK']].sum(); print('TNFSF4 T/NK share', T4[['CD4 T', 'CD8 T', 'Other T/NK']].round(1).to_dict(), round(TSH4, 1))
    assert TSH4 > 50 and TOPC['TNFSF4'] == 'CD4 T', TSH4
    r0 = lambda v: str(Decimal(repr(round(float(v), 6))).quantize(Decimal(1), rounding=ROUND_HALF_UP))   # half-up (36.5 -> 37)
    assert r0(36.5) == '37' and r0(10.1) == '10'
    prog_txt = '; '.join(f'{nm} [{t}]' for t, nm, _, _ in PROG)
    LEG = (f"Supplemental Figure 10. Cell sources of the ligands of the Fig 6A Treg activation programs in TAURUS Crohn's disease scRNA-seq before and "
       f"after anti-TNF (Remission {NP['rem_pre']} patients, {NB['rem_pre']} pre-/{NB['rem_post']} post-treatment biopsies; Non-remission "
       f"{NP['nonrem_pre']} patients pre ({NB['nonrem_pre']} biopsies), {NP['nonrem_post']} post ({NB['nonrem_post']} biopsies)). "
       f"Rows, the {len(GENES)} cognate ligands of the four Treg Activation Programs of Fig 6A ({prog_txt}); grey text after each ligand, its receptor(s) on Tregs "
       f"and program. Transcripts of TNFSF4 (OX40L), the OX40 ligand, come mostly from T cells ({r0(TSH4)}% of pre-treatment transcripts: CD4 T ({r0(T4['CD4 T'])}%), CD8 T ({r0(T4['CD8 T'])}%) and other T/NK cells ({r0(T4['Other T/NK'])}%)). "
       f"Rows ordered by the sender group contributing most to the remission − non-remission difference (dot colour of largest magnitude; "
       f"groups separated by grey lines, in column order), then by difference P. "
       f"Dot size, % of the ligand's pre-treatment transcripts from each cell group; dot colour, the group's contribution to the difference in "
       f"change (remission − non-remission), as % of each remission group's pre-treatment all-cell mean; red, relatively higher in non-remission "
       f"after treatment. Myeloid groups left of the dashed line; Stromal, endothelial cells, fibroblasts and pericytes. "
           + ("Other, granulocytes/neutrophils, other T/NK cells (γδ T, NK, resting T), epithelial cells and neurons. " if OTH_SHOW else
              f"Remaining groups (granulocytes/neutrophils, other T/NK, epithelial cells, neurons) are omitted (together < 5% of each "
          f"ligand's pre-treatment transcripts; at most {OTH_MAX:.1f}%). ")
           + f"TNFSF18 had {int(TOT113)} pre-treatment transcripts. "
       f"Right, all-cell change P within each remission group and difference P (bold, P < 0.05); biopsy-level mixed model, value ~ timepoint × "
       f"remission + (1 | patient) + (1 | patient:timepoint).")
    nw = len(LEG.split()); print('legend words', nw); assert 150 <= nw <= 280, nw
    for bad in ['S6', 'Supplementary Figure', 'Nominal', 'territor', 'blunt', 'SELENOP', 'D6', 'quiescent', 'TNFSF9', 'TNFRSF9', '4-1BB', '-derived']: assert bad not in LEG, bad
    with open(f'{STEM}_legend.txt', 'w') as fh: fh.write(LEG + '\n')
    print('wrote', f'{STEM}_legend.txt', flush=True)
    print('\ntop-contributing sender per ligand:', {g: CLAB.get(TOPC[g], TOPC[g]) for g in GENES})
