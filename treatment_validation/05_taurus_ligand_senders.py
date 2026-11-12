#!/usr/bin/env python
"""TAURUS CD scRNA-seq (anti-TNF): senders of the ligands of the four Treg activation programs.

Panels: Supp 10 (per-ligand sender contributions to the Remission - Non-remission difference in change).
Inputs: TAURUS raw_norm_counts_anno.h5ad.gz; cell cache of 04 (Fig6_treg_antitnf_taurus/supplementary/cache/cn10_cytokine_cells.parquet).
Outputs (SuppFig10_treg_activation_ligand_senders/): cache/treg_activation_ligand_cells.parquet, treg_activation_ligand_sources.csv,
  SuppFig10_treg_activation_ligand_senders.png.
Sections (run order, after 04): S10_ligand_sources -> S10_ligand_senders.
Run: RUN_SECTIONS=<id>[,<id>...] python <file> runs selected sections only (file order).
"""
import os
import time
import shutil
import tempfile
import subprocess
import warnings
import logging
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
logging.getLogger('fontTools').setLevel(logging.WARNING)
_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'; R = '/share/fsmresfiles/UC/revision_jenny'
F6S = f'{R}/Fig6_treg_antitnf_taurus/supplementary'; HERE = f'{R}/SuppFig10_treg_activation_ligand_senders'
CACHE = f'{F6S}/cache/cn10_cytokine_cells.parquet'; CACHE2 = f'{HERE}/cache/treg_activation_ligand_cells.parquet'
CSV = f'{HERE}/treg_activation_ligand_sources.csv'; PNG = f'{HERE}/SuppFig10_treg_activation_ligand_senders.png'
t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)


# ======================================================================================================================
# Section S10_ligand_sources: senders of the Treg activation program ligands; all-cell and main-source models; decomposition of
# the change (coarse groups, eligibility rule, model and decomposition as 04 F6B_coarse_screen). Output: treg_activation_ligand_sources.csv
# ======================================================================================================================
if _run('S10_ligand_sources'):
    KEYS = ['Patient', 'Remission_status', 'Treatment']; MIN_N = 10; MIN_PRE_PAT = 3; MAIN_SHARE = 20; DRV_PCT = 25
    FAM = {'IL-2/STAT5': ['IL2'], 'TNFR2 non-canonical NF-kB': ['TNF', 'LTA'],
           'Co-stimulatory/suppressive receptor module': ['CD80', 'CD86', 'ICOSLG', 'TNFSF4', 'TNFSF18', 'CD70', 'CD274', 'PDCD1LG2', 'CD40LG', 'IL10'],
           'TNFRSF9 / 4-1BB': ['TNFSF9']}
    ALIAS = {'TNFSF4': 'OX40L', 'TNFSF18': 'GITRL', 'CD274': 'PD-L1', 'PDCD1LG2': 'PD-L2', 'TNFSF9': '4-1BBL', 'CD40LG': 'CD40L'}
    GENES = [g for v in FAM.values() for g in v]; G2F = {g: f for f, v in FAM.items() for g in v}
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')

    # cells (cytokine cell cache of 04) + ligand columns (log-norm + raw counts; own cache, same cell order)
    C = pd.read_parquet(CACHE)[['sample_id'] + KEYS + ['aucell_cell_type_short', 'aucell_cell_category_from_type']]; log('cells', C.shape)
    if not os.path.exists(CACHE2):
        import scanpy as sc
        ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); log('loaded', ad.shape)   # X = log1p(CP10k); layers['counts'] = raw UMIs
        gi = [ad.var_names.get_loc(g) for g in GENES]
        cnt = ad.layers['counts'][:, gi].toarray(); lg = ad.X[:, gi].toarray()
        L = pd.DataFrame(index=ad.obs_names.copy())
        for j, g in enumerate(GENES): L[g] = lg[:, j]; L[f'{g}_counts'] = cnt[:, j]
        del ad
        with open(CACHE2, 'wb') as fh: L.to_parquet(fh)
    L = pd.read_parquet(CACHE2)
    C = pd.concat([C, L[[c for g in GENES for c in (g, f'{g}_counts')]].astype(np.float64)], axis=1)   # float64: exact decomposition sums
    C['biopsy'] = C.sample_id

    # coarse groups (as 04); Cycling joins Epithelial; Neuron is in no group
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
    T2G = {t: g for g, ts in COARSE.items() for t in ts}
    C['grp'] = C.aucell_cell_type_short.map(T2G).fillna(NOGRP)
    GROUPS = list(COARSE)
    MEMBERS = {g: '; '.join(v) for g, v in COARSE.items()}; MEMBERS[ALL] = 'every cell'; MEMBERS[NOGRP] = 'Neuron'
    GRP = [('rem_pre', 'Remission', 'Pre'), ('rem_post', 'Remission', 'Post'), ('nonrem_pre', 'Non_Remission', 'Pre'), ('nonrem_post', 'Non_Remission', 'Post')]
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}
    BM = C.groupby('biopsy')[KEYS].first().assign(n_all=C.groupby('biopsy').size())

    # part 1: sources (share of pre-treatment transcripts, % expressing, mean) per ligand x coarse group; top 3 fine types
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
    P1['source_rank'] = P1[P1.cell_group.isin(GROUPS)].groupby('feature').share_of_transcripts_pre.rank(ascending=False, method='first')
    P1['main_source'] = P1.cell_group.isin(GROUPS) & ((P1.share_of_transcripts_pre >= MAIN_SHARE) | (P1.source_rank == 1))
    MAIN = {g: list(P1[(P1.feature == g) & P1.main_source].sort_values('source_rank').cell_group) for g in GENES}
    FT = []
    Cp = C[pre]
    for g in GENES:
        s = Cp.groupby('aucell_cell_type_short')[f'{g}_counts'].sum() / TOT_PRE[g] * 100; s = s.sort_values(ascending=False).head(3)
        pe = Cp.groupby('aucell_cell_type_short')[g].apply(lambda x: 100 * (x > 0).mean())
        for i, (t, v) in enumerate(s.items()):
            FT.append(dict(part='1b top fine types', program=G2F[g], feature=g, cell_group=T2G.get(t, NOGRP), fine_type=t, source_rank=i + 1,
                           share_of_transcripts_pre=v, pct_cells_expressing_pre=pe[t], n_cells_pre=int((Cp.aucell_cell_type_short == t).sum())))
    FT = pd.DataFrame(FT)

    # part 3: decomposition of the all-cell change (abundance x expression, two-factor midpoint split)
    DEC = []
    cnt = C.groupby(['biopsy', 'grp']).size().unstack(fill_value=0).reindex(BM.index)
    frac = cnt.div(BM.n_all, axis=0)
    for f in GENES:
        mg = C.groupby(['biopsy', 'grp'])[f].mean().unstack().reindex(BM.index)
        contrib = (frac * mg).fillna(0)
        allm = C.groupby('biopsy')[f].mean().reindex(BM.index)
        for rs, rlab in [('Remission', 'remission'), ('Non_Remission', 'non-remission')]:
            b0 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Pre')]; b1 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Post')]
            tot0, tot1 = allm[b0].mean(), allm[b1].mean(); dtot = tot1 - tot0
            for g in GROUPS + [NOGRP]:
                C0, C1 = contrib.loc[b0, g].mean(), contrib.loc[b1, g].mean(); F0, F1 = frac.loc[b0, g].mean(), frac.loc[b1, g].mean()
                W0 = C0 / F0 if F0 > 0 else 0.0; W1 = C1 / F1 if F1 > 0 else 0.0
                DEC.append(dict(feature=f, program=G2F[f], cell_group=g, remission_group=rlab, contribution_pre=C0, contribution_post=C1, contribution_to_change=C1 - C0,
                                pct_of_all_cell_change=100 * (C1 - C0) / dtot if dtot != 0 else np.nan, abundance_part=(F1 - F0) * (W0 + W1) / 2,
                                expression_part=(W1 - W0) * (F0 + F1) / 2, mean_fraction_pre=100 * F0, mean_fraction_post=100 * F1,
                                weighted_group_mean_pre=W0, weighted_group_mean_post=W1, all_cell_mean_pre=tot0, all_cell_mean_post=tot1,
                                all_cell_change=dtot, n_biopsies_pre=len(b0), n_biopsies_post=len(b1)))
    DEC = pd.DataFrame(DEC)
    W = DEC.pivot_table(index=['feature', 'program', 'cell_group'], columns='remission_group',
                        values=['contribution_to_change', 'abundance_part', 'expression_part', 'all_cell_change'])
    DD = pd.DataFrame({k: W[k]['remission'] - W[k]['non-remission'] for k in ['contribution_to_change', 'abundance_part', 'expression_part', 'all_cell_change']}).reset_index()
    DD['pct_of_all_cell_change'] = 100 * DD.contribution_to_change / DD.all_cell_change
    DD['remission_group'] = 'differential (remission − non-remission)'
    DEC = pd.concat([DEC, DD], ignore_index=True)
    DEC['part'] = '3 decomposition'
    # main sources = transcript sources (>= 20% of transcripts, or the top group) + differential drivers (>= 25% of the gross differential change)
    DIFF = DEC[DEC.remission_group.str.startswith('differential') & DEC.cell_group.isin(GROUPS)].copy()
    DIFF['pct_of_gross_differential'] = 100 * DIFF.contribution_to_change.abs() / DIFF.groupby('feature').contribution_to_change.transform(lambda x: x.abs().sum())
    DRIVER = {g: list(DIFF[(DIFF.feature == g) & (DIFF.pct_of_gross_differential >= DRV_PCT)].sort_values('pct_of_gross_differential', ascending=False).cell_group) for g in GENES}
    TSRC = {g: list(MAIN[g]) for g in GENES}
    MAIN = {g: TSRC[g] + [x for x in DRIVER[g] if x not in TSRC[g]] for g in GENES}
    ROLE = lambda g, cg: '; '.join(r for r, ok in [('transcript source', cg in TSRC[g]), ('differential driver', cg in DRIVER[g])] if ok)

    # part 2: model input (all-cell mean; main-source within-group expression; main-source abundance)
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

    # biopsy-level model: y ~ post * remission + (1 | patient) + (1 | patient:timepoint); lme4 + lmerTest REML, Satterthwaite df;
    # y standardised (estimate back-transformed; estimate_sd in SD units); bobyqa retry if the default optimizer errors
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
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)', 'baseline_difference': 'Baseline difference (remission − non-remission, pre)'}
    RR['contrast'] = RR.contrast_key.map(CON)
    tq = stats.t.ppf(0.975, RR.df); RR['ci_low'] = RR.estimate - tq * RR.se; RR['ci_high'] = RR.estimate + tq * RR.se
    RR = DESC.merge(RR, on='feature_key', how='left')
    # BH within the ligand family per contrast (all-cell; main-source expression) and across the tested groups (abundance)
    RR['q_bh'] = np.nan; RR['bh_family'] = ''
    for (pt_, pr, k), i in RR[RR.p.notna()].groupby(['part', 'program', 'contrast_key']).groups.items():
        RR.loc[i, 'q_bh'] = multipletests(RR.loc[i, 'p'], method='fdr_bh')[1]
        RR.loc[i, 'bh_family'] = (f'{pr}: ligands in the family, this contrast' if pt_ == '2 all-cell model' else
                                  f'{pr}: ligand x main-source tests in the family, this contrast' if pt_.startswith('2b') else 'main-source groups, this contrast')
    RR['n_tests_in_bh_family'] = RR.groupby(['part', 'program', 'contrast_key']).q_bh.transform('count')

    DEC['main_source'] = [cg in MAIN[f] for f, cg in zip(DEC.feature, DEC.cell_group)]
    DEC['source_role'] = [ROLE(f, cg) if cg in GROUPS else '' for f, cg in zip(DEC.feature, DEC.cell_group)]
    DEC = DEC.merge(DIFF[['feature', 'cell_group', 'pct_of_gross_differential']], on=['feature', 'cell_group'], how='left')
    DEC['pct_abundance_of_group_change'] = 100 * DEC.abundance_part / DEC.contribution_to_change
    DEC['unit'] = 'log-normalised expression, log1p(CP10k), all-cell biopsy mean scale; mean_fraction in % of all cells'
    DEC['method'] = ('observed means over biopsies (no model): contribution = mean over biopsies of fraction × within-group mean; '
                     'abundance / expression parts = two-factor midpoint split, add up to the contribution')

    # CSV
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
    RES = ALLR.reindex(columns=COLS)
    with open(CSV, 'w') as fh:
        RES.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', CSV, RES.shape)


# ======================================================================================================================
# Section S10_ligand_senders: Supp 10 dot plot (plot only). Rows = 13 ligands of the four Fig 6A programs; columns = sender groups;
# size = % of pre-treatment transcripts; colour = group contribution to the difference in change; right strip = all-cell P
# ======================================================================================================================
if _run('S10_ligand_senders'):
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize

    COSTIM_CSV = 'Co-stimulatory/suppressive receptor module'
    PROG = [('IL-2/STAT5', ['IL2'], 'IL-2/STAT5'), ('TNFR2', ['TNF', 'LTA'], 'TNFR2 non-canonical NF-kB'),
            ('co-stim', ['CD80', 'CD86', 'CD274', 'PDCD1LG2', 'IL10', 'TNFSF18', 'CD40LG', 'CD70', 'ICOSLG'], COSTIM_CSV),
            ('OX40', ['TNFSF4'], COSTIM_CSV)]
    DROPPED = ['TNFSF9']                                         # 4-1BBL: TNFRSF9 / 4-1BB is not a Fig 6A program
    TAG = {g: t for t, gs, _ in PROG for g in gs}
    ALIAS = {'TNFSF4': 'OX40L', 'TNFSF18': 'GITRL', 'CD274': 'PD-L1', 'PDCD1LG2': 'PD-L2', 'CD40LG': 'CD40L'}
    # receptor(s) on Tregs (CellChatDB receptors of each ligand)
    REC = {'IL2': ['IL2RA', 'IL2RB', 'IL2RG'], 'TNF': ['TNFRSF1B'], 'LTA': ['TNFRSF1B'], 'CD80': ['CD28', 'CTLA4'], 'CD86': ['CD28', 'CTLA4'],
           'CD274': ['PDCD1'], 'PDCD1LG2': ['PDCD1'], 'IL10': ['IL10RA', 'IL10RB'], 'TNFSF18': ['TNFRSF18'], 'CD40LG': ['CD40'],
           'TNFSF4': ['TNFRSF4'], 'CD70': ['CD27'], 'ICOSLG': ['ICOS']}
    GENES = [g for _, gs, _ in PROG for g in gs]

    def rlab(rs):   # IL2RA/IL2RB/IL2RG -> IL2RA/B/G when the receptors differ only by the last letter
        cp = os.path.commonprefix(rs)
        return rs[0] + '/' + '/'.join(r[len(cp):] for r in rs[1:]) if len(rs) > 1 and cp and all(len(r) - len(cp) == 1 for r in rs) else '/'.join(rs)

    def sig2(p):    # 2 significant figures, half-up on the CSV's 3-significant-figure value, trailing zero kept (0.998 -> 1.0)
        d_ = Decimal(repr(float(p))); r_ = d_.quantize(Decimal(1).scaleb(d_.adjusted() - 1), rounding=ROUND_HALF_UP)
        if r_.adjusted() > d_.adjusted(): r_ = r_.quantize(Decimal(1).scaleb(r_.adjusted() - 1))
        return str(r_)
    fmtv = lambda p: '< 0.001' if p < 0.001 else sig2(p)
    COLS = ['DC', 'Macrophage', 'Monocyte-derived Mac', 'B cell', 'Plasma', 'CD4 T', 'CD8 T', 'Stromal', 'Other']
    CLAB = {'Macrophage': 'Mac', 'Monocyte-derived Mac': 'Mo-Mac'}
    STROMAL = ['Endothelial', 'Fibroblast', 'Pericyte']
    OTHER = ['Granulocyte/Neutrophil', 'Other T/NK', 'Epithelial', 'Neuron (no coarse group)']
    MYELOID = 3
    OTHER_MIN = 5.0                                              # % of a ligand's transcripts needed to plot the Other column

    # values from the S10_ligand_sources CSV (no refit)
    d = pd.read_csv(CSV)
    d = d[~d.feature.isin(DROPPED)].copy()
    def col(s): return s.where(~s.isin(STROMAL), 'Stromal').where(~s.isin(OTHER), 'Other')
    S = d[(d.part == '1 source') & (d.cell_group != 'All cells')].copy(); S['col'] = col(S.cell_group)
    SHARE = S.groupby(['feature', 'col']).share_of_transcripts_pre.sum().unstack().reindex(index=GENES, columns=COLS)
    PCOLS = COLS if (SHARE['Other'] > OTHER_MIN).any() else COLS[:-1]
    DC_ = d[d.part == '3 decomposition'].copy(); DC_['col'] = col(DC_.cell_group)

    def contrib(rg):
        x = DC_[DC_.remission_group == rg]
        pre = x.groupby('feature').all_cell_mean_pre.first()
        c = x.groupby(['feature', 'col']).contribution_to_change.sum().unstack().reindex(index=GENES, columns=COLS)
        return c, pre.reindex(GENES)
    cR, preR = contrib('remission'); cN, preN = contrib('non-remission')
    VAL = 100 * (cR.div(preR, axis=0) - cN.div(preN, axis=0))   # contribution to the difference in change, % of pre-treatment all-cell mean
    A = d[d.part == '2 all-cell model'].pivot_table(index='feature', columns='contrast', values='p').reindex(GENES)
    PC = {'rem': 'Change in remission (post − pre)', 'non': 'Change in non-remission (post − pre)', 'diff': 'Differential change (remission − non-remission)'}

    # row order: sender column with the largest |contribution to the difference|, then difference P
    TOPC = VAL[PCOLS].abs().idxmax(axis=1)
    ORD = pd.DataFrame({'top': TOPC, 'ci': TOPC.map(PCOLS.index), 'p': A[PC['diff']]}).rename_axis('g').reset_index()
    ORD = ORD.sort_values(['ci', 'p', 'g'], kind='mergesort')
    GENES = ORD.g.tolist()

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.linewidth': 1.0})
    BLUE, RED, K, GR = '#2166ac', '#b2182b', '#111111', '#8c8c8c'
    VMAX = 75.0
    CMAP = plt.get_cmap('RdBu'); NORM = Normalize(-VMAX, VMAX)
    FS_HD, FS_LIG, FS_REC, FS_P, FS_LEG, FS_AT = 13, 13, 10, 12.5, 11.5, 13.5
    SMAX = 500.0; size = lambda s: SMAX * np.asarray(s) / 100.0    # marker area (pt^2) at 100%
    GAP = 0.5                                                      # extra space (rows) between groups of rows
    ypos, seps, y = {}, [], 0.0
    for i, g in enumerate(GENES):
        if i and TOPC[g] != TOPC[GENES[i - 1]]: seps.append(y - 0.5 + GAP / 2); y += GAP
        ypos[g] = y; y += 1.0
    YMAX = y - 1.0
    UNIT, COLW, PGAP = 0.32, 0.35, 0.13                            # inches per row, per dot column; gap between P columns
    WIN, HIN, DPI = 16.0, 9.0, 300
    NCOL = len(PCOLS); XLO, XHI = -0.5, NCOL - 0.45
    AX0, BOTA = 4.4, 1.0; AXW = COLW * (XHI - XLO); AXH = UNIT * (YMAX + 1.1); TOPA = BOTA + AXH
    fig = plt.figure(figsize=(WIN, HIN), dpi=DPI)
    axin = lambda x, y_, w, h: fig.add_axes([x / WIN, y_ / HIN, w / WIN, h / HIN])
    axA = axin(AX0, BOTA, AXW, AXH)
    SX0 = AX0 + AXW + 0.05
    axS = axin(SX0, BOTA, 2.0, AXH); axS.sharey(axA)               # width reset below from the widest P value per column
    for g in GENES:
        for j, c in enumerate(PCOLS):
            s = SHARE.loc[g, c]
            if s > 0: axA.scatter(j, ypos[g], s=size(s), c=[CMAP(NORM(VAL.loc[g, c]))], edgecolors='#444444', linewidths=0.6, zorder=3)
    axA.set_xlim(XLO, XHI); axA.set_ylim(YMAX + 0.55, -0.55); axA.set_xticks([]); axA.set_yticks([])
    for sp in axA.spines.values(): sp.set_visible(False)
    for j in range(NCOL): axA.axvline(j, color='#eeeeee', lw=0.8, zorder=0)
    for g in GENES: axA.axhline(ypos[g], color='#eeeeee', lw=0.8, zorder=0)
    axA.axvline(MYELOID - 0.5, color=GR, lw=1.0, ls=(0, (3, 3)), zorder=1)
    HOFF = (-3, 4)
    HD = [axA.annotate(CLAB.get(c, c), xy=(j, 1), xycoords=('data', 'axes fraction'), xytext=HOFF, textcoords='offset points', rotation=45,
                       ha='left', va='bottom', rotation_mode='anchor', fontsize=FS_HD, color=K, annotation_clip=False) for j, c in enumerate(PCOLS)]
    # ligand (black) and receptor + program tag (grey) label columns, left of the dots
    yt_ = axA.get_yaxis_transform()
    LIG = {g: axA.text(0, ypos[g], f'{g} ({ALIAS[g]})' if g in ALIAS else g, transform=yt_, ha='left', va='center', fontsize=FS_LIG, color=K) for g in GENES}
    RCP = {g: axA.text(0, ypos[g], f'→ {rlab(REC[g])} [{TAG[g]}]', transform=yt_, ha='left', va='center', fontsize=FS_REC, color='#7a7a7a') for g in GENES}
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); axw = axA.get_window_extent(rend).width
    wl = max(t.get_window_extent(rend).width for t in LIG.values()); wr = max(t.get_window_extent(rend).width for t in RCP.values())
    PT = DPI / 72 / axw                                            # 1 pt in axes fraction
    XR = -(wr / axw + 6 * PT); XL = XR - 5 * PT
    for g in GENES: LIG[g].set_x(XL); LIG[g].set_ha('right'); RCP[g].set_x(XR)
    # P strip: change P in remission, in non-remission and difference P (bold if P < 0.05)
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
    # keys right of the P strip: colour bar, then dot-size scale
    LX = SX0 + SW + 0.2
    CBT = fig.text(LX / WIN, TOPA / HIN, 'Contribution to\ndifference in\nchange (% of\npre-treatment\nmean;\nremission −\nnon-remission)',
                   ha='left', va='top', fontsize=FS_LEG, linespacing=1.05)
    fig.canvas.draw(); cbt_y0 = CBT.get_window_extent(rend).y0 / DPI
    CBH = 1.2; cax = axin(LX, cbt_y0 - 0.12 - CBH, 0.16, CBH)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation='vertical', ticks=[-75, -50, -25, 0, 25, 50, 75])
    cb.ax.tick_params(labelsize=FS_LEG, length=3); cb.outline.set_linewidth(0.6)
    fig.canvas.draw(); cax_y0 = cax.get_tightbbox(rend).y0 / DPI
    SZT = fig.text(LX / WIN, (cax_y0 - 0.22) / HIN, '% of pre-\ntreatment\ntranscripts', ha='left', va='top', fontsize=FS_LEG, linespacing=1.05)
    fig.canvas.draw(); szt_y0 = SZT.get_window_extent(rend).y0 / DPI
    SZS = [10, 25, 50, 75]; SSP = 0.28
    sax = axin(LX, szt_y0 - 0.08 - SSP * len(SZS), 1.0, SSP * len(SZS)); sax.axis('off'); sax.set_xlim(0, 1.0); sax.set_ylim(SSP * len(SZS), 0)
    for i, s in enumerate(SZS):
        sax.scatter(0.16, SSP * (i + 0.5), s=size(s), facecolors='white', edgecolors='#444444', linewidths=0.8, clip_on=False)
        sax.text(0.40, SSP * (i + 0.5), f'{s}%', ha='left', va='center', fontsize=FS_LEG)
    # axis titles
    axA.text(XS0 - 5 * PT, 0.5, 'Treg Activation Program (Ligand → Treg Receptor)', transform=axA.transAxes, rotation=90, ha='right', va='center', fontsize=FS_AT)
    fig.canvas.draw()
    ytop = max(h.get_window_extent(rend).y1 for h in HD + SH)
    yT = fig.transFigure.inverted().transform((0, ytop + 6 * DPI / 72))[1]
    xT = fig.transFigure.inverted().transform((axA.transData.transform(((NCOL - 1) / 2, 0))[0], 0))[0]
    xP = fig.transFigure.inverted().transform((axS.transData.transform((SW / 2, 0))[0], 0))[0]
    fig.text(xT, yT, 'Sender cell group', ha='center', va='bottom', fontsize=FS_AT)
    fig.text(xP, yT, 'P (biopsy-level mixed model)', ha='center', va='bottom', fontsize=FS_AT)
    fig.canvas.draw()
    fig.savefig(PNG, dpi=DPI, bbox_inches='tight', pad_inches=0.08, facecolor='white')
    log('wrote', PNG)
