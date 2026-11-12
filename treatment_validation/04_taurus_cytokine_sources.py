#!/usr/bin/env python
"""TAURUS CD scRNA-seq (anti-TNF): cell-group sources of the change in OSM, IL1B and IL6 after treatment.

Panels: Fig 6B; Supp Table 11 (cytokine models).
Inputs: TAURUS raw_norm_counts_anno.h5ad.gz (01, 02).
Outputs (Fig6_treg_antitnf_taurus/): supplementary/cache/cn10_cytokine_cells.parquet (also read by 05),
  supplementary/cache/cn10_cycling_lineage_check.csv, supplementary/cn10_cytokine_coarse_screen.csv,
  supplementary/cn10_cytokine_contrib_bootstrap.csv, fig6_cn10_cytokine_combined.png.
Sections (run order): cytokine_cell_cache -> F6B_coarse_screen -> F6B_contribution_se -> F6B_cytokine_sources.
Run: RUN_SECTIONS=<id>[,<id>...] python <file> runs selected sections only (file order).
"""
import os
import re
import math
import time
import shutil
import tempfile
import subprocess
import warnings
import logging

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


# ---- paths and shared definitions
T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'; R = '/share/fsmresfiles/UC/revision_jenny'
TOP = f'{R}/Fig6_treg_antitnf_taurus'; HERE = f'{TOP}/supplementary'
CACHE = f'{HERE}/cache/cn10_cytokine_cells.parquet'; CYC = f'{HERE}/cache/cn10_cycling_lineage_check.csv'
CSV = f'{HERE}/cn10_cytokine_coarse_screen.csv'; BOOT = f'{HERE}/cn10_cytokine_contrib_bootstrap.csv'
PNG = f'{TOP}/fig6_cn10_cytokine_combined.png'
KEYS = ['Patient', 'Remission_status', 'Treatment']; MIN_N = 10; MIN_PRE_PAT = 3
GENES = ['OSM', 'IL1B', 'IL6']; AB = 'Abundance'; FEATS = GENES + [AB]
GRP = [('rem_pre', 'Remission', 'Pre'), ('rem_post', 'Remission', 'Post'), ('nonrem_pre', 'Non_Remission', 'Pre'), ('nonrem_post', 'Non_Remission', 'Post')]
t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')


def p2(p):   # 2 significant figures, trailing zeros kept (0.0020)
    d = max(1 - math.floor(math.log10(abs(p))), 0); return f'{p:.{d}f}'


# coarse cell groups (aucell_cell_type_short -> group); Cycling joins Epithelial (cycling lineage check); Neuron is in no group
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
T2G = {t: g for g, ts in COARSE.items() for t in ts}
NOGRP = 'Neuron (no coarse group)'; ALL = 'All cells'; GROUPS = list(COARSE)
NAME = {'Monocyte-derived Mac': 'Mo-Mac'}                                       # display name only (CSV keeps the analysis name)
disp = lambda s: re.sub(r'\+(?=[\s)A-Z]|$)', '⁺', NAME.get(s, s))
CHG = {'remission': 'Change in remission (post − pre)', 'non-remission': 'Change in non-remission (post − pre)'}
DIFC = 'Differential change (remission − non-remission)'
CON = {'baseline_difference': 'Baseline difference (remission − non-remission, pre)', 'differential_change': DIFC,
       'change_remission': CHG['remission'], 'change_nonremission': CHG['non-remission']}


def load_cells(float64=False):
    """Per-cell cache with biopsy id and coarse group."""
    C = pd.read_parquet(CACHE); C['biopsy'] = C.sample_id
    if float64: C[GENES] = C[GENES].astype('float64')
    C['grp'] = C.aucell_cell_type_short.map(T2G).fillna(NOGRP)
    return C


def decompose(C, BM):
    """All-cell biopsy mean = sum over groups of fraction x within-group mean; per remission group, each group's change in mean
    contribution, split into abundance and expression parts (two-factor midpoint split)."""
    DEC, ALLM = [], {}
    cnt = C.groupby(['biopsy', 'grp']).size().unstack(fill_value=0).reindex(BM.index)
    frac = cnt.div(BM.n_all, axis=0)
    for f in GENES:
        mg = C.groupby(['biopsy', 'grp'])[f].mean().unstack().reindex(BM.index)
        contrib = (frac * mg).fillna(0)
        allm = C.groupby('biopsy')[f].mean().reindex(BM.index); ALLM[f] = allm
        for rs, rlab in [('Remission', 'remission'), ('Non_Remission', 'non-remission')]:
            b0 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Pre')]; b1 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Post')]
            tot0, tot1 = allm[b0].mean(), allm[b1].mean(); dtot = tot1 - tot0
            for g in GROUPS + [NOGRP]:
                C0, C1 = contrib.loc[b0, g].mean(), contrib.loc[b1, g].mean(); F0, F1 = frac.loc[b0, g].mean(), frac.loc[b1, g].mean()
                W0 = C0 / F0 if F0 > 0 else 0.0; W1 = C1 / F1 if F1 > 0 else 0.0
                DEC.append(dict(feature=f, cell_group=g, remission_group=rlab, contribution_pre=C0, contribution_post=C1, contribution_to_change=C1 - C0,
                                pct_of_all_cell_change=100 * (C1 - C0) / dtot, abundance_part=(F1 - F0) * (W0 + W1) / 2,
                                expression_part=(W1 - W0) * (F0 + F1) / 2, mean_fraction_pre=100 * F0, mean_fraction_post=100 * F1,
                                weighted_group_mean_pre=W0, weighted_group_mean_post=W1, all_cell_mean_pre=tot0, all_cell_mean_post=tot1,
                                all_cell_change=dtot, n_biopsies_pre=len(b0), n_biopsies_post=len(b1)))
    return pd.DataFrame(DEC), ALLM


def biopsy_meta(C):
    return C.groupby('biopsy')[KEYS].first().assign(n_all=C.groupby('biopsy').size())


# biopsy-level model: y ~ post * remission + (1 | patient) + (1 | patient:timepoint); lme4 + lmerTest REML, Satterthwaite df;
# y standardised before fitting (estimate back-transformed; estimate_sd in SD units); bobyqa retry if the default optimizer errors
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


def fit_lmm(LONG):
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        LONG.to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/r_out.csv'], check=True)
        return pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ======================================================================================================================
# Section cytokine_cell_cache: per-cell OSM / IL1B / IL6 (log-norm and raw counts) for all cells (read by the Fig 6B sections and 05)
# ======================================================================================================================
if _run('cytokine_cell_cache') and not os.path.exists(CACHE):
    import scanpy as sc
    ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); log('loaded', ad.shape)   # X = log1p(CP10k); layers['counts'] = raw UMIs
    gi = [ad.var_names.get_loc(g) for g in GENES]
    cnt = ad.layers['counts'][:, gi].toarray(); lg = ad.X[:, gi].toarray()
    C = ad.obs[['sample_id'] + KEYS + ['aucell_cell_type_short', 'aucell_cell_category_from_type']].astype(str).copy()
    for j, g in enumerate(GENES): C[g] = lg[:, j]; C[f'{g}_counts'] = cnt[:, j]
    del ad
    os.makedirs(os.path.dirname(CACHE), exist_ok=True); C.to_parquet(CACHE); log('cached', CACHE, C.shape)


# ======================================================================================================================
# Section F6B_coarse_screen: coarse cell-group models of OSM, IL1B, IL6 (within-group expression and abundance) and the
# decomposition of the all-cell change. Panels: Fig 6B dots and P values; Supp Table 11
# ======================================================================================================================
if _run('F6B_coarse_screen'):
    # Cycling lineage check on the data (Cycling joins Epithelial); numbers quoted in the Fig 6B legend
    if not os.path.exists(CYC):
        import anndata as an
        ad = an.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz', backed='r')
        cy = (ad.obs.aucell_cell_type_short == 'Cycling').values; idx = np.where(cy)[0]
        MK = ['CD3D', 'CD3E', 'PTPRC', 'NKG7', 'EPCAM', 'KRT8', 'MKI67']
        X = ad.X[idx][:, [ad.var_names.get_loc(g) for g in MK]]; X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
        rows = [dict(item=f'% of Cycling cells expressing {g}', value=100 * (X[:, j] > 0).mean()) for j, g in enumerate(MK)]
        vc = ad.obs.loc[cy, 'major'].astype(str).value_counts(normalize=True) * 100
        rows += [dict(item=f'% of Cycling cells with TAURUS study annotation major = {k}', value=v) for k, v in vc.items() if v > 0]
        with open(CYC, 'w') as fh: pd.DataFrame(rows).assign(n_cycling_cells=len(idx)).to_csv(fh, index=False)
        del ad

    C = load_cells(); log('cells', C.shape)
    MEMBERS = {g: '; '.join(v) for g, v in COARSE.items()}; MEMBERS[ALL] = 'every cell'; MEMBERS[NOGRP] = 'Neuron'
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}
    BM = biopsy_meta(C)

    def pseudobulk(D, col):
        return D.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'),
                                       Treatment=('Treatment', 'first'), n_cells=(col, 'size'), y=(col, 'mean')).reset_index()

    # parts 1 (within-group expression; biopsies with >= 10 cells of the group) and 2 (abundance, % of all biopsy cells; all biopsies)
    GROWS, EROWS, LONG = [], [], []
    for name in GROUPS + [ALL]:
        m = (C.grp == name).values if name != ALL else np.ones(len(C), bool)
        D = C[m]
        nb = D.groupby('biopsy').size().reindex(BM.index, fill_value=0); bmeta = BM.assign(n=nb)
        keep = bmeta[bmeta.n >= MIN_N]
        row = dict(cell_group=name, n_cells_total=len(D), n_biopsies_with_cells=int((nb > 0).sum()), n_biopsies_kept=len(keep),
                   n_patients_kept=keep.Patient.nunique(), pct_of_all_cells_pre=100 * m[pre].sum() / pre.sum())
        for nm, r_, t_ in GRP:
            s = keep[(keep.Remission_status == r_) & (keep.Treatment == t_)]
            row.update({f'n_patients_{nm}': s.Patient.nunique(), f'n_biopsies_{nm}': len(s), f'n_cells_{nm}': int(s.n.sum())})
        # eligibility: >= 3 patients per remission group with a pre-treatment biopsy of >= 10 group cells; >= 1 such post biopsy per group
        elig = row['n_patients_rem_pre'] >= MIN_PRE_PAT and row['n_patients_nonrem_pre'] >= MIN_PRE_PAT
        post_ok = row['n_biopsies_rem_post'] >= 1 and row['n_biopsies_nonrem_post'] >= 1
        row['eligible'] = elig and post_ok
        row['exclusion_reason'] = '' if row['eligible'] else (f"< {MIN_PRE_PAT} patients with a pre-treatment biopsy with >= {MIN_N} cells in a remission group" if not elig else 'no eligible post-treatment biopsy in a remission group')
        GROWS.append(row)
        Dp = D[D.Treatment == 'Pre']; Dk = D[D.biopsy.isin(keep.index)]
        for f in FEATS:
            e = dict(cell_group=name, feature=f)
            if f == AB:
                PB = BM.assign(y=100 * nb / BM.n_all).reset_index()[['biopsy'] + KEYS + ['y']]
                PB['n_cells'] = BM.n_all.values
                e.update(share_of_transcripts_pre=np.nan, pct_cells_expressing_pre=np.nan, mean_expression_pre=np.nan,
                         biopsies_used='all biopsies (fraction of all cells; no drop)')
                e['estimable'] = row['eligible'] and name != ALL
                if name == ALL: e['exclusion_reason_feature'] = 'not applicable (fraction = 100% by definition)'
            else:
                PB = pseudobulk(Dk, f) if len(Dk) else None
                e.update(share_of_transcripts_pre=100 * Dp[f'{f}_counts'].sum() / TOT_PRE[f], pct_cells_expressing_pre=100 * (Dp[f] > 0).mean(),
                         mean_expression_pre=Dp[f].mean(), biopsies_used=f'biopsies with >= {MIN_N} cells of the group')
                e['estimable'] = bool(row['eligible'] and PB is not None and PB.y.std() > 0 and (Dk[f] > 0).any())
            if row['eligible']:
                for nm, r_, t_ in GRP:
                    s = PB[(PB.Remission_status == r_) & (PB.Treatment == t_)].y
                    e[f'mean_of_biopsy_means_{nm}'] = s.mean(); e[f'sem_of_biopsy_means_{nm}'] = s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else np.nan
                    e[f'n_biopsies_model_{nm}'] = len(s)
            if e['estimable']:
                k = f'{name}||{f}'
                LONG.append(PB.assign(feature_key=k, pt=PB.Patient + '_' + PB.Treatment, post=(PB.Treatment == 'Post').astype(float),
                                      rem=(PB.Remission_status == 'Remission').astype(float))[['feature_key', 'Patient', 'pt', 'post', 'rem', 'y']])
            EROWS.append(e)
    for f in GENES:   # Neuron share of transcripts (decomposition term only)
        EROWS.append(dict(cell_group=NOGRP, feature=f, share_of_transcripts_pre=100 * C.loc[pre & (C.grp == NOGRP).values, f'{f}_counts'].sum() / TOT_PRE[f],
                          pct_cells_expressing_pre=100 * (C.loc[pre & (C.grp == NOGRP).values, f] > 0).mean(), estimable=False,
                          exclusion_reason_feature='not a coarse group (decomposition term only)'))
    G = pd.DataFrame(GROWS); E = pd.DataFrame(EROWS)
    LONG = pd.concat(LONG, ignore_index=True)

    log(f'fitting {LONG.feature_key.nunique()} group x feature models in lme4 / lmerTest ...')
    RR = fit_lmm(LONG)
    RR['cell_group'] = RR.feature_key.str.split(r'\|\|').str[0]; RR['feature'] = RR.feature_key.str.split(r'\|\|').str[1]
    tq = stats.t.ppf(0.975, RR.df); RR['ci_low'] = RR.estimate - tq * RR.se; RR['ci_high'] = RR.estimate + tq * RR.se
    base = E[E.cell_group != NOGRP].merge(G, on='cell_group')
    full = base.assign(_k=1).merge(pd.DataFrame({'contrast_key': list(CON), '_k': 1}), on='_k').drop(columns='_k')
    RR = full.merge(RR.drop(columns=['feature_key']), on=['cell_group', 'feature', 'contrast_key'], how='left')
    RR['contrast'] = RR.contrast_key.map(CON)
    # BH across the 13 coarse groups, per feature and contrast (All cells not in the family)
    RR['q_bh'] = np.nan
    for (f, k), i in RR[RR.p.notna() & (RR.cell_group != ALL)].groupby(['feature', 'contrast_key']).groups.items():
        RR.loc[i, 'q_bh'] = multipletests(RR.loc[i, 'p'], method='fdr_bh')[1]
    RR['n_tests_in_bh_family'] = RR.groupby(['feature', 'contrast_key']).q_bh.transform('count').where(RR.cell_group != ALL)
    RR['part'] = np.where(RR.feature == AB, '2 abundance', '1 expression')
    RR['unit'] = np.where(RR.feature == AB, '% of all cells in the biopsy (estimate in percentage points)', 'log-normalised expression, log1p(CP10k)')
    ORD = E[E.feature.isin(GENES) & E.cell_group.isin(GROUPS)].groupby('cell_group').share_of_transcripts_pre.sum().sort_values(ascending=False)
    ROWS = list(ORD.index) + [ALL]
    RR['_o'] = RR.cell_group.map({g: i for i, g in enumerate(ROWS)}); RR['_f'] = RR.feature.map({f: i for i, f in enumerate(FEATS)})
    RR['_c'] = RR.contrast_key.map({k: i for i, k in enumerate(CON)})
    RR = RR.sort_values(['_f', '_o', '_c']).reset_index(drop=True)

    # part 3: decomposition of the all-cell change (observed biopsy values), plus remission - non-remission rows
    DEC, _ = decompose(C, BM)
    W = DEC.pivot_table(index=['feature', 'cell_group'], columns='remission_group',
                        values=['contribution_to_change', 'abundance_part', 'expression_part', 'all_cell_change'])
    DD = pd.DataFrame({k: W[k]['remission'] - W[k]['non-remission']
                       for k in ['contribution_to_change', 'abundance_part', 'expression_part', 'all_cell_change']}).reset_index()
    DD['pct_of_all_cell_change'] = 100 * DD.contribution_to_change / DD.all_cell_change
    DD['remission_group'] = 'differential (remission − non-remission)'
    DEC = pd.concat([DEC, DD], ignore_index=True)
    DEC['part'] = '3 decomposition'

    # CSV
    RR['mapping'] = RR.cell_group.map(MEMBERS); DEC['mapping'] = DEC.cell_group.map(MEMBERS)
    RR['model_formula'] = 'biopsy value ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df; no covariates'
    RR['eligibility_rule'] = (f'group tested if >= {MIN_PRE_PAT} patients per remission group have a pre-treatment biopsy with >= {MIN_N} cells of the group '
                              'and each remission group keeps >= 1 such post-treatment biopsy')
    RR['exclusion_reason'] = RR.exclusion_reason.where(RR.exclusion_reason != '', RR.get('exclusion_reason_feature')).fillna('')
    DEC['unit'] = 'log-normalised expression, log1p(CP10k), all-cell biopsy mean scale; mean_fraction in % of all cells'
    DEC['method'] = ('observed means over biopsies (no model): contribution = mean over biopsies of fraction × within-group mean; '
                     'abundance / expression parts = two-factor midpoint split, add up to the contribution')
    NEU = E[E.cell_group == NOGRP][['cell_group', 'feature', 'share_of_transcripts_pre', 'pct_cells_expressing_pre', 'exclusion_reason_feature']].rename(columns={'exclusion_reason_feature': 'exclusion_reason'})
    NEU['part'] = '1 expression'; NEU['mapping'] = MEMBERS[NOGRP]
    COLS = (['part', 'cell_group', 'mapping', 'feature', 'contrast', 'remission_group', 'estimate', 'ci_low', 'ci_high', 'se', 'estimate_sd', 'df', 'statistic',
             'p', 'q_bh', 'n_tests_in_bh_family', 'eligible', 'estimable', 'exclusion_reason', 'share_of_transcripts_pre', 'pct_cells_expressing_pre',
             'mean_expression_pre', 'pct_of_all_cells_pre']
            + [f'{a}_{nm}' for a in ['mean_of_biopsy_means', 'sem_of_biopsy_means', 'n_biopsies_model'] for nm, _, _ in GRP]
            + ['contribution_pre', 'contribution_post', 'contribution_to_change', 'pct_of_all_cell_change', 'abundance_part', 'expression_part',
               'mean_fraction_pre', 'mean_fraction_post', 'weighted_group_mean_pre', 'weighted_group_mean_post', 'all_cell_mean_pre', 'all_cell_mean_post',
               'all_cell_change', 'n_biopsies_pre', 'n_biopsies_post', 'biopsies_used', 'n_cells_total', 'n_biopsies_kept', 'n_patients_kept']
            + [f'{a}_{nm}' for nm, _, _ in GRP for a in ['n_patients', 'n_biopsies', 'n_cells']]
            + ['var_patient', 'var_patient_timepoint', 'var_residual', 'singular', 'converged', 'optimizer', 'unit', 'model_formula', 'eligibility_rule', 'method'])
    RES = pd.concat([RR, NEU, DEC], ignore_index=True).reindex(columns=COLS)
    RES['estimable'] = RES.estimable.astype('boolean')
    with open(CSV, 'w') as fh:
        RES.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', CSV, RES.shape)


# ======================================================================================================================
# Section F6B_contribution_se: net % change of OSM, IL1B, IL6 with model-based SE and 95% CI (Fig 6B bars), plus a
# patient-level cluster bootstrap (CSV only). Output: cn10_cytokine_contrib_bootstrap.csv (Fig 6B source data)
# ======================================================================================================================
if _run('F6B_contribution_se'):
    N_BOOT = 2000; SEED = 20261105
    RES = pd.read_csv(CSV)
    A = RES[(RES.part == '1 expression') & (RES.cell_group == ALL) & RES.feature.isin(GENES)].set_index(['feature', 'contrast'])
    PV = {(f, k): float(A.loc[(f, c), 'p']) for f in GENES for k, c in list(CHG.items()) + [('difference', DIFC)]}
    C = load_cells(float64=True)
    BM = biopsy_meta(C)
    DEC, ALLM = decompose(C, BM)
    DEC['pct_net'] = 100 * DEC.all_cell_change / DEC.all_cell_mean_pre

    # patient-level cluster bootstrap of the net % change: patients resampled within remission group; all biopsies used;
    # a replicate without a pre- or post-treatment biopsy is redrawn; same patient draws for all three genes; percentile CI
    rng = np.random.default_rng(SEED); BS = []
    for rs, rlab in [('Remission', 'remission'), ('Non_Remission', 'non-remission')]:
        bm = BM[BM.Remission_status == rs]; pats = np.array(sorted(bm.Patient.unique()))
        is0 = (bm.Treatment == 'Pre').values; is1 = (bm.Treatment == 'Post').values
        pidx = pd.Index(pats).get_indexer(bm.Patient)
        n0 = np.bincount(pidx, weights=is0, minlength=len(pats)); n1 = np.bincount(pidx, weights=is1, minlength=len(pats))
        S0 = {f: np.bincount(pidx, weights=ALLM[f][bm.index].values * is0, minlength=len(pats)) for f in GENES}
        S1 = {f: np.bincount(pidx, weights=ALLM[f][bm.index].values * is1, minlength=len(pats)) for f in GENES}
        Wb = []; redraw = 0
        while len(Wb) < N_BOOT:
            w = np.bincount(rng.integers(0, len(pats), len(pats)), minlength=len(pats))
            if (w @ n0) == 0 or (w @ n1) == 0: redraw += 1; continue
            Wb.append(w)
        Wb = np.array(Wb, dtype=float)
        for f in GENES:
            obs = 100 * ((S1[f].sum() / n1.sum()) - (S0[f].sum() / n0.sum())) / (S0[f].sum() / n0.sum())
            m0 = (Wb @ S0[f]) / (Wb @ n0); m1 = (Wb @ S1[f]) / (Wb @ n1); b = 100 * (m1 - m0) / m0
            lo, hi = np.percentile(b, [2.5, 97.5])
            BS.append(dict(gene=f, remission_group=rlab, net_pct_change=obs, bootstrap_ratio_ci95_low=lo, bootstrap_ratio_ci95_high=hi,
                           bootstrap_ratio_se=b.std(ddof=1), change_p=PV[(f, rlab)], difference_p=PV[(f, 'difference')],
                           n_patients=len(pats), n_patients_pre=int((n0 > 0).sum()), n_patients_post=int((n1 > 0).sum()),
                           n_biopsies_pre=int(n0.sum()), n_biopsies_post=int(n1.sum()), all_cell_mean_pre=S0[f].sum() / n0.sum(),
                           bootstrap_ratio_n_boot=N_BOOT, bootstrap_ratio_seed=SEED, bootstrap_ratio_n_redrawn_no_pre_or_post=redraw,
                           bootstrap_ratio_method='patient-level cluster bootstrap within remission group; percentile 95% CI; all biopsies, not paired-only',
                           p_source='All cells, biopsy-level mixed model: biopsy mean ~ timepoint x remission + (1|patient) + (1|patient:timepoint)'))
    BS = pd.DataFrame(BS)

    # model-based estimate, SE and 95% CI in %: All cells change estimate / observed pre-treatment all-cell mean (the diamond's denominator)
    MOD = []
    for f in GENES:
        for rl in ['remission', 'non-remission']:
            r_ = A.loc[(f, CHG[rl])]
            dr = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].iloc[0]
            den, obs = dr.all_cell_mean_pre, dr.pct_net
            e, lo, hi = (100 * r_[c] / den for c in ['estimate', 'ci_low', 'ci_high'])
            MOD.append(dict(gene=f, remission_group=rl, model_pct_change=e, model_ci95_low=lo, model_ci95_high=hi,
                            model_estimate_lognorm=r_.estimate, model_ci95_low_lognorm=r_.ci_low, model_ci95_high_lognorm=r_.ci_high,
                            model_se_lognorm=r_.se, model_se_pct=100 * r_.se / den, model_minus_observed_pct_points=e - obs,
                            ci_excludes_0=(lo > 0) or (hi < 0)))
    BS = pd.DataFrame(MOD).merge(BS, on=['gene', 'remission_group'], validate='1:1')
    BS['model_ci_source'] = ('All cells change (post - pre) from the biopsy-level mixed model (cn10_cytokine_coarse_screen.csv), log-normalised '
                             'units / observed pre-treatment all-cell mean of the group (constant denominator) x 100')
    BS['model_minus_se_pct'] = BS.model_pct_change - BS.model_se_pct; BS['model_plus_se_pct'] = BS.model_pct_change + BS.model_se_pct
    BS = BS[['gene', 'remission_group', 'net_pct_change', 'model_pct_change', 'model_se_pct', 'model_minus_se_pct', 'model_plus_se_pct', 'model_ci95_low',
             'model_ci95_high', 'model_minus_observed_pct_points', 'ci_excludes_0', 'change_p', 'difference_p', 'all_cell_mean_pre', 'model_estimate_lognorm',
             'model_ci95_low_lognorm', 'model_ci95_high_lognorm', 'model_se_lognorm', 'n_patients', 'n_patients_pre', 'n_patients_post', 'n_biopsies_pre',
             'n_biopsies_post', 'model_ci_source', 'p_source'] + [c for c in BS.columns if c.startswith('bootstrap_ratio')]]
    for c in BS.columns:
        if BS[c].dtype.kind == 'f': BS[c] = BS[c].map(sig3)
    with open(BOOT, 'w') as fh: BS.to_csv(fh, index=False)
    log('wrote', BOOT)


# ======================================================================================================================
# Section F6B_cytokine_sources: Fig 6B (plot only). Stacked cell-group contributions to the % change in the all-cell mean with
# net change and model SE; abundance dots for main-source and contrast groups; strip of transcript and cell shares
# ======================================================================================================================
if _run('F6B_cytokine_sources'):
    MAIN_TX = 15.0          # main source: >= 15% of the pre-treatment transcripts of OSM, IL1B or IL6
    CONTRAST_CELLS = 8.0    # contrast: remaining groups with >= 8% of pre-treatment cells
    RES = pd.read_csv(CSV)
    A = RES[(RES.part == '1 expression') & (RES.cell_group == ALL) & RES.feature.isin(GENES)].set_index(['feature', 'contrast'])
    PV = {(f, k): float(A.loc[(f, c), 'p']) for f in GENES for k, c in list(CHG.items()) + [('difference', DIFC)]}
    RR = RES[RES.part.isin(['1 expression', '2 abundance']) & RES.contrast.notna()].copy()
    RR['contrast_key'] = RR.contrast.map({v: k for k, v in CON.items()})

    # strip shares and decomposition from the cell cache
    C = load_cells(float64=True)
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}
    SHX = {(g, f): 100 * C.loc[pre & (C.grp == g).values, f'{f}_counts'].sum() / TOT_PRE[f] for g in GROUPS for f in GENES}
    PCX = {g: 100 * (pre & (C.grp == g).values).sum() / pre.sum() for g in GROUPS}
    MAXTX = {g: max(SHX[(g, f)] for f in GENES) for g in GROUPS}
    MAIN = sorted([g for g in GROUPS if MAXTX[g] >= MAIN_TX], key=lambda g: -MAXTX[g])
    CONTRAST = sorted([g for g in GROUPS if g not in MAIN and PCX[g] >= CONTRAST_CELLS], key=lambda g: -PCX[g])
    ROWS = MAIN + CONTRAST
    DEC, _ = decompose(C, biopsy_meta(C))
    DEC['pct_contribution'] = 100 * DEC.contribution_to_change / DEC.all_cell_mean_pre
    DEC['pct_net'] = 100 * DEC.all_cell_change / DEC.all_cell_mean_pre

    # model estimate +/- SE in % (no refit)
    EST = {}; SEP = {}
    for f in GENES:
        for rl in ['remission', 'non-remission']:
            r_ = A.loc[(f, CHG[rl])]
            den = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].all_cell_mean_pre.iloc[0]
            EST[(f, rl)] = 100 * r_.estimate / den; SEP[(f, rl)] = 100 * r_.se / den

    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.legend_handler import HandlerTuple
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14, 'xtick.labelsize': 13, 'ytick.labelsize': 13,
                         'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0})
    K, GR = '#111111', '#A6A6A6'
    # coloured groups: >= 10% of a bar's net change or >= 10% of a gene's pre-treatment transcripts; the rest pooled as Other groups
    big = DEC.assign(r=DEC.contribution_to_change.abs() / DEC.all_cell_change.abs())
    SRCG = {g for (g, f), v in SHX.items() if v >= 10}
    COLG = [g for g in GROUPS if (big[big.cell_group == g].r >= 0.10).any() or g in SRCG]
    PAL = {'Macrophage': '#D55E00', 'Monocyte-derived Mac': '#E69F00', 'DC': '#CC79A7', 'Granulocyte/Neutrophil': '#F0E442', 'B cell': '#7B3294',
           'Plasma': '#C2A5CF', 'CD4 T': '#0072B2', 'CD8 T': '#56B4E9', 'Other T/NK': '#009E73', 'Fibroblast': '#8C6D31',
           'Pericyte': '#C9B37E', 'Endothelial': '#4D4D4D', 'Epithelial': '#80CDC1'}
    OTH = 'Other groups'; PAL[OTH] = '#D9D9D9'
    BARCOL = {g: (PAL[g] if g in COLG else PAL[OTH]) for g in ROWS}

    # one gridspec row: [3 bar panels][row labels][abundance dots][strip][keys]
    fig = plt.figure(figsize=(20, 5.8))
    gs = fig.add_gridspec(1, 7, width_ratios=[1.85, 1.85, 1.85, 1.9, 2.1, 3.9, 3.2], wspace=0.32, left=0.05, right=0.99, bottom=0.1, top=0.83)

    # bars: stacked contributions, net change (diamond), model estimate +/- SE (clipped with arrow and value if beyond the cap)
    CAPF = 0.6
    for i, f in enumerate(GENES):
        a = fig.add_subplot(gs[0, i]); ext = []
        stk = {}
        for rl in ['remission', 'non-remission']:
            s_ = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].pct_contribution
            stk[rl] = (s_[s_ > 0].sum(), s_[s_ < 0].sum())
        bmax = max([0] + [u for u, _ in stk.values()]); bmin = min([0] + [d_ for _, d_ in stk.values()]); span = bmax - bmin
        capH, capL = bmax + CAPF * span, bmin - CAPF * span
        for xi, rl in enumerate(['remission', 'non-remission']):
            d = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].set_index('cell_group')
            s = d.pct_contribution
            vals = [(g, s[g]) for g in COLG] + [(OTH, s.drop(COLG).sum())]
            up = dn = 0.0
            for g, v in vals:
                if v >= 0: a.bar(xi, v, bottom=up, width=0.62, color=PAL[g], edgecolor='white', lw=0.6, zorder=2); up += v
                else: a.bar(xi, v, bottom=dn, width=0.62, color=PAL[g], edgecolor='white', lw=0.6, zorder=2); dn += v
            net = sum(v for _, v in vals)
            e_ = EST[(f, rl)]; lo, hi = e_ - SEP[(f, rl)], e_ + SEP[(f, rl)]
            lo_d, hi_d = max(lo, capL), min(hi, capH)
            a.plot([xi, xi], [lo_d, hi_d], color=K, lw=1.2, solid_capstyle='butt', zorder=3.5)
            for v_, vd_, mk in [(hi, hi_d, '^'), (lo, lo_d, 'v')]:
                if v_ != vd_:
                    a.scatter(xi, vd_, marker=mk, s=45, color=K, lw=0, zorder=3.6)
                    a.annotate(f'{v_:.0f}%', (xi, vd_), xytext=(7, -6 if mk == '^' else 6), textcoords='offset points', ha='left', va='top' if mk == '^' else 'bottom', fontsize=12, zorder=5)
            a.scatter(xi, net, marker='D', s=70, color='white', edgecolor=K, linewidths=1.8, zorder=4)
            ext += [max(up, hi_d), min(dn, lo_d)]
        lo_, hi_ = min(ext + [0]), max(ext + [0]); rg = hi_ - lo_
        a.set_ylim(lo_ - 0.08 * rg, hi_ + 0.08 * rg)
        a.axhline(0, color=K, lw=1.0, zorder=3)
        a.set_xticks([0, 1]); a.set_xticklabels(['Remission', 'Non-\nremission'], fontsize=13); a.set_xlim(-0.55, 1.55)
        a.spines[['top', 'right']].set_visible(False); a.tick_params(axis='y', labelsize=12)
        a.yaxis.set_major_locator(plt.MaxNLocator(5)); a.ticklabel_format(axis='y', style='plain')
        a.set_title(f'Difference P = {p2(PV[(f, "difference")])}', fontsize=12.5, pad=10)
        a.annotate(f, (0.5, 1), xycoords=a.title, xytext=(0, 8), textcoords='offset points', ha='center', va='bottom', fontsize=15, fontstyle='italic')
        if i == 0: a.set_ylabel('% change in all-cell mean\n(relative to pre-treatment)', fontsize=13.5)

    # abundance dots (difference in change, baseline): colour = estimate in SD units, size = 4 P bins, black ring = P < 0.05
    VM = 1.5; norm = TwoSlopeNorm(vcenter=0, vmin=-VM, vmax=VM); cmap = plt.get_cmap('RdBu_r')
    P_SIZES = [520, 330, 140, 60]                     # P < 0.01 | 0.01-0.05 | 0.05-0.1 | P >= 0.1
    P_LABELS = ['P ≥ 0.1', '0.05–0.1', '0.01–0.05', '< 0.01']
    psize = lambda p: P_SIZES[0 if p < 0.01 else 1 if p < 0.05 else 2 if p < 0.1 else 3]
    COLSP = [(AB, k) for k in ['differential_change', 'baseline_difference']]
    ypos = {r: i for i, r in enumerate(ROWS)}
    axl = fig.add_subplot(gs[0, 3]); ax = fig.add_subplot(gs[0, 4], sharey=axl); axs = fig.add_subplot(gs[0, 5], sharey=axl)
    EL = RR[RR.contrast_key.isin(['baseline_difference', 'differential_change'])].set_index(['cell_group', 'feature', 'contrast_key'])
    for j, (f, k) in enumerate(COLSP):
        for r, y in ypos.items():
            x = EL.loc[(r, f, k)]
            ax.scatter(j, y, s=psize(x.p), c=[cmap(norm(np.clip(x.estimate_sd, -VM, VM)))], edgecolors=K if x.p < 0.05 else '#9A9A9A',
                       linewidths=2.0 if x.p < 0.05 else 0.6, zorder=3)
    ax.set_xlim(-0.6, len(COLSP) - 0.4); ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_yticks(list(ypos.values())); ax.tick_params(axis='y', labelleft=False)
    ax.set_xticks(range(len(COLSP))); ax.set_xticklabels(['Difference\nin change', 'Base-\nline'], fontsize=12.5, linespacing=1.0)
    ax.xaxis.tick_top(); ax.tick_params(axis='x', length=0, pad=2)
    ax.set_title('Abundance', fontsize=15, pad=6)
    ax.annotate('% of all cells', (0.5, 1), xycoords=ax.title, xytext=(0, 2), textcoords='offset points', ha='center', va='bottom', fontsize=14, color='#444444')
    ax.axhline(len(MAIN) - 0.5, color=GR, lw=0.8)                                  # separator: main sources | contrast
    for a_ in [ax, axs]:
        a_.spines[['top', 'right', 'bottom']].set_visible(False)
        for r, yy in ypos.items():
            if yy % 2 == 0: a_.axhspan(yy - 0.5, yy + 0.5, color='#F2F2F2', zorder=0, lw=0)
    # row names with a colour swatch (bar colour), grouped as main sources | contrast
    axl.axis('off'); axl.set_xlim(0, 1)
    for r, y in ypos.items():
        axl.scatter(0.45, y, marker='s', s=140, color=BARCOL[r], edgecolor='white', lw=0.6)
        axl.text(0.53, y, disp(r), ha='left', va='center', fontsize=14)
    for lab, blk in [('Main sources', MAIN), ('Contrast', CONTRAST)]:
        yt, yb = ypos[blk[0]] - 0.4, ypos[blk[-1]] + 0.4
        axl.plot([0.35, 0.35], [yb, yt], color=GR, lw=0.8)
        axl.text(0.32, (yt + yb) / 2, lab, rotation=90, ha='right', va='center', fontsize=11.5, fontstyle='italic', color='#7F7F7F')
    # strip: % of pre-treatment transcripts per gene and % of pre-treatment cells
    PCM = plt.get_cmap('Purples'); PVMAX = 40
    for j, f in enumerate(GENES + ['cells']):
        for r, y in ypos.items():
            v = SHX[(r, f)] if f != 'cells' else PCX[r]
            axs.add_patch(plt.Rectangle((j - 0.47, y - 0.45), 0.94, 0.9, color=PCM(0.08 + 0.85 * min(v, PVMAX) / PVMAX), lw=0, zorder=1))
            lab = f'{sig3(round(v, 1)):.1f}' if v < 10 else f'{v:.0f}'
            axs.text(j, y, lab, ha='center', va='center', fontsize=12, color='white' if v > 0.6 * PVMAX else K, zorder=2)
    axs.axvline(2.5, color='white', lw=4, zorder=3)
    axs.set_xlim(-0.55, 3.55); axs.tick_params(axis='y', left=False, labelleft=False)
    axs.set_xticks(range(4)); axs.set_xticklabels(GENES + ['Cells'], fontsize=13); axs.xaxis.tick_top(); axs.tick_params(axis='x', length=0)
    for tl in axs.get_xticklabels()[:3]: tl.set_fontstyle('italic')
    axs.set_title('% of pre-treatment\ntranscripts | cells', fontsize=14, linespacing=1.05)
    axs.spines['left'].set_visible(False)
    # keys: colour scale | P size bins + ring note
    axk = fig.add_subplot(gs[0, 6]); axk.axis('off')
    axk.text(0, 1, 'Remission −\nnon-remission\n(SD of biopsy\nvalues)', ha='left', va='top', fontsize=12.5, linespacing=1.05, transform=axk.transAxes)
    cax = axk.inset_axes([0.02, 0.24, 0.06, 0.48])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[-1.5, 0, 1.5]); cb.ax.set_yticklabels(['−1.5', '0', '1.5'], fontsize=12)
    hd = [axk.scatter([], [], s=sz, c='#BBBBBB', edgecolors='#9A9A9A', linewidths=0.6) for sz in P_SIZES[::-1]]
    axk.legend(hd, P_LABELS, title='P', loc='upper left', bbox_to_anchor=(0.5, 1.02), frameon=False,
               fontsize=12, title_fontsize=13, labelspacing=0.9, borderpad=0.2, handletextpad=0.6)
    axk.text(0.52, 0.42, '○ black ring:\n   P < 0.05', ha='left', va='top', fontsize=12, linespacing=1.25, transform=axk.transAxes)
    # bar colours not keyed by a row swatch, and the net-change / SE symbol
    EXTRA = [g for g in COLG if g not in ROWS] + [OTH]
    hx = [Patch(facecolor=PAL[g], edgecolor='white') for g in EXTRA]
    dia = (Line2D([], [], marker='|', ls='', ms=22, mew=1.2, color=K), Line2D([], [], marker='D', ls='', ms=8, mfc='white', mec=K, mew=1.8))
    axl.legend(hx + [dia], [disp(g) for g in EXTRA] + ['Net change (sum of groups); bar = ± SE (mixed model)'], handler_map={tuple: HandlerTuple(ndivide=1)},
               loc='upper left', bbox_to_anchor=(0.37, -0.01), ncol=len(EXTRA) + 1, frameon=False, fontsize=12.5,
               handlelength=1.1, handletextpad=0.5, columnspacing=1.2, borderpad=0.1)
    fig.savefig(PNG, dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    log('wrote', PNG)
