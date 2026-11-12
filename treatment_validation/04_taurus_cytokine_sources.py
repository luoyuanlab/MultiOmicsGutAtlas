#!/usr/bin/env python
"""TAURUS CD scRNA-seq (anti-TNF): cell-group sources of the change in OSM, IL1B and IL6 after treatment.

Panels made
  Fig 6B  change in the all-cell mean of OSM, IL1B and IL6 split into cell-group contributions; cell-group abundance contrasts
  Supp Table 11 (cytokine models)

Inputs
  /share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny/raw_norm_counts_anno.h5ad.gz (01_taurus_preprocessing.ipynb, 02_taurus_aucell_annotation.R)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/archive/csv/msigdb_human_snapshot.csv.gz (MSigDB snapshot)

Outputs (all in /share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/)
  supplementary/nonremission_persistent_signals.csv (F6B_inflammatory_signals_biopsy; All Cells OSM / IL1B / IL6 models)
  supplementary/cache/cn10_cytokine_cells.parquet (per-cell OSM / IL1B / IL6 values; also read by 05_taurus_ligand_senders.py)
  supplementary/cn10_cytokine_celltype_screen.csv, cn10_cytokine_coarse_screen.csv, cache/cn10_cycling_lineage_check.csv,
    cn10_cytokine_contrib_bootstrap.csv (Supp Table 11 and Fig 6B source data)
  fig6_cn10_cytokine_combined.png and legend (Fig 6B)

Run order of sections
  F6B_inflammatory_signals_biopsy -> cytokine_celltype_screen -> F6B_coarse_screen -> F6B_contribution_se -> F6B_cytokine_sources

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
# Section F6B_inflammatory_signals_biopsy: biopsy-level models of inflammatory ligands, receptors and APC programs by cell group
# Panels: Fig 6B All Cells OSM / IL1B / IL6 P values (Results: difference in change P = 0.035, 0.018, 0.029); its CSV is read by the Fig 6B sections
# ======================================================================================================================
if _run('F6B_inflammatory_signals_biopsy'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/56_nonremission_persistent_signals.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: which DC / OGN+RSPO3+ fibroblast interaction signals of the CD CosMx niches (Fig 5) fall in
remitters but persist or rise in non-remitters? Numbers only, no figure.

Design: the biopsy-level 'patient + patient:timepoint' model of 03_taurus_treg_programs.py (F6A_programs_biopsy).
  Pseudobulk: per biopsy (sample_id), the mean of the per-cell feature over the cells of the feature's cell group. Biopsies with < 10 cells
    of that group are dropped; patients / biopsies / cells kept per remission group x timepoint are reported per feature.
  Model: biopsy mean ~ post * remission + (1 | patient) + (1 | patient:timepoint); R lme4 + lmerTest, REML, Satterthwaite df; no covariates;
    y standardised before fitting (estimates back-transformed to y units).
  Contrasts: change in remission = b_post + b_post:rem; change in non-remission = b_post; differential change = b_post:rem;
    baseline difference (remission - non-remission, pre) = b_rem. 95% CI = estimate +/- t(0.975, Satterthwaite df) x SE.
  BH within family, per contrast: ligands; receptors; programs.
Cell groups (aucell_cell_type_short): DC = cDC1 + cDC2 + CCR7+ DC; SELENOP+ Mac = Mac S+XS-, Mac S+M+S+, Mac S+SG+, Mac S+M+P+; Mo-Mac;
  Inf Mo-Mac; OGN+RSPO3+ Fib; all fibroblasts = 9 fibroblast types (FIB9; SELENOP+ Fib is not among them); All APCs (pooled) = DC +
  SELENOP+ Mac + CD9+ Mac + Mo-Mac + B cell; All cells = every cell of the biopsy.
Features (pre-specified from the CN10 signalling of the CD CosMx niches and the DC co-stimulation axis). y = log-normalised expression
  (X = log1p CP10k):
  Ligands: OSM and IL1B in Mo-Mac, Inf Mo-Mac, SELENOP+ Mac, DC; IL6 in all fibroblasts and OGN+RSPO3+ Fib; TNF in Mo-Mac, Inf Mo-Mac, DC
    (positive control for the drug target); CD80, CD86, CD40, ICOSLG in DC; each ligand also as the biopsy mean over All cells.
  Receptors: OSMR, IL6ST, IL1R1, IL6R, TNFRSF1A in the receiver. Receiver = OGN+RSPO3+ Fib if, after the 10-cell drop, >= 3 remitters and
    >= 3 non-remitters keep both a pre and a post biopsy (pre-specified rule); otherwise all fibroblasts. The non-chosen receiver is fitted
    as a sensitivity (not in a BH family); IL6 in OGN+RSPO3+ Fib likewise leaves the ligand BH family when OGN+RSPO3+ Fib fails the rule.
  Programs (within the cell group): 'Inflammatory Response (Hallmark)' in All APCs (decoupler ULM on HALLMARK_INFLAMMATORY_RESPONSE and
    REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION, MSigDB snapshot); 'Co-Stimulatory Ligands' (CD80, CD86, CD40, ICOSLG, TNFSF4, TNFSF18,
    TNFSF9, CD70) and 'MHC-II Antigen Presentation (10 Core Genes)' in DC (scanpy score_genes, random_state 0).
Descriptives: % of cells with y > 0 and the mean of biopsy means, per remission group x timepoint (kept biopsies).
Output Fig6_treg_antitnf_taurus/supplementary/nonremission_persistent_signals.csv (3 significant figures; Supp Table 11).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time, subprocess, tempfile, shutil
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, decoupler as dc
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'; R = '/share/fsmresfiles/UC/revision_jenny'
    HERE = f'{R}/Fig6_treg_antitnf_taurus/supplementary'; OUT = f'{HERE}/nonremission_persistent_signals.csv'
    MSIG = f'{R}/cd_treg_rebuild_rowfix/archive/csv/msigdb_human_snapshot.csv.gz'
    KEYS = ['Patient', 'Remission_status', 'Treatment']; MIN_N = 10; MIN_PAIRED = 3
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)

    # ================================================================ definitions
    SELM = ['Mac S+XS-', 'Mac S+M+S+', 'Mac S+SG+', 'Mac S+M+P+']; DCT = ['cDC1', 'cDC2', 'CCR7+ DC']
    FIB9 = ['VSTM2A+ Crypt Top Fib', 'ADAMDEC1+ Activ Fib', 'ADAMDEC1+ Fib', 'OGN+RSPO3+ Fib', 'T-interact Fib', 'Activ Fib', 'HHIP+ Myofib',
            'PLA2G2A+ ECM Fib', 'GREM2+ Myofib']
    MHC2 = ['HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DMA', 'HLA-DMB', 'CD74', 'CIITA']
    COSTIM = ['CD80', 'CD86', 'CD40', 'ICOSLG', 'TNFSF4', 'TNFSF18', 'TNFSF9', 'CD70']
    S_MHC2 = 'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION'; S_INFL = 'HALLMARK_INFLAMMATORY_RESPONSE'
    OGNL, FIBL = 'OGN⁺RSPO3⁺ Fib', 'Fibroblasts (All 9 Types)'
    GROUPS = {'Mo-Mac': ['Mo-Mac'], 'Inf Mo-Mac': ['Inf Mo-Mac'], 'SELENOP⁺ Mac': SELM, 'DC': DCT, FIBL: FIB9, OGNL: ['OGN+RSPO3+ Fib'],
              'All APCs (Pooled)': DCT + SELM + ['CD9+ Mac', 'Mo-Mac', 'B cell'], 'All Cells': None}
    MYE = ['Mo-Mac', 'Inf Mo-Mac', 'SELENOP⁺ Mac', 'DC']
    LIG_SPEC = [('OSM', MYE), ('IL1B', MYE), ('IL6', [FIBL, OGNL]), ('TNF', ['Mo-Mac', 'Inf Mo-Mac', 'DC']),
                ('CD80', ['DC']), ('CD86', ['DC']), ('CD40', ['DC']), ('ICOSLG', ['DC'])]
    REC = ['OSMR', 'IL6ST', 'IL1R1', 'IL6R', 'TNFRSF1A']
    PROGS = [('hallmark_inflammatory', 'Inflammatory Response (Hallmark)', 'All APCs (Pooled)'),
             ('costim', 'Co-Stimulatory Ligands', 'DC'), ('mhc2_core', 'MHC-II Antigen Presentation (10 Core Genes)', 'DC')]
    FAM_L, FAM_R, FAM_P = 'Ligands', 'Receptors', 'Programs'
    GENES = list(dict.fromkeys([g for g, _ in LIG_SPEC] + REC))

    # ================================================================ per-cell values
    ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); log('loaded', ad.shape)
    assert ad.raw is None and np.allclose(np.expm1(ad.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)            # X = log1p(CP10k)
    assert ad.obs_names.is_unique
    ct = ad.obs['aucell_cell_type_short'].astype(str).values
    for g, ts in GROUPS.items(): assert ts is None or all(t in set(ct) for t in ts), g
    assert all(x in ad.var_names for x in GENES + MHC2 + COSTIM)
    meta = ad.obs[KEYS + ['sample_id']].astype(str).copy(); meta['biopsy'] = meta.sample_id
    assert meta.groupby('biopsy').Patient.nunique().max() == 1 and meta.groupby('biopsy').Treatment.nunique().max() == 1
    assert meta.groupby('Patient').Remission_status.nunique().max() == 1
    assert set(meta.Treatment) == {'Pre', 'Post'} and set(meta.Remission_status) == {'Remission', 'Non_Remission'}
    Xg = ad[:, GENES].X; Xg = pd.DataFrame(Xg.toarray(), columns=GENES, index=ad.obs_names)
    mask_of = lambda g: np.ones(ad.n_obs, bool) if GROUPS[g] is None else np.isin(ct, GROUPS[g])

    # programs: subset copy, score_genes (random_state 0), ULM on the group's ULM sets
    msig = pd.read_csv(MSIG, usecols=['geneset', 'genesymbol'])
    PC = {}
    sub = ad[mask_of('DC')].copy(); sub.layers['lognorm'] = sub.X.copy()
    sc.tl.score_genes(sub, MHC2, score_name='mhc2_core', random_state=0); sc.tl.score_genes(sub, COSTIM, score_name='costim', random_state=0)
    PC['mhc2_core'] = sub.obs['mhc2_core'].copy(); PC['costim'] = sub.obs['costim'].copy(); del sub
    sub = ad[mask_of('All APCs (Pooled)')].copy(); sub.layers['lognorm'] = sub.X.copy()
    net = msig[msig.geneset.isin([S_MHC2, S_INFL])].drop_duplicates().rename(columns={'geneset': 'source', 'genesymbol': 'target'})
    dc.mt.ulm(data=sub, net=net, layer='lognorm'); A = dc.pp.get_obsm(sub, key='score_ulm')
    PC['hallmark_inflammatory'] = pd.Series(np.asarray(A[:, S_INFL].X).ravel(), index=sub.obs_names); del sub, A
    log('scored programs', {k: len(v) for k, v in PC.items()})

    # ================================================================ feature list
    def cells(g, col):
        mk = mask_of(g); v = Xg.loc[mk, col] if col in GENES else PC[col]
        return meta.loc[v.index].assign(v=v.values)

    def paired(PB):
        k = PB.groupby(['Patient', 'Remission_status']).Treatment.nunique().reset_index(); k = k[k.Treatment == 2]
        return int((k.Remission_status == 'Remission').sum()), int((k.Remission_status == 'Non_Remission').sum())

    def pseudobulk(C):
        PB = C.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                     n_cells=('v', 'size'), y=('v', 'mean'), pct=('v', lambda x: 100 * (x > 0).mean())).reset_index()
        return PB[PB.n_cells >= MIN_N].reset_index(drop=True), PB[PB.n_cells < MIN_N]

    # receiver rule (pre-specified)
    pbO, _ = pseudobulk(cells(OGNL, 'OSMR')); OGN_PAIRS = paired(pbO)
    OGN_OK = OGN_PAIRS[0] >= MIN_PAIRED and OGN_PAIRS[1] >= MIN_PAIRED
    RECV, RECV_SENS = (OGNL, FIBL) if OGN_OK else (FIBL, OGNL)
    log(f'OGN+RSPO3+ Fib, biopsies with >= {MIN_N} cells: {len(pbO)}; paired remitters / non-remitters = {OGN_PAIRS} -> receiver = {RECV}')

    FEATS = []   # family, key, label, group, var, kind, role
    for gene, gs in LIG_SPEC:
        for g in gs + ['All Cells']:
            sens = (g == OGNL and not OGN_OK)
            FEATS.append((FAM_L, f'lig_{gene}_{g}', f'{gene} in {g}', g, gene, 'gene',
                          'sensitivity (OGN⁺RSPO3⁺ Fib below pairing rule; not in BH family)' if sens else 'primary'))
    for gene in REC:
        FEATS.append((FAM_R, f'rec_{gene}_{RECV}', f'{gene} in {RECV}', RECV, gene, 'gene', 'primary'))
    for gene in REC:
        FEATS.append((FAM_R, f'rec_{gene}_{RECV_SENS}', f'{gene} in {RECV_SENS}', RECV_SENS, gene, 'gene', 'sensitivity (alternative receiver; not in BH family)'))
    for k, lab, g in PROGS:
        FEATS.append((FAM_P, f'prog_{k}', f'{lab} in {g}', g, k, 'program', 'primary'))
    FD = pd.DataFrame(FEATS, columns=['family', 'feature_key', 'feature', 'cell_group', 'gene_or_program', 'kind', 'role'])
    assert FD.feature_key.is_unique

    LONG, NROWS = [], []
    GRP = [('rem_pre', 'Remission', 'Pre'), ('rem_post', 'Remission', 'Post'), ('nonrem_pre', 'Non_Remission', 'Pre'), ('nonrem_post', 'Non_Remission', 'Post')]
    CC = {}
    for f in FD.itertuples():
        key = (f.cell_group, f.gene_or_program)
        if key not in CC: CC[key] = cells(f.cell_group, f.gene_or_program)
        C = CC[key]; PB, DROP = pseudobulk(C)
        kept = C[C.biopsy.isin(PB.biopsy)]
        row = dict(feature_key=f.feature_key, n_biopsies_dropped=len(DROP), n_cells_dropped=int(DROP.n_cells.sum()),
                   n_patients=PB.Patient.nunique(), n_biopsies=len(PB), n_cells=int(PB.n_cells.sum()))
        row['n_paired_patients_rem'], row['n_paired_patients_nonrem'] = paired(PB)
        for nm, r_, t_ in GRP:
            s = PB[(PB.Remission_status == r_) & (PB.Treatment == t_)]; kc = kept[(kept.Remission_status == r_) & (kept.Treatment == t_)]
            row.update({f'n_patients_{nm}': s.Patient.nunique(), f'n_biopsies_{nm}': len(s), f'n_cells_{nm}': int(s.n_cells.sum()),
                        f'mean_of_biopsy_means_{nm}': s.y.mean() if len(s) else np.nan,
                        f'pct_cells_expressing_{nm}': (100 * (kc.v > 0).mean() if len(kc) else np.nan) if f.kind == 'gene' else np.nan})
        row['pct_cells_expressing_all'] = 100 * (kept.v > 0).mean() if f.kind == 'gene' else np.nan
        est = all(row[f'n_biopsies_{nm}'] >= 1 for nm, _, _ in GRP) and PB.y.std() > 0
        row['estimable'] = est
        NROWS.append(row)
        if est:
            LONG.append(PB.assign(feature_key=f.feature_key, pt=PB.Patient + '_' + PB.Treatment, post=(PB.Treatment == 'Post').astype(float),
                                  rem=(PB.Remission_status == 'Remission').astype(float))[['feature_key', 'Patient', 'pt', 'post', 'rem', 'y']])
    NN = pd.DataFrame(NROWS); LONG = pd.concat(LONG, ignore_index=True)
    log(f'{len(FD)} features; estimable {int(NN.estimable.sum())}; not estimable: {NN[~NN.estimable].feature_key.tolist()}')

    # ================================================================ lme4 / lmerTest, Satterthwaite (patient + patient:timepoint model)
    RCODE = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); D <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); out <- list()
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1),
           baseline_difference = c(rem = 1))
for (cc in unique(D$feature_key)) {
  d <- D[D$feature_key == cc, ]; d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
  y <- d$y; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  m <- lmer(ys ~ post * rem + (1 | Patient) + (1 | pt), data = d, REML = TRUE)
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m)); stopifnot(length(fe) == 4)
  for (k in names(Ls)) {
    L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
    ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
    out[[length(out) + 1]] <- data.frame(feature_key = cc, contrast_key = k, estimate = ct$Estimate * sy, se = ct$`Std. Error` * sy,
      df = ct$df, statistic = ct$`t value`, p = ct$`Pr(>|t|)`, var_patient = vc$vcov[vc$grp == "Patient"] * sy^2,
      var_patient_timepoint = vc$vcov[vc$grp == "pt"] * sy^2, var_residual = vc$vcov[vc$grp == "Residual"] * sy^2, singular = isSingular(m),
      converged = !any(grepl("converge|gradient|Hessian|eigenvalue", c(m@optinfo$conv$lme4$messages, m@optinfo$warnings))))
  }
}
write.csv(do.call(rbind, out), a[2], row.names = FALSE)
"""
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        LONG.to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        log(f'fitting {LONG.feature_key.nunique()} features in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/r_out.csv'], check=True)
        RR = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(RR) == 4 * LONG.feature_key.nunique() and RR.p.notna().all()
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'baseline_difference': 'Baseline difference (remission − non-remission, pre)'}
    RR['contrast'] = RR.contrast_key.map(CON)
    tq = stats.t.ppf(0.975, RR.df); RR['ci_low'] = RR.estimate - tq * RR.se; RR['ci_high'] = RR.estimate + tq * RR.se
    RR = FD.merge(RR, on='feature_key', how='left').merge(NN, on='feature_key')
    RR['q_bh'] = np.nan
    P1 = RR[(RR.role == 'primary') & RR.p.notna()]
    for _, i in P1.groupby(['family', 'contrast_key']).groups.items():
        RR.loc[i, 'q_bh'] = multipletests(RR.loc[i, 'p'], method='fdr_bh')[1]
    RR['bh_family'] = np.where(RR.role == 'primary', RR.family + ' (n = ' + RR.family.map(FD[FD.role == 'primary'].family.value_counts()).astype(str)
                               + '), BH within family, per contrast', 'not in a BH family')
    RR['model_formula'] = 'biopsy mean ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df; no covariates'
    RR['unit'] = np.where(RR.kind == 'gene', 'log-normalised expression, log1p(CP10k)', np.where(RR.gene_or_program == 'hallmark_inflammatory', 'ULM score', 'score_genes score'))
    RR['min_cells_per_biopsy'] = MIN_N; RR['receiver_rule'] = (f'OGN⁺RSPO3⁺ Fib receiver if >= {MIN_PAIRED} paired remitters and non-remitters after the drop; '
                                                         f'observed {OGN_PAIRS[0]} / {OGN_PAIRS[1]} -> {RECV}')

    COLS = (['family', 'feature', 'cell_group', 'gene_or_program', 'role', 'contrast', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh',
             'bh_family', 'unit', 'model_formula', 'var_patient', 'var_patient_timepoint',
             'var_residual', 'singular', 'converged', 'estimable', 'min_cells_per_biopsy', 'n_patients', 'n_biopsies', 'n_cells',
             'n_biopsies_dropped', 'n_cells_dropped', 'n_paired_patients_rem', 'n_paired_patients_nonrem']
            + [f'{a}_{nm}' for nm, _, _ in GRP for a in ['n_patients', 'n_biopsies', 'n_cells']]
            + [f'mean_of_biopsy_means_{nm}' for nm, _, _ in GRP] + [f'pct_cells_expressing_{nm}' for nm, _, _ in GRP] + ['pct_cells_expressing_all',
               'receiver_rule', 'feature_key'])
    RR['_f'] = RR.feature_key.map({k: i for i, k in enumerate(FD.feature_key)}); RR['_c'] = RR.contrast_key.map({k: i for i, k in enumerate(CON)})
    FO = {FAM_P: 0, FAM_L: 1, FAM_R: 2}; RR['_fam'] = RR.family.map(FO)
    RR = RR.sort_values(['_fam', '_f', '_c'])[COLS].reset_index(drop=True)
    with open(OUT, 'w') as fh:
        RR.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, RR.shape)

    # ================================================================ console summary
    S = RR.copy(); S['feature'] = S.feature.str[:52]; S['contrast'] = S.contrast.str.split(' \\(').str[0]
    print(S[['feature', 'role', 'contrast', 'estimate', 'ci_low', 'ci_high', 'df', 'p', 'q_bh', 'singular']]
          .apply(lambda c: c.map(sig3)).to_string(index=False))
    N = RR.drop_duplicates('feature_key')
    print('\nCounts per feature (patients / biopsies after drop):')
    print(N[['feature', 'n_patients', 'n_biopsies', 'n_biopsies_dropped', 'n_paired_patients_rem', 'n_paired_patients_nonrem']
            + [f'n_patients_{nm}' for nm, _, _ in GRP] + [f'n_biopsies_{nm}' for nm, _, _ in GRP] + ['pct_cells_expressing_all']].to_string(index=False))
    for k, lab in CON.items():
        s = RR[(RR.contrast_key == k) & RR.p.notna()]
        print(f'\n{lab}: P<0.05: {[(f, sig3(e), sig3(p), sig3(q), r[:3]) for f, e, p, q, r in s[s.p < 0.05][["feature", "estimate", "p", "q_bh", "role"]].values]}'
          f' | q<0.05: {s[s.q_bh < 0.05].feature.tolist()}')


# ======================================================================================================================
# Section cytokine_celltype_screen: fine cell-type screen of OSM, IL1B and IL6; builds the per-cell cytokine cache
# Panels: Fig 6B inputs (cache cn10_cytokine_cells.parquet; Supp Table 11)
# ======================================================================================================================
if _run('cytokine_celltype_screen'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/57_cn10_cytokine_celltype_screen.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: cell-type screen for the remission vs non-remission difference in the CN10 cytokines OSM, IL1B
and IL6 (the cytokines of the inflammatory CN10 niche of the CD CosMx data, Fig 5). Every cell type of the annotation is tested and the
sender is read off the data.

Loading, pseudobulk and model as F6B_inflammatory_signals_biopsy.
  y = per biopsy (sample_id), the mean over the biopsy's cells of that cell type of log-normalised expression (X = log1p CP10k).
  Model: y ~ post * remission + (1 | patient) + (1 | patient:timepoint); R lme4 + lmerTest, REML, Satterthwaite df; no covariates;
    y standardised before fitting; if the default optimizer errors, the fit is retried with bobyqa (column optimizer); estimates
    back-transformed to log-normalised units (column estimate) and also kept in SD units of the biopsy means (column estimate_sd).
  Contrasts: change in remission = b_post + b_post:rem; change in non-remission = b_post; differential change = b_post:rem;
    baseline difference (remission - non-remission, pre-treatment) = b_rem. 95% CI = estimate +/- t(0.975, df) x SE.
Cell types: aucell_cell_type_short, all types, with one merge: SELENOP+ Mac = Mac S+XS-, Mac S+M+S+, Mac S+SG+, Mac S+M+P+.
  Lineage groups (aucell_cell_category_from_type): Myeloid (All) = category Myeloid + Granulocyte; Fibroblasts (All) = category Connective
  without the two pericyte types (10 fibroblast types, SELENOP+ Fib included); T Cells (All) = category T; B / Plasma (All) = B + Plasma;
  Epithelial (All) = category Epithelial; Endothelial (lineage = the single Endothelial type). All Cells = every cell (reference row).
Eligibility (set before any model was fitted; per cell type, gene independent): biopsies with < 10 cells of the type are dropped; the type
  is eligible if >= 3 Remission and >= 3 Non-remission patients keep an eligible pre-treatment biopsy. A type must also keep >= 1
  post-treatment biopsy in each group for the interaction model to be estimable; a gene with no expressing cell in the kept biopsies is not
  estimable. Every excluded type is listed with its reason and counts.
Multiple testing: BH across cell types, per gene and contrast, separately within the fine-type family and the lineage family (lineages
  overlap the fine types). All Cells is not in a BH family. P and q are both reported.
Descriptives per cell type x gene (pre-treatment cells of all biopsies, before the 10-cell drop): % of cells expressing (count > 0), mean
  log-normalised expression, and share of total transcripts = the type's sum of raw UMI counts of the gene / sum over all pre-treatment
  cells (fine-type shares sum to 100%). Same % and mean over all cells (both timepoints) are given too. main_source = share >= 1%.
Output: supplementary/cn10_cytokine_celltype_screen.csv (3 significant figures; Supp Table 11).
Cache: supplementary/cache/cn10_cytokine_cells.parquet (per-cell OSM / IL1B / IL6 log-norm + counts, Hallmark Inflammatory Response score
  in All APCs); built if absent; read by the other Fig 6B sections and by 05_taurus_ligand_senders.py.
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, re, warnings, time, subprocess, tempfile, shutil
    warnings.filterwarnings('ignore')
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import numpy as np, pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'; R = '/share/fsmresfiles/UC/revision_jenny'
    TOP = f'{R}/Fig6_treg_antitnf_taurus'; HERE = f'{TOP}/supplementary'
    OUT = f'{HERE}/cn10_cytokine_celltype_screen.csv'; CACHE = f'{HERE}/cache/cn10_cytokine_cells.parquet'
    MSIG = f'{R}/cd_treg_rebuild_rowfix/archive/csv/msigdb_human_snapshot.csv.gz'
    KEYS = ['Patient', 'Remission_status', 'Treatment']; MIN_N = 10; MIN_PRE_PAT = 3
    GENES = ['OSM', 'IL1B', 'IL6']
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500); pd.set_option('display.max_columns', 50)

    SELM = ['Mac S+XS-', 'Mac S+M+S+', 'Mac S+SG+', 'Mac S+M+P+']; SELL = 'SELENOP+ Mac'
    DCT = ['cDC1', 'cDC2', 'CCR7+ DC']; APC = DCT + SELM + ['CD9+ Mac', 'Mo-Mac', 'B cell']       # 'All APCs (Pooled)' of F6B_inflammatory_signals_biopsy
    S_MHC2 = 'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION'; S_INFL = 'HALLMARK_INFLAMMATORY_RESPONSE'
    PERI = ['CD36+ Peri', 'RERGL+ Contr Peri']

    # ================================================================ per-cell values (cache)
    if not os.path.exists(CACHE):
        import scanpy as sc, decoupler as dc
        ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); log('loaded', ad.shape)
        assert ad.raw is None and np.allclose(np.expm1(ad.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)        # X = log1p(CP10k)
        assert ad.obs_names.is_unique and 'counts' in ad.layers and all(g in ad.var_names for g in GENES)
        gi = [ad.var_names.get_loc(g) for g in GENES]
        cnt = ad.layers['counts'][:, gi].toarray(); lg = ad.X[:, gi].toarray()
        assert np.allclose(cnt[:2000].astype(float), np.round(cnt[:2000])) and ((cnt > 0) == (lg > 0)).all()          # raw integer UMIs, same zeros
        assert np.allclose(np.asarray(ad.layers['counts'][:500].sum(1, dtype=np.float64)).ravel(), ad.obs.total_counts.values[:500], rtol=1e-3)
        C = ad.obs[['sample_id'] + KEYS + ['aucell_cell_type_short', 'aucell_cell_category_from_type']].astype(str).copy()
        for j, g in enumerate(GENES): C[g] = lg[:, j]; C[f'{g}_counts'] = cnt[:, j]
        # Hallmark Inflammatory Response in All APCs, as F6B_inflammatory_signals_biopsy
        msig = pd.read_csv(MSIG, usecols=['geneset', 'genesymbol'])
        sub = ad[np.isin(C.aucell_cell_type_short.values, APC)].copy(); sub.layers['lognorm'] = sub.X.copy()
        net = msig[msig.geneset.isin([S_MHC2, S_INFL])].drop_duplicates().rename(columns={'geneset': 'source', 'genesymbol': 'target'})
        dc.mt.ulm(data=sub, net=net, layer='lognorm'); A = dc.pp.get_obsm(sub, key='score_ulm')
        C['hallmark_inflammatory'] = np.nan
        C.loc[sub.obs_names, 'hallmark_inflammatory'] = np.asarray(A[:, S_INFL].X).ravel(); del sub, A, ad
        os.makedirs(os.path.dirname(CACHE), exist_ok=True); C.to_parquet(CACHE); log('cached', CACHE, C.shape)
    C = pd.read_parquet(CACHE); log('cells', C.shape)
    C['biopsy'] = C.sample_id
    assert C.groupby('biopsy').Patient.nunique().max() == 1 and C.groupby('biopsy').Treatment.nunique().max() == 1
    assert C.groupby('Patient').Remission_status.nunique().max() == 1
    assert set(C.Treatment) == {'Pre', 'Post'} and set(C.Remission_status) == {'Remission', 'Non_Remission'}
    C['ct'] = np.where(C.aucell_cell_type_short.isin(SELM), SELL, C.aucell_cell_type_short)
    cat = C.aucell_cell_category_from_type

    # ================================================================ groups
    FINE = sorted(C.ct.unique()); assert len(FINE) == 54 - 4 + 1, len(FINE)
    LIN = {'Myeloid (All)': cat.isin(['Myeloid', 'Granulocyte']).values,
           'Fibroblasts (All)': ((cat == 'Connective') & ~C.ct.isin(PERI)).values,
           'T Cells (All)': (cat == 'T').values, 'B / Plasma (All)': cat.isin(['B', 'Plasma']).values,
           'Epithelial (All)': (cat == 'Epithelial').values, 'Endothelial (All)': (C.ct == 'Endothelial').values}
    assert sorted(C.ct[LIN['Fibroblasts (All)']].unique()) == sorted(t for t in FINE if 'fib' in t.lower()), C.ct[LIN['Fibroblasts (All)']].unique()
    GDEF = [(t, 'Fine cell type', (C.ct == t).values) for t in FINE] + [(k, 'Lineage', m) for k, m in LIN.items()] \
           + [('All Cells', 'All cells (reference)', np.ones(len(C), bool))]
    LIN_MEMBERS = {k: '; '.join(sorted(C.ct[m].unique())) for k, m in LIN.items()}
    for k, v in LIN_MEMBERS.items(): log(k, ':', v)

    GRP = [('rem_pre', 'Remission', 'Pre'), ('rem_post', 'Remission', 'Post'), ('nonrem_pre', 'Non_Remission', 'Pre'), ('nonrem_post', 'Non_Remission', 'Post')]
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}

    def pseudobulk(D, col):
        return D.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'),
                                       Treatment=('Treatment', 'first'), n_cells=(col, 'size'), y=(col, 'mean')).reset_index()

    GROWS, LONG, PBS = [], [], {}
    for name, level, m in GDEF:
        D = C[m]
        nb = D.groupby('biopsy').size(); bmeta = D.groupby('biopsy')[['Patient', 'Remission_status', 'Treatment']].first().assign(n=nb)
        keep = bmeta[bmeta.n >= MIN_N]
        row = dict(cell_group=name, level=level, n_cells_total=len(D), n_biopsies_total=len(bmeta), n_biopsies_kept=len(keep),
                   n_cells_kept=int(keep.n.sum()), n_biopsies_dropped=len(bmeta) - len(keep), n_patients_kept=keep.Patient.nunique())
        for nm, r_, t_ in GRP:
            s = keep[(keep.Remission_status == r_) & (keep.Treatment == t_)]
            row.update({f'n_patients_{nm}': s.Patient.nunique(), f'n_biopsies_{nm}': len(s), f'n_cells_{nm}': int(s.n.sum())})
        elig = row['n_patients_rem_pre'] >= MIN_PRE_PAT and row['n_patients_nonrem_pre'] >= MIN_PRE_PAT
        post_ok = row['n_biopsies_rem_post'] >= 1 and row['n_biopsies_nonrem_post'] >= 1
        if not elig:
            reason = (f"< {MIN_PRE_PAT} patients with an eligible (>= {MIN_N} cells) pre-treatment biopsy in "
                      + ' and '.join([lab for lab, k in [('Remission', 'n_patients_rem_pre'), ('Non-remission', 'n_patients_nonrem_pre')] if row[k] < MIN_PRE_PAT])
                      + f" (Remission {row['n_patients_rem_pre']}, Non-remission {row['n_patients_nonrem_pre']})")
        elif not post_ok:
            reason = 'no eligible post-treatment biopsy in ' + ' and '.join([lab for lab, k in [('Remission', 'n_biopsies_rem_post'), ('Non-remission', 'n_biopsies_nonrem_post')] if row[k] < 1])
        else:
            reason = ''
        row['eligible'] = elig and post_ok; row['exclusion_reason'] = reason
        GROWS.append(row)
        Dp = D[D.Treatment == 'Pre']; Dk = D[D.biopsy.isin(keep.index)]
        for g in GENES:
            e = dict(cell_group=name, gene=g, pct_cells_expressing_pre=100 * (Dp[g] > 0).mean() if len(Dp) else np.nan,
                     mean_expression_pre=Dp[g].mean() if len(Dp) else np.nan, pct_cells_expressing_all=100 * (D[g] > 0).mean(),
                     mean_expression_all=D[g].mean(), share_of_transcripts_pre=100 * Dp[f'{g}_counts'].sum() / TOT_PRE[g],
                     pct_cells_expressing_kept=100 * (Dk[g] > 0).mean() if len(Dk) else np.nan)
            if row['eligible']:
                PB = pseudobulk(Dk, g); assert len(PB) == len(keep)
                for nm, r_, t_ in GRP:
                    s = PB[(PB.Remission_status == r_) & (PB.Treatment == t_)].y
                    e[f'mean_of_biopsy_means_{nm}'] = s.mean(); e[f'sem_of_biopsy_means_{nm}'] = s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else np.nan
                e['estimable'] = bool(PB.y.std() > 0 and (Dk[g] > 0).any())
                if e['estimable']:
                    k = f'{name}||{g}'; PBS[k] = PB
                    LONG.append(PB.assign(feature_key=k, pt=PB.Patient + '_' + PB.Treatment, post=(PB.Treatment == 'Post').astype(float),
                                          rem=(PB.Remission_status == 'Remission').astype(float))[['feature_key', 'Patient', 'pt', 'post', 'rem', 'y']])
                else:
                    e['exclusion_reason_gene'] = f'{g} not expressed in the kept biopsies'
            else:
                e['estimable'] = False
            GROWS[-1].setdefault('_genes', []).append(e)
    G = pd.DataFrame([{k: v for k, v in r.items() if k != '_genes'} for r in GROWS])
    E = pd.DataFrame([e for r in GROWS for e in r['_genes']])
    # main source of a gene (descriptive flag): >= 1% of the gene's pre-treatment transcripts
    SRC_MIN = 1.0
    E['main_source'] = E.share_of_transcripts_pre >= SRC_MIN
    fine_tot = E[E.cell_group.isin(FINE)].groupby('gene').share_of_transcripts_pre.sum()
    assert np.allclose(fine_tot.values, 100), fine_tot
    LONG = pd.concat(LONG, ignore_index=True)
    log(f"eligible fine types {int(G[(G.level == 'Fine cell type')].eligible.sum())} / {len(FINE)}; lineages {int(G[G.level == 'Lineage'].eligible.sum())} / {len(LIN)}")
    print(G[~G.eligible][['cell_group', 'level', 'n_cells_total', 'n_biopsies_kept', 'n_patients_rem_pre', 'n_patients_nonrem_pre',
                          'n_biopsies_rem_post', 'n_biopsies_nonrem_post', 'exclusion_reason']].to_string(index=False), flush=True)
    print(G[G.cell_group == 'OGN+RSPO3+ Fib'].T.to_string(), flush=True)

    # ================================================================ lme4 / lmerTest, Satterthwaite (+ estimate_sd)
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
        log(f'fitting {LONG.feature_key.nunique()} cell type x gene models in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/r_out.csv'], check=True)
        RR = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(RR) == 4 * LONG.feature_key.nunique() and (RR.p.notna() | (RR.optimizer == 'fit failed')).all()
    log('optimizer use:', RR.drop_duplicates('feature_key').optimizer.value_counts().to_dict(), 'fit failed:', RR[RR.optimizer == 'fit failed'].feature_key.unique().tolist())
    log('fitted')
    CON = {'baseline_difference': 'Baseline difference (remission − non-remission, pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)'}
    RR['cell_group'] = RR.feature_key.str.split(r'\|\|').str[0]; RR['gene'] = RR.feature_key.str.split(r'\|\|').str[1]
    RR['contrast'] = RR.contrast_key.map(CON)
    tq = stats.t.ppf(0.975, RR.df); RR['ci_low'] = RR.estimate - tq * RR.se; RR['ci_high'] = RR.estimate + tq * RR.se

    # full table: every cell group x gene x contrast (excluded / not estimable groups keep one row per contrast with empty statistics)
    base = E.merge(G, on='cell_group')
    full = base.assign(_k=1).merge(pd.DataFrame({'contrast_key': list(CON), '_k': 1}), on='_k').drop(columns='_k')
    full['contrast'] = full.contrast_key.map(CON)
    RR = full.merge(RR.drop(columns=['contrast', 'feature_key']), on=['cell_group', 'gene', 'contrast_key'], how='left')
    RR['bh_family'] = RR.level.map({'Fine cell type': 'fine cell types', 'Lineage': 'lineage groups', 'All cells (reference)': 'none'})
    RR['q_bh'] = np.nan
    for (fam, g, k), i in RR[RR.p.notna() & (RR.bh_family != 'none')].groupby(['bh_family', 'gene', 'contrast_key']).groups.items():
        RR.loc[i, 'q_bh'] = multipletests(RR.loc[i, 'p'], method='fdr_bh')[1]
    nfam = RR[RR.p.notna()].groupby(['bh_family', 'gene', 'contrast_key']).size()
    RR['n_tests_in_bh_family'] = [nfam.get((f, g, k), np.nan) if f != 'none' else np.nan for f, g, k in RR[['bh_family', 'gene', 'contrast_key']].values]
    RR['bh_family'] = np.where(RR.bh_family == 'none', 'not in a BH family (reference)', 'BH across ' + RR.bh_family + ', per gene and contrast')
    RR['lineage_members'] = RR.cell_group.map(LIN_MEMBERS).fillna(RR.cell_group.map({SELL: '; '.join(SELM), 'All Cells': 'every cell'}))
    RR['model_formula'] = 'biopsy mean ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df; no covariates'
    RR['unit'] = 'log-normalised expression, log1p(CP10k); estimate_sd = SD units of the biopsy means'
    RR['eligibility_rule'] = (f'biopsies with < {MIN_N} cells of the type dropped; >= {MIN_PRE_PAT} patients with an eligible pre-treatment biopsy in each '
                          'remission group; >= 1 eligible post-treatment biopsy per group (model estimability)')

    # row order: main sources first = sum over the 3 genes of the share of pre-treatment transcripts
    ORD = E.groupby('cell_group').share_of_transcripts_pre.sum()
    LV = {'Fine cell type': 0, 'Lineage': 1, 'All cells (reference)': 2}
    RR['_lv'] = RR.level.map(LV); RR['_o'] = -RR.cell_group.map(ORD); RR['_g'] = RR.gene.map({g: i for i, g in enumerate(GENES)})
    RR['_c'] = RR.contrast_key.map({k: i for i, k in enumerate(CON)})
    RR['source_rank_score'] = RR.cell_group.map(ORD)
    RR = RR.sort_values(['_lv', 'eligible', '_o', '_g', '_c'], ascending=[True, False, True, True, True]).reset_index(drop=True)
    RR['estimable'] = RR.estimable.astype(bool)

    COLS = (['level', 'cell_group', 'lineage_members', 'gene', 'contrast', 'estimate', 'ci_low', 'ci_high', 'se', 'estimate_sd', 'df', 'statistic',
             'p', 'q_bh', 'bh_family', 'n_tests_in_bh_family', 'eligible', 'estimable', 'exclusion_reason',
             'pct_cells_expressing_pre', 'mean_expression_pre', 'share_of_transcripts_pre', 'pct_cells_expressing_all', 'mean_expression_all',
             'pct_cells_expressing_kept', 'main_source', 'source_rank_score']
            + [f'{a}_{nm}' for a in ['mean_of_biopsy_means', 'sem_of_biopsy_means'] for nm, _, _ in GRP]
            + ['n_cells_total', 'n_biopsies_total', 'n_biopsies_kept', 'n_cells_kept', 'n_biopsies_dropped', 'n_patients_kept']
            + [f'{a}_{nm}' for nm, _, _ in GRP for a in ['n_patients', 'n_biopsies', 'n_cells']]
            + ['var_patient', 'var_patient_timepoint', 'var_residual', 'singular', 'converged', 'optimizer', 'unit', 'model_formula', 'eligibility_rule'])
    if 'exclusion_reason_gene' in RR:
        RR['exclusion_reason'] = np.where(RR.exclusion_reason_gene.notna() & (RR.exclusion_reason == ''), RR.exclusion_reason_gene, RR.exclusion_reason)
    RR['exclusion_reason'] = np.where((RR.optimizer == 'fit failed') & (RR.exclusion_reason == ''), 'model fit failed (both optimizers)', RR.exclusion_reason)
    RES = RR[COLS].copy()
    with open(OUT, 'w') as fh:
        RES.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, RES.shape)

    # ================================================================ console summary
    S = RR[RR.p.notna()]
    for g in GENES:
        src = E[(E.gene == g) & E.cell_group.isin(FINE)].sort_values('share_of_transcripts_pre', ascending=False).head(6)
        print(f'\n{g} main sources (fine types, pre-treatment): ' + '; '.join(f'{r.cell_group} share {r.share_of_transcripts_pre:.1f}% pct {r.pct_cells_expressing_pre:.1f}% mean {r.mean_expression_pre:.3f}' for r in src.itertuples()))
        for lvl in ['Fine cell type', 'Lineage', 'All cells (reference)']:
            for k in ['baseline_difference', 'differential_change']:
                s = S[(S.gene == g) & (S.level == lvl) & (S.contrast_key == k)].sort_values('p').head(4)
                print(f'  {lvl:22s} {k:20s}: ' + '; '.join(f'{r.cell_group} est {r.estimate:+.3g} (sd {r.estimate_sd:+.2f}) P {r.p:.3g} q {r.q_bh:.3g}' for r in s.itertuples()))
    TOPB, TOPS = {}, {}
    for g in GENES:
        s = S[(S.gene == g) & (S.level == 'Fine cell type') & (S.contrast_key == 'baseline_difference')].sort_values('p')
        TOPB[g] = s.iloc[0].cell_group; TOPS[g] = s[s.main_source.astype(bool)].iloc[0].cell_group
    log('top baseline fine type per gene (any type):', TOPB, '| among main sources (>= 1% of transcripts):', TOPS)


# ======================================================================================================================
# Section F6B_coarse_screen: coarse cell-group screen of OSM, IL1B and IL6 with abundance, and decomposition of the all-cell change
# Panels: Fig 6B dots and P values (Supp Table 11)
# ======================================================================================================================
if _run('F6B_coarse_screen'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/57b_cn10_cytokine_coarse_screen.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: coarse-group screen of the CN10 cytokines OSM, IL1B and IL6, with abundance and a decomposition of
the all-cell change.
Question: is the all-cell persistence of OSM / IL1B / IL6 in non-remitters carried by more of the source cells (abundance) or by more
expression per cell?

Loading, cache, eligibility rule and model (R code) as cytokine_celltype_screen.
Coarse groups (aucell_cell_type_short -> group; mapping table COARSE below, written into the CSV and legend):
  DC = cDC1, cDC2, CCR7+ DC, Trans cDC2/Mac; Macrophage = Mac S+XS-, Mac S+M+S+, Mac S+SG+, Mac S+M+P+, CD9+ Mac;
  Monocyte-derived Mac = Mo-Mac, Inf Mo-Mac; Granulocyte/Neutrophil = Granulocyte, Neutrophil; B cell; Plasma;
  CD4 T = CD4 Effector, CD4 Trm(Prolif), CD4 Trm+D42, Exhaus CD4 Tfh, Tfh, Treg; CD8 T = CD8 GZMK+ Trm, CD8 GZMK- Trm, CD8 GZMK- Trm(Prolif);
  Other T/NK = Gd T, Gd T (Vd2g9+), NK, Resting T; Fibroblast = the 10 fibroblast / myofibroblast types; Pericyte = CD36+ Peri,
  RERGL+ Contr Peri; Endothelial; Epithelial = the 11 epithelial types + Cycling; All cells (reference).
  'Cycling' was to join Other T/NK only if T: it is checked on the data (cache/cn10_cycling_lineage_check.csv, built once from the h5ad):
  if < 5% of Cycling cells express CD3D/CD3E and > 50% express EPCAM, Cycling is epithelial and joins Epithelial (asserted below).
  Neuron is in no coarse group (All cells only; kept as its own term in the decomposition so the terms add up exactly).
Part 1 (expression): y = per biopsy, mean log1p(CP10k) of the gene over the group's cells; biopsies with < 10 cells of the group dropped.
Part 2 (abundance): y = per biopsy, % of all the biopsy's cells that belong to the group; all biopsies kept (a low count is the readout).
  Groups are tested in both parts only if eligible by the cytokine_celltype_screen rule: >= 3 Remission and >= 3 Non-remission patients with a pre-treatment
  biopsy holding >= 10 cells of the group, and >= 1 such post-treatment biopsy per remission group.
Model (both parts): y ~ post * remission + (1 | patient) + (1 | patient:timepoint); lme4 + lmerTest REML, Satterthwaite df; no covariates;
  y standardised before fitting (estimate_sd = SD units of the biopsy means; estimate = back-transformed units); bobyqa retry if the default
  optimizer errors. Contrasts: baseline difference = b_rem; change in remission = b_post + b_post:rem; change in non-remission = b_post;
  differential change = b_post:rem. BH across the 13 coarse groups, per feature (gene or abundance) and contrast; All cells not in the family.
Share of transcripts (pre-treatment): the group's sum of raw UMI counts of the gene / the sum over all pre-treatment cells (13 groups +
  Neuron = 100%).
Part 3 (decomposition, observed values, all 96 biopsies): per biopsy b, all-cell mean M_b = sum_g f_gb * m_gb exactly (f = group fraction of
  the biopsy's cells, m = within-group mean; f * m = 0 when the biopsy has no cell of g). Per remission group and timepoint, the mean over
  biopsies of f_gb * m_gb is group g's contribution C_g; sum_g C_g = mean of the biopsy all-cell means. Contribution to the change =
  C_g(post) - C_g(pre); these add up to the observed all-cell change. Each group's change is split exactly into an abundance part and an
  expression part (two-factor midpoint split): C = F * W, F = mean fraction, W = C / F (fraction-weighted within-group mean);
  dC = dF * (W_pre + W_post)/2 + dW * (F_pre + F_post)/2.
Outputs: supplementary/cn10_cytokine_coarse_screen.csv (3 significant figures; Supp Table 11, Fig 6B source data) and
  supplementary/cache/cn10_cycling_lineage_check.csv.
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, re, warnings, time, subprocess, tempfile, shutil
    warnings.filterwarnings('ignore')
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import numpy as np, pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'; R = '/share/fsmresfiles/UC/revision_jenny'
    TOP = f'{R}/Fig6_treg_antitnf_taurus'; HERE = f'{TOP}/supplementary'
    OUT = f'{HERE}/cn10_cytokine_coarse_screen.csv'; CACHE = f'{HERE}/cache/cn10_cytokine_cells.parquet'
    CYC = f'{HERE}/cache/cn10_cycling_lineage_check.csv'
    KEYS = ['Patient', 'Remission_status', 'Treatment']; MIN_N = 10; MIN_PRE_PAT = 3
    GENES = ['OSM', 'IL1B', 'IL6']; AB = 'Abundance'; FEATS = GENES + [AB]
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500); pd.set_option('display.max_columns', 50)
    assert os.path.exists(CACHE), f'{CACHE} missing: run cytokine_celltype_screen first (it builds the per-cell cache)'

    # ================================================================ coarse mapping
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

    # ---- Cycling lineage check (data, not assumption); built once from the h5ad, then cached
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
    CY = pd.read_csv(CYC).set_index('item').value
    cyc_t = max(CY['% of Cycling cells expressing CD3D'], CY['% of Cycling cells expressing CD3E']); cyc_e = CY['% of Cycling cells expressing EPCAM']
    cyc_epi_study = CY[[k for k in CY.index if 'epithelium' in k]].sum()
    assert cyc_t < 5 and cyc_e > 50, (cyc_t, cyc_e)
    log(f'Cycling: CD3D/CD3E+ {cyc_t:.1f}%, EPCAM+ {cyc_e:.1f}%, TAURUS study annotation epithelium {cyc_epi_study:.1f}% -> Epithelial')

    C = pd.read_parquet(CACHE); log('cells', C.shape)
    C['biopsy'] = C.sample_id
    assert C.groupby('biopsy').Patient.nunique().max() == 1 and C.groupby('biopsy').Treatment.nunique().max() == 1
    assert C.groupby('Patient').Remission_status.nunique().max() == 1
    assert set(C.Treatment) == {'Pre', 'Post'} and set(C.Remission_status) == {'Remission', 'Non_Remission'}
    T2G = {t: g for g, ts in COARSE.items() for t in ts}
    assert len(T2G) == sum(len(v) for v in COARSE.values()), 'a type is in two groups'
    types = set(C.aucell_cell_type_short)
    assert set(T2G) <= types, set(T2G) - types
    assert types - set(T2G) == {'Neuron'}, types - set(T2G)
    assert sorted(FIBS) == sorted(t for t in types if 'fib' in t.lower())
    assert sorted(EPIS) == sorted(C.aucell_cell_type_short[C.aucell_cell_category_from_type == 'Epithelial'].unique())
    C['grp'] = C.aucell_cell_type_short.map(T2G).fillna(NOGRP)
    GROUPS = list(COARSE)
    MEMBERS = {g: '; '.join(v) for g, v in COARSE.items()}; MEMBERS[ALL] = 'every cell'; MEMBERS[NOGRP] = 'Neuron'

    GRP = [('rem_pre', 'Remission', 'Pre'), ('rem_post', 'Remission', 'Post'), ('nonrem_pre', 'Non_Remission', 'Pre'), ('nonrem_post', 'Non_Remission', 'Post')]
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}
    BM = C.groupby('biopsy')[KEYS].first().assign(n_all=C.groupby('biopsy').size())
    log('biopsies', len(BM), 'cells per biopsy min', BM.n_all.min())
    assert BM.n_all.min() >= MIN_N

    def pseudobulk(D, col):
        return D.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'),
                                       Treatment=('Treatment', 'first'), n_cells=(col, 'size'), y=(col, 'mean')).reset_index()

    # ================================================================ parts 1 + 2: per-group tables and model input
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
    # Neuron share of transcripts (for the 100% check and the CSV)
    for f in GENES:
        EROWS.append(dict(cell_group=NOGRP, feature=f, share_of_transcripts_pre=100 * C.loc[pre & (C.grp == NOGRP).values, f'{f}_counts'].sum() / TOT_PRE[f],
                          pct_cells_expressing_pre=100 * (C.loc[pre & (C.grp == NOGRP).values, f] > 0).mean(), estimable=False,
                          exclusion_reason_feature='not a coarse group (decomposition term only)'))
    G = pd.DataFrame(GROWS); E = pd.DataFrame(EROWS)
    tot = E[E.cell_group.isin(GROUPS + [NOGRP]) & E.feature.isin(GENES)].groupby('feature').share_of_transcripts_pre.sum()
    assert np.allclose(tot.values, 100), tot
    LONG = pd.concat(LONG, ignore_index=True)
    print(G[['cell_group', 'n_cells_total', 'pct_of_all_cells_pre', 'n_biopsies_kept', 'n_patients_rem_pre', 'n_patients_nonrem_pre',
             'n_biopsies_rem_post', 'n_biopsies_nonrem_post', 'eligible']].to_string(index=False), flush=True)

    # ================================================================ lme4 / lmerTest (as cytokine_celltype_screen)
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
        log(f'fitting {LONG.feature_key.nunique()} group x feature models in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/r_out.csv'], check=True)
        RR = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(RR) == 4 * LONG.feature_key.nunique() and (RR.p.notna() | (RR.optimizer == 'fit failed')).all()
    log('optimizer use:', RR.drop_duplicates('feature_key').optimizer.value_counts().to_dict())
    CON = {'baseline_difference': 'Baseline difference (remission − non-remission, pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)'}
    RR['cell_group'] = RR.feature_key.str.split(r'\|\|').str[0]; RR['feature'] = RR.feature_key.str.split(r'\|\|').str[1]
    tq = stats.t.ppf(0.975, RR.df); RR['ci_low'] = RR.estimate - tq * RR.se; RR['ci_high'] = RR.estimate + tq * RR.se
    base = E[E.cell_group != NOGRP].merge(G, on='cell_group')
    full = base.assign(_k=1).merge(pd.DataFrame({'contrast_key': list(CON), '_k': 1}), on='_k').drop(columns='_k')
    RR = full.merge(RR.drop(columns=['feature_key']), on=['cell_group', 'feature', 'contrast_key'], how='left')
    RR['contrast'] = RR.contrast_key.map(CON)
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
    S = RR[RR.p.notna()]
    for f in FEATS:
        for k in ['differential_change', 'baseline_difference', 'change_remission', 'change_nonremission']:
            s = S[(S.feature == f) & (S.contrast_key == k)].sort_values('p').head(5)
            print(f'{f:9s} {k:20s}: ' + '; '.join(f'{r.cell_group} {r.estimate:+.3g} (sd {r.estimate_sd:+.2f}) P {r.p:.3g} q {r.q_bh:.3g}' for r in s.itertuples()))

    # ================================================================ part 3: decomposition of the all-cell change (observed biopsy values)
    DEC = []
    cnt = C.groupby(['biopsy', 'grp']).size().unstack(fill_value=0).reindex(BM.index)
    frac = cnt.div(BM.n_all, axis=0)                                         # f_gb
    for f in GENES:
        mg = C.groupby(['biopsy', 'grp'])[f].mean().unstack().reindex(BM.index)  # m_gb (NaN if no cell)
        contrib = (frac * mg).fillna(0)                                        # f_gb * m_gb
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
                DEC.append(dict(feature=f, cell_group=g, remission_group=rlab, contribution_pre=C0, contribution_post=C1, contribution_to_change=C1 - C0,
                                pct_of_all_cell_change=100 * (C1 - C0) / dtot, abundance_part=dab, expression_part=dex,
                                mean_fraction_pre=100 * F0, mean_fraction_post=100 * F1, weighted_group_mean_pre=W0, weighted_group_mean_post=W1,
                                all_cell_mean_pre=tot0, all_cell_mean_post=tot1, all_cell_change=dtot, n_biopsies_pre=len(b0), n_biopsies_post=len(b1)))
    DEC = pd.DataFrame(DEC)
    chk = DEC.groupby(['feature', 'remission_group']).agg(s=('contribution_to_change', 'sum'), t=('all_cell_change', 'first'))
    assert np.allclose(chk.s, chk.t, atol=1e-12), chk
    # all-cell observed means reproduce the expression table (All cells row)
    for f in GENES:
        a_ = RR[(RR.cell_group == ALL) & (RR.feature == f)].iloc[0]
        for rl, nm in [('remission', 'rem'), ('non-remission', 'nonrem')]:
            d_ = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].iloc[0]
            assert np.isclose(d_.all_cell_mean_pre, a_[f'mean_of_biopsy_means_{nm}_pre']) and np.isclose(d_.all_cell_mean_post, a_[f'mean_of_biopsy_means_{nm}_post'])
    # differential (remission - non-remission) rows
    W = DEC.pivot_table(index=['feature', 'cell_group'], columns='remission_group',
                        values=['contribution_to_change', 'abundance_part', 'expression_part', 'all_cell_change'])
    DD = pd.DataFrame({'contribution_to_change': W['contribution_to_change']['remission'] - W['contribution_to_change']['non-remission'],
                       'abundance_part': W['abundance_part']['remission'] - W['abundance_part']['non-remission'],
                       'expression_part': W['expression_part']['remission'] - W['expression_part']['non-remission'],
                       'all_cell_change': W['all_cell_change']['remission'] - W['all_cell_change']['non-remission']}).reset_index()
    DD['pct_of_all_cell_change'] = 100 * DD.contribution_to_change / DD.all_cell_change
    DD['remission_group'] = 'differential (remission − non-remission)'
    DEC = pd.concat([DEC, DD], ignore_index=True)
    DEC['part'] = '3 decomposition'
    for f in GENES:
        print(f'\n{f} decomposition of the all-cell change (observed):')
        s = DEC[DEC.feature == f].pivot_table(index='cell_group', columns='remission_group', values='contribution_to_change')
        s = s.reindex(s['differential (remission − non-remission)'].abs().sort_values(ascending=False).index)
        print(s.round(5).to_string())
        print(DEC[(DEC.feature == f)].groupby('remission_group')[['abundance_part', 'expression_part', 'contribution_to_change']].sum().round(5).to_string())

    # ================================================================ CSV
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
    with open(OUT, 'w') as fh:
        RES.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, RES.shape)


# ======================================================================================================================
# Section F6B_contribution_se: model-based SE and bootstrap CI of the net all-cell % change of OSM, IL1B and IL6 (Fig 6B bars)
# Panels: Fig 6B bars (source data cn10_cytokine_contrib_bootstrap.csv)
# ======================================================================================================================
if _run('F6B_contribution_se'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/57e_fig6_cn10_cytokine_contrib.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: error bars and source data of the Fig 6B stacked bars (no model is refitted).
Inputs:
  supplementary/cn10_cytokine_coarse_screen.csv   (F6B_coarse_screen: All cells change / difference P from the biopsy-level mixed model)
  supplementary/nonremission_persistent_signals.csv (the same All Cells P values, estimates and CIs; both CSVs must agree)
  supplementary/cache/cn10_cytokine_cells.parquet  (per-cell cache; stacked segments as in F6B_coarse_screen)
  supplementary/cache/cn10_cycling_lineage_check.csv
Net % change (diamond): change in the all-cell mean (post - pre, mean over biopsies) / the group's observed pre-treatment all-cell mean.
Error bar = +/- 1 SE of the biopsy-level mixed model: All cells 'change in remission' / 'change in non-remission' estimate and SE
  (log-normalised units) divided by the group's observed pre-treatment all-cell mean (constant denominator, the diamond's own), centred on
  the model estimate % (|model - observed| <= 3 percentage points asserted). The model 95% CI is in the CSV (model_ci95_*); it excludes 0
  exactly when change P < 0.05 (asserted).
Bootstrap (CSV columns 'bootstrap_ratio_*' only, not plotted): patient-level cluster bootstrap, patients resampled with replacement within
  each remission group (2,000 iterations, seed 20261105); each replicate recomputes the diamond: per-biopsy all-cell mean, mean over the
  replicate's pre and post biopsies (a patient drawn k times contributes its biopsies k times), 100 * (post - pre) / pre. All biopsies of each
  remission group are used (not paired-only): Remission 10 patients (all paired); Non-remission 6 patients pre, 5 post. A replicate with no
  pre- or no post-treatment biopsy is redrawn (counted in the CSV). The same patient draws are used for all three genes. Percentile CI.
Output: supplementary/cn10_cytokine_contrib_bootstrap.csv (read by F6B_cytokine_sources).
Run: OPENBLAS_CORETYPE=SkylakeX /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, re, math, warnings
    warnings.filterwarnings('ignore')
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import numpy as np, pandas as pd

    R = '/share/fsmresfiles/UC/revision_jenny'; TOP = f'{R}/Fig6_treg_antitnf_taurus'; HERE = f'{TOP}/supplementary'
    CSV = f'{HERE}/cn10_cytokine_coarse_screen.csv'; NRP = f'{HERE}/nonremission_persistent_signals.csv'
    CACHE = f'{HERE}/cache/cn10_cytokine_cells.parquet'; CYC = f'{HERE}/cache/cn10_cycling_lineage_check.csv'
    BOOT = f'{HERE}/cn10_cytokine_contrib_bootstrap.csv'
    GENES = ['OSM', 'IL1B', 'IL6']; N_BOOT = 2000; SEED = 20261105
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
    def p2(p):                                                                    # 2 significant figures, trailing zeros kept (0.0020)
        d = max(1 - math.floor(math.log10(abs(p))), 0); return f'{p:.{d}f}'

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
    NOGRP = 'Neuron (no coarse group)'; ALL = 'All cells'; GROUPS = list(COARSE)
    NAME = {'Monocyte-derived Mac': 'Mo-Mac'}                                   # display name only (CSV keeps the analysis name)
    disp = lambda s: re.sub(r'\+(?=[\s)A-Z]|$)', '⁺', NAME.get(s, s))
    CHG = {'remission': 'Change in remission (post − pre)', 'non-remission': 'Change in non-remission (post − pre)'}
    DIFC = 'Differential change (remission − non-remission)'

    # ================================================================ P values (All cells, biopsy-level mixed model; F6B_coarse_screen CSV)
    RES = pd.read_csv(CSV)
    A = RES[(RES.part == '1 expression') & (RES.cell_group == ALL) & RES.feature.isin(GENES)].set_index(['feature', 'contrast'])
    NP = pd.read_csv(NRP); NP = NP[(NP.cell_group == 'All Cells') & NP.gene_or_program.isin(GENES)].set_index(['gene_or_program', 'contrast'])
    PV = {}
    for f in GENES:
        for k, c in list(CHG.items()) + [('difference', DIFC)]:
            p = float(A.loc[(f, c), 'p']); assert np.isclose(p, NP.loc[(f, c), 'p'], rtol=1e-6), (f, c)   # both CSVs agree
            PV[(f, k)] = p
    print('P values', {k: p2(v) for k, v in PV.items()})

    # ================================================================ decomposition (as F6B_coarse_screen)
    G = RES[RES.part.isin(['1 expression', '2 abundance']) & RES.contrast.notna()].drop_duplicates('cell_group').set_index('cell_group')
    CY = pd.read_csv(CYC).set_index('item').value
    cyc_t = max(CY['% of Cycling cells expressing CD3D'], CY['% of Cycling cells expressing CD3E']); cyc_e = CY['% of Cycling cells expressing EPCAM']
    C = pd.read_parquet(CACHE); C['biopsy'] = C.sample_id
    C[GENES] = C[GENES].astype('float64')
    T2G = {t: g for g, ts in COARSE.items() for t in ts}
    assert set(C.aucell_cell_type_short) - set(T2G) == {'Neuron'}
    C['grp'] = C.aucell_cell_type_short.map(T2G).fillna(NOGRP)
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}
    SHX = {(g, f): 100 * C.loc[pre & (C.grp == g).values, f'{f}_counts'].sum() / TOT_PRE[f] for g in GROUPS for f in GENES}
    BM = C.groupby('biopsy')[['Patient', 'Remission_status', 'Treatment']].first().assign(n_all=C.groupby('biopsy').size())
    cnt = C.groupby(['biopsy', 'grp']).size().unstack(fill_value=0).reindex(BM.index)
    frac = cnt.div(BM.n_all, axis=0)
    DEC = []; ALLM = {}
    for f in GENES:
        mg = C.groupby(['biopsy', 'grp'])[f].mean().unstack().reindex(BM.index)
        contrib = (frac * mg).fillna(0); allm = C.groupby('biopsy')[f].mean().reindex(BM.index); ALLM[f] = allm
        assert np.allclose(contrib.sum(1), allm, atol=1e-10)
        for rs, rlab in [('Remission', 'remission'), ('Non_Remission', 'non-remission')]:
            b0 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Pre')]; b1 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Post')]
            tot0, tot1 = allm[b0].mean(), allm[b1].mean()
            for g in GROUPS + [NOGRP]:
                C0, C1 = contrib.loc[b0, g].mean(), contrib.loc[b1, g].mean()
                DEC.append(dict(feature=f, cell_group=g, remission_group=rlab, contribution_to_change=C1 - C0,
                                all_cell_mean_pre=tot0, all_cell_change=tot1 - tot0))
    DEC = pd.DataFrame(DEC)
    DEC['pct_contribution'] = 100 * DEC.contribution_to_change / DEC.all_cell_mean_pre
    DEC['pct_net'] = 100 * DEC.all_cell_change / DEC.all_cell_mean_pre
    for (f, rl), d in DEC.groupby(['feature', 'remission_group']):
        assert np.isclose(d.pct_contribution.sum(), d.pct_net.iloc[0], rtol=1e-12, atol=1e-10), (f, rl)
    CD = RES[(RES.part == '3 decomposition') & RES.remission_group.isin(['remission', 'non-remission'])]
    chk = DEC.merge(CD[['feature', 'cell_group', 'remission_group', 'contribution_to_change', 'all_cell_mean_pre', 'all_cell_change']],
                    on=['feature', 'cell_group', 'remission_group'], suffixes=('', '_csv'))
    assert len(chk) == len(DEC) == len(CD)
    for c in ['contribution_to_change', 'all_cell_mean_pre', 'all_cell_change']: assert np.allclose(chk[c], chk[f'{c}_csv'], rtol=5e-3, atol=1e-12), c
    DIFF = RES[(RES.part == '3 decomposition') & RES.remission_group.astype(str).str.startswith('differential')]
    SPLIT = {f: 100 * DIFF[DIFF.feature == f].abundance_part.sum() / DIFF[DIFF.feature == f].contribution_to_change.sum() for f in GENES}
    print('decomposition matches CSV; split', {f: round(v, 1) for f, v in SPLIT.items()})

    # ================================================================ patient-level cluster bootstrap of the net % change (diamond)
    rng = np.random.default_rng(SEED); BS = []; CI = {}
    for rs, rlab in [('Remission', 'remission'), ('Non_Remission', 'non-remission')]:
        bm = BM[BM.Remission_status == rs]; pats = np.array(sorted(bm.Patient.unique()))
        is0 = (bm.Treatment == 'Pre').values; is1 = (bm.Treatment == 'Post').values
        pidx = pd.Index(pats).get_indexer(bm.Patient); assert (pidx >= 0).all()
        n0 = np.bincount(pidx, weights=is0, minlength=len(pats)); n1 = np.bincount(pidx, weights=is1, minlength=len(pats))
        S0 = {f: np.bincount(pidx, weights=ALLM[f][bm.index].values * is0, minlength=len(pats)) for f in GENES}
        S1 = {f: np.bincount(pidx, weights=ALLM[f][bm.index].values * is1, minlength=len(pats)) for f in GENES}
        W = []; redraw = 0
        while len(W) < N_BOOT:
            w = np.bincount(rng.integers(0, len(pats), len(pats)), minlength=len(pats))
            if (w @ n0) == 0 or (w @ n1) == 0: redraw += 1; continue
            W.append(w)
        W = np.array(W, dtype=float)
        for f in GENES:
            obs = 100 * ((S1[f].sum() / n1.sum()) - (S0[f].sum() / n0.sum())) / (S0[f].sum() / n0.sum())
            net_d = DEC[(DEC.feature == f) & (DEC.remission_group == rlab)].pct_net.iloc[0]
            assert np.isclose(obs, net_d, rtol=1e-12, atol=1e-10), (f, rlab, obs, net_d)     # bootstrap statistic == diamond
            m0 = (W @ S0[f]) / (W @ n0); m1 = (W @ S1[f]) / (W @ n1); b = 100 * (m1 - m0) / m0
            lo, hi = np.percentile(b, [2.5, 97.5]); CI[(f, rlab)] = (obs, lo, hi)
            BS.append(dict(gene=f, remission_group=rlab, net_pct_change=obs, ci95_low=lo, ci95_high=hi, bootstrap_se=b.std(ddof=1),
                           change_p=PV[(f, rlab)], difference_p=PV[(f, 'difference')],
                           n_patients=len(pats), n_patients_pre=int((n0 > 0).sum()), n_patients_post=int((n1 > 0).sum()),
                           n_biopsies_pre=int(n0.sum()), n_biopsies_post=int(n1.sum()), all_cell_mean_pre=S0[f].sum() / n0.sum(),
                           n_boot=N_BOOT, seed=SEED, n_redrawn_no_pre_or_post=redraw,
                           method='patient-level cluster bootstrap within remission group; percentile 95% CI; all biopsies, not paired-only',
                           p_source='All cells, biopsy-level mixed model: biopsy mean ~ timepoint x remission + (1|patient) + (1|patient:timepoint)'))
    BS = pd.DataFrame(BS).rename(columns={'ci95_low': 'bootstrap_ratio_ci95_low', 'ci95_high': 'bootstrap_ratio_ci95_high',
                                           'bootstrap_se': 'bootstrap_ratio_se', 'n_boot': 'bootstrap_ratio_n_boot', 'seed': 'bootstrap_ratio_seed',
                                           'n_redrawn_no_pre_or_post': 'bootstrap_ratio_n_redrawn_no_pre_or_post', 'method': 'bootstrap_ratio_method'})
    BCI = dict(CI)                                                                 # bootstrap (kept for the CSV only)
    # ---------------------------------------------------------------- model-based SE and 95% CI
    CI = {}; SEP = {}; MOD = []
    for f in GENES:
        for rl in ['remission', 'non-remission']:
            r_ = A.loc[(f, CHG[rl])]; n_ = NP.loc[(f, CHG[rl])]
            for c in ['estimate', 'ci_low', 'ci_high']: assert np.isclose(r_[c], n_[c], rtol=1e-6), (f, rl, c)   # both CSVs agree
            den = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].all_cell_mean_pre.iloc[0]
            obs = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].pct_net.iloc[0]
            e, lo, hi = (100 * r_[c] / den for c in ['estimate', 'ci_low', 'ci_high'])
            assert lo < e < hi
            excl = (lo > 0) or (hi < 0); assert excl == (PV[(f, rl)] < 0.05), (f, rl, lo, hi, PV[(f, rl)])        # CI excludes 0 <=> P < 0.05
            gap = e - obs; assert abs(gap) <= 3, ('model vs observed differ by > 3 percentage points', f, rl, e, obs)
            CI[(f, rl)] = (e, lo, hi); sep = 100 * r_.se / den; SEP[(f, rl)] = sep
            assert np.isclose((hi - lo) / 2 / sep, 1.986, atol=0.02), (f, rl, (hi - lo) / 2 / sep)        # CI half-width = t(0.975, df 92) x SE
            MOD.append(dict(gene=f, remission_group=rl, model_pct_change=e, model_ci95_low=lo, model_ci95_high=hi,
                            model_estimate_lognorm=r_.estimate, model_ci95_low_lognorm=r_.ci_low, model_ci95_high_lognorm=r_.ci_high,
                            model_se_lognorm=r_.se, model_se_pct=sep, model_minus_observed_pct_points=gap, ci_excludes_0=excl))
            print(f'MODEL CI {f:5s} {rl:14s} est {e:7.2f}%  CI {lo:7.2f} to {hi:7.2f}%  observed {obs:7.2f}%  gap {gap:+.2f} pp  P {p2(PV[(f, rl)])}  excl0 {excl}')
    BS = pd.DataFrame(MOD).merge(BS, on=['gene', 'remission_group'], validate='1:1')
    BS['model_ci_source'] = ('All cells change (post - pre) from the biopsy-level mixed model (cn10_cytokine_coarse_screen.csv), log-normalised '
                         'units / observed pre-treatment all-cell mean of the group (constant denominator) x 100')
    BS['model_minus_se_pct'] = BS.model_pct_change - BS.model_se_pct; BS['model_plus_se_pct'] = BS.model_pct_change + BS.model_se_pct
    BS = BS[['gene', 'remission_group', 'net_pct_change', 'model_pct_change', 'model_se_pct', 'model_minus_se_pct', 'model_plus_se_pct', 'model_ci95_low', 'model_ci95_high', 'model_minus_observed_pct_points',
             'ci_excludes_0', 'change_p', 'difference_p', 'all_cell_mean_pre', 'model_estimate_lognorm', 'model_ci95_low_lognorm',
             'model_ci95_high_lognorm', 'model_se_lognorm', 'n_patients', 'n_patients_pre', 'n_patients_post', 'n_biopsies_pre', 'n_biopsies_post',
             'model_ci_source', 'p_source'] + [c for c in BS.columns if c.startswith('bootstrap_ratio')]]
    out = BS.copy()
    for c in out.columns:
        if out[c].dtype.kind == 'f': out[c] = out[c].map(sig3)
    with open(BOOT, 'w') as fh: out.to_csv(fh, index=False)
    print(out[['gene', 'remission_group', 'net_pct_change', 'model_pct_change', 'model_se_pct', 'model_minus_se_pct', 'model_plus_se_pct', 'model_ci95_low', 'model_ci95_high', 'change_p',
               'bootstrap_ratio_ci95_low', 'bootstrap_ratio_ci95_high']].to_string())


# ======================================================================================================================
# Section F6B_cytokine_sources: cell-group sources of the change in OSM, IL1B and IL6 after anti-TNF (figure and legend; plot only)
# Panels: Fig 6B
# ======================================================================================================================
if _run('F6B_cytokine_sources'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/57f_fig6_cn10_cytokine_combined.py'   # working folder of this step (relative paths resolve from here)
    """Fig 6B, TAURUS anti-TNF, CD scRNA-seq. Plot only: no model is refitted.
One row, left to right: 3 bar panels | abundance dots | strip | keys.
Bars: stacked cell-group contributions to the % change in the all-cell mean of OSM, IL1B and IL6 (Remission, Non-remission), net change
  (diamond) with +/- 1 SE of the biopsy-level mixed model, 'Difference P' (All cells, differential change) above each gene.
Dots: abundance (% of all cells) 'Difference in change' and 'Baseline' contrasts of the coarse groups, plus the strip (% of pre-treatment
  transcripts of OSM, IL1B, IL6 and % of cells). Rows selected by a stated rule applied to the strip values (cell cache, = CSV):
  Main sources = groups holding >= 15% of the pre-treatment transcripts of OSM, IL1B or IL6, ordered by their largest share |
  thin separator | Contrast = the remaining groups with >= 8% of pre-treatment cells, ordered by % of cells.
  A colour swatch before each row name = that group's colour in the stacked bars (groups pooled in 'Other groups' get the Other-groups
  grey). Dot size = 4 discrete P bins; black ring = P < 0.05.
Inputs: supplementary/cn10_cytokine_coarse_screen.csv, supplementary/nonremission_persistent_signals.csv (same All Cells P values),
  supplementary/cache/cn10_cytokine_cells.parquet, supplementary/cache/cn10_cycling_lineage_check.csv,
  supplementary/cn10_cytokine_contrib_bootstrap.csv (F6B_contribution_se; same model estimate and SE).
Outputs: fig6_cn10_cytokine_combined.png, fig6_cn10_cytokine_combined_legend.txt
Run: OPENBLAS_CORETYPE=SkylakeX /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, re, math, warnings
    warnings.filterwarnings('ignore')
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import numpy as np, pandas as pd

    R = '/share/fsmresfiles/UC/revision_jenny'; TOP = f'{R}/Fig6_treg_antitnf_taurus'; HERE = f'{TOP}/supplementary'
    CSV = f'{HERE}/cn10_cytokine_coarse_screen.csv'; NRP = f'{HERE}/nonremission_persistent_signals.csv'
    CACHE = f'{HERE}/cache/cn10_cytokine_cells.parquet'; CYC = f'{HERE}/cache/cn10_cycling_lineage_check.csv'
    BOOT = f'{HERE}/cn10_cytokine_contrib_bootstrap.csv'
    PNG = f'{TOP}/fig6_cn10_cytokine_combined.png'; LEG = f'{TOP}/fig6_cn10_cytokine_combined_legend.txt'
    MIN_N = 10; MIN_PRE_PAT = 3
    GENES = ['OSM', 'IL1B', 'IL6']; AB = 'Abundance'; FEATS = GENES + [AB]
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
    def p2(p):                                                                    # 2 significant figures, trailing zeros kept (0.0020)
        d = max(1 - math.floor(math.log10(abs(p))), 0); return f'{p:.{d}f}'

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
    NOGRP = 'Neuron (no coarse group)'; ALL = 'All cells'; GROUPS = list(COARSE)
    NAME = {'Monocyte-derived Mac': 'Mo-Mac'}                                   # display name only (CSV keeps the analysis name)
    disp = lambda s: re.sub(r'\+(?=[\s)A-Z]|$)', '⁺', NAME.get(s, s))
    CHG = {'remission': 'Change in remission (post − pre)', 'non-remission': 'Change in non-remission (post − pre)'}
    DIFC = 'Differential change (remission − non-remission)'
    CON = {'baseline_difference': 'Baseline difference (remission − non-remission, pre)',
           'differential_change': DIFC, 'change_remission': CHG['remission'], 'change_nonremission': CHG['non-remission']}
    MAIN_TX = 15.0    # main source: >= 15% of the pre-treatment transcripts of OSM, IL1B or IL6
    CONTRAST_CELLS = 8.0   # contrast: remaining groups with >= 8% of pre-treatment cells

    # ================================================================ P values (All cells, biopsy-level mixed model; F6B_coarse_screen CSV)
    RES = pd.read_csv(CSV)
    A = RES[(RES.part == '1 expression') & (RES.cell_group == ALL) & RES.feature.isin(GENES)].set_index(['feature', 'contrast'])
    NP = pd.read_csv(NRP); NP = NP[(NP.cell_group == 'All Cells') & NP.gene_or_program.isin(GENES)].set_index(['gene_or_program', 'contrast'])
    PV = {}
    for f in GENES:
        for k, c in list(CHG.items()) + [('difference', DIFC)]:
            p = float(A.loc[(f, c), 'p']); assert np.isclose(p, NP.loc[(f, c), 'p'], rtol=1e-6), (f, c)
            PV[(f, k)] = p

    # ================================================================ dot-plot table
    RR = RES[RES.part.isin(['1 expression', '2 abundance']) & RES.contrast.notna()].copy()
    RR['contrast_key'] = RR.contrast.map({v: k for k, v in CON.items()}); assert RR.contrast_key.notna().all()
    assert len(RR) == 4 * len(FEATS) * (len(GROUPS) + 1) and not RR.duplicated(['cell_group', 'feature', 'contrast_key']).any()
    E = RES[RES.part == '1 expression'].drop_duplicates(['cell_group', 'feature'])[['cell_group', 'feature', 'share_of_transcripts_pre']]
    G = RR.drop_duplicates('cell_group').set_index('cell_group')
    CY = pd.read_csv(CYC).set_index('item').value
    cyc_t = max(CY['% of Cycling cells expressing CD3D'], CY['% of Cycling cells expressing CD3E']); cyc_e = CY['% of Cycling cells expressing EPCAM']

    # ================================================================ decomposition + strip shares from the cell cache
    C = pd.read_parquet(CACHE); C['biopsy'] = C.sample_id
    C[GENES] = C[GENES].astype('float64')
    T2G = {t: g for g, ts in COARSE.items() for t in ts}
    assert set(C.aucell_cell_type_short) - set(T2G) == {'Neuron'}
    C['grp'] = C.aucell_cell_type_short.map(T2G).fillna(NOGRP)
    pre = (C.Treatment == 'Pre').values
    TOT_PRE = {g: C.loc[pre, f'{g}_counts'].sum() for g in GENES}
    SHX = {(g, f): 100 * C.loc[pre & (C.grp == g).values, f'{f}_counts'].sum() / TOT_PRE[f] for g in GROUPS for f in GENES}
    SHX.update({(ALL, f): 100.0 for f in GENES})
    PCX = {g: 100 * (pre & (C.grp == g).values).sum() / pre.sum() for g in GROUPS}; PCX[ALL] = 100.0
    for (g, f), v in SHX.items(): assert np.isclose(v, E.set_index(['cell_group', 'feature']).share_of_transcripts_pre[(g, f)], rtol=5e-3, atol=1e-3), (g, f)
    for g, v in PCX.items(): assert np.isclose(v, G.pct_of_all_cells_pre[g], rtol=5e-3, atol=1e-3), g
    # ---- row selection rule (exact strip values)
    MAXTX = {g: max(SHX[(g, f)] for f in GENES) for g in GROUPS}
    MAIN = sorted([g for g in GROUPS if MAXTX[g] >= MAIN_TX], key=lambda g: -MAXTX[g])
    CONTRAST = sorted([g for g in GROUPS if g not in MAIN and PCX[g] >= CONTRAST_CELLS], key=lambda g: -PCX[g])
    ROWS = MAIN + CONTRAST
    for g in GROUPS:
        print(f'RULE {disp(g):24s} max transcript share {MAXTX[g]:6.2f}%  cells {PCX[g]:6.2f}%  ->',
              'main' if g in MAIN else 'contrast' if g in CONTRAST else 'not shown')
    print('MAIN', [disp(g) for g in MAIN]); print('CONTRAST', [disp(g) for g in CONTRAST])
    BM = C.groupby('biopsy')[['Patient', 'Remission_status', 'Treatment']].first().assign(n_all=C.groupby('biopsy').size())
    cnt = C.groupby(['biopsy', 'grp']).size().unstack(fill_value=0).reindex(BM.index)
    frac = cnt.div(BM.n_all, axis=0)
    DEC = []
    for f in GENES:
        mg = C.groupby(['biopsy', 'grp'])[f].mean().unstack().reindex(BM.index)
        contrib = (frac * mg).fillna(0); allm = C.groupby('biopsy')[f].mean().reindex(BM.index)
        assert np.allclose(contrib.sum(1), allm, atol=1e-10)
        for rs, rlab in [('Remission', 'remission'), ('Non_Remission', 'non-remission')]:
            b0 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Pre')]; b1 = BM.index[(BM.Remission_status == rs) & (BM.Treatment == 'Post')]
            tot0, tot1 = allm[b0].mean(), allm[b1].mean()
            for g in GROUPS + [NOGRP]:
                C0, C1 = contrib.loc[b0, g].mean(), contrib.loc[b1, g].mean()
                DEC.append(dict(feature=f, cell_group=g, remission_group=rlab, contribution_to_change=C1 - C0,
                                all_cell_mean_pre=tot0, all_cell_change=tot1 - tot0))
    DEC = pd.DataFrame(DEC)
    DEC['pct_contribution'] = 100 * DEC.contribution_to_change / DEC.all_cell_mean_pre
    DEC['pct_net'] = 100 * DEC.all_cell_change / DEC.all_cell_mean_pre
    for (f, rl), d in DEC.groupby(['feature', 'remission_group']):
        assert np.isclose(d.pct_contribution.sum(), d.pct_net.iloc[0], rtol=1e-12, atol=1e-10), (f, rl)
    CD = RES[(RES.part == '3 decomposition') & RES.remission_group.isin(['remission', 'non-remission'])]
    chk = DEC.merge(CD[['feature', 'cell_group', 'remission_group', 'contribution_to_change', 'all_cell_mean_pre', 'all_cell_change']],
                    on=['feature', 'cell_group', 'remission_group'], suffixes=('', '_csv'))
    assert len(chk) == len(DEC) == len(CD)
    for c in ['contribution_to_change', 'all_cell_mean_pre', 'all_cell_change']: assert np.allclose(chk[c], chk[f'{c}_csv'], rtol=5e-3, atol=1e-12), c
    DIFF = RES[(RES.part == '3 decomposition') & RES.remission_group.astype(str).str.startswith('differential')]
    SPLIT = {f: 100 * DIFF[DIFF.feature == f].abundance_part.sum() / DIFF[DIFF.feature == f].contribution_to_change.sum() for f in GENES}
    print('decomposition matches CSV; split', {f: round(v, 1) for f, v in SPLIT.items()})

    # ================================================================ model estimate +/- SE in % (from the CSVs, no refit)
    EST = {}; SEP = {}
    BS = pd.read_csv(BOOT).set_index(['gene', 'remission_group'])
    for f in GENES:
        for rl in ['remission', 'non-remission']:
            r_ = A.loc[(f, CHG[rl])]
            den = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].all_cell_mean_pre.iloc[0]
            obs = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].pct_net.iloc[0]
            EST[(f, rl)] = 100 * r_.estimate / den; SEP[(f, rl)] = 100 * r_.se / den
            assert abs(EST[(f, rl)] - obs) <= 3
            for col, v in [('model_pct_change', EST[(f, rl)]), ('model_se_pct', SEP[(f, rl)]), ('net_pct_change', obs)]:
                assert np.isclose(BS.loc[(f, rl), col], sig3(v), rtol=2e-3), (f, rl, col)          # same values as the F6B_contribution_se CSV

    # ================================================================ figure
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.legend_handler import HandlerTuple
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14, 'xtick.labelsize': 13, 'ytick.labelsize': 13,
                         'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0})
    K, GR = '#111111', '#A6A6A6'
    big = DEC.assign(r=DEC.contribution_to_change.abs() / DEC.all_cell_change.abs())
    SRCG = {g for (g, f), v in SHX.items() if g in GROUPS and v >= 10}
    COLG = [g for g in GROUPS if (big[big.cell_group == g].r >= 0.10).any() or g in SRCG]
    PAL = {'Macrophage': '#D55E00', 'Monocyte-derived Mac': '#E69F00', 'DC': '#CC79A7', 'Granulocyte/Neutrophil': '#F0E442', 'B cell': '#7B3294',
           'Plasma': '#C2A5CF', 'CD4 T': '#0072B2', 'CD8 T': '#56B4E9', 'Other T/NK': '#009E73', 'Fibroblast': '#8C6D31',
           'Pericyte': '#C9B37E', 'Endothelial': '#4D4D4D', 'Epithelial': '#80CDC1'}
    OTH = 'Other groups'; PAL[OTH] = '#D9D9D9'
    BARCOL = {g: (PAL[g] if g in COLG else PAL[OTH]) for g in MAIN + CONTRAST}     # colour of the group in the stacked bars
    print('coloured groups', COLG); print('row swatches', {disp(g): c for g, c in BARCOL.items()})
    assert set(MAIN) <= set(COLG)

    # ---------------------------------------------------------------- one row, explicit layout in inches
    # [3 compact bar panels][group label | swatch | row name][abundance dots][strip][colour key][P key]
    PH = 4.2; TOPM = 1.0; BOTM = 0.6; H = TOPM + PH + BOTM                         # panel height shared by bars, dots and strip
    ny = len(ROWS); RU = PH / ny                                                   # rows fill the bar-panel height
    BX0 = 0.95; BW = 1.85; BG = 0.72                                              # bar panels: left, width, gap
    XA = BX0 + 3 * BW + 2 * BG + 2.05; U = 0.95                                     # dot axes start after the row-label margin
    WA = 2.2 * U; XS = XA + WA + 0.2; WS = 4.1 * U; XK = XS + WS + 0.45; XP = XK + 1.65
    W = XP + 1.75
    fig = plt.figure(figsize=(W, H))
    yb_ = BOTM / H; hA = PH / H; ytop = (BOTM + PH) / H
    in2y = lambda v: v / H                                                         # inches -> figure fraction (vertical)

    # ---------------------------------------------------------------- bars (compact panels, shared styling)
    CAPF = 0.6
    BARV = []
    for i, f in enumerate(GENES):
        a = fig.add_axes([(BX0 + i * (BW + BG)) / W, yb_, BW / W, hA]); ext = []                 # own y-range per gene (not shared); tick labels on every panel; y title on the first only
        stk = {}
        for rl in ['remission', 'non-remission']:
            s_ = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].pct_contribution
            stk[rl] = (s_[s_ > 0].sum(), s_[s_ < 0].sum())
        bmax = max([0] + [u for u, _ in stk.values()]); bmin = min([0] + [d_ for _, d_ in stk.values()]); span = bmax - bmin
        capH, capL = bmax + CAPF * span, bmin - CAPF * span
        for xi, rl in enumerate(['remission', 'non-remission']):
            d = DEC[(DEC.feature == f) & (DEC.remission_group == rl)].set_index('cell_group')
            s = d.pct_contribution; net_true = d.pct_net.iloc[0]
            vals = [(g, s[g]) for g in COLG] + [(OTH, s.drop(COLG).sum())]
            up = dn = 0.0
            for g, v in vals:
                if v >= 0: a.bar(xi, v, bottom=up, width=0.62, color=PAL[g], edgecolor='white', lw=0.6, zorder=2); up += v
                else: a.bar(xi, v, bottom=dn, width=0.62, color=PAL[g], edgecolor='white', lw=0.6, zorder=2); dn += v
            net = sum(v for _, v in vals)
            assert abs((up + dn) - net) < 1e-9 and abs(net - net_true) < 1e-9, (f, rl, up + dn, net, net_true)
            e_ = EST[(f, rl)]; lo, hi = e_ - SEP[(f, rl)], e_ + SEP[(f, rl)]
            lo_d, hi_d = max(lo, capL), min(hi, capH)
            a.plot([xi, xi], [lo_d, hi_d], color=K, lw=1.2, solid_capstyle='butt', zorder=3.5)
            for v_, vd_, mk in [(hi, hi_d, '^'), (lo, lo_d, 'v')]:
                if v_ != vd_:
                    a.scatter(xi, vd_, marker=mk, s=45, color=K, lw=0, zorder=3.6)
                    a.annotate(f'{v_:.0f}%', (xi, vd_), xytext=(7, -6 if mk == '^' else 6), textcoords='offset points', ha='left', va='top' if mk == '^' else 'bottom', fontsize=12, zorder=5)
            a.scatter(xi, net, marker='D', s=70, color='white', edgecolor=K, linewidths=1.8, zorder=4)
            top_, bot_ = max(up, hi_d), min(dn, lo_d); ext += [top_, bot_]
            BARV.append((f, rl, {disp(g): sig3(v) for g, v in vals}, 'net %', sig3(net), '±SE', (sig3(lo), sig3(hi)), 'P', p2(PV[(f, rl)])))
        lo_, hi_ = min(ext + [0]), max(ext + [0]); rg = hi_ - lo_
        a.set_ylim(lo_ - 0.08 * rg, hi_ + 0.08 * rg)
        a.axhline(0, color=K, lw=1.0, zorder=3)
        a.set_xticks([0, 1]); a.set_xticklabels(['Remission', 'Non-\nremission'], fontsize=13); a.set_xlim(-0.55, 1.55)
        a.spines[['top', 'right']].set_visible(False); a.tick_params(axis='y', labelsize=12)
        a.yaxis.set_major_locator(plt.MaxNLocator(5)); a.ticklabel_format(axis='y', style='plain')
        xc = (BX0 + i * (BW + BG) + BW / 2) / W
        fig.text(xc, ytop + in2y(0.44), f, ha='center', va='bottom', fontsize=15, fontstyle='italic')
        fig.text(xc, ytop + in2y(0.08), f'Difference P = {p2(PV[(f, "difference")])}', ha='center', va='bottom', fontsize=12.5)
        if i == 0: a.set_ylabel('% change in all-cell mean\n(relative to pre-treatment)', fontsize=13.5)

    # ---------------------------------------------------------------- dot plot (abundance) + strip
    VM = 1.5; norm = TwoSlopeNorm(vcenter=0, vmin=-VM, vmax=VM); cmap = plt.get_cmap('RdBu_r')
    # 4 DISCRETE size bins matching the key exactly; every ringed dot (P < 0.05) is in one of the two largest bins
    P_BINS = [0.01, 0.05, 0.1]; P_SIZES = [520, 330, 140, 60]                     # P < 0.01 | 0.01-0.05 | 0.05-0.1 | P >= 0.1
    P_LABELS = ['P ≥ 0.1', '0.05–0.1', '0.01–0.05', '< 0.01']
    pbin = lambda p: 0 if p < 0.01 else 1 if p < 0.05 else 2 if p < 0.1 else 3
    psize = lambda p: P_SIZES[pbin(p)]
    COLSP = [(AB, k) for k in ['differential_change', 'baseline_difference']]
    ypos = {r: i for i, r in enumerate(ROWS)}
    ax = fig.add_axes([XA / W, yb_, WA / W, hA]); axs = fig.add_axes([XS / W, yb_, WS / W, hA], sharey=ax)
    EL = RR[RR.contrast_key.isin(['baseline_difference', 'differential_change'])].set_index(['cell_group', 'feature', 'contrast_key'])
    ELp = EL.loc[[(r, f, k) for r in ROWS for f, k in COLSP]]
    assert not (ELp.q_bh < 0.05).any(), 'a plotted dot reaches BH q < 0.05: the key needs the asterisk back'
    assert not ELp.p.between(0.0495, 0.0505).any(), 'a P value sits on the 0.05 ring threshold at 3 s.f.'
    PLOTTED = []
    for j, (f, k) in enumerate(COLSP):
        for r, y in ypos.items():
            x = EL.loc[(r, f, k)]; assert np.isfinite(x.p)
            ax.scatter(j, y, s=psize(x.p), c=[cmap(norm(np.clip(x.estimate_sd, -VM, VM)))], edgecolors=K if x.p < 0.05 else '#9A9A9A',
                       linewidths=2.0 if x.p < 0.05 else 0.6, zorder=3)
            assert (pbin(x.p) <= 1) == (x.p < 0.05)                                 # ring <=> one of the two largest bins
            PLOTTED.append((disp(r), f, k, round(x.estimate_sd, 2), sig3(x.p), round(float(psize(x.p))), bool(x.p < 0.05)))
    ax.set_xlim(-0.6, len(COLSP) - 0.4); ax.set_ylim(ny - 0.5, -0.5)
    ax.set_yticks(list(ypos.values())); ax.set_yticklabels([disp(r) for r in ypos], fontsize=14)
    ax.set_xticks(range(len(COLSP))); ax.set_xticklabels(['Difference\nin change', 'Base-\nline'], fontsize=12.5, linespacing=1.0)
    ax.xaxis.tick_top(); ax.tick_params(axis='x', length=0, pad=2)
    ax.text(0.5, ytop + in2y(0.44), 'Abundance', transform=fig.transFigure, ha='center', va='bottom', fontsize=15)
    ax.text(0.5, ytop + in2y(0.70), '% of all cells', transform=fig.transFigure, ha='center', va='bottom', fontsize=14, color='#444444')
    ax.texts[-2].set_x((XA + WA / 2) / W); ax.texts[-1].set_x((XA + WA / 2) / W)
    ax.axhline(len(MAIN) - 0.5, color=GR, lw=0.8)                                  # thin separator: main sources | contrast
    for a_ in [ax, axs]:
        a_.spines[['top', 'right', 'bottom']].set_visible(False)
        for r, yy in ypos.items():
            if yy % 2 == 0: a_.axhspan(yy - 0.5, yy + 0.5, color='#F2F2F2', zorder=0, lw=0)
    PCM = plt.get_cmap('Purples'); PVMAX = 40
    STRIP = []
    for j, f in enumerate(GENES + ['cells']):
        for r, y in ypos.items():
            v = SHX[(r, f)] if f != 'cells' else PCX[r]
            axs.add_patch(plt.Rectangle((j - 0.47, y - 0.45), 0.94, 0.9, color=PCM(0.08 + 0.85 * min(v, PVMAX) / PVMAX), lw=0, zorder=1))
            lab = f'{sig3(round(v, 1)):.1f}' if v < 10 else f'{v:.0f}'
            axs.text(j, y, lab, ha='center', va='center', fontsize=12, color='white' if v > 0.6 * PVMAX else K, zorder=2)
            STRIP.append((disp(r), f, lab))
    axs.axvline(2.5, color='white', lw=4, zorder=3)
    axs.set_xlim(-0.55, 3.55); axs.set_xticks([]); axs.tick_params(axis='y', left=False, labelleft=False)
    for j, f in enumerate(GENES + ['Cells']): axs.text(j, -0.62, f, ha='center', va='bottom', fontsize=13, fontstyle='italic' if f in GENES else 'normal', clip_on=False)
    fig.text((XS + WS / 2) / W, ytop + in2y(0.44), '% of pre-treatment\ntranscripts | cells', ha='center', va='bottom', fontsize=14, linespacing=1.05)
    axs.spines['left'].set_visible(False)
    # keys, far right, top-aligned with the panels: colour scale | P size bins + ring note
    fig.text(XK / W, ytop, 'Remission −\nnon-remission\n(SD of biopsy\nvalues)', ha='left', va='top', fontsize=12.5, linespacing=1.05)
    cax = fig.add_axes([(XK + 0.05) / W, ytop - in2y(3.2), 0.2 / W, in2y(2.0)])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[-1.5, 0, 1.5]); cb.ax.set_yticklabels(['−1.5', '0', '1.5'], fontsize=12)
    hd = [plt.scatter([], [], s=sz, c='#BBBBBB', edgecolors='#9A9A9A', linewidths=0.6) for sz in P_SIZES[::-1]]   # key = the 4 plotted sizes
    fig.legend(hd, P_LABELS, title='P', loc='upper left', bbox_to_anchor=(XP / W, ytop + in2y(0.05)), frameon=False,
               fontsize=12, title_fontsize=13, labelspacing=0.9, borderpad=0.2, handletextpad=0.6)
    fig.text((XP + 0.08) / W, ytop - in2y(2.45), '○ black ring:\n   P < 0.05', ha='left', va='top', fontsize=12, linespacing=1.25)

    # ---- colour swatches before each row name + grey italic group labels at the left margin
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    SW = 0.17; SWG = 0.07                                                          # swatch side and gap to the name (inches)
    sx0 = {}
    for tl, r in zip(ax.get_yticklabels(), ROWS):
        if r == ALL: continue
        bb = tl.get_window_extent(rend).transformed(inv)
        x1 = bb.x0 - SWG / W; x0 = x1 - SW / W; yc = (bb.y0 + bb.y1) / 2
        fig.add_artist(Rectangle((x0, yc - SW / H / 2), SW / W, SW / H, transform=fig.transFigure, facecolor=BARCOL[r], edgecolor='white', lw=0.6, clip_on=False))
        sx0[r] = x0
    xmin = min(sx0.values())
    ylab = lambda yy: inv.transform(ax.transData.transform((0, yy)))[1]
    for lab, blk in [('Main sources', MAIN), ('Contrast', CONTRAST)]:
        yt, yb = ylab(ypos[blk[0]] - 0.4), ylab(ypos[blk[-1]] + 0.4)
        xl = xmin - 0.10 / W
        fig.add_artist(Line2D([xl, xl], [yb, yt], transform=fig.transFigure, color=GR, lw=0.8))
        fig.text(xl - 0.06 / W, (yt + yb) / 2, lab, rotation=90, ha='right', va='center', fontsize=11.5, fontstyle='italic', color='#7F7F7F')
    EXTRA = [g for g in COLG if g not in ROWS] + [OTH]                            # bar colours not keyed by a row swatch
    print('mini key', [disp(g) for g in EXTRA])
    hx = [Patch(facecolor=PAL[g], edgecolor='white') for g in EXTRA]
    dia = (Line2D([], [], marker='|', ls='', ms=22, mew=1.2, color=K), Line2D([], [], marker='D', ls='', ms=8, mfc='white', mec=K, mew=1.8))
    fig.legend(hx + [dia], [disp(g) for g in EXTRA] + ['Net change (sum of groups); bar = ± SE (mixed model)'], handler_map={tuple: HandlerTuple(ndivide=1)},
               loc='upper left', bbox_to_anchor=(xmin - 0.02 / W, yb_ - in2y(0.08)), ncol=len(EXTRA) + 1, frameon=False, fontsize=12.5,
               handlelength=1.1, handletextpad=0.5, columnspacing=1.2, borderpad=0.1)
    fig.savefig(PNG, dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    print('wrote', PNG)

    # ================================================================ legend
    allc = G.loc[ALL]; ne = int(G.loc[GROUPS].eligible.astype(str).eq('True').sum())
    mp = '; '.join(f"{disp(g)} = {', '.join(disp(t) for t in v)}" for g, v in COARSE.items())
    n3 = lambda v: f'{v:#.3g}'.rstrip('.').replace('-', '−')
    OBS = {(f, rl): DEC[(DEC.feature == f) & (DEC.remission_group == rl)].pct_net.iloc[0] for f in GENES for rl in ['remission', 'non-remission']}
    fmt = lambda f, rl: f"{n3(OBS[(f, rl)])}% (model {n3(EST[(f, rl)])} ± {n3(SEP[(f, rl)])}%)"
    gapmax = np.ceil(max(abs(EST[k] - OBS[k]) for k in OBS) * 10) / 10
    notshown = [disp(g) for g in GROUPS if g not in ROWS]
    grey = [disp(g) for g in MAIN + CONTRAST if g not in COLG]
    # expression sentence, derived from the data: no main-source expression contrast reaches BH q < 0.05
    XQ = RES[(RES.part == '1 expression') & RES.cell_group.isin(MAIN) & RES.contrast.notna()]
    assert len(XQ) == 4 * len(GENES) * len(MAIN) and XQ.q_bh.min() > 0.05, XQ.q_bh.min()
    print('main-source expression contrasts: min q', XQ.q_bh.min(), '; P < 0.05:',
          XQ[XQ.p < 0.05][['cell_group', 'feature', 'contrast', 'p']].values.tolist())
    EXTRA_ = [disp(g) for g in COLG if g not in ROWS] + ['Other groups']
    swt = (f"The swatch before each row name is the group's colour in the stacked bars and serves as the bar key"
           + (f" ({' and '.join(grey) if len(grey) <= 2 else ', '.join(grey[:-1]) + ' and ' + grey[-1]} {'is' if len(grey) == 1 else 'are'} within Other groups)" if grey else "")
           + f"; {', '.join(EXTRA_[:-1])} and {EXTRA_[-1]}, which appear in the bars but not as rows, are keyed below the dot plot. ")
    txt = (f"Cell-group sources of the change in OSM, IL1B and IL6 after anti-TNF in TAURUS Crohn's disease scRNA-seq (Remission "
       f"{int(allc.n_patients_rem_pre)} patients, {int(allc.n_biopsies_rem_pre)} pre- and {int(allc.n_biopsies_rem_post)} post-treatment biopsies; "
       f"Non-remission {int(allc.n_patients_nonrem_pre)} patients pre-treatment ({int(allc.n_biopsies_nonrem_pre)} biopsies) and "
       f"{int(allc.n_patients_nonrem_post)} post-treatment ({int(allc.n_biopsies_nonrem_post)} biopsies); all biopsies used, not only paired patients). "
       f"Groups were defined from the cell annotation (aucell_cell_type_short): {mp}. Cycling cells were placed in Epithelial because {cyc_e:.0f}% of them "
       f"express EPCAM and {cyc_t:.1f}% express CD3D/CD3E; Neuron cells are in no group. "
       f"Left, cell-group contributions to the change in the all-cell mean of each gene (log-normalised expression; post- minus pre-treatment, mean over "
       f"biopsies), as % of that remission group's pre-treatment all-cell mean. Per biopsy, the all-cell mean equals the sum over cell groups of (group "
       f"fraction × within-group mean); each segment is the change in a group's mean contribution divided by the pre-treatment all-cell mean, so the "
       f"segments add up to the net % change (diamond, Net change (sum of groups)). Stacked bars and diamonds are observed values. Groups holding ≥ 10% of "
       f"a bar's net change or ≥ 10% of a gene's pre-treatment transcripts are coloured; the rest are pooled (Other groups, which include Neuron). "
       f"Error bars (bar = ± SE (mixed model)): model estimate ± 1 standard error of the change after treatment within each remission group from the "
       f"biopsy-level mixed model that gives the P values, converted to % by dividing by that group's observed pre-treatment all-cell mean (the diamond's "
       f"denominator); each bar is centred on the model estimate, which lies within {gapmax:g} percentage points of the observed net change. Net % change, "
       f"observed (model estimate ± SE): OSM, Remission {fmt('OSM', 'remission')}, Non-remission {fmt('OSM', 'non-remission')}; IL1B, Remission "
       f"{fmt('IL1B', 'remission')}, Non-remission {fmt('IL1B', 'non-remission')}; IL6, Remission {fmt('IL6', 'remission')}, Non-remission "
       f"{fmt('IL6', 'non-remission')}. Difference P above each gene, change in Remission minus change in Non-remission (All cells); the P values for the "
       f"change after treatment within each remission group are in the source data. Y-axes differ between genes. "
       f"Right, cell-group abundance before and after treatment. Rows: main sources, groups holding ≥ {MAIN_TX:g}% of the pre-treatment "
       f"transcripts of OSM, IL1B or IL6 (ordered by largest share); contrast, the remaining groups with ≥ {CONTRAST_CELLS:g}% of pre-treatment cells "
       f"(ordered by % of cells). {len(ROWS)} of the {ne} tested groups are shown ({', '.join(notshown)} are in the source data). " + swt + f"Abundance, per "
       f"biopsy, the group's percentage of all cells (all biopsies). A group was tested if "
       f"≥ {MIN_PRE_PAT} Remission and ≥ {MIN_PRE_PAT} Non-remission patients had a pre-treatment biopsy with ≥ {MIN_N} cells of the group (all {ne} groups "
       f"qualified). Difference in change, change after treatment in Remission minus that in Non-remission; Baseline, Remission − Non-remission before "
       f"treatment. Dot colour, estimate in SD units of the biopsy values (red, higher in Remission; blue, higher in Non-remission; clipped at ±{VM}); dot "
       f"size, P in four bins (P ≥ 0.1, 0.05–0.1, 0.01–0.05 and < 0.01; ringed dots are always in the two largest bins); black ring, P < 0.05; no group reached Benjamini-Hochberg q < 0.05 across "
       f"the {len(GROUPS)} groups per column. Side strip, % of each gene's pre-treatment transcripts "
       f"(raw UMI counts) held by the group, and the group's % of pre-treatment cells. Within-group expression per cell did not differ between remission "
       f"groups for any main source after Benjamini-Hochberg correction (all q > 0.05; source data). "
       f"All P values (bars and dots) are from a mixed model of the per-biopsy value, value ~ timepoint × remission + (1 | patient) + (1 | patient:timepoint) "
       f"(lme4/lmerTest REML, Satterthwaite df; no covariates). Values, SE, model 95% CIs, P and q values, including the groups not shown, are in the source data.")
    assert not re.search(r'Nominal|Monocyte-derived|Diff\.|\bD6\b|Supplementary (Fig|figure)|supplementary table', txt, re.I)
    with open(LEG, 'w') as fh: fh.write(txt + '\n')
    print('\n' + txt)
