#!/usr/bin/env python
"""TAURUS CD scRNA-seq (anti-TNF): Treg Suppression Score and Treg activation programs before and after treatment.

Panels: Fig 6A; Supp Table 11 (biopsy-level models).
Inputs: TAURUS raw_norm_counts_anno.h5ad.gz (01, 02); Treg Suppression Score gene sets (cosmx_cd/10); Treg TF and pathway tables
  of the CD atlas (cd_scrna/12, cd_scrna/14); CellChatDB, CollecTRI and MSigDB snapshots.
Outputs: panels/fig6c_taurus/ (Treg pseudotime, cell-level LMM); Fig6_treg_antitnf_taurus/supplementary/ (biopsy-level model CSVs,
  per-cell cache/); Fig6_treg_antitnf_taurus/fig6_treg_state_lines.png.
Sections (run order): treg_pseudotime, cell_level_lmm (--cache <supplementary>/cache; copy its cell_level_lmm_results.csv to
  supplementary/), F6A_programs_biopsy, F6A_suppression_score_biopsy, F6A_post_treatment_difference (run once after the two
  previous sections; appends rows), F6A_pathways_biopsy, F6A_tnfrsf4_biopsy, F6A_treg_state_lines.
Run: RUN_SECTIONS=<id>[,<id>...] python <file>. Model sections call Rscript (lme4, lmerTest).
"""
import os, sys, glob, time, shutil, tempfile, subprocess, warnings
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
_RUN_SECTIONS = [_s.strip() for _s in os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# ---------------------------------------------------------------- shared paths, constants and helpers
RJ = '/share/fsmresfiles/UC/revision_jenny'; D = f'{RJ}/cd_treg_rebuild_rowfix'; T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
H5AD = f'{T}/raw_norm_counts_anno.h5ad.gz'; P6 = f'{D}/panels/fig6c_taurus'
HERE = f'{RJ}/Fig6_treg_antitnf_taurus/supplementary'; CACHE = f'{HERE}/cache/Treg.parquet'
O54 = f'{HERE}/biopsy_level_prepost.csv'; O55 = f'{HERE}/treg_suppression_biopsy_level.csv'
O65 = f'{HERE}/fig3G_pathways_taurus_prepost.csv'; O66 = f'{HERE}/tnfrsf4_tnfrsf18_prepost.csv'
PCACHE = f'{HERE}/cache/fig3G_pathways_treg.parquet'
MIN_TREG = 10; MOD = 'patient + patient:timepoint'
IL2, TNFR2 = 'HALLMARK_IL2_STAT5_SIGNALING', 'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY'
CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
       'differential_change': 'Differential change (remission − non-remission)',
       'baseline_difference': 'Baseline difference (remission − non-remission, pre)',
       'post_difference': 'Post-treatment difference (remission − non-remission, post)'}
CON4 = ['change_remission', 'change_nonremission', 'differential_change', 'baseline_difference']
GRP = [('rem_pre', 1.0, 0.0), ('rem_post', 1.0, 1.0), ('nonrem_pre', 0.0, 0.0), ('nonrem_post', 0.0, 1.0)]
FORMULA = 'biopsy mean ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df'
sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
_t0 = time.time(); log = lambda *a: print(f'[{time.time() - _t0:5.0f}s]', *a, flush=True)


def write_sig3(R, path):
    with open(path, 'w') as fh:
        R.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', path, R.shape)


def load_tregs():
    import scanpy as sc
    ad = sc.read_h5ad(H5AD); ct = ad.obs['aucell_cell_type_short'].astype(str).values
    treg = ad[ct == 'Treg'].copy(); del ad; log('Tregs', treg.n_obs)
    return treg


def suppression_score(treg):
    """Treg Suppression Score = score_genes(early) - score_genes(late), biopsy (WTX) gene sets present in TAURUS."""
    import scanpy as sc
    SG = pd.read_csv(glob.glob(f'{D}/panels/fig4_and_related/D_shared_inputs/treg_*_signature_genes.csv')[0]); s_ = SG[SG.cohort == 'biopsy']
    GS = {side: [g for g in s_.loc[s_.set.str.startswith(side), 'gene'] if g in treg.var_names] for side in ['early', 'late']}
    sc.tl.score_genes(treg, GS['late'], score_name='l', random_state=0)
    sc.tl.score_genes(treg, GS['early'], score_name='e', random_state=0)
    return (treg.obs.e - treg.obs.l).values


def gene_values(treg, g):
    x = treg[:, g].X
    return np.asarray(x.todense() if hasattr(x, 'todense') else x).ravel()


def pseudobulk(F, FCOLS):
    """Per-biopsy mean of each per-cell feature over its Tregs; biopsies with < MIN_TREG Tregs dropped."""
    F['biopsy'] = F.sample_id.astype(str)
    F['post'] = (F.Treatment == 'Post').astype(float); F['rem'] = (F.Remission_status == 'Remission').astype(float)
    PB = (F.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                  post=('post', 'first'), rem=('rem', 'first'), n_treg=('post', 'size'))
          .join(F.groupby('biopsy')[FCOLS].mean()).reset_index())
    DROP = PB[PB.n_treg < MIN_TREG]; PB = PB[PB.n_treg >= MIN_TREG].reset_index(drop=True)
    log(f'{len(PB)} biopsies, {int(PB.n_treg.sum())} Tregs, {PB.Patient.nunique()} patients; dropped {len(DROP)}')
    return PB, DROP


def group_n(PB, tregs=False):
    NN = {}
    for nm, r_, p_ in GRP:
        s = PB[(PB.rem == r_) & (PB.post == p_)]
        NN.update({f'n_biopsies_{nm}': len(s), f'n_patients_{nm}': s.Patient.nunique()})
        if tregs: NN[f'n_tregs_{nm}'] = int(s.n_treg.sum())
    return NN


def paired_n(PB):
    pairs = PB.groupby(['Patient', 'rem']).post.nunique().reset_index()
    return {'n_paired_patients_rem': int(((pairs.post == 2) & (pairs.rem == 1)).sum()),
            'n_paired_patients_nonrem': int(((pairs.post == 2) & (pairs.rem == 0)).sum())}


RCODE_BIOPSY = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); d <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); cols <- readLines(a[2]); out <- list()
mods <- strsplit(a[4], ",")[[1]]; cons <- strsplit(a[5], ",")[[1]]
d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
FORMS <- list(patient_only = ys ~ post * rem + (1 | Patient), patient_timepoint = ys ~ post * rem + (1 | Patient) + (1 | pt))
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1),
           baseline_difference = c(rem = 1), post_difference = c(rem = 1, `post:rem` = 1))
for (mk in mods) for (cc in cols) {
  y <- d[[cc]]; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  m <- lmer(FORMS[[mk]], data = d, REML = TRUE)
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m))
  for (k in cons) {
    if (mk == "patient_only" && k == "baseline_difference") next
    L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
    ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
    out[[length(out) + 1]] <- data.frame(model_key = mk, feature_key = cc, contrast_key = k, estimate = ct$Estimate * sy, se = ct$`Std. Error` * sy,
      df = ct$df, statistic = ct$`t value`, p = ct$`Pr(>|t|)`, var_patient = vc$vcov[vc$grp == "Patient"] * sy^2,
      var_patient_timepoint = if (any(vc$grp == "pt")) vc$vcov[vc$grp == "pt"] * sy^2 else NA,
      var_residual = vc$vcov[vc$grp == "Residual"] * sy^2, singular = isSingular(m),
      converged = !any(grepl("converge|gradient|Hessian|eigenvalue", c(m@optinfo$conv$lme4$messages, m@optinfo$warnings))))
  }
}
write.csv(do.call(rbind, out), a[3], row.names = FALSE)
"""


def lmer_biopsy(PB, FCOLS, mods=('patient_timepoint',), cons=tuple(CON4)):
    """Biopsy-level LMM per feature (y standardised, estimates back-transformed); adds contrast label and 95% CI (Satterthwaite t)."""
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        PB.assign(pt=PB.Patient + '_' + PB.Treatment)[['Patient', 'pt', 'post', 'rem'] + FCOLS].to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/cols.txt', 'w') as fh: fh.write('\n'.join(FCOLS) + '\n')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE_BIOPSY)
        log(f'fitting {len(FCOLS)} features x {len(mods)} models in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/cols.txt', f'{tmp}/r_out.csv', ','.join(mods), ','.join(cons)], check=True)
        R = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    R['contrast'] = R.contrast_key.map(CON)
    tq = stats.t.ppf(0.975, R.df); R['ci_low'] = R.estimate - tq * R.se; R['ci_high'] = R.estimate + tq * R.se
    return R


def order_rows(R, FCOLS, cols, by=('_c', '_f')):
    R['_f'] = R.feature_key.map({k: i for i, k in enumerate(FCOLS)}); R['_c'] = R.contrast.map({v: i for i, v in enumerate(CON.values())})
    return R.sort_values(list(by))[cols].reset_index(drop=True)


# ======================================================================================================================
# treg_pseudotime: diffusion pseudotime of TAURUS Tregs (per cell and per patient x timepoint)
# ======================================================================================================================
if _run('treg_pseudotime'):
    import scanpy as sc, scanpy.external as sce
    cd = sc.read_h5ad(H5AD)
    treg = cd[cd.obs['aucell_cell_type_short'] == 'Treg'].copy(); del cd
    sc.pp.highly_variable_genes(treg, n_top_genes=3000, flavor='seurat_v3', layer='counts')
    sc.pp.scale(treg, max_value=10)
    sc.tl.pca(treg, n_comps=50, use_highly_variable=True)
    sce.pp.harmony_integrate(treg, key='Patient')
    sc.pp.neighbors(treg, use_rep='X_pca_harmony', n_neighbors=15)
    sc.tl.leiden(treg, resolution=0.4, key_added='leiden_treg_full', random_state=42)
    sc.tl.score_genes(treg, ['FOXP3', 'CTLA4', 'BATF', 'IL2RA', 'CD27'], score_name='suppressive_score')
    root_cluster = treg.obs.groupby('leiden_treg_full')['suppressive_score'].mean().idxmax()
    root_cell = treg.obs[treg.obs['leiden_treg_full'] == root_cluster].index[0]
    treg.uns['iroot'] = treg.obs_names.get_loc(root_cell)
    sc.tl.diffmap(treg); sc.tl.dpt(treg)
    dmin, dmax = treg.obs['dpt_pseudotime'].min(), treg.obs['dpt_pseudotime'].max()
    treg.obs['dpt_pseudotime_norm'] = (treg.obs['dpt_pseudotime'] - dmin) / (dmax - dmin)
    treg.obs[['Patient', 'Remission_status', 'Treatment', 'leiden_treg_full', 'dpt_pseudotime', 'dpt_pseudotime_norm']].to_csv(f'{P6}/fig6B_treg_pseudotime_per_cell.csv')


# ======================================================================================================================
# cell_level_lmm: per-cell Treg / inflammatory-axis features and cell-level LMMs (patient + sample random intercepts);
# writes the per-cell cache read by the biopsy-level sections (--cache <dir>)
# ======================================================================================================================
if _run('cell_level_lmm'):
    import statsmodels.api as sm, statsmodels.formula.api as smf
    OUT = f'{P6}/cell_level_lmm'; os.makedirs(OUT, exist_ok=True)
    KEYS = ['Patient', 'Remission_status', 'Treatment']; OGN = 'OGN+RSPO3+ Fib'
    CDIR = sys.argv[sys.argv.index('--cache') + 1] if '--cache' in sys.argv else None

    # features: (family, figure, column, label, compartment table, type)
    FAM_6 = 'Treg programs (5)'; FAM_TF = 'Treg TFs (28)'; FAM_INF = 'Inflammatory axis (18)'
    SET_B = ['STAT5A', 'STAT5B', 'FOXP3', 'RUNX1', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL', 'PRDM1', 'IKZF2', 'BATF']
    E8 = pd.read_csv(f'{D}/panels/extended_data_fig8/ext8_tf_activity_along_pseudotime.csv')
    TFS = list(dict.fromkeys(SET_B + list(E8.tf)))
    MYE = ['Mo-Mac', 'Inf Mo-Mac', 'Trans cDC2/Mac', 'Mac S+XS-', 'Mac S+M+S+', 'Mac S+SG+', 'Mac S+M+P+', 'CD9+ Mac', 'Neutrophil', 'cDC1',
           'cDC2', 'CCR7+ DC']
    FIB = ['VSTM2A+ Crypt Top Fib', 'ADAMDEC1+ Activ Fib', 'ADAMDEC1+ Fib', 'OGN+RSPO3+ Fib', 'T-interact Fib', 'Activ Fib', 'HHIP+ Myofib',
           'PLA2G2A+ ECM Fib', 'GREM2+ Myofib']                                                                           # no SELENOP+ Fib
    A_FEAT = [('OSM', 'Myeloid'), ('IL6', 'Stromal (Fib + Endo)'), ('IL6', 'Myeloid'), ('IL1A', 'Myeloid'), ('IL1B', 'Myeloid'),
              ('IL11', 'Fibroblasts'), ('TNF', 'Myeloid'), ('TNF', 'T Cells')]
    B_FEAT = [(g, c) for c in ['Fibroblasts', 'OGN⁺RSPO3⁺ Fib'] for g in ['OSMR', 'IL11RA', 'IL6ST', 'IL1R1', 'TNFRSF1A']]
    fid = lambda g, c: f'{g} | {c}'
    SAFE = lambda s: ''.join(ch if ch.isalnum() else '_' for ch in s)
    FEATURES = [(FAM_6, 'Treg program', 'pseudotime', 'Treg pseudotime', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'il2_stat5', 'IL-2/STAT5 signaling (Hallmark, ULM)', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'costim_module', 'Co-stimulation / suppressive Treg receptor module', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'tnfr2_pathway', 'TNFR2 non-canonical NF-κB (Reactome, ULM)', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'gene_TNFRSF9', 'TNFRSF9 (4-1BB)', 'Treg', 'gene')]
    FEATURES += [(FAM_TF, 'Treg TF', f'tf_{tf}', f'{tf} activity (CollecTRI, ULM)', 'Treg', 'score') for tf in TFS]
    FEATURES += [(FAM_INF, 'inflammatory ligand', fid(g, c), fid(g, c), c, 'gene') for g, c in A_FEAT]
    FEATURES += [(FAM_INF, 'inflammatory receptor', fid(g, c), fid(g, c), c, 'gene') for g, c in B_FEAT]

    def build_tables():
        import scanpy as sc, decoupler as dc
        ad = sc.read_h5ad(H5AD); log('loaded', ad.shape)
        ct = ad.obs['aucell_cell_type_short'].astype(str); cat = ad.obs['aucell_cell_category_from_type'].astype(str)
        meta = ad.obs[KEYS + ['sample_id', 'Site']].astype(str).assign(total_counts=ad.obs['total_counts'].astype(float).values)
        TB = {}
        # Tregs: IL-2/STAT5 and TNFR2 (no PSM genes) ULM, TNFRSF9, receptor module, 28 TFs, pseudotime
        treg = ad[(ct == 'Treg').values].copy(); treg.layers['lognorm'] = treg.X.copy(); F = meta.loc[treg.obs_names].copy()
        msig = pd.read_csv(f'{D}/archive/csv/msigdb_human_snapshot.csv.gz', usecols=['geneset', 'genesymbol'])
        SETS = [IL2, 'REACTOME_RUNX1_AND_FOXP3_CONTROL_THE_DEVELOPMENT_OF_REGULATORY_T_LYMPHOCYTES_TREGS',
                'GOBP_T_CELL_TOLERANCE_INDUCTION', 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION']
        net = msig[msig.geneset.isin(SETS)].drop_duplicates().rename(columns={'geneset': 'source', 'genesymbol': 'target'})
        dc.mt.ulm(data=treg, net=net, layer='lognorm'); F['il2_stat5'] = np.asarray(dc.pp.get_obsm(treg, key='score_ulm')[:, IL2].X).ravel()
        tn = msig[msig.geneset == TNFR2].drop_duplicates().rename(columns={'geneset': 'source', 'genesymbol': 'target'})
        tn = tn[~tn.target.str.startswith('PSM')]
        dc.mt.ulm(data=treg, net=tn, layer='lognorm'); F['tnfr2_pathway'] = np.asarray(dc.pp.get_obsm(treg, key='score_ulm').X).ravel()
        x = treg[:, 'TNFRSF9'].layers['lognorm']; F['gene_TNFRSF9'] = np.asarray(x.todense() if hasattr(x, 'todense') else x).ravel()
        cc = pd.read_csv(f'{D}/objects/cellchatdb_human_interactions_snapshot.csv')
        gs = lambda s: [g.strip() for g in str(s).replace(';', ',').split(',') if g.strip() and g.strip() != 'nan']
        det = pd.Series(np.asarray((treg.layers['lognorm'] > 0).mean(0)).ravel(), index=treg.var_names)
        pws = ['CD80', 'CD86', 'PD-L1', 'PDL2', 'ICOS', 'OX40', 'GITRL', 'CD70', 'CD40', 'IL2', 'IL10']
        recs = sorted({g for r in cc.loc[cc.pathway_name.isin(pws), 'receptor.symbol'] for g in gs(r)} & set(treg.var_names))
        recs = [g for g in recs if det[g] >= 0.05]
        sc.tl.score_genes(treg, recs, score_name='costim_module'); F['costim_module'] = treg.obs['costim_module'].values
        F.attrs['costim_receptors'] = ';'.join(recs)
        net_all = pd.read_csv(f'{D}/objects/collectri_human_snapshot.csv.gz'); netf = net_all[net_all.source.isin(TFS)]
        dc.mt.ulm(data=treg, net=netf, layer='lognorm', tmin=5); A = dc.pp.get_obsm(treg, key='score_ulm')
        for tf in TFS: F[f'tf_{tf}'] = np.asarray(A[:, tf].X).ravel()
        pt = pd.read_csv(f'{P6}/fig6B_treg_pseudotime_per_cell.csv', index_col=0)
        F['pseudotime'] = pt.loc[F.index, 'dpt_pseudotime_norm'].values
        TB['Treg'] = F
        # inflammatory-axis compartments
        TCELL = sorted(set(ct[cat == 'T']))
        COMP = {'Myeloid': MYE, 'T Cells': TCELL, 'Stromal (Fib + Endo)': FIB + ['Endothelial'], 'Fibroblasts': FIB, 'OGN⁺RSPO3⁺ Fib': [OGN]}
        for c, types in COMP.items():
            genes = sorted({g for g, cc_ in A_FEAT + B_FEAT if cc_ == c}); m = ct.isin(types).values
            X = ad[m][:, genes].X; X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
            tab = meta.loc[ad.obs_names[m]].copy()
            for i, g in enumerate(genes): tab[fid(g, c)] = X[:, i]
            TB[c] = tab
        return TB

    if CDIR and os.path.exists(f'{CDIR}/tables.csv'):
        names = pd.read_csv(f'{CDIR}/tables.csv')
        TB = {r['name']: pd.read_parquet(f"{CDIR}/{r['file']}") for _, r in names.iterrows()}
    else:
        TB = build_tables()
        if CDIR:
            os.makedirs(CDIR, exist_ok=True)
            for k, v in TB.items():
                v = v.copy(); v.attrs = {}; v.to_parquet(f'{CDIR}/{SAFE(k)}.parquet')
            pd.DataFrame({'name': list(TB), 'file': [f'{SAFE(k)}.parquet' for k in TB]}).to_csv(f'{CDIR}/tables.csv', index=False)
            open(f'{CDIR}/costim_receptors.txt', 'w').write(TB['Treg'].attrs['costim_receptors'])

    # models: A nested (primary), B patient only, D + biopsy, E site-matched, F no libsize, G paired-t diagnostic
    CONTRASTS = {CON['change_remission']: {'post': 1.0, 'post:rem': 1.0}, CON['change_nonremission']: {'post': 1.0},
                 CON['differential_change']: {'post:rem': 1.0}}
    CKEY = {k: s for k, s in zip(CONTRASTS, ['change_remission', 'change_nonremission', 'differential_change'])}
    M_NEST = 'A nested LMM: patient + sample random intercepts (primary)'
    M_PAT = 'B patient random intercept only (sensitivity a; liberal)'
    M_BIO = 'D three-level LMM: + biopsy random intercept (sensitivity c)'
    M_SITE = 'E three-level LMM, site-matched biopsies only (sensitivity d)'
    M_NOLIB = 'F nested LMM without libsize_z (sensitivity e)'
    M_DIAG = 'G paired t on per-patient change in sample means (diagnostic)'
    METHODS = [['lbfgs'], ['powell', 'lbfgs'], ['nm', 'lbfgs'], ['bfgs'], ['cg']]
    VC = {'nested': {'sample': '0 + C(sample)'}, 'patient': None, 'biopsy': {'sample': '0 + C(sample)', 'biopsy': '0 + C(biopsy)'}}

    def prep(tab, col, site_matched=False):
        d = tab[['Patient', 'Remission_status', 'Treatment', 'sample_id', 'Site', 'total_counts', col]].rename(columns={col: 'y', 'sample_id': 'biopsy'}).copy()
        if site_matched:   # sites sampled both before and after treatment in that patient
            both = d.groupby(['Patient', 'Site']).Treatment.nunique(); both = both[both == 2].index
            d = d[pd.MultiIndex.from_frame(d[['Patient', 'Site']]).isin(both)]
        d['post'] = (d.Treatment == 'Post').astype(float); d['rem'] = (d.Remission_status == 'Remission').astype(float)
        d['sample'] = d.Patient + '_' + d.Treatment
        ll = np.log(d.total_counts); d['lib_z'] = (ll - ll.mean()) / ll.std()
        mu, sy = d.y.mean(), d.y.std(); d['ys'] = (d.y - mu) / sy
        return d.reset_index(drop=True), sy

    def fit_lmm(d, sy, level, lib=True):
        """REML fit from each optimiser sequence; keep the converged fit with the highest REML log-likelihood."""
        md = smf.mixedlm('ys ~ post * rem' + (' + lib_z' if lib else ''), d, groups='Patient', re_formula='1', vc_formula=VC[level])
        fits, notes, wmsgs = [], [], set()
        for meth in METHODS:
            with warnings.catch_warnings(record=True) as W:
                warnings.simplefilter('always')
                try:
                    r = md.fit(reml=True, method=meth)
                except Exception as e:
                    notes.append(f"{'+'.join(meth)}: {type(e).__name__}"); continue
            fe = r.fe_params; V = r.cov_params().loc[fe.index, fe.index]; dg = np.diag(V.values)
            if r.converged and np.isfinite(r.llf) and np.all(np.isfinite(dg)) and np.all(dg > 0):
                fits.append(('+'.join(meth), r, fe, V))
                wmsgs |= {str(w.message).split('\n')[0][:90] for w in W if not issubclass(w.category, (FutureWarning, DeprecationWarning))}
            else:
                notes.append(f"{'+'.join(meth)}: converged={r.converged}, finite SE={np.all(np.isfinite(dg))}")
        if not fits:
            return dict(ok=False, retries='; '.join(notes))
        llf = np.array([f[1].llf for f in fits]); meth, r, fe, V = fits[int(np.argmax(llf))]
        vcs = dict(zip(md.exog_vc.names, np.asarray(r.vcomp) * sy ** 2)) if VC[level] else {}
        out = dict(ok=True, optimizer=f'{meth} (best REML logLik of {len(fits)} converged; range {llf.max() - llf.min():.3g})',
                   warnings='; '.join(sorted(wmsgs)), retries='; '.join(notes), llf=float(llf.max()),
                   var_patient=float(r.cov_re.iloc[0, 0]) * sy ** 2, var_sample=vcs.get('sample', np.nan), var_biopsy=vcs.get('biopsy', np.nan),
                   var_residual=float(r.scale) * sy ** 2, lib_coef=float(fe['lib_z']) * sy if lib else np.nan, contrasts={})
        for k, w in CONTRASTS.items():
            L = pd.Series(0.0, index=fe.index)
            for n, v in w.items(): L[n] = v
            out['contrasts'][k] = (float(L @ fe) * sy, float(np.sqrt(L.values @ V.values @ L.values)) * sy)
        return out

    def diag_paired_t(d):
        s = d.groupby(['Patient', 'rem', 'post']).y.mean().unstack('post').dropna(); s['chg'] = s[1.0] - s[0.0]; s = s.reset_index()
        r = sm.OLS(s.chg.values, sm.add_constant(s.rem.values)).fit(); df = int(r.df_resid); out = {}
        for k, L in zip(CONTRASTS, [np.array([1.0, 1.0]), np.array([1.0, 0.0]), np.array([0.0, 1.0])]):
            tt = r.t_test(L); out[k] = (float(np.squeeze(tt.effect)), float(np.squeeze(tt.sd)), df, float(np.squeeze(tt.pvalue)))
        return out, len(s)

    def counts(dd, ftype):
        o = {}
        for nm, r_, p_ in GRP:
            x = dd[(dd.rem == r_) & (dd.post == p_)]
            o[f'n_cells_{nm}'] = len(x); o[f'n_samples_{nm}'] = x['sample'].nunique()
            o[f'pct_cells_expressing_{nm}'] = 100 * (x.y > 0).mean() if ftype == 'gene' else np.nan
        return o

    def lmm_rows(fit, base, model, df_bw, flag=''):
        rows = []
        for k in CONTRASTS:
            if not fit['ok']:
                rows.append({**base, 'model': model, 'contrast': k, 'fit_flag': 'FIT FAILED: ' + fit['retries'], 'test': 'LMM REML'}); continue
            est, se = fit['contrasts'][k]; tz = est / se
            p = 2 * (stats.t.sf(abs(tz), df_bw) if np.isfinite(df_bw) else stats.norm.sf(abs(tz)))
            rows.append({**base, 'model': model, 'contrast': k, 'estimate': est, 'se': se, 'df': df_bw, 'statistic': tz, 'p': p,
                         'p_wald_z': 2 * stats.norm.sf(abs(tz)), 'var_patient': fit['var_patient'], 'var_sample': fit['var_sample'],
                         'var_biopsy': fit['var_biopsy'], 'var_residual': fit['var_residual'], 'libsize_coef': fit['lib_coef'], 'reml_loglik': fit['llf'],
                         'optimizer': fit['optimizer'], 'fit_warnings': fit['warnings'], 'fit_flag': flag,
                         'test': 'LMM REML, t with between-within df (sample-level contrast)' if np.isfinite(df_bw) else 'LMM REML, Wald z'})
        return rows

    def run_feature(spec, tab):
        from threadpoolctl import threadpool_limits
        fam, figp, col, lab, comp, ftype = spec
        with threadpool_limits(2), warnings.catch_warnings():
            warnings.simplefilter('ignore')
            d, sy = prep(tab, col); n_pat, n_smp = d.Patient.nunique(), d['sample'].nunique(); df_bw = n_smp - n_pat - 2
            base = dict(family=fam, figure_panel=figp, feature=lab, feature_key=col, compartment=comp, feature_type=ftype,
                        n_patients=n_pat, n_samples=n_smp, n_biopsies=d.biopsy.nunique(), n_cells=len(d), **counts(d, ftype))
            fn = fit_lmm(d, sy, 'nested'); fp = fit_lmm(d, sy, 'patient')
            flag = ''
            if not fn['ok']:
                flag = 'NESTED FIT FAILED: primary row = patient-only random intercept (fallback); ' + fn['retries']
            elif fn['var_sample'] < 1e-6 * fn['var_residual']:
                flag = 'sample variance estimated at ~0 (boundary)'
            rows = lmm_rows(fn if fn['ok'] else fp, base, M_NEST, df_bw if fn['ok'] else np.inf, flag)
            rows += lmm_rows(fp, base, M_PAT, np.inf)
            rows += lmm_rows(fit_lmm(d, sy, 'biopsy'), base, M_BIO, df_bw)
            ds, sys_ = prep(tab, col, site_matched=True); np_, ns_ = ds.Patient.nunique(), ds['sample'].nunique()
            bs = {**base, 'n_patients': np_, 'n_samples': ns_, 'n_biopsies': ds.biopsy.nunique(), 'n_cells': len(ds), **counts(ds, ftype)}
            rows += lmm_rows(fit_lmm(ds, sys_, 'biopsy'), bs, M_SITE, ns_ - np_ - 2)
            rows += lmm_rows(fit_lmm(d, sy, 'nested', lib=False), base, M_NOLIB, df_bw)
            dg, n_pair = diag_paired_t(d)
            for k, (est, se, dfd, p) in dg.items():
                rows.append({**base, 'model': M_DIAG, 'contrast': k, 'estimate': est, 'se': se, 'df': dfd, 'statistic': est / se, 'p': p,
                             'test': f'OLS of per-patient change in sample mean on remission ({n_pair} paired patients), t'})
        return rows

    from joblib import Parallel, delayed
    cols_needed = lambda comp, col: TB[comp][['Patient', 'Remission_status', 'Treatment', 'sample_id', 'Site', 'total_counts', col]]
    OUTS = Parallel(n_jobs=min(len(FEATURES), 26), backend='loky', verbose=0)(delayed(run_feature)(f, cols_needed(f[4], f[2])) for f in FEATURES)
    R = pd.DataFrame([r for rows in OUTS for r in rows])

    # primary model in R lme4 / lmerTest (Satterthwaite df and P)
    RCODE = r'''
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); spec <- read.csv(a[1], stringsAsFactors = FALSE); out <- list()
for (f in unique(spec$file)) {
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE); d$Patient <- factor(d$Patient); d$sample <- factor(d$sample)
  for (i in which(spec$file == f)) {
    y <- d[[spec$column[i]]]; sy <- sd(y); d$ys <- (y - mean(y)) / sy
    m <- lmer(ys ~ post * rem + lib_z + (1 | Patient) + (1 | sample), data = d, REML = TRUE)
    fe <- names(fixef(m))
    Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1))
    for (k in names(Ls)) {
      L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
      ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
      out[[length(out) + 1]] <- data.frame(feature_key = spec$column[i], compartment = spec$compartment[i], contrast_key = k,
        r_estimate = ct$Estimate * sy, r_se = ct$`Std. Error` * sy, df_satterthwaite = ct$df, p_satterthwaite = ct$`Pr(>|t|)`,
        r_singular = isSingular(m))
    }
  }
}
write.csv(do.call(rbind, out), a[2], row.names = FALSE)
'''
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        spec = []
        for comp in sorted({f[4] for f in FEATURES}):
            cols = [f[2] for f in FEATURES if f[4] == comp]; d0, _ = prep(TB[comp], cols[0])
            tab = d0[['Patient', 'sample', 'post', 'rem', 'lib_z']].copy()
            for c in cols: tab[c] = TB[comp][c].values
            fn_ = f'{tmp}/{SAFE(comp)}.csv'; tab.to_csv(fn_, index=False, float_format='%.9g')
            spec += [dict(file=fn_, column=c, compartment=comp) for c in cols]
        pd.DataFrame(spec).to_csv(f'{tmp}/spec.csv', index=False); open(f'{tmp}/x.R', 'w').write(RCODE)
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/spec.csv', f'{tmp}/r_out.csv'], check=True)
        RR = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    RR['contrast'] = RR.contrast_key.map({v: k for k, v in CKEY.items()})
    R = R.merge(RR[['feature_key', 'compartment', 'contrast', 'df_satterthwaite', 'p_satterthwaite', 'r_singular']]
                .assign(model=M_NEST), on=['feature_key', 'compartment', 'contrast', 'model'], how='left')

    # BH within family, per model and contrast
    R['q_bh'] = np.nan; R['q_bh_wald_z'] = np.nan; R['q_satterthwaite'] = np.nan
    for (fam, mod, con), i in R.groupby(['family', 'model', 'contrast']).groups.items():
        i = i[R.loc[i, 'p'].notna().values]
        R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
        if mod == M_NEST:
            R.loc[i, 'q_bh_wald_z'] = multipletests(R.loc[i, 'p_wald_z'], method='fdr_bh')[1]
            R.loc[i, 'q_satterthwaite'] = multipletests(R.loc[i, 'p_satterthwaite'], method='fdr_bh')[1]
    R['estimate_definition'] = np.where(R.model == M_DIAG, 'per-patient change in sample mean (OLS contrast)', 'fixed-effect contrast (units of y)')
    R['bh_family'] = R.family + ' (BH within family, per model and contrast)'
    R['model'] = pd.Categorical(R.model, [M_NEST, M_PAT, M_BIO, M_SITE, M_NOLIB, M_DIAG], ordered=True)
    R['contrast'] = pd.Categorical(R.contrast, list(CONTRASTS), ordered=True)
    feat_order = {(f[2], f[4]): i for i, f in enumerate(FEATURES)}
    R['_f'] = [feat_order[(a, b)] for a, b in zip(R.feature_key, R.compartment)]
    R = R.sort_values(['_f', 'contrast', 'model']).reset_index(drop=True)
    COLS = ['family', 'figure_panel', 'feature', 'compartment', 'feature_type', 'model', 'contrast', 'estimate', 'se', 'df', 'statistic', 'p', 'q_bh',
            'p_wald_z', 'q_bh_wald_z', 'df_satterthwaite', 'p_satterthwaite', 'q_satterthwaite', 'test', 'estimate_definition', 'bh_family',
            'n_patients', 'n_samples', 'n_biopsies', 'n_cells'] + [f'n_cells_{g[0]}' for g in GRP] + [f'n_samples_{g[0]}' for g in GRP] + \
           [f'pct_cells_expressing_{g[0]}' for g in GRP] + \
           ['var_patient', 'var_sample', 'var_biopsy', 'var_residual', 'libsize_coef', 'optimizer', 'fit_warnings', 'fit_flag', 'r_singular', 'feature_key']
    write_sig3(R[COLS], f'{OUT}/cell_level_lmm_results.csv')


# ======================================================================================================================
# F6A_programs_biopsy: biopsy-level LMMs of the Treg programs, pseudotime and 28 Treg TFs (Supp Table 11; Fig 6A receptor module)
# ======================================================================================================================
if _run('F6A_programs_biopsy'):
    FAM_P, FAM_TF = 'Programs (4 + pseudotime)', 'Treg TFs (28)'
    PROG = [('il2_stat5', 'IL-2/STAT5 Signaling (Hallmark, ULM)'), ('costim_module', 'Co-Stimulatory/Suppressive Receptor Module'),
            ('gene_TNFRSF9', 'TNFRSF9 (4-1BB)'), ('tnfr2_pathway', 'TNFR2 Non-Canonical NF-κB (Reactome, ULM)'), ('pseudotime', 'Treg Pseudotime')]
    F = pd.read_parquet(CACHE)
    TFS = [c[3:] for c in F.columns if c.startswith('tf_')]
    FEATS = [(FAM_P, k, l) for k, l in PROG] + [(FAM_TF, f'tf_{t}', f'{t} activity (CollecTRI, ULM)') for t in TFS]
    FCOLS = [k for _, k, _ in FEATS]
    PB, DROP = pseudobulk(F, FCOLS)
    NN = {**group_n(PB, tregs=True), **paired_n(PB), 'n_biopsies_dropped': len(DROP)}

    R = lmer_biopsy(PB, FCOLS, mods=('patient_only', 'patient_timepoint'))
    MODS = {'patient_only': 'patient only', 'patient_timepoint': MOD}
    R['model'] = R.model_key.map(MODS)
    R['model_formula'] = R.model.map({'patient only': 'biopsy mean ~ post × remission + (1 | patient); lme4/lmerTest REML, Satterthwaite df',
                                      MOD: FORMULA})
    R = R.merge(pd.DataFrame(FEATS, columns=['family', 'feature_key', 'feature']), on='feature_key')
    R['q_bh'] = np.nan
    for _, i in R.groupby(['family', 'model', 'contrast']).groups.items():
        R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
    R['bh_family'] = R.family + ', BH within family, per model and contrast'

    # cell-level model F of cell_level_lmm, side by side
    C = pd.read_csv(f'{HERE}/cell_level_lmm_results.csv'); C = C[C.model.str.startswith('F nested LMM without libsize_z')]
    C = C[['feature_key', 'contrast', 'estimate', 'p', 'q_bh']].rename(columns={'estimate': 'cell_F_estimate', 'p': 'cell_F_p', 'q_bh': 'cell_F_q_bh'})
    R = R.merge(C, on=['feature_key', 'contrast'], how='left')
    R['sign_agrees_with_cell_F'] = (np.sign(R.estimate) == np.sign(R.cell_F_estimate)).where(R.cell_F_estimate.notna())
    R['sig_biopsy'] = R.q_bh < 0.05; R['sig_cell_F'] = R.cell_F_q_bh < 0.05
    R['min_tregs_per_biopsy'] = MIN_TREG
    for k, v in NN.items(): R[k] = v
    R['n_patients'] = PB.Patient.nunique(); R['n_biopsies'] = len(PB); R['n_tregs'] = int(PB.n_treg.sum())
    COLS = (['family', 'feature', 'contrast', 'model', 'model_formula', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'bh_family',
             'cell_F_estimate', 'cell_F_p', 'cell_F_q_bh', 'sign_agrees_with_cell_F', 'sig_biopsy', 'sig_cell_F', 'var_patient', 'var_patient_timepoint', 'var_residual',
             'singular', 'converged', 'min_tregs_per_biopsy', 'n_patients', 'n_biopsies', 'n_tregs'] + list(NN) + ['feature_key'])
    R['_m'] = R.model.map({v: i for i, v in enumerate(MODS.values())})
    write_sig3(order_rows(R, FCOLS, COLS, by=('_m', '_c', '_f')), O54)


# ======================================================================================================================
# F6A_suppression_score_biopsy: biopsy-level LMMs of the Treg Suppression Score and suppressive effector genes (Fig 6A; Supp Table 11)
# ======================================================================================================================
if _run('F6A_suppression_score_biopsy'):
    EFFECTORS = ['CTLA4', 'IL10', 'TGFB1', 'EBI3', 'LAG3', 'ENTPD1', 'NT5E', 'LRRC32', 'TIGIT', 'GZMB', 'IL2RA', 'HAVCR2']
    NOTE = {'IL2RA': 'IL2RA (CD25) is also an IL-2/STAT5 target and reflects IL-2 signalling, not only IL-2 consumption as a suppression mechanism'}
    treg = load_tregs()
    F = pd.read_parquet(CACHE)[['Patient', 'Remission_status', 'Treatment', 'sample_id', 'pseudotime']].loc[treg.obs_names].copy()
    F['supp_score'] = suppression_score(treg)
    GENES = [g for g in EFFECTORS if g in treg.var_names]; ABSENT = [g for g in EFFECTORS if g not in treg.var_names]
    Xg = treg[:, GENES].X; Xg = Xg.toarray() if hasattr(Xg, 'toarray') else np.asarray(Xg)
    for i, g in enumerate(GENES): F[f'gene_{g}'] = Xg[:, i]
    del treg

    FAM_S, FAM_G = 'Treg Suppression Score', f'Suppressive effector genes ({len(GENES)})'
    FEATS = [(FAM_S, 'supp_score', 'Treg Suppression Score (early − late, score_genes, biopsy WTX gene sets)')] + \
            [(FAM_G, f'gene_{g}', f'{g} (log-normalised expression)') for g in GENES]
    FCOLS = [k for _, k, _ in FEATS]
    PB, DROP = pseudobulk(F, FCOLS)
    NN = {**group_n(PB, tregs=True), **paired_n(PB), 'n_biopsies_dropped': len(DROP)}
    KEPT = F.biopsy.isin(PB.biopsy)
    PCT = {f'gene_{g}': (100 * (F[f'gene_{g}'] > 0).mean(), 100 * (F.loc[KEPT, f'gene_{g}'] > 0).mean()) for g in GENES}
    GM = {k: {nm: PB[(PB.rem == r_) & (PB.post == p_)][k].mean() for nm, r_, p_ in GRP} for k in FCOLS}

    R = lmer_biopsy(PB, FCOLS)
    R['model'] = MOD; R['model_formula'] = FORMULA + '; no covariates'
    R = R.merge(pd.DataFrame(FEATS, columns=['family', 'feature_key', 'feature']), on='feature_key')
    R['q_bh'] = np.nan
    for _, i in R.groupby(['family', 'contrast']).groups.items():
        R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
    R['bh_family'] = R.family.map({FAM_S: 'Treg Suppression Score alone (1 test; q = P), per contrast',
                                   FAM_G: f'Suppressive effector genes ({len(GENES)} present), BH per contrast'})
    R['units'] = np.where(R.feature_key == 'supp_score', 'score units (early − late score_genes)', 'log1p(CP10k)')
    R['pct_tregs_expressing_all'] = R.feature_key.map(lambda k: PCT[k][0] if k in PCT else np.nan)
    R['pct_tregs_expressing_analysed_biopsies'] = R.feature_key.map(lambda k: PCT[k][1] if k in PCT else np.nan)
    for nm, _, _ in GRP: R[f'mean_{nm}'] = R.feature_key.map(lambda k: GM[k][nm])
    R['note'] = R.feature_key.map(lambda k: NOTE.get(k[5:], '') if k.startswith('gene_') else
                                  'primary Treg Suppression Score; high = Activated, low = Low-suppressive Tregs')
    R['effector_genes_absent'] = ', '.join(ABSENT) if ABSENT else 'none'
    R['min_tregs_per_biopsy'] = MIN_TREG
    R['n_patients'] = PB.Patient.nunique(); R['n_biopsies'] = len(PB); R['n_tregs'] = int(PB.n_treg.sum())
    for k, v in NN.items(): R[k] = v
    COLS = (['family', 'feature', 'contrast', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'bh_family', 'units',
             'pct_tregs_expressing_all', 'pct_tregs_expressing_analysed_biopsies'] + [f'mean_{nm}' for nm, _, _ in GRP] +
            ['note', 'effector_genes_absent', 'model', 'model_formula', 'var_patient', 'var_patient_timepoint', 'var_residual', 'singular', 'converged',
             'min_tregs_per_biopsy', 'n_patients', 'n_biopsies', 'n_tregs'] + list(NN) + ['feature_key'])
    write_sig3(order_rows(R, FCOLS, COLS), O55)


# ======================================================================================================================
# F6A_post_treatment_difference: post-treatment Remission vs Non-remission contrast (b_rem + b_post:rem) of the same model,
# appended to the O54 / O55 CSVs (run once after each run of the two sections above)
# ======================================================================================================================
if _run('F6A_post_treatment_difference'):
    NEWC = CON['post_difference']
    A54, A55 = pd.read_csv(O54), pd.read_csv(O55)
    if NEWC in set(A54.contrast) or NEWC in set(A55.contrast): sys.exit('post-treatment difference rows already present')
    PROGK = ['il2_stat5', 'costim_module', 'gene_TNFRSF9', 'tnfr2_pathway', 'pseudotime']          # programs family of F6A_programs_biopsy
    FCOLS = ['supp_score'] + PROGK
    treg = load_tregs()
    F = pd.read_parquet(CACHE).loc[treg.obs_names].copy(); F['supp_score'] = suppression_score(treg); del treg
    PB, _ = pseudobulk(F, FCOLS)
    R = lmer_biopsy(PB, FCOLS, cons=('post_difference',))
    R['q_bh'] = np.nan   # BH: programs family of 5; Treg Suppression Score alone (q = P)
    i = R.index[R.feature_key.isin(PROGK)]; R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
    i = R.index[R.feature_key == 'supp_score']; R.loc[i, 'q_bh'] = R.loc[i, 'p']

    # new rows: other columns copied from the feature's baseline row
    NR = R.set_index('feature_key')
    STATC = ['estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'var_patient', 'var_patient_timepoint', 'var_residual', 'singular', 'converged']

    def new_rows(A, keys):
        base = A[(A.model == MOD) & (A.contrast == CON['baseline_difference']) & A.feature_key.isin(keys)].copy()
        base['contrast'] = NEWC
        for c in STATC: base[c] = base.feature_key.map(NR[c]).map(sig3)
        if 'sig_biopsy' in base: base['sig_biopsy'] = base.q_bh < 0.05
        for c in ['cell_F_estimate', 'cell_F_p', 'cell_F_q_bh', 'sign_agrees_with_cell_F']:
            if c in base: base[c] = np.nan
        if 'sig_cell_F' in base: base['sig_cell_F'] = False
        return base
    with open(O54, 'w') as fh: pd.concat([A54, new_rows(A54, PROGK)], ignore_index=True).to_csv(fh, index=False)
    with open(O55, 'w') as fh: pd.concat([A55, new_rows(A55, ['supp_score'])], ignore_index=True).to_csv(fh, index=False)


# ======================================================================================================================
# F6A_pathways_biopsy: the six Treg suppression pathways of Fig 4G scored in TAURUS Tregs, biopsy-level LMM (Fig 6A IL-2/STAT5, TNFR2)
# ======================================================================================================================
if _run('F6A_pathways_biopsy'):
    G3 = pd.read_csv(f'{RJ}/Fig3D-I_treg_pseudotime_signals/Fig3G_pathways_along_pseudotime_stats.csv')
    G3 = G3[G3.resource == 'Suppression_panel'].reset_index(drop=True)   # the 6 Treg suppression rows of Fig 4G
    SETS = list(G3.gene_set); PATH = dict(zip(G3.gene_set, G3.pathway))
    msig = pd.read_csv(f'{D}/archive/csv/msigdb_human_snapshot.csv.gz', usecols=['geneset', 'genesymbol'])
    net = (msig[msig.geneset.isin(SETS)].drop_duplicates(['geneset', 'genesymbol'])
           .rename(columns={'geneset': 'source', 'genesymbol': 'target'})[['source', 'target']])
    tn = net[(net.source == TNFR2) & ~net.target.str.startswith('PSM')].assign(source='TNFR2_NO_PSM')           # as cell_level_lmm

    if os.path.exists(PCACHE):
        S = pd.read_parquet(PCACHE); NG = pd.read_csv(PCACHE.replace('.parquet', '_ngenes.csv'), index_col=0).n_genes_scored.to_dict()
    else:   # per-cell ULM scores on log-normalised expression (default tmin)
        import decoupler as dc
        treg = load_tregs(); treg.layers['lognorm'] = treg.X.copy()
        S = treg.obs[['Patient', 'Remission_status', 'Treatment', 'sample_id']].astype(str).copy()
        NG = {s: int(net[net.source == s].target.isin(treg.var_names).sum()) for s in SETS}
        NG['TNFR2_NO_PSM'] = int(tn.target.isin(treg.var_names).sum())
        dc.mt.ulm(data=treg, net=net, layer='lognorm'); A = dc.pp.get_obsm(treg, key='score_ulm')
        for s in SETS: S[s] = np.asarray(A[:, s].X).ravel()
        dc.mt.ulm(data=treg, net=tn, layer='lognorm'); S['TNFR2_NO_PSM'] = np.asarray(dc.pp.get_obsm(treg, key='score_ulm').X).ravel()
        S.attrs = {}; S.to_parquet(PCACHE)
        with open(PCACHE.replace('.parquet', '_ngenes.csv'), 'w') as fh:
            pd.Series(NG, name='n_genes_scored').rename_axis('gene_set').to_csv(fh)
    NSET = net.groupby('source').size().to_dict(); NSET['TNFR2_NO_PSM'] = len(tn)

    FCOLS = SETS + ['TNFR2_NO_PSM']
    PB, DROP = pseudobulk(S.copy(), FCOLS)
    R = lmer_biopsy(PB, FCOLS)
    R['model'] = MOD; R['model_formula'] = FORMULA
    FAM3G, FAMCHK = 'Fig 3G Treg suppression pathways (6)', 'biopsy_level_prepost.csv definitions (not in BH family)'
    R['family'] = np.where(R.feature_key.isin(SETS), FAM3G, FAMCHK)
    R['pathway'] = R.feature_key.map(PATH).fillna('TNFR2 non-canonical NF-κB (Reactome, without PSM genes; as biopsy_level_prepost.csv)')
    R['q_bh'] = np.nan
    for c, i in R[R.family == FAM3G].groupby('contrast').groups.items():
        R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
    R['bh_family'] = np.where(R.family == FAM3G, 'BH over the 6 Fig 3G pathways, per contrast', 'none')
    R['n_genes_in_set'] = R.feature_key.map(NSET); R['n_genes_scored'] = R.feature_key.map(NG)
    # P of the matching features in biopsy_level_prepost.csv, side by side
    B = pd.read_csv(O54); B = B[B.model == MOD]
    R['fig6A_existing_p'] = np.nan
    for mine, theirs in [(IL2, 'il2_stat5'), (TNFR2, 'tnfr2_pathway'), ('TNFR2_NO_PSM', 'tnfr2_pathway')]:
        m = B[B.feature_key == theirs].set_index('contrast').p
        i = R.feature_key == mine; R.loc[i, 'fig6A_existing_p'] = R.loc[i, 'contrast'].map(m)
    R['fig6A_existing_feature'] = R.feature_key.map({IL2: 'IL-2/STAT5 (identical gene set)', TNFR2: 'TNFR2 without PSM genes (gene set differs)',
                                                     'TNFR2_NO_PSM': 'TNFR2 without PSM genes (identical gene set)'})
    # CD Treg pseudotime result (Fig 4G) side by side
    g = G3.set_index('gene_set')
    R['fig3G_direction_late_vs_early'] = R.feature_key.map(lambda s: ('down' if g.loc[s, 'median_delta_last_vs_first_quintile'] < 0 else 'up') if s in g.index else np.nan)
    R['fig3G_median_delta_late_vs_early'] = R.feature_key.map(g.median_delta_last_vs_first_quintile)
    R['fig3G_q_late_vs_early'] = R.feature_key.map(g.q_last_vs_first_quintile_BH)
    R['fig3G_q_bin_effect'] = R.feature_key.map(g.q_bin_effect_mixed_model_BH)
    R['fig3G_trajectory_shape'] = R.feature_key.map(g.trajectory_shape)
    R['fig3G_test'] = np.where(R.family == FAM3G, 'last vs first pseudotime quintile, patient-paired Wilcoxon (25 CD patients), BH over the panel', '')
    NN = {'n_patients': PB.Patient.nunique(), 'n_biopsies': len(PB), 'n_tregs': int(PB.n_treg.sum()), 'n_biopsies_dropped': len(DROP), **group_n(PB)}
    for k, v in NN.items(): R[k] = v
    R['min_tregs_per_biopsy'] = MIN_TREG
    COLS = (['family', 'pathway', 'contrast', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'bh_family',
             'fig3G_direction_late_vs_early', 'fig3G_median_delta_late_vs_early', 'fig3G_q_late_vs_early', 'fig3G_q_bin_effect', 'fig3G_trajectory_shape',
             'fig3G_test', 'fig6A_existing_p', 'fig6A_existing_feature', 'model', 'model_formula', 'var_patient', 'var_patient_timepoint', 'var_residual',
             'singular', 'converged', 'n_genes_in_set', 'n_genes_scored', 'min_tregs_per_biopsy'] + list(NN) + ['feature_key'])
    write_sig3(order_rows(R, FCOLS, COLS), O65)


# ======================================================================================================================
# F6A_tnfrsf4_biopsy: TNFRSF4 (OX40; Fig 6A panel) and TNFRSF18 (GITR; CSV only) in Tregs, biopsy-level LMM
# BH for TNFRSF4: the 4 Fig 6A activation programs (IL-2/STAT5, TNFR2 full set, receptor module, TNFRSF4), per contrast
# ======================================================================================================================
if _run('F6A_tnfrsf4_biopsy'):
    GENES = {'TNFRSF4': 'TNFRSF4 (OX40)', 'TNFRSF18': 'TNFRSF18 (GITR)'}
    treg = load_tregs()
    F = pd.read_parquet(CACHE).loc[treg.obs_names].copy()
    for g in GENES: F[f'gene_{g}'] = gene_values(treg, g)
    del treg
    FCOLS = [f'gene_{g}' for g in GENES]
    PB, DROP = pseudobulk(F, FCOLS)
    FK = F[F.biopsy.isin(PB.biopsy)]
    NN = group_n(PB)
    R = lmer_biopsy(PB, FCOLS, cons=tuple(CON))

    A54 = pd.read_csv(O54)
    G = pd.read_csv(O65); G = G[(G.model == MOD) & (G.family == 'Fig 3G Treg suppression pathways (6)') & G.feature_key.isin([IL2, TNFR2])]
    CM = A54[(A54.model == MOD) & (A54.feature_key == 'costim_module')]
    R['q_bh'] = np.nan; R['bh_family'] = ''
    for k, cname in CON.items():
        if k == 'post_difference': continue
        others = pd.concat([G[G.contrast == cname][['feature_key', 'p']], CM[CM.contrast == cname][['feature_key', 'p']]])
        i = R.index[(R.feature_key == 'gene_TNFRSF4') & (R.contrast_key == k)]
        ps = np.r_[others.p.values, R.loc[i, 'p'].values]; R.loc[i, 'q_bh'] = multipletests(ps, method='fdr_bh')[1][-1]
        R.loc[i, 'bh_family'] = 'Fig 6A activation programs (4: IL-2/STAT5, TNFR2 full set, receptor module, TNFRSF4), per contrast'
    R.loc[(R.feature_key == 'gene_TNFRSF4') & (R.contrast_key == 'post_difference'), 'bh_family'] = 'none (no post-treatment difference rows for the other 3 programs)'
    R.loc[R.feature_key == 'gene_TNFRSF18', 'bh_family'] = 'none (CSV-only companion, unadjusted P)'

    # % of Tregs expressing, observed group means of the biopsy means
    EX = {}
    for g in GENES:
        c = f'gene_{g}'; e = {'pct_tregs_expressing_all': 100 * (FK[c] > 0).mean()}
        for nm, r_, p_ in GRP:
            s = FK[(FK.rem == r_) & (FK.post == p_)]; e[f'pct_tregs_expressing_{nm}'] = 100 * (s[c] > 0).mean(); e[f'n_tregs_{nm}'] = len(s)
            b = PB[(PB.rem == r_) & (PB.post == p_)][c]; e[f'mean_biopsy_mean_{nm}'] = b.mean(); e[f'sem_biopsy_mean_{nm}'] = b.std(ddof=1) / np.sqrt(len(b))
        EX[c] = e
    R = R.join(pd.DataFrame.from_dict(EX, orient='index'), on='feature_key')
    R['feature'] = R.feature_key.str[5:].map(GENES)
    R['model'] = MOD; R['model_formula'] = FORMULA
    R['value'] = 'per-Treg log-normalised expression (log1p CP10k), averaged per biopsy'
    R['min_tregs_per_biopsy'] = MIN_TREG; R['n_patients'] = PB.Patient.nunique(); R['n_biopsies'] = len(PB); R['n_tregs'] = len(FK)
    for k, v in NN.items(): R[k] = v
    R['n_biopsies_dropped'] = len(DROP)
    EXC = [c for c in R.columns if c.startswith(('pct_tregs', 'mean_biopsy', 'sem_biopsy', 'n_tregs_'))]
    COLS = (['feature', 'contrast', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'bh_family', 'model', 'model_formula', 'value']
            + EXC + ['var_patient', 'var_patient_timepoint', 'var_residual', 'singular', 'converged', 'min_tregs_per_biopsy', 'n_patients', 'n_biopsies',
                     'n_tregs'] + list(NN) + ['n_biopsies_dropped', 'feature_key'])
    write_sig3(order_rows(R, FCOLS, COLS, by=('_f', '_c')), O66)


# ======================================================================================================================
# F6A_treg_state_lines: Fig 6A figure (statistics read from the biopsy-level CSVs, not refitted)
# ======================================================================================================================
if _run('F6A_treg_state_lines'):
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from decimal import Decimal, ROUND_HALF_UP
    STEM = f'{RJ}/Fig6_treg_antitnf_taurus/fig6_treg_state_lines'
    CK = {'rem': 'change_remission', 'non': 'change_nonremission', 'diff': 'differential_change', 'base': 'baseline_difference', 'post': 'post_difference'}

    # panels: (key, source, panel title, y label)
    PANELS = [('supp_score', 's55', 'Treg Suppression Score', 'Treg Suppression Score'),
              (IL2, 's65', 'IL-2/STAT5 Signaling', 'Activity (ULM score)'),
              (TNFR2, 's65', 'TNFR2 Non-Canonical\nNF-κB', 'Activity (ULM score)'),
              ('costim_module', 's54', 'Co-Stimulatory/Suppressive\nReceptor Module', 'Module score'),
              ('gene_TNFRSF4', 's66', 'TNFRSF4 (OX40)', 'Log-normalised expression')]
    PK = [p[0] for p in PANELS]; AK = PK[1:]                               # plotted panels; the 4 activation programs

    S5 = pd.read_csv(O55); S5 = S5[(S5.model == MOD) & (S5.feature_key == 'supp_score')]
    G = pd.read_csv(O65); G = G[(G.model == MOD) & (G.family == 'Fig 3G Treg suppression pathways (6)') & G.feature_key.isin([IL2, TNFR2])]
    B = pd.read_csv(O54); B = B[(B.model == MOD) & (B.feature_key == 'costim_module')]
    S6 = pd.read_csv(O66); S6 = S6[(S6.model == MOD) & (S6.feature_key == 'gene_TNFRSF4')]
    W = pd.concat([S5, G, B, S6], ignore_index=True).pivot(index='feature_key', columns='contrast', values=['estimate', 'p'])

    # per-cell values -> biopsy means -> observed group means
    F = pd.read_parquet(CACHE)
    PF = pd.read_parquet(PCACHE).loc[F.index]
    F[IL2] = PF[IL2].values; F[TNFR2] = PF[TNFR2].values
    treg = load_tregs()
    F = F.loc[treg.obs_names].copy()
    F['supp_score'] = suppression_score(treg); F['gene_TNFRSF4'] = gene_values(treg, 'TNFRSF4'); del treg
    PB, _ = pseudobulk(F, PK)
    GK = {('rem', 'Pre'): 'Remission', ('rem', 'Post'): 'Remission', ('non', 'Pre'): 'Non_Remission', ('non', 'Post'): 'Non_Remission'}
    NP = {g: PB[(PB.Remission_status == rs) & (PB.Treatment == g[1])].Patient.nunique() for g, rs in GK.items()}
    rows = []
    for k in PK:
        r = dict(key=k)
        for g, rs in GK.items():
            v = PB.loc[(PB.Remission_status == rs) & (PB.Treatment == g[1]), k]
            r[f'{g[0]}_{g[1]}'] = v.mean(); r[f'{g[0]}_{g[1]}_se'] = v.std(ddof=1) / np.sqrt(len(v))
        for c, ck in CK.items():
            has = ('p', CON[ck]) in W.columns and pd.notna(W.loc[k, ('p', CON[ck])])
            r[f'{c}_est'] = W.loc[k, ('estimate', CON[ck])] if has else np.nan; r[f'{c}_p'] = W.loc[k, ('p', CON[ck])] if has else np.nan
        rows.append(r)
    P = pd.DataFrame(rows).set_index('key')
    for c in ['rem', 'non', 'diff', 'base']:   # BH q over the 4 activation programs, per contrast (reported in the legend)
        P.loc[AK, f'{c}_q'] = multipletests(P.loc[AK, f'{c}_p'].values, method='fdr_bh')[1]
    print(P[[f'{c}_{s}' for c in ['rem', 'non', 'diff', 'base'] for s in ['est', 'p', 'q']]].to_string(), flush=True)

    def sig2(p):   # 2 significant figures, half-up on the CSV's 3-significant-figure value
        d = Decimal(repr(float(p))); q = Decimal(1).scaleb(d.adjusted() - 1)
        return str(d.quantize(q, rounding=ROUND_HALF_UP))

    def fmt_p(p, s='P'):
        return f'{s} < 0.001' if p < 0.001 else f'{s} = {sig2(float(f"{p:.3g}"))}'

    # figure: Suppression State (row 1), Activation Programs (row 2); fixed grid, saved with a tight bounding box
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14, 'xtick.labelsize': 13, 'ytick.labelsize': 13,
                         'axes.titleweight': 'normal', 'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0})
    BLUE, RED, K, GR = '#2166ac', '#b2182b', '#111111', '#8c8c8c'
    DX = 0.06; FS_P = 12; DPI = 300; XLIM = [-0.45, 1.45]
    WIN, HIN = 15.4, 9.4; fy = lambda inch: inch / HIN
    fig = plt.figure(figsize=(WIN, HIN), dpi=DPI)
    gs = fig.add_gridspec(2, 4, left=1.05 / WIN, right=1 - 0.25 / WIN, bottom=fy(0.80), top=fy(8.00), wspace=0.44, hspace=0.67)
    axs = [fig.add_subplot(gs[0, 0])] + [fig.add_subplot(gs[1, j]) for j in range(4)]
    for ax, (k, src, title, ylab) in zip(axs, PANELS):
        r = P.loc[k]
        for g, col, dx in [('non', RED, DX), ('rem', BLUE, -DX)]:
            ax.errorbar(np.array([0, 1]) + dx, [r[f'{g}_Pre'], r[f'{g}_Post']], yerr=[r[f'{g}_Pre_se'], r[f'{g}_Post_se']], color=col, lw=2.0,
                        marker='o', ms=6.5, capsize=3, elinewidth=1.4, capthick=1.4, zorder=3 if g == 'rem' else 2)
        ax.set_xlim(*XLIM); ax.set_xticks([0, 1]); ax.set_xticklabels(['Pre-\ntreatment', 'Post-\ntreatment'], linespacing=1.0)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylabel(ylab, fontsize=14, labelpad=4, linespacing=1.1)
        ax.set_title(title, fontsize=15, pad=8, linespacing=1.1)
        lo_ = min(r[f'{g}_{t}'] - r[f'{g}_{t}_se'] for g in ['rem', 'non'] for t in ['Pre', 'Post'])
        hi_ = max(r[f'{g}_{t}'] + r[f'{g}_{t}_se'] for g in ['rem', 'non'] for t in ['Pre', 'Post'])
        ax.set_ylim(lo_ - 0.10 * (hi_ - lo_), hi_ + 0.45 * (hi_ - lo_))   # headroom for the change P values
        for i, (g, col) in enumerate([('rem', BLUE), ('non', RED)]):
            ax.text(0.5, 0.98 - 0.09 * i, fmt_p(r[f'{g}_p']), transform=ax.transAxes, ha='center', va='top', fontsize=FS_P, color=col)

    # group labels + brackets above each row (offset scaled by the number of title lines)
    for i0, i1, lab in [(0, 0, 'Suppression State'), (1, 4, 'Activation Programs')]:
        b0, b1 = axs[i0].get_position(), axs[i1].get_position()
        nl = 1 + max(PANELS[i][2].count('\n') for i in range(i0, i1 + 1))
        y = b0.y1 + fy(0.33 + 0.25 * nl); x0, x1 = b0.x0 - 0.68 / WIN, b1.x1
        fig.add_artist(Line2D([x0, x0, x1, x1], [y - fy(0.13), y, y, y - fy(0.13)], transform=fig.transFigure, color=GR, lw=1.2))
        fig.text((x0 + x1) / 2, y + fy(0.05), lab, ha='center', va='bottom', fontsize=16, color=K)

    # colour key + one-line P key right of the Suppression State panel
    HD = [Line2D([], [], color=BLUE, lw=2.0, marker='o', ms=6.5, label=f"Remission (n = {NP[('rem', 'Pre')]} patients)"),
          Line2D([], [], color=RED, lw=2.0, marker='o', ms=6.5, label=f"Non-remission (n = {NP[('non', 'Pre')]} pre-, {NP[('non', 'Post')]} post-treatment)")]
    kx = axs[0].get_position().x1 + 0.75 / WIN; ktop = axs[0].get_position().y1 + fy(0.05)
    fig.legend(handles=HD, loc='upper left', bbox_to_anchor=(kx, ktop), ncol=1, frameon=False, fontsize=14, handlelength=1.8,
               borderaxespad=0, labelspacing=0.6)
    fig.text(kx + 0.004, ktop - fy(0.95), 'Blue / red: change after treatment within group (biopsy-level mixed model)', ha='left', va='top', fontsize=13, color=K)
    fig.savefig(f'{STEM}.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08); plt.close(fig)
    log('wrote', f'{STEM}.png')
