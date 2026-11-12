#!/usr/bin/env python
"""TAURUS CD scRNA-seq (anti-TNF): Treg Suppression Score and Treg activation programs before and after treatment.

Panels made
  Fig 6A  Treg Suppression Score, IL-2/STAT5 signaling, TNFR2 non-canonical NF-kB, co-stimulatory/suppressive receptor module and
          TNFRSF4, Remission vs Non-remission, pre vs post (biopsy-level LMM with patient and patient x timepoint random effects)
  Supp Table 11 (biopsy-level models)

Inputs
  /share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny/raw_norm_counts_anno.h5ad.gz (01_taurus_preprocessing.ipynb, 02_taurus_aucell_annotation.R)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/D_shared_inputs/treg_*_signature_genes.csv
    (Treg Suppression Score gene sets, cosmx_cd/10_treg_suppression_score.py)
  Treg pseudotime and TF tables of the CD scRNA-seq atlas (cd_scrna/12_treg_pseudotime.py, cd_scrna/14_treg_pathway_tf_activity.py):
    cd_treg_rebuild_rowfix/panels/extended_data_fig8/ext8_tf_activity_along_pseudotime.csv (TF list),
    Fig3D-I_treg_pseudotime_signals/Fig3G_pathways_along_pseudotime_stats.csv and fig3F_marker_dotplot.csv,
    SuppFig8A-I_treg_pseudotime_robustness/ext8_tf_activity_along_pseudotime_drivers.csv
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/ (CellChatDB, CollecTRI snapshots) and
    /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/archive/csv/msigdb_human_snapshot.csv.gz (MSigDB snapshot)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig6c_taurus/ (Treg pseudotime, cell-level LMM, patient-level
    Treg Suppression Score tests)
  /share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/ (biopsy-level model CSVs = Supp Table 11; per-cell cache/)
  /share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/fig6_treg_state_lines.png and legend (Fig 6A)

Run order of sections
  treg_pseudotime -> cell_level_lmm -> suppression_score_patient_level -> F6A_programs_biopsy -> F6A_suppression_score_biopsy
  -> F6A_post_treatment_difference -> F6A_pathways_biopsy -> F6A_tnfrsf4_biopsy -> F6A_treg_state_lines
  The biopsy-level sections read the per-cell Treg table supplementary/cache/Treg.parquet and supplementary/cell_level_lmm_results.csv:
  run cell_level_lmm with `--cache <supplementary>/cache` and place its cell_level_lmm_results.csv in supplementary/.

Notes
  Model sections call Rscript (lme4 / lmerTest) through temporary R scripts.

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
# Section treg_pseudotime: diffusion pseudotime of TAURUS Tregs, per cell and per patient x timepoint
# Panels: Fig 6A inputs (TAURUS Treg diffusion pseudotime)
# ======================================================================================================================
if _run('treg_pseudotime'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/taurus_fig6c/treg_pseudotime.py'   # working folder of this step (relative paths resolve from here)
    """Diffusion pseudotime of TAURUS Tregs (aucell_cell_type_short == 'Treg').
Treg subset; 3,000 HVGs (seurat_v3 on counts); scale; PCA (50 PCs); Harmony on Patient; neighbour graph (15 neighbours on all
Harmony PCs); Leiden 0.4; suppressive score (FOXP3, CTLA4, BATF, IL2RA, CD27); root = first cell of the cluster with the highest
mean suppressive score; diffusion map and DPT; pseudotime scaled 0-1.
Outputs (panels/fig6c_taurus/): fig6B_treg_pseudotime_per_cell.csv (per-cell pseudotime, read by cell_level_lmm and
suppression_score_patient_level) and fig6B_patient_pseudotime.csv (per-patient pre / post means and their change)."""
    import scanpy as sc, scanpy.external as sce, pandas as pd, numpy as np
    T='/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'; P='/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig6c_taurus'
    cd=sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz')
    treg=cd[cd.obs['aucell_cell_type_short']=='Treg'].copy(); del cd
    sc.pp.highly_variable_genes(treg, n_top_genes=3000, flavor='seurat_v3', layer='counts')
    sc.pp.scale(treg, max_value=10)
    sc.tl.pca(treg, n_comps=50, use_highly_variable=True)
    sce.pp.harmony_integrate(treg, key='Patient')
    # neighbour graph for Leiden / DPT: 15 neighbours on all Harmony PCs
    sc.pp.neighbors(treg, use_rep='X_pca_harmony', n_neighbors=15)
    sc.tl.leiden(treg, resolution=0.4, key_added='leiden_treg_full', random_state=42)
    sc.tl.score_genes(treg, ['FOXP3','CTLA4','BATF','IL2RA','CD27'], score_name='suppressive_score')
    root_cluster=treg.obs.groupby('leiden_treg_full')['suppressive_score'].mean().idxmax()
    root_cell=treg.obs[treg.obs['leiden_treg_full']==root_cluster].index[0]
    treg.uns['iroot']=treg.obs_names.get_loc(root_cell)
    sc.tl.diffmap(treg); sc.tl.dpt(treg)
    dmin,dmax=treg.obs['dpt_pseudotime'].min(),treg.obs['dpt_pseudotime'].max()
    treg.obs['dpt_pseudotime_norm']=(treg.obs['dpt_pseudotime']-dmin)/(dmax-dmin)
    print(treg.obs['dpt_pseudotime_norm'].describe().round(3))
    treg.obs[['Patient','Remission_status','Treatment','leiden_treg_full','dpt_pseudotime','dpt_pseudotime_norm']].to_csv(f'{P}/fig6B_treg_pseudotime_per_cell.csv')
    pt_stats=treg.obs.groupby(['Patient','Remission_status','Treatment'],observed=True).agg(mean_dpt=('dpt_pseudotime_norm','mean')).reset_index().dropna()
    pt_pivot=pt_stats.pivot_table(index=['Patient','Remission_status'],columns='Treatment',values='mean_dpt').reset_index().dropna()
    pt_pivot['delta']=pt_pivot['Post']-pt_pivot['Pre']; pt_pivot.to_csv(f'{P}/fig6B_patient_pseudotime.csv',index=False)


# ======================================================================================================================
# Section cell_level_lmm: per-cell Treg and inflammatory-axis features and cell-level mixed models of their pre/post change
# Panels: Fig 6A inputs (per-cell Treg feature table; cell-level model columns of Supp Table 11)
# Note: optional argument --cache <dir>: write (or reuse) the per-cell feature tables there; the biopsy-level sections read <dir>/Treg.parquet
# ======================================================================================================================
if _run('cell_level_lmm'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/29_taurus_cell_level_lmm.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: per-cell feature tables and cell-level linear mixed models with patient and sample random
effects for the Treg programs, 28 Treg TFs and the IL-6 / IL-1 / OSM / IL-11 / TNF inflammatory axis.

Every cell is modelled with random intercepts for patient and for sample (sample = patient x timepoint), so the pre/post contrasts
are tested against sample-to-sample variation, not against cell-to-cell variation (as in the TAURUS study's MAST models with a
treatment x remission interaction; Thomas et al. 2024 Nat Immunol, PMID 39438660).

Features (per cell, X = log1p(CP10k)):
  Treg programs (5 features, Tregs = aucell_cell_type_short == 'Treg'):
    Treg diffusion pseudotime: dpt_pseudotime_norm of treg_pseudotime (fig6B_treg_pseudotime_per_cell.csv);
    Hallmark IL-2/STAT5: decoupler ULM on log-normalised expression with a 6-set MSigDB net;
    co-stimulation / suppressive Treg receptor module: 16 CellChatDB receptors (>= 5% of Tregs), scanpy score_genes;
    TNFR2 non-canonical NF-kB: Reactome without PSM genes, ULM; TNFRSF9 log-normalised expression.
  TF family (28 features): CollecTRI ULM, tmin = 5: STAT5A, STAT5B, FOXP3, RUNX1, NFKB1, NFKB2, RELA, RELB, REL, PRDM1, IKZF2,
    BATF and the 20 TFs of the CD atlas Treg pseudotime TF table (ext8_tf_activity_along_pseudotime.csv).
  Inflammatory-axis family (18 features): ligands OSM | Myeloid, IL6 | Stromal, IL6 | Myeloid, IL1A | Myeloid, IL1B | Myeloid,
    IL11 | Fibroblasts, TNF | Myeloid, TNF | T Cells, and receptors OSMR, IL11RA, IL6ST, IL1R1, TNFRSF1A in fibroblasts and in
    OGN+RSPO3+ Fib; per-cell log-normalised expression in the feature's compartment. SELENOP+ Fib excluded from fibroblasts; the
    SELENOP+ Mac subtypes are in myeloid.
Model (per feature; all cells of the feature's compartment; all samples including the non-remitter sampled at baseline only; no
minimum cell count):
  y ~ post * remission + libsize_z + (1 | patient) + (1 | sample)
  statsmodels MixedLM, REML: groups = patient, re_formula = '1', vc_formula = {'sample': '0 + C(sample)'}; post and remission 0/1;
  libsize_z = log total counts z-scored across the cells analysed; y is standardised before fitting (estimates back-transformed).
  Contrasts: change in remission = b_post + b_post:remission; change in non-remission = b_post; differential change (remission minus
  non-remission) = b_post:remission.
  Primary P: t test on the contrast with between-within df (sample-level contrasts: n_samples - n_patients - 2 = 13). The Wald z P
  is also reported. The same model is fitted in R lme4 / lmerTest with Satterthwaite df (estimates must agree with statsmodels).
  BH within each family (Treg programs; TFs; inflammatory axis), separately for each model and contrast.
  Convergence: every fit is run from 5 optimiser sequences (lbfgs; powell+lbfgs; nm+lbfgs; bfgs; cg) and the converged fit with the
  highest REML log-likelihood is kept (lbfgs alone can report convergence at an inferior optimum).
  If no sequence converges for the nested fit, the primary row falls back to the patient-only random-intercept fit and is flagged.
Sensitivity / diagnostics:
  (a) patient random intercept only (no sample term; cells independent within patient), Wald z;
  (c) + biopsy random intercept (patient / sample / biopsy; sample_id): each TAURUS CD patient x timepoint holds 1-5 biopsies from
      different sites (terminal ileum, ascending / descending / sigmoid colon, rectum);
  (d) as (c) on site-matched biopsies only (sites sampled both before and after treatment in that patient; the TAURUS longitudinal
      rule). Drops the baseline-only non-remitter and the unmatched sites of CD2, CD4 and CD5;
  (e) the nested model without libsize_z (read by the biopsy-level sections as the cell-level model F);
  (f) diagnostic: pooled-variance paired t test on the per-patient change in sample means (OLS of change on remission; paired
      patients only; df = n_paired - 2 = 13).
  Model labels in the CSV: A primary, B (a), D (c), E (d), F (e), G (f).
Output panels/fig6c_taurus/cell_level_lmm/cell_level_lmm_results.csv (3 significant figures).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>
     optional: --cache <dir>  write / reuse the per-cell feature tables (Treg.parquet etc.; all statistics are recomputed)."""
    import os, sys, warnings, subprocess, tempfile, shutil, time
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import numpy as np, pandas as pd
    from scipy import stats
    import statsmodels.api as sm, statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests

    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
    P = f'{D}/panels/fig6c_taurus'; OUT = f'{P}/cell_level_lmm'; os.makedirs(OUT, exist_ok=True)
    KEYS = ['Patient', 'Remission_status', 'Treatment']; OGN = 'OGN+RSPO3+ Fib'
    CACHE = sys.argv[sys.argv.index('--cache') + 1] if '--cache' in sys.argv else None
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 400)
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:6.0f}s]', *a, flush=True)

    # ================================================================ feature definitions
    FAM_6 = 'Treg programs (5)'; FAM_TF = 'Treg TFs (28)'; FAM_INF = 'Inflammatory axis (18)'
    SET_B = ['STAT5A', 'STAT5B', 'FOXP3', 'RUNX1', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL', 'PRDM1', 'IKZF2', 'BATF']
    E8 = pd.read_csv(f'{D}/panels/extended_data_fig8/ext8_tf_activity_along_pseudotime.csv'); assert len(E8) == 20
    TFS = list(dict.fromkeys(SET_B + list(E8.tf))); assert len(TFS) == 28
    MYE = ['Mo-Mac', 'Inf Mo-Mac', 'Trans cDC2/Mac', 'Mac S+XS-', 'Mac S+M+S+', 'Mac S+SG+', 'Mac S+M+P+', 'CD9+ Mac', 'Neutrophil', 'cDC1',
           'cDC2', 'CCR7+ DC']
    FIB = ['VSTM2A+ Crypt Top Fib', 'ADAMDEC1+ Activ Fib', 'ADAMDEC1+ Fib', 'OGN+RSPO3+ Fib', 'T-interact Fib', 'Activ Fib', 'HHIP+ Myofib',
           'PLA2G2A+ ECM Fib', 'GREM2+ Myofib']                                                                           # no SELENOP+ Fib
    A_FEAT = [('OSM', 'Myeloid'), ('IL6', 'Stromal (Fib + Endo)'), ('IL6', 'Myeloid'), ('IL1A', 'Myeloid'), ('IL1B', 'Myeloid'),
              ('IL11', 'Fibroblasts'), ('TNF', 'Myeloid'), ('TNF', 'T Cells')]
    B_FEAT = [(g, c) for c in ['Fibroblasts', 'OGN⁺RSPO3⁺ Fib'] for g in ['OSMR', 'IL11RA', 'IL6ST', 'IL1R1', 'TNFRSF1A']]
    fid = lambda g, c: f'{g} | {c}'
    SAFE = lambda s: ''.join(ch if ch.isalnum() else '_' for ch in s)

    # FEATURES: (family, figure, column, label, compartment table, type)
    FEATURES = [(FAM_6, 'Treg program', 'pseudotime', 'Treg pseudotime', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'il2_stat5', 'IL-2/STAT5 signaling (Hallmark, ULM)', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'costim_module', 'Co-stimulation / suppressive Treg receptor module', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'tnfr2_pathway', 'TNFR2 non-canonical NF-κB (Reactome, ULM)', 'Treg', 'score'),
                (FAM_6, 'Treg program', 'gene_TNFRSF9', 'TNFRSF9 (4-1BB)', 'Treg', 'gene')]
    FEATURES += [(FAM_TF, 'Treg TF', f'tf_{tf}', f'{tf} activity (CollecTRI, ULM)', 'Treg', 'score') for tf in TFS]
    FEATURES += [(FAM_INF, 'inflammatory ligand', fid(g, c), fid(g, c), c, 'gene') for g, c in A_FEAT]
    FEATURES += [(FAM_INF, 'inflammatory receptor', fid(g, c), fid(g, c), c, 'gene') for g, c in B_FEAT]
    assert len(FEATURES) == 5 + 28 + 18


    # ================================================================ per-cell features
    def build_tables():
        import scanpy as sc, decoupler as dc
        ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); log('loaded', ad.shape)
        ct = ad.obs['aucell_cell_type_short'].astype(str); cat = ad.obs['aucell_cell_category_from_type'].astype(str); minor = ad.obs['minor'].astype(str)
        assert np.allclose(np.expm1(ad.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)                       # X = log1p(CP10k)
        # library size = obs total_counts (QC total over all genes before gene filtering; equals n_counts); the stored counts layer
        # (28,972 genes) sums to >= 99.8% of it
        rat = np.asarray(ad.layers['counts'][:20000].sum(1)).ravel() / ad.obs['total_counts'].values[:20000]
        assert np.allclose(ad.obs['n_counts'], ad.obs['total_counts']) and rat.min() > 0.99 and rat.max() <= 1 + 1e-9, (rat.min(), rat.max())
        meta = ad.obs[KEYS + ['sample_id', 'Site']].astype(str).assign(total_counts=ad.obs['total_counts'].astype(float).values)
        TB = {}
        # ---- Tregs
        treg = ad[(ct == 'Treg').values].copy(); treg.layers['lognorm'] = treg.X.copy(); F = meta.loc[treg.obs_names].copy()
        msig = pd.read_csv(f'{D}/archive/csv/msigdb_human_snapshot.csv.gz', usecols=['geneset', 'genesymbol'])
        SETS = ['HALLMARK_IL2_STAT5_SIGNALING', 'REACTOME_RUNX1_AND_FOXP3_CONTROL_THE_DEVELOPMENT_OF_REGULATORY_T_LYMPHOCYTES_TREGS',
                'GOBP_T_CELL_TOLERANCE_INDUCTION', 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION']
        net = msig[msig.geneset.isin(SETS)].drop_duplicates().rename(columns={'geneset': 'source', 'genesymbol': 'target'})
        dc.mt.ulm(data=treg, net=net, layer='lognorm'); F['il2_stat5'] = np.asarray(dc.pp.get_obsm(treg, key='score_ulm')[:, 'HALLMARK_IL2_STAT5_SIGNALING'].X).ravel()
        tn = msig[msig.geneset == 'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY'].drop_duplicates().rename(columns={'geneset': 'source', 'genesymbol': 'target'})
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
        assert set(A.var_names) == set(TFS), sorted(set(TFS) - set(A.var_names))
        for tf in TFS: F[f'tf_{tf}'] = np.asarray(A[:, tf].X).ravel()
        pt = pd.read_csv(f'{P}/fig6B_treg_pseudotime_per_cell.csv', index_col=0)
        assert set(pt.index) == set(F.index) and (pt.loc[F.index, 'Patient'].astype(str).values == F.Patient.values).all()
        F['pseudotime'] = pt.loc[F.index, 'dpt_pseudotime_norm'].values
        TB['Treg'] = F; log('Treg features', F.shape)
        # ---- inflammatory-axis compartments
        TCELL = sorted(set(ct[cat == 'T']))
        assert set(MYE) == set(ct[cat == 'Myeloid']) and 'SELENOP+ Fib' in set(ct) and 'SELENOP+ Fib' not in FIB and set(FIB) <= set(ct)
        assert (minor[ct == 'Granulocyte'] == 'Mast').mean() > 0.95
        COMP = {'Myeloid': MYE, 'T Cells': TCELL, 'Stromal (Fib + Endo)': FIB + ['Endothelial'], 'Fibroblasts': FIB, 'OGN⁺RSPO3⁺ Fib': [OGN]}
        for c, types in COMP.items():
            genes = sorted({g for g, cc_ in A_FEAT + B_FEAT if cc_ == c}); m = ct.isin(types).values
            X = ad[m][:, genes].X; X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
            tab = meta.loc[ad.obs_names[m]].copy()
            for i, g in enumerate(genes): tab[fid(g, c)] = X[:, i]
            TB[c] = tab; log(c, tab.shape, 'cell types:', len(types))
        return TB


    if CACHE and os.path.exists(f'{CACHE}/tables.csv'):
        names = pd.read_csv(f'{CACHE}/tables.csv')
        TB = {r['name']: pd.read_parquet(f"{CACHE}/{r['file']}") for _, r in names.iterrows()}
        RECEPTORS = open(f'{CACHE}/costim_receptors.txt').read().strip().split(';'); log('loaded cache', {k: v.shape for k, v in TB.items()})
    else:
        TB = build_tables(); RECEPTORS = TB['Treg'].attrs['costim_receptors'].split(';')
        if CACHE:
            os.makedirs(CACHE, exist_ok=True)
            for k, v in TB.items():
                v = v.copy(); v.attrs = {}; v.to_parquet(f'{CACHE}/{SAFE(k)}.parquet')
            pd.DataFrame({'name': list(TB), 'file': [f'{SAFE(k)}.parquet' for k in TB]}).to_csv(f'{CACHE}/tables.csv', index=False)
            open(f'{CACHE}/costim_receptors.txt', 'w').write(';'.join(RECEPTORS))

    # ================================================================ models
    CONTRASTS = {'Change in remission (post − pre)': {'post': 1.0, 'post:rem': 1.0},
                 'Change in non-remission (post − pre)': {'post': 1.0},
                 'Differential change (remission − non-remission)': {'post:rem': 1.0}}
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
        if site_matched:   # TAURUS longitudinal rule: keep biopsies from sites sampled both before and after treatment in that patient
            both = d.groupby(['Patient', 'Site']).Treatment.nunique(); both = both[both == 2].index
            d = d[pd.MultiIndex.from_frame(d[['Patient', 'Site']]).isin(both)]
        d['post'] = (d.Treatment == 'Post').astype(float); d['rem'] = (d.Remission_status == 'Remission').astype(float)
        d['sample'] = d.Patient + '_' + d.Treatment
        ll = np.log(d.total_counts); d['lib_z'] = (ll - ll.mean()) / ll.std()
        mu, sy = d.y.mean(), d.y.std(); d['ys'] = (d.y - mu) / sy
        return d.reset_index(drop=True), sy


    def fit_lmm(d, sy, level, lib=True):
        """REML fit from every optimiser sequence in METHODS; keep the converged fit (finite, positive fixed-effect variances) with the
    highest REML log-likelihood. statsmodels' default lbfgs alone can stop at an inferior optimum while reporting convergence."""
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


    GT = [('rem_pre', 1.0, 0.0), ('rem_post', 1.0, 1.0), ('nonrem_pre', 0.0, 0.0), ('nonrem_post', 0.0, 1.0)]


    def counts(dd, ftype):
        o = {}
        for nm, r_, p_ in GT:
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
                         'var_biopsy': fit['var_biopsy'], 'var_residual': fit['var_residual'], 'libsize_coef': fit['lib_coef'], 'reml_loglik': fit['llf'], 'optimizer': fit['optimizer'],
                         'fit_warnings': fit['warnings'], 'fit_flag': flag,
                         'test': 'LMM REML, t with between-within df (sample-level contrast)' if np.isfinite(df_bw) else 'LMM REML, Wald z'})
        return rows


    def run_feature(spec, tab):
        from threadpoolctl import threadpool_limits
        fam, figp, col, lab, comp, ftype = spec; t1 = time.time()
        with threadpool_limits(2), warnings.catch_warnings():
            warnings.simplefilter('ignore')
            d, sy = prep(tab, col); n_pat, n_smp = d.Patient.nunique(), d['sample'].nunique(); df_bw = n_smp - n_pat - 2
            base = dict(family=fam, figure_panel=figp, feature=lab, feature_key=col, compartment=comp, feature_type=ftype,
                        n_patients=n_pat, n_samples=n_smp, n_biopsies=d.biopsy.nunique(), n_cells=len(d), **counts(d, ftype))
            fn = fit_lmm(d, sy, 'nested'); fp = fit_lmm(d, sy, 'patient'); assert fp['ok'], f'{col}: patient-only fit failed'
            flag = ''
            if not fn['ok']:
                flag = 'NESTED FIT FAILED: primary row = patient-only random intercept (fallback); ' + fn['retries']
            elif fn['var_sample'] < 1e-6 * fn['var_residual']:
                flag = 'sample variance estimated at ~0 (boundary)'
            rows = lmm_rows(fn if fn['ok'] else fp, base, M_NEST, df_bw if fn['ok'] else np.inf, flag)
            rows += lmm_rows(fp, base, M_PAT, np.inf)
            # sensitivities (c) + biopsy random intercept, (d) the same on site-matched biopsies; no fallback, failures flagged
            fb = fit_lmm(d, sy, 'biopsy'); rows += lmm_rows(fb, base, M_BIO, df_bw)
            ds, sys_ = prep(tab, col, site_matched=True); np_, ns_ = ds.Patient.nunique(), ds['sample'].nunique()
            bs = {**base, 'n_patients': np_, 'n_samples': ns_, 'n_biopsies': ds.biopsy.nunique(), 'n_cells': len(ds), **counts(ds, ftype)}
            fs = fit_lmm(ds, sys_, 'biopsy'); rows += lmm_rows(fs, bs, M_SITE, ns_ - np_ - 2)
            fl = fit_lmm(d, sy, 'nested', lib=False); rows += lmm_rows(fl, base, M_NOLIB, df_bw)
            dg, n_pair = diag_paired_t(d)
            for k, (est, se, dfd, p) in dg.items():
                rows.append({**base, 'model': M_DIAG, 'contrast': k, 'estimate': est, 'se': se, 'df': dfd, 'statistic': est / se, 'p': p,
                             'test': f'OLS of per-patient change in sample mean on remission ({n_pair} paired patients), t'})
        msg = (f'{col:30s} {comp:22s} cells {len(d):6d} samples {n_smp} | nested {"ok" if fn["ok"] else "FAILED"} {fn.get("optimizer", "-")} '
           f'| var_sample {fn.get("var_sample", np.nan):.3g} var_patient {fn.get("var_patient", np.nan):.3g} | biopsy {"ok" if fb["ok"] else "FAILED"} '
           f'site {"ok" if fs["ok"] else "FAILED"} | {time.time() - t1:.0f}s')
        return rows, msg


    from joblib import Parallel, delayed
    log(f'fitting {len(FEATURES)} features x 5 LMMs x {len(METHODS)} optimiser sequences (parallel)')
    cols_needed = lambda comp, col: TB[comp][['Patient', 'Remission_status', 'Treatment', 'sample_id', 'Site', 'total_counts', col]]
    OUTS = Parallel(n_jobs=min(len(FEATURES), 26), backend='loky', verbose=0)(
        delayed(run_feature)(f, cols_needed(f[4], f[2])) for f in FEATURES)
    ROWS = [r for rows, _ in OUTS for r in rows]
    for _, m in OUTS: log(m)
    R = pd.DataFrame(ROWS)


    # ================================================================ cross-check in R (lme4 / lmerTest, Satterthwaite df)
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
        log('running lme4 / lmerTest cross-check ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/spec.csv', f'{tmp}/r_out.csv'], check=True)
        RR = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    RR['contrast'] = RR.contrast_key.map({v: k for k, v in CKEY.items()})
    R = R.merge(RR[['feature_key', 'compartment', 'contrast', 'r_estimate', 'r_se', 'df_satterthwaite', 'p_satterthwaite', 'r_singular']]
                .assign(model=M_NEST), on=['feature_key', 'compartment', 'contrast', 'model'], how='left')
    nn = R[(R.model == M_NEST) & (R.fit_flag.fillna('').str.startswith('NESTED FIT FAILED') == False)]
    dev_est = ((nn.estimate - nn.r_estimate).abs() / nn.se).max(); dev_se = ((nn.se - nn.r_se).abs() / nn.se).max(); dev_se_smaller = ((nn.r_se - nn.se) / nn.r_se).max()
    log(f'statsmodels vs lme4: max |estimate difference| = {dev_est:.2g} SE, max relative SE difference = {dev_se:.2g}')
    assert dev_est < 0.01 and dev_se < 0.15, 'statsmodels and lme4 disagree'
    R = R.drop(columns=['r_estimate', 'r_se'])


    # ================================================================ BH within family, per model and contrast
    R['q_bh'] = np.nan; R['q_bh_wald_z'] = np.nan; R['q_satterthwaite'] = np.nan
    for (fam, mod, con), i in R.groupby(['family', 'model', 'contrast']).groups.items():
        i = i[R.loc[i, 'p'].notna().values]
        R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
        if mod == M_NEST:
            R.loc[i, 'q_bh_wald_z'] = multipletests(R.loc[i, 'p_wald_z'], method='fdr_bh')[1]
            R.loc[i, 'q_satterthwaite'] = multipletests(R.loc[i, 'p_satterthwaite'], method='fdr_bh')[1]
    R['estimate_definition'] = np.where(R.model == M_DIAG, 'per-patient change in sample mean (OLS contrast)', 'fixed-effect contrast (units of y)')
    R['bh_family'] = R.family + ' (BH within family, per model and contrast)'
    ALL = R.copy()
    ALL['model'] = pd.Categorical(ALL.model, [M_NEST, M_PAT, M_BIO, M_SITE, M_NOLIB, M_DIAG], ordered=True)
    ALL['contrast'] = pd.Categorical(ALL.contrast, list(CONTRASTS), ordered=True)
    fam_order = {FAM_6: 0, FAM_TF: 1, FAM_INF: 2}; feat_order = {(f[2], f[4]): i for i, f in enumerate(FEATURES)}
    ALL['_f'] = [feat_order[(a, b)] for a, b in zip(ALL.feature_key, ALL.compartment)]
    ALL = ALL.sort_values(['_f', 'contrast', 'model']).drop(columns='_f').reset_index(drop=True)
    COLS = ['family', 'figure_panel', 'feature', 'compartment', 'feature_type', 'model', 'contrast', 'estimate', 'se', 'df', 'statistic', 'p', 'q_bh',
            'p_wald_z', 'q_bh_wald_z', 'df_satterthwaite', 'p_satterthwaite', 'q_satterthwaite', 'test', 'estimate_definition', 'bh_family',
            'n_patients', 'n_samples', 'n_biopsies', 'n_cells', 'n_cells_rem_pre', 'n_cells_rem_post', 'n_cells_nonrem_pre', 'n_cells_nonrem_post',
            'n_samples_rem_pre', 'n_samples_rem_post', 'n_samples_nonrem_pre', 'n_samples_nonrem_post',
            'pct_cells_expressing_rem_pre', 'pct_cells_expressing_rem_post', 'pct_cells_expressing_nonrem_pre', 'pct_cells_expressing_nonrem_post',
            'var_patient', 'var_sample', 'var_biopsy', 'var_residual', 'libsize_coef', 'optimizer', 'fit_warnings', 'fit_flag', 'r_singular', 'feature_key']
    ALL = ALL[COLS]
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')
    ALL.apply(lambda c: c.map(sig3)).to_csv(f'{OUT}/cell_level_lmm_results.csv', index=False)
    log('wrote cell_level_lmm_results.csv', ALL.shape)


# ======================================================================================================================
# Section suppression_score_patient_level: Treg Suppression Score in TAURUS Tregs, per-patient tests and cell-level model
# Panels: Fig 6A inputs (Treg Suppression Score per patient; Supp Table 11)
# ======================================================================================================================
if _run('suppression_score_patient_level'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/40_taurus_treg_suppression_score.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF CD scRNA-seq: does the Treg Suppression Score change pre -> post by remission, and does it move with the
Treg pseudotime and the Treg activation programs? Numbers only, no figure.

Score: gene sets = the early and late Treg signatures of the CD scRNA-seq Treg pseudotime (bin 4 vs bin 0;
  panels/fig4_and_related/D_shared_inputs/treg_D6_signature_genes.csv, cosmx_cd/10_treg_suppression_score.py). Primary: the biopsy
  (CosMx WTX) gene sets, 50 early + 50 late; genes present in TAURUS kept.
  Per TAURUS Treg (aucell_cell_type_short == 'Treg', all 10,246) on X = log1p(CP10k): scanpy score_genes(late set) and
  score_genes(early set) (random_state 0, all TAURUS genes as control pool, scored over all Tregs).
  Treg Suppression Score = early-set score - late-set score (high = Activated, low = Low-suppressive Tregs).
  Not z-scored within sample: a within-sample z sets every sample mean to 0 and would erase the pre/post difference tested here. The
  per-cell score is used in its raw units; the LMM standardises it over all Tregs (estimates back-transformed).
  Sensitivity scores (patient-level tests only): (s1) mean z of early genes - mean z of late genes, z per gene over all TAURUS Tregs
  (cohort-wide); (s2) score_genes with the CD resection (CosMx 6K) gene sets (50 early, 35 late).
Tests
  1 per-patient mean score, pre vs post: Wilcoxon signed-rank within remitters (10 pairs) and non-remitters (5 pairs); Mann-Whitney U of
    the change (remission vs non-remission). Nested cell-level LMM (no covariate): y ~ post * remission + (1 | patient) + (1 | sample),
    statsmodels MixedLM REML, all 31 samples / 16 patients, t with between-within df 31 - 16 - 2 = 13; same model + libsize_z as
    sensitivity.
  2 Spearman of change in score vs change in inflammation (TAURUS Inflammation_score, Treg-weighted sample mean), all 15 paired and
    the 10 remitters.
  3 Spearman of change in score vs change in pseudotime (fig6B_patient_pseudotime.csv), vs change in IL-2/STAT5 (Hallmark ULM,
    per-patient median) and vs change in the co-stimulatory/suppressive receptor module (16 receptors, score_genes, per-patient mean);
    within-sample Spearman of score vs pseudotime (samples with >= 20 Tregs) for direction.
  4 baseline (pre) score, remitters (10) vs non-remitters (6): Mann-Whitney U.
Output panels/fig6c_taurus/treg_suppression_score/taurus_treg_suppression_score.csv (3 significant figures).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, decoupler as dc
    from scipy import stats
    from scipy.stats import spearmanr, wilcoxon, mannwhitneyu
    import statsmodels.formula.api as smf

    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
    P = f'{D}/panels/fig6c_taurus'; OUT = f'{P}/treg_suppression_score'; os.makedirs(OUT, exist_ok=True)
    KEYS = ['Patient', 'Remission_status', 'Treatment']
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')

    # ================================================================ Tregs and per-cell features
    ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); ct = ad.obs['aucell_cell_type_short'].astype(str).values
    treg = ad[ct == 'Treg'].copy(); del ad; treg.layers['lognorm'] = treg.X.copy(); log('Tregs', treg.n_obs)
    assert np.allclose(np.expm1(treg.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)                      # X = log1p(CP10k)
    F = treg.obs[KEYS + ['sample_id']].astype(str).copy(); F['infl'] = treg.obs['Inflammation_score'].astype(float).values
    F['total_counts'] = treg.obs['total_counts'].astype(float).values
    pc = pd.read_csv(f'{P}/fig6B_treg_pseudotime_per_cell.csv', index_col=0)
    assert set(pc.index) == set(F.index); F['pseudotime'] = pc.loc[F.index, 'dpt_pseudotime_norm'].values

    SG = pd.read_csv(f'{D}/panels/fig4_and_related/D_shared_inputs/treg_D6_signature_genes.csv')
    GS, NG = {}, {}
    for coh in ['biopsy', 'resection']:
        s_ = SG[SG.cohort == coh]
        GS[coh] = {side: [g for g in s_.loc[s_.set.str.startswith(side), 'gene'] if g in treg.var_names] for side in ['early', 'late']}
        NG[coh] = {side: (len(GS[coh][side]), int(s_.set.str.startswith(side).sum())) for side in ['early', 'late']}
    log('genes present (TAURUS / signature):', NG)
    for coh in ['biopsy', 'resection']:
        sc.tl.score_genes(treg, GS[coh]['late'], score_name='l', random_state=0)
        sc.tl.score_genes(treg, GS[coh]['early'], score_name='e', random_state=0)
        F[f'supp_{coh}'] = (treg.obs.e - treg.obs.l).values
    Xz = treg[:, GS['biopsy']['early'] + GS['biopsy']['late']].X; Xz = Xz.toarray() if hasattr(Xz, 'toarray') else np.asarray(Xz)
    Z = (Xz - Xz.mean(0)) / Xz.std(0, ddof=1); ne = len(GS['biopsy']['early'])
    F['supp_meanz'] = Z[:, :ne].mean(1) - Z[:, ne:].mean(1)
    # activation programs
    msig = pd.read_csv(f'{D}/archive/csv/msigdb_human_snapshot.csv.gz', usecols=['geneset', 'genesymbol'])
    SETS = ['HALLMARK_IL2_STAT5_SIGNALING', 'REACTOME_RUNX1_AND_FOXP3_CONTROL_THE_DEVELOPMENT_OF_REGULATORY_T_LYMPHOCYTES_TREGS',
            'GOBP_T_CELL_TOLERANCE_INDUCTION', 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION']
    net6 = msig[msig.geneset.isin(SETS)].drop_duplicates().rename(columns={'geneset': 'source', 'genesymbol': 'target'})
    dc.mt.ulm(data=treg, net=net6, layer='lognorm'); F['il2_stat5'] = np.asarray(dc.pp.get_obsm(treg, key='score_ulm')[:, 'HALLMARK_IL2_STAT5_SIGNALING'].X).ravel()
    cc = pd.read_csv(f'{D}/objects/cellchatdb_human_interactions_snapshot.csv')
    gs = lambda s: [g.strip() for g in str(s).replace(';', ',').split(',') if g.strip() and g.strip() != 'nan']
    det = pd.Series(np.asarray((treg.layers['lognorm'] > 0).mean(0)).ravel(), index=treg.var_names)
    pws = ['CD80', 'CD86', 'PD-L1', 'PDL2', 'ICOS', 'OX40', 'GITRL', 'CD70', 'CD40', 'IL2', 'IL10']
    recs = sorted({g for r in cc.loc[cc.pathway_name.isin(pws), 'receptor.symbol'] for g in gs(r)} & set(treg.var_names)); recs = [g for g in recs if det[g] >= 0.05]
    assert len(recs) == 16
    sc.tl.score_genes(treg, recs, score_name='costim_module'); F['costim_module'] = treg.obs['costim_module'].values
    del treg; log('per-cell features done')

    # ================================================================ per-sample / per-patient values
    SCORES = {'supp_biopsy': 'primary: score_genes, biopsy WTX gene sets', 'supp_meanz': 'sensitivity: mean z (cohort-wide), biopsy WTX gene sets',
              'supp_resection': 'sensitivity: score_genes, resection 6K gene sets'}
    S = F.groupby(KEYS).agg(**{k: (k, 'mean') for k in SCORES}, pt_mean=('pseudotime', 'mean'), il2_stat5=('il2_stat5', 'median'),
                            costim_module=('costim_module', 'mean'), inflammation=('infl', 'mean'), n_treg=('infl', 'size')).reset_index()
    assert len(S) == 31
    pv = lambda col: S.pivot_table(index=['Patient', 'Remission_status'], columns='Treatment', values=col)
    DL = pd.DataFrame({c: pv(c).Post - pv(c).Pre for c in list(SCORES) + ['pt_mean', 'il2_stat5', 'costim_module', 'inflammation']}).dropna().reset_index()
    assert len(DL) == 15
    PP = pd.read_csv(f'{P}/fig6B_patient_pseudotime.csv').set_index('Patient'); DL['pt_linechart'] = DL.Patient.map(PP.delta)
    PRE = S[S.Treatment == 'Pre']; assert len(PRE) == 16
    rem, non = DL.Remission_status == 'Remission', DL.Remission_status == 'Non_Remission'

    rows = []
    def add(score, section, test, group, n, estimate, estimate_def, statistic, p, **kw):
        rows.append(dict(score=SCORES.get(score, score), section=section, test=test, group=group, n=n, estimate=estimate,
                         estimate_definition=estimate_def, statistic=statistic, p=p, **kw))

    for sc_ in SCORES:
        pvs = pv(sc_)
        for gname, m in [('Remission', rem), ('Non_Remission', non)]:
            d_ = DL[m]; pp_ = pvs.loc[d_.set_index(['Patient', 'Remission_status']).index]; w = wilcoxon(pp_.Post, pp_.Pre)
            add(sc_, '1 pre vs post', 'Wilcoxon signed-rank (per-patient mean)', gname, len(d_), d_[sc_].median(), 'median change post - pre (score units)',
                w.statistic, w.pvalue, n_patients_down=int((d_[sc_] < 0).sum()), pre_median=pp_.Pre.median(), post_median=pp_.Post.median())
        u = mannwhitneyu(DL[rem][sc_], DL[non][sc_])
        add(sc_, '1 pre vs post', 'Mann-Whitney U on change', 'Remission vs Non_Remission', f'{rem.sum()} vs {non.sum()}',
            DL[rem][sc_].median() - DL[non][sc_].median(), 'difference of median changes (remission - non-remission)', u.statistic, u.pvalue)
        for gname, m in [('All paired', rem | non), ('Remission', rem)]:
            r = spearmanr(DL[m].inflammation, DL[m][sc_])
            add(sc_, '2 inflammation coupling', 'Spearman, change in score vs change in inflammation', gname, int(m.sum()), r.statistic, 'rho', np.nan, r.pvalue)
        for col, lab in [('pt_linechart', 'change in Treg pseudotime (line chart)'), ('il2_stat5', 'change in IL-2/STAT5'),
                         ('costim_module', 'change in co-stimulatory/suppressive receptor module')]:
            r = spearmanr(DL[col], DL[sc_])
            add(sc_, '3 agreement', f'Spearman, change in score vs {lab}', 'All paired', 15, r.statistic, 'rho', np.nan, r.pvalue)
        u = mannwhitneyu(PRE[PRE.Remission_status == 'Remission'][sc_], PRE[PRE.Remission_status == 'Non_Remission'][sc_])
        add(sc_, '4 baseline', 'Mann-Whitney U, pre-treatment per-patient mean', 'Remission vs Non_Remission', '10 vs 6',
            PRE[PRE.Remission_status == 'Remission'][sc_].median() - PRE[PRE.Remission_status == 'Non_Remission'][sc_].median(),
            'difference of medians (remission - non-remission)', u.statistic, u.pvalue)
    # direction: within-sample score vs pseudotime (primary score)
    rh = np.array([spearmanr(g.pseudotime, g.supp_biopsy).statistic for _, g in F.groupby(['Patient', 'Treatment']) if len(g) >= 20])
    add('supp_biopsy', '3 agreement', 'within-sample Spearman, score vs pseudotime (samples >= 20 Tregs)', 'All samples', len(rh), np.median(rh),
        'median within-sample rho', int((rh < 0).sum()), stats.binomtest(int((rh < 0).sum()), len(rh)).pvalue, n_samples_negative=int((rh < 0).sum()))
    log('patient-level tests done')

    # ================================================================ nested cell-level LMM (no covariate; + libsize_z as sensitivity)
    METHODS = [['lbfgs'], ['powell', 'lbfgs'], ['nm', 'lbfgs'], ['bfgs'], ['cg']]
    CONTRASTS = {'Change in remission (post − pre)': {'post': 1.0, 'post:rem': 1.0}, 'Change in non-remission (post − pre)': {'post': 1.0},
                 'Differential change (remission − non-remission)': {'post:rem': 1.0}, 'Baseline (remission − non-remission, pre)': {'rem': 1.0}}
    d = F[['Patient', 'Remission_status', 'Treatment', 'total_counts', 'supp_biopsy']].rename(columns={'supp_biopsy': 'y'}).copy()
    d['post'] = (d.Treatment == 'Post').astype(float); d['rem'] = (d.Remission_status == 'Remission').astype(float); d['sample'] = d.Patient + '_' + d.Treatment
    ll = np.log(d.total_counts); d['lib_z'] = (ll - ll.mean()) / ll.std(); sy = d.y.std(); d['ys'] = (d.y - d.y.mean()) / sy; d = d.reset_index(drop=True)
    df_bw = d['sample'].nunique() - d.Patient.nunique() - 2; assert df_bw == 13
    for lib in [False, True]:
        md = smf.mixedlm('ys ~ post * rem' + (' + lib_z' if lib else ''), d, groups='Patient', re_formula='1', vc_formula={'sample': '0 + C(sample)'})
        fits = []
        for meth in METHODS:
            try: r = md.fit(reml=True, method=meth)
            except Exception: continue
            dg = np.diag(r.cov_params().loc[r.fe_params.index, r.fe_params.index].values)
            if r.converged and np.isfinite(r.llf) and np.all(np.isfinite(dg)) and np.all(dg > 0): fits.append(r)
        assert fits, 'no converged fit'
        r = fits[int(np.argmax([f.llf for f in fits]))]; fe = r.fe_params; V = r.cov_params().loc[fe.index, fe.index]
        for k, w in CONTRASTS.items():
            L = pd.Series(0.0, index=fe.index)
            for n_, v_ in w.items(): L[n_] = v_
            est, se = float(L @ fe), float(np.sqrt(L.values @ V.values @ L.values)); t_ = est / se
            add('supp_biopsy', '1 pre vs post' if 'Baseline' not in k else '4 baseline',
                'nested LMM y ~ post*remission' + (' + libsize_z' if lib else '') + ' + (1|patient) + (1|sample), REML, t 13 df' + (' [sensitivity]' if lib else ' [primary]'),
                k, f'{len(d)} cells / 31 samples / 16 patients', est * sy, 'contrast estimate (score units)', t_, 2 * stats.t.sf(abs(t_), df_bw),
                estimate_sd_units=est, se=se * sy, df=df_bw)
    log('LMM done')

    R = pd.DataFrame(rows)
    R['genes_present'] = R.score.map({v: ('early %d/%d, late %d/%d' % (*NG['resection']['early'], *NG['resection']['late']) if 'resection' in v else
                                          'early %d/%d, late %d/%d' % (*NG['biopsy']['early'], *NG['biopsy']['late'])) for v in SCORES.values()})
    R = R.apply(lambda c: c.map(sig3))
    with open(f'{OUT}/taurus_treg_suppression_score.csv', 'w') as fh: R.to_csv(fh, index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70); pd.set_option('display.max_rows', 200)
    print(R.drop(columns=['estimate_definition', 'genes_present']).to_string(index=False))
    print('per-patient changes (primary score):'); print(DL[['Patient', 'Remission_status', 'supp_biopsy', 'pt_linechart', 'il2_stat5', 'costim_module', 'inflammation']].round(3).to_string(index=False))
    log('done')


# ======================================================================================================================
# Section F6A_programs_biopsy: biopsy-level mixed models of the Treg programs, pseudotime and 28 Treg TFs before and after anti-TNF
# Panels: Fig 6A receptor module (Supp Table 11)
# ======================================================================================================================
if _run('F6A_programs_biopsy'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/54_biopsy_level_prepost.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: pre/post Treg comparison with the BIOPSY as the unit. No covariates (no library size,
inflammation or dissociation).
  Pseudobulk: per biopsy (sample_id), the mean of each per-cell feature over its Tregs. Biopsies with < 10 Tregs are dropped
    (count reported).
  Models (column 'model'): every biopsy weighted equally (unweighted lmer); R lme4 + lmerTest, REML, Satterthwaite df; y standardised
    before fitting (estimates back-transformed).
    'patient only':                biopsy mean ~ post * remission + (1 | patient)
    'patient + patient:timepoint': biopsy mean ~ post * remission + (1 | patient) + (1 | patient:timepoint)
      (the patient x timepoint level = the 'sample' of cell_level_lmm: biopsies of one patient at one timepoint are not independent
       replicates of the pre/post contrast).
  Contrasts: change in remission = b_post + b_post:rem; change in non-remission = b_post; difference (remission - non-remission) = b_post:rem;
    'patient + patient:timepoint' model only: baseline difference (remission - non-remission, pre) = b_rem.
  95% CI: estimate +/- t(0.975, Satterthwaite df) x SE.  BH within family (TFs 28; programs 4 + pseudotime 5), per model and contrast.
Features (per-cell values from supplementary/cache/Treg.parquet, written by cell_level_lmm):
  TF family (28): pre-specified Treg TFs, decoupler ULM on CollecTRI (tmin 5).
  Program family (5): Hallmark IL-2/STAT5 (ULM), co-stimulatory/suppressive receptor module (16 CellChatDB receptors, score_genes),
    TNFRSF9 (log-normalised), TNFR2 non-canonical NF-kB (Reactome without PSM genes, ULM), and Treg diffusion pseudotime.
Side by side: cell-level model F of cell_level_lmm (y ~ post * remission + (1 | patient) + (1 | sample), no covariates; between-within
  df): estimate, P and BH q from supplementary/cell_level_lmm_results.csv.
Output: supplementary/biopsy_level_prepost.csv (3 significant figures; Supp Table 11). No figure.
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time, subprocess, tempfile, shutil
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd
    from statsmodels.stats.multitest import multipletests

    HERE = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary'
    CACHE = f'{HERE}/cache/Treg.parquet'; CELL = f'{HERE}/cell_level_lmm_results.csv'; OUT = f'{HERE}/biopsy_level_prepost.csv'
    MIN_TREG = 10
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')

    # ================================================================ features (cell_level_lmm definitions)
    FAM_P, FAM_TF = 'Programs (4 + pseudotime)', 'Treg TFs (28)'
    PROG = [('il2_stat5', 'IL-2/STAT5 Signaling (Hallmark, ULM)'), ('costim_module', 'Co-Stimulatory/Suppressive Receptor Module'),
            ('gene_TNFRSF9', 'TNFRSF9 (4-1BB)'), ('tnfr2_pathway', 'TNFR2 Non-Canonical NF-κB (Reactome, ULM)'), ('pseudotime', 'Treg Pseudotime')]
    F = pd.read_parquet(CACHE)
    TFS = [c[3:] for c in F.columns if c.startswith('tf_')]; assert len(TFS) == 28, TFS
    FEATS = [(FAM_P, k, l) for k, l in PROG] + [(FAM_TF, f'tf_{t}', f'{t} activity (CollecTRI, ULM)') for t in TFS]
    FCOLS = [k for _, k, _ in FEATS]; assert set(FCOLS) <= set(F.columns) and F[FCOLS].notna().all().all()

    # ================================================================ biopsy pseudobulk
    F['biopsy'] = F.sample_id.astype(str)
    assert F.groupby('biopsy').Patient.nunique().max() == 1 and F.groupby('biopsy').Treatment.nunique().max() == 1
    F['post'] = (F.Treatment == 'Post').astype(float); F['rem'] = (F.Remission_status == 'Remission').astype(float)
    assert set(F.Treatment) == {'Pre', 'Post'} and set(F.Remission_status) == {'Remission', 'Non_Remission'}
    PB = (F.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                  post=('post', 'first'), rem=('rem', 'first'), n_treg=('post', 'size'))
          .join(F.groupby('biopsy')[FCOLS].mean()).reset_index())
    x = F[F.biopsy == PB.biopsy.iloc[0]]; assert np.allclose(PB.loc[0, FCOLS].astype(float).values, x[FCOLS].mean().values)
    assert PB.n_treg.sum() == len(F)
    DROP = PB[PB.n_treg < MIN_TREG]; PB = PB[PB.n_treg >= MIN_TREG].reset_index(drop=True)
    log(f'{len(F)} Tregs in {len(PB) + len(DROP)} biopsies from {F.Patient.nunique()} patients; dropped {len(DROP)} biopsies with < {MIN_TREG} Tregs '
    f'({int(DROP.n_treg.sum())} Tregs): ' + ', '.join(f'{r.Patient} {r.Treatment} ({r.n_treg})' for r in DROP.itertuples()))
    GRP = [('rem_pre', 1.0, 0.0), ('rem_post', 1.0, 1.0), ('nonrem_pre', 0.0, 0.0), ('nonrem_post', 0.0, 1.0)]
    NN = {}
    for nm, r_, p_ in GRP:
        s = PB[(PB.rem == r_) & (PB.post == p_)]
        NN.update({f'n_biopsies_{nm}': len(s), f'n_patients_{nm}': s.Patient.nunique(), f'n_tregs_{nm}': int(s.n_treg.sum())})
    pairs = PB.groupby(['Patient', 'rem']).post.nunique().reset_index()
    NN['n_paired_patients_rem'] = int(((pairs.post == 2) & (pairs.rem == 1)).sum()); NN['n_paired_patients_nonrem'] = int(((pairs.post == 2) & (pairs.rem == 0)).sum())
    NN['n_biopsies_dropped'] = len(DROP)
    log('after drop: ' + ', '.join(f'{k}={v}' for k, v in NN.items()))

    # ================================================================ lme4 / lmerTest, Satterthwaite
    RCODE = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); d <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); cols <- readLines(a[2]); out <- list()
d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
FORMS <- list(patient_only = ys ~ post * rem + (1 | Patient), patient_timepoint = ys ~ post * rem + (1 | Patient) + (1 | pt))
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1))
for (mk in names(FORMS)) for (cc in cols) {
  y <- d[[cc]]; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  m <- lmer(FORMS[[mk]], data = d, REML = TRUE)
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m))
  for (k in c(names(Ls), if (mk == "patient_timepoint") "baseline_difference")) {
    Lk <- if (k == "baseline_difference") c(rem = 1) else Ls[[k]]
    L <- setNames(rep(0, length(fe)), fe); L[names(Lk)] <- Lk
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
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        PB.assign(pt=PB.Patient + '_' + PB.Treatment)[['Patient', 'pt', 'post', 'rem'] + FCOLS].to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/cols.txt', 'w') as fh: fh.write('\n'.join(FCOLS) + '\n')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        log(f'fitting {len(FCOLS)} features x 2 models in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/cols.txt', f'{tmp}/r_out.csv'], check=True)
        R = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(R) == (3 + 4) * len(FCOLS) and R.p.notna().all()
    MODS = {'patient_only': 'patient only', 'patient_timepoint': 'patient + patient:timepoint'}
    FORM = {'patient only': 'biopsy mean ~ post × remission + (1 | patient); lme4/lmerTest REML, Satterthwaite df',
            'patient + patient:timepoint': 'biopsy mean ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df'}
    R['model'] = R.model_key.map(MODS); R['model_formula'] = R.model.map(FORM)
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'baseline_difference': 'Baseline difference (remission − non-remission, pre)'}
    R['contrast'] = R.contrast_key.map(CON)
    from scipy import stats
    tq = stats.t.ppf(0.975, R.df); R['ci_low'] = R.estimate - tq * R.se; R['ci_high'] = R.estimate + tq * R.se
    R = R.merge(pd.DataFrame(FEATS, columns=['family', 'feature_key', 'feature']), on='feature_key')
    R['q_bh'] = np.nan
    for _, i in R.groupby(['family', 'model', 'contrast']).groups.items():
        R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
    R['bh_family'] = R.family + ', BH within family, per model and contrast'

    # ================================================================ cell-level model F (cell_level_lmm)
    C = pd.read_csv(CELL); C = C[C.model.str.startswith('F nested LMM without libsize_z')]
    C = C[['feature_key', 'contrast', 'estimate', 'p', 'q_bh']].rename(columns={'estimate': 'cell_F_estimate', 'p': 'cell_F_p', 'q_bh': 'cell_F_q_bh'})
    R = R.merge(C, on=['feature_key', 'contrast'], how='left'); assert (R.cell_F_q_bh.notna() == (R.contrast_key != 'baseline_difference')).all()
    R['sign_agrees_with_cell_F'] = (np.sign(R.estimate) == np.sign(R.cell_F_estimate)).where(R.cell_F_estimate.notna())
    R['sig_biopsy'] = R.q_bh < 0.05; R['sig_cell_F'] = R.cell_F_q_bh < 0.05
    R['min_tregs_per_biopsy'] = MIN_TREG
    for k, v in NN.items(): R[k] = v
    R['n_patients'] = PB.Patient.nunique(); R['n_biopsies'] = len(PB); R['n_tregs'] = int(PB.n_treg.sum())

    COLS = (['family', 'feature', 'contrast', 'model', 'model_formula', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'bh_family',
             'cell_F_estimate', 'cell_F_p', 'cell_F_q_bh', 'sign_agrees_with_cell_F', 'sig_biopsy', 'sig_cell_F', 'var_patient', 'var_patient_timepoint', 'var_residual',
             'singular', 'converged', 'min_tregs_per_biopsy', 'n_patients', 'n_biopsies', 'n_tregs'] + list(NN) + ['feature_key'])
    R['_f'] = R.feature_key.map({k: i for i, k in enumerate(FCOLS)}); R['_c'] = R.contrast.map({v: i for i, v in enumerate(CON.values())})
    R['_m'] = R.model.map({v: i for i, v in enumerate(MODS.values())})
    R = R.sort_values(['_m', '_c', '_f'])[COLS].reset_index(drop=True)
    with open(OUT, 'w') as fh:
        R.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, R.shape)

    show = R.copy(); show['feature'] = show.feature.str[:28]; show['contrast'] = show.contrast.str[:18]
    show = show[show.model == 'patient + patient:timepoint']
    print(show[['contrast', 'feature', 'estimate', 'ci_low', 'ci_high', 'df', 'p', 'q_bh', 'cell_F_estimate', 'cell_F_q_bh', 'singular']]
          .apply(lambda c: c.map(sig3)).to_string(index=False))
    for mod in MODS.values():
      M = R[R.model == mod]
      print(f'\n== {mod}: singular fits: {sorted(M[M.singular].feature.str.split(" ").str[0].unique())}; non-converged: {int((~M.converged).sum())}; '
        f'Satterthwaite df range {M.df.min():.3g}-{M.df.max():.3g} (by contrast: ' +
            ', '.join(f"{k.split(' (')[0]} {g.df.min():.3g}-{g.df.max():.3g}" for k, g in M.groupby('contrast', sort=False)) + ')')
      for k in CON.values():
        s = M[M.contrast == k]
        for fam in [FAM_P, FAM_TF]:
            t = s[s.family == fam]
            if t.empty: continue
            print(f'{k} | {fam}: biopsy P<0.05: {[(f.split(" ")[0], sig3(e), sig3(q)) for f, e, q in t[t.p < 0.05][["feature", "estimate", "q_bh"]].values]} | biopsy q<0.05: {t[t.sig_biopsy].feature.str.split(" ").str[0].tolist()} | cell F q<0.05: '
              f'{t[t.sig_cell_F].feature.str.split(" ").str[0].tolist()} | sign disagreements: {int((t.sign_agrees_with_cell_F == False).sum())}')


# ======================================================================================================================
# Section F6A_suppression_score_biopsy: biopsy-level mixed models of the Treg Suppression Score and suppressive effector genes
# Panels: Fig 6A Treg Suppression Score (Supp Table 11)
# ======================================================================================================================
if _run('F6A_suppression_score_biopsy'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/55_treg_suppression_biopsy_level.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: Treg suppressive capacity in remitters vs non-remitters, at baseline and in its pre/post change,
with the BIOPSY as the unit. Numbers only, no figure.

Model (as F6A_programs_biopsy, 'patient + patient:timepoint' model only; no covariates):
  Pseudobulk: per biopsy (sample_id), the mean of each per-cell feature over its Tregs; biopsies with < 10 Tregs dropped.
  biopsy mean ~ post * remission + (1 | patient) + (1 | patient:timepoint); R lme4 + lmerTest, REML, Satterthwaite df; every biopsy weighted
  equally; y standardised before fitting (estimates back-transformed to feature units).
  Contrasts: change in remission = b_post + b_post:rem; change in non-remission = b_post; differential change (remission - non-remission)
  = b_post:rem; baseline difference (remission - non-remission, pre) = b_rem.  95% CI = estimate +/- t(0.975, Satterthwaite df) x SE.
  BH: family (a) the Treg Suppression Score alone (1 test, q = P); family (b) the effector genes present; each per contrast.
Features (per TAURUS Treg: aucell_cell_type_short == 'Treg', 10,246 Tregs; X = log1p(CP10k)):
  (a) Treg Suppression Score, primary version of suppression_score_patient_level: scanpy score_genes(late set) and score_genes(early set)
      with the biopsy (CosMx WTX) early/late gene sets of the CD scRNA-seq Treg pseudotime signature (genes present in TAURUS),
      random_state 0, all TAURUS genes as control pool, scored over all Tregs; score = early - late (high = Activated,
      low = Low-suppressive). Raw units.
  (b) Suppressive effector genes (log-normalised expression), pre-specified canonical Treg suppression mechanisms:
      CTLA4, IL10, TGFB1, EBI3, LAG3, ENTPD1, NT5E, LRRC32, TIGIT, GZMB, IL2RA, HAVCR2. Absent genes skipped; % of Tregs expressing
      (count > 0) reported over all Tregs and over the Tregs of the analysed biopsies. IL2RA (CD25) also reflects IL-2 signalling.
Output: supplementary/treg_suppression_biopsy_level.csv (3 significant figures; Supp Table 11).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time, subprocess, tempfile, shutil, glob
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    HERE = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary'
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
    CACHE = f'{HERE}/cache/Treg.parquet'; OUT = f'{HERE}/treg_suppression_biopsy_level.csv'
    MIN_TREG = 10
    EFFECTORS = ['CTLA4', 'IL10', 'TGFB1', 'EBI3', 'LAG3', 'ENTPD1', 'NT5E', 'LRRC32', 'TIGIT', 'GZMB', 'IL2RA', 'HAVCR2']
    NOTE = {'IL2RA': 'IL2RA (CD25) is also an IL-2/STAT5 target and reflects IL-2 signalling, not only IL-2 consumption as a suppression mechanism'}
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')

    # ================================================================ Tregs and per-cell features
    ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); ct = ad.obs['aucell_cell_type_short'].astype(str).values
    treg = ad[ct == 'Treg'].copy(); del ad; treg.layers['lognorm'] = treg.X.copy(); log('Tregs', treg.n_obs)
    assert np.allclose(np.expm1(treg.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)                      # X = log1p(CP10k)
    SGF = glob.glob(f'{D}/panels/fig4_and_related/D_shared_inputs/treg_*_signature_genes.csv'); assert len(SGF) == 1, SGF
    SG = pd.read_csv(SGF[0]); s_ = SG[SG.cohort == 'biopsy']
    GS = {side: [g for g in s_.loc[s_.set.str.startswith(side), 'gene'] if g in treg.var_names] for side in ['early', 'late']}
    NG = {side: (len(GS[side]), int(s_.set.str.startswith(side).sum())) for side in ['early', 'late']}
    assert NG == {'early': (49, 50), 'late': (50, 50)}, NG                                              # genes present in TAURUS
    sc.tl.score_genes(treg, GS['late'], score_name='l', random_state=0)
    sc.tl.score_genes(treg, GS['early'], score_name='e', random_state=0)
    F = pd.read_parquet(CACHE)[['Patient', 'Remission_status', 'Treatment', 'sample_id', 'pseudotime']]
    assert set(F.index) == set(treg.obs_names) and len(F) == treg.n_obs
    F = F.loc[treg.obs_names].copy()
    assert (F.Patient.astype(str).values == treg.obs.Patient.astype(str).values).all() and (F.Treatment.astype(str).values == treg.obs.Treatment.astype(str).values).all()
    F['supp_score'] = (treg.obs.e - treg.obs.l).values
    GENES = [g for g in EFFECTORS if g in treg.var_names]; ABSENT = [g for g in EFFECTORS if g not in treg.var_names]
    Xg = treg[:, GENES].layers['lognorm']; Xg = Xg.toarray() if hasattr(Xg, 'toarray') else np.asarray(Xg)
    for i, g in enumerate(GENES): F[f'gene_{g}'] = Xg[:, i]
    del treg; log(f'effector genes present {len(GENES)}/{len(EFFECTORS)}; absent: {ABSENT}')

    FAM_S, FAM_G = 'Treg Suppression Score', f'Suppressive effector genes ({len(GENES)})'
    FEATS = [(FAM_S, 'supp_score', 'Treg Suppression Score (early − late, score_genes, biopsy WTX gene sets)')] + \
            [(FAM_G, f'gene_{g}', f'{g} (log-normalised expression)') for g in GENES]
    FCOLS = [k for _, k, _ in FEATS]; assert F[FCOLS].notna().all().all()

    # ================================================================ biopsy pseudobulk
    F['biopsy'] = F.sample_id.astype(str)
    assert F.groupby('biopsy').Patient.nunique().max() == 1 and F.groupby('biopsy').Treatment.nunique().max() == 1
    F['post'] = (F.Treatment == 'Post').astype(float); F['rem'] = (F.Remission_status == 'Remission').astype(float)
    assert set(F.Treatment) == {'Pre', 'Post'} and set(F.Remission_status) == {'Remission', 'Non_Remission'}
    PB = (F.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                  post=('post', 'first'), rem=('rem', 'first'), n_treg=('post', 'size'))
          .join(F.groupby('biopsy')[FCOLS].mean()).reset_index())
    x = F[F.biopsy == PB.biopsy.iloc[0]]; assert np.allclose(PB.loc[0, FCOLS].astype(float).values, x[FCOLS].mean().values)
    assert PB.n_treg.sum() == len(F)
    DROP = PB[PB.n_treg < MIN_TREG]; PB = PB[PB.n_treg >= MIN_TREG].reset_index(drop=True)
    log(f'{len(F)} Tregs in {len(PB) + len(DROP)} biopsies from {F.Patient.nunique()} patients; dropped {len(DROP)} biopsies with < {MIN_TREG} Tregs '
    f'({int(DROP.n_treg.sum())} Tregs)')
    GRP = [('rem_pre', 1.0, 0.0), ('rem_post', 1.0, 1.0), ('nonrem_pre', 0.0, 0.0), ('nonrem_post', 0.0, 1.0)]
    NN = {}
    for nm, r_, p_ in GRP:
        s = PB[(PB.rem == r_) & (PB.post == p_)]
        NN.update({f'n_biopsies_{nm}': len(s), f'n_patients_{nm}': s.Patient.nunique(), f'n_tregs_{nm}': int(s.n_treg.sum())})
    pairs = PB.groupby(['Patient', 'rem']).post.nunique().reset_index()
    NN['n_paired_patients_rem'] = int(((pairs.post == 2) & (pairs.rem == 1)).sum()); NN['n_paired_patients_nonrem'] = int(((pairs.post == 2) & (pairs.rem == 0)).sum())
    NN['n_biopsies_dropped'] = len(DROP)
    log('after drop: ' + ', '.join(f'{k}={v}' for k, v in NN.items()))
    KEPT = F.biopsy.isin(PB.biopsy)
    PCT = {f'gene_{g}': (100 * (F[f'gene_{g}'] > 0).mean(), 100 * (F.loc[KEPT, f'gene_{g}'] > 0).mean()) for g in GENES}
    # group means (biopsy-level, of biopsy means) for direction reading
    GM = {k: {nm: PB[(PB.rem == r_) & (PB.post == p_)][k].mean() for nm, r_, p_ in GRP} for k in FCOLS}

    # ================================================================ lme4 / lmerTest, Satterthwaite (patient + patient:timepoint model)
    RCODE = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); d <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); cols <- readLines(a[2]); out <- list()
d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
FORMS <- list(patient_timepoint = ys ~ post * rem + (1 | Patient) + (1 | pt))
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1),
           baseline_difference = c(rem = 1))
for (mk in names(FORMS)) for (cc in cols) {
  y <- d[[cc]]; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  m <- lmer(FORMS[[mk]], data = d, REML = TRUE)
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m))
  for (k in names(Ls)) {
    L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
    ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
    out[[length(out) + 1]] <- data.frame(model_key = mk, feature_key = cc, contrast_key = k, estimate = ct$Estimate * sy, se = ct$`Std. Error` * sy,
      df = ct$df, statistic = ct$`t value`, p = ct$`Pr(>|t|)`, var_patient = vc$vcov[vc$grp == "Patient"] * sy^2,
      var_patient_timepoint = vc$vcov[vc$grp == "pt"] * sy^2, var_residual = vc$vcov[vc$grp == "Residual"] * sy^2, singular = isSingular(m),
      converged = !any(grepl("converge|gradient|Hessian|eigenvalue", c(m@optinfo$conv$lme4$messages, m@optinfo$warnings))))
  }
}
write.csv(do.call(rbind, out), a[3], row.names = FALSE)
"""
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        PB.assign(pt=PB.Patient + '_' + PB.Treatment)[['Patient', 'pt', 'post', 'rem'] + FCOLS].to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/cols.txt', 'w') as fh: fh.write('\n'.join(FCOLS) + '\n')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        log(f'fitting {len(FCOLS)} features in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/cols.txt', f'{tmp}/r_out.csv'], check=True)
        R = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(R) == 4 * len(FCOLS) and R.p.notna().all()
    R['model'] = 'patient + patient:timepoint'
    R['model_formula'] = 'biopsy mean ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df; no covariates'
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'baseline_difference': 'Baseline difference (remission − non-remission, pre)'}
    R['contrast'] = R.contrast_key.map(CON)
    tq = stats.t.ppf(0.975, R.df); R['ci_low'] = R.estimate - tq * R.se; R['ci_high'] = R.estimate + tq * R.se
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
    R['_f'] = R.feature_key.map({k: i for i, k in enumerate(FCOLS)}); R['_c'] = R.contrast.map({v: i for i, v in enumerate(CON.values())})
    R = R.sort_values(['_c', '_f'])[COLS].reset_index(drop=True)
    with open(OUT, 'w') as fh:
        R.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, R.shape)

    pd.set_option('display.width', 250)
    show = R.copy(); show['feature'] = show.feature.str.split(' ').str[0]; show['contrast'] = show.contrast.str[:24]
    print(show[['contrast', 'feature', 'estimate', 'ci_low', 'ci_high', 'df', 'p', 'q_bh', 'pct_tregs_expressing_analysed_biopsies', 'singular']]
          .apply(lambda c: c.map(sig3)).to_string(index=False))
    print(f'singular fits: {sorted(R[R.singular].feature.str.split(" ").str[0].unique())}; non-converged: {int((~R.converged).sum())}')
    for k in CON.values():
        t = R[(R.contrast == k) & (R.family == FAM_G)]
        print(f'{k}: genes P<0.05: {[(f.split(" ")[0], sig3(e), sig3(p), sig3(q)) for f, e, p, q in t[t.p < 0.05][["feature", "estimate", "p", "q_bh"]].values]} | q<0.05: {t[t.q_bh < 0.05].feature.str.split(" ").str[0].tolist()}')


# ======================================================================================================================
# Section F6A_post_treatment_difference: post-treatment Remission vs Non-remission contrast added to the biopsy-level models
# Panels: Fig 6A (adds post-treatment difference rows to the F6A_programs_biopsy / F6A_suppression_score_biopsy CSVs; Supp Table 11)
# ======================================================================================================================
if _run('F6A_post_treatment_difference'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/60_post_treatment_difference.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: post-treatment difference between remitters and non-remitters, from the same biopsy-level model
as F6A_programs_biopsy / F6A_suppression_score_biopsy. Numbers only, no figure.

Model ('patient + patient:timepoint' model only; no covariates), fitted for the 6 Fig 6A features:
  Pseudobulk: per biopsy (sample_id), the mean of each per-cell feature over its Tregs; biopsies with < 10 Tregs dropped.
  biopsy mean ~ post * remission + (1 | patient) + (1 | patient:timepoint); R lme4 + lmerTest, REML, Satterthwaite df; every biopsy weighted
  equally; y standardised before fitting (estimates back-transformed to feature units).
  Post-treatment difference (remission - non-remission, post) = b_rem + b_post:rem (linear contrast, lmerTest contest, Satterthwaite df).
  95% CI = estimate +/- t(0.975, Satterthwaite df) x SE.
  The 4 other contrasts (change in remission, change in non-remission, differential change, baseline difference) are fitted too and
  must equal (3 significant figures) the rows written by the two upstream sections, so the new rows come from the same model.
  BH for the new rows: programs family (4 programs + pseudotime, 5 tests) within the new contrast; Treg Suppression Score alone
  (1 test, q = P).
Features: Treg Suppression Score (score_genes, biopsy WTX early / late gene sets, random_state 0, early - late); Treg Pseudotime,
  IL-2/STAT5 (Hallmark, ULM), TNFR2 non-canonical NF-kB (Reactome, ULM), TNFRSF9 (4-1BB), co-stimulatory/suppressive receptor module
  (supplementary/cache/Treg.parquet).
Output: the new contrast rows are appended to supplementary/biopsy_level_prepost.csv (5 rows) and
  supplementary/treg_suppression_biopsy_level.csv (1 row), contrast = 'Post-treatment difference (remission − non-remission, post)'.
  The cell-level model F columns are empty for the new rows (no cell-level counterpart), as for the baseline rows.
  Run it once after each run of the two upstream sections (it stops if the rows already exist).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time, subprocess, tempfile, shutil, glob
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    HERE = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary'
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
    CACHE = f'{HERE}/cache/Treg.parquet'; O54 = f'{HERE}/biopsy_level_prepost.csv'; O55 = f'{HERE}/treg_suppression_biopsy_level.csv'
    MIN_TREG = 10; MOD = 'patient + patient:timepoint'
    NEWC = 'Post-treatment difference (remission − non-remission, post)'
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'baseline_difference': 'Baseline difference (remission − non-remission, pre)', 'post_difference': NEWC}
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')

    A54, A55 = pd.read_csv(O54), pd.read_csv(O55)
    assert NEWC not in set(A54.contrast) and NEWC not in set(A55.contrast), 'post-treatment difference rows already present'
    PROGK = ['il2_stat5', 'costim_module', 'gene_TNFRSF9', 'tnfr2_pathway', 'pseudotime']          # programs family of F6A_programs_biopsy
    FCOLS = ['supp_score'] + PROGK

    # ================================================================ per-cell features (Treg Suppression Score as F6A_suppression_score_biopsy)
    ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); ct = ad.obs['aucell_cell_type_short'].astype(str).values
    treg = ad[ct == 'Treg'].copy(); del ad; log('Tregs', treg.n_obs)
    assert np.allclose(np.expm1(treg.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)                      # X = log1p(CP10k)
    SGF = glob.glob(f'{D}/panels/fig4_and_related/D_shared_inputs/treg_*_signature_genes.csv'); assert len(SGF) == 1, SGF
    SG = pd.read_csv(SGF[0]); s_ = SG[SG.cohort == 'biopsy']
    GS = {side: [g for g in s_.loc[s_.set.str.startswith(side), 'gene'] if g in treg.var_names] for side in ['early', 'late']}
    NG = {side: (len(GS[side]), int(s_.set.str.startswith(side).sum())) for side in ['early', 'late']}
    assert NG == {'early': (49, 50), 'late': (50, 50)}, NG                                              # genes present in TAURUS
    sc.tl.score_genes(treg, GS['late'], score_name='l', random_state=0)
    sc.tl.score_genes(treg, GS['early'], score_name='e', random_state=0)
    F = pd.read_parquet(CACHE)
    assert set(F.index) == set(treg.obs_names) and len(F) == treg.n_obs
    F = F.loc[treg.obs_names].copy()
    assert (F.Patient.astype(str).values == treg.obs.Patient.astype(str).values).all() and (F.Treatment.astype(str).values == treg.obs.Treatment.astype(str).values).all()
    F['supp_score'] = (treg.obs.e - treg.obs.l).values; del treg
    assert F[FCOLS].notna().all().all()

    # ================================================================ biopsy pseudobulk
    F['biopsy'] = F.sample_id.astype(str)
    assert F.groupby('biopsy').Patient.nunique().max() == 1 and F.groupby('biopsy').Treatment.nunique().max() == 1
    F['post'] = (F.Treatment == 'Post').astype(float); F['rem'] = (F.Remission_status == 'Remission').astype(float)
    assert set(F.Treatment) == {'Pre', 'Post'} and set(F.Remission_status) == {'Remission', 'Non_Remission'}
    PB = (F.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                  post=('post', 'first'), rem=('rem', 'first'), n_treg=('post', 'size'))
          .join(F.groupby('biopsy')[FCOLS].mean()).reset_index())
    assert PB.n_treg.sum() == len(F)
    DROP = PB[PB.n_treg < MIN_TREG]; PB = PB[PB.n_treg >= MIN_TREG].reset_index(drop=True)
    ref = A54[A54.model == MOD].iloc[0]
    assert len(PB) == ref.n_biopsies and int(PB.n_treg.sum()) == ref.n_tregs and PB.Patient.nunique() == ref.n_patients and len(DROP) == ref.n_biopsies_dropped
    log(f'{len(PB)} biopsies, {int(PB.n_treg.sum())} Tregs, {PB.Patient.nunique()} patients; dropped {len(DROP)}')

    # ================================================================ lme4 / lmerTest, Satterthwaite (+ the post-treatment contrast)
    RCODE = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); d <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); cols <- readLines(a[2]); out <- list()
d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
FORMS <- list(patient_timepoint = ys ~ post * rem + (1 | Patient) + (1 | pt))
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1),
           baseline_difference = c(rem = 1), post_difference = c(rem = 1, `post:rem` = 1))
for (mk in names(FORMS)) for (cc in cols) {
  y <- d[[cc]]; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  m <- lmer(FORMS[[mk]], data = d, REML = TRUE)
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m))
  for (k in names(Ls)) {
    L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
    ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
    out[[length(out) + 1]] <- data.frame(model_key = mk, feature_key = cc, contrast_key = k, estimate = ct$Estimate * sy, se = ct$`Std. Error` * sy,
      df = ct$df, statistic = ct$`t value`, p = ct$`Pr(>|t|)`, var_patient = vc$vcov[vc$grp == "Patient"] * sy^2,
      var_patient_timepoint = vc$vcov[vc$grp == "pt"] * sy^2, var_residual = vc$vcov[vc$grp == "Residual"] * sy^2, singular = isSingular(m),
      converged = !any(grepl("converge|gradient|Hessian|eigenvalue", c(m@optinfo$conv$lme4$messages, m@optinfo$warnings))))
  }
}
write.csv(do.call(rbind, out), a[3], row.names = FALSE)
"""
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        PB.assign(pt=PB.Patient + '_' + PB.Treatment)[['Patient', 'pt', 'post', 'rem'] + FCOLS].to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/cols.txt', 'w') as fh: fh.write('\n'.join(FCOLS) + '\n')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        log(f'fitting {len(FCOLS)} features in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/cols.txt', f'{tmp}/r_out.csv'], check=True)
        R = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(R) == 5 * len(FCOLS) and R.p.notna().all()
    R['contrast'] = R.contrast_key.map(CON); R['model'] = MOD
    tq = stats.t.ppf(0.975, R.df); R['ci_low'] = R.estimate - tq * R.se; R['ci_high'] = R.estimate + tq * R.se
    R['q_bh'] = np.nan
    for k in CON:   # BH: programs family of 5; Treg Suppression Score alone
        i = R.index[(R.contrast_key == k) & R.feature_key.isin(PROGK)]; R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
        i = R.index[(R.contrast_key == k) & (R.feature_key == 'supp_score')]; R.loc[i, 'q_bh'] = R.loc[i, 'p']

    # ================================================================ the 4 other contrasts must equal the upstream rows (same model and biopsies)
    CHK = ['estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'var_patient', 'var_residual']
    up = pd.concat([A54[A54.model == MOD], A55[A55.model == MOD]]); up = up[up.feature_key.isin(FCOLS)].set_index(['feature_key', 'contrast'])
    fit = R[R.contrast_key != 'post_difference'].set_index(['feature_key', 'contrast']).reindex(up.index)
    assert len(up) == 4 * len(FCOLS) and fit.estimate.notna().all()
    for c in CHK:
        a_, b_ = fit[c].map(sig3).astype(float).values, up[c].astype(float).values
        assert np.allclose(a_, b_, rtol=2e-3, atol=1e-12), (c, pd.DataFrame({'fit': a_, 'upstream': b_}, index=up.index)[~np.isclose(a_, b_, rtol=2e-3)])
    assert (fit.singular.values == up.singular.astype(bool).values).all() and (fit.converged.values == up.converged.astype(bool).values).all()
    log(f'fit matches all {len(up)} upstream rows of the 6 features (estimate, CI, SE, df, t, P, q, variance components; 3 significant figures)')

    # ================================================================ append the new contrast rows (other columns copied from the feature's baseline row)
    NR = R[R.contrast_key == 'post_difference'].set_index('feature_key')
    STATC = ['estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'var_patient', 'var_patient_timepoint', 'var_residual', 'singular', 'converged']
    def new_rows(A, keys):
        base = A[(A.model == MOD) & (A.contrast == CON['baseline_difference']) & A.feature_key.isin(keys)].copy()
        assert len(base) == len(keys)
        base['contrast'] = NEWC
        for c in STATC: base[c] = base.feature_key.map(NR[c]).map(sig3)
        if 'sig_biopsy' in base: base['sig_biopsy'] = base.q_bh < 0.05
        for c in ['cell_F_estimate', 'cell_F_p', 'cell_F_q_bh', 'sign_agrees_with_cell_F']:
            if c in base: base[c] = np.nan
        if 'sig_cell_F' in base: base['sig_cell_F'] = False
        return base
    N54, N55 = new_rows(A54, PROGK), new_rows(A55, ['supp_score'])
    with open(O54, 'w') as fh: pd.concat([A54, N54], ignore_index=True).to_csv(fh, index=False)
    with open(O55, 'w') as fh: pd.concat([A55, N55], ignore_index=True).to_csv(fh, index=False)
    C54, C55 = pd.read_csv(O54), pd.read_csv(O55)
    assert len(C54) == len(A54) + 5 and len(C55) == len(A55) + 1
    pd.testing.assert_frame_equal(C54.iloc[:len(A54)], A54); pd.testing.assert_frame_equal(C55.iloc[:len(A55)], A55)
    log(f'appended {len(N54)} rows to {os.path.basename(O54)} and {len(N55)} row to {os.path.basename(O55)}')

    pd.set_option('display.width', 250)
    show = R[['feature_key', 'contrast_key', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'p', 'q_bh', 'singular']].copy()
    print(show.apply(lambda c: c.map(sig3)).to_string(index=False))


# ======================================================================================================================
# Section F6A_pathways_biopsy: the six Treg suppression pathways of Fig 4G scored in TAURUS Tregs, biopsy-level model
# Panels: Fig 6A IL-2/STAT5 and TNFR2 (Supp Table 11)
# ======================================================================================================================
if _run('F6A_pathways_biopsy'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/65_fig3G_pathways_taurus.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: the six 'Treg suppression' pathways of the CD Treg pseudotime (Fig 4G) scored in TAURUS Tregs and
tested before/after anti-TNF with the Fig 6A biopsy-level model. Numbers only, no figure.

Gene sets: the six 'Treg suppression' rows of the Fig 4G pathway table (Fig3D-I_treg_pseudotime_signals/
  Fig3G_pathways_along_pseudotime_stats.csv), taken from the MSigDB snapshot (cd_treg_rebuild_rowfix/archive/csv/
  msigdb_human_snapshot.csv.gz), full sets as in Fig 4G (TNFR2 Reactome WITH the PSM genes).
Scoring: TAURUS raw_norm_counts_anno.h5ad.gz, Tregs = aucell_cell_type_short == 'Treg', decoupler dc.mt.ulm on X = log1p(CP10k)
  ('lognorm' layer), default tmin.
  The Treg cells must equal those of supplementary/cache/Treg.parquet; IL-2/STAT5 (identical gene set) must equal its 'il2_stat5'
  column and TNFR2 without PSM genes must equal its 'tnfr2_pathway' column. The full TNFR2 set (with PSM genes) is the one plotted.
Test ('patient + patient:timepoint' model; no covariates):
  per biopsy (sample_id) the mean of each per-cell score over its Tregs; biopsies with < 10 Tregs dropped;
  biopsy mean ~ post * remission + (1 | patient) + (1 | patient:timepoint); R lme4 + lmerTest, REML, Satterthwaite df; every biopsy
  weighted equally; y standardised before fitting (estimates back-transformed to score units). 95% CI = estimate +/- t(0.975, df) x SE.
  Contrasts: change in remission = b_post + b_post:rem; change in non-remission = b_post; difference in change = b_post:rem;
  pre-treatment difference (remission - non-remission) = b_rem.
  BH over the 6 pathways, per contrast.
  Rows of family 'biopsy_level_prepost.csv definitions' (not in the BH family): IL-2/STAT5 and TNFR2 (no PSM) from the per-cell table; their
  estimate / P must equal supplementary/biopsy_level_prepost.csv.
Side by side: the CD Treg pseudotime result of each pathway (last vs first pseudotime quintile, patient-paired Wilcoxon, BH;
  direction = sign of median change; plus mixed-model bin-effect q and shape).
Output: supplementary/fig3G_pathways_taurus_prepost.csv (3 significant figures; Supp Table 11). Per-cell scores cached in
  supplementary/cache/fig3G_pathways_treg.parquet (re-used if present; delete it to re-score).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time, subprocess, tempfile, shutil
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    HERE = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary'
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
    F3G = '/share/fsmresfiles/UC/revision_jenny/Fig3D-I_treg_pseudotime_signals/Fig3G_pathways_along_pseudotime_stats.csv'
    CACHE29 = f'{HERE}/cache/Treg.parquet'; MYCACHE = f'{HERE}/cache/fig3G_pathways_treg.parquet'
    B54 = f'{HERE}/biopsy_level_prepost.csv'; OUT = f'{HERE}/fig3G_pathways_taurus_prepost.csv'
    MIN_TREG = 10; MOD = 'patient + patient:timepoint'
    KEYS = ['Patient', 'Remission_status', 'Treatment', 'sample_id']
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')

    # ================================================================ 'Treg suppression' pathway block of Fig 4G
    G3 = pd.read_csv(F3G); G3 = G3[G3.group == 'Treg suppression'].reset_index(drop=True)
    assert len(G3) == 6 and (G3.resource == 'Suppression_panel').all()
    SETS = list(G3.gene_set)
    PATH = dict(zip(G3.gene_set, G3.pathway))
    TNFR2 = 'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY'; IL2 = 'HALLMARK_IL2_STAT5_SIGNALING'
    assert {TNFR2, IL2} <= set(SETS)
    msig = pd.read_csv(f'{D}/archive/csv/msigdb_human_snapshot.csv.gz', usecols=['geneset', 'genesymbol'])
    net = (msig[msig.geneset.isin(SETS)].drop_duplicates(['geneset', 'genesymbol'])
           .rename(columns={'geneset': 'source', 'genesymbol': 'target'})[['source', 'target']])
    assert set(net.source) == set(SETS)
    tn = net[(net.source == TNFR2) & ~net.target.str.startswith('PSM')].assign(source='TNFR2_NO_PSM')           # as cell_level_lmm


    # ================================================================ per-cell scores (as cell_level_lmm build_tables)
    def score():
        import scanpy as sc, decoupler as dc
        ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); log('loaded', ad.shape)
        ct = ad.obs['aucell_cell_type_short'].astype(str)
        assert np.allclose(np.expm1(ad.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)                       # X = log1p(CP10k)
        treg = ad[(ct == 'Treg').values].copy(); del ad
        treg.layers['lognorm'] = treg.X.copy()
        S = treg.obs[KEYS].astype(str).copy()
        ngen = {s: int(net[net.source == s].target.isin(treg.var_names).sum()) for s in SETS}
        ngen['TNFR2_NO_PSM'] = int(tn.target.isin(treg.var_names).sum())
        dc.mt.ulm(data=treg, net=net, layer='lognorm'); A = dc.pp.get_obsm(treg, key='score_ulm')
        assert set(A.var_names) == set(SETS), set(SETS) - set(A.var_names)
        for s in SETS: S[s] = np.asarray(A[:, s].X).ravel()
        dc.mt.ulm(data=treg, net=tn, layer='lognorm'); S['TNFR2_NO_PSM'] = np.asarray(dc.pp.get_obsm(treg, key='score_ulm').X).ravel()
        S.attrs = {}
        return S, ngen


    if os.path.exists(MYCACHE):
        S = pd.read_parquet(MYCACHE); NG = pd.read_csv(MYCACHE.replace('.parquet', '_ngenes.csv'), index_col=0).n_genes_scored.to_dict()
        log('loaded per-cell score cache', S.shape)
    else:
        S, NG = score()
        S.to_parquet(MYCACHE)
        with open(MYCACHE.replace('.parquet', '_ngenes.csv'), 'w') as fh:
            pd.Series(NG, name='n_genes_scored').rename_axis('gene_set').to_csv(fh)
        log('scored and cached', S.shape)
    NSET = net.groupby('source').size().to_dict(); NSET['TNFR2_NO_PSM'] = len(tn)

    # ================================================================ same Tregs and scores as the per-cell table of cell_level_lmm
    C = pd.read_parquet(CACHE29)
    assert set(C.index) == set(S.index) and len(C) == len(S), 'Treg cells differ from the cell_level_lmm table'
    C = C.loc[S.index]
    for k in ['Patient', 'Remission_status', 'Treatment', 'sample_id']:
        assert (C[k].astype(str).values == S[k].values).all(), k
    for mine, theirs in [(IL2, 'il2_stat5'), ('TNFR2_NO_PSM', 'tnfr2_pathway')]:
        d = np.abs(S[mine].values - C[theirs].values).max()
        assert np.allclose(S[mine].values, C[theirs].values, rtol=1e-5, atol=1e-6), (mine, d)
        log(f'per-cell {mine} == cache {theirs} (max |diff| {d:.2g}, {len(S)} Tregs)')
    r_t = np.corrcoef(S[TNFR2], S['TNFR2_NO_PSM'])[0, 1]
    log(f'TNFR2 full set (with PSM, {NSET[TNFR2]} genes) vs without PSM ({NSET["TNFR2_NO_PSM"]} genes): per-cell r = {r_t:.3f}')

    # ================================================================ biopsy pseudobulk
    FCOLS = SETS + ['TNFR2_NO_PSM']
    F = S.copy(); F['biopsy'] = F.sample_id.astype(str)
    assert F.groupby('biopsy').Patient.nunique().max() == 1 and F.groupby('biopsy').Treatment.nunique().max() == 1
    F['post'] = (F.Treatment == 'Post').astype(float); F['rem'] = (F.Remission_status == 'Remission').astype(float)
    assert set(F.Treatment) == {'Pre', 'Post'} and set(F.Remission_status) == {'Remission', 'Non_Remission'}
    PB = (F.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                  post=('post', 'first'), rem=('rem', 'first'), n_treg=('post', 'size'))
          .join(F.groupby('biopsy')[FCOLS].mean()).reset_index())
    assert PB.n_treg.sum() == len(F)
    DROP = PB[PB.n_treg < MIN_TREG]; PB = PB[PB.n_treg >= MIN_TREG].reset_index(drop=True)
    log(f'{len(F)} Tregs; kept {len(PB)} biopsies ({int(PB.n_treg.sum())} Tregs, {PB.Patient.nunique()} patients); dropped {len(DROP)} with < {MIN_TREG} Tregs')
    assert (len(PB), int(PB.n_treg.sum()), PB.Patient.nunique()) == (88, 10205, 16), 'biopsy set differs from Fig 6A (88 biopsies, 10,205 Tregs, 16 patients)'

    RCODE = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); d <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); cols <- readLines(a[2]); out <- list()
d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1),
           baseline_difference = c(rem = 1))
for (cc in cols) {
  y <- d[[cc]]; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  m <- lmer(ys ~ post * rem + (1 | Patient) + (1 | pt), data = d, REML = TRUE)
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m))
  for (k in names(Ls)) {
    L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
    ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
    out[[length(out) + 1]] <- data.frame(feature_key = cc, contrast_key = k, estimate = ct$Estimate * sy, se = ct$`Std. Error` * sy,
      df = ct$df, statistic = ct$`t value`, p = ct$`Pr(>|t|)`, var_patient = vc$vcov[vc$grp == "Patient"] * sy^2,
      var_patient_timepoint = vc$vcov[vc$grp == "pt"] * sy^2, var_residual = vc$vcov[vc$grp == "Residual"] * sy^2, singular = isSingular(m),
      converged = !any(grepl("converge|gradient|Hessian|eigenvalue", c(m@optinfo$conv$lme4$messages, m@optinfo$warnings))))
  }
}
write.csv(do.call(rbind, out), a[3], row.names = FALSE)
"""
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        PB.assign(pt=PB.Patient + '_' + PB.Treatment)[['Patient', 'pt', 'post', 'rem'] + FCOLS].to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/cols.txt', 'w') as fh: fh.write('\n'.join(FCOLS) + '\n')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        log(f'fitting {len(FCOLS)} features in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/cols.txt', f'{tmp}/r_out.csv'], check=True)
        R = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(R) == 4 * len(FCOLS) and R.p.notna().all()
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'baseline_difference': 'Baseline difference (remission − non-remission, pre)'}
    R['contrast'] = R.contrast_key.map(CON); R['model'] = MOD
    R['model_formula'] = 'biopsy mean ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df'
    tq = stats.t.ppf(0.975, R.df); R['ci_low'] = R.estimate - tq * R.se; R['ci_high'] = R.estimate + tq * R.se
    FAM3G, FAMCHK = 'Fig 3G Treg suppression pathways (6)', 'biopsy_level_prepost.csv definitions (not in BH family)'
    R['family'] = np.where(R.feature_key.isin(SETS), FAM3G, FAMCHK)
    R['pathway'] = R.feature_key.map(PATH).fillna('TNFR2 non-canonical NF-κB (Reactome, without PSM genes; as biopsy_level_prepost.csv)')
    R['q_bh'] = np.nan
    for c, i in R[R.family == FAM3G].groupby('contrast').groups.items():
        R.loc[i, 'q_bh'] = multipletests(R.loc[i, 'p'], method='fdr_bh')[1]
    R['bh_family'] = np.where(R.family == FAM3G, 'BH over the 6 Fig 3G pathways, per contrast', 'none')
    R['n_genes_in_set'] = R.feature_key.map(NSET); R['n_genes_scored'] = R.feature_key.map(NG)

    # ================================================================ IL-2/STAT5 and TNFR2 (no PSM): same estimate / P as biopsy_level_prepost.csv
    B = pd.read_csv(B54); B = B[B.model == MOD]
    for mine, theirs in [(IL2, 'il2_stat5'), ('TNFR2_NO_PSM', 'tnfr2_pathway')]:
        for c in CON.values():
            a = R[(R.feature_key == mine) & (R.contrast == c)].iloc[0]; b = B[(B.feature_key == theirs) & (B.contrast == c)].iloc[0]
            for k in ['estimate', 'se', 'df', 'p']:
                assert np.isclose(sig3(a[k]), b[k], rtol=2e-3), (mine, c, k, a[k], b[k])
        log(f'{mine}: estimate, SE, df and P of all 4 contrasts equal {theirs} in {os.path.basename(B54)}')
    R['fig6A_existing_p'] = np.nan
    for mine, theirs in [(IL2, 'il2_stat5'), (TNFR2, 'tnfr2_pathway'), ('TNFR2_NO_PSM', 'tnfr2_pathway')]:
        m = B[B.feature_key == theirs].set_index('contrast').p
        i = R.feature_key == mine; R.loc[i, 'fig6A_existing_p'] = R.loc[i, 'contrast'].map(m)
    R['fig6A_existing_feature'] = R.feature_key.map({IL2: 'IL-2/STAT5 (identical gene set)', TNFR2: 'TNFR2 without PSM genes (gene set differs)',
                                                     'TNFR2_NO_PSM': 'TNFR2 without PSM genes (identical gene set)'})

    # ================================================================ CD Treg pseudotime result side by side
    g = G3.set_index('gene_set')
    R['fig3G_direction_late_vs_early'] = R.feature_key.map(lambda s: ('down' if g.loc[s, 'median_delta_last_vs_first_quintile'] < 0 else 'up') if s in g.index else np.nan)
    R['fig3G_median_delta_late_vs_early'] = R.feature_key.map(g.median_delta_last_vs_first_quintile)
    R['fig3G_q_late_vs_early'] = R.feature_key.map(g.q_last_vs_first_quintile_BH)
    R['fig3G_q_bin_effect'] = R.feature_key.map(g.q_bin_effect_mixed_model_BH)
    R['fig3G_trajectory_shape'] = R.feature_key.map(g.trajectory_shape)
    R['fig3G_test'] = np.where(R.family == FAM3G, 'last vs first pseudotime quintile, patient-paired Wilcoxon (25 CD patients), BH over the panel', '')

    # ================================================================ counts and output
    NN = {'n_patients': PB.Patient.nunique(), 'n_biopsies': len(PB), 'n_tregs': int(PB.n_treg.sum()), 'n_biopsies_dropped': len(DROP)}
    for nm, r_, p_ in [('rem_pre', 1.0, 0.0), ('rem_post', 1.0, 1.0), ('nonrem_pre', 0.0, 0.0), ('nonrem_post', 0.0, 1.0)]:
        s = PB[(PB.rem == r_) & (PB.post == p_)]
        NN.update({f'n_biopsies_{nm}': len(s), f'n_patients_{nm}': s.Patient.nunique()})
    for k, v in NN.items(): R[k] = v
    R['min_tregs_per_biopsy'] = MIN_TREG
    COLS = (['family', 'pathway', 'contrast', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'bh_family',
             'fig3G_direction_late_vs_early', 'fig3G_median_delta_late_vs_early', 'fig3G_q_late_vs_early', 'fig3G_q_bin_effect', 'fig3G_trajectory_shape',
             'fig3G_test', 'fig6A_existing_p', 'fig6A_existing_feature', 'model', 'model_formula', 'var_patient', 'var_patient_timepoint', 'var_residual',
             'singular', 'converged', 'n_genes_in_set', 'n_genes_scored', 'min_tregs_per_biopsy'] + list(NN) + ['feature_key'])
    R['_f'] = R.feature_key.map({k: i for i, k in enumerate(FCOLS)}); R['_c'] = R.contrast.map({v: i for i, v in enumerate(CON.values())})
    R = R.sort_values(['_c', '_f'])[COLS].reset_index(drop=True)
    with open(OUT, 'w') as fh:
        R.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, R.shape)

    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 45)
    show = R.copy(); show['pathway'] = show.pathway.str[:40]; show['contrast'] = show.contrast.str[:24]
    print(show[['contrast', 'pathway', 'estimate', 'ci_low', 'ci_high', 'df', 'p', 'q_bh', 'fig3G_direction_late_vs_early', 'fig3G_q_late_vs_early',
                'fig6A_existing_p', 'singular']].apply(lambda c: c.map(sig3)).to_string(index=False))


# ======================================================================================================================
# Section F6A_tnfrsf4_biopsy: TNFRSF4 (OX40) and TNFRSF18 (GITR) in Tregs before and after anti-TNF, biopsy-level model
# Panels: Fig 6A (TNFRSF4 values; Supp Table 11)
# ======================================================================================================================
if _run('F6A_tnfrsf4_biopsy'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary/66_tnfrsf4_prepost.py'   # working folder of this step (relative paths resolve from here)
    """TAURUS anti-TNF, CD scRNA-seq: TNFRSF4 (OX40) and TNFRSF18 (GITR) in Tregs before and after anti-TNF, remitters vs
non-remitters, with the biopsy as the unit. Numbers only, no figure (TNFRSF4 is drawn by F6A_treg_state_lines).
Why these genes: TNFRSF4 and TNFRSF18 are the top early-bin (Activated) Treg genes of Fig 4F (fig3F_marker_dotplot.csv); TNFRSF4 is the
  top driver target of YY1 (Supp Fig 8F, ext8_tf_activity_along_pseudotime_drivers.csv); OX40L -> OX40 is a Fig 5C co-stimulatory pair.
  TNFRSF4 = Fig 6A panel; TNFRSF18 = CSV-only companion.
Method: as TNFRSF9 in cell_level_lmm and the biopsy-level sections.
  Per-cell value: log-normalised expression (X of the TAURUS h5ad = log1p(CP10k)) in Tregs (aucell_cell_type_short == 'Treg').
    TNFRSF9 extracted the same way must equal the gene_TNFRSF9 column of supplementary/cache/Treg.parquet.
  Pseudobulk: per biopsy (sample_id), mean over its Tregs; biopsies with < 10 Tregs dropped (the same 88 biopsies as the other
    biopsy-level sections).
  Model: biopsy mean ~ post * remission + (1 | patient) + (1 | patient:timepoint); R lme4 + lmerTest, REML, Satterthwaite df; every biopsy
    weighted equally; y standardised before fitting (estimates back-transformed to log-normalised expression units).
  Contrasts (lmerTest contest): change in remission = b_post + b_post:rem; change in non-remission = b_post; differential change = b_post:rem;
    baseline (pre-treatment) difference = b_rem; post-treatment difference = b_rem + b_post:rem. 95% CI = estimate +/- t(0.975, df) x SE.
  TNFRSF9 is fitted with the same code and must equal (3 s.f.) its 5 rows of supplementary/biopsy_level_prepost.csv.
  BH (q_bh): TNFRSF4 joins the Fig 6A activation-program family, 4 tests per contrast: IL-2/STAT5 and TNFR2 non-canonical NF-kB (full
    Reactome set) from supplementary/fig3G_pathways_taurus_prepost.csv (F6A_pathways_biopsy), the co-stimulatory/suppressive receptor
    module from supplementary/biopsy_level_prepost.csv (F6A_programs_biopsy), and TNFRSF4 (here); statsmodels fdr_bh per contrast.
    F6A_pathways_biopsy has no post-treatment difference rows, so that contrast has no q. TNFRSF18 is not in a BH family (CSV-only
    companion): q_bh empty.
  % of Tregs expressing (log-normalised expression > 0): over all Tregs of the 88 biopsies, and per group x timepoint (pooled Tregs);
    also the observed mean of the biopsy means per group x timepoint.
Output: supplementary/tnfrsf4_tnfrsf18_prepost.csv (3 significant figures; Supp Table 11).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file>"""
    import os, sys, warnings, time, subprocess, tempfile, shutil
    warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    HERE = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/supplementary'
    T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
    CACHE = f'{HERE}/cache/Treg.parquet'; O54 = f'{HERE}/biopsy_level_prepost.csv'; O65 = f'{HERE}/fig3G_pathways_taurus_prepost.csv'
    OUT = f'{HERE}/tnfrsf4_tnfrsf18_prepost.csv'
    MIN_TREG = 10; MOD = 'patient + patient:timepoint'
    GENES = {'TNFRSF4': 'TNFRSF4 (OX40)', 'TNFRSF18': 'TNFRSF18 (GITR)', 'TNFRSF9': 'TNFRSF9 (4-1BB)'}
    CON = {'change_remission': 'Change in remission (post − pre)', 'change_nonremission': 'Change in non-remission (post − pre)',
           'differential_change': 'Differential change (remission − non-remission)',
           'baseline_difference': 'Baseline difference (remission − non-remission, pre)',
           'post_difference': 'Post-treatment difference (remission − non-remission, post)'}
    IL2, TNFR2 = 'HALLMARK_IL2_STAT5_SIGNALING', 'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY'
    t0 = time.time(); log = lambda *a: print(f'[{time.time() - t0:5.0f}s]', *a, flush=True)
    sig3 = lambda v: v if not isinstance(v, (float, np.floating)) or not np.isfinite(v) else float(f'{v:.3g}')

    # ================================================================ per-cell log-normalised expression
    ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); ct = ad.obs['aucell_cell_type_short'].astype(str).values
    treg = ad[ct == 'Treg'].copy(); del ad; log('Tregs', treg.n_obs)
    assert np.allclose(np.expm1(treg.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)                      # X = log1p(CP10k)
    F = pd.read_parquet(CACHE); assert set(F.index) == set(treg.obs_names) and len(F) == treg.n_obs
    F = F.loc[treg.obs_names].copy()
    assert (F.Patient.astype(str).values == treg.obs.Patient.astype(str).values).all() and (F.Treatment.astype(str).values == treg.obs.Treatment.astype(str).values).all()
    for g in GENES:
        assert g in treg.var_names, g
        x = treg[:, g].X; F[f'gene_{g}'] = np.asarray(x.todense() if hasattr(x, 'todense') else x).ravel()
    assert np.allclose(F['gene_TNFRSF9'].values, pd.read_parquet(CACHE).loc[F.index, 'gene_TNFRSF9'].values, rtol=0, atol=1e-6), 'TNFRSF9 differs from the cell_level_lmm table'
    del treg
    FCOLS = [f'gene_{g}' for g in GENES]

    # ================================================================ biopsy pseudobulk
    F['biopsy'] = F.sample_id.astype(str)
    assert F.groupby('biopsy').Patient.nunique().max() == 1 and F.groupby('biopsy').Treatment.nunique().max() == 1
    F['post'] = (F.Treatment == 'Post').astype(float); F['rem'] = (F.Remission_status == 'Remission').astype(float)
    PB = (F.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                  post=('post', 'first'), rem=('rem', 'first'), n_treg=('post', 'size'))
          .join(F.groupby('biopsy')[FCOLS].mean()).reset_index())
    DROP = PB[PB.n_treg < MIN_TREG]; PB = PB[PB.n_treg >= MIN_TREG].reset_index(drop=True)
    A54 = pd.read_csv(O54); ref = A54[A54.model == MOD].iloc[0]
    assert len(PB) == ref.n_biopsies == 88 and int(PB.n_treg.sum()) == ref.n_tregs and PB.Patient.nunique() == ref.n_patients and len(DROP) == ref.n_biopsies_dropped
    FK = F[F.biopsy.isin(PB.biopsy)]; assert len(FK) == ref.n_tregs
    log(f'{len(PB)} biopsies, {len(FK)} Tregs, {PB.Patient.nunique()} patients; dropped {len(DROP)}')
    GRP = [('rem_pre', 1.0, 0.0), ('rem_post', 1.0, 1.0), ('nonrem_pre', 0.0, 0.0), ('nonrem_post', 0.0, 1.0)]
    NN = {}
    for nm, r_, p_ in GRP:
        s = PB[(PB.rem == r_) & (PB.post == p_)]
        NN.update({f'n_biopsies_{nm}': len(s), f'n_patients_{nm}': s.Patient.nunique()})
        assert NN[f'n_biopsies_{nm}'] == ref[f'n_biopsies_{nm}'] and NN[f'n_patients_{nm}'] == ref[f'n_patients_{nm}'], nm

    # ================================================================ lme4 / lmerTest, Satterthwaite
    RCODE = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE); d <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE); cols <- readLines(a[2]); out <- list()
d$Patient <- factor(d$Patient); d$pt <- factor(d$pt)
Ls <- list(change_remission = c(post = 1, `post:rem` = 1), change_nonremission = c(post = 1), differential_change = c(`post:rem` = 1),
           baseline_difference = c(rem = 1), post_difference = c(rem = 1, `post:rem` = 1))
for (cc in cols) {
  y <- d[[cc]]; sy <- sd(y); d$ys <- (y - mean(y)) / sy
  m <- lmer(ys ~ post * rem + (1 | Patient) + (1 | pt), data = d, REML = TRUE)
  fe <- names(fixef(m)); vc <- as.data.frame(VarCorr(m))
  for (k in names(Ls)) {
    L <- setNames(rep(0, length(fe)), fe); L[names(Ls[[k]])] <- Ls[[k]]
    ct <- contest(m, L = matrix(L, nrow = 1), joint = FALSE, ddf = "Satterthwaite")
    out[[length(out) + 1]] <- data.frame(feature_key = cc, contrast_key = k, estimate = ct$Estimate * sy, se = ct$`Std. Error` * sy,
      df = ct$df, statistic = ct$`t value`, p = ct$`Pr(>|t|)`, var_patient = vc$vcov[vc$grp == "Patient"] * sy^2,
      var_patient_timepoint = vc$vcov[vc$grp == "pt"] * sy^2, var_residual = vc$vcov[vc$grp == "Residual"] * sy^2, singular = isSingular(m),
      converged = !any(grepl("converge|gradient|Hessian|eigenvalue", c(m@optinfo$conv$lme4$messages, m@optinfo$warnings))))
  }
}
write.csv(do.call(rbind, out), a[3], row.names = FALSE)
"""
    tmp = tempfile.mkdtemp(prefix='taurus_lmm_')
    try:
        PB.assign(pt=PB.Patient + '_' + PB.Treatment)[['Patient', 'pt', 'post', 'rem'] + FCOLS].to_csv(f'{tmp}/pb.csv', index=False, float_format='%.12g')
        with open(f'{tmp}/cols.txt', 'w') as fh: fh.write('\n'.join(FCOLS) + '\n')
        with open(f'{tmp}/x.R', 'w') as fh: fh.write(RCODE)
        log(f'fitting {len(FCOLS)} genes in lme4 / lmerTest ...')
        subprocess.run(['Rscript', f'{tmp}/x.R', f'{tmp}/pb.csv', f'{tmp}/cols.txt', f'{tmp}/r_out.csv'], check=True)
        R = pd.read_csv(f'{tmp}/r_out.csv')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    assert len(R) == 5 * len(FCOLS) and R.p.notna().all()
    R['contrast'] = R.contrast_key.map(CON)
    tq = stats.t.ppf(0.975, R.df); R['ci_low'] = R.estimate - tq * R.se; R['ci_high'] = R.estimate + tq * R.se

    # ================================================================ TNFRSF9 must equal its 5 rows of biopsy_level_prepost.csv (same model and biopsies)
    up = A54[(A54.model == MOD) & (A54.feature_key == 'gene_TNFRSF9')].set_index('contrast'); assert len(up) == 5
    fit = R[R.feature_key == 'gene_TNFRSF9'].set_index('contrast').reindex(up.index)
    for c in ['estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'var_patient', 'var_residual']:
        assert np.allclose(fit[c].map(sig3).astype(float).values, up[c].astype(float).values, rtol=2e-3, atol=1e-12), c
    log('TNFRSF9 matches its 5 rows of biopsy_level_prepost.csv (estimate, CI, SE, df, t, P, variance components)')
    R = R[R.feature_key != 'gene_TNFRSF9'].copy()

    # ================================================================ BH: Fig 6A activation programs (4) per contrast, TNFRSF4 only
    G = pd.read_csv(O65); G = G[(G.model == MOD) & (G.family == 'Fig 3G Treg suppression pathways (6)') & G.feature_key.isin([IL2, TNFR2])]
    CM = A54[(A54.model == MOD) & (A54.feature_key == 'costim_module')]
    R['q_bh'] = np.nan; R['bh_family'] = ''
    for k, cname in CON.items():
        if k == 'post_difference': continue
        others = pd.concat([G[G.contrast == cname][['feature_key', 'p']], CM[CM.contrast == cname][['feature_key', 'p']]])
        assert len(others) == 3, (k, len(others))
        i = R.index[(R.feature_key == 'gene_TNFRSF4') & (R.contrast_key == k)]; assert len(i) == 1
        ps = np.r_[others.p.values, R.loc[i, 'p'].values]; R.loc[i, 'q_bh'] = multipletests(ps, method='fdr_bh')[1][-1]
        R.loc[i, 'bh_family'] = 'Fig 6A activation programs (4: IL-2/STAT5, TNFR2 full set, receptor module, TNFRSF4), per contrast'
    R.loc[(R.feature_key == 'gene_TNFRSF4') & (R.contrast_key == 'post_difference'), 'bh_family'] = 'none (no post-treatment difference rows for the other 3 programs)'
    R.loc[R.feature_key == 'gene_TNFRSF18', 'bh_family'] = 'none (CSV-only companion, unadjusted P)'

    # ================================================================ % of Tregs expressing, observed group means
    EX = {}
    for g in ['TNFRSF4', 'TNFRSF18']:
        c = f'gene_{g}'; e = {'pct_tregs_expressing_all': 100 * (FK[c] > 0).mean()}
        for nm, r_, p_ in GRP:
            s = FK[(FK.rem == r_) & (FK.post == p_)]; e[f'pct_tregs_expressing_{nm}'] = 100 * (s[c] > 0).mean(); e[f'n_tregs_{nm}'] = len(s)
            b = PB[(PB.rem == r_) & (PB.post == p_)][c]; e[f'mean_biopsy_mean_{nm}'] = b.mean(); e[f'sem_biopsy_mean_{nm}'] = b.std(ddof=1) / np.sqrt(len(b))
        EX[c] = e
    R = R.join(pd.DataFrame.from_dict(EX, orient='index'), on='feature_key')
    R['feature'] = R.feature_key.str[5:].map(GENES)
    R['model'] = MOD
    R['model_formula'] = 'biopsy mean ~ post × remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df'
    R['value'] = 'per-Treg log-normalised expression (log1p CP10k), averaged per biopsy'
    R['min_tregs_per_biopsy'] = MIN_TREG; R['n_patients'] = PB.Patient.nunique(); R['n_biopsies'] = len(PB); R['n_tregs'] = len(FK)
    for k, v in NN.items(): R[k] = v
    R['n_biopsies_dropped'] = len(DROP)
    EXC = [c for c in R.columns if c.startswith(('pct_tregs', 'mean_biopsy', 'sem_biopsy', 'n_tregs_'))]
    COLS = (['feature', 'contrast', 'estimate', 'ci_low', 'ci_high', 'se', 'df', 'statistic', 'p', 'q_bh', 'bh_family', 'model', 'model_formula', 'value']
            + EXC + ['var_patient', 'var_patient_timepoint', 'var_residual', 'singular', 'converged', 'min_tregs_per_biopsy', 'n_patients', 'n_biopsies',
                     'n_tregs'] + list(NN) + ['n_biopsies_dropped', 'feature_key'])
    R['_f'] = R.feature_key.map({'gene_TNFRSF4': 0, 'gene_TNFRSF18': 1}); R['_c'] = R.contrast_key.map({k: i for i, k in enumerate(CON)})
    R = R.sort_values(['_f', '_c'])[COLS].reset_index(drop=True)
    assert len(R) == 10 and R.converged.all()
    with open(OUT, 'w') as fh:
        R.apply(lambda c: c.map(sig3)).to_csv(fh, index=False)
    log('wrote', OUT, R.shape)
    pd.set_option('display.width', 250)
    print(R[['feature', 'contrast', 'estimate', 'ci_low', 'ci_high', 'df', 'p', 'q_bh', 'singular']].apply(lambda c: c.map(sig3)).to_string(index=False))
    print(R.drop_duplicates('feature')[['feature'] + EXC].apply(lambda c: c.map(sig3)).T.to_string())


# ======================================================================================================================
# Section F6A_treg_state_lines: Treg suppression state and Treg activation programs before and after anti-TNF (figure and legend)
# Panels: Fig 6A
# Note: refuses to overwrite an existing PNG / legend; pass --check <dir> to write them to another folder
# ======================================================================================================================
if _run('F6A_treg_state_lines'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig6_treg_antitnf_taurus/58g_fig6_treg_state_lines.py'   # working folder of this step (relative paths resolve from here)
    """Fig 6A: Treg suppression state and Treg activation programs before and after anti-TNF in TAURUS Crohn's disease scRNA-seq,
remitters vs non-remitters. Plot only: statistics are read from the biopsy-level CSVs, not refitted.
Layout (~15.4 in wide, two rows; fonts ~14 pt): row 1 'Suppression State': Treg Suppression Score, colour key + one-line P key to its
  right. Row 2 'Activation Programs': IL-2/STAT5 Signaling, TNFR2 Non-Canonical NF-κB, Co-Stimulatory/Suppressive Receptor Module,
  TNFRSF4 (OX40). TNFRSF4: with TNFRSF18 (GITR) the top early-bin (Activated) Treg gene in Fig 4F; the top driver target of YY1, the
  top falling TF in Supplementary Fig 8F; OX40L -> OX40 is one of the Fig 5C co-stimulatory pairs.
Statistics (read, not refitted), biopsy-level model 'patient + patient:timepoint':
  biopsy mean ~ post x remission + (1 | patient) + (1 | patient:timepoint); lme4/lmerTest REML, Satterthwaite df; biopsies >= 10 Tregs.
  Treg Suppression Score: supplementary/treg_suppression_biopsy_level.csv (F6A_suppression_score_biopsy). IL-2/STAT5 and TNFR2 (full
  set): supplementary/fig3G_pathways_taurus_prepost.csv (F6A_pathways_biopsy; IL-2/STAT5 must equal the F6A_programs_biopsy row,
  identical gene set). Co-stimulatory module and the pseudotime source-data note: supplementary/biopsy_level_prepost.csv
  (F6A_programs_biopsy, F6A_post_treatment_difference). TNFRSF4: supplementary/tnfrsf4_tnfrsf18_prepost.csv (F6A_tnfrsf4_biopsy;
  its q must equal the BH computed here).
  On the plot: the change after treatment within each group (2 s.f., 'P < 0.001' below 0.001; no q).
  Legend only: difference in change and pre-treatment difference, P and BH q over the 4 activation programs (IL-2/STAT5, TNFR2 full
  set, receptor module, TNFRSF4; statsmodels fdr_bh on the CSV P values), per contrast.
Points: observed biopsy means (each biopsy's mean per-cell value over its Tregs, >= 10 Tregs), averaged over the biopsies of each group x
  timepoint; error bars, SEM across biopsies. Biopsy/patient/Treg counts must equal all four CSVs.
Per-cell values: pseudotime and costim_module from supplementary/cache/Treg.parquet (cell_level_lmm); IL-2/STAT5 and TNFR2 (full) from
  supplementary/cache/fig3G_pathways_treg.parquet (F6A_pathways_biopsy); TNFRSF4 from X of the TAURUS h5ad (log1p CP10k). The Treg
  Suppression Score is computed per cell as in F6A_suppression_score_biopsy; group means must equal the CSVs.
Outputs: fig6_treg_state_lines.png + fig6_treg_state_lines_legend.txt (never overwritten; --check <dir> writes them to <dir> instead).
Run: OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=32 /data/ydn4687/geneformer_env/bin/python <this file> [--check <outdir>]"""
    import os, sys, glob, warnings
    warnings.filterwarnings('ignore')
    import logging; logging.getLogger('fontTools').setLevel(logging.WARNING)
    import numpy as np, pandas as pd, matplotlib, scanpy as sc
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import matplotlib.transforms as mtrans
    from decimal import Decimal, ROUND_HALF_UP
    from statsmodels.stats.multitest import multipletests

    R = '/share/fsmresfiles/UC/revision_jenny'; D = f'{R}/cd_treg_rebuild_rowfix'; T = '/share/fsmresfiles/UC/NatureImmunology_TAURUS_forJenny'
    OUT = f'{R}/Fig6_treg_antitnf_taurus'; STEM = f'{OUT}/fig6_treg_state_lines'; SUP = f'{OUT}/supplementary'
    if '--check' in sys.argv: STEM = f"{sys.argv[sys.argv.index('--check') + 1]}/fig6_treg_state_lines"; os.makedirs(os.path.dirname(STEM), exist_ok=True)
    for f in [f'{STEM}.png', f'{STEM}_legend.txt']: assert not os.path.exists(f), f'refusing to overwrite {f}'
    CACHE = f'{SUP}/cache/Treg.parquet'; PCACHE = f'{SUP}/cache/fig3G_pathways_treg.parquet'; MIN_TREG = 10
    MOD = 'patient + patient:timepoint'
    CON = {'rem': 'Change in remission (post − pre)', 'non': 'Change in non-remission (post − pre)',
           'diff': 'Differential change (remission − non-remission)', 'base': 'Baseline difference (remission − non-remission, pre)',
           'post': 'Post-treatment difference (remission − non-remission, post)'}
    IL2, TNFR2 = 'HALLMARK_IL2_STAT5_SIGNALING', 'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY'

    # ---------------------------------------------------------------- panels: (key, source, CSV feature name, panel title, y label)
    PANELS = [('supp_score', 's55', 'Treg Suppression Score (early − late, score_genes, biopsy WTX gene sets)', 'Treg Suppression Score', 'Treg Suppression Score'),
              (IL2, 's65', 'IL-2/STAT5 signaling (Hallmark)', 'IL-2/STAT5 Signaling', 'Activity (ULM score)'),
              (TNFR2, 's65', 'TNFR2 non-canonical NF-κB (Reactome)', 'TNFR2 Non-Canonical\nNF-κB', 'Activity (ULM score)'),
              ('costim_module', 's54', 'Co-Stimulatory/Suppressive Receptor Module', 'Co-Stimulatory/Suppressive\nReceptor Module', 'Module score'),
              ('gene_TNFRSF4', 's66', 'TNFRSF4 (OX40)', 'TNFRSF4 (OX40)', 'Log-normalised expression')]
    NOTE = ('pseudotime', 's54', 'Treg Pseudotime')                       # legend source-data note only
    PK = [p[0] for p in PANELS]; AK = PK[1:]                               # plotted panels; the 4 activation programs
    AP3 = [IL2, TNFR2, 'costim_module']                                    # the 3 pathway/module programs
    KEYS = PK + ['pseudotime']

    # ---------------------------------------------------------------- biopsy-level statistics
    S5 = pd.read_csv(f'{SUP}/treg_suppression_biopsy_level.csv'); S5 = S5[(S5.model == MOD) & (S5.feature_key == 'supp_score')].copy()
    G = pd.read_csv(f'{SUP}/fig3G_pathways_taurus_prepost.csv')
    G = G[(G.model == MOD) & (G.family == 'Fig 3G Treg suppression pathways (6)') & G.feature_key.isin([IL2, TNFR2])].copy()
    B = pd.read_csv(f'{SUP}/biopsy_level_prepost.csv'); B = B[(B.model == MOD) & B.feature_key.isin(['costim_module', 'pseudotime', 'il2_stat5'])].copy()
    S6 = pd.read_csv(f'{SUP}/tnfrsf4_tnfrsf18_prepost.csv'); S6 = S6[(S6.model == MOD) & (S6.feature_key == 'gene_TNFRSF4')].copy()
    SRC = {'s55': S5, 's65': G, 's54': B, 's66': S6}
    for k, src, feat, *_ in PANELS + [NOTE]:
        s = SRC[src]; s = s[s.feature_key == k]; ncon = 4 if src == 's65' else 5
        assert len(s) == ncon and set(s.contrast) <= set(CON.values()) and (s[('pathway' if src == 's65' else 'feature')] == feat).all(), (k, len(s))
    # the 2 pathways are the IL-2 and TNF pathways that fall along the CD Treg pseudotime (Fig 4G; BH q < 0.05, direction down)
    for k in [IL2, TNFR2]:
        g = G[G.feature_key == k].iloc[0]; assert g.fig3G_direction_late_vs_early == 'down' and g.fig3G_q_late_vs_early < 0.05, k
    assert int(G[G.feature_key == TNFR2].n_genes_in_set.iloc[0]) == 188, 'TNFR2 must be the full Reactome set (with PSM genes)'
    # IL-2/STAT5: identical gene set in F6A_programs_biopsy and F6A_pathways_biopsy -> identical model results
    for c in ['rem', 'non', 'diff', 'base']:
        a_ = G[(G.feature_key == IL2) & (G.contrast == CON[c])].iloc[0]; b_ = B[(B.feature_key == 'il2_stat5') & (B.contrast == CON[c])].iloc[0]
        assert a_.p == b_.p and a_.estimate == b_.estimate, ('IL-2/STAT5 differs between F6A_programs_biopsy and F6A_pathways_biopsy', c)
    ST = pd.concat([S5, G.rename(columns={'pathway': 'feature'}), B[B.feature_key.isin(['costim_module', 'pseudotime'])], S6], ignore_index=True)
    NCOL = ['n_patients', 'n_biopsies', 'n_tregs'] + [f'n_{a}_{nm}' for nm in ['rem_pre', 'rem_post', 'nonrem_pre', 'nonrem_post'] for a in ['biopsies', 'patients']]
    assert (ST[NCOL].nunique() == 1).all(), 'the biopsy-level sections analysed different biopsies'
    first = ST.iloc[0]
    W = ST.pivot(index='feature_key', columns='contrast', values=['estimate', 'p'])

    # ---------------------------------------------------------------- per-cell values
    F = pd.read_parquet(CACHE)
    PF = pd.read_parquet(PCACHE)
    assert set(PF.index) == set(F.index) and len(PF) == len(F)
    PF = PF.loc[F.index]
    for c in ['Patient', 'Remission_status', 'Treatment', 'sample_id']: assert (PF[c].astype(str).values == F[c].astype(str).values).all(), c
    assert np.allclose(PF[IL2].values, F.il2_stat5.values, rtol=0, atol=1e-6), 'IL-2/STAT5 scores differ between caches'
    F[IL2] = PF[IL2].values; F[TNFR2] = PF[TNFR2].values
    ad = sc.read_h5ad(f'{T}/raw_norm_counts_anno.h5ad.gz'); ct = ad.obs['aucell_cell_type_short'].astype(str).values
    treg = ad[ct == 'Treg'].copy(); del ad
    assert np.allclose(np.expm1(treg.X[:200].toarray()).sum(1), 1e4, rtol=1e-3)                      # X = log1p(CP10k)
    SGF = glob.glob(f'{D}/panels/fig4_and_related/D_shared_inputs/treg_*_signature_genes.csv'); assert len(SGF) == 1, SGF
    SG = pd.read_csv(SGF[0]); s_ = SG[SG.cohort == 'biopsy']
    GS = {side: [g for g in s_.loc[s_.set.str.startswith(side), 'gene'] if g in treg.var_names] for side in ['early', 'late']}
    NG = {side: (len(GS[side]), int(s_.set.str.startswith(side).sum())) for side in ['early', 'late']}
    assert NG == {'early': (49, 50), 'late': (50, 50)}, NG                                              # genes present in TAURUS
    sc.tl.score_genes(treg, GS['late'], score_name='l', random_state=0)
    sc.tl.score_genes(treg, GS['early'], score_name='e', random_state=0)
    assert set(F.index) == set(treg.obs_names) and len(F) == treg.n_obs
    F = F.loc[treg.obs_names].copy()
    assert (F.Patient.astype(str).values == treg.obs.Patient.astype(str).values).all() and (F.Treatment.astype(str).values == treg.obs.Treatment.astype(str).values).all()
    F['supp_score'] = (treg.obs.e - treg.obs.l).values
    x_ = treg[:, 'TNFRSF4'].X; F['gene_TNFRSF4'] = np.asarray(x_.todense() if hasattr(x_, 'todense') else x_).ravel(); del treg
    assert F[KEYS].notna().all().all()
    print(f'{len(F)} Tregs; Treg Suppression Score computed per cell', flush=True)

    # ---------------------------------------------------------------- biopsy pseudobulk -> observed group means
    F['biopsy'] = F.sample_id.astype(str)
    assert F.groupby('biopsy').Patient.nunique().max() == 1 and F.groupby('biopsy').Treatment.nunique().max() == 1
    PB = (F.groupby('biopsy').agg(Patient=('Patient', 'first'), Remission_status=('Remission_status', 'first'), Treatment=('Treatment', 'first'),
                                  n_treg=('Patient', 'size')).join(F.groupby('biopsy')[KEYS].mean()).reset_index())
    PB = PB[PB.n_treg >= MIN_TREG].reset_index(drop=True)
    assert len(PB) == int(first.n_biopsies) and int(PB.n_treg.sum()) == int(first.n_tregs) and PB.Patient.nunique() == int(first.n_patients)
    GK = {('rem', 'Pre'): ('Remission', 'Pre', 'rem_pre'), ('rem', 'Post'): ('Remission', 'Post', 'rem_post'),
          ('non', 'Pre'): ('Non_Remission', 'Pre', 'nonrem_pre'), ('non', 'Post'): ('Non_Remission', 'Post', 'nonrem_post')}
    NB, NP = {}, {}
    for g, (rs, tp, nm) in GK.items():
        s = PB[(PB.Remission_status == rs) & (PB.Treatment == tp)]; NB[g] = len(s); NP[g] = s.Patient.nunique()
        assert NB[g] == int(first[f'n_biopsies_{nm}']) and NP[g] == int(first[f'n_patients_{nm}']), (g, NB[g], NP[g])
        m = s.supp_score.mean(); ref = float(S5[f'mean_{nm}'].iloc[0])
        assert float(f'{m:.3g}') == ref, ('Treg Suppression Score group mean differs from F6A_suppression_score_biopsy', nm, m, ref)
        m = s.gene_TNFRSF4.mean(); ref = float(S6[f'mean_biopsy_mean_{nm}'].iloc[0])
        assert float(f'{m:.3g}') == ref, ('TNFRSF4 group mean differs from F6A_tnfrsf4_biopsy', nm, m, ref)
    rows = []
    for k in KEYS:
        r = dict(key=k)
        for g, (rs, tp, _) in GK.items():
            v = PB.loc[(PB.Remission_status == rs) & (PB.Treatment == tp), k]
            r[f'{g[0]}_{g[1]}'] = v.mean(); r[f'{g[0]}_{g[1]}_se'] = v.std(ddof=1) / np.sqrt(len(v))
        for c in CON:
            has = ('p', CON[c]) in W.columns and pd.notna(W.loc[k, ('p', CON[c])])
            r[f'{c}_est'] = W.loc[k, ('estimate', CON[c])] if has else np.nan; r[f'{c}_p'] = W.loc[k, ('p', CON[c])] if has else np.nan
        r['obs_rem_change'] = r['rem_Post'] - r['rem_Pre']; r['obs_non_change'] = r['non_Post'] - r['non_Pre']
        r['obs_diff'] = r['obs_rem_change'] - r['obs_non_change']; r['obs_base'] = r['rem_Pre'] - r['non_Pre']
        rows.append(r)
    P = pd.DataFrame(rows).set_index('key')
    # BH q over the 4 activation programs, per contrast
    assert len(AK) == 4
    for c in ['rem', 'non', 'diff', 'base']:
        P.loc[AK, f'{c}_q'] = multipletests(P.loc[AK, f'{c}_p'].values, method='fdr_bh')[1]
        q6 = float(S6[S6.contrast == CON[c]].q_bh.iloc[0])
        assert np.isclose(P.loc['gene_TNFRSF4', f'{c}_q'], q6, rtol=5e-3), ('TNFRSF4 q differs from F6A_tnfrsf4_biopsy', c)   # its BH used unrounded P
    assert P.loc[PK, ['rem_p', 'non_p', 'diff_p', 'base_p']].notna().all().all()
    for k in PK:   # the model change agrees in sign with the observed change wherever the change is significant
        for g in ['rem', 'non']:
            if P.loc[k, f'{g}_p'] < 0.05: assert np.sign(P.loc[k, f'{g}_est']) == np.sign(P.loc[k, f'obs_{g}_change']), (k, g)
    pd.set_option('display.width', 250)
    print(P[['rem_Pre', 'rem_Post', 'non_Pre', 'non_Post', 'obs_rem_change', 'rem_est', 'rem_p', 'obs_non_change', 'non_est', 'non_p',
             'obs_base', 'base_est', 'base_p', 'base_q', 'obs_diff', 'diff_est', 'diff_p', 'diff_q']].round(4).to_string(), flush=True)

    def sig2(p):   # 2 significant figures, half-up on the CSV's 3-significant-figure value, trailing zero kept
        d = Decimal(repr(float(p))); q = Decimal(1).scaleb(d.adjusted() - 1)
        return str(d.quantize(q, rounding=ROUND_HALF_UP))
    def fmt_p(p, s='P'):
        return f'{s} < 0.001' if p < 0.001 else f'{s} = {sig2(float(f"{p:.3g}"))}'

    # ---------------------------------------------------------------- figure (fonts ~14 pt; two rows)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14, 'xtick.labelsize': 13, 'ytick.labelsize': 13,
                         'axes.titleweight': 'normal', 'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0})
    BLUE, RED, K, GR = '#2166ac', '#b2182b', '#111111', '#8c8c8c'
    DX = 0.06; FS_P = 12; DPI = 300
    XLIM = [-0.45, 1.45]
    CLR = 5.0
    TIERS = [(20.0, 8.0, 60.0), (20.0, 0.5, 60.0), (40.0, 0.5, 5.0)]
    NEAR = max(t[0] for t in TIERS)

    WIN, HIN = 15.4, 9.4
    fig = plt.figure(figsize=(WIN, HIN), dpi=DPI)
    L, Rr, WS = 1.05 / WIN, 1 - 0.25 / WIN, 0.44
    def fy(inch): return inch / HIN
    gs2 = fig.add_gridspec(1, 4, left=L, right=Rr, bottom=fy(0.80), top=fy(3.50), wspace=WS)
    gs1 = fig.add_gridspec(1, 4, left=L, right=Rr, bottom=fy(5.30), top=fy(8.00), wspace=WS)
    axs = [fig.add_subplot(gs1[0, 0])] + [fig.add_subplot(gs2[0, j]) for j in range(4)]
    ANN, DATA = {}, {}
    for ax, (k, src, feat, title, ylab) in zip(axs, PANELS):
        r = P.loc[k]; DATA[k] = []
        for g, col, dx in [('non', RED, DX), ('rem', BLUE, -DX)]:
            x = np.array([0, 1]) + dx; y = np.array([r[f'{g}_Pre'], r[f'{g}_Post']]); e = np.array([r[f'{g}_Pre_se'], r[f'{g}_Post_se']])
            ax.errorbar(x, y, yerr=e, color=col, lw=2.0, marker='o', ms=6.5, capsize=3, elinewidth=1.4, capthick=1.4,
                        zorder=3 if g == 'rem' else 2)
            DATA[k].append((g, x, y, e))
        ax.set_xlim(*XLIM); ax.set_xticks([0, 1]); ax.set_xticklabels(['Pre-\ntreatment', 'Post-\ntreatment'], linespacing=1.0)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylabel(ylab, fontsize=14, labelpad=4, linespacing=1.1)
        a = {}
        a['rem'] = ax.text(0.5, r.rem_Pre, fmt_p(r.rem_p), ha='center', va='center', fontsize=FS_P, color=BLUE, zorder=5)
        a['non'] = ax.text(0.5, r.non_Pre, fmt_p(r.non_p), ha='center', va='center', fontsize=FS_P, color=RED, zorder=5)
        a['title'] = ax.annotate(title, xy=(0.5, 1), xycoords='axes fraction', xytext=(0, 8), textcoords='offset points',
                                 ha='center', va='bottom', fontsize=15, linespacing=1.1)
        ANN[k] = a
        lo_ = min(r[f'{g}_{t}'] - r[f'{g}_{t}_se'] for g in ['rem', 'non'] for t in ['Pre', 'Post'])
        hi_ = max(r[f'{g}_{t}'] + r[f'{g}_{t}_se'] for g in ['rem', 'non'] for t in ['Pre', 'Post'])
        ax._lohi = (lo_, hi_); ax._pad = [0.10, 0.10]

    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    def bb(t, pad=0.0):
        e = t.get_window_extent(rend); return mtrans.Bbox.from_extents(e.x0 - pad, e.y0 - pad, e.x1 + pad, e.y1 + pad)
    def ov(a, b, pad=1.0):
        return not (a.x1 + pad <= b.x0 or b.x1 + pad <= a.x0 or a.y1 + pad <= b.y0 or b.y1 + pad <= a.y0)
    CAP = 3 * DPI / 72
    MS = 6.5 * DPI / 72 / 2
    def data_boxes(ax, k):
        boxes, pts = [], {}
        for g, x, y, e in DATA[k]:
            for xi, yi, ei in zip(x, y, e):
                (px, py0), (_, py1) = ax.transData.transform([(xi, yi - ei), (xi, yi + ei)])
                boxes.append(mtrans.Bbox.from_extents(px - max(CAP, MS), py0, px + max(CAP, MS), py1))
                pm = ax.transData.transform((xi, yi)); boxes.append(mtrans.Bbox.from_extents(pm[0] - MS, pm[1] - MS, pm[0] + MS, pm[1] + MS))
            t_ = np.linspace(0, 1, 400)[:, None]
            pts[g] = ax.transData.transform(np.column_stack([x[0] + t_[:, 0] * (x[1] - x[0]), y[0] + t_[:, 0] * (y[1] - y[0])]))
        return boxes, pts

    def box_dist(cx, cy, w, h, pts):
        dx = np.maximum(np.abs(pts[:, 0] - cx) - w / 2, 0); dy = np.maximum(np.abs(pts[:, 1] - cy) - h / 2, 0)
        return np.sqrt(dx ** 2 + dy ** 2).min()

    def place_change(ax, k):
        """each change P at the free position nearest its own line's midpoint (tiers strict -> relaxed); False if none fits"""
        a = ANN[k]; fig.canvas.draw()
        abb = ax.get_window_extent(rend); boxes, pts = data_boxes(ax, k)
        inv = ax.transData.inverted(); WH = {}
        for g in ['rem', 'non']: a[g].set_position((0.5, 0))
        fig.canvas.draw()
        for g in ['rem', 'non']: e0 = bb(a[g]); WH[g] = (e0.width, e0.height)
        Gd = np.arange(-300, 301, 3.0); GX, GY = np.meshgrid(Gd, Gd); D2 = (GX ** 2 + GY ** 2).ravel(); o_ = np.argsort(D2, kind='stable')
        GX, GY, D2 = GX.ravel()[o_], GY.ravel()[o_], D2[o_]
        def first_free(g, placed):
            w, h = WH[g]; own, other = pts[g], pts['non' if g == 'rem' else 'rem']; mx, my = own[len(own) // 2]
            for near_, sep_, gap_ in TIERS:
                for dx_, dy_, d2 in zip(GX, GY, D2):
                    cx, cy = mx + dx_, my + dy_
                    bx = mtrans.Bbox.from_extents(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
                    if bx.x0 < abb.x0 + 4 or bx.x1 > abb.x1 - 2 or bx.y0 < abb.y0 + 2 or bx.y1 > abb.y1 - 2: continue
                    do = box_dist(cx, cy, w, h, own)
                    if do < CLR or do > near_: continue
                    if any(ov(bx, b_, pad=CLR) for b_ in boxes): continue
                    if any(ov(bx, b_, pad=max(gap_, do + sep_)) for b_ in placed): continue
                    if box_dist(cx, cy, w, h, other) < do + sep_: continue
                    return cx, cy, bx, np.sqrt(d2) + 1e4 * TIERS.index((near_, sep_, gap_))
            return None
        best = None
        for order in [('rem', 'non'), ('non', 'rem')]:
            res, placed = {}, []
            for g in order:
                f_ = first_free(g, placed)
                if f_ is None: break
                res[g] = f_; placed.append(f_[2])
            if len(res) == 2 and (best is None or res['rem'][3] + res['non'][3] < best['rem'][3] + best['non'][3]): best = res
        if best is None: return False
        for g in ['rem', 'non']: a[g].set_position(tuple(inv.transform((best[g][0], best[g][1]))))
        return True

    for ax, (k, *_) in zip(axs, PANELS):
        lo_, hi_ = ax._lohi; rng = hi_ - lo_
        for it in range(200):
            ax.set_ylim(lo_ - ax._pad[0] * rng, hi_ + ax._pad[1] * rng); fig.canvas.draw()
            grow = not place_change(ax, k)
            if grow: ax._pad[0] += 0.03; ax._pad[1] += 0.03
            else: break
        assert not grow, ('y-limits did not converge', k)
        print(f'  {k[:30]:30s}: y padding {ax._pad[0]:.2f} / {ax._pad[1]:.2f} of the data range', flush=True)

    # group labels + light brackets above the titles of each row
    fig.canvas.draw()
    GL, GLB = [], []
    for idx, lab in [((0, 0), 'Suppression State'), ((1, 4), 'Activation Programs')]:
        a0, a1 = axs[idx[0]], axs[idx[1]]
        ttl = max(bb(ANN[k]['title']).transformed(fig.transFigure.inverted()).y1 for ax_, (k, *_) in zip(axs, PANELS) if ax_ in axs[idx[0]:idx[1] + 1])
        b0, b1 = a0.get_position(), a1.get_position(); y = ttl + fy(0.22); x0, x1 = b0.x0 - 0.68 / WIN, b1.x1
        GLB.append(fig.add_artist(Line2D([x0, x0, x1, x1], [y - fy(0.13), y, y, y - fy(0.13)], transform=fig.transFigure, color=GR, lw=1.2)))
        GL.append(fig.text((x0 + x1) / 2, y + fy(0.05), lab, ha='center', va='bottom', fontsize=16, color=K))

    # colour key + one-line P key, right of the Suppression State panel
    HD = [Line2D([], [], color=BLUE, lw=2.0, marker='o', ms=6.5, label=f"Remission (n = {NP[('rem', 'Pre')]} patients)"),
          Line2D([], [], color=RED, lw=2.0, marker='o', ms=6.5,
                 label=f"Non-remission (n = {NP[('non', 'Pre')]} pre-, {NP[('non', 'Post')]} post-treatment)")]
    kx = axs[0].get_position().x1 + 0.75 / WIN
    ktop = axs[0].get_position().y1 + fy(0.05)
    LG = fig.legend(handles=HD, loc='upper left', bbox_to_anchor=(kx, ktop), ncol=1, frameon=False, fontsize=14, handlelength=1.8,
                    borderaxespad=0, labelspacing=0.6)
    fig.canvas.draw()
    lgb = LG.get_window_extent(rend).transformed(fig.transFigure.inverted())
    PKEY = 'Blue / red: change after treatment within group (biopsy-level mixed model)'
    KT = fig.text(kx + 0.004, lgb.y0 - fy(0.30), PKEY, ha='left', va='top', fontsize=13, color=K)

    # ---------------------------------------------------------------- layout and value checks
    fig.canvas.draw()
    for ax, (k, *_) in zip(axs, PANELS):
        a = ANN[k]; r = P.loc[k]
        txt = {nm: bb(a[nm]) for nm in ['rem', 'non', 'title']}
        names = list(txt)
        for i in range(len(names)):
            for j in range(i + 1, len(names)): assert not ov(txt[names[i]], txt[names[j]]), (k, names[i], names[j])
        boxes, pts = data_boxes(ax, k); allpts = np.vstack(list(pts.values()))
        for nm in ['rem', 'non', 'title']:
            for b_ in boxes: assert not ov(txt[nm], b_, pad=2.0), (k, nm, 'overlaps a point/error bar')
            e = txt[nm]; inside = (allpts[:, 0] >= e.x0 - 2) & (allpts[:, 0] <= e.x1 + 2) & (allpts[:, 1] >= e.y0 - 2) & (allpts[:, 1] <= e.y1 + 2)
            assert not inside.any(), (k, nm, 'overlaps a line')
        abb = ax.get_window_extent(rend)
        for nm in ['rem', 'non']:
            e = txt[nm]; assert e.x0 >= abb.x0 + 2 and e.x1 <= abb.x1 and e.y0 >= abb.y0 and e.y1 <= abb.y1, (k, nm, 'outside panel')
        assert txt['title'].y0 >= abb.y1, (k, 'title must sit above the panel')
        for g in ['rem', 'non']:
            e = txt[g]; cx, cy = (e.x0 + e.x1) / 2, (e.y0 + e.y1) / 2
            own, oth = pts[g], pts['non' if g == 'rem' else 'rem']
            do, dt = box_dist(cx, cy, e.width, e.height, own), box_dist(cx, cy, e.width, e.height, oth)
            assert CLR - 0.5 <= do <= NEAR + 0.5 and dt > do, (k, g, do, dt)
            mx, my = own[len(own) // 2]; dm = np.hypot(cx - mx, cy - my)
            assert dm < 0.45 * abb.width, (k, g, 'change P far from its line midpoint', dm)
            print(f'  {k[:30]:30s} {g}: {a[g].get_text():10s} {do:5.1f} px from own line, {dt:6.1f} px from other line, centre {dm:5.0f} px from midpoint')
        # every P printed equals the CSV
        csvp = {g: float(ST[(ST.feature_key == k) & (ST.contrast == CON[g])].p.iloc[0]) for g in ['rem', 'non']}
        assert [a['rem'].get_text(), a['non'].get_text()] == [fmt_p(csvp['rem']), fmt_p(csvp['non'])], k
        assert sum(isinstance(c_, matplotlib.text.Text) and c_.get_text().startswith('P ') for c_ in ax.get_children()) == 2, (k, 'only the 2 change P drawn')
        assert a['rem'].get_color() == BLUE and a['non'].get_color() == RED
        tl = [t.get_window_extent(rend) for t in ax.get_xticklabels()]; assert not ov(tl[0], tl[1], pad=4), (k, 'tick labels overlap')
    tb = [a.get_tightbbox(rend) for a in axs]
    for i, j in [(1, 2), (2, 3), (3, 4)]: assert tb[i].x1 < tb[j].x0, ('neighbouring panels overlap', i, j)
    titles = [bb(ANN[k]['title']) for k in PK]
    for i, j in [(1, 2), (2, 3), (3, 4)]: assert not ov(titles[i], titles[j], pad=6), ('titles overlap', i, j)
    for j in range(1, 5): assert tb[0].y0 > tb[j].y1, ('rows overlap', j)
    keys = [LG.get_window_extent(rend), bb(KT)]
    for b_ in tb:
        for kb in keys: assert not ov(b_, kb), 'key overlaps a panel'
    assert not ov(keys[0], keys[1]), 'key items overlap'
    for gl, gb in zip(GL, GLB):
        for b_ in tb: assert not ov(bb(gl), b_) and not ov(gb.get_window_extent(rend), b_, pad=0), 'group label/bracket overlaps a panel'
        for kb in keys: assert not ov(bb(gl), kb) and not ov(gb.get_window_extent(rend), kb, pad=0), 'group label/bracket overlaps key'
    fbb = fig.bbox
    for e in keys + [bb(g) for g in GL] + tb: assert e.x0 >= fbb.x0 - 1 and e.x1 <= fbb.x1 + 1 and e.y0 >= fbb.y0 - 1 and e.y1 <= fbb.y1 + 1, ('outside figure', e, fbb)
    alltxt = ' '.join(t.get_text() for t in fig.findobj(matplotlib.text.Text))
    for w in ['Nominal', 'nominal', 'D6', 'quiescent', 'territory', 'blunted', 'robust', 'colocali', 'Difference', 'bracket', 'TNFRSF9']: assert w not in alltxt, w
    assert axs[4].get_ylabel() == 'Log-normalised expression' and ANN['gene_TNFRSF4']['title'].get_text() == 'TNFRSF4 (OX40)'
    print('layout checks passed (annotations vs points/error bars/lines, panels, rows, keys, group labels; P values = CSV)', flush=True)
    fig.savefig(f'{STEM}.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    print('figure size (in):', [round(v / DPI, 2) for v in fig.get_tightbbox(rend).size * DPI], flush=True); plt.close(fig)
    print('wrote', f'{STEM}.png', flush=True)

    # ---------------------------------------------------------------- P table (printed for the report)
    NAME = {'supp_score': 'Treg Suppression Score', IL2: 'IL-2/STAT5', TNFR2: 'TNFR2', 'costim_module': 'receptor module', 'gene_TNFRSF4': 'TNFRSF4',
            'pseudotime': 'Treg Pseudotime'}
    print('\nP table (drawn: change within group; not drawn: difference in change, pre-treatment difference; q = BH over the 4 activation programs)')
    for k in PK:
        r = P.loc[k]; qq = lambda c: (', ' + fmt_p(r[f'{c}_q'], 'q')) if k in AK else ''
        print(f"  {NAME[k]:24s} rem {fmt_p(r.rem_p):10s}{qq('rem'):14s} | non {fmt_p(r.non_p):10s}{qq('non'):14s} | "
          f"diff {fmt_p(r.diff_p):10s}{qq('diff'):14s} | pre {fmt_p(r.base_p):10s}{qq('base')}")

    # ---------------------------------------------------------------- legend
    pf = lambda k, c: fmt_p(P.loc[k, f'{c}_p'])
    pq = lambda k, c: f"{fmt_p(P.loc[k, f'{c}_p'])}, {fmt_p(P.loc[k, f'{c}_q'], 'q')}"
    assert (P.loc[['supp_score', 'pseudotime'], ['rem_p', 'non_p']] >= 0.05).all().all()
    assert (P.loc[AP3, 'rem_p'] < 0.05).all() and (P.loc[AP3, 'rem_est'] < 0).all() and (P.loc[AP3, 'non_p'] >= 0.05).all()
    assert (P.loc[AP3, 'rem_q'] < 0.05).all()
    X4 = 'gene_TNFRSF4'
    REC = 'CD27, CD28, CTLA4, ICOS, IL10RA, IL21R, IL2RA, IL2RB, IL2RG, IL4R, IL7R, ITGB1, ITGB2, PDCD1, TNFRSF18, TNFRSF4'   # the 16 module receptors of cell_level_lmm
    assert len(REC.split(', ')) == 16
    NGI = {k: (int(G[G.feature_key == k].n_genes_scored.iloc[0]), int(G[G.feature_key == k].n_genes_in_set.iloc[0])) for k in [IL2, TNFR2]}
    QMAX = fmt_p(P.loc[AP3, 'rem_q'].max(), 'q')
    PCT4 = float(S6.pct_tregs_expressing_all.iloc[0])
    # rationale for TNFRSF4, asserted from the Fig 4F and Supplementary Fig 8F source tables
    F3F = pd.read_csv(f'{R}/Fig3D-I_treg_pseudotime_signals/fig3F_marker_dotplot.csv').drop_duplicates('names').sort_values('log2FC_bin0', ascending=False)
    assert F3F.names.iloc[:2].tolist() == ['TNFRSF18', 'TNFRSF4'], F3F.names.iloc[:3].tolist()
    DRV = pd.read_csv(f'{R}/SuppFig8A-I_treg_pseudotime_robustness/ext8_tf_activity_along_pseudotime_drivers.csv')
    DRV = DRV[(DRV.tf == 'YY1') & DRV.is_driver].sort_values('contribution'); assert DRV.gene.iloc[0] == 'TNFRSF4' and DRV.contribution.iloc[0] < 0
    def chg(k, g):   # within-group change of TNFRSF4, worded from the model estimate and P
        p_ = P.loc[k, f'{g}_p']
        if p_ >= 0.05: return f"did not change in {'Remission' if g == 'rem' else 'Non-remission'} ({pq(k, g)})"
        return f"{'fell' if P.loc[k, f'{g}_est'] < 0 else 'rose'} in {'Remission' if g == 'rem' else 'Non-remission'} ({pq(k, g)})"
    def base_sentence():
        sig = [k for k in PK if P.loc[k, 'base_p'] < 0.05]; ns = [k for k in PK if k not in sig]
        one = lambda k: f"{NAME[k]} {pq(k, 'base') if k in AK else pf(k, 'base')}"
        s = '; '.join(f"{NAME[k]} was {'lower' if P.loc[k, 'base_est'] < 0 else 'higher'} in Remission than in Non-remission "
                  f"({pq(k, 'base') if k in AK else pf(k, 'base')})" for k in sig)
        return (f"Before treatment, {s}; " if sig else 'Before treatment, ') + \
               ('the other panels did not differ (' if sig else 'Remission and Non-remission did not differ in any panel (') + '; '.join(one(k) for k in ns) + ').'
    txt = (f"Treg suppression state and Treg activation programs before and after anti-TNF, by remission. TAURUS Crohn's disease scRNA-seq: "
       f"{int(first.n_tregs):,} Tregs in {int(first.n_biopsies)} biopsies (each with ≥ {MIN_TREG} Tregs) from {int(first.n_patients)} patients "
       f"(Remission: {NP[('rem', 'Pre')]} patients, {NB[('rem', 'Pre')]} pre-treatment and {NB[('rem', 'Post')]} post-treatment biopsies; "
       f"Non-remission: {NP[('non', 'Pre')]} patients pre-treatment ({NB[('non', 'Pre')]} biopsies) and {NP[('non', 'Post')]} post-treatment "
       f"({NB[('non', 'Post')]} biopsies)). Top row, Suppression State: Treg Suppression Score (early minus late Treg pseudotime gene-set "
       f"score; higher values, Activated, more suppressive Tregs; lower values, Low-suppressive Tregs). Bottom row, Activation Programs: "
       f"IL-2/STAT5 signaling and TNFR2 non-canonical NF-κB, the activation pathways of the Treg Activation & Suppressive Function block in "
       f"Fig 4G (IL-2/STAT5 Signaling, Hallmark, {NGI[IL2][0]} of {NGI[IL2][1]} genes scored; TNFR2 Non-Canonical NF-κB, Reactome, full gene "
       f"set as in Fig 4G, {NGI[TNFR2][0]} of {NGI[TNFR2][1]} genes scored), which fall along the Treg pseudotime, each scored per Treg by "
       f"univariate linear model (ULM, decoupler) on log-normalised expression; the Co-Stimulatory/Suppressive Receptor Module, the "
       f"receptors of the Fig 5C co-stimulatory signals (16 receptors from CellChatDB co-stimulatory/checkpoint pathways, expressed in ≥ 5% "
       f"of Tregs: {REC}), scored per Treg with scanpy score_genes; and TNFRSF4 (OX40), per-Treg log-normalised expression (detected in "
       f"{PCT4:.0f}% of Tregs). TNFRSF4 is, with TNFRSF18 (GITR), the top early-bin (Activated) Treg gene in Fig 4F, and the top driver "
       f"target of YY1, the top falling TF in Supplementary Fig 8F; OX40L→OX40 is one of the Fig 5C co-stimulatory pairs. Each feature was "
       f"averaged over each biopsy's Tregs. Biopsy means were modelled as feature ~ timepoint × remission + (1 | patient) + "
       f"(1 | patient:timepoint), with patient and patient-by-timepoint random effects and every biopsy weighted equally (lme4/lmerTest REML, "
       f"Satterthwaite df); all P are linear contrasts of this one model, fitted per feature; q, Benjamini–Hochberg over the four activation "
       f"programs, per contrast. Points, mean of the biopsy means of each group before (Pre-treatment) and after (Post-treatment) treatment; "
       f"error bars, SEM across biopsies; blue, Remission; red, Non-remission. P shown alongside each line, in the line's colour, are for the "
       f"change after treatment (post- vs pre-treatment) within that group. The Treg Suppression Score did not change after treatment in "
       f"either group (Remission {pf('supp_score', 'rem')}; Non-remission {pf('supp_score', 'non')}). IL-2/STAT5, TNFR2 and the receptor "
       f"module fell in Remission (IL-2/STAT5 {pf(IL2, 'rem')}; TNFR2 {pf(TNFR2, 'rem')}; receptor module {pf('costim_module', 'rem')}; "
       f"all {QMAX}) and did not change in Non-remission (P = {sig2(P.loc[AP3, 'non_p'].min())} to {sig2(P.loc[AP3, 'non_p'].max())}); "
       f"TNFRSF4 {chg(X4, 'rem')} and {chg(X4, 'non')}. Contrasts from the same model (not drawn): the difference in change between groups "
       f"(change in Remission vs change in Non-remission) was Treg Suppression Score {pf('supp_score', 'diff')}; "
           + '; '.join(f"{NAME[k]} {pq(k, 'diff')}" for k in AK) + f". {base_sentence()} After treatment, TNFRSF4 was "
       f"{'lower' if P.loc[X4, 'post_est'] < 0 else 'higher'} in Remission than in Non-remission ({pf(X4, 'post')}). Y-axes differ between panels. "
       f"Source data note (not shown): Treg Pseudotime (diffusion pseudotime), analysed with the same model, did not change after treatment "
       f"in either group (Remission {pf('pseudotime', 'rem')}; Non-remission {pf('pseudotime', 'non')}); Remission vs Non-remission before "
       f"treatment {pf('pseudotime', 'base')}, after treatment {pf('pseudotime', 'post')}; difference in change {pf('pseudotime', 'diff')}.")
    assert 'the activation pathways of the Treg Activation & Suppressive Function block in Fig 4G' in txt
    assert txt.lower().count('supplementary') == 1 and 'Supplementary Fig 8F' in txt, 'only the Supplementary Fig 8F citation'
    assert all(w not in txt.lower() for w in ['supplementary/', 'nominal', 'd6', 'quiescent', 'territory', 'blunted', 'robust', 'colocali', 'tnfrsf9']), 'legend wording'
    with open(f'{STEM}_legend.txt', 'w') as fh: fh.write(txt + '\n')
    print(txt)
