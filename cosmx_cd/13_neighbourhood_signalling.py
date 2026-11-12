#!/usr/bin/env python
"""CD CosMx 6K resections: ligand-receptor signalling per cellular neighbourhood (CN).

Ligand-receptor (LR) score per cell = LIANA bivariate local score (distance-weighted ligand-receptor co-expression, Gaussian kernel, 200 px bandwidth,
NicheNet prior-network pairs), grouped into CellChatDB-based modules and pre-specified inflammatory axes and summarised per CN and section.

Panels made
  Fig 5C      Treg Suppression Score per CN (left; forest table from 12_neighbourhood_treg_score.py) with CN signalling (right), one figure
  Supp 9A     CN characterisation: cell-type log2 enrichment, LR modules and inflammatory axes per CN
  Supp 9B     signalling of each CN vs the other CNs (linear mixed model, 95% CI)
  Supp Table 10 (LR module membership)
  Also writes the Treg receiver pair table and per-Treg module scores read by 17_neighbourhood_apc_mhcii.py.

Inputs
  /share/fsmresfiles/UC/AtoMx/CD_surgical_resection/ (per-section Processed/norm_log1p.h5ad; combined/all_samples.parquet; combined/spot_lr_df.parquet from
    07_resection_liana_merge.ipynb; combined/raw_counts_strict_tregs.csv from 05_resection_treg_counts.ipynb)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/cellchatdb_human_interactions_snapshot.csv (CellChatDB human interactions);
    /share/fsmresfiles/ydn4687/nichenet/lr_network_human_21122021.rds (NicheNet prior network)
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/: D_shared_inputs/treg_definition_comparison_per_cell.csv
    (10_treg_suppression_score.py); B_resection_analyses/cellular_neighborhoods_v2/ (cn_per_cell.parquet, cn_treg_state_enrichment_all.csv) and
    B_resection_analyses/ligand_receptor_state/treg_local_ligand_receptor.csv.gz (11_neighbourhoods_local_tables.py); A_fig4_panels/fig4_cn_composition_heatmap.csv
    and A_fig4_panels/fig4_cn_treg_state_forest.csv + _legend.txt (12_neighbourhood_treg_score.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/panels/fig4_and_related/A_fig4_panels/ (CN-level panels and tables),
    B_resection_analyses/lr_modules_cellchat/ and B_resection_analyses/inflammatory_axes/ (caches and tables), METHODS_inflammatory_axes.md
  /share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/ (Fig 5C) and /share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/ (Supp 9A-B)

Run order of sections
  cn_module_relative_scores -> treg_receiver_lr_modules -> cn_lr_modules_celltype_restricted -> cn_lr_modules_allcells -> cn10_lr_pairs
  -> lr_costim_inhibitory_modules -> cn_lr_modules_immunomodulatory -> cn_lr_modules_costim -> cn_inflammatory_axes -> S9A_cn_characterisation
  -> cn_signalling_table -> F5C_score_signalling -> S9B_cn_signalling
  cn_lr_modules_allcells writes the module membership and per-cell module cache used by the later module sections; cn10_lr_pairs and
  lr_costim_inhibitory_modules write the per-pair caches; S9A_cn_characterisation combines the cn_lr_modules_costim and cn_inflammatory_axes tables with the
  CN composition; cn_signalling_table and F5C_score_signalling read the S9A_cn_characterisation table; S9B_cn_signalling reads the cn_signalling_table table.
  Fig 5C is drawn with VM=1.2 WIDE=1 RUN_SECTIONS=F5C_score_signalling python <this file>.

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its upstream outputs exist. Command-line
  arguments and environment variables used by a section are listed in its own docstring.
  Paths are the lab's absolute paths (to be replaced by placeholders before public release).
  Environment: see the README section "Environment" (OPENBLAS_CORETYPE=SkylakeX on Intel Sapphire Rapids nodes; OPENBLAS_NUM_THREADS=32 for the
  per-cell sections; the environment asserts are kept as comments).
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
# Section cn_module_relative_scores: Relative LR module score per cell, summarised per CN (pooled, section-adjusted, section consistency)
# Panels: LR module inputs (relative-score reference table)
# ======================================================================================================================
if _run('cn_module_relative_scores'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/63_module_cn_adjusted.py'   # working folder of this step (relative paths resolve from here)
    """Relative ligand-receptor module score per cell, summarised per cellular neighbourhood (CN) three ways on the same cells:
  (a) pooled      mean relative module score per CN, all sections pooled;
  (b) adjusted    each cell's module score z-standardised within its section, then averaged per CN (sections differ in baseline signal and
                  contribute unevenly to each CN);
  (c) consistency for each module x CN, the number of sections in which the CN sits above that section's own median (sections contributing
                  >= 200 cells of the CN).
Cells: random 300,000 cells of the 11 resection sections (seed 0). Modules (CellChatDB pathway_name, objects/cellchatdb_human_interactions_snapshot.csv):
  Antigen presentation = MHC-II; Co-stimulation / suppressive Treg = CD80, CD86, PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2, IL10;
  Inflammatory cytokine = IL1, IL6, TNF, IFN-II, IL12, IL17, OSM; Matrix-integrin = all ECM-Receptor pathways. A LIANA pair of combined/spot_lr_df.parquet
  joins a module when a CellChatDB interaction of the module has the same ligand set and an overlapping receptor set (module_of()). Relative score = mean
  local score over the module's pairs divided by the mean over all pairs of the cell.
CN = CN_k10 of cellular_neighborhoods_v2/cn_per_cell.parquet (11_neighbourhoods_local_tables.py).
Outputs: B_resection_analyses/lr_modules_cellchat/module_scores_per_cell_sample.csv.gz and module_by_CN_pooled_vs_adjusted.{csv,png} (the CSV gives the CN
  row labels and the relative-score reference used by cn_lr_modules_celltype_restricted and cn_lr_modules_allcells)."""
    import os, gc, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/B_resection_analyses/lr_modules_cellchat'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; SPOT = f'{R}/combined/spot_lr_df.parquet'
    MOD = {'Antigen presentation': ['MHC-II'], 'Co-stimulation / suppressive Treg': ['CD80', 'CD86', 'PD-L1', 'PDL2', 'ICOS', 'OX40', 'GITRL', 'CD70', 'CD40', 'IL2', 'IL10'],
           'Inflammatory cytokine': ['IL1', 'IL6', 'TNF', 'IFN-II', 'IL12', 'IL17', 'OSM']}
    cc = pd.read_csv(f'{D}/objects/cellchatdb_human_interactions_snapshot.csv')
    MOD['Matrix-integrin'] = sorted(cc.loc[cc.annotation == 'ECM-Receptor', 'pathway_name'].unique())
    gs = lambda x: frozenset(y.strip() for y in str(x).replace(';', ',').split(',') if y.strip() and y.strip() != 'nan')
    cc['L'] = cc['ligand.symbol'].map(gs); cc['R'] = cc['receptor.symbol'].map(gs)
    cols = [c for c in pq.read_schema(SPOT).names if '^' in c]
    def module_of(c):
        l, r = c.split('^'); L, Rr = frozenset(l.split('_')), frozenset(r.split('_')); hits = []
        for m, pws in MOD.items():
            sub = cc[cc.pathway_name.isin(pws)]
            if any((L == a) and (Rr == b or Rr <= b or b <= Rr) for a, b in zip(sub.L, sub.R)): hits.append(m)
        return hits
    MEM = {c: module_of(c) for c in cols}; MEM = {c: m for c, m in MEM.items() if m}
    cn = pq.read_table(f'{P4}/B_resection_analyses/cellular_neighborhoods_v2/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
    NM = pd.read_csv(f'{P4}/B_resection_analyses/cellular_neighborhoods_v2/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN\u207aRSPO3\u207a Fib').replace('PLA2G2A Fib', 'PLA2G2A\u207a Fib') for k, v in NM.items()}
    idx_all = pq.read_table(SPOT, columns=['__index_level_0__']).column(0).to_pylist(); pos = pd.Series(np.arange(len(idx_all)), index=idx_all)
    rng = np.random.default_rng(0); pool = np.asarray(cn.index[cn.index.isin(pos.index)])
    samp = pd.Index(rng.choice(pool, 300000, replace=False)); take = np.sort(pos.loc[samp].values); rows_idx = np.asarray(idx_all)[take]
    tot = np.zeros(len(take)); acc = {m: np.zeros(len(take)) for m in MOD}; npair = {m: 0 for m in MOD}
    for i in range(0, len(cols), 250):
        ch = cols[i:i + 250]; A = pq.read_table(SPOT, columns=ch).take(take).to_pandas().astype(np.float32).fillna(0.0)
        tot += A.values.sum(1)
        for m in MOD:
            ps = [c for c in ch if c in MEM and m in MEM[c]]
            if ps: acc[m] += A[ps].values.sum(1); npair[m] += len(ps)
        del A; gc.collect()
    den = tot / len(cols); den[den == 0] = np.nan
    S = pd.DataFrame({m: (acc[m] / max(npair[m], 1)) / den for m in MOD}, index=rows_idx).fillna(0.0)
    S['CN'] = cn.CN_k10.reindex(S.index).values; S['section'] = cn.section.reindex(S.index).values
    S.to_csv(f'{OUT}/module_scores_per_cell_sample.csv.gz', float_format='%.4g')
    mods = list(MOD)
    pooled = S.groupby('CN')[mods].mean()
    Z = S.groupby('section')[mods].transform(lambda x: (x - x.mean()) / (x.std() or 1))
    adjusted = Z.assign(CN=S.CN.values).groupby('CN')[mods].mean()
    cons = pd.DataFrame(0, index=sorted(S.CN.unique()), columns=mods)
    ntest = pd.Series(0, index=sorted(S.CN.unique()))
    for sec, g in S.groupby('section'):
        med = g[mods].median()
        for c, gg in g.groupby('CN'):
            if len(gg) < 200: continue
            ntest[c] += 1
            for m in mods: cons.loc[c, m] += int(gg[m].mean() > med[m])
    pooled.index = [NM[i] for i in pooled.index]; adjusted.index = [NM[i] for i in adjusted.index]; cons.index = [NM[i] for i in cons.index]; ntest.index = cons.index
    out = pd.concat({'pooled_mean': pooled, 'section_adjusted_z': adjusted, 'sections_above_median': cons}, axis=1)
    out['n_sections_tested'] = ntest; out.to_csv(f'{OUT}/module_by_CN_pooled_vs_adjusted.csv', float_format='%.3g')
    plt.rcParams.update({'font.size': 14})
    fig, axes = plt.subplots(1, 2, figsize=(17, 5.4), gridspec_kw={'wspace': 0.35})
    Zp = (pooled - pooled.mean()) / pooled.std()
    ann = np.array([[f'{pooled.iloc[r, c]:.2f}' for c in range(len(mods))] for r in range(len(pooled))])
    sns.heatmap(Zp, cmap='RdBu_r', center=0, vmin=-2, vmax=2, annot=ann, fmt='', ax=axes[0], linewidths=0.4, cbar_kws={'label': 'z across CN', 'shrink': 0.8})
    axes[0].set_title('(a) Pooled across sections', fontsize=15, loc='left')
    ann2 = np.array([[f'{adjusted.iloc[r, c]:+.2f}\n{int(cons.iloc[r, c])}/{int(ntest.iloc[r])}' for c in range(len(mods))] for r in range(len(adjusted))])
    sns.heatmap(adjusted, cmap='RdBu_r', center=0, vmin=-0.5, vmax=0.5, annot=ann2, fmt='', annot_kws={'fontsize': 10}, ax=axes[1], linewidths=0.4, cbar_kws={'label': 'mean within-section z', 'shrink': 0.8})
    axes[1].set_title('(b) Section-adjusted  \u2014  value and sections above median', fontsize=15, loc='left')
    for ax in axes:
        ax.set_xticklabels(['Antigen\npresentation', 'Co-stim /\nsuppr. Treg', 'Inflammatory\ncytokine', 'Matrix-\nintegrin'], rotation=0, fontsize=12)
        ax.set_ylabel(''); ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=11)
    fig.savefig(f'{OUT}/module_by_CN_pooled_vs_adjusted.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    pd.set_option('display.width', 220); print('POOLED'); print(pooled.round(3).to_string()); print('ADJUSTED'); print(adjusted.round(3).to_string()); print('SECTIONS ABOVE MEDIAN (of n tested)'); print(cons.assign(n=ntest).to_string())


# ======================================================================================================================
# Section treg_receiver_lr_modules: LR modules delivered to strict Tregs from the canonical sender cell types, per CN
# Panels: LR module inputs (Treg receiver pairs and per-Treg module scores, read by cn_lr_modules_celltype_restricted and 17_neighbourhood_apc_mhcii.py)
# ======================================================================================================================
if _run('treg_receiver_lr_modules'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/74_lr_modules_celltype_restricted.py'   # working folder of this step (relative paths resolve from here)
    """CellChatDB-labelled ligand-receptor modules delivered to strict Tregs (receivers) from the main canonical sender cell types of each signal, by
cellular neighbourhood, with the section-to-section spread.
Module labels (CellChatDB.human pathway_name, objects/cellchatdb_human_interactions_snapshot.csv):
  Antigen presentation = MHC-II; Co-stimulation / suppressive Treg = CD80, CD86, PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2, IL10;
  Inflammatory cytokine = IL1, IL6, TNF, IFN-II, IL12, IL17, OSM; Matrix-integrin = all 'ECM-Receptor' pathways.
Pairs: the NicheNet prior network (/share/fsmresfiles/ydn4687/nichenet/lr_network_human_21122021.rds; the pairs of the resection LIANA run), assigned to
  modules with CellChatDB (pair in module if a CellChatDB interaction of the module's pathways has the NicheNet ligand among its ligand genes and the
  NicheNet receptor among its receptor genes; unmatched pairs dropped); ligand and receptor both on the CosMx 6K panel. NicheNet pairs are single genes;
  expression z within section over all cells. pair_source_nichenet_vs_cellchatdb.csv gives the pair counts of both pair sources per module.
Senders (cell_type_general of all_samples.parquet; Tregs are never senders), by CellChatDB pathway:
  Antigen presentation  MHC-II                                          -> DC, B cell, CD9+ M2 Mac, Mo-Mac
  Co-stimulation        CD80, CD86, ICOS(L), PD-L1, PDL2, CD70, OX40, GITRL (APC ligands) -> DC, B cell
                        IL2 (IL2-family), CD40 (CD40LG)                 -> CD4 Effector, CD4 Trm, Tfh
                        IL10 (IL10-family)                              -> DC, CD9+ M2 Mac, Mo-Mac
  Inflammatory cytokine TNF, IL1 (incl. IL18), IL6, IL12 (incl. IL23, IL27), OSM -> Mo-Mac, CD9+ M2 Mac, DC, Neutrophil
                        IFN-II (IFNG), IL17                             -> CD4 Effector, CD8 Trm, NK, Gd T
  Matrix-integrin       all ECM-Receptor                                -> OGN+RSPO3+ Fib, ADAMDEC1+ Fib, Activated Fib, Crypt Top Fib,
                                                                           PLA2G2A+ ECM Fib, Myofibroblast
Per strict Treg (6,554 Tregs, 11 CD resection sections), per pair:
  L = mean within-section z of the ligand over the pair's own sender group's cells within 50 um (pair counts if >= 1 such cell);
  R = the Treg's own within-section z of the receptor;
  pair score = max(L, 0) x max(R, 0) (positive part; raw L.R kept as a sensitivity metric).
  Module score = mean over the pairs that count; a Treg is scored for a module if >= 1 pair counts. Sensitivity (CSV only): the same per-pair rule with
  >= 3 sender cells per pair (metrics 'L min3', 'LR min3', 'LxR raw min3').
Metrics per module: 'L' = ligand supply = mean L over the module's unique (ligand, sender set) units; 'LR' = mean pair score over the module's pairs;
  'LxR raw' = mean L.R; 'Sender density' = log1p number of the module's sender cells within 50 um.
Neighbourhood = the Treg's own CN_k10 (cellular_neighborhoods_v2/cn_per_cell.parquet; CN names from cn_treg_state_enrichment_all.csv, compact labels).
Stats, per metric x module x CN: score z-standardised within section over the scored Tregs; sections with >= 10 scored Tregs in the CN (and >= 10
  outside it); per section the mean in the CN minus the mean of the other Tregs (the dots); linear mixed models over all scored Tregs of those sections,
  z ~ in_CN + library size, with (ri) a section random intercept (figure model) and (rs) a random intercept and a random in_CN slope by section
  (sensitivity, CSV only). Library size = log total counts per Treg from combined/raw_counts_strict_tregs.csv, z within section. 95% Wald CI; BH over
  the CNs within each module (separately per metric and model); sections agreeing = sections whose difference has the sign of the estimate;
  section_ttest_p = one-sample t-test of the per-section differences. CN6 and CN9 (>= 10 Tregs in < 3 sections) are not drawn.
Outputs: A_fig4_panels/supp_lr_modules_treg_ligand_supply.{png,csv} + _legend.txt (ligand supply; the CSV holds all metrics and both models);
  A_fig4_panels/ext_lr_modules_treg_ligand_x_receptor.png + _legend.txt (ligand x receptor);
  B_resection_analyses/lr_modules_cellchat/: treg_receiver_pairs_kept.csv, pair_source_nichenet_vs_cellchatdb.csv,
  treg_receiver_module_scores_per_treg.csv.gz (per-Treg table), treg_receiver_modules_sender_density.png.
Environment variables: SKIP_COMPUTE=1 re-plots from the per-Treg table; FIG_MODEL=ri|rs selects the plotted model (default ri);
  TOTAL_SUPPLY=1 adds metric 'T' = total ligand supply per Treg (per pair, the sum over the pair's sender cells within 50 um of the ligand's
  log1p-normalised expression, 0 if no sender, so no Treg is dropped; module score = mean over the module's pairs; the per-sender mean 'L' becomes a
  sensitivity metric) and writes only the *_total files: A_fig4_panels/supp_lr_modules_treg_ligand_supply_total.{png,csv} + _legend.txt and
  B_resection_analyses/lr_modules_cellchat/treg_receiver_module_scores_per_treg_total.csv.gz."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
    from anndata.io import read_elem
    from scipy import sparse
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    STEM = f'{OUTA}/supp_lr_modules_treg_ligand_supply'; STEM_EXT = f'{OUTA}/ext_lr_modules_treg_ligand_x_receptor'
    TOTAL = os.environ.get('TOTAL_SUPPLY') == '1'; SFX = '_total' if TOTAL else ''
    if TOTAL: STEM = f'{STEM}_total'
    CACHE = f'{OUTB}/treg_receiver_module_scores_per_treg{SFX}.csv.gz'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50; MIN_SEND = 1; MIN_SEND_SENS = 3; MIN_TREG_CN = 10
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    # ---------------- modules and sender sets (main canonical senders per pathway)
    MODS = ['Antigen presentation', 'Co-stimulation / suppressive Treg', 'Inflammatory cytokine', 'Matrix-integrin']
    MOD = {'Antigen presentation': ['MHC-II'], 'Co-stimulation / suppressive Treg': ['CD80', 'CD86', 'PD-L1', 'PDL2', 'ICOS', 'OX40', 'GITRL', 'CD70', 'CD40', 'IL2', 'IL10'],
           'Inflammatory cytokine': ['IL1', 'IL6', 'TNF', 'IFN-II', 'IL12', 'IL17', 'OSM']}
    cc = pd.read_csv(f'{D}/objects/cellchatdb_human_interactions_snapshot.csv')
    MOD['Matrix-integrin'] = sorted(cc.loc[cc.annotation == 'ECM-Receptor', 'pathway_name'].unique())
    SETS = {'MHCII_APC': ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac'],
            'DC_B': ['DC', 'B cell'],
            'CD4_T': ['CD4 Effector', 'CD4 Trm', 'Tfh'],
            'IL10_myeloid': ['DC', 'CD9+ M2 Mac', 'Mo-Mac'],
            'Inflam_myeloid': ['Mo-Mac', 'CD9+ M2 Mac', 'DC', 'Neutrophil'],
            'Inflam_T_NK': ['CD4 Effector', 'CD8 Trm', 'NK', 'Gd T'],
            'Fibroblast': ['OGN+RSPO3+ Fib', 'ADAMDEC1+ Fib', 'Activated Fib', 'Crypt Top Fib', 'PLA2G2A+ ECM Fib', 'Myofibroblast']}
    COSTIM_APC = ('CD80', 'CD86', 'ICOS', 'PD-L1', 'PDL2', 'CD70', 'OX40', 'GITRL')
    def sender_set(module, pathway):
        if module == 'Antigen presentation': assert pathway == 'MHC-II'; return 'MHCII_APC'
        if module == 'Co-stimulation / suppressive Treg':
            if pathway in ('IL2', 'CD40'): return 'CD4_T'
            if pathway == 'IL10': return 'IL10_myeloid'
            assert pathway in COSTIM_APC, pathway; return 'DC_B'
        if module == 'Inflammatory cytokine':
            if pathway in ('IFN-II', 'IL17'): return 'Inflam_T_NK'
            assert pathway in ('IL1', 'IL6', 'TNF', 'IL12', 'OSM'), pathway; return 'Inflam_myeloid'
        return 'Fibroblast'
    gs = lambda s: tuple(sorted(x.strip() for x in str(s).replace(';', ',').split(',') if x.strip() and x.strip() != 'nan'))
    with h5py.File(f'{R}/{FOLDER["43563_A8"]}/Processed/norm_log1p.h5ad', 'r') as h: PANEL = set(read_elem(h['var']).index.astype(str))
    # Pair source = NicheNet prior network (the LR pairs of the resection LIANA run, 06_resection_liana.ipynb);
    # CellChatDB is used only to assign pairs to modules (complexes split on '_'): a NicheNet pair (ligand l, receptor r)
    # belongs to module m if a CellChatDB interaction of m's pathways has l among its ligand gene(s) and a receptor gene set that contains r.
    # (Requiring an equal ligand SET would drop every single-gene NicheNet ligand whose CellChatDB ligand is a complex - all IL12/IL23/IL27
    # pairs - so the ligand side is matched by gene membership.)
    import pyreadr
    NN = list(pyreadr.read_r('/share/fsmresfiles/ydn4687/nichenet/lr_network_human_21122021.rds').values())[0][['from', 'to']].drop_duplicates()
    ccL = cc['ligand.symbol'].map(lambda x: frozenset(gs(x))); ccR = cc['receptor.symbol'].map(lambda x: frozenset(gs(x)))
    PR = []
    for m in MODS:
        sub = cc.pathway_name.isin(MOD[m])
        for l, r in zip(NN['from'], NN['to']):
            hit = sorted({pw for pw, a, b in zip(cc.pathway_name[sub], ccL[sub], ccR[sub]) if l in a and r in b})
            if not hit: continue
            ss = {sender_set(m, pw) for pw in hit}; assert len(ss) == 1, (m, l, r, hit)
            miss = [g for g in (l, r) if g not in PANEL]
            PR.append(dict(module=m, pathway='/'.join(hit), interaction=f'{l}_{r}', ligand=l, receptor=r, sender_set=ss.pop(), kept=not miss,
                           missing_genes=','.join(miss), pair_source='NicheNet v2 lr_network_human_21122021'))
    PR = pd.DataFrame(PR); PR['sender_cell_types'] = PR.sender_set.map(lambda s: ', '.join(SETS[s]))
    dup = PR.duplicated(['module', 'ligand', 'receptor']); assert not dup.any(), PR[dup]
    if not TOTAL: PR.to_csv(f'{OUTB}/treg_receiver_pairs_kept.csv', index=False)
    # pair counts per module: CellChatDB-direct pairs vs NicheNet pairs
    CCD = []
    for m in MODS:
        for _, r in cc[cc.pathway_name.isin(MOD[m])].iterrows():
            L_, R_ = gs(r['ligand.symbol']), gs(r['receptor.symbol'])
            CCD.append(dict(module=m, ligand='_'.join(L_), receptor='_'.join(R_), on_panel=all(g in PANEL for g in L_ + R_)))
    CCD = pd.DataFrame(CCD).drop_duplicates(['module', 'ligand', 'receptor'])
    PCMP = pd.DataFrame({'n_pairs_cellchatdb_direct_on_panel': CCD[CCD.on_panel].groupby('module').size(),
                         'n_ligands_cellchatdb_direct_on_panel': CCD[CCD.on_panel].groupby('module').ligand.nunique(),
                         'n_pairs_nichenet_matched': PR.groupby('module').size(), 'n_pairs_nichenet_on_panel': PR[PR.kept].groupby('module').size(),
                         'n_ligands_nichenet_on_panel': PR[PR.kept].groupby('module').ligand.nunique(),
                         'n_pairs_nichenet_exact_ligand_set_rule_of_24': PR.assign(e=[any(frozenset([l]) == a and r in b for a, b in zip(ccL[cc.pathway_name.isin(MOD[m])], ccR[cc.pathway_name.isin(MOD[m])])) for m, l, r in zip(PR.module, PR.ligand, PR.receptor)])[lambda x: x.kept & x.e].groupby('module').size()}).reindex(MODS).fillna(0).astype(int)
    PCMP.index.name = 'module'
    if not TOTAL: PCMP.to_csv(f'{OUTB}/pair_source_nichenet_vs_cellchatdb.csv')
    print(PCMP.to_string(), flush=True)
    KP = PR[PR.kept].reset_index(drop=True)
    print(KP.groupby('module').size().to_string(), flush=True)
    LUN = sorted(set(KP.ligand)); RUN = sorted(set(KP.receptor)); GENES = sorted({g for u in LUN + RUN for g in u.split('_')})
    LSET = KP[['ligand', 'sender_set']].drop_duplicates()
    NEED = {s: sorted(LSET[LSET.sender_set == s].ligand) for s in LSET.sender_set.unique()}
    MODSETS = {m: sorted(KP[KP.module == m].sender_set.unique()) for m in MODS}
    DENS = {m: sorted({c for s in MODSETS[m] for c in SETS[s]}) for m in MODS}
    print('sender sets per module:', MODSETS, flush=True)
    # ---------------- per-Treg computation
    PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); pcr = PC[PC.cohort == 'resection']
    assert len(pcr) == 6554
    if os.environ.get('SKIP_COMPUTE') != '1':
        sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
        sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
        sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
        allsend = set(c for v in SETS.values() for c in v); assert allsend <= set(sa.ct.unique()), allsend - set(sa.ct.unique())
        assert 'Treg' not in allsend and not any(c.startswith('SELENOP') for c in allsend)
        T = sa.loc[pcr.index, ['section']].copy(); assert (T.section.values == pcr.section.astype(str).values).all()
        LCOL = {(s, u): f'L|{s}|{u}' for s in NEED for u in NEED[s]}; RCOL = {u: f'R|{u}' for u in RUN}
        blocks = []
        for s_, f in FOLDER.items():
            sub = sa[sa.section == s_]
            with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
                ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
            assert set(GENES) <= set(v), set(GENES) - set(v)
            nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
            keep = np.where(nm.isin(sub.index))[0]
            E = pd.DataFrame(X[keep][:, [v.get_loc(g) for g in GENES]].toarray().astype(np.float32), index=nm[keep], columns=GENES); del X; gc.collect()
            cs = sub.loc[E.index]
            U = {u: E[u.split('_')].min(1).values for u in sorted(set(LUN) | set(RUN))}       # complex = minimum of subunits
            Z = {}
            for u, x in U.items():
                sd = x.std(); Z[u] = (x - x.mean()) / sd if sd > 0 else np.full_like(x, np.nan)
            tg = T[T.section == s_]; assert tg.index.isin(E.index).all()
            nb = cKDTree(cs[['x', 'y']].values).query_ball_point(sa.loc[tg.index, ['x', 'y']].values, r=RAD)
            ctv = cs.ct.values; out = pd.DataFrame(index=tg.index)
            out['nb_n'] = [len(i) for i in nb]
            rows_ = np.repeat(np.arange(len(nb)), [len(i) for i in nb]); cols_ = np.concatenate([np.asarray(i, dtype=int) for i in nb])
            for sname, types in SETS.items():
                msk = np.isin(ctv, types)
                sel = msk[cols_]; Nm = sparse.csr_matrix((np.ones(sel.sum(), np.float32), (rows_[sel], cols_[sel])), shape=(len(nb), len(cs)))
                n = np.asarray(Nm.sum(1)).ravel(); out[f'n|{sname}'] = n
                if sname in NEED:
                    ZL = np.column_stack([Z[u] for u in NEED[sname]]); ZL = np.nan_to_num(ZL, nan=0.0)
                    Lm = np.asarray(Nm @ ZL) / np.where(n > 0, n, 1)[:, None]; Lm[n < 1] = np.nan   # per-pair sender threshold applied in the module scores
                    for j, u in enumerate(NEED[sname]): out[LCOL[(sname, u)]] = Lm[:, j] if np.isfinite(Z[u]).any() else np.nan
                    if TOTAL:   # total supply: SUM over the group's cells within 50 um of log1p expression (>= 0; 0 if no sender)
                        UL = np.column_stack([U[u] for u in NEED[sname]]); assert (UL >= 0).all()
                        Tm = np.asarray(Nm @ UL)
                        for j, u in enumerate(NEED[sname]): out[f'T|{sname}|{u}'] = Tm[:, j]
            for m_, types in DENS.items(): out[f'dens|{m_}'] = np.asarray([np.isin(ctv[np.asarray(i, dtype=int)], types).sum() for i in nb])
            pos = E.index.get_indexer(tg.index)
            for u in RUN: out[RCOL[u]] = Z[u][pos]
            blocks.append(out); print(s_, 'Tregs', len(tg), 'cells', len(cs), flush=True)
            del E, U, Z; gc.collect()
        W = pd.concat(blocks).loc[T.index]; del sa; gc.collect()
        W.insert(0, 'section', T.section)
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
        assert (cn.loc[W.index, 'section'].astype(str).values == W.section.values).all()
        W.insert(1, 'CN', cn.loc[W.index, 'CN_k10'].astype(int).values)
        cnt = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0); assert W.index.isin(cnt.index).all()
        W.insert(2, 'log_total_counts', np.log1p(cnt.loc[W.index].sum(1).values)); del cnt; gc.collect()
        # module scores
        # each pair is scored from ITS OWN sender group and counts for a Treg if that group has >= MIN_SEND cells
        # within 50 um; module score = mean over the pairs that count; a Treg is scored if >= 1 pair counts (no all-groups requirement).
        # Sensitivity (CSV only): the same per-pair rule with >= MIN_SEND_SENS sender cells per pair (metrics 'L min3', 'LR min3', 'LxR raw min3').
        def Lmat(units, thr):
            return np.column_stack([np.where(W[f'n|{s}'].values >= thr, W[f'L|{s}|{l}'].values, np.nan) for l, s in units])
        for m in MODS:
            kp = KP[KP.module == m]; lu = kp[['ligand', 'sender_set']].drop_duplicates()
            Rp = np.column_stack([W[f'R|{r}'].values for r in kp.receptor])
            for thr, sfx in [(MIN_SEND, ''), (MIN_SEND_SENS, ' min3')]:
                Lp = Lmat(zip(kp.ligand, kp.sender_set), thr); cnt_ = np.isfinite(Lp)
                W[f'{m}|n_pairs_counting{sfx}'] = cnt_.sum(1)
                W[f'{m}|LR{sfx}'] = np.nanmean(np.where(cnt_, np.clip(Lp, 0, None) * np.clip(Rp, 0, None), np.nan), 1)
                W[f'{m}|LxR raw{sfx}'] = np.nanmean(np.where(cnt_, Lp * Rp, np.nan), 1)
                W[f'{m}|L{sfx}'] = np.nanmean(Lmat(zip(lu.ligand, lu.sender_set), thr), 1)
            W[f'{m}|Sender density'] = np.log1p(W[f'dens|{m}'].values)
            if TOTAL: W[f'{m}|T'] = np.mean(np.column_stack([W[f'T|{s}|{l}'].values for l, s in zip(kp.ligand, kp.sender_set)]), 1)   # mean over pairs
        W.to_csv(CACHE, float_format='%.4g')
    W = pd.read_csv(CACHE, index_col=0)
    assert all(f'n|{s}' in W for s in SETS) and f'{MODS[0]}|L min3' in W, 'cache is from a different sender/scoring definition - recompute without SKIP_COMPUTE'
    METRICS = ['L', 'LR', 'LxR raw', 'Sender density', 'L min3', 'LR min3', 'LxR raw min3']
    if TOTAL: assert f'{MODS[0]}|T' in W and W[[f'{m}|T' for m in MODS]].notna().all().all(); METRICS = ['T'] + METRICS
    print('Tregs scored per module (>=1 sender/pair):', {m: int(W[f'{m}|L'].notna().sum()) for m in MODS}, '| >=3 sender/pair:', {m: int(W[f'{m}|L min3'].notna().sum()) for m in MODS}, flush=True)
    # ---------------- CN names (as in 12_neighbourhood_treg_score.py)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [('OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
            ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
            ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
            ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    NMS = {k: short_cn(v) for k, v in NM.items()}
    # figure-only compact tick labels (narrow panel); dominant cell types of each CN, full names in the CSV (CN_name)
    FIGCN = {1: 'CN1: Epithelial', 2: 'CN2: Plasma', 3: 'CN3: Myofib / Neuron', 4: 'CN4: DC / Mac', 5: 'CN5: InfMo-Mac',
             7: 'CN7: Activ Fib / Neut', 8: 'CN8: B / Tfh', 10: 'CN10: OGN⁺RSPO3⁺ Fib'}
    for k_, v_ in FIGCN.items(): assert NMS[k_].startswith(v_.split(':')[0] + ':') and v_.split(': ')[1].split(' / ')[0][:5] in NMS[k_].replace('Myofibroblast', 'Myofib'), (k_, NMS[k_])
    CNS = list(range(1, 11))
    # ---------------- stats
    zs = lambda v: (v - v.mean()) / v.std()
    W['lib_z'] = W.groupby('section').log_total_counts.transform(zs)
    FIG_MODEL = os.environ.get('FIG_MODEL', 'ri')   # 'ri' random intercept (figure); 'rs' random slope = sensitivity (CSV)
    from scipy import stats as sst
    def fit_mm(dd, slope):
        """y ~ in_cn + lib_z, section random intercept (ri) or random intercept + random in_cn slope (rs); optimiser fallbacks if the Hessian fails."""
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = (smf.mixedlm('y ~ in_cn + lib_z', dd, groups=dd.section, re_formula='~in_cn') if slope else
                     smf.mixedlm('y ~ in_cn + lib_z', dd, groups=dd.section)).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], f.params['lib_z'], meth
            except Exception: pass
        return (np.nan,) * 4 + ('failed',)
    SEC, ST = [], []
    for met in METRICS:
        for m in MODS:
            col = f'{m}|{met}'; d = W[['section', 'CN', 'lib_z', col]].dropna().copy()
            d['y'] = d.groupby('section')[col].transform(zs)
            for c in CNS:
                secs = []
                for s_, g in d.groupby('section'):
                    inn, oth = g[g.CN == c], g[g.CN != c]
                    if len(inn) >= MIN_TREG_CN and len(oth) >= MIN_TREG_CN:
                        secs.append(s_)
                        SEC.append(dict(metric=met, module=m, CN=c, CN_name=NM[c], section=s_, n_tregs_in_CN=len(inn), n_tregs_other=len(oth),
                                        mean_score_in_CN=inn[col].mean(), mean_score_other=oth[col].mean(), mean_z_in_CN=inn.y.mean(),
                                        mean_z_other=oth.y.mean(), delta_z=inn.y.mean() - oth.y.mean()))
                o = dict(metric=met, module=m, CN=c, CN_name=NM[c], n_sections=len(secs), n_tregs_in_CN_all_sections=int((d.CN == c).sum()))
                if len(secs) >= 3:
                    dd = d[d.section.isin(secs)].assign(in_cn=lambda x: (x.CN == c).astype(float))
                    dl = np.array([r['delta_z'] for r in SEC if r['metric'] == met and r['module'] == m and r['CN'] == c])
                    o.update(n_tregs_in_CN=int(dd.in_cn.sum()), n_tregs_total=len(dd), median_section_delta=float(np.median(dl)),
                             section_ttest_p=sst.ttest_1samp(dl, 0).pvalue)
                    for tag, sl in [('ri', False), ('rs', True)]:
                        b, se, p_, lb, meth = fit_mm(dd, sl)
                        o.update({f'estimate_{tag}': b, f'ci_low_{tag}': b - 1.96 * se, f'ci_high_{tag}': b + 1.96 * se, f'p_{tag}': p_,
                                  f'lib_size_coef_{tag}': lb, f'optimiser_{tag}': meth,
                                  f'sections_agreeing_{tag}': int((np.sign(dl) == np.sign(b)).sum()) if np.isfinite(b) else np.nan})
                ST.append(o)
    SEC = pd.DataFrame(SEC); ST = pd.DataFrame(ST)
    for tag in ['ri', 'rs']:
        ST[f'q_{tag}'] = np.nan
        for (met, m), g in ST.groupby(['metric', 'module']):
            ok = g[f'p_{tag}'].notna(); ST.loc[g.index[ok], f'q_{tag}'] = multipletests(g[f'p_{tag}'][ok], method='fdr_bh')[1]
        ST[f'sections_agreeing_str_{tag}'] = [f'{int(a)}/{int(n)}' if pd.notna(a) else '' for a, n in zip(ST[f'sections_agreeing_{tag}'], ST.n_sections)]
    for k in ['estimate', 'ci_low', 'ci_high', 'p', 'q', 'sections_agreeing', 'sections_agreeing_str']: ST[k] = ST[f'{k}_{FIG_MODEL}']
    ST['figure_model'] = {'rs': 'random intercept + random in_CN slope by section', 'ri': 'random intercept by section'}[FIG_MODEL]
    print('failed fits:', ST[[c for c in ST if c.startswith('optimiser_')]].stack().value_counts().to_dict())
    pd.set_option('display.width', 260); pd.set_option('display.max_rows', 200)
    for met in METRICS:
        print('====', met); print(ST[ST.metric == met][['module', 'CN', 'n_sections', 'n_tregs_in_CN', 'estimate_ri', 'ci_low_ri', 'ci_high_ri', 'q_ri', 'sections_agreeing_str_ri', 'estimate_rs', 'q_rs', 'sections_agreeing_str_rs', 'section_ttest_p']].round(3).to_string(index=False))
    # ---------------- CSV for the panel (per-section rows + stats rows, 3 sig figs)
    sig3 = lambda x: x.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    A_ = SEC.assign(row_type='section').copy(); B_ = ST.assign(row_type='mixed_model', section='all (mixed model)').copy()
    OUTCSV = pd.concat([B_, A_], ignore_index=True)
    OUTCSV['metric_description'] = OUTCSV.metric.map({'T': 'TOTAL ligand supply: sum over sender cells within 50 um of log1p ligand expression, mean over pairs (main figure, _total version)',
                                                      'L': 'ligand supply from the specified senders (main figure)' if not TOTAL else 'per-sender mean ligand supply (sensitivity)', 'LR': 'ligand x receptor, positive parts (Ext figure)',
                                                      'LxR raw': 'raw ligand x receptor product (sensitivity)', 'Sender density': 'log1p sender cells within 50 um',
                                                      'L min3': 'ligand supply, sensitivity: pair counts only with >= 3 sender cells', 'LR min3': 'ligand x receptor, sensitivity: >= 3 sender cells per pair',
                                                      'LxR raw min3': 'raw product, sensitivity: >= 3 sender cells per pair'})
    OUTCSV['CN_label'] = OUTCSV.CN.map(NMS)
    OUTCSV['module_senders'] = OUTCSV.module.map(lambda m: '; '.join(f'{"/".join(sorted(KP[(KP.module == m) & (KP.sender_set == s)].pathway.unique()))}: {", ".join(SETS[s])}' for s in MODSETS[m]))
    OUTCSV['n_pairs'] = OUTCSV.module.map(KP.groupby('module').size())
    cols = ['row_type', 'metric', 'metric_description', 'module', 'CN', 'CN_label', 'CN_name', 'section', 'n_sections', 'n_tregs_in_CN', 'n_tregs_in_CN_all_sections', 'n_tregs_other', 'n_tregs_total',
            'figure_model', 'estimate', 'ci_low', 'ci_high', 'p', 'q', 'sections_agreeing_str',
            'estimate_ri', 'ci_low_ri', 'ci_high_ri', 'p_ri', 'q_ri', 'sections_agreeing_str_ri', 'lib_size_coef_ri', 'optimiser_ri',
            'estimate_rs', 'ci_low_rs', 'ci_high_rs', 'p_rs', 'q_rs', 'sections_agreeing_str_rs', 'lib_size_coef_rs', 'optimiser_rs',
            'median_section_delta', 'section_ttest_p', 'delta_z', 'mean_z_in_CN', 'mean_z_other', 'mean_score_in_CN', 'mean_score_other', 'n_pairs', 'module_senders']
    OUTCSV = OUTCSV[cols].apply(sig3)
    OUTCSV.to_csv(f'{STEM}.csv', index=False)
    # ---------------- figures
    SHORTM = {'Antigen presentation': 'Antigen Presentation', 'Co-stimulation / suppressive Treg': 'Co-Stimulation / Suppressive Treg',
              'Inflammatory cytokine': 'Inflammatory Cytokine', 'Matrix-integrin': 'Matrix–Integrin'}
    COL = {'Antigen presentation': ('#9ecae1', '#08519c'), 'Co-stimulation / suppressive Treg': ('#bcbddc', '#54278f'),
           'Inflammatory cytokine': ('#fcae91', '#a50f15'), 'Matrix-integrin': ('#a1d99b', '#006d2c')}
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    YL = {'T': 'Total Ligand Supply to Tregs (z vs Other Tregs)', 'LR': 'Ligand × Receptor Score in Tregs (z vs Other Tregs)', 'L': 'Ligand Supply to Tregs (z vs Other Tregs)',
          'Sender density': 'Sender Cells near Tregs (z vs Other Tregs)'}
    EXCL = [6, 9]   # too few Tregs (>= 10 Tregs in < 3 sections)
    PCNS = [c for c in CNS if c not in EXCL]
    for met in ['L', 'LR']:
        for c in EXCL: assert (W[W.CN == c].groupby('section').size() >= MIN_TREG_CN).sum() < 3, c   # < 3 sections with >= 10 Tregs, scoring aside
    NT = ST[ST.metric.isin(['T', 'L', 'LR']) & ST.CN.isin(PCNS) & ST.estimate.isna()][['metric', 'module', 'CN', 'n_sections']]
    print('plotted CNs:', PCNS, '| not testable within plotted CNs:'); print(NT.to_string(index=False))
    def draw(met, path_stem):   # PNG only
        CNS = PCNS
        plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
        fig, axes = plt.subplots(4, 1, figsize=(6.0, 12), sharex=True, gridspec_kw={'hspace': 0.24})
        rng = np.random.default_rng(1); NTX = []
        for ax, m in zip(axes, MODS):
            lc, dc = COL[m]; g = ST[(ST.metric == met) & (ST.module == m)].set_index('CN'); s = SEC[(SEC.metric == met) & (SEC.module == m)]
            ax.axhline(0, color='#888888', lw=1, zorder=0)
            for i, c in enumerate(CNS):
                tested = c in g.index and pd.notna(g.loc[c, 'estimate'])
                v = s[s.CN == c].delta_z.values if tested else np.array([])   # untested cells: no dots (1-2 sections; values in the CSV)
                ax.scatter(i + rng.uniform(-0.14, 0.14, len(v)), v, s=45, color=lc, edgecolor='white', linewidth=0.5, alpha=0.95, zorder=2)
                if c in g.index and pd.notna(g.loc[c, 'estimate']):
                    e, lo, hi, q = g.loc[c, ['estimate', 'ci_low', 'ci_high', 'q']]
                    ax.errorbar(i, e, yerr=[[e - lo], [hi - e]], fmt='o', ms=10, color=dc, mfc=dc if q < 0.05 else 'white', mew=2.2, elinewidth=2.6, capsize=0, zorder=3)
                else: NTX.append((ax, i))
            vals = np.concatenate([s[s.CN.isin([c for c in CNS if c in g.index and pd.notna(g.loc[c, 'estimate'])])].delta_z.values, g.loc[g.index.isin(CNS), 'ci_low'].dropna().values, g.loc[g.index.isin(CNS), 'ci_high'].dropna().values])
            ymin, ymax = vals.min(), vals.max(); span = ymax - ymin
            stop = [ymax]
            for i, c in enumerate(CNS):
                if c in g.index and pd.notna(g.loc[c, 'q']) and g.loc[c, 'q'] < 0.05:
                    v = s[s.CN == c].delta_z.values; top = max(v.max() if len(v) else -np.inf, g.loc[c, 'ci_high'])
                    ax.text(i, top + 0.012 * span, stars(g.loc[c, 'q']), ha='center', va='bottom', fontsize=17, color=dc, fontweight='bold')
                    stop.append(top + 0.2 * span)   # room for the star text below the panel title
            ax.set_ylim(ymin - 0.05 * span, max(ymax + 0.05 * span, max(stop)))
            for a_, i in NTX:
                if a_ is ax: ax.text(i, 0, 'n.t.', ha='center', va='center', fontsize=14, color='#969696', style='italic')
            ax.set_title(SHORTM[m], fontsize=17, loc='left', pad=3, fontweight='bold')
            ax.tick_params(axis='y', labelsize=15)
            ax.spines[['top', 'right']].set_visible(False); ax.set_xlim(-0.55, len(CNS) - 0.45)
        axes[-1].set_xticks(range(len(CNS))); axes[-1].set_xticklabels([FIGCN[c] for c in CNS], rotation=50, ha='right', rotation_mode='anchor', fontsize=14)
        axes[-1].set_xlabel('Neighbourhood of the Treg', fontsize=16)
        fig.supylabel(YL[met], fontsize=16, x=-0.045)
        h = [Line2D([], [], marker='o', ls='', color='#bdbdbd', ms=9, label='One Section'),
             Line2D([], [], marker='o', ls='-', color='#252525', mfc='#252525', ms=10, lw=2.4, label='Estimate ± 95% CI, q < 0.05'),
             Line2D([], [], ls='', label=' '),
             Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=10, mew=2, lw=2.4, label='q ≥ 0.05')]   # rows: section | blank ; filled | open
        axes[0].legend(handles=h, loc='lower left', bbox_to_anchor=(-0.02, 1.12), ncol=2, frameon=False, fontsize=14, handletextpad=0.3, columnspacing=0.8, borderaxespad=0)
        fig.savefig(f'{path_stem}.png', dpi=300, bbox_inches='tight', pad_inches=0.04)
        plt.close(fig)
    if TOTAL: draw('T', STEM)
    else:
        draw('L', STEM)
        draw('LR', STEM_EXT)
        draw('Sender density', f'{OUTB}/treg_receiver_modules_sender_density')
    # ---------------- legends
    npairs = KP.groupby('module').size(); nlig = KP[['module', 'ligand', 'sender_set']].drop_duplicates().groupby('module').size()
    nsc = {m: int(W[f'{m}|L'].notna().sum()) for m in MODS}
    ex_txt = (f'{" and ".join(f"CN{c}" for c in EXCL)} are omitted because they contain too few Tregs (' + ' and '.join(f'{int((W.CN == c).sum())}' for c in EXCL)
              + ' Tregs in total, respectively; at least 10 Tregs in fewer than 3 sections). ')
    common = ('CD resections (CosMx 6K; 11 sections, 4 patients; '
          f'{len(W):,} strict Tregs as receivers). Ligand-receptor pairs were taken from the NicheNet v2 prior network, matching Fig 3, and grouped '
          'into functional modules using CellChatDB pathway annotations (a pair joins a module when a CellChatDB interaction of the module\'s '
          'pathways shares its ligand gene and contains its receptor gene; pairs matching no module are not used). Modules: antigen presentation = MHC-II; co-stimulation / '
          'suppressive Treg = CD80, CD86, PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2, IL10; inflammatory cytokine = IL1, IL6, TNF, IFN-II, IL12, '
          'IL17, OSM; matrix–integrin = all ECM-Receptor pathways; pairs with both genes on the 6K panel: '
          f'{npairs["Antigen presentation"]}, {npairs["Co-stimulation / suppressive Treg"]}, {npairs["Inflammatory cytokine"]} and {npairs["Matrix-integrin"]}. Each ligand is counted only from its main canonical sender cell '
          'types: MHC-II from DC, B cells, CD9+ M2 macrophages and Mo-Macs; the APC co-stimulatory and co-inhibitory ligands (CD80, CD86, ICOSLG, '
          'CD274, PDCD1LG2, CD70, TNFSF4, TNFSF18) from DC and B cells; IL2-family cytokines and CD40LG from CD4 effector, CD4 Trm and Tfh cells; '
          'IL10-family cytokines from DC, CD9+ M2 macrophages and Mo-Macs; TNF, IL1-family (including IL18), IL6-family, IL12-family (including '
          'IL23 and IL27) and OSM from Mo-Macs, CD9+ M2 macrophages, DC and neutrophils; IFNG and IL17-family from CD4 effector, CD8 Trm, NK and '
          'γδ T cells; extracellular-matrix ligands from fibroblasts (OGN+RSPO3+, ADAMDEC1+, activated, crypt-top and PLA2G2A+ ECM) and '
          'myofibroblasts. For each Treg and pair, ligand supply is the mean within-section z of the ligand over that pair\'s sender cells within '
          '50 µm; each pair is scored from its own sender group and counts for a Treg when that group has at least one cell within 50 µm, '
          'the module score is the mean over the pairs that count, and a Treg is scored when at least one pair counts '
          f'(Tregs scored: {nsc["Antigen presentation"]:,}, {nsc["Co-stimulation / suppressive Treg"]:,}, {nsc["Inflammatory cytokine"]:,} and '
          f'{nsc["Matrix-integrin"]:,}). Each Treg is assigned to its own neighbourhood. ')
    stats_txt = ('the score is z-standardised within section over the scored Tregs, and each value is the mean of Tregs in the neighbourhood minus the '
             'mean of the other Tregs of the same section (0 = no difference). Light dots: one section each (sections with at least 10 Tregs in the '
             'neighbourhood). Dark points and bars: linear mixed model, score ~ in neighbourhood + Treg library size (log total counts, within-section z) '
             '+ (1 | section), estimate and 95% CI. Filled points: BH q < 0.05 over neighbourhoods within each module; * q < 0.05, ** q < 0.01, '
             '*** q < 0.001. ' + ex_txt + '{NT}')
    legL = ('Ligand supply to Tregs from the canonical sender cells of each ligand-receptor module, across spatial cellular neighbourhoods, ' + common
            + f'The module\'s ligand supply is the mean over its unique ligands ({nlig["Antigen presentation"]}, {nlig["Co-stimulation / suppressive Treg"]}, '
        f'{nlig["Inflammatory cytokine"]} and {nlig["Matrix-integrin"]}); receptor expression is not used, so the score measures the signal reaching the Treg. '
        'Y axis ("Ligand Supply to Tregs (z vs Other Tregs)"): ' + stats_txt
            + 'The ligand × receptor version is in Extended Data. The source data hold the per-section values, the random-intercept model shown and a '
        'random-intercept plus random-slope model (1 + in neighbourhood | section) as a sensitivity check, the same analysis requiring at least '
        '3 sender cells per pair (sensitivity check), and the ligand × receptor, raw-product and sender-density metrics. Neighbourhood names are abbreviated; full names are in the source data.')
    legLR = ('Ligand × receptor signalling received by Tregs from the canonical sender cells of each ligand-receptor module, across spatial cellular '
         'neighbourhoods (Extended Data companion to the Fig 4 ligand-supply panel; same Tregs, pairs, senders and statistics), ' + common
             + 'The receptor term is the Treg\'s own within-section z of the receptor; the pair score is the product of the positive parts of ligand supply '
         'and receptor, and the module score is the mean over pairs. Y axis ("Ligand × Receptor Score in Tregs (z vs Other Tregs)"): ' + stats_txt
             + 'Per-section values and both mixed models (random intercept shown; random intercept plus random slope as a sensitivity check) are in the '
         'Fig 4 ligand-supply source data (metric = LR). Neighbourhood names are abbreviated; full names are in the source data.')
    def nt_txt(met):
        x = NT[NT.metric == met]
        if not len(x): return ''
        parts = [f'{m.replace("Matrix-integrin", "matrix–integrin").lower().replace("treg", "Treg")} in {", ".join(f"CN{c}" for c in g.CN)}' for m, g in x.groupby('module', sort=False)]
        return ('n.t., not tested: ' + '; '.join(parts) + ' (fewer than 3 sections with at least 10 scored Tregs in the neighbourhood). ')
    legL = legL.replace('{NT}', nt_txt('L')); legLR = legLR.replace('{NT}', nt_txt('LR'))
    for t in (legL, legLR): assert 'activation' not in t.lower() and 'quiescent' not in t.lower()
    if TOTAL:
        legT = ('Total ligand supply to Tregs from the canonical sender cells of each ligand-receptor module, across spatial cellular neighbourhoods, '
                + common.split('For each Treg and pair, ligand supply is')[0]
                + 'For each Treg and pair, total ligand supply is the sum, over that pair\'s sender cells within 50 µm, of the ligand\'s log-normalised '
            'expression, so it reflects both the number of sender cells and their expression (a pair with no sender cell nearby contributes 0, so '
            f'every Treg is scored; {len(W):,} Tregs); the module score is the mean over the module\'s pairs. Each Treg is assigned to its own neighbourhood. '
            'Receptor expression is not used. Y axis ("Total Ligand Supply to Tregs (z vs Other Tregs)"): ' + stats_txt.replace('{NT}', nt_txt('T'))
                + 'The source data hold the per-section values, the random-intercept model shown, a random-intercept plus random-slope model '
            '(1 + in neighbourhood | section) and the per-sender mean ligand supply (average over sender cells, which does not scale with the number '
            'of senders) as sensitivity analyses. Neighbourhood names are abbreviated; full names are in the source data.')
        assert 'activation' not in legT.lower() and 'quiescent' not in legT.lower()
        open(f'{STEM}_legend.txt', 'w').write(legT + '\n')
    else:
        open(f'{STEM}_legend.txt', 'w').write(legL + '\n')
        open(f'{STEM_EXT}_legend.txt', 'w').write(legLR + '\n')
    print('done', flush=True)


# ======================================================================================================================
# Section cn_lr_modules_celltype_restricted: LR module activity per CN restricted to the sender and receiver cell types of each pair
# Panels: LR module inputs (cell-type-restricted reference columns of cn_lr_modules_allcells)
# ======================================================================================================================
if _run('cn_lr_modules_celltype_restricted'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/79_lr_modules_niche.py'   # working folder of this step (relative paths resolve from here)
    """Ligand-receptor module activity per cellular neighbourhood, restricted to the cell types that send and receive each signal. This characterises the
niche; it is not signalling to Tregs.
Pairs, module labels and pathway-level sender sets as in treg_receiver_lr_modules (the pair table is asserted identical to treg_receiver_pairs_kept.csv).
Senders (cell_type_general): MHC-II -> DC, B cell, CD9+ M2 Mac, Mo-Mac, SELENOP+ Mac; APC co-stimulatory ligands -> DC, B cell; IL2-family + CD40LG -> CD4
  Effector, CD4 Trm, Tfh; IL10-family -> DC, CD9+ M2 Mac, Mo-Mac, SELENOP+ Mac; TNF / IL1 / IL6 / IL12 / OSM -> Mo-Mac, CD9+ M2 Mac, SELENOP+ Mac, DC,
  Neutrophil; IFNG / IL17 -> CD4 Effector, CD8 Trm, NK, Gd T; ECM -> OGN+RSPO3+, ADAMDEC1+, Activated, Crypt Top, PLA2G2A+ ECM Fib, Myofibroblast.
  Per-ligand additions from established biology, erring wide (LIG_SENDERS: epithelial / fibroblast / endothelial / myeloid / T sources for IL33,
  IL36A/B/G, IL11, IL17B/C, IL25, IL7, IL15, TSLP, IL6, TNF; epithelial includes Cycling).
Receivers (per pair):
  Antigen presentation (HLA-* -> CD4)          -> CD4 T: CD4 Effector, CD4 Trm, Tfh, Exhausted CD4 Tfh, Resting T, Treg
  Co-stim APC ligands (CD80/CD86/ICOSLG/PD-L1/PD-L2/CD70/OX40L/GITRL) and IL2-family -> all T + NK: CD4 T above, CD8 Trm, Gd T, NK
  CD40LG -> CD40                               -> DC, B cell, CD9+ M2 Mac, Mo-Mac, SELENOP+ Mac
  CD40LG -> integrins                          -> myeloid (DC, CD9+ M2 Mac, Mo-Mac, SELENOP+ Mac, Neutrophil, Granulocyte) + all T + NK
  IL10                                         -> all immune cells (myeloid, B, plasma, T, NK)
  IL10-family IL19 / IL20 / IL22 / IL24        -> epithelial (Enterocyte, Absorptive Epi, Secretory Epi, Microfold, Stem, Cycling) + fibroblasts
  IL10-family IL26                             -> epithelial + myeloid + T / NK
  Inflammatory cytokine, Matrix-integrin       -> all cell types except SELENOP+ Fib (receptor expression gates the score)
  SELENOP+ Fib is never a sender or receiver.
Score (per section): neighbours = all cells within 50 um (excluding the cell itself). For receiver cell r of the pair's receiver set, ligand = mean
  log1p-normalised expression over the pair's sender cells within 50 um (0 if none); pair activity = ligand x r's own receptor log1p-normalised expression.
  Main ('receiver_mean_zero'): cell score = mean over the module's pairs it receives; CN value = mean over all receiver cells of the CN; CN x section kept
  if >= 100 cells and >= 20 receiver cells.
Sensitivities (CSV only): 'per_receiver_mean' = receivers with >= 1 sender within 50 um only (>= 20 such receivers); 'total' = abundance-weighted (summed
  sender ligand, summed over the CN's cells and divided by its number of cells, non-receivers 0); 'min3' = 'total' with pairs counted only where the
  sender group has >= 3 cells within 50 um.
Heatmap value ('Mean Within-Section z'): CN values z-standardised across CNs within each section, then averaged over sections. Stars: per module x CN,
  linear mixed model on the per-section z values, z ~ in_CN + (1 | section), BH over the 10 CNs within module (Wilcoxon signed-rank q in the CSV).
  All 10 CNs are shown. The CSV also carries the section-adjusted relative score of cn_module_relative_scores as a reference column.
Outputs: A_fig4_panels/fig4_lr_modules_niche_heatmap.{png,csv} + _legend.txt; A_fig4_panels/supp_lr_modules_niche.{png,csv} + _legend.txt (read by
  cn_lr_modules_allcells); B_resection_analyses/lr_modules_cellchat/lr_pairs_senders_receivers.csv and lr_modules_niche_per_section_cn.csv (cache).
SKIP_COMPUTE=1 redraws from the cache."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX' and os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq, pyreadr
    from anndata.io import read_elem
    from scipy import sparse, stats as sst
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.lines import Line2D
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    HSTEM = f'{OUTA}/fig4_lr_modules_niche_heatmap'; SSTEM = f'{OUTA}/supp_lr_modules_niche'; CACHE = f'{OUTB}/lr_modules_niche_per_section_cn.csv'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; UM = 0.12028; RAD = 50; MIN_CELLS_CN = 100; MIN_REC_CN = 20
    MAIN = 'receiver_mean_zero'   # receivers without a sender within 50 um contribute 0 (no border-cell selection);
                                  # sensitivities: 'per_receiver_mean' (with-sender receivers only), 'total' / 'min3' (abundance-weighted)
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    # ---------------- modules and sender sets (treg_receiver_lr_modules sets, SELENOP+ Mac added to the myeloid sets)
    MODS = ['Antigen presentation', 'Co-stimulation / suppressive Treg', 'Inflammatory cytokine', 'Matrix-integrin']
    MOD = {'Antigen presentation': ['MHC-II'], 'Co-stimulation / suppressive Treg': ['CD80', 'CD86', 'PD-L1', 'PDL2', 'ICOS', 'OX40', 'GITRL', 'CD70', 'CD40', 'IL2', 'IL10'],
           'Inflammatory cytokine': ['IL1', 'IL6', 'TNF', 'IFN-II', 'IL12', 'IL17', 'OSM']}
    cc = pd.read_csv(f'{D}/objects/cellchatdb_human_interactions_snapshot.csv')
    MOD['Matrix-integrin'] = sorted(cc.loc[cc.annotation == 'ECM-Receptor', 'pathway_name'].unique())
    SETS = {'MHCII_APC': ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac'], 'DC_B': ['DC', 'B cell'], 'CD4_T': ['CD4 Effector', 'CD4 Trm', 'Tfh'],
            'IL10_myeloid': ['DC', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac'], 'Inflam_myeloid': ['Mo-Mac', 'CD9+ M2 Mac', 'SELENOP+ Mac', 'DC', 'Neutrophil'],
            'Inflam_T_NK': ['CD4 Effector', 'CD8 Trm', 'NK', 'Gd T'],
            'Fibroblast': ['OGN+RSPO3+ Fib', 'ADAMDEC1+ Fib', 'Activated Fib', 'Crypt Top Fib', 'PLA2G2A+ ECM Fib', 'Myofibroblast']}
    # senders and receivers from established biology, erring wide; per-ligand additions to the pathway-level sender sets above
    EPI = ['Enterocyte', 'Absorptive Epi', 'Absoprtive Epi', 'Secretory Epi', 'Microfold', 'Stem', 'Cycling']   # epithelial incl. Cycling
    FIB = SETS['Fibroblast']; MYE3 = ['DC', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac']; ENDO = ['Endothelial']
    LIG_SENDERS = {'IL33': SETS['Inflam_myeloid'] + EPI + ENDO + FIB, 'IL36A': SETS['Inflam_myeloid'] + EPI, 'IL36B': SETS['Inflam_myeloid'] + EPI,
                   'IL36G': SETS['Inflam_myeloid'] + EPI, 'IL11': SETS['Inflam_myeloid'] + FIB, 'IL17B': SETS['Inflam_T_NK'] + EPI,
                   'IL17C': SETS['Inflam_T_NK'] + EPI, 'IL25': SETS['Inflam_T_NK'] + EPI, 'IL7': SETS['CD4_T'] + EPI + FIB + ENDO,
                   'IL15': SETS['CD4_T'] + MYE3 + EPI + FIB, 'TSLP': SETS['CD4_T'] + EPI + FIB, 'IL6': SETS['Inflam_myeloid'] + FIB + ENDO,
                   'TNF': SETS['Inflam_myeloid'] + SETS['Inflam_T_NK']}
    for l_, t_ in LIG_SENDERS.items(): SETS[f'S_{l_}'] = list(dict.fromkeys(t_))
    assert not any(c == 'SELENOP+ Fib' for v_ in SETS.values() for c in v_)   # SELENOP+ Mac allowed, SELENOP+ Fib excluded
    COSTIM_APC = ('CD80', 'CD86', 'ICOS', 'PD-L1', 'PDL2', 'CD70', 'OX40', 'GITRL')
    def sender_set(module, pathway):
        if module == 'Antigen presentation': assert pathway == 'MHC-II'; return 'MHCII_APC'
        if module == 'Co-stimulation / suppressive Treg':
            if pathway in ('IL2', 'CD40'): return 'CD4_T'
            if pathway == 'IL10': return 'IL10_myeloid'
            assert pathway in COSTIM_APC, pathway; return 'DC_B'
        if module == 'Inflammatory cytokine':
            if pathway in ('IFN-II', 'IL17'): return 'Inflam_T_NK'
            assert pathway in ('IL1', 'IL6', 'TNF', 'IL12', 'OSM'), pathway; return 'Inflam_myeloid'
        return 'Fibroblast'
    # ---------------- receiver sets
    ALL_TYPES = ['Myofibroblast', 'Enterocyte', 'Plasma', 'Secretory Epi', 'DC', 'OGN+RSPO3+ Fib', 'Activated Fib', 'Microfold', 'NK', 'SELENOP+ Mac',
                 'Absorptive Epi', 'CD8 Trm', 'Stem', 'Pericyte', 'Endothelial', 'B cell', 'CD4 Trm', 'ADAMDEC1+ Fib', 'Crypt Top Fib', 'Granulocyte',
                 'PLA2G2A+ ECM Fib', 'Neuron', 'Tfh', 'Gd T', 'Absoprtive Epi', 'Treg', 'Exhausted CD4 Tfh', 'Resting T', 'CD9+ M2 Mac', 'SELENOP+ Fib',
                 'Neutrophil', 'CD4 Effector', 'Mo-Mac', 'Cycling']   # cell_type_general of all_samples.parquet ('Absoprtive Epi' = source typo)
    CD4T = ['CD4 Effector', 'CD4 Trm', 'Tfh', 'Exhausted CD4 Tfh', 'Resting T', 'Treg']
    TALL = CD4T + ['CD8 Trm', 'Gd T', 'NK']
    MYE5 = ['DC', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac', 'Neutrophil', 'Granulocyte']
    RSETS = {'CD4_T_all': CD4T, 'T_NK_all': TALL, 'APC_CD40': ['DC', 'B cell', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac'],
             'Myeloid_T_NK': MYE5 + TALL, 'IL10_myeloid_T_NK': ['DC', 'CD9+ M2 Mac', 'Mo-Mac', 'SELENOP+ Mac'] + TALL,
             'Immune_all': MYE5 + ['B cell', 'Plasma'] + TALL,
             'Epithelial_Fib': EPI + FIB, 'Epithelial_Myeloid_T_NK': EPI + MYE5 + TALL,
             'All_non_SELENOP_Fib': [c for c in ALL_TYPES if c != 'SELENOP+ Fib']}
    assert not any(c == 'SELENOP+ Fib' for v_ in RSETS.values() for c in v_)
    def receiver_set(module, pathway, ligand, receptor):
        if module == 'Antigen presentation': return 'CD4_T_all'
        if module == 'Co-stimulation / suppressive Treg':
            if pathway in COSTIM_APC or pathway == 'IL2': return 'T_NK_all'
            if pathway == 'CD40': return 'APC_CD40' if receptor == 'CD40' else 'Myeloid_T_NK'
            assert pathway == 'IL10', pathway
            if ligand == 'IL10': return 'Immune_all'                                   # IL10R broadly on leukocytes
            if ligand in ('IL19', 'IL20', 'IL22', 'IL24'): return 'Epithelial_Fib'      # IL20RA/IL20RB/IL22RA1: epithelium + some stroma
            assert ligand == 'IL26', ligand; return 'Epithelial_Myeloid_T_NK'          # IL20RA/IL10RB on epithelium; IL26 also acts on DC / T
        return 'All_non_SELENOP_Fib'
    # ---------------- pairs (as treg_receiver_lr_modules)
    gs = lambda s: tuple(sorted(x.strip() for x in str(s).replace(';', ',').split(',') if x.strip() and x.strip() != 'nan'))
    with h5py.File(f'{R}/{FOLDER["43563_A8"]}/Processed/norm_log1p.h5ad', 'r') as h: PANEL = set(read_elem(h['var']).index.astype(str))
    NN = list(pyreadr.read_r('/share/fsmresfiles/ydn4687/nichenet/lr_network_human_21122021.rds').values())[0][['from', 'to']].drop_duplicates()
    ccL = cc['ligand.symbol'].map(lambda x: frozenset(gs(x))); ccR = cc['receptor.symbol'].map(lambda x: frozenset(gs(x)))
    PR = []
    for m in MODS:
        sub = cc.pathway_name.isin(MOD[m])
        for l, r in zip(NN['from'], NN['to']):
            hit = sorted({pw for pw, a, b in zip(cc.pathway_name[sub], ccL[sub], ccR[sub]) if l in a and r in b})
            if not hit: continue
            ss = {sender_set(m, pw) for pw in hit}; assert len(ss) == 1, (m, l, r, hit)
            if l in LIG_SENDERS: ss = {f'S_{l}'}
            rs = {receiver_set(m, pw, l, r) for pw in hit}; assert len(rs) == 1, (m, l, r, hit)
            miss = [g for g in (l, r) if g not in PANEL]
            PR.append(dict(module=m, pathway='/'.join(hit), interaction=f'{l}_{r}', ligand=l, receptor=r, sender_set=ss.pop(), receiver_set=rs.pop(),
                           kept=not miss, missing_genes=','.join(miss), pair_source='NicheNet v2 lr_network_human_21122021'))
    PR = pd.DataFrame(PR); PR['sender_cell_types'] = PR.sender_set.map(lambda s: ', '.join(SETS[s])); PR['receiver_cell_types'] = PR.receiver_set.map(lambda s: ', '.join(RSETS[s]))
    assert not PR.duplicated(['module', 'ligand', 'receptor']).any()
    kp74 = pd.read_csv(f'{OUTB}/treg_receiver_pairs_kept.csv')   # treg_receiver_lr_modules pair table: the pairs must be identical (only sender sets differ)
    assert (PR[['module', 'ligand', 'receptor', 'kept']].reset_index(drop=True) == kp74[['module', 'ligand', 'receptor', 'kept']]).all().all()
    assert set(LIG_SENDERS) <= set(PR[PR.kept].ligand), set(LIG_SENDERS) - set(PR[PR.kept].ligand)
    PR.to_csv(f'{OUTB}/lr_pairs_senders_receivers.csv', index=False)
    KP = PR[PR.kept].reset_index(drop=True)
    print(KP.groupby(['module', 'sender_set', 'receiver_set']).size().to_string(), flush=True)
    LSET = KP[['ligand', 'sender_set']].drop_duplicates(); NEED = {s: sorted(LSET[LSET.sender_set == s].ligand) for s in LSET.sender_set.unique()}
    GENES = sorted(set(KP.ligand) | set(KP.receptor)); CNS = list(range(1, 11))
    VARS = ['receiver_mean_zero', 'per_receiver_mean', 'total', 'min3']
    # ---------------- per-section computation -> per (section, CN, module) values
    if os.environ.get('SKIP_COMPUTE') != '1':
        sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
        sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]
        sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
        assert set(sa.ct.unique()) == set(ALL_TYPES), set(sa.ct.unique()) ^ set(ALL_TYPES)
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
        assert cn.index.isin(sa.index).all()
        sa['CN'] = cn.CN_k10.reindex(sa.index).values   # NaN = no neighbourhood assigned (0.4% of cells); they still send
        ROWS = []
        for s_, f in FOLDER.items():
            sub = sa[sa.section == s_]
            with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
                ob = read_elem(h['obs']); v = read_elem(h['var']).index.astype(str); X = read_elem(h['X'])
            nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
            keep = np.where(nm.isin(sub.index))[0]
            E = pd.DataFrame(X[keep][:, [v.get_loc(g) for g in GENES]].toarray().astype(np.float32), index=nm[keep], columns=GENES); del X; gc.collect()
            assert (E.values >= 0).all()
            cs = sub.loc[E.index]; n = len(cs); ct = cs.ct.values; cnv = cs.CN.values
            pr = cKDTree(cs[['x', 'y']].values).query_pairs(RAD, output_type='ndarray')   # neighbours within 50 um, self excluded
            A = sparse.csr_matrix((np.ones(2 * len(pr), np.float32), (np.r_[pr[:, 0], pr[:, 1]], np.r_[pr[:, 1], pr[:, 0]])), shape=(n, n)); del pr
            NS, LS = {}, {}
            for sname, lig in NEED.items():
                msk = np.isin(ct, SETS[sname]).astype(np.float32)
                NS[sname] = np.asarray(A @ msk).ravel(); LS[sname] = np.asarray(A @ (E[lig].values * msk[:, None]))
            del A; gc.collect()
            RM = {k: np.isin(ct, t) for k, t in RSETS.items()}
            for m in MODS:
                kp = KP[KP.module == m]; P = len(kp)
                tot = np.zeros(n); m3 = np.zeros(n); pnum = np.zeros(n); pden = np.zeros(n); znum = np.zeros(n); zden = np.zeros(n)
                for l, r, ss, rs in zip(kp.ligand, kp.receptor, kp.sender_set, kp.receiver_set):
                    Ls = LS[ss][:, NEED[ss].index(l)]; ns = NS[ss]; rm = RM[rs]; Rv = E[r].values
                    act = rm * Ls * Rv; tot += act; m3 += act * (ns >= 3)
                    ok = rm & (ns >= 1); pnum += np.where(ok, Ls / np.maximum(ns, 1) * Rv, 0.0); pden += ok
                    znum += np.where(rm, Ls / np.maximum(ns, 1) * Rv, 0.0); zden += rm   # MAIN: mean sender ligand, 0 if no sender; every receiver counted
                tot /= P; m3 /= P; pm = np.where(pden > 0, pnum / np.maximum(pden, 1), np.nan); rz = np.where(zden > 0, znum / np.maximum(zden, 1), np.nan)
                isrec = np.zeros(n, bool)
                for rs in kp.receiver_set.unique(): isrec |= RM[rs]
                for c in CNS:
                    ix = cnv == c; nc = int(ix.sum())
                    if nc == 0: continue
                    ROWS.append(dict(section=s_, CN=c, module=m, n_cells=nc, n_receiver_cells=int(isrec[ix].sum()), n_receivers_with_sender=int((pden[ix] > 0).sum()),
                                     total=tot[ix].sum() / nc, min3=m3[ix].sum() / nc, per_receiver_mean=np.nanmean(pm[ix]) if (pden[ix] > 0).any() else np.nan,
                                     receiver_mean_zero=np.nanmean(rz[ix]) if isrec[ix].any() else np.nan))
            print(s_, 'cells', n, flush=True); del E, LS, NS, RM; gc.collect()
        pd.DataFrame(ROWS).to_csv(CACHE, index=False, float_format='%.6g')
    V = pd.read_csv(CACHE)
    # ---------------- CN names (long labels, as module_by_CN_pooled_vs_adjusted.csv)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [(' Fib / PLA2G2A⁺ Fib / Endothelial', '/PLA2G2A⁺ Fib'), (' / Crypt Top Fib / Cycling', ' / CryptTop / Cycl'),
            (' / ADAMDEC1 Fib', ' / ADAMDEC1'), (' / Neuron / Pericyte', ' / Neuron / Peri'), (' / Mac S_M_P / CD9 Mac', ' / Mac / CD9 Mac'),
            ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'), ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'),
            ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {c: short_cn(NM[c]) for c in CNS}
    rel63 = pd.read_csv(f'{OUTB}/module_by_CN_pooled_vs_adjusted.csv', header=[0, 1], index_col=0)
    assert [short_cn(i) for i in rel63.index] == [LAB[c] for c in CNS], 'row labels must equal module_by_CN_pooled_vs_adjusted.csv'
    RELZ = rel63['section_adjusted_z'][MODS].set_axis(CNS)
    FIGCN = {1: 'CN1: Epithelial', 2: 'CN2: Plasma', 3: 'CN3: Myofib / Neuron', 4: 'CN4: DC / Mac', 5: 'CN5: InfMo-Mac', 6: 'CN6: Epithelial',
             7: 'CN7: Activ Fib / Neut', 8: 'CN8: B / Tfh', 9: 'CN9: Myofibroblast', 10: 'CN10: OGN⁺RSPO3⁺ Fib'}
    # ---------------- within-section z across CNs, mean over sections; LMM stars
    zs = lambda x: (x - x.mean()) / x.std()
    V = V[(V.n_cells >= MIN_CELLS_CN) & (V.n_receiver_cells >= MIN_REC_CN)].copy()
    V.loc[V.n_receivers_with_sender < MIN_REC_CN, 'per_receiver_mean'] = np.nan   # with-sender sensitivity keeps its own >= 20 rule
    for var in VARS: V[f'z_{var}'] = V.groupby(['section', 'module'])[var].transform(zs)
    def fit(dd):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm('y ~ in_cn', dd, groups=dd.section).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], meth
            except Exception: pass
        return np.nan, np.nan, np.nan, 'failed'
    ST = []
    for var in VARS:
        for m in MODS:
            d = V[V.module == m].dropna(subset=[f'z_{var}'])
            for c in CNS:
                dd = d.assign(y=d[f'z_{var}'], in_cn=(d.CN == c).astype(float)); zc = dd[dd.CN == c].y.values
                b, se, p, meth = fit(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
                t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
                ST.append(dict(variant=var, module=m, CN=c, CN_label=LAB[c], CN_name=NM[c], n_sections=len(zc), mean_within_section_z=zc.mean(),
                               mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)), mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)),
                               sections_above_0=int((zc > 0).sum()), wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan, lmm_estimate=b, lmm_ci_low=b - 1.96 * se, lmm_ci_high=b + 1.96 * se, lmm_p=p, optimiser=meth))
    ST = pd.DataFrame(ST); ST['lmm_q'] = np.nan; ST['wilcoxon_q'] = np.nan
    for (var, m), g in ST.groupby(['variant', 'module']):
        for k_ in ['lmm', 'wilcoxon']:
            ok = g[f'{k_}_p'].notna(); ST.loc[g.index[ok], f'{k_}_q'] = multipletests(g[f'{k_}_p'][ok], method='fdr_bh')[1]
    ST['original_63_mean_within_section_z'] = [RELZ.loc[c, m] for c, m in zip(ST.CN, ST.module)]   # relative-score reference (column name kept)
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
    print('optimisers:', ST.optimiser.value_counts().to_dict())
    for var in VARS:
        x = ST[ST.variant == var]
        print('====', var, '(value, LMM stars | Wilcoxon stars)'); print(x.pivot(index='CN', columns='module', values='mean_within_section_z').round(2).astype(str).add(
            x.pivot(index='CN', columns='module', values='lmm_q').applymap(stars)).add(' | ').add(x.pivot(index='CN', columns='module', values='wilcoxon_q').applymap(stars)).loc[CNS, MODS].to_string())
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    SHORTM = {'Antigen presentation': 'Antigen Presentation', 'Co-stimulation / suppressive Treg': 'Co-Stimulation / Suppressive Treg',
              'Inflammatory cytokine': 'Inflammatory Cytokine', 'Matrix-integrin': 'Matrix–Integrin'}
    S = ST[ST.variant == MAIN].copy()
    # ---------------- heatmap
    E = S.pivot(index='CN', columns='module', values='mean_within_section_z').loc[CNS, MODS]; Q = S.pivot(index='CN', columns='module', values='lmm_q').loc[CNS, MODS]
    assert E.notna().all().all() and Q.notna().all().all()
    VM = float(np.ceil(np.abs(E.values).max() * 10) / 10)
    ann = np.array([[f'{E.iloc[r, c]:+.2f}{stars(Q.iloc[r, c])}' for c in range(len(MODS))] for r in range(len(CNS))])
    SHORT = ['Antigen\nPresent.', 'Co-stim /\nSuppr. Treg', 'Inflamm.\nCytokine', 'Matrix-\nIntegrin']
    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(figsize=(6.9, 5.0))
    sns.heatmap(E.set_axis([LAB[c] for c in CNS]), cmap='RdBu_r', center=0, vmin=-VM, vmax=VM, annot=ann, fmt='', annot_kws={'fontsize': 15}, ax=ax,
                linewidths=0.5, linecolor='white', cbar_kws={'label': 'Mean Within-Section z', 'shrink': 0.85, 'pad': 0.012})
    ax.set_xticklabels(SHORT, rotation=0, fontsize=14); ax.set_yticklabels([LAB[c] for c in CNS], rotation=0, fontsize=13)
    ax.set_xlabel(''); ax.set_ylabel('Spatial Neighborhood', fontsize=15)
    fig.savefig(f'{HSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    H = S.assign(module=S.module.map(SHORTM), cell_text=[ann[CNS.index(c), MODS.index(m)] for c, m in zip(S.CN, S.module)],
                 metric='niche LR signalling intensity: mean over ALL receiver cells of the CN of the mean over the cell\'s pairs of (mean sender ligand within 50 um, 0 if no sender, x own receptor); z across CNs within section, mean over sections',
                 model='z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within module')
    H.drop(columns=['variant']).apply(sig3).to_csv(f'{HSTEM}.csv', index=False)
    # ---------------- per-section dot plot
    COL = {'Antigen presentation': ('#9ecae1', '#08519c'), 'Co-stimulation / suppressive Treg': ('#bcbddc', '#54278f'),
           'Inflammatory cytokine': ('#fcae91', '#a50f15'), 'Matrix-integrin': ('#a1d99b', '#006d2c')}
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(4, 1, figsize=(7.6, 12), sharex=True, gridspec_kw={'hspace': 0.24})
    rng = np.random.default_rng(1)
    for ax, m in zip(axes, MODS):
        lc, dc = COL[m]; g = S[S.module == m].set_index('CN'); vv = V[V.module == m]
        ax.axhline(0, color='#888888', lw=1, zorder=0)
        for i, c in enumerate(CNS):
            z = vv[vv.CN == c][f'z_{MAIN}'].values
            ax.scatter(i + rng.uniform(-0.14, 0.14, len(z)), z, s=45, color=lc, edgecolor='white', linewidth=0.5, alpha=0.95, zorder=2)
            e, lo, hi, q = g.loc[c, ['mean_within_section_z', 'mean_ci_low', 'mean_ci_high', 'lmm_q']]
            ax.errorbar(i, e, yerr=[[e - lo], [hi - e]], fmt='o', ms=10, color=dc, mfc=dc if q < 0.05 else 'white', mew=2.2, elinewidth=2.6, capsize=0, zorder=3)
        vals = np.concatenate([vv[f'z_{MAIN}'].dropna().values, g.mean_ci_low.values, g.mean_ci_high.values]); ymin, ymax = vals.min(), vals.max(); span = ymax - ymin; stop = [ymax]
        for i, c in enumerate(CNS):
            q = g.loc[c, 'lmm_q']
            if q < 0.05:
                top = max(vv[vv.CN == c][f'z_{MAIN}'].max(), g.loc[c, 'mean_ci_high'])
                ax.text(i, top + 0.012 * span, stars(q), ha='center', va='bottom', fontsize=17, color=dc, fontweight='bold'); stop.append(top + 0.2 * span)
        ax.set_ylim(ymin - 0.05 * span, max(ymax + 0.05 * span, max(stop)))
        ax.set_title(SHORTM[m], fontsize=17, loc='left', pad=3, fontweight='bold'); ax.tick_params(axis='y', labelsize=15)
        ax.spines[['top', 'right']].set_visible(False); ax.set_xlim(-0.55, len(CNS) - 0.45)
    axes[-1].set_xticks(range(len(CNS))); axes[-1].set_xticklabels([LAB[c] for c in CNS], rotation=45, ha='right', rotation_mode='anchor', fontsize=13)   # the heatmap's long labels
    axes[-1].set_xlabel('Spatial Neighborhood', fontsize=16)
    fig.supylabel('Signalling Intensity (Within-Section z)', fontsize=16, x=0.04)
    h = [Line2D([], [], marker='o', ls='', color='#bdbdbd', ms=9, label='One Section'),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='#252525', ms=10, lw=2.4, label='Mean ± 95% CI, q < 0.05'),
         Line2D([], [], ls='', label=' '),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=10, mew=2, lw=2.4, label='q ≥ 0.05')]
    axes[0].legend(handles=h, loc='lower left', bbox_to_anchor=(-0.02, 1.12), ncol=2, frameon=False, fontsize=14, handletextpad=0.3, columnspacing=0.8, borderaxespad=0)
    fig.savefig(f'{SSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    SUP = pd.concat([ST.assign(row_type='summary'), V.melt(id_vars=['section', 'CN', 'module', 'n_cells', 'n_receiver_cells', 'n_receivers_with_sender'],
                     value_vars=[f'z_{v}' for v in VARS] + VARS, var_name='quantity', value_name='value').assign(row_type='section')], ignore_index=True)
    SUP.apply(sig3).to_csv(f'{SSTEM}.csv', index=False)
    # ---------------- legends
    npairs = KP.groupby('module').size()
    common = ('CD resections (CosMx 6K; 11 sections, 4 patients; all cells). Ligand-receptor pairs were taken from the NicheNet v2 prior network, '
          'matching Fig 3, and grouped into functional modules using CellChatDB pathway annotations (antigen presentation = MHC-II; co-stimulation / '
          'suppressive Treg = CD80, CD86, PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2, IL10; inflammatory cytokine = IL1, IL6, TNF, IFN-II, IL12, '
          f'IL17, OSM; matrix–integrin = all ECM-Receptor pathways; {npairs[MODS[0]]}, {npairs[MODS[1]]}, {npairs[MODS[2]]} and {npairs[MODS[3]]} pairs '
          'with both genes on the 6K panel). Each pair is restricted to the cell types that send and receive it, assigned from established expression '
          'and erring wide (fibroblasts = OGN+RSPO3+, ADAMDEC1+, activated, crypt-top and PLA2G2A+ ECM fibroblasts and myofibroblasts; epithelial = '
          'enterocytes, absorptive, secretory, microfold, stem and cycling cells; macrophages = CD9+ M2, SELENOP+ and Mo-Macs; SELENOP+ fibroblasts never used). Senders: MHC-II from '
          'DC, B cells and macrophages; APC co-stimulatory ligands (CD80, CD86, ICOSLG, CD274, PDCD1LG2, CD70, TNFSF4, TNFSF18) from DC '
          'and B cells; IL2, IL4, IL9, IL21 and CD40LG from CD4 effector, CD4 Trm and Tfh cells, with IL7 also from epithelial cells, fibroblasts and '
          'endothelial cells, IL15 also from DC, macrophages, epithelial cells and fibroblasts, and TSLP also from epithelial cells and '
          'fibroblasts; IL10-family cytokines from DC and macrophages; IL1-family, IL12-family and OSM from macrophages, '
          'DC and neutrophils, with IL33 also from epithelial cells, endothelial cells and fibroblasts, IL36A/B/G also from epithelial cells, '
          'IL6 also from fibroblasts and endothelial cells, IL11 also from fibroblasts, and TNF also from T and NK cells; IFNG and IL17-family from CD4 '
          'effector, CD8 Trm, NK and γδ T cells, with IL17B, IL17C and IL25 also from epithelial cells; extracellular-matrix ligands from fibroblasts. '
          'Receivers: CD4 T cells (CD4 effector, CD4 Trm, Tfh, exhausted CD4 Tfh, resting T, Treg) for antigen presentation; all T cells and NK cells '
          'for APC co-stimulatory ligands and IL2-family cytokines; DC, B cells and macrophages for CD40LG–CD40, and myeloid, T and NK '
          'cells for CD40LG–integrins; all immune cells (myeloid, B, plasma, T, NK) for IL10; epithelial cells and fibroblasts for IL19, IL20, IL22 and '
          'IL24; epithelial, myeloid, T and NK cells for IL26; all cell types except SELENOP+ fibroblasts for inflammatory cytokines and matrix–integrin '
          '(receptor expression gates the score). '
          'For each receiver cell and pair, pair activity is the mean log-normalised ligand expression of the pair\'s sender cells within 50 µm '
          '(0 when no sender cell is within 50 µm) multiplied by the receiver\'s own log-normalised receptor expression; a receiver cell\'s module '
          'intensity is the mean over the module\'s pairs it receives, and the neighbourhood value in a section is the mean over all its receiver cells '
          '(non-receivers excluded). Receivers without a sender within 50 µm contribute zero, avoiding selection of border cells in sender-poor '
          'neighbourhoods; sender presence counts but sender number does not. The score therefore reflects per-cell signalling intensity (mean sender ligand '
          '× receiver receptor), independent of how many sender or receiver cells a neighbourhood contains (neighbourhood × section combinations with '
          'at least 100 cells and at least 20 receiver cells). Values are z-standardised across neighbourhoods within each section and averaged over '
          'sections ("Mean Within-Section z"; 0 = neighbourhood average). Significance: linear mixed model on the per-section values, '
          'z ~ in neighbourhood + (1 | section), BH-corrected over the 10 neighbourhoods within each module; * q < 0.05, ** q < 0.01, *** q < 0.001; '
          'no star, q ≥ 0.05. ')
    legH = ('Ligand-receptor module activity across spatial cellular neighbourhoods (rows) by module (columns), ' + common
            + f'Colour scale symmetric, ±{VM:.1f}. Section-level values are shown in the companion dot plot. The source data hold the per-section values, '
        'the mixed-model estimates and q values, a Wilcoxon signed-rank test of the per-section values (BH within module), and, as sensitivity '
        'analyses, a version restricted to receivers with at least one sender within 50 µm (and at least 20 such receivers), and an abundance-weighted version (summed sender ligand, summed over the neighbourhood\'s cells and divided by its number of cells, '
        'non-receivers counted as 0) with and without a minimum of 3 sender cells per pair. Neighbourhood names are abbreviated by their dominant cell types; full names are in the source data.')
    legS = ('Section-level ligand-receptor module activity across spatial cellular neighbourhoods (companion to the Fig 4 niche heatmap), ' + common
            + 'Light dots: one section each. Dark points and bars: mean over sections ± 95% CI (t distribution); filled, q < 0.05. Neighbourhood names are '
        'abbreviated by their dominant cell types; full names are in the source data.')
    for t in (legH, legS): assert 'activation score' not in t.lower() and 'quiescent' not in t.lower() and 'to tregs' not in t.lower()
    open(f'{HSTEM}_legend.txt', 'w').write(legH + '\n'); open(f'{SSTEM}_legend.txt', 'w').write(legS + '\n')
    print('colour limits +/-', VM); print('done', flush=True)


# ======================================================================================================================
# Section cn_lr_modules_allcells: Cell-type-agnostic LR module scores over all cells, per CN and section
# Panels: Fig 5C (right) and Supp 9A-B inputs (module membership, per-cell cache, per-section CN values); Supp Table 10 (membership)
# ======================================================================================================================
if _run('cn_lr_modules_allcells'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/82_lr_modules_allcells.py'   # working folder of this step (relative paths resolve from here)
    """Ligand-receptor module activity per cellular neighbourhood from cell-type-agnostic LIANA spatial scores over all cells, with abundance included (raw
local scores, not the density-normalised relative score of cn_module_relative_scores).
Data: per-cell LIANA bivariate local scores, combined/spot_lr_df.parquet (NicheNet lr_network pairs; Gaussian kernel, 200 px bandwidth);
  3,002 pairs x 2,948,012 cells; values in [0, 1], NaN -> 0.
Modules: module_of() as in cn_module_relative_scores (CellChatDB pathway labels: Antigen presentation = MHC-II; Co-stimulation / suppressive Treg = CD80,
  CD86, PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2, IL10; Inflammatory cytokine = IL1, IL6, TNF, IFN-II, IL12, IL17, OSM; Matrix-integrin =
  ECM-Receptor pathways).
Main ('raw'): all cells of the parquet with a CN (no subsampling). Cell module score = mean raw local score over the module's pairs; CN x section value =
  mean over all cells of that CN in that section (cells without interacting partners contribute their low scores, so abundance is included); CN x section
  kept if >= 100 cells; z across CNs within section, mean over sections ('Mean Within-Section z'); stars: linear mixed model on the per-section z,
  z ~ in_CN + (1 | section), BH over the 10 CNs within module (Wilcoxon signed-rank q in the CSV).
Sensitivity and reference columns (CSV only):
  'relative_63'      = the relative score of cn_module_relative_scores (random 300,000 cells, seed 0; module mean / mean over all 3,002 pairs of the cell;
                       z within section over cells; mean per CN), recomputed from the per-cell cache and asserted equal to module_by_CN_pooled_vs_adjusted.csv;
  'relative_cnmean'  = relative score, all cells, aggregated like the main metric (isolates raw vs relative);
  niche_*            = the cell-type-restricted values of cn_lr_modules_celltype_restricted.
Outputs: A_fig4_panels/supp_lr_modules_allcells_heatmap.{png,csv} + _legend.txt (read by the later module sections); A_fig4_panels/fig4_lr_modules_allcells_dotplot.{png,csv}
  + _legend.txt; B_resection_analyses/lr_modules_cellchat/lr_modules_allcells_per_cell.parquet (per-cell cache), lr_modules_allcells_membership.csv
  (module membership) and lr_modules_allcells_per_section_cn.csv (CN x section values).
SKIP_COMPUTE=1 redraws from the per-cell cache."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.lines import Line2D
    from scipy import stats as sst
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; SPOT = f'{R}/combined/spot_lr_df.parquet'
    HSTEM = f'{OUTA}/supp_lr_modules_allcells_heatmap'; SSTEM = f'{OUTA}/fig4_lr_modules_allcells_dotplot'; CACHE = f'{OUTB}/lr_modules_allcells_per_cell.parquet'
    MIN_CELLS_CN = 100; MAIN = 'raw'
    # ---------------- modules as in cn_module_relative_scores
    MOD = {'Antigen presentation': ['MHC-II'], 'Co-stimulation / suppressive Treg': ['CD80', 'CD86', 'PD-L1', 'PDL2', 'ICOS', 'OX40', 'GITRL', 'CD70', 'CD40', 'IL2', 'IL10'],
           'Inflammatory cytokine': ['IL1', 'IL6', 'TNF', 'IFN-II', 'IL12', 'IL17', 'OSM']}
    cc = pd.read_csv(f'{D}/objects/cellchatdb_human_interactions_snapshot.csv')
    MOD['Matrix-integrin'] = sorted(cc.loc[cc.annotation == 'ECM-Receptor', 'pathway_name'].unique())
    MODS = list(MOD)
    gs = lambda x: frozenset(y.strip() for y in str(x).replace(';', ',').split(',') if y.strip() and y.strip() != 'nan')
    cc['L'] = cc['ligand.symbol'].map(gs); cc['R'] = cc['receptor.symbol'].map(gs)
    cols = [c for c in pq.read_schema(SPOT).names if '^' in c]
    def module_of(c):
        l, r = c.split('^'); L, Rr = frozenset(l.split('_')), frozenset(r.split('_')); hits = []
        for m, pws in MOD.items():
            sub = cc[cc.pathway_name.isin(pws)]
            if any((L == a) and (Rr == b or Rr <= b or b <= Rr) for a, b in zip(sub.L, sub.R)): hits.append(m)
        return hits
    MEM = {c: module_of(c) for c in cols}; MEM = {c: m for c, m in MEM.items() if m}
    memb = pd.DataFrame([(c, m) for c, ms in MEM.items() for m in ms], columns=['lr_pair', 'module']); memb.to_csv(f'{OUTB}/lr_modules_allcells_membership.csv', index=False)
    NPAIR = memb.module.value_counts().reindex(MODS); print(NPAIR.to_string(), flush=True)
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
    idx_all = pq.read_table(SPOT, columns=['__index_level_0__']).column(0).to_pylist(); pos = pd.Series(np.arange(len(idx_all)), index=idx_all)
    # ---------------- per-cell module scores over ALL cells (raw and relative), cached
    if os.environ.get('SKIP_COMPUTE') != '1':
        n = len(idx_all); tot = np.zeros(n); acc = {m: np.zeros(n) for m in MOD}; npair = {m: 0 for m in MOD}
        for i in range(0, len(cols), 250):
            ch = cols[i:i + 250]; A = pq.read_table(SPOT, columns=ch).to_pandas().astype(np.float32).fillna(0.0)
            tot += A.values.sum(1)
            for m in MOD:
                ps = [c for c in ch if c in MEM and m in MEM[c]]
                if ps: acc[m] += A[ps].values.sum(1); npair[m] += len(ps)
            del A; gc.collect(); print('columns', i + len(ch), '/', len(cols), flush=True)
        assert all(npair[m] == NPAIR[m] for m in MOD)
        den = tot / len(cols); den[den == 0] = np.nan
        C = pd.DataFrame({**{f'raw|{m}': acc[m] / npair[m] for m in MOD}, **{f'rel|{m}': (acc[m] / npair[m]) / den for m in MOD}}, index=pd.Index(idx_all, name='cell'))
        C[[f'rel|{m}' for m in MOD]] = C[[f'rel|{m}' for m in MOD]].fillna(0.0)
        C.astype(np.float32).to_parquet(CACHE)
    C = pd.read_parquet(CACHE)
    # ---------------- (1) relative score of cn_module_relative_scores from the per-cell cache (random 300,000, seed 0, per-cell z within section, mean per CN)
    rng = np.random.default_rng(0); pool = np.asarray(cn.index[cn.index.isin(pos.index)])
    samp = pd.Index(rng.choice(pool, 300000, replace=False)); take = np.sort(pos.loc[samp].values); rows_idx = np.asarray(idx_all)[take]
    S = C.loc[rows_idx, [f'rel|{m}' for m in MODS]].set_axis(MODS, axis=1).astype(np.float64)
    S['CN'] = cn.CN_k10.reindex(S.index).values; S['section'] = cn.section.reindex(S.index).values
    Z = S.groupby('section')[MODS].transform(lambda x: (x - x.mean()) / (x.std() or 1))
    REP = Z.assign(CN=S.CN.values).groupby('CN')[MODS].mean()
    REF63 = pd.read_csv(f'{OUTB}/module_by_CN_pooled_vs_adjusted.csv', header=[0, 1], index_col=0)['section_adjusted_z'][MODS]
    REF63.index = [int(i.split(':')[0][2:]) for i in REF63.index]
    print('relative_63 recomputed vs module_by_CN_pooled_vs_adjusted.csv:'); print(pd.concat({'recomputed': REP.round(3), 'table': REF63}, axis=1).to_string())
    assert np.allclose(REP.loc[REF63.index].values, REF63.values, rtol=5e-3, atol=1e-3), 'relative_63 differs from module_by_CN_pooled_vs_adjusted.csv'
    for (c, m, v) in [(4, 'Antigen presentation', 0.43), (8, 'Co-stimulation / suppressive Treg', 0.743), (9, 'Matrix-integrin', 0.844), (10, 'Antigen presentation', 0.186)]:
        assert abs(REP.loc[c, m] - v) < 2e-3, (c, m, REP.loc[c, m], v)
    # ---------------- (2) CN x section means over ALL cells (main raw; relative_cnmean)
    X = C.join(cn[['CN_k10', 'section']], how='inner').rename(columns={'CN_k10': 'CN'})
    print('cells with LIANA scores:', len(C), '| with a CN:', len(X), '| CN cells without LIANA scores:', int((~cn.index.isin(C.index)).sum()), flush=True)
    G = X.groupby(['section', 'CN'])
    V = (G[[f'raw|{m}' for m in MODS]].mean().set_axis(MODS, axis=1).reset_index().melt(id_vars=['section', 'CN'], var_name='module', value_name='raw')
         .merge(G[[f'rel|{m}' for m in MODS]].mean().set_axis(MODS, axis=1).reset_index().melt(id_vars=['section', 'CN'], var_name='module', value_name='relative_cnmean'),
                on=['section', 'CN', 'module'])
         .merge(G.size().rename('n_cells').reset_index(), on=['section', 'CN']))
    V.to_csv(f'{OUTB}/lr_modules_allcells_per_section_cn.csv', index=False, float_format='%.6g')
    CNS = list(range(1, 11)); VARS = ['raw', 'relative_cnmean']
    V = V[V.n_cells >= MIN_CELLS_CN].copy()
    zs = lambda x: (x - x.mean()) / x.std()
    for var in VARS: V[f'z_{var}'] = V.groupby(['section', 'module'])[var].transform(zs)
    def fit(dd):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm('y ~ in_cn', dd, groups=dd.section).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], meth
            except Exception: pass
        return np.nan, np.nan, np.nan, 'failed'
    # ---------------- CN names / labels (as module_by_CN_pooled_vs_adjusted.csv)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [(' Fib / PLA2G2A⁺ Fib / Endothelial', '/PLA2G2A⁺ Fib'), (' / Crypt Top Fib / Cycling', ' / CryptTop / Cycl'),
            (' / ADAMDEC1 Fib', ' / ADAMDEC1'), (' / Neuron / Pericyte', ' / Neuron / Peri'), (' / Mac S_M_P / CD9 Mac', ' / Mac / CD9 Mac'),
            ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'), ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'),
            ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {c: short_cn(NM[c]) for c in CNS}
    assert [LAB[c] for c in CNS] == [short_cn(i) for i in pd.read_csv(f'{OUTB}/module_by_CN_pooled_vs_adjusted.csv', header=[0, 1], index_col=0).index]
    ST = []
    for var in VARS:
        for m in MODS:
            d = V[V.module == m].dropna(subset=[f'z_{var}'])
            for c in CNS:
                dd = d.assign(y=d[f'z_{var}'], in_cn=(d.CN == c).astype(float)); zc = dd[dd.CN == c].y.values
                b, se, p, meth = fit(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
                t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
                ST.append(dict(variant=var, module=m, CN=c, CN_label=LAB[c], CN_name=NM[c], n_sections=len(zc), mean_within_section_z=zc.mean(),
                               mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)), mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)),
                               sections_above_0=int((zc > 0).sum()), wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan,
                               lmm_estimate=b, lmm_ci_low=b - 1.96 * se, lmm_ci_high=b + 1.96 * se, lmm_p=p, optimiser=meth))
    ST = pd.DataFrame(ST); ST['lmm_q'] = np.nan; ST['wilcoxon_q'] = np.nan
    for (var, m), g in ST.groupby(['variant', 'module']):
        for k_ in ['lmm', 'wilcoxon']:
            ok = g[f'{k_}_p'].notna(); ST.loc[g.index[ok], f'{k_}_q'] = multipletests(g[f'{k_}_p'][ok], method='fdr_bh')[1]
    ST['relative_63_original_reproduced'] = [REP.loc[c, m] for c, m in zip(ST.CN, ST.module)]
    ST['relative_63_original_archived'] = [REF63.loc[c, m] for c, m in zip(ST.CN, ST.module)]   # column names kept
    # reference: cell-type-restricted values (cn_lr_modules_celltype_restricted)
    NI = pd.read_csv(f'{OUTA}/supp_lr_modules_niche.csv'); NI = NI[NI.row_type == 'summary']
    for v_ in NI.variant.unique():
        x = NI[NI.variant == v_].set_index(['module', 'CN'])
        ST[f'niche_{v_}_z'] = [x.loc[(m, c), 'mean_within_section_z'] for m, c in zip(ST.module, ST.CN)]
        ST[f'niche_{v_}_q'] = [x.loc[(m, c), 'lmm_q'] for m, c in zip(ST.module, ST.CN)]
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
    print('optimisers:', ST.optimiser.value_counts().to_dict())
    for var in VARS:
        x = ST[ST.variant == var]
        print('====', var, '(value, LMM stars | Wilcoxon stars); significant (LMM):', int((x.lmm_q < 0.05).sum()), '/ 40')
        print(x.pivot(index='CN', columns='module', values='mean_within_section_z').round(2).astype(str).add(
            x.pivot(index='CN', columns='module', values='lmm_q').applymap(stars)).add(' | ').add(x.pivot(index='CN', columns='module', values='wilcoxon_q').applymap(stars)).loc[CNS, MODS].to_string())
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    SHORTM = {'Antigen presentation': 'Antigen Presentation', 'Co-stimulation / suppressive Treg': 'Co-Stimulation / Suppressive Treg',
              'Inflammatory cytokine': 'Inflammatory Cytokine', 'Matrix-integrin': 'Matrix–Integrin'}
    S_ = ST[ST.variant == MAIN].copy()
    # ---------------- heatmap
    E = S_.pivot(index='CN', columns='module', values='mean_within_section_z').loc[CNS, MODS]; Q = S_.pivot(index='CN', columns='module', values='lmm_q').loc[CNS, MODS]
    assert E.notna().all().all() and Q.notna().all().all()
    VM = float(np.ceil(np.abs(E.values).max() * 10) / 10)
    ann = np.array([[f'{E.iloc[r, c]:+.2f}{stars(Q.iloc[r, c])}' for c in range(len(MODS))] for r in range(len(CNS))])
    SHORT = ['Antigen\nPresent.', 'Co-stim /\nSuppr. Treg', 'Inflamm.\nCytokine', 'Matrix-\nIntegrin']
    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(figsize=(6.9, 5.0))
    sns.heatmap(E.set_axis([LAB[c] for c in CNS]), cmap='RdBu_r', center=0, vmin=-VM, vmax=VM, annot=ann, fmt='', annot_kws={'fontsize': 15}, ax=ax,
                linewidths=0.5, linecolor='white', cbar_kws={'label': 'Mean Within-Section z', 'shrink': 0.85, 'pad': 0.012})
    ax.set_xticklabels(SHORT, rotation=0, fontsize=14); ax.set_yticklabels([LAB[c] for c in CNS], rotation=0, fontsize=13)
    ax.set_xlabel(''); ax.set_ylabel('Spatial Neighborhood', fontsize=15)
    fig.savefig(f'{HSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    H = S_.assign(module=S_.module.map(SHORTM), cell_text=[ann[CNS.index(c), MODS.index(m)] for c, m in zip(S_.CN, S_.module)],
                  metric='LIANA bivariate local score (raw), mean over module pairs per cell, mean over ALL cells of the CN per section; z across CNs within section, mean over sections',
                  model='z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within module')
    H.drop(columns=['variant']).apply(sig3).to_csv(f'{HSTEM}.csv', index=False)
    # ---------------- per-section dot plot (long labels)
    COL = {'Antigen presentation': ('#9ecae1', '#08519c'), 'Co-stimulation / suppressive Treg': ('#bcbddc', '#54278f'),
           'Inflammatory cytokine': ('#fcae91', '#a50f15'), 'Matrix-integrin': ('#a1d99b', '#006d2c')}
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(4, 1, figsize=(7.6, 12), sharex=True, gridspec_kw={'hspace': 0.24})
    rng2 = np.random.default_rng(1)
    for ax, m in zip(axes, MODS):
        lc, dc = COL[m]; g = S_[S_.module == m].set_index('CN'); vv = V[V.module == m]
        ax.axhline(0, color='#888888', lw=1, zorder=0)
        for i, c in enumerate(CNS):
            z = vv[vv.CN == c][f'z_{MAIN}'].values
            ax.scatter(i + rng2.uniform(-0.14, 0.14, len(z)), z, s=45, color=lc, edgecolor='white', linewidth=0.5, alpha=0.95, zorder=2)
            e, lo, hi, q = g.loc[c, ['mean_within_section_z', 'mean_ci_low', 'mean_ci_high', 'lmm_q']]
            ax.errorbar(i, e, yerr=[[e - lo], [hi - e]], fmt='o', ms=10, color=dc, mfc=dc if q < 0.05 else 'white', mew=2.2, elinewidth=2.6, capsize=0, zorder=3)
        vals = np.concatenate([vv[f'z_{MAIN}'].dropna().values, g.mean_ci_low.values, g.mean_ci_high.values]); ymin, ymax = vals.min(), vals.max(); span = ymax - ymin; stop = [ymax]
        for i, c in enumerate(CNS):
            q = g.loc[c, 'lmm_q']
            if q < 0.05:
                top = max(vv[vv.CN == c][f'z_{MAIN}'].max(), g.loc[c, 'mean_ci_high'])
                ax.text(i, top + 0.012 * span, stars(q), ha='center', va='bottom', fontsize=17, color=dc, fontweight='bold'); stop.append(top + 0.2 * span)
        ax.set_ylim(ymin - 0.05 * span, max(ymax + 0.05 * span, max(stop)))
        ax.set_title(SHORTM[m], fontsize=17, loc='left', pad=3, fontweight='bold'); ax.tick_params(axis='y', labelsize=15)
        ax.spines[['top', 'right']].set_visible(False); ax.set_xlim(-0.55, len(CNS) - 0.45)
    axes[-1].set_xticks(range(len(CNS))); axes[-1].set_xticklabels([LAB[c] for c in CNS], rotation=45, ha='right', rotation_mode='anchor', fontsize=13)
    axes[-1].set_xlabel('Spatial Neighborhood', fontsize=16)
    fig.supylabel('LR Module Score (Within-Section z)', fontsize=16, x=0.04)
    h = [Line2D([], [], marker='o', ls='', color='#bdbdbd', ms=9, label='One Section'),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='#252525', ms=10, lw=2.4, label='Mean ± 95% CI, q < 0.05'),
         Line2D([], [], ls='', label=' '),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=10, mew=2, lw=2.4, label='q ≥ 0.05')]
    axes[0].legend(handles=h, loc='lower left', bbox_to_anchor=(-0.02, 1.12), ncol=2, frameon=False, fontsize=14, handletextpad=0.3, columnspacing=0.8, borderaxespad=0)
    fig.savefig(f'{SSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    SUP = pd.concat([ST.assign(row_type='summary'), V.melt(id_vars=['section', 'CN', 'module', 'n_cells'], value_vars=[f'z_{v}' for v in VARS] + VARS,
                     var_name='quantity', value_name='value').assign(row_type='section')], ignore_index=True)
    SUP.apply(sig3).to_csv(f'{SSTEM}.csv', index=False)
    # ---------------- legends
    common = (f'CD resections (CosMx 6K; 11 sections, 4 patients; all {len(X):,} cells). Ligand-receptor (LR) scores are cell-type-agnostic spatial LIANA '
          'bivariate local scores (Gaussian kernel, 200 px bandwidth) for the NicheNet v2 prior-network pairs (as Fig 3), computed for every cell; pairs '
          'are grouped into modules by CellChatDB pathway annotation (antigen presentation = MHC-II; co-stimulation / suppressive Treg = CD80, CD86, '
          'PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2, IL10; inflammatory cytokine = IL1, IL6, TNF, IFN-II, IL12, IL17, OSM; matrix–integrin = '
          f'ECM-Receptor pathways; {NPAIR[MODS[0]]}, {NPAIR[MODS[1]]}, {NPAIR[MODS[2]]} and {NPAIR[MODS[3]]} LR pairs). No sender or receiver cell types '
          'are imposed. A cell\'s module score is the mean raw local score over the module\'s pairs, and a neighbourhood\'s value in a section is the '
          'mean over all its cells, so the score reflects both ligand–receptor expression and the abundance of interacting cells (neighbourhood × '
          'section combinations with at least 100 cells). Values are z-standardised across neighbourhoods within each section and averaged over sections '
          '("Mean Within-Section z"; 0 = neighbourhood average). Significance: linear mixed model on the per-section values, z ~ in neighbourhood + '
          '(1 | section), BH-corrected over the 10 neighbourhoods within each module; * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. ')
    legH = ('Spatial ligand-receptor module activity across cellular neighbourhoods (rows) by module (columns), ' + common
            + f'Colour scale symmetric, ±{VM:.1f}. Section-level values are shown in the companion dot plot. The source data hold the per-section values, '
        'the mixed-model estimates and q values, a Wilcoxon signed-rank test (BH within module), and, as sensitivity analyses, the density-normalised '
        '(relative) score and the cell-type-restricted versions. Neighbourhood names are abbreviated by their dominant cell '
        'types; full names are in the source data.')
    legS = ('Section-level spatial ligand-receptor module activity across cellular neighbourhoods (companion to the Fig 4 heatmap), ' + common
            + 'Light dots: one section each. Dark points and bars: mean over sections ± 95% CI (t distribution); filled, q < 0.05. Neighbourhood names are '
        'abbreviated by their dominant cell types; full names are in the source data.')
    for t in (legH, legS): assert 'activation score' not in t.lower() and 'quiescent' not in t.lower() and 'to tregs' not in t.lower()
    open(f'{HSTEM}_legend.txt', 'w').write(legH + '\n'); open(f'{SSTEM}_legend.txt', 'w').write(legS + '\n')
    print('colour limits +/-', VM); print('done', flush=True)


# ======================================================================================================================
# Section cn10_lr_pairs: Pair-level LR scores per CN and section, with CN10 statistics per pair
# Panels: per-pair cache read by lr_costim_inhibitory_modules, cn_lr_modules_immunomodulatory, cn_lr_modules_costim and cn_inflammatory_axes
# ======================================================================================================================
if _run('cn10_lr_pairs'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/83_cn10_lr_pairs.py'   # working folder of this step (relative paths resolve from here)
    """Pair-level breakdown of the four ligand-receptor modules for CN10 (OGN+RSPO3+ / PLA2G2A+ Fib). Module averages (e.g. Co-Stimulation, 50 pairs) can
hide single pairs that are higher or lower in CN10.
Statistic as in cn_lr_modules_allcells, applied to each pair instead of the module mean:
  cell-type-agnostic LIANA bivariate local score (raw, NaN -> 0) of the pair for every cell of combined/spot_lr_df.parquet;
  CN x section value = mean over all cells of that CN in that section (CN x section kept if >= 100 cells);
  z across the 10 CNs within section (sections where the pair is constant, e.g. all zero in 18888_A10, give undefined z and drop out);
  CN10: mean z over sections (+ t 95% CI), n sections, sections above 0, Wilcoxon signed-rank p, and the linear mixed model
  z ~ in_CN10 + (1 | section) (estimate, 95% CI, p; optimiser fallback); BH within module (LMM and Wilcoxon separately).
  CN4 and CN8 mean z given for context.
Module membership: lr_modules_allcells_membership.csv (cn_lr_modules_allcells); CN: cn_per_cell.parquet, CN_k10.
The mean over a module's pairs of the per-pair CN x section raw means equals the module values of cn_lr_modules_allcells (linearity; asserted).
Outputs: A_fig4_panels/supp_cn10_lr_pairs.{png,csv} + _legend.txt; B_resection_analyses/lr_modules_cellchat/lr_pairs_allcells_per_section_cn.csv (per-pair cache).
SKIP_COMPUTE=1 redraws from the cache."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scipy import stats as sst, sparse
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; SPOT = f'{R}/combined/spot_lr_df.parquet'
    STEM = f'{OUTA}/supp_cn10_lr_pairs'; CACHE = f'{OUTB}/lr_pairs_allcells_per_section_cn.csv'
    MIN_CELLS_CN = 100; TARGET = 10; CONTEXT = [4, 8]; CNS = list(range(1, 11))
    MODS = ['Antigen presentation', 'Co-stimulation / suppressive Treg', 'Inflammatory cytokine', 'Matrix-integrin']
    SHORTM = {'Antigen presentation': 'Antigen Presentation', 'Co-stimulation / suppressive Treg': 'Co-Stimulation / Suppressive Treg',
              'Inflammatory cytokine': 'Inflammatory Cytokine', 'Matrix-integrin': 'Matrix–Integrin'}
    memb = pd.read_csv(f'{OUTB}/lr_modules_allcells_membership.csv')
    assert not memb.lr_pair.duplicated().any()
    NPAIR = memb.module.value_counts().reindex(MODS); assert NPAIR.tolist() == [9, 50, 34, 227], NPAIR
    PAIRS = memb.lr_pair.tolist(); MODOF = memb.set_index('lr_pair').module.to_dict()
    # ---------------- per pair x section x CN means over ALL cells (as cn_lr_modules_allcells: inner join of LIANA cells with CN cells)
    if os.environ.get('SKIP_COMPUTE') != '1':
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas()
        assert not cn.index.duplicated().any()
        idx_all = pq.read_table(SPOT, columns=['__index_level_0__']).column(0).to_pylist()
        rows = pd.Index(idx_all); assert not rows.duplicated().any()
        keep = rows.isin(cn.index); print('LIANA cells:', len(rows), '| with a CN:', int(keep.sum()), flush=True)
        sec = cn.section.reindex(rows[keep]).values; cnn = cn.CN_k10.reindex(rows[keep]).values.astype(int)
        grp = pd.MultiIndex.from_arrays([sec, cnn], names=['section', 'CN']); codes, uniq = pd.factorize(grp)
        M = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
        ncell = np.asarray(M.sum(1)).ravel(); kidx = np.where(keep)[0]; out = []
        for i in range(0, len(PAIRS), 64):
            ch = PAIRS[i:i + 64]; A = pq.read_table(SPOT, columns=ch).to_pandas().astype(np.float64).fillna(0.0).values[kidx]
            out.append(pd.DataFrame(np.asarray(M @ A) / ncell[:, None], columns=ch)); del A; gc.collect(); print('pairs', i + len(ch), '/', len(PAIRS), flush=True)
        W = pd.concat(out, axis=1); W.index = pd.MultiIndex.from_tuples(list(uniq), names=['section', 'CN']); W['n_cells'] = ncell
        L = W.reset_index().melt(id_vars=['section', 'CN', 'n_cells'], var_name='lr_pair', value_name='raw')
        L['module'] = L.lr_pair.map(MODOF); L.to_csv(CACHE, index=False, float_format='%.8g')
    L = pd.read_csv(CACHE)
    # check against cn_lr_modules_allcells: module value = mean over pairs (linearity)
    V82 = pd.read_csv(f'{OUTB}/lr_modules_allcells_per_section_cn.csv')
    chk = L.groupby(['section', 'CN', 'module']).agg(raw=('raw', 'mean'), n=('n_cells', 'first')).reset_index().merge(V82, on=['section', 'CN', 'module'], suffixes=('', '_82'))
    assert len(chk) == len(V82) and (chk.n == chk.n_cells).all() and np.allclose(chk.raw, chk.raw_82, rtol=1e-4, atol=1e-7), 'per-pair means do not reproduce cn_lr_modules_allcells'
    print('matches cn_lr_modules_allcells module values:', len(chk), 'section x CN x module rows; max abs diff', float((chk.raw - chk.raw_82).abs().max()), flush=True)
    # ---------------- z across CNs within section, per pair
    L = L[L.n_cells >= MIN_CELLS_CN].copy()
    L['z'] = L.groupby(['section', 'lr_pair']).raw.transform(lambda x: (x - x.mean()) / x.std())
    L.loc[~np.isfinite(L.z), 'z'] = np.nan
    def fit(dd):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm('y ~ in_cn', dd, groups=dd.section).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], meth
            except Exception: pass
        return np.nan, np.nan, np.nan, 'failed'
    ST = []
    for p, d in L.dropna(subset=['z']).groupby('lr_pair'):
        dd = d.assign(y=d.z, in_cn=(d.CN == TARGET).astype(float)); zc = dd[dd.CN == TARGET].y.values
        b, se, pv, meth = fit(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
        t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
        lig, rec = p.split('^')
        r = dict(module=MODOF[p], lr_pair=p, label=f"{lig.replace('_', '+')} → {rec.replace('_', '+')}", ligand=lig.replace('_', '+'), receptor=rec.replace('_', '+'),
                 n_sections_total=d.section.nunique(), cn10_n_sections=len(zc), cn10_mean_within_section_z=zc.mean() if len(zc) else np.nan,
                 cn10_mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)) if len(zc) > 1 else np.nan,
                 cn10_mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)) if len(zc) > 1 else np.nan,
                 cn10_sections_above_0=int((zc > 0).sum()), cn10_wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan,
                 cn10_lmm_estimate=b, cn10_lmm_se=se, cn10_lmm_ci_low=b - 1.96 * se, cn10_lmm_ci_high=b + 1.96 * se, cn10_lmm_p=pv, optimiser=meth)
        for c in CONTEXT:
            zz = d[d.CN == c].z.values; r[f'cn{c}_mean_within_section_z'] = zz.mean() if len(zz) else np.nan; r[f'cn{c}_n_sections'] = len(zz)
        r['cn10_mean_raw_score'] = L[(L.lr_pair == p) & (L.CN == TARGET)].raw.mean(); r['all_cn_mean_raw_score'] = L[L.lr_pair == p].raw.mean()
        ST.append(r)
    ST = pd.DataFrame(ST)
    missing = sorted(set(PAIRS) - set(ST.lr_pair)); print('pairs with no defined z in any section:', missing)
    if missing: ST = pd.concat([ST, pd.DataFrame([dict(module=MODOF[p], lr_pair=p, label=f"{p.split('^')[0].replace('_', '+')} → {p.split('^')[1].replace('_', '+')}",
                                                       cn10_n_sections=0, optimiser='not fitted') for p in missing])], ignore_index=True)
    for k_ in ['lmm', 'wilcoxon']:
        ST[f'cn10_{k_}_q'] = np.nan
        for m, g in ST.groupby('module'):
            ok = g[f'cn10_{k_}_p'].notna(); ST.loc[g.index[ok], f'cn10_{k_}_q'] = multipletests(g[f'cn10_{k_}_p'][ok], method='fdr_bh')[1]
    ST['module'] = pd.Categorical(ST.module, MODS); ST = ST.sort_values(['module', 'cn10_lmm_estimate'], ascending=[True, False]).reset_index(drop=True)
    ST['module_label'] = ST.module.astype(str).map(SHORTM)
    # module-level reference (cn_lr_modules_allcells, CN10)
    M82 = pd.read_csv(f'{OUTA}/supp_lr_modules_allcells_heatmap.csv'); M82 = M82[M82.CN == TARGET].set_index('module')
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 400)
    print('optimisers:', ST.optimiser.value_counts().to_dict())
    SHOWC = ['label', 'cn10_lmm_estimate', 'cn10_lmm_ci_low', 'cn10_lmm_ci_high', 'cn10_lmm_q', 'cn10_wilcoxon_q', 'cn10_mean_within_section_z', 'cn10_n_sections',
             'cn10_sections_above_0', 'cn4_mean_within_section_z', 'cn8_mean_within_section_z']
    for m in MODS:
        g = ST[ST.module == m]; mm = M82.loc[SHORTM[m]]
        print(f'==== {m}: {len(g)} pairs; up q<0.05: {int(((g.cn10_lmm_q < 0.05) & (g.cn10_lmm_estimate > 0)).sum())}, down q<0.05: '
          f'{int(((g.cn10_lmm_q < 0.05) & (g.cn10_lmm_estimate < 0)).sum())}; q<0.1: {int((g.cn10_lmm_q < 0.1).sum())} | module-level CN10 z {mm.mean_within_section_z:+.3f}, q {mm.lmm_q:.3g}')
        print(g[SHOWC].round(3).to_string(index=False) if m != 'Matrix-integrin' else g[g.cn10_lmm_q < 0.05][SHOWC].round(3).to_string(index=False))
    # ---------------- figure selection
    def select(g):
        if g.module.iloc[0] == 'Co-stimulation / suppressive Treg': return g, 'all pairs'
        if len(g) <= 12: return g, 'all pairs'
        s = g[g.cn10_lmm_q < 0.05]
        if len(s) <= 12: return s, 'all pairs with q < 0.05'
        up = s[s.cn10_lmm_estimate > 0].nlargest(8, 'cn10_lmm_estimate'); dn = s[s.cn10_lmm_estimate < 0].nsmallest(8, 'cn10_lmm_estimate')
        return pd.concat([up, dn]).sort_values('cn10_lmm_estimate', ascending=False), (f'the {len(up)} highest' + (f' and {len(dn)} lowest' if len(dn) else '')
                                                                                        + f' of the {len(s)} pairs with q < 0.05' + ('' if len(dn) else ' (none lower)'))
    SEL = {m: select(ST[(ST.module == m) & ST.cn10_lmm_estimate.notna()]) for m in MODS}
    for m in MODS: print(m, '->', SEL[m][1], len(SEL[m][0]))
    ST['shown_in_figure'] = False
    for m in MODS: ST.loc[ST.lr_pair.isin(SEL[m][0].lr_pair), 'shown_in_figure'] = True
    COL = {'Antigen presentation': '#08519c', 'Co-stimulation / suppressive Treg': '#54278f', 'Inflammatory cytokine': '#a50f15', 'Matrix-integrin': '#006d2c'}
    plt.rcParams.update({'font.size': 14, 'font.family': 'DejaVu Sans'})
    ROWH = 0.25; TG = 0.42; XT = 0.42; XLAB = 0.32; TOP = 0.55; LEG = 0.45
    right = ['Antigen presentation', 'Inflammatory cytokine', 'Matrix-integrin']
    nL = len(SEL['Co-stimulation / suppressive Treg'][0]); nR = [len(SEL[m][0]) for m in right]
    hL = TOP + TG + nL * ROWH + XT + XLAB; hR = TOP + sum(TG + n * ROWH + XT for n in nR) + XLAB + LEG
    H = max(hL, hR); W = 14.0
    fig = plt.figure(figsize=(W, H))
    def draw(ax, m):
        g = SEL[m][0].sort_values('cn10_lmm_estimate', ascending=True).reset_index(drop=True); c = COL[m]; y = np.arange(len(g))
        ax.axvline(0, color='#888888', lw=1, zorder=0)
        sig = (g.cn10_lmm_q < 0.05).values
        ax.hlines(y, g.cn10_lmm_ci_low, g.cn10_lmm_ci_high, color=c, lw=2.2, zorder=2)
        ax.scatter(g.cn10_lmm_estimate[sig], y[sig], s=70, color=c, edgecolor=c, linewidth=1.8, zorder=3)
        ax.scatter(g.cn10_lmm_estimate[~sig], y[~sig], s=70, color='white', edgecolor=c, linewidth=1.8, zorder=3)
        lo, hi = g.cn10_lmm_ci_low.min(), g.cn10_lmm_ci_high.max(); span = hi - lo
        for yi, (h_, q) in enumerate(zip(g.cn10_lmm_ci_high, g.cn10_lmm_q)):
            s = stars(q) if np.isfinite(q) else ''
            if s: ax.text(h_ + 0.015 * span, yi, s, ha='left', va='center', fontsize=14, color=c)
        ax.set_xlim(min(lo, 0) - 0.04 * span, max(hi, 0) + 0.13 * span)
        ax.set_yticks(y); ax.set_yticklabels(g.label, fontsize=13); ax.set_ylim(-0.6, len(g) - 0.4)
        ax.set_title(SHORTM[m], fontsize=15, loc='left', pad=5)
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(axis='x', labelsize=13)
    XLABT = 'CN10 Estimate (Within-Section z)'
    xL0, xL1, xR0, xR1 = 2.35, 6.45, 9.0, 13.75
    top = H - TOP - TG; ax = fig.add_axes([xL0 / W, (top - nL * ROWH) / H, (xL1 - xL0) / W, nL * ROWH / H])
    draw(ax, 'Co-stimulation / suppressive Treg'); ax.set_xlabel(XLABT, fontsize=14)
    top = H - TOP
    for i, (m, n) in enumerate(zip(right, nR)):
        top -= TG; ax = fig.add_axes([xR0 / W, (top - n * ROWH) / H, (xR1 - xR0) / W, n * ROWH / H]); draw(ax, m)
        if i == len(right) - 1: ax.set_xlabel(XLABT, fontsize=14)
        top -= n * ROWH + XT
    fig.suptitle('CN10 Ligand–Receptor Pairs', fontsize=16, x=0.5, y=1 - 0.1 / H, va='top')
    h = [Line2D([], [], marker='o', ls='-', color='#252525', mfc='#252525', ms=8, lw=2, label='Estimate ± 95% CI, q < 0.05'),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=8, mew=1.8, lw=2, label='q ≥ 0.05')]
    fig.legend(handles=h, loc='upper left', bbox_to_anchor=((xR0 - 0.1) / W, (top - XLAB - 0.15) / H), ncol=1, frameon=False, fontsize=13, handletextpad=0.4)
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)
    # ---------------- CSV (all pairs, all stats) + legend
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    OUT = ST.drop(columns=['module']).rename(columns={'module_label': 'module'})
    OUT = OUT[['module'] + [c for c in OUT.columns if c != 'module']]
    OUT['module_cn10_mean_within_section_z_82'] = OUT.module.map(M82.mean_within_section_z); OUT['module_cn10_lmm_q_82'] = OUT.module.map(M82.lmm_q)
    OUT['metric'] = 'LIANA bivariate local score (raw) of the pair, mean over ALL cells of the CN per section (>= 100 cells); z across CNs within section'
    OUT['model'] = 'z ~ in_CN10 + (1 | section) on per-section z; BH over the pairs of the module'
    OUT.apply(sig3).to_csv(f'{STEM}.csv', index=False)
    nsec = L.section.nunique()
    def upn(m, sign): g = ST[(ST.module == m) & (ST.cn10_lmm_q < 0.05)]; return int((np.sign(g.cn10_lmm_estimate) == sign).sum())
    cnt = '; '.join(f'{SHORTM[m]}, {upn(m, 1)} higher and {upn(m, -1)} lower of {NPAIR[m]}' for m in MODS)
    leg = (f'Individual ligand-receptor (LR) pairs of the Fig 4D modules in neighbourhood CN10 (OGN⁺RSPO3⁺ Fib / PLA2G2A⁺ Fib / Endothelial), '
       f'CD resections (CosMx 6K; {nsec} sections, 4 patients; all cells). Scores are the same cell-type-agnostic spatial LIANA bivariate local scores '
       '(Gaussian kernel, 200 px bandwidth; NicheNet v2 prior-network pairs) and the same pair-to-module assignment (CellChatDB pathway annotation) as '
       'Fig 4D, analysed per pair rather than averaged over the module: for each pair, the neighbourhood value in a section is the mean raw local score over '
       'all cells of the neighbourhood (neighbourhood x section combinations with at least 100 cells), z-standardised across the 10 neighbourhoods within the '
       'section (sections in which a pair has no variation, for example no signal, are left out for that pair). Points show the CN10 estimate of a linear '
       'mixed model on the per-section values, z ~ in CN10 + (1 | section), with 95% confidence intervals; P values are Benjamini–Hochberg-corrected over '
       'the pairs of each module; filled points and * q < 0.05, ** q < 0.01, *** q < 0.001; open points, q ≥ 0.05. Pairs are ordered by estimate; '
       'positive values mean higher in CN10 than in the average neighbourhood of the same section. Co-Stimulation / Suppressive Treg shows all '
       f'{NPAIR["Co-stimulation / suppressive Treg"]} pairs and Antigen Presentation all {NPAIR["Antigen presentation"]}; Inflammatory Cytokine and Matrix–Integrin show '
       f'{SEL["Inflammatory cytokine"][1]} and {SEL["Matrix-integrin"][1]}, respectively. Pairs significant at q < 0.05: {cnt}. Labels give ligand → receptor (all pairs are single-gene ligand and receptor). The source data hold every pair with the CN10 mean within-section z (95% CI, number of sections, '
       'sections above 0), the mixed-model estimate, confidence interval and q value, a Wilcoxon signed-rank test (BH within module), and the CN4 and CN8 '
       'mean within-section z for comparison.')
    for t in (leg,): assert 'activation score' not in t.lower() and 'quiescent' not in t.lower() and 'to tregs' not in t.lower()
    open(f'{STEM}_legend.txt', 'w').write(leg + '\n')
    print(leg); print('done', flush=True)


# ======================================================================================================================
# Section lr_costim_inhibitory_modules: Pre-specified co-stimulatory, immunoregulatory and gamma-c cytokine LR modules per CN
# Panels: Supp 9A-B inputs (per-pair cache lr_pairs_v2_per_section_cn.csv, membership lr_modules_v2_membership.csv); Supp Table 10 (membership)
# ======================================================================================================================
if _run('lr_costim_inhibitory_modules'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/84_lr_modules_costim_inhibitory.py'   # working folder of this step (relative paths resolve from here)
    """Pre-specified co-stimulatory, immunoregulatory and gamma-c cytokine ligand-receptor modules, which split the broad 'Co-Stimulation / Suppressive Treg'
module of cn_lr_modules_allcells into biologically coherent modules. The definitions below are set from biology (CellChatDB pathway_name where available,
gene names otherwise) before any CN value of these modules is computed; they are not tuned on the result.
  Co-Stimulatory (Activating)          CD80 and CD86 (receptor CD28 only); ICOS (ICOSLG -> ICOS, CD28); OX40 (TNFSF4 -> TNFRSF4); GITRL (TNFSF18 -> TNFRSF18);
                                       4-1BB (CellChatDB 'CD137': TNFSF9 -> TNFRSF9); CD70 (CD70 -> CD27); CD40 (CD40LG -> CD40 only).
  Immunoregulatory (Inhibitory /       inhibitory checkpoints: PD-L1 (CD274 -> PDCD1); PDL2 (PDCD1LG2 -> PDCD1); CD80 / CD86 -> CTLA4; NECTIN2, PVR -> TIGIT, CD96;
   Suppressive)                        GALECTIN (LGALS9 -> HAVCR2 only); BTLA / CD160 (TNFRSF14 -> BTLA, CD160); CD200 (CD200 -> CD200R1); HLA-E -> KLRC1 (NKG2A);
                                       suppressive cytokines: IL10 (IL10 -> IL10RA, IL10RB); TGFb (TGFB1/2/3 -> TGFBR1, TGFBR2).
  Treg-Supportive gamma-c Cytokines    IL2 pathway restricted to IL2, IL7, IL15 -> their CellChatDB receptors (reported separately).
  Pre-specified exclusions (biologically ambiguous or outside the definition): ICOSLG -> CTLA4 (inhibitory receptor, not in the checkpoint list);
  CD40LG -> integrins (CellChatDB CD40 pathway, not co-stimulation); CD80 -> CD274 (cis PD-L1 interaction); NECTIN/PVR -> CD226 (activating);
  LGALS9 -> PTPRC / CD44 / P4HB; TGFb -> TGFBR3 or ACVR receptors; IL10-family cytokines other than IL10; HLA-E -> KLRD1 (CD94 pairs with both inhibitory
  NKG2A and activating NKG2C, so it is not specific without KLRC1); reverse-direction NicheNet pairs (receptor -> ligand).
Matching: a NicheNet pair joins a module when a reference interaction has the same ligand and an overlapping receptor (the module_of rule of
  cn_lr_modules_allcells); CellChatDB-derived references keep CellChatDB's heteromeric receptors, gene-name references are single genes.
Antigen Presentation, Inflammatory Cytokine and Matrix-Integrin are taken unchanged from lr_modules_allcells_membership.csv.
Scores / aggregation / statistics as cn_lr_modules_allcells (main 'raw' variant): per-cell raw LIANA bivariate local score (NaN -> 0), module = mean over the
  module's pairs per cell; CN x section value = mean over all cells of that CN in that section (>= 100 cells); z across CNs within section; mean over
  sections; linear mixed model z ~ in_CN + (1 | section), BH over the 10 CNs within module (Wilcoxon signed-rank q in the CSV). Computed per pair x section
  x CN as in cn10_lr_pairs (module value = mean over its pairs, identical by linearity); the unchanged modules and the broad module recomputed this way are
  asserted equal to cn_lr_modules_allcells.
Outputs: A_fig4_panels/fig4_lr_modules_allcells_v2_dotplot.{png,csv} + _legend.txt; A_fig4_panels/supp_lr_modules_allcells_v2_heatmap.{png,csv} + _legend.txt;
  B_resection_analyses/lr_modules_cellchat/: lr_modules_v2_membership.csv (pair list + sub-family), lr_modules_v2_wanted_pairs.csv (pre-specified vs present
  in the LIANA output), lr_pairs_v2_per_section_cn.csv (per-pair cache), lr_pairs_v2_detection.csv (fraction of cells with a non-zero score per pair),
  lr_modules_v2_cn10_pairs.csv (per-pair CN10 statistics, as cn10_lr_pairs).
SKIP_COMPUTE=1 redraws from the cache."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib, h5py
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.lines import Line2D
    from anndata.io import read_elem
    from scipy import stats as sst, sparse
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; SPOT = f'{R}/combined/spot_lr_df.parquet'
    HSTEM = f'{OUTA}/supp_lr_modules_allcells_v2_heatmap'; SSTEM = f'{OUTA}/fig4_lr_modules_allcells_v2_dotplot'
    CACHE = f'{OUTB}/lr_pairs_v2_per_section_cn.csv'; CACHE83 = f'{OUTB}/lr_pairs_allcells_per_section_cn.csv'
    MIN_CELLS_CN = 100; CNS = list(range(1, 11)); TARGET = 10
    CS, IR, GC = 'Co-Stimulatory (Activating)', 'Immunoregulatory (Inhibitory / Suppressive)', 'Treg-Supportive γc Cytokines'
    AP, IC, MI = 'Antigen presentation', 'Inflammatory cytokine', 'Matrix-integrin'
    MODS = [AP, CS, IR, IC, MI, GC]          # gamma-c module last, shown separated
    BROAD = 'Co-stimulation / suppressive Treg'
    # ---------------- PRE-SPECIFIED module definitions (set before any CN value of these modules is computed)
    cc = pd.read_csv(f'{D}/objects/cellchatdb_human_interactions_snapshot.csv')
    gs = lambda x: frozenset(y.strip() for y in str(x).replace(';', ',').split(',') if y.strip() and y.strip() != 'nan')
    cc['L'] = cc['ligand.symbol'].map(gs); cc['R'] = cc['receptor.symbol'].map(gs)
    F = frozenset
    def ccref(pw, lig=None, rec_any=None, rec_eq=None):
        s = cc[cc.pathway_name == pw]; assert len(s), pw
        out = []
        for a, b in zip(s.L, s.R):
            if lig is not None and not (a <= F(lig)): continue
            if rec_any is not None and not (b & F(rec_any)) or rec_any is not None and not (b <= F(rec_any)): continue
            if rec_eq is not None and b != F(rec_eq): continue
            out.append((a, b))
        assert out, (pw, lig, rec_any, rec_eq); return out
    gn = lambda ligs, recs: [(F([l]), F([r])) for l in ligs for r in recs]
    SPEC = {CS: [('CD80 -> CD28', 'CellChatDB CD80 (receptor CD28 only)', ccref('CD80', rec_eq=['CD28'])),
                 ('CD86 -> CD28', 'CellChatDB CD86 (receptor CD28 only)', ccref('CD86', rec_eq=['CD28'])),
                 ('ICOS', 'CellChatDB ICOS (receptors ICOS, CD28)', ccref('ICOS', rec_any=['ICOS', 'CD28'])),
                 ('OX40', 'CellChatDB OX40', ccref('OX40')), ('GITRL', 'CellChatDB GITRL', ccref('GITRL')),
                 ('4-1BB', 'CellChatDB CD137', ccref('CD137', lig=['TNFSF9'], rec_eq=['TNFRSF9'])),
                 ('CD70', 'CellChatDB CD70', ccref('CD70')), ('CD40', 'CellChatDB CD40 (CD40LG -> CD40 only)', ccref('CD40', rec_eq=['CD40']))],
            IR: [('PD-L1', 'CellChatDB PD-L1', ccref('PD-L1')), ('PDL2', 'CellChatDB PDL2', ccref('PDL2')),
                 ('CD80 / CD86 -> CTLA4', 'CellChatDB CD80, CD86 (receptor CTLA4 only)', ccref('CD80', rec_eq=['CTLA4']) + ccref('CD86', rec_eq=['CTLA4'])),
                 ('NECTIN / PVR -> TIGIT, CD96', 'gene names (CellChatDB lists CD96 as ligand)', gn(['NECTIN2', 'PVR'], ['TIGIT', 'CD96'])),
                 ('GALECTIN-9 -> TIM-3', 'CellChatDB GALECTIN (receptor HAVCR2 only)', ccref('GALECTIN', rec_eq=['HAVCR2'])),
                 ('HVEM -> BTLA, CD160', 'gene names (CellChatDB lists BTLA / CD160 as ligand)', gn(['TNFRSF14'], ['BTLA', 'CD160'])),
                 ('CD200 -> CD200R1', 'CellChatDB CD200', ccref('CD200')),
                 ('HLA-E -> NKG2A', 'gene names', gn(['HLA-E'], ['KLRC1'])),
                 ('IL10', 'CellChatDB IL10 (ligand IL10 only)', ccref('IL10', lig=['IL10'])),
                 ('TGFb', 'CellChatDB TGFb (receptors TGFBR1, TGFBR2)', ccref('TGFb', rec_eq=['TGFBR1', 'TGFBR2']))],
            GC: [('IL2', 'CellChatDB IL2 (ligand IL2)', ccref('IL2', lig=['IL2'])), ('IL7', 'CellChatDB IL2 (ligand IL7)', ccref('IL2', lig=['IL7'])),
                 ('IL15', 'CellChatDB IL2 (ligand IL15)', ccref('IL2', lig=['IL15']))]}
    cols = [c for c in pq.read_schema(SPOT).names if '^' in c]
    def match(c, refs):
        l, r = c.split('^'); L, Rr = F(l.split('_')), F(r.split('_'))
        return any((L == a) and (Rr == b or Rr <= b or b <= Rr) for a, b in refs)
    rows = []
    for m, fams in SPEC.items():
        for fam, src_, refs in fams:
            for c in cols:
                if match(c, refs): rows.append(dict(module=m, subfamily=fam, source=src_, lr_pair=c))
    NEW = pd.DataFrame(rows); assert not NEW.lr_pair.duplicated().any(), NEW[NEW.lr_pair.duplicated(keep=False)]
    m82 = pd.read_csv(f'{OUTB}/lr_modules_allcells_membership.csv'); assert m82.module.value_counts().reindex([AP, BROAD, IC, MI]).tolist() == [9, 50, 34, 227]
    BASE = m82[m82.module.isin([AP, IC, MI])].assign(subfamily='', source='cn_lr_modules_allcells')
    MEMB = pd.concat([BASE, NEW], ignore_index=True)[['module', 'subfamily', 'source', 'lr_pair']]
    MEMB['in_costim_suppressive_treg_module'] = MEMB.lr_pair.isin(m82[m82.module == BROAD].lr_pair)
    MEMB.to_csv(f'{OUTB}/lr_modules_v2_membership.csv', index=False)
    NPAIR = MEMB.module.value_counts().reindex(MODS); print(NPAIR.to_string(), flush=True)
    # wanted (expanded to ligand x receptor subunit) vs present in the LIANA output, and why missing
    with h5py.File(f'{R}/Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437/Processed/norm_log1p.h5ad', 'r') as h:
        PANEL = set(read_elem(h['var']).index.astype(str))
    colset = set(cols); WANT = []
    for m, fams in SPEC.items():
        for fam, src_, refs in fams:
            for a, b in refs:
                for l in sorted(a):
                    for r in sorted(b):
                        p = f'{l}^{r}'; pres = p in colset
                        why = '' if pres else ('; '.join([f'{g} not on 6K panel' for g in (l, r) if g not in PANEL]) or 'not in NicheNet / LIANA output')
                        WANT.append(dict(module=m, subfamily=fam, wanted_pair=p, present_in_liana=pres, reason_missing=why))
    WANT = pd.DataFrame(WANT).drop_duplicates()
    # pairs of the broad co-stimulation module that are in none of the pre-specified modules
    dropped = sorted(set(m82[m82.module == BROAD].lr_pair) - set(NEW.lr_pair))
    # ---------------- per pair x section x CN means over ALL cells (as cn10_lr_pairs) + per-pair detection
    NP = NEW.lr_pair.tolist()
    if os.environ.get('SKIP_COMPUTE') != '1':
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas(); assert not cn.index.duplicated().any()
        idx_all = pq.read_table(SPOT, columns=['__index_level_0__']).column(0).to_pylist(); rows_ = pd.Index(idx_all); assert not rows_.duplicated().any()
        keep = rows_.isin(cn.index); print('LIANA cells:', len(rows_), '| with a CN:', int(keep.sum()), flush=True)
        sec = cn.section.reindex(rows_[keep]).values; cnn = cn.CN_k10.reindex(rows_[keep]).values.astype(int)
        codes, uniq = pd.factorize(pd.MultiIndex.from_arrays([sec, cnn], names=['section', 'CN']))
        M = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
        ncell = np.asarray(M.sum(1)).ravel(); kidx = np.where(keep)[0]; out = []; det = []
        for i in range(0, len(NP), 64):
            ch = NP[i:i + 64]; A = pq.read_table(SPOT, columns=ch).to_pandas().astype(np.float64).fillna(0.0).values[kidx]
            out.append(pd.DataFrame(np.asarray(M @ A) / ncell[:, None], columns=ch))
            det.append(pd.DataFrame({'lr_pair': ch, 'frac_cells_nonzero': (A > 0).mean(0), 'mean_raw_all_cells': A.mean(0)})); del A; gc.collect()
        W = pd.concat(out, axis=1); W.index = pd.MultiIndex.from_tuples(list(uniq), names=['section', 'CN']); W['n_cells'] = ncell
        L = W.reset_index().melt(id_vars=['section', 'CN', 'n_cells'], var_name='lr_pair', value_name='raw')
        L.to_csv(CACHE, index=False, float_format='%.8g'); pd.concat(det).to_csv(f'{OUTB}/lr_pairs_v2_detection.csv', index=False, float_format='%.4g')
    L = pd.read_csv(CACHE); DET = pd.read_csv(f'{OUTB}/lr_pairs_v2_detection.csv').set_index('lr_pair')
    L83 = pd.read_csv(CACHE83)
    LL = pd.concat([L83[L83.lr_pair.isin(BASE.lr_pair)][['section', 'CN', 'n_cells', 'lr_pair', 'raw']], L[['section', 'CN', 'n_cells', 'lr_pair', 'raw']]], ignore_index=True)
    LL = LL.merge(MEMB[['lr_pair', 'module']], on='lr_pair')
    assert set(LL.lr_pair) == set(MEMB.lr_pair)
    V = LL.groupby(['section', 'CN', 'module']).agg(raw=('raw', 'mean'), n_cells=('n_cells', 'first')).reset_index()
    V82 = pd.read_csv(f'{OUTB}/lr_modules_allcells_per_section_cn.csv')
    chk = V[V.module.isin([AP, IC, MI])].merge(V82, on=['section', 'CN', 'module'], suffixes=('', '_82'))
    assert len(chk) == (V82.module.isin([AP, IC, MI])).sum() and (chk.n_cells == chk.n_cells_82).all() and np.allclose(chk.raw, chk.raw_82, rtol=1e-4, atol=1e-7)
    nchk = V[V.module == CS][['section', 'CN', 'n_cells']].merge(V82[V82.module == BROAD][['section', 'CN', 'n_cells']], on=['section', 'CN'], suffixes=('', '_82'))
    assert len(nchk) == (V82.module == BROAD).sum() and (nchk.n_cells == nchk.n_cells_82).all()
    print('unchanged modules reproduce cn_lr_modules_allcells: max abs diff', float((chk.raw - chk.raw_82).abs().max()), flush=True)
    BROADV = V82[V82.module == BROAD][['section', 'CN', 'n_cells', 'raw']].assign(module=BROAD)
    V = pd.concat([V, BROADV], ignore_index=True)
    V = V[V.n_cells >= MIN_CELLS_CN].copy()
    zs = lambda x: (x - x.mean()) / x.std()
    V['z_raw'] = V.groupby(['section', 'module']).raw.transform(zs); V.loc[~np.isfinite(V.z_raw), 'z_raw'] = np.nan
    def fit(dd):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm('y ~ in_cn', dd, groups=dd.section).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], meth
            except Exception: pass
        return np.nan, np.nan, np.nan, 'failed'
    # ---------------- CN names / labels (as cn_lr_modules_allcells)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [(' Fib / PLA2G2A⁺ Fib / Endothelial', '/PLA2G2A⁺ Fib'), (' / Crypt Top Fib / Cycling', ' / CryptTop / Cycl'),
            (' / ADAMDEC1 Fib', ' / ADAMDEC1'), (' / Neuron / Pericyte', ' / Neuron / Peri'), (' / Mac S_M_P / CD9 Mac', ' / Mac / CD9 Mac'),
            ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'), ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'),
            ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {c: short_cn(NM[c]) for c in CNS}
    ST = []
    for m in MODS + [BROAD]:
        d = V[V.module == m].dropna(subset=['z_raw'])
        for c in CNS:
            dd = d.assign(y=d.z_raw, in_cn=(d.CN == c).astype(float)); zc = dd[dd.CN == c].y.values
            b, se, p, meth = fit(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
            t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
            ST.append(dict(module=m, CN=c, CN_label=LAB[c], CN_name=NM[c], n_pairs=int((MEMB.module == m).sum()) if m != BROAD else 50, n_sections=len(zc),
                           mean_within_section_z=zc.mean(), mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)),
                           mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)), sections_above_0=int((zc > 0).sum()),
                           wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan,
                           lmm_estimate=b, lmm_ci_low=b - 1.96 * se, lmm_ci_high=b + 1.96 * se, lmm_p=p, optimiser=meth))
    ST = pd.DataFrame(ST); ST['lmm_q'] = np.nan; ST['wilcoxon_q'] = np.nan
    for m, g in ST.groupby('module'):
        for k_ in ['lmm', 'wilcoxon']:
            ok = g[f'{k_}_p'].notna(); ST.loc[g.index[ok], f'{k_}_q'] = multipletests(g[f'{k_}_p'][ok], method='fdr_bh')[1]
    # unchanged modules and the broad module must equal the cn_lr_modules_allcells heatmap values
    H82 = pd.read_csv(f'{OUTA}/supp_lr_modules_allcells_heatmap.csv')
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    SHORTM = {AP: 'Antigen Presentation', CS: 'Co-Stimulatory (Activating)', IR: 'Immunoregulatory (Inhibitory / Suppressive)', IC: 'Inflammatory Cytokine',
              MI: 'Matrix–Integrin', GC: 'Treg-Supportive γc Cytokines', BROAD: 'Co-Stimulation / Suppressive Treg'}
    for m in [AP, IC, MI, BROAD]:
        a = ST[ST.module == m].set_index('CN'); b = H82[H82.module == SHORTM[m]].set_index('CN')
        assert np.allclose(a.loc[CNS, 'mean_within_section_z'], b.loc[CNS, 'mean_within_section_z'], rtol=5e-3, atol=1e-3), m   # the heatmap CSV is rounded to 3 sig figs
        assert (a.loc[CNS, 'n_sections'].values == b.loc[CNS, 'n_sections'].values).all(), m
        # q: same significance category as cn_lr_modules_allcells (mixed-model optimiser convergence differs slightly on inputs re-read from the CSV caches)
        cat = lambda q: np.digitize(q, [0.001, 0.01, 0.05]); assert (cat(a.loc[CNS, 'lmm_q'].values) == cat(b.loc[CNS, 'lmm_q'].values)).all(), m
        print(m, 'max |log10 q - log10 q82|', float(np.abs(np.log10(a.loc[CNS, 'lmm_q'].values) - np.log10(b.loc[CNS, 'lmm_q'].values)).max()), flush=True)
    print('unchanged modules (incl. the broad co-stimulation module) reproduce the cn_lr_modules_allcells heatmap values', flush=True)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
    print(ST.pivot(index='CN', columns='module', values='mean_within_section_z').round(2).astype(str).add(
          ST.pivot(index='CN', columns='module', values='lmm_q').applymap(stars)).loc[CNS, MODS + [BROAD]].to_string())
    # ---------------- per-pair CN10 statistics for the pre-specified modules (as cn10_lr_pairs; reported, not plotted)
    LP = L.merge(MEMB[['lr_pair', 'module', 'subfamily']], on='lr_pair'); LP = LP[LP.n_cells >= MIN_CELLS_CN].copy()
    LP['z'] = LP.groupby(['section', 'lr_pair']).raw.transform(zs); LP.loc[~np.isfinite(LP.z), 'z'] = np.nan
    PR = []
    for p, d in LP.groupby('lr_pair'):
        d2 = d.dropna(subset=['z']); dd = d2.assign(y=d2.z, in_cn=(d2.CN == TARGET).astype(float)); zc = dd[dd.CN == TARGET].y.values
        b, se, pv, meth = fit(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
        PR.append(dict(module=SHORTM[d.module.iloc[0]], subfamily=d.subfamily.iloc[0], lr_pair=p, label=p.replace('^', ' → ').replace('_', '+'),
                       frac_cells_nonzero=DET.loc[p, 'frac_cells_nonzero'], mean_raw_all_cells=DET.loc[p, 'mean_raw_all_cells'],
                       cn10_n_sections=len(zc), cn10_mean_within_section_z=zc.mean() if len(zc) else np.nan, cn10_sections_above_0=int((zc > 0).sum()),
                       cn10_lmm_estimate=b, cn10_lmm_ci_low=b - 1.96 * se, cn10_lmm_ci_high=b + 1.96 * se, cn10_lmm_p=pv,
                       cn8_mean_within_section_z=d2[d2.CN == 8].z.mean(), cn4_mean_within_section_z=d2[d2.CN == 4].z.mean()))
    PR = pd.DataFrame(PR); PR['cn10_lmm_q'] = np.nan
    for m, g in PR.groupby('module'):
        ok = g.cn10_lmm_p.notna(); PR.loc[g.index[ok], 'cn10_lmm_q'] = multipletests(g.cn10_lmm_p[ok], method='fdr_bh')[1]
    PR = PR.sort_values(['module', 'cn10_lmm_estimate'], ascending=[True, False])
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    PR.apply(sig3).to_csv(f'{OUTB}/lr_modules_v2_cn10_pairs.csv', index=False)
    print(PR[['module', 'label', 'frac_cells_nonzero', 'cn10_n_sections', 'cn10_lmm_estimate', 'cn10_lmm_q', 'cn8_mean_within_section_z']].round(4).to_string(index=False))
    WANT = WANT.merge(DET[['frac_cells_nonzero', 'mean_raw_all_cells']], left_on='wanted_pair', right_index=True, how='left')
    WANT['module'] = WANT.module.map(SHORTM); WANT.apply(sig3).to_csv(f'{OUTB}/lr_modules_v2_wanted_pairs.csv', index=False)
    print(WANT.to_string(index=False)); print('broad co-stimulation pairs not in a pre-specified module:', dropped)
    # ---------------- heatmap (gamma-c column separated)
    S_ = ST[ST.module.isin(MODS)].copy()
    E = S_.pivot(index='CN', columns='module', values='mean_within_section_z').loc[CNS, MODS]; Q = S_.pivot(index='CN', columns='module', values='lmm_q').loc[CNS, MODS]
    assert E.notna().all().all() and Q.notna().all().all()
    VM = float(np.ceil(np.abs(E.values).max() * 10) / 10)
    ann = np.array([[f'{E.iloc[r, c]:+.2f}{stars(Q.iloc[r, c])}' for c in range(len(MODS))] for r in range(len(CNS))])
    SHORT = ['Antigen\nPresent.', 'Co-Stim.\n(Activ.)', 'Immuno-\nregulatory', 'Inflamm.\nCytokine', 'Matrix-\nIntegrin', 'γc\nCytokine']
    plt.rcParams.update({'font.size': 17, 'font.family': 'DejaVu Sans'})
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    sns.heatmap(E.set_axis([LAB[c] for c in CNS]), cmap='RdBu_r', center=0, vmin=-VM, vmax=VM, annot=ann, fmt='', annot_kws={'fontsize': 14}, ax=ax,
                linewidths=0.5, linecolor='white', cbar_kws={'label': 'Mean Within-Section z', 'shrink': 0.85, 'pad': 0.012})
    ax.axvline(len(MODS) - 1, color='white', lw=6)
    ax.set_xticklabels(SHORT, rotation=0, fontsize=14); ax.set_yticklabels([LAB[c] for c in CNS], rotation=0, fontsize=13)
    ax.set_xlabel(''); ax.set_ylabel('Spatial Neighborhood', fontsize=15)
    fig.savefig(f'{HSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    H = ST.assign(module_label=ST.module.map(SHORTM), cell_text=[ann[CNS.index(c), MODS.index(m)] if m in MODS else '' for c, m in zip(ST.CN, ST.module)],
                  shown_in_figure=ST.module.isin(MODS),
                  metric='LIANA bivariate local score (raw), mean over module pairs per cell, mean over ALL cells of the CN per section; z across CNs within section, mean over sections',
                  model='z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within module')
    H = H.drop(columns=['module']).rename(columns={'module_label': 'module'}); H = H[['module'] + [c for c in H.columns if c != 'module']]
    H.apply(sig3).to_csv(f'{HSTEM}.csv', index=False)
    # ---------------- per-section dot plot (long labels, no bold)
    COL = {AP: ('#9ecae1', '#08519c'), CS: ('#fdae6b', '#d94801'), IR: ('#bcbddc', '#54278f'), IC: ('#fcae91', '#a50f15'), MI: ('#a1d99b', '#006d2c'),
           GC: ('#dfc27d', '#8c510a')}
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(len(MODS), 1, figsize=(7.6, 17.6), sharex=True, gridspec_kw={'hspace': 0.27})
    rng2 = np.random.default_rng(1)
    for ax, m in zip(axes, MODS):
        lc, dc = COL[m]; g = S_[S_.module == m].set_index('CN'); vv = V[V.module == m]
        ax.axhline(0, color='#888888', lw=1, zorder=0)
        for i, c in enumerate(CNS):
            z = vv[vv.CN == c].z_raw.dropna().values
            ax.scatter(i + rng2.uniform(-0.14, 0.14, len(z)), z, s=45, color=lc, edgecolor='white', linewidth=0.5, alpha=0.95, zorder=2)
            e, lo, hi, q = g.loc[c, ['mean_within_section_z', 'mean_ci_low', 'mean_ci_high', 'lmm_q']]
            ax.errorbar(i, e, yerr=[[e - lo], [hi - e]], fmt='o', ms=10, color=dc, mfc=dc if q < 0.05 else 'white', mew=2.2, elinewidth=2.6, capsize=0, zorder=3)
        vals = np.concatenate([vv.z_raw.dropna().values, g.mean_ci_low.values, g.mean_ci_high.values]); ymin, ymax = vals.min(), vals.max(); span = ymax - ymin; stop = [ymax]
        for i, c in enumerate(CNS):
            q = g.loc[c, 'lmm_q']
            if q < 0.05:
                top = max(vv[vv.CN == c].z_raw.max(), g.loc[c, 'mean_ci_high'])
                ax.text(i, top + 0.012 * span, stars(q), ha='center', va='bottom', fontsize=17, color=dc); stop.append(top + 0.2 * span)
        ax.set_ylim(ymin - 0.05 * span, max(ymax + 0.05 * span, max(stop)))
        ax.set_title(SHORTM[m], fontsize=17, loc='left', pad=3); ax.tick_params(axis='y', labelsize=15)
        ax.spines[['top', 'right']].set_visible(False); ax.set_xlim(-0.55, len(CNS) - 0.45)
    axes[-1].set_xticks(range(len(CNS))); axes[-1].set_xticklabels([LAB[c] for c in CNS], rotation=45, ha='right', rotation_mode='anchor', fontsize=13)
    axes[-1].set_xlabel('Spatial Neighborhood', fontsize=16)
    fig.supylabel('LR Module Score (Within-Section z)', fontsize=16, x=0.04)
    h = [Line2D([], [], marker='o', ls='', color='#bdbdbd', ms=9, label='One Section'),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='#252525', ms=10, lw=2.4, label='Mean ± 95% CI, q < 0.05'),
         Line2D([], [], ls='', label=' '),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=10, mew=2, lw=2.4, label='q ≥ 0.05')]
    axes[0].legend(handles=h, loc='lower left', bbox_to_anchor=(-0.02, 1.12), ncol=2, frameon=False, fontsize=14, handletextpad=0.3, columnspacing=0.8, borderaxespad=0)
    fig.savefig(f'{SSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    SUP = pd.concat([H.assign(row_type='summary'), V.assign(module=V.module.map(SHORTM)).rename(columns={'z_raw': 'within_section_z'}).assign(row_type='section')],
                    ignore_index=True)
    SUP.apply(sig3).to_csv(f'{SSTEM}.csv', index=False)
    # ---------------- legends
    def plist(m): return ', '.join(f"{p.split('^')[0]} → {p.split('^')[1].replace('_', '+')}" for p in MEMB[MEMB.module == m].lr_pair)
    miss = WANT[~WANT.present_in_liana]
    missing_txt = '; '.join(f"{r.wanted_pair.replace('^', ' → ')} ({r.reason_missing})" for r in miss.itertuples())
    ncell_txt = '2,936,007'
    common = (f'CD resections (CosMx 6K; 11 sections, 4 patients; all {ncell_txt} cells with a neighbourhood). Ligand-receptor (LR) scores are cell-type-agnostic '
          'spatial LIANA bivariate local scores (Gaussian kernel, 200 px bandwidth) for the NicheNet v2 prior-network pairs (as Fig 3), computed for every cell. '
          'Modules were defined from biology before analysis (CellChatDB pathway annotation where available, gene names otherwise; a pair joins a module '
          'when a reference interaction has the same ligand and an overlapping receptor): Antigen Presentation = CellChatDB MHC-II '
          f'({NPAIR[AP]} pairs); Co-Stimulatory (Activating) = CD80 and CD86 with receptor CD28, ICOS, OX40, GITRL, 4-1BB, CD70 and CD40 ({NPAIR[CS]} pairs: '
          f'{plist(CS)}); Immunoregulatory (Inhibitory / Suppressive) = the inhibitory checkpoints PD-L1, PD-L2, CD80 / CD86 → CTLA4, NECTIN2 / PVR → TIGIT '
          f'and CD96, galectin-9 → TIM-3, HVEM → BTLA and CD160, CD200 → CD200R1 and HLA-E → NKG2A plus the suppressive cytokines IL-10 and TGF-β ({NPAIR[IR]} '
          f'pairs: {plist(IR)}); Inflammatory Cytokine = IL1, IL6, TNF, IFN-II, IL12, IL17, OSM ({NPAIR[IC]} pairs); Matrix–Integrin = ECM-Receptor pathways '
          f'({NPAIR[MI]} pairs); Treg-Supportive γc Cytokines = IL-2, IL-7 and IL-15 with their receptors ({NPAIR[GC]} pairs: {plist(GC)}; shown separately). '
          f'Pre-specified interactions without a score: {missing_txt}. No sender or receiver cell types are imposed. A cell\'s module score is the mean raw '
          'local score over the module\'s pairs, and a neighbourhood\'s value in a section is the mean over all its cells, so the score reflects both '
          'ligand–receptor expression and the abundance of interacting cells (neighbourhood × section combinations with at least 100 cells). Values are '
          'z-standardised across neighbourhoods within each section and averaged over sections ("Mean Within-Section z"; 0 = neighbourhood average). '
          'Significance: linear mixed model on the per-section values, z ~ in neighbourhood + (1 | section), BH-corrected over the 10 neighbourhoods within '
          'each module; * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. ')
    legH = ('Spatial ligand-receptor module activity across cellular neighbourhoods (rows) by module (columns), ' + common
            + f'Colour scale symmetric, ±{VM:.1f}; the γc cytokine column is set apart. Section-level values are shown in the companion dot plot. The source data '
        'hold the per-section values, mixed-model estimates and q values, a Wilcoxon signed-rank test (BH within module) and, for reference, the '
        'combined Co-Stimulation / Suppressive Treg module. Neighbourhood names are abbreviated by their dominant cell types; full names are in the source data.')
    legS = ('Section-level spatial ligand-receptor module activity across cellular neighbourhoods, ' + common
            + 'Light dots: one section each. Dark points and bars: mean over sections ± 95% CI (t distribution); filled, q < 0.05. Neighbourhood names are '
        'abbreviated by their dominant cell types; full names are in the source data.')
    for t in (legH, legS): assert 'activation score' not in t.lower() and 'quiescent' not in t.lower() and 'territory' not in t.lower() and 'selenop+ fib' not in t.lower()
    open(f'{HSTEM}_legend.txt', 'w').write(legH + '\n'); open(f'{SSTEM}_legend.txt', 'w').write(legS + '\n')
    print('colour limits +/-', VM); print('done', flush=True)


# ======================================================================================================================
# Section cn_lr_modules_immunomodulatory: LR module heatmap and dot plot with the merged 28-pair Immunomodulatory module per CN
# Panels: Supp 9A input (cell text of supp_lr_modules_allcells_clean_heatmap.csv)
# ======================================================================================================================
if _run('cn_lr_modules_immunomodulatory'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/82c_lr_modules_allcells_clean.py'   # working folder of this step (relative paths resolve from here)
    """All-cells ligand-receptor module heatmap and dot plot in which the broad 'Co-Stimulation / Suppressive Treg' module of cn_lr_modules_allcells is replaced
by a pre-specified 'Immunomodulatory' module.
The broad module (CellChatDB pathways CD80, CD86, PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2, IL10; 50 NicheNet pairs) also contains the IL-20-family
  cytokines (IL19, IL20, IL22, IL24, IL26 -> IL20RA / IL20RB / IL22RA1 / IL10RB), IL4, IL9, IL21, TSLP and the gamma-c cytokines IL2 / IL7 / IL15 (all
  through the CellChatDB IL2 / IL10 pathways), CD40LG -> integrins (ITGA2B, ITGA5, ITGAM, ITGB1/2/3), CD80 -> CD274 (cis) and ICOSLG -> CTLA4.
Immunomodulatory (28 pairs) = the co-stimulatory and immunoregulatory modules of lr_costim_inhibitory_modules merged (membership read from
  lr_modules_v2_membership.csv, asserted):
  Co-Stimulatory (Activating), 9: CD80 -> CD28, CD86 -> CD28, ICOSLG -> CD28, ICOSLG -> ICOS, TNFSF4 -> TNFRSF4, TNFSF18 -> TNFRSF18, TNFSF9 -> TNFRSF9,
    CD70 -> CD27, CD40LG -> CD40.
  Immunoregulatory (Inhibitory / Suppressive), 19: CD274 -> PDCD1, PDCD1LG2 -> PDCD1, CD80 -> CTLA4, CD86 -> CTLA4, NECTIN2 / PVR -> TIGIT, CD96,
    LGALS9 -> HAVCR2, TNFRSF14 -> BTLA, CD160, IL10 -> IL10RA, IL10RB, TGFB1/2/3 -> TGFBR1, TGFBR2.
  The gamma-c cytokines (IL2 / IL7 / IL15) are not included. Antigen Presentation, Inflammatory Cytokine and Matrix-Integrin as in cn_lr_modules_allcells.
Metric / statistics as cn_lr_modules_allcells (main 'raw' variant, incl. its float32 precision path): per-cell raw LIANA bivariate local score
  (spot_lr_df.parquet, NicheNet pairs, NaN -> 0); cell module score = mean over the module's pairs; CN x section value = mean over all cells of the CN in
  the section (>= 100 cells); z across CNs within section; mean over sections ('Mean Within-Section z'); linear mixed model z ~ in_CN + (1 | section), BH
  over the 10 CNs within module (Wilcoxon signed-rank q in the CSV).
Computation:
  Antigen Presentation, Inflammatory Cytokine, Matrix-Integrin (and the broad module, CSV reference) = the cn_lr_modules_allcells computation on its
  per-cell cache.
  Immunomodulatory (and its two sub-modules, CSV only): every cell has every pair (NaN -> 0) and the cell set is the same for all pairs, so the CN x section
  module value (mean over cells of the mean over pairs) = mean over pairs of the per-pair CN x section means (lr_pairs_v2_per_section_cn.csv). Asserted
  before use: this shortcut reproduces the CN x section values of all four cn_lr_modules_allcells modules (per-pair cache of cn10_lr_pairs, max abs diff
  < 1e-6), the two per-pair caches agree for the pairs they share, and the unchanged modules and the broad module reproduce the cn_lr_modules_allcells
  heatmap values.
Outputs: A_fig4_panels/fig4_lr_modules_allcells_clean_dotplot.{png,csv} + _legend.txt; A_fig4_panels/supp_lr_modules_allcells_clean_heatmap.{png,csv}
  + _legend.txt (cell text used by S9A_cn_characterisation); A_fig4_panels/supp_lr_modules_allcells_clean_submodules.csv (Co-Stimulatory and
  Immunoregulatory separately, CN10 / CN4 / CN8)."""
    import os, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.lines import Line2D
    from scipy import stats as sst
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    warnings.simplefilter('ignore')
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    HSTEM = f'{OUTA}/supp_lr_modules_allcells_clean_heatmap'; SSTEM = f'{OUTA}/fig4_lr_modules_allcells_clean_dotplot'
    SUBCSV = f'{OUTA}/supp_lr_modules_allcells_clean_submodules.csv'
    CACHE82 = f'{OUTB}/lr_modules_allcells_per_cell.parquet'; PAIR83 = f'{OUTB}/lr_pairs_allcells_per_section_cn.csv'; PAIR84 = f'{OUTB}/lr_pairs_v2_per_section_cn.csv'
    MIN_CELLS_CN = 100; CNS = list(range(1, 11)); REPORT_CNS = [10, 4, 8]
    AP, BROAD, IC, MI = 'Antigen presentation', 'Co-stimulation / suppressive Treg', 'Inflammatory cytokine', 'Matrix-integrin'
    CS, IR, IMM = 'Co-Stimulatory (Activating)', 'Immunoregulatory (Inhibitory / Suppressive)', 'Immunomodulatory'
    MODS = [AP, IMM, IC, MI]          # shown, in the cn_lr_modules_allcells order (Immunomodulatory in place of the broad module)
    REF = [BROAD, CS, IR]             # CSV only
    MOD4 = [AP, BROAD, IC, MI]
    # ---------------- membership: cn_lr_modules_allcells (unchanged modules) and lr_costim_inhibitory_modules (pre-specified sub-modules), asserted against the definition
    m82 = pd.read_csv(f'{OUTB}/lr_modules_allcells_membership.csv'); assert not m82.lr_pair.duplicated().any()
    assert m82.module.value_counts().reindex(MOD4).tolist() == [9, 50, 34, 227]
    vmem = pd.read_csv(f'{OUTB}/lr_modules_v2_membership.csv')
    CS_P = vmem[vmem.module == CS].lr_pair.tolist(); IR_P = vmem[vmem.module == IR].lr_pair.tolist()
    CS_EXP = ['CD80^CD28', 'CD86^CD28', 'ICOSLG^CD28', 'ICOSLG^ICOS', 'TNFSF4^TNFRSF4', 'TNFSF18^TNFRSF18', 'TNFSF9^TNFRSF9', 'CD70^CD27', 'CD40LG^CD40']
    CHK_EXP = ['CD274^PDCD1', 'PDCD1LG2^PDCD1', 'CD80^CTLA4', 'CD86^CTLA4', 'NECTIN2^TIGIT', 'NECTIN2^CD96', 'PVR^TIGIT', 'PVR^CD96', 'LGALS9^HAVCR2',
               'TNFRSF14^BTLA', 'TNFRSF14^CD160']
    CYT_EXP = ['IL10^IL10RA', 'IL10^IL10RB'] + [f'TGFB{i}^TGFBR{j}' for i in (1, 2, 3) for j in (1, 2)]
    assert sorted(CS_P) == sorted(CS_EXP) and len(CS_P) == 9
    assert sorted(IR_P) == sorted(CHK_EXP + CYT_EXP) and len(IR_P) == 19
    IMM_P = CS_P + IR_P; assert len(set(IMM_P)) == 28
    assert not any(p.split('^')[0] in ('IL2', 'IL7', 'IL15') for p in IMM_P)          # gamma-c cytokines not included
    PAIRS = {IMM: IMM_P, CS: CS_P, IR: IR_P}
    NPAIR = {AP: 9, BROAD: 50, IC: 34, MI: 227, IMM: 28, CS: 9, IR: 19}
    dropped = sorted(set(m82[m82.module == BROAD].lr_pair) - set(IMM_P)); added = sorted(set(IMM_P) - set(m82[m82.module == BROAD].lr_pair))
    print('broad co-stimulation pairs not in Immunomodulatory (%d):' % len(dropped), dropped); print('pairs added (%d):' % len(added), added, flush=True)
    # ---------------- (1) the four modules of cn_lr_modules_allcells, from its per-cell cache
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas(); assert not cn.index.duplicated().any()
    C = pd.read_parquet(CACHE82, columns=[f'raw|{m}' for m in MOD4])
    X = C.join(cn[['CN_k10', 'section']], how='inner').rename(columns={'CN_k10': 'CN'})
    NCELL = len(X); print('cells with LIANA scores:', len(C), '| with a CN:', NCELL, flush=True)
    G = X.groupby(['section', 'CN'])
    V82 = (G[[f'raw|{m}' for m in MOD4]].mean().set_axis(MOD4, axis=1).reset_index().melt(id_vars=['section', 'CN'], var_name='module', value_name='raw')
           .merge(G.size().rename('n_cells').reset_index(), on=['section', 'CN']))
    del C, X; assert V82.raw.dtype == np.float32   # precision path of cn_lr_modules_allcells (float32 per-cell cache -> float32 CN x section means)
    f82 = pd.read_csv(f'{OUTB}/lr_modules_allcells_per_section_cn.csv').merge(V82, on=['section', 'CN', 'module'], suffixes=('_file', ''))
    assert len(f82) == len(V82) and (f82.n_cells == f82.n_cells_file).all() and np.allclose(f82.raw, f82.raw_file, rtol=1e-5, atol=1e-9)   # file: 6 sig figs
    # ---------------- (2) shortcut check: mean over pairs of per-pair CN x section means == the module values of cn_lr_modules_allcells (all four modules)
    L83 = pd.read_csv(PAIR83); assert set(L83.lr_pair) == set(m82.lr_pair)
    for m in MOD4:
        s = L83[L83.module == m].groupby(['section', 'CN']).agg(raw=('raw', 'mean'), n=('lr_pair', 'size'), n_cells=('n_cells', 'first'))
        assert (s.n == NPAIR[m]).all()
        r = V82[V82.module == m].set_index(['section', 'CN']).reindex(s.index)
        d = float((s.raw - r.raw).abs().max()); assert d < 1e-6 and (s.n_cells == r.n_cells).all(), (m, d)
        print(f'shortcut reproduces cn_lr_modules_allcells for {m}: max abs diff {d:.2e} over {len(s)} CN x section values', flush=True)
    L = pd.read_csv(PAIR84); assert set(IMM_P) <= set(L.lr_pair)
    ov = L.merge(L83, on=['section', 'CN', 'lr_pair'], suffixes=('', '_83'))
    assert ov.lr_pair.nunique() >= 14 and (ov.n_cells == ov.n_cells_83).all() and np.allclose(ov.raw, ov.raw_83, rtol=1e-7, atol=1e-10)
    print('lr_costim_inhibitory_modules pair cache == cn10_lr_pairs pair cache for', ov.lr_pair.nunique(), 'shared pairs', flush=True)
    # ---------------- (3) Immunomodulatory and sub-modules: mean over pairs per CN x section
    VN = []
    for m, ps in PAIRS.items():
        s = L[L.lr_pair.isin(ps)].groupby(['section', 'CN']).agg(raw=('raw', 'mean'), n=('lr_pair', 'size'), n_cells=('n_cells', 'first')).reset_index()
        assert (s.n == len(ps)).all(); VN.append(s.drop(columns='n').assign(module=m))
    VN = pd.concat(VN, ignore_index=True)
    chk = VN[VN.module == IMM].merge(V82[V82.module == AP], on=['section', 'CN'], suffixes=('', '_ap')); assert len(chk) == (V82.module == AP).sum() and (chk.n_cells == chk.n_cells_ap).all()
    VN['raw'] = VN.raw.astype(np.float32)   # same float32 precision path as cn_lr_modules_allcells (the LMM optimiser fallback is sensitive to it at the variance boundary)
    V = pd.concat([V82[['section', 'CN', 'module', 'raw', 'n_cells']], VN[['section', 'CN', 'module', 'raw', 'n_cells']]], ignore_index=True)
    assert V.raw.dtype == np.float32
    V = V[V.n_cells >= MIN_CELLS_CN].copy()
    zs = lambda x: (x - x.mean()) / x.std()
    V['z_raw'] = V.groupby(['section', 'module']).raw.transform(zs); V.loc[~np.isfinite(V.z_raw), 'z_raw'] = np.nan
    print('sections without variation (z undefined) per module:', V[V.z_raw.isna()].groupby('module').section.unique().to_dict(), flush=True)
    def fit(dd):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm('y ~ in_cn', dd, groups=dd.section).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], meth
            except Exception: pass
        return np.nan, np.nan, np.nan, 'failed'
    # ---------------- CN names / labels (as cn_lr_modules_allcells)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [(' Fib / PLA2G2A⁺ Fib / Endothelial', '/PLA2G2A⁺ Fib'), (' / Crypt Top Fib / Cycling', ' / CryptTop / Cycl'),
            (' / ADAMDEC1 Fib', ' / ADAMDEC1'), (' / Neuron / Pericyte', ' / Neuron / Peri'), (' / Mac S_M_P / CD9 Mac', ' / Mac / CD9 Mac'),
            ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'), ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'),
            ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {c: short_cn(NM[c]) for c in CNS}
    ST = []
    for m in MODS + REF:
        d = V[V.module == m].dropna(subset=['z_raw'])
        for c in CNS:
            dd = d.assign(y=d.z_raw, in_cn=(d.CN == c).astype(float)); zc = dd[dd.CN == c].y.values
            b, se, p, meth = fit(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
            t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
            ST.append(dict(module=m, CN=c, CN_label=LAB[c], CN_name=NM[c], n_pairs=NPAIR[m], n_sections=len(zc), mean_within_section_z=zc.mean(),
                           mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)), mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)),
                           sections_above_0=int((zc > 0).sum()), wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan,
                           lmm_estimate=b, lmm_ci_low=b - 1.96 * se, lmm_ci_high=b + 1.96 * se, lmm_p=p, optimiser=meth))
    ST = pd.DataFrame(ST); ST['lmm_q'] = np.nan; ST['wilcoxon_q'] = np.nan
    for m, g in ST.groupby('module'):
        for k_ in ['lmm', 'wilcoxon']:
            ok = g[f'{k_}_p'].notna(); ST.loc[g.index[ok], f'{k_}_q'] = multipletests(g[f'{k_}_p'][ok], method='fdr_bh')[1]
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    cat = lambda q: np.digitize(q, [0.001, 0.01, 0.05])
    SHORTM = {AP: 'Antigen Presentation', IMM: 'Immunomodulatory', IC: 'Inflammatory Cytokine', MI: 'Matrix–Integrin',
              BROAD: 'Co-Stimulation / Suppressive Treg', CS: 'Co-Stimulatory (Activating)', IR: 'Immunoregulatory (Inhibitory / Suppressive)'}
    # ---------------- check against the cn_lr_modules_allcells heatmap (3 sig figs)
    H82 = pd.read_csv(f'{OUTA}/supp_lr_modules_allcells_heatmap.csv')
    for m in MOD4:
        a = ST[ST.module == m].set_index('CN').loc[CNS]; b = H82[H82.module == SHORTM[m]].set_index('CN').loc[CNS]
        assert np.allclose(a.mean_within_section_z, b.mean_within_section_z, rtol=5e-3, atol=1e-3), m
        assert (a.n_sections.values == b.n_sections.values).all() and (cat(a.lmm_q.values) == cat(b.lmm_q.values)).all(), m
        assert np.allclose(a.lmm_q, b.lmm_q, rtol=5e-3, atol=0), (m, float(np.abs(np.log10(a.lmm_q.values) - np.log10(b.lmm_q.values)).max()))
    print('Antigen Presentation, Inflammatory Cytokine, Matrix-Integrin and the broad module reproduce cn_lr_modules_allcells (z and q, 3 sig figs)', flush=True)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
    print(ST.pivot(index='CN', columns='module', values='mean_within_section_z').round(2).astype(str).add(
          ST.pivot(index='CN', columns='module', values='lmm_q').applymap(stars)).loc[CNS, MODS + REF].to_string(), flush=True)
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    # ---------------- heatmap
    S_ = ST[ST.module.isin(MODS)].copy()
    E = S_.pivot(index='CN', columns='module', values='mean_within_section_z').loc[CNS, MODS]; Q = S_.pivot(index='CN', columns='module', values='lmm_q').loc[CNS, MODS]
    assert E.notna().all().all() and Q.notna().all().all()
    VM = float(np.ceil(np.abs(E.values).max() * 10) / 10)
    ann = np.array([[f'{E.iloc[r, c]:+.2f}{stars(Q.iloc[r, c])}' for c in range(len(MODS))] for r in range(len(CNS))])
    SHORT = ['Antigen\nPresent.', 'Immuno-\nmodulatory', 'Inflamm.\nCytokine', 'Matrix-\nIntegrin']
    plt.rcParams.update({'font.size': 17, 'font.family': 'DejaVu Sans'})
    fig, ax = plt.subplots(figsize=(6.9, 5.0))
    sns.heatmap(E.set_axis([LAB[c] for c in CNS]), cmap='RdBu_r', center=0, vmin=-VM, vmax=VM, annot=ann, fmt='', annot_kws={'fontsize': 15}, ax=ax,
                linewidths=0.5, linecolor='white', cbar_kws={'label': 'Mean Within-Section z', 'shrink': 0.85, 'pad': 0.012})
    ax.set_xticklabels(SHORT, rotation=0, fontsize=14); ax.set_yticklabels([LAB[c] for c in CNS], rotation=0, fontsize=13)
    ax.set_xlabel(''); ax.set_ylabel('Spatial Neighborhood', fontsize=15)
    fig.savefig(f'{HSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    METRIC = 'LIANA bivariate local score (raw), mean over module pairs per cell, mean over ALL cells of the CN per section; z across CNs within section, mean over sections'
    MODEL = 'z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within module'
    H = ST[ST.module.isin(MODS + [BROAD])].copy()
    H['cell_text'] = [ann[CNS.index(c), MODS.index(m)] if m in MODS else '' for c, m in zip(H.CN, H.module)]
    H['shown_in_figure'] = H.module.isin(MODS); H['module'] = H.module.map(SHORTM).where(H.module != BROAD, 'Co-Stimulation / Suppressive Treg (broad module, reference)')
    H = H.assign(metric=METRIC, model=MODEL); H.apply(sig3).to_csv(f'{HSTEM}.csv', index=False)
    # ---------------- sub-modules for CN10 / CN4 / CN8 (CSV only)
    SB = ST[ST.module.isin([CS, IR, IMM]) & ST.CN.isin(REPORT_CNS)].copy()
    SB['module'] = pd.Categorical(SB.module.map({CS: 'Co-Stimulatory (Activating)', IR: 'Immunoregulatory (Inhibitory / Suppressive)', IMM: 'Immunomodulatory (merged)'}),
                                  ['Co-Stimulatory (Activating)', 'Immunoregulatory (Inhibitory / Suppressive)', 'Immunomodulatory (merged)'])
    SB['CN'] = pd.Categorical(SB.CN, REPORT_CNS); SB = SB.sort_values(['CN', 'module'])
    SB['stars'] = SB.lmm_q.map(stars)
    SB = SB[['CN', 'CN_label', 'module', 'n_pairs', 'n_sections', 'mean_within_section_z', 'mean_ci_low', 'mean_ci_high', 'sections_above_0', 'lmm_estimate',
             'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q', 'stars', 'wilcoxon_q']].assign(
        note='q = BH over all 10 CNs within each (sub-)module; z and model exactly as the figure; sub-modules not shown in the figure')
    SB.apply(sig3).to_csv(SUBCSV, index=False)
    # ---------------- per-section dot plot (no bold)
    COL = {AP: ('#9ecae1', '#08519c'), IMM: ('#bcbddc', '#54278f'), IC: ('#fcae91', '#a50f15'), MI: ('#a1d99b', '#006d2c')}
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(4, 1, figsize=(7.6, 12), sharex=True, gridspec_kw={'hspace': 0.24})
    # title of the Immunomodulatory panel: long form only if it fits within the axes width
    LONG = 'Immunomodulatory (Co-Stimulatory + Inhibitory)'
    tt = axes[1].set_title(LONG, fontsize=17, loc='left', pad=3); fig.canvas.draw()
    fits = tt.get_window_extent().width <= axes[1].get_window_extent().width
    TITLE = {**{m: SHORTM[m] for m in MODS}, IMM: LONG if fits else 'Immunomodulatory'}
    print('Immunomodulatory title:', TITLE[IMM], '| long form fits:', fits, flush=True)
    rng2 = np.random.default_rng(1)
    for ax, m in zip(axes, MODS):
        lc, dc = COL[m]; g = S_[S_.module == m].set_index('CN'); vv = V[V.module == m]
        ax.axhline(0, color='#888888', lw=1, zorder=0)
        for i, c in enumerate(CNS):
            z = vv[vv.CN == c].z_raw.dropna().values
            ax.scatter(i + rng2.uniform(-0.14, 0.14, len(z)), z, s=45, color=lc, edgecolor='white', linewidth=0.5, alpha=0.95, zorder=2)
            e, lo, hi, q = g.loc[c, ['mean_within_section_z', 'mean_ci_low', 'mean_ci_high', 'lmm_q']]
            ax.errorbar(i, e, yerr=[[e - lo], [hi - e]], fmt='o', ms=10, color=dc, mfc=dc if q < 0.05 else 'white', mew=2.2, elinewidth=2.6, capsize=0, zorder=3)
        vals = np.concatenate([vv.z_raw.dropna().values, g.mean_ci_low.values, g.mean_ci_high.values]); ymin, ymax = vals.min(), vals.max(); span = ymax - ymin; stop = [ymax]
        for i, c in enumerate(CNS):
            q = g.loc[c, 'lmm_q']
            if q < 0.05:
                top = max(vv[vv.CN == c].z_raw.max(), g.loc[c, 'mean_ci_high'])
                ax.text(i, top + 0.012 * span, stars(q), ha='center', va='bottom', fontsize=17, color=dc); stop.append(top + 0.2 * span)
        ax.set_ylim(ymin - 0.05 * span, max(ymax + 0.05 * span, max(stop)))
        ax.set_title(TITLE[m], fontsize=17, loc='left', pad=3); ax.tick_params(axis='y', labelsize=15)
        ax.spines[['top', 'right']].set_visible(False); ax.set_xlim(-0.55, len(CNS) - 0.45)
    axes[-1].set_xticks(range(len(CNS))); axes[-1].set_xticklabels([LAB[c] for c in CNS], rotation=45, ha='right', rotation_mode='anchor', fontsize=13)
    axes[-1].set_xlabel('Spatial Neighborhood', fontsize=16)
    fig.supylabel('LR Module Score (Within-Section z)', fontsize=16, x=0.04)
    h = [Line2D([], [], marker='o', ls='', color='#bdbdbd', ms=9, label='One Section'),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='#252525', ms=10, lw=2.4, label='Mean ± 95% CI, q < 0.05'),
         Line2D([], [], ls='', label=' '),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=10, mew=2, lw=2.4, label='q ≥ 0.05')]
    axes[0].legend(handles=h, loc='lower left', bbox_to_anchor=(-0.02, 1.12), ncol=2, frameon=False, fontsize=14, handletextpad=0.3, columnspacing=0.8, borderaxespad=0)
    fig.savefig(f'{SSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    SUM = ST.assign(shown_in_figure=ST.module.isin(MODS), module=ST.module.map(SHORTM), row_type='summary')
    SEC = V.assign(module=V.module.map(SHORTM), shown_in_figure=V.module.isin(MODS), row_type='section').rename(columns={'z_raw': 'within_section_z'})
    OUTS = pd.concat([SUM, SEC], ignore_index=True).assign(metric=METRIC, model=MODEL)
    for c_ in ['n_cells', 'n_pairs', 'n_sections', 'sections_above_0']: OUTS[c_] = OUTS[c_].round().astype('Int64')   # counts stay exact
    OUTS.apply(sig3).to_csv(f'{SSTEM}.csv', index=False)
    # ---------------- legends
    DET = pd.read_csv(f'{OUTB}/lr_pairs_v2_detection.csv').set_index('lr_pair').frac_cells_nonzero
    pct = lambda ps: (100 * DET.loc[ps].min(), 100 * DET.loc[ps].max())
    cd8086 = pct(['CD80^CD28', 'CD86^CD28', 'CD80^CTLA4', 'CD86^CTLA4']); ox40 = pct(['TNFSF4^TNFRSF4']); il2 = pct(['IL2^IL2RA', 'IL2^IL2RB', 'IL2^IL2RG'])
    low = pct(['TGFB2^TGFBR1', 'TGFB2^TGFBR2', 'TNFRSF14^CD160'])
    def grp(ps):
        out = {}
        for p in ps: l, r = p.split('^'); out.setdefault(l, []).append(r)
        return '; '.join(f'{l} → {", ".join(rs)}' for l, rs in out.items())
    common = (f'CD resections (CosMx 6K; 11 sections, 4 patients; all {NCELL:,} cells with a neighbourhood). Ligand-receptor (LR) scores are cell-type-agnostic '
          'spatial LIANA bivariate local scores (Gaussian kernel, 200 px bandwidth) for the NicheNet v2 prior-network pairs (as Fig 3), computed for every cell. '
          'Modules (a pair joins a module when a reference interaction has the same ligand and an overlapping receptor): Antigen Presentation = CellChatDB '
          f'MHC-II ({NPAIR[AP]} pairs); Immunomodulatory = {NPAIR[IMM]} pre-specified co-stimulatory and inhibitory pairs, namely co-stimulatory ({len(CS_EXP)}): '
          f'{grp(CS_EXP)}; inhibitory checkpoints ({len(CHK_EXP)}): {grp(CHK_EXP)}; suppressive cytokines ({len(CYT_EXP)}): {grp(CYT_EXP)}; Inflammatory '
          f'Cytokine = CellChatDB IL1, IL6, TNF, IFN-II, IL12, IL17, OSM ({NPAIR[IC]} pairs); Matrix–Integrin = CellChatDB ECM-Receptor pathways ({NPAIR[MI]} '
          'pairs). The γc cytokines IL-2, IL-7 and IL-15 are not included in the Immunomodulatory module; CD200 → CD200R1 and HLA-E → KLRC1 (NKG2A) were '
          'pre-specified but cannot be scored because CD200R1 and KLRC1 are not on the 6K panel. Some pairs are sparsely detected on the 6K panel: CD80 and '
          f'CD86 pairs have a non-zero local score in only {cd8086[0]:.0f}–{cd8086[1]:.0f}% of cells, OX40 (TNFSF4 → TNFRSF4) in {ox40[0]:.0f}% and IL-2 pairs '
          f'(not included) in {il2[0]:.0f}–{il2[1]:.0f}% (TGFB2 and CD160 pairs, {low[0]:.0f}–{low[1]:.0f}%). The Immunomodulatory module is used instead of the '
          'broad Co-Stimulation / Suppressive Treg module (CellChatDB CD80, CD86, PD-L1, PDL2, ICOS, OX40, GITRL, CD70, CD40, IL2 and IL10 pathways; '
          f'{NPAIR[BROAD]} pairs) because that module also included non-immunomodulatory cytokine pairs (the IL-20-family cytokines IL-19, IL-20, IL-22, '
          'IL-24 and IL-26; the γc-family cytokines IL-2, IL-4, IL-7, IL-9, IL-15 and IL-21; and TSLP) and CD40LG → integrin pairs. No sender or receiver '
          'cell types are imposed. '
          'A cell\'s module score is the mean raw local score over the module\'s pairs, and a neighbourhood\'s value in a section is the mean over all its '
          'cells, so the score reflects both ligand–receptor expression and the abundance of interacting cells (neighbourhood × section combinations with at '
          'least 100 cells). Values are z-standardised across neighbourhoods within each section and averaged over sections ("Mean Within-Section z"; '
          '0 = neighbourhood average). Significance: linear mixed model on the per-section values, z ~ in neighbourhood + (1 | section), BH-corrected over '
          'the 10 neighbourhoods within each module; * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. ')
    legH = ('Spatial ligand-receptor module activity across cellular neighbourhoods (rows) by module (columns), ' + common
            + f'Colour scale symmetric, ±{VM:.1f}. Section-level values are shown in the companion dot plot. The source data hold the mixed-model estimates '
        'and q values, a Wilcoxon signed-rank test (BH within module) and, for reference, the broad Co-Stimulation / Suppressive Treg module; the '
        'co-stimulatory and inhibitory sub-modules are reported separately for CN10, CN4 and CN8 (supp_lr_modules_allcells_clean_submodules.csv). '
        'Neighbourhood names are abbreviated by their dominant cell types; full names are in the source data.')
    legS = ('Section-level spatial ligand-receptor module activity across cellular neighbourhoods, ' + common
            + 'Light dots: one section each. Dark points and bars: mean over sections ± 95% CI (t distribution); filled, q < 0.05. Neighbourhood names are '
        'abbreviated by their dominant cell types; full names are in the source data.')
    for t in (legH, legS):
        tl = t.lower(); assert 'activation score' not in tl and 'quiescent' not in tl and 'territor' not in tl and 'to tregs' not in tl
    open(f'{HSTEM}_legend.txt', 'w').write(legH + '\n'); open(f'{SSTEM}_legend.txt', 'w').write(legS + '\n')
    print('colour limits +/-', VM); print('done', flush=True)


# ======================================================================================================================
# Section cn_lr_modules_costim: LR module dot plot with the 7-pair APC -> T-cell Co-Stimulatory module per CN
# Panels: Fig 5C, Supp 9A-B inputs (fig4_lr_modules_allcells_costim_dotplot.csv + _legend.txt)
# ======================================================================================================================
if _run('cn_lr_modules_costim'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/82d_lr_modules_allcells_costim.py'   # working folder of this step (relative paths resolve from here)
    """Ligand-receptor module dot plot (four stacked modules) with the canonical APC -> T-cell Co-Stimulatory module. Its CSV is the module source of
S9A_cn_characterisation, cn_signalling_table, F5C_score_signalling and S9B_cn_signalling.
Method as cn_lr_modules_immunomodulatory (same LIANA scores, per-CN x section means over all cells, >= 100 cells, within-section z, mixed model, BH over
the 10 CNs within each module):
  1. Co-Stimulatory = canonical APC -> T-cell co-stimulation (ligand on the antigen-presenting cell, receptor on the T cell / Treg), 7 pairs: CD80 -> CD28,
     CD86 -> CD28, ICOSLG -> ICOS, TNFSF4 -> TNFRSF4, TNFSF18 -> TNFRSF18, TNFSF9 -> TNFRSF9, CD70 -> CD27. Two pairs of the pre-specified 9-pair
     co-stimulatory set of lr_costim_inhibitory_modules are excluded: CD40LG -> CD40 (T cell -> APC direction; licenses APCs, does not co-stimulate
     T cells) and ICOSLG -> CD28 (non-canonical). Not plotted, source data only: Inhibitory / Suppressive (19 pairs), the 9-pair Co-Stimulatory set, the
     merged Immunomodulatory module and the broad Co-Stimulation / Suppressive Treg module.
  2. Panel order: Antigen Presentation, Inflammatory Cytokine, Co-Stimulatory, Matrix-Integrin. No heatmap is produced.
  3. Significance stars sit on the side of 0 matching the direction of the estimate: positive -> just above the highest element (section dot, CI bar, mean
     marker), negative -> just below the lowest element; exact glyph extents (TextPath) are used, y-limits are widened only as far as needed, and the
     placement is checked (inside axes, no overlap with dots / CI bars / mean markers / titles / legend, correct side of 0).
  4. Per-pair co-stimulatory values (per pair, z across CNs within section, mean over sections) for the 7 kept and 2 excluded pairs, with the fraction of
     cells with a non-zero local score, are added to the CSV (row_type 'pair_summary') and asserted against lr_modules_v2_cn10_pairs.csv.
Consistency asserts (before any output is written): the 9-pair Co-Stimulatory, Inhibitory / Suppressive and merged Immunomodulatory values and q equal
  supp_lr_modules_allcells_clean_submodules.csv (CN4 / CN8 / CN10) and the cn_lr_modules_immunomodulatory dot-plot CSV (all CNs, summary and section rows);
  Antigen Presentation, Inflammatory Cytokine, Matrix-Integrin and the broad module equal the cn_lr_modules_immunomodulatory dot-plot and heatmap CSVs
  (3 significant figures, identical after the CSV round trip); plus the checks against cn_lr_modules_allcells and the per-pair caches.
Outputs: A_fig4_panels/fig4_lr_modules_allcells_costim_dotplot.{png,csv} + _legend.txt"""
    import os, sys, io, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    from scipy import stats as sst
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    warnings.simplefilter('ignore')
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'
    OUTA = f'{P4}/A_fig4_panels'; OUTB = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    SSTEM = f'{OUTA}/fig4_lr_modules_allcells_costim_dotplot'
    # cn_lr_modules_immunomodulatory outputs (read only, for the consistency asserts)
    DOT82C = f'{OUTA}/fig4_lr_modules_allcells_clean_dotplot.csv'; HEAT82C = f'{OUTA}/supp_lr_modules_allcells_clean_heatmap.csv'
    SUB82C = f'{OUTA}/supp_lr_modules_allcells_clean_submodules.csv'; PR84 = f'{OUTB}/lr_modules_v2_cn10_pairs.csv'
    CACHE82 = f'{OUTB}/lr_modules_allcells_per_cell.parquet'; PAIR83 = f'{OUTB}/lr_pairs_allcells_per_section_cn.csv'; PAIR84 = f'{OUTB}/lr_pairs_v2_per_section_cn.csv'
    OUTFILES = [f'{SSTEM}.png', f'{SSTEM}.csv', f'{SSTEM}_legend.txt']
    assert all('_costim_' in os.path.basename(f) for f in OUTFILES) and not set(OUTFILES) & {DOT82C, HEAT82C, SUB82C, PR84}   # outputs never overwrite an input
    MIN_CELLS_CN = 100; CNS = list(range(1, 11)); REPORT_CNS = [10, 4, 8]
    AP, BROAD, IC, MI = 'Antigen presentation', 'Co-stimulation / suppressive Treg', 'Inflammatory cytokine', 'Matrix-integrin'
    CS9, IR, IMM = 'Co-Stimulatory (Activating)', 'Immunoregulatory (Inhibitory / Suppressive)', 'Immunomodulatory'
    CS7 = 'Co-Stimulatory (APC -> T cell)'
    MODS = [AP, IC, CS7, MI]          # shown, in this order (Inflammatory Cytokine before Co-Stimulatory)
    REF = [IR, CS9, IMM, BROAD]       # CSV only (not plotted)
    MOD4 = [AP, BROAD, IC, MI]
    # ---------------- membership: cn_lr_modules_allcells (unchanged modules) and lr_costim_inhibitory_modules (pre-specified sub-modules), asserted against the definition
    m82 = pd.read_csv(f'{OUTB}/lr_modules_allcells_membership.csv'); assert not m82.lr_pair.duplicated().any()
    assert m82.module.value_counts().reindex(MOD4).tolist() == [9, 50, 34, 227]
    vmem = pd.read_csv(f'{OUTB}/lr_modules_v2_membership.csv')
    CS_P = vmem[vmem.module == CS9].lr_pair.tolist(); IR_P = vmem[vmem.module == IR].lr_pair.tolist()
    CS_EXP = ['CD80^CD28', 'CD86^CD28', 'ICOSLG^CD28', 'ICOSLG^ICOS', 'TNFSF4^TNFRSF4', 'TNFSF18^TNFRSF18', 'TNFSF9^TNFRSF9', 'CD70^CD27', 'CD40LG^CD40']
    CHK_EXP = ['CD274^PDCD1', 'PDCD1LG2^PDCD1', 'CD80^CTLA4', 'CD86^CTLA4', 'NECTIN2^TIGIT', 'NECTIN2^CD96', 'PVR^TIGIT', 'PVR^CD96', 'LGALS9^HAVCR2',
               'TNFRSF14^BTLA', 'TNFRSF14^CD160']
    CYT_EXP = ['IL10^IL10RA', 'IL10^IL10RB'] + [f'TGFB{i}^TGFBR{j}' for i in (1, 2, 3) for j in (1, 2)]
    assert sorted(CS_P) == sorted(CS_EXP) and len(CS_P) == 9
    assert sorted(IR_P) == sorted(CHK_EXP + CYT_EXP) and len(IR_P) == 19
    IMM_P = CS_P + IR_P; assert len(set(IMM_P)) == 28
    assert not any(p.split('^')[0] in ('IL2', 'IL7', 'IL15') for p in IMM_P)          # gamma-c cytokines not included
    # canonical APC -> T-cell co-stimulation: ligand on the APC, receptor on the T cell / Treg
    CS7_P = ['CD80^CD28', 'CD86^CD28', 'ICOSLG^ICOS', 'TNFSF4^TNFRSF4', 'TNFSF18^TNFRSF18', 'TNFSF9^TNFRSF9', 'CD70^CD27']
    CS_DROP = {'CD40LG^CD40': 'T cell -> APC direction (CD40L licenses APCs; does not co-stimulate T cells)', 'ICOSLG^CD28': 'non-canonical interaction'}
    assert len(CS7_P) == 7 and sorted(CS7_P + list(CS_DROP)) == sorted(CS_P) and not set(CS7_P) & set(CS_DROP)
    assert all(p.split('^')[1] in ('CD28', 'ICOS', 'TNFRSF4', 'TNFRSF18', 'TNFRSF9', 'CD27') for p in CS7_P)   # T-cell-side receptors only
    PAIRS = {IMM: IMM_P, CS9: CS_P, IR: IR_P, CS7: CS7_P}
    NPAIR = {AP: 9, BROAD: 50, IC: 34, MI: 227, IMM: 28, CS9: 9, IR: 19, CS7: 7}
    dropped = sorted(set(m82[m82.module == BROAD].lr_pair) - set(IMM_P)); added = sorted(set(IMM_P) - set(m82[m82.module == BROAD].lr_pair))
    print('broad co-stimulation pairs not in Immunomodulatory (%d):' % len(dropped), dropped); print('pairs added (%d):' % len(added), added, flush=True)
    # ---------------- (1) the four modules of cn_lr_modules_allcells, from its per-cell cache
    cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas(); assert not cn.index.duplicated().any()
    C = pd.read_parquet(CACHE82, columns=[f'raw|{m}' for m in MOD4])
    X = C.join(cn[['CN_k10', 'section']], how='inner').rename(columns={'CN_k10': 'CN'})
    NCELL = len(X); print('cells with LIANA scores:', len(C), '| with a CN:', NCELL, flush=True)
    G = X.groupby(['section', 'CN'])
    V82 = (G[[f'raw|{m}' for m in MOD4]].mean().set_axis(MOD4, axis=1).reset_index().melt(id_vars=['section', 'CN'], var_name='module', value_name='raw')
           .merge(G.size().rename('n_cells').reset_index(), on=['section', 'CN']))
    del C, X; assert V82.raw.dtype == np.float32   # precision path of cn_lr_modules_allcells (float32 per-cell cache -> float32 CN x section means)
    f82 = pd.read_csv(f'{OUTB}/lr_modules_allcells_per_section_cn.csv').merge(V82, on=['section', 'CN', 'module'], suffixes=('_file', ''))
    assert len(f82) == len(V82) and (f82.n_cells == f82.n_cells_file).all() and np.allclose(f82.raw, f82.raw_file, rtol=1e-5, atol=1e-9)   # file: 6 sig figs
    # ---------------- (2) shortcut check: mean over pairs of per-pair CN x section means == the module values of cn_lr_modules_allcells (all four modules)
    L83 = pd.read_csv(PAIR83); assert set(L83.lr_pair) == set(m82.lr_pair)
    for m in MOD4:
        s = L83[L83.module == m].groupby(['section', 'CN']).agg(raw=('raw', 'mean'), n=('lr_pair', 'size'), n_cells=('n_cells', 'first'))
        assert (s.n == NPAIR[m]).all()
        r = V82[V82.module == m].set_index(['section', 'CN']).reindex(s.index)
        d = float((s.raw - r.raw).abs().max()); assert d < 1e-6 and (s.n_cells == r.n_cells).all(), (m, d)
        print(f'shortcut reproduces cn_lr_modules_allcells for {m}: max abs diff {d:.2e} over {len(s)} CN x section values', flush=True)
    L = pd.read_csv(PAIR84); assert set(IMM_P) <= set(L.lr_pair)
    ov = L.merge(L83, on=['section', 'CN', 'lr_pair'], suffixes=('', '_83'))
    assert ov.lr_pair.nunique() >= 14 and (ov.n_cells == ov.n_cells_83).all() and np.allclose(ov.raw, ov.raw_83, rtol=1e-7, atol=1e-10)
    print('lr_costim_inhibitory_modules pair cache == cn10_lr_pairs pair cache for', ov.lr_pair.nunique(), 'shared pairs', flush=True)
    # ---------------- (3) Immunomodulatory, its sub-modules and the 7-pair Co-Stimulatory module: mean over pairs per CN x section
    VN = []
    for m, ps in PAIRS.items():
        s = L[L.lr_pair.isin(ps)].groupby(['section', 'CN']).agg(raw=('raw', 'mean'), n=('lr_pair', 'size'), n_cells=('n_cells', 'first')).reset_index()
        assert (s.n == len(ps)).all(); VN.append(s.drop(columns='n').assign(module=m))
    VN = pd.concat(VN, ignore_index=True)
    chk = VN[VN.module == IMM].merge(V82[V82.module == AP], on=['section', 'CN'], suffixes=('', '_ap')); assert len(chk) == (V82.module == AP).sum() and (chk.n_cells == chk.n_cells_ap).all()
    VN['raw'] = VN.raw.astype(np.float32)   # same float32 precision path as cn_lr_modules_allcells (the LMM optimiser fallback is sensitive to it at the variance boundary)
    V = pd.concat([V82[['section', 'CN', 'module', 'raw', 'n_cells']], VN[['section', 'CN', 'module', 'raw', 'n_cells']]], ignore_index=True)
    assert V.raw.dtype == np.float32
    V = V[V.n_cells >= MIN_CELLS_CN].copy()
    zs = lambda x: (x - x.mean()) / x.std()
    V['z_raw'] = V.groupby(['section', 'module']).raw.transform(zs); V.loc[~np.isfinite(V.z_raw), 'z_raw'] = np.nan
    print('sections without variation (z undefined) per module:', V[V.z_raw.isna()].groupby('module').section.unique().to_dict(), flush=True)
    def fit(dd):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm('y ~ in_cn', dd, groups=dd.section).fit(reml=True, method=meth)
                if np.isfinite(f.bse['in_cn']) and f.bse['in_cn'] > 0: return f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn'], meth
            except Exception: pass
        return np.nan, np.nan, np.nan, 'failed'
    # ---------------- CN names / labels (as cn_lr_modules_allcells)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [(' Fib / PLA2G2A⁺ Fib / Endothelial', '/PLA2G2A⁺ Fib'), (' / Crypt Top Fib / Cycling', ' / CryptTop / Cycl'),
            (' / ADAMDEC1 Fib', ' / ADAMDEC1'), (' / Neuron / Pericyte', ' / Neuron / Peri'), (' / Mac S_M_P / CD9 Mac', ' / Mac / CD9 Mac'),
            ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'), ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'),
            ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {c: short_cn(NM[c]) for c in CNS}
    ST = []
    for m in MODS + REF:
        d = V[V.module == m].dropna(subset=['z_raw'])
        for c in CNS:
            dd = d.assign(y=d.z_raw, in_cn=(d.CN == c).astype(float)); zc = dd[dd.CN == c].y.values
            b, se, p, meth = fit(dd) if len(zc) >= 3 else (np.nan,) * 3 + ('not fitted',)
            t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
            ST.append(dict(module=m, CN=c, CN_label=LAB[c], CN_name=NM[c], n_pairs=NPAIR[m], n_sections=len(zc), mean_within_section_z=zc.mean(),
                           mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)), mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)),
                           sections_above_0=int((zc > 0).sum()), wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan,
                           lmm_estimate=b, lmm_ci_low=b - 1.96 * se, lmm_ci_high=b + 1.96 * se, lmm_p=p, optimiser=meth))
    ST = pd.DataFrame(ST); ST['lmm_q'] = np.nan; ST['wilcoxon_q'] = np.nan
    for m, g in ST.groupby('module'):
        for k_ in ['lmm', 'wilcoxon']:
            ok = g[f'{k_}_p'].notna(); ST.loc[g.index[ok], f'{k_}_q'] = multipletests(g[f'{k_}_p'][ok], method='fdr_bh')[1]
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    cat = lambda q: np.digitize(q, [0.001, 0.01, 0.05])
    # module names in the cn_lr_modules_immunomodulatory outputs (for the asserts) and in this section's outputs
    SHORTM_IMM = {AP: 'Antigen Presentation', IMM: 'Immunomodulatory', IC: 'Inflammatory Cytokine', MI: 'Matrix–Integrin',
                  BROAD: 'Co-Stimulation / Suppressive Treg', CS9: 'Co-Stimulatory (Activating)', IR: 'Immunoregulatory (Inhibitory / Suppressive)'}
    SHORTM = {AP: 'Antigen Presentation', CS7: 'Co-Stimulatory', IC: 'Inflammatory Cytokine', MI: 'Matrix–Integrin', IR: 'Inhibitory / Suppressive',
              CS9: 'Co-Stimulatory (9 pairs, reference)', IMM: 'Immunomodulatory (28 pairs, reference)', BROAD: 'Co-Stimulation / Suppressive Treg (50 pairs, reference)'}
    # ---------------- checks against the cn_lr_modules_allcells heatmap (3 sig figs), as in cn_lr_modules_immunomodulatory
    H82 = pd.read_csv(f'{OUTA}/supp_lr_modules_allcells_heatmap.csv')
    for m in MOD4:
        a = ST[ST.module == m].set_index('CN').loc[CNS]; b = H82[H82.module == SHORTM_IMM[m]].set_index('CN').loc[CNS]
        assert np.allclose(a.mean_within_section_z, b.mean_within_section_z, rtol=5e-3, atol=1e-3), m
        assert (a.n_sections.values == b.n_sections.values).all() and (cat(a.lmm_q.values) == cat(b.lmm_q.values)).all(), m
        assert np.allclose(a.lmm_q, b.lmm_q, rtol=5e-3, atol=0), (m, float(np.abs(np.log10(a.lmm_q.values) - np.log10(b.lmm_q.values)).max()))
    print('Antigen Presentation, Inflammatory Cytokine, Matrix-Integrin and the broad module reproduce cn_lr_modules_allcells (z and q, 3 sig figs)', flush=True)
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    roundtrip = lambda df: pd.read_csv(io.StringIO(df.apply(sig3).to_csv(index=False)))   # exactly what a CSV written by cn_lr_modules_immunomodulatory holds
    def same(a, b, cols, what):
        assert len(a) == len(b) and len(a) > 0, (what, len(a), len(b))
        for c in cols:
            x, y = a[c].reset_index(drop=True), b[c].reset_index(drop=True)
            if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
                ok = np.array_equal(x.astype(float).values, y.astype(float).values, equal_nan=True)
            else:
                ok = (x.fillna('').astype(str).values == y.fillna('').astype(str).values).all()
            assert ok, (what, c, pd.concat([x, y], axis=1)[x.astype(str) != y.astype(str)].head())
        print(f'  identical to cn_lr_modules_immunomodulatory: {what} ({len(a)} rows x {len(cols)} columns)', flush=True)
    # ---------------- (A) 9-pair Co-Stimulatory, Inhibitory / Suppressive and merged Immunomodulatory == supp_lr_modules_allcells_clean_submodules.csv (CN4 / CN8 / CN10)
    S82 = pd.read_csv(SUB82C)
    SB = ST[ST.module.isin([CS9, IR, IMM]) & ST.CN.isin(REPORT_CNS)].copy()
    SB['module'] = pd.Categorical(SB.module.map({CS9: 'Co-Stimulatory (Activating)', IR: 'Immunoregulatory (Inhibitory / Suppressive)', IMM: 'Immunomodulatory (merged)'}),
                                  ['Co-Stimulatory (Activating)', 'Immunoregulatory (Inhibitory / Suppressive)', 'Immunomodulatory (merged)'])
    SB['CN'] = pd.Categorical(SB.CN, REPORT_CNS); SB = SB.sort_values(['CN', 'module']); SB['stars'] = SB.lmm_q.map(stars)
    SB = roundtrip(SB[[c for c in S82.columns if c != 'note']])
    print('assert vs', os.path.basename(SUB82C), flush=True)
    same(SB, S82, [c for c in S82.columns if c != 'note'], '9-pair Co-Stimulatory / Inhibitory / merged, CN4 CN8 CN10, all columns incl. lmm_q and stars')
    # ---------------- (B) every module of cn_lr_modules_immunomodulatory == its dot-plot CSV (summary rows, all CNs; section rows); unchanged modules == its heatmap CSV
    O82 = pd.read_csv(DOT82C)
    SUMC = [c for c in O82.columns if c not in ('module', 'shown_in_figure', 'row_type', 'section', 'raw', 'n_cells', 'within_section_z', 'metric', 'model')]
    SECC = ['section', 'CN', 'raw', 'n_cells', 'within_section_z']
    print('assert vs', os.path.basename(DOT82C), flush=True)
    for m in [AP, IC, MI, CS9, IR, IMM, BROAD]:
        a = roundtrip(ST[ST.module == m].sort_values('CN')[SUMC]); b = O82[(O82.row_type == 'summary') & (O82.module == SHORTM_IMM[m])].sort_values('CN')
        same(a, b, SUMC, f'{SHORTM_IMM[m]} summary')
        a = V[V.module == m].rename(columns={'z_raw': 'within_section_z'})[SECC].sort_values(['section', 'CN'])
        for c_ in ['n_cells']: a[c_] = a[c_].round().astype('Int64')
        a = roundtrip(a); b = O82[(O82.row_type == 'section') & (O82.module == SHORTM_IMM[m])].sort_values(['section', 'CN'])
        same(a, b, SECC, f'{SHORTM_IMM[m]} section rows')
    HM82 = pd.read_csv(HEAT82C)
    print('assert vs', os.path.basename(HEAT82C), flush=True)
    for m in [AP, IC, MI]:
        a = ST[ST.module == m].sort_values('CN').copy(); a['cell_text'] = [f'{e:+.2f}{stars(q)}' for e, q in zip(a.mean_within_section_z, a.lmm_q)]
        cols = [c for c in HM82.columns if c not in ('module', 'shown_in_figure', 'metric', 'model')]
        same(roundtrip(a[cols]), HM82[HM82.module == SHORTM_IMM[m]].sort_values('CN'), cols, f'{SHORTM_IMM[m]} heatmap rows')
    print('all asserts vs cn_lr_modules_immunomodulatory passed before the 7-pair Co-Stimulatory set is used', flush=True)
    # ---------------- per-pair co-stimulatory values (per-pair method of lr_costim_inhibitory_modules; reported in the CSV, not plotted)
    DET = pd.read_csv(f'{OUTB}/lr_pairs_v2_detection.csv').set_index('lr_pair')
    LP = L[L.lr_pair.isin(CS_P)].copy(); LP = LP[LP.n_cells >= MIN_CELLS_CN].copy()
    LP['z'] = LP.groupby(['section', 'lr_pair']).raw.transform(zs); LP.loc[~np.isfinite(LP.z), 'z'] = np.nan
    PP = []
    for p in CS_P:
        d2 = LP[LP.lr_pair == p].dropna(subset=['z'])
        for c in CNS:
            zc = d2[d2.CN == c].z.values
            PP.append(dict(module='Co-Stimulatory pair (not plotted)', lr_pair=p, pair_label=p.replace('^', ' → '),
                           pair_status='kept (in the 7-pair module)' if p in CS7_P else 'removed: ' + CS_DROP[p].replace('->', '→'),
                           frac_cells_nonzero=DET.loc[p, 'frac_cells_nonzero'], CN=c, CN_label=LAB[c], CN_name=NM[c], n_sections=len(zc),
                           mean_within_section_z=zc.mean() if len(zc) else np.nan, sections_above_0=int((zc > 0).sum())))
    PP = pd.DataFrame(PP)
    P84 = pd.read_csv(PR84).set_index('lr_pair')
    for p in CS_P:
        a = roundtrip(PP[PP.lr_pair == p].set_index('CN')[['mean_within_section_z', 'n_sections']].T.reset_index(drop=True))
        for c, col in [(10, 'cn10_mean_within_section_z'), (8, 'cn8_mean_within_section_z'), (4, 'cn4_mean_within_section_z')]:
            assert float(a.loc[0, str(c)]) == float(P84.loc[p, col]), (p, c, float(a.loc[0, str(c)]), float(P84.loc[p, col]))
        assert int(a.loc[1, '10']) == int(P84.loc[p, 'cn10_n_sections']) and float(f"{DET.loc[p, 'frac_cells_nonzero']:.3g}") == float(P84.loc[p, 'frac_cells_nonzero'])
    print('per-pair co-stimulatory values reproduce lr_modules_v2_cn10_pairs.csv (CN4 / CN8 / CN10 z, CN10 n sections, detection)', flush=True)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
    PT = PP[PP.CN.isin([4, 8, 10])].pivot(index='lr_pair', columns='CN', values='mean_within_section_z').loc[CS_P]
    PT.insert(0, 'pct_cells_nonzero', [100 * DET.loc[p, 'frac_cells_nonzero'] for p in CS_P]); PT.insert(0, 'kept', [p in CS7_P for p in CS_P])
    print(PT.round(3).to_string(), flush=True)
    print(ST.pivot(index='CN', columns='module', values='mean_within_section_z').round(2).astype(str).add(
          ST.pivot(index='CN', columns='module', values='lmm_q').applymap(stars)).loc[CNS, MODS + REF].to_string(), flush=True)
    print(ST[ST.module.isin([CS7, IR, CS9])].assign(stars=lambda d: d.lmm_q.map(stars))[
          ['module', 'CN', 'n_sections', 'mean_within_section_z', 'lmm_estimate', 'lmm_p', 'lmm_q', 'stars', 'optimiser']].round(4).to_string(index=False), flush=True)
    # ---------------- per-section dot plot (no bold); stars on the side of 0 matching the direction
    S_ = ST[ST.module.isin(MODS)].copy()
    assert S_.mean_within_section_z.notna().all() and S_.lmm_q.notna().all()
    sg = S_[S_.lmm_q < 0.05]
    assert (np.sign(sg.mean_within_section_z) == np.sign(sg.lmm_estimate)).all(), sg   # direction is unambiguous for every star
    COL = {AP: ('#9ecae1', '#08519c'), CS7: ('#bcbddc', '#54278f'), IC: ('#fcae91', '#a50f15'), MI: ('#a1d99b', '#006d2c')}
    TITLE = {AP: 'Antigen Presentation', CS7: 'Co-Stimulatory', IC: 'Inflammatory Cytokine', MI: 'Matrix–Integrin'}
    plt.rcParams.update({'font.size': 15, 'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(4, 1, figsize=(7.6, 12), sharex=True, gridspec_kw={'hspace': 0.24})
    STAR_FS = 17; FP = FontProperties(family='DejaVu Sans', size=STAR_FS)
    GLY = {s: TextPath((0, 0), s, size=STAR_FS, prop=FP).get_extents() for s in ['*', '**', '***']}   # glyph extents, pt relative to the baseline
    DOT_S, MEAN_MS, MEAN_MEW, ELW = 45, 10, 2.2, 2.6
    R_DOT = np.sqrt(DOT_S) / 2 + 0.25; R_MEAN = MEAN_MS / 2 + MEAN_MEW / 2; R_CI = ELW / 2       # pt beyond the data value (dot / mean marker / CI end)
    GAP, MARG = 2.0, 2.0                                                                         # pt: star-to-element gap; star-to-axes-edge margin
    rng2 = np.random.default_rng(1); ART = []
    for ax, m in zip(axes, MODS):
        lc, dc = COL[m]; g = S_[S_.module == m].set_index('CN'); vv = V[V.module == m]
        ax.axhline(0, color='#888888', lw=1, zorder=0)
        pts = {}
        for i, c in enumerate(CNS):
            z = vv[vv.CN == c].z_raw.dropna().values; xj = i + rng2.uniform(-0.14, 0.14, len(z))
            ax.scatter(xj, z, s=DOT_S, color=lc, edgecolor='white', linewidth=0.5, alpha=0.95, zorder=2)
            e, lo, hi, q = g.loc[c, ['mean_within_section_z', 'mean_ci_low', 'mean_ci_high', 'lmm_q']]
            ax.errorbar(i, e, yerr=[[e - lo], [hi - e]], fmt='o', ms=MEAN_MS, color=dc, mfc=dc if q < 0.05 else 'white', mew=MEAN_MEW, elinewidth=ELW, capsize=0, zorder=3)
            pts[c] = (xj, z, e, lo, hi, q)
        vals = np.concatenate([vv.z_raw.dropna().values, g.mean_ci_low.values, g.mean_ci_high.values]); ymin, ymax = vals.min(), vals.max(); span = ymax - ymin
        b0, b1 = ymin - 0.05 * span, ymax + 0.05 * span            # base limits; widened below only as far as the stars need
        H_pt = ax.get_window_extent().height * 72 / fig.dpi
        sigc = [c for c in CNS if pts[c][5] < 0.05]
        def place(y0, y1):
            ppd = H_pt / (y1 - y0); out = {}                        # points per data unit
            for c in sigc:
                xj, z, e, lo, hi, q = pts[c]; s = stars(q); gb = GLY[s]
                if e > 0:
                    edge = max(z.max() + R_DOT / ppd, hi + R_CI / ppd, e + R_MEAN / ppd); off = GAP - gb.y0
                    out[c] = (s, edge, off, edge + (GAP) / ppd, edge + (GAP + gb.y1 - gb.y0) / ppd)
                else:
                    edge = min(z.min() - R_DOT / ppd, lo - R_CI / ppd, e - R_MEAN / ppd); off = -GAP - gb.y1
                    out[c] = (s, edge, off, edge - (GAP + gb.y1 - gb.y0) / ppd, edge - GAP / ppd)
            return out, ppd
        y0, y1 = b0, b1
        for it in range(100):
            out, ppd = place(y0, y1)
            n0 = min([b0] + [o[3] - MARG / ppd for o in out.values()]); n1 = max([b1] + [o[4] + MARG / ppd for o in out.values()])
            if abs(n0 - y0) < 1e-9 and abs(n1 - y1) < 1e-9: break
            y0, y1 = n0, n1
        assert it < 99
        ax.set_ylim(y0, y1); out, ppd = place(y0, y1)
        for i, c in enumerate(CNS):
            if c in out:
                s, edge, off, gy_lo, gy_hi = out[c]
                t = ax.annotate(s, xy=(i, edge), xytext=(0, off), textcoords='offset points', ha='center', va='baseline', fontsize=STAR_FS, color=COL[m][1],
                                annotation_clip=False)
                ART.append(dict(ax=ax, module=m, CN=c, i=i, text=t, s=s, edge=edge, off=off, sign=np.sign(pts[c][2]), pts=pts[c]))
        ax.set_title(TITLE[m], fontsize=17, loc='left', pad=3); ax.tick_params(axis='y', labelsize=15)
        ax.spines[['top', 'right']].set_visible(False); ax.set_xlim(-0.55, len(CNS) - 0.45)
        print(f'{TITLE[m]}: y-limits {y0:.3f} .. {y1:.3f} (base {b0:.3f} .. {b1:.3f}); stars at CN', {c: (out[c][0], '+' if pts[c][2] > 0 else '-') for c in sigc}, flush=True)
    axes[-1].set_xticks(range(len(CNS))); axes[-1].set_xticklabels([LAB[c] for c in CNS], rotation=45, ha='right', rotation_mode='anchor', fontsize=13)
    axes[-1].set_xlabel('Spatial Neighborhood', fontsize=16)
    sy = fig.supylabel('LR Module Score (Within-Section z)', fontsize=16, x=0.04)
    h = [Line2D([], [], marker='o', ls='', color='#bdbdbd', ms=9, label='One Section'),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='#252525', ms=10, lw=2.4, label='Mean ± 95% CI, q < 0.05'),
         Line2D([], [], ls='', label=' '),
         Line2D([], [], marker='o', ls='-', color='#252525', mfc='white', ms=10, mew=2, lw=2.4, label='q ≥ 0.05')]
    leg = axes[0].legend(handles=h, loc='lower left', bbox_to_anchor=(-0.02, 1.12), ncol=2, frameon=False, fontsize=14, handletextpad=0.3, columnspacing=0.8, borderaxespad=0)
    # ---------------- check of the star placement and text layout (display coordinates)
    from matplotlib.transforms import Bbox
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); px = fig.dpi / 72
    ov_ = lambda a, b: a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1
    circ_hits = lambda cx, cy, r, bb: (max(bb.x0 - cx, 0, cx - bb.x1) ** 2 + max(bb.y0 - cy, 0, cy - bb.y1) ** 2) < r ** 2
    TITLES = [a.title.get_window_extent(rd) for a in axes]; LEGB = leg.get_window_extent(rd); AXB = [a.get_window_extent(rd) for a in axes]
    YTL = [t.get_window_extent(rd) for a in axes for t in a.get_yticklabels() if t.get_visible() and t.get_text()]
    XTL = [t.get_window_extent(rd) for t in axes[-1].get_xticklabels()]
    CLEAR = []
    for r in ART:
        ax = r['ax']; gb = GLY[r['s']]; tb = r['text'].get_window_extent(rd); AB = ax.get_window_extent(rd)
        cx, cy = ax.transData.transform((r['i'], r['edge'])); base = cy + r['off'] * px; hw = (gb.x1 - gb.x0) / 2 * px
        assert abs((tb.x0 + tb.x1) / 2 - cx) < 1.0, ('star not centred on its CN', r['module'], r['CN'])
        G_ = Bbox([[cx - hw, base + gb.y0 * px], [cx + hw, base + gb.y1 * px]])          # the asterisk glyphs themselves
        assert G_.y0 >= AB.y0 + 1 and G_.y1 <= AB.y1 - 1 and G_.x0 >= AB.x0 and G_.x1 <= AB.x1, ('star outside axes', r['module'], r['CN'])
        y0d = ax.transData.transform((0, 0))[1]
        assert (G_.y0 > y0d) if r['sign'] > 0 else (G_.y1 < y0d), ('star on wrong side of 0', r['module'], r['CN'])
        hits = []; dmin = np.inf
        for coll in ax.collections:                                                       # section dots (PathCollections) and CI bars (LineCollection)
            if hasattr(coll, 'get_segments'):
                for seg in coll.get_segments():
                    (sx, sy0), (_, sy1) = ax.transData.transform(np.asarray(seg)[[0, -1]])
                    SB_ = Bbox([[sx - R_CI * px, min(sy0, sy1) - R_CI * px], [sx + R_CI * px, max(sy0, sy1) + R_CI * px]])
                    if ov_(SB_, G_): hits.append(('CI bar', sx))
                    if abs(sx - cx) < hw: dmin = min(dmin, (G_.y0 - SB_.y1) if r['sign'] > 0 else (SB_.y0 - G_.y1))
            else:
                for (x_, y_) in ax.transData.transform(coll.get_offsets()):
                    if circ_hits(x_, y_, R_DOT * px, G_): hits.append(('section dot', x_, y_))
                    if abs(x_ - cx) < hw + R_DOT * px: dmin = min(dmin, (G_.y0 - (y_ + R_DOT * px)) if r['sign'] > 0 else ((y_ - R_DOT * px) - G_.y1))
        for ln in ax.lines:                                                               # mean markers
            if ln.get_marker() == 'o' and len(ln.get_xdata()) == 1:
                x_, y_ = ax.transData.transform((ln.get_xdata()[0], ln.get_ydata()[0]))
                if circ_hits(x_, y_, R_MEAN * px, G_): hits.append(('mean marker', x_, y_))
        assert not hits, ('star overlaps', r['module'], r['CN'], hits)
        assert not any(ov_(G_, b) for b in TITLES + [LEGB] + YTL + XTL), ('star overlaps a title / legend / tick label', r['module'], r['CN'])
        CLEAR.append(dmin / px)
    for k in range(1, len(axes)):                                                         # titles sit between panels, clear of the panel above
        assert TITLES[k].y0 >= AXB[k].y1 - 0.5 and TITLES[k].y1 < AXB[k - 1].y0, ('title overlaps the panel above', k)
    assert LEGB.y0 > TITLES[0].y1 and not ov_(LEGB, TITLES[0]), 'legend overlaps the first title'
    assert not any(ov_(sy.get_window_extent(rd), b) for b in YTL), 'y label overlaps tick labels'
    print(f'star placement check passed: {len(ART)} stars inside axes, on the correct side of 0, no overlap with section dots / mean markers / CI bars / titles / '
      f'legend / tick labels; clearance to the nearest element in the column {min(CLEAR):.1f}-{max(CLEAR):.1f} pt', flush=True)
    fig.savefig(f'{SSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    # ---------------- source data (CSV): plotted modules + reference modules (not plotted) + per-pair co-stimulatory values
    METRIC = 'LIANA bivariate local score (raw), mean over module pairs per cell, mean over ALL cells of the CN per section; z across CNs within section, mean over sections'
    MODEL = 'z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within module'
    PAIRS_TXT = {m: '; '.join(p.replace('^', ' → ') for p in ps) for m, ps in [(CS7, CS7_P), (CS9, CS_P), (IR, IR_P), (IMM, IMM_P)]}
    SUM = ST.assign(shown_in_figure=ST.module.isin(MODS), stars=ST.lmm_q.map(stars), lr_pairs=ST.module.map(PAIRS_TXT), module=ST.module.map(SHORTM), row_type='summary')
    SEC = V.assign(module=V.module.map(SHORTM), shown_in_figure=V.module.isin(MODS), row_type='section').rename(columns={'z_raw': 'within_section_z'})
    PPO = PP.assign(shown_in_figure=False, row_type='pair_summary', n_pairs=1,
                    metric='per pair: LIANA bivariate local score (raw), mean over ALL cells of the CN per section; z across CNs within section, mean over sections; '
                       'frac_cells_nonzero = fraction of all cells with a non-zero local score', model='')
    OUTS = pd.concat([SUM.assign(metric=METRIC, model=MODEL), SEC.assign(metric=METRIC, model=MODEL), PPO], ignore_index=True)
    for c_ in ['n_cells', 'n_pairs', 'n_sections', 'sections_above_0']: OUTS[c_] = OUTS[c_].round().astype('Int64')   # counts stay exact
    OUTS.apply(sig3).to_csv(f'{SSTEM}.csv', index=False)
    # ---------------- legend
    DET_ = DET.frac_cells_nonzero
    pct = lambda ps: (100 * DET_.loc[ps].min(), 100 * DET_.loc[ps].max())
    cd8086 = pct(['CD80^CD28', 'CD86^CD28', 'CD80^CTLA4', 'CD86^CTLA4']); ox40 = pct(['TNFSF4^TNFRSF4']); il2 = pct(['IL2^IL2RA', 'IL2^IL2RB', 'IL2^IL2RG'])
    low = pct(['TGFB2^TGFBR1', 'TGFB2^TGFBR2', 'TNFRSF14^CD160'])
    def grp(ps):
        out = {}
        for p in ps: l, r = p.split('^'); out.setdefault(l, []).append(r)
        return '; '.join(f'{l} → {", ".join(rs)}' for l, rs in out.items())
    legS = ('Section-level spatial ligand-receptor module activity across cellular neighbourhoods, '
        f'CD resections (CosMx 6K; 11 sections, 4 patients; all {NCELL:,} cells with a neighbourhood). Ligand-receptor (LR) scores are cell-type-agnostic '
        'spatial LIANA bivariate local scores (Gaussian kernel, 200 px bandwidth) for the NicheNet v2 prior-network pairs, computed for every cell. '
        f'Modules (in panel order): Antigen Presentation = CellChatDB MHC-II ({NPAIR[AP]} pairs); Inflammatory Cytokine = CellChatDB IL1, IL6, TNF, '
        f'IFN-II, IL12, IL17, OSM ({NPAIR[IC]} pairs); Co-Stimulatory = canonical antigen-presenting cell (APC) → T-cell '
        'co-stimulation, i.e. pre-specified pairs with the ligand on the APC and the receptor on the T cell / Treg '
        f'({NPAIR[CS7]} pairs): {grp(CS7_P)}; '
        f'Matrix–Integrin = CellChatDB ECM-Receptor pathways ({NPAIR[MI]} pairs) (for the CellChatDB modules, a pair joins a module when a reference '
        'interaction has the same ligand and an overlapping receptor). Two pairs of the pre-specified 9-pair co-stimulatory set are excluded: CD40LG → CD40, '
        'because it signals from the T cell to the APC (CD40L licenses APCs and does not co-stimulate T cells), and ICOSLG → CD28, a non-canonical '
        'interaction; values for the 9-pair set are in the source data for comparison. Co-stimulatory and inhibitory / suppressive pairs, which together '
        f'form the "Immunomodulatory" module ({NPAIR[IMM]} pairs), are shown separately, because co-stimulation (signal 2 for T-cell and Treg '
        'activation) and inhibitory checkpoints / suppressive cytokines are biologically distinct. The Inhibitory / Suppressive module '
        f'({NPAIR[IR]} pre-specified pairs; inhibitory checkpoints ({len(CHK_EXP)}): {grp(CHK_EXP)}; suppressive cytokines ({len(CYT_EXP)}): '
        f'{grp(CYT_EXP)}) is not plotted; its values are in the source data. The γc cytokines IL-2, IL-7 and IL-15 are in neither module; CD200 → CD200R1 '
        'and HLA-E → KLRC1 (NKG2A) were pre-specified as inhibitory but cannot be scored because CD200R1 and KLRC1 are not on the 6K panel. Some pairs '
        f'are sparsely detected on the 6K panel: CD80 and CD86 pairs have a non-zero local score in only {cd8086[0]:.0f}–{cd8086[1]:.0f}% of cells and '
        f'OX40 (TNFSF4 → TNFRSF4) in {ox40[0]:.0f}%; IL-2 pairs (not included) in {il2[0]:.0f}–{il2[1]:.0f}%, and the TGFB2 and CD160 pairs '
        f'(Inhibitory / Suppressive) in {low[0]:.0f}–{low[1]:.0f}%. The per-pair detection and neighbourhood values of the nine co-stimulatory pairs are '
        'in the source data. These modules are used instead of the broad Co-Stimulation / Suppressive Treg module (CellChatDB CD80, CD86, PD-L1, PDL2, ICOS, OX40, '
        f'GITRL, CD70, CD40, IL2 and IL10 pathways; {NPAIR[BROAD]} pairs) because that module also included non-immunomodulatory cytokine pairs (the '
        'IL-20-family cytokines IL-19, IL-20, IL-22, IL-24 and IL-26; the γc-family cytokines IL-2, IL-4, IL-7, IL-9, IL-15 and IL-21; and TSLP) and '
        'CD40LG → integrin pairs. No sender or receiver cell types are imposed. '
        'A cell\'s module score is the mean raw local score over the module\'s pairs, and a neighbourhood\'s value in a section is the mean over all its '
        'cells, so the score reflects both ligand–receptor expression and the abundance of interacting cells (neighbourhood × section combinations with at '
        'least 100 cells). Values are z-standardised across neighbourhoods within each section and averaged over sections ("Mean Within-Section z"; '
        '0 = neighbourhood average). Significance: linear mixed model on the per-section values, z ~ in neighbourhood + (1 | section), BH-corrected over '
        'the 10 neighbourhoods within each module; * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. Stars sit on the side of 0 that matches '
        'the direction of the estimate: above the highest section value or upper CI bound for positive estimates, below the lowest section value or '
        'lower CI bound for negative estimates. Light dots: one section each. Dark points and bars: mean over sections ± 95% CI (t distribution); filled, '
        'q < 0.05. Neighbourhood names are abbreviated by their dominant cell types; full names are in the source data.')
    tl = legS.lower()
    assert not any(w in tl for w in ['activation score', 'quiescent', 'territor', 'to tregs', 'blunt', 'selenop⁺ fib', 'selenop+ fib'])
    open(f'{SSTEM}_legend.txt', 'w').write(legS + '\n')
    print('wrote', *OUTFILES, sep='\n  '); print('done', flush=True)


# ======================================================================================================================
# Section cn_inflammatory_axes: Pre-specified inflammatory signalling axes per CN, CN10 ligand sources and the link to the Treg Suppression Score
# Panels: Fig 5C, Supp 9A-B inputs (supp_inflammatory_axes_by_neighbourhood.csv + _legend.txt)
# ======================================================================================================================
if _run('cn_inflammatory_axes'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/86_cn10_inflammatory_axes.py'   # working folder of this step (relative paths resolve from here)
    """CN10 (OGN+RSPO3+ / PLA2G2A+ Fib neighbourhood) as a TNF-independent inflammatory stroma, and a pre-specified link to the Treg Suppression Score.
CD resections (CosMx 6K; 11 sections, 4 patients).
Motivation: in the pair-level breakdown of cn10_lr_pairs (A_fig4_panels/supp_cn10_lr_pairs.csv) IL6, IL1A, OSM and IL11 pairs are higher in CN10, but not
TNF (high in CN4 / CN8). These axes are linked to anti-TNF non-response (West 2017 PMID 28368383; Friedrich 2021 PMID 34675383; Smillie 2019 PMID 31348891;
Atreya 2000 PMID 10802717). All axes, ligands, receptors and tests below are set in the PRE-SPECIFIED block before any data are read.
PART 1  Axis heatmap. Axes = LIANA pairs present in combined/spot_lr_df.parquet (all members of the 'Inflammatory cytokine' module of cn_lr_modules_allcells):
          IL-6 = IL6->IL6R, IL6->IL6ST; IL-1 = IL1A->IL1R1, IL1A->IL1RAP, IL1B->IL1R1, IL1B->IL1RAP; OSM = OSM->OSMR, OSM->LIFR, OSM->IL6ST;
          IL-11 = IL11->IL11RA, IL11->IL6ST; TNF (contrast) = TNF->TNFRSF1A, TNF->TNFRSF1B.
        Computed like a module of cn_lr_modules_allcells: raw per-cell LIANA bivariate local score (NaN -> 0), mean over the axis pairs per cell; CN x section
        mean over all cells (inner join with cn_per_cell CN_k10; >= 100 cells); z across CNs within section (a section where the axis is constant gives
        undefined z and drops out); mean over sections; LMM z ~ in_CN + (1 | section) (optimiser fallback); BH over the 10 CNs within axis (Wilcoxon
        signed-rank q in the CSV).
        Asserted: (a) the per-pair cache of cn10_lr_pairs reproduces the module values of cn_lr_modules_allcells; (b) the direct per-cell axis computation
        from the parquet equals the mean of the per-pair CN x section means (linearity).
        Figure: A_fig4_panels/supp_inflammatory_axes_by_neighbourhood.{png,csv} + _legend.txt; TNF set apart by a gap.
PART 2  Ligand sources (CSV only): IL6, IL1A, IL1B, OSM, IL11, TNF; % of cells expressing (norm_log1p > 0) per cell_type_general, CN10 vs all
        other CNs, per section mean (sections with >= 20 cells of the type in both); share of the ligand-positive cells (and of the summed ligand
        expression) contributed by each cell type in CN10, per section mean. Expression: per-section Processed/norm_log1p.h5ad; cell types:
        combined/all_samples.parquet ('Absoprtive Epi' shown as 'BEST4+ Epi'); CN: cellular_neighborhoods_v2/cn_per_cell.parquet CN_k10.
        SELENOP+ Fib has no CN assignment and is therefore absent; SELENOP+ Mac is included.
        Output: A_fig4_panels/supp_inflammatory_ligand_sources_cn10.csv.
PART 3  Pre-specified link to the Treg Suppression Score (-D6 signature, z within section), 6,554 strict Tregs.
        Local availability of IL6, IL1A, OSM, IL11, TNF (TNF = control) = mean within-section z of the ligand over all cells within 50 um (self
        included, >= 5 cells), as in the per-Treg local table of 11_neighbourhoods_local_tables.py (asserted to reproduce its lig_IL6, lig_TNF and
        Treg Suppression Score).
        libsize_z = log(1 + total counts) per Treg (combined/raw_counts_strict_tregs.csv), z within section over Tregs.
        PRIMARY (all Tregs): score ~ local ligand + libsize_z + (1 | section); BH over the 5 ligands; per-section OLS sign agreement.
        CN10 (descriptive): same model in Tregs with CN_k10 == 10 only if >= 5 sections have >= 10 CN10 Tregs (otherwise reported as underpowered and
          not fitted); n Tregs, sections, estimate and 95% CI; no BH-based claims.
        EXPLORATORY: score ~ local ligand x in_CN10 + libsize_z + (1 | section), all Tregs (interaction and CN10-specific slope).
        Own receptor expression of the Tregs (descriptive): % expressing IL6R, IL6ST, IL1R1, OSMR, IL11RA, TNFRSF1B (additional co-receptors IL1RAP,
          LIFR, TNFRSF1A), CN10 vs other CNs, pooled and per-section mean.
        Outputs: B_resection_analyses/inflammatory_axes/; a forest plot A_fig4_panels/supp_inflammatory_ligands_treg_suppression.png only if a
        primary q < 0.05.
Methods text: A_fig4_panels/supp_inflammatory_axes_by_neighbourhood_legend.txt and panels/fig4_and_related/METHODS_inflammatory_axes.md.
SKIP_COMPUTE=1 redraws from the caches in B_resection_analyses/inflammatory_axes/."""
    import os, gc, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_NUM_THREADS') == '32'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker, seaborn as sns
    from anndata.io import read_elem
    from scipy import sparse, stats as sst
    from scipy.spatial import cKDTree
    from statsmodels.stats.multitest import multipletests
    import statsmodels.formula.api as smf
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUTA = f'{P4}/A_fig4_panels'
    LRC = f'{P4}/B_resection_analyses/lr_modules_cellchat'; CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'
    LRS = f'{P4}/B_resection_analyses/ligand_receptor_state'; OUTB = f'{P4}/B_resection_analyses/inflammatory_axes'; os.makedirs(OUTB, exist_ok=True)
    R = '/share/fsmresfiles/UC/AtoMx/CD_surgical_resection'; SPOT = f'{R}/combined/spot_lr_df.parquet'
    HSTEM = f'{OUTA}/supp_inflammatory_axes_by_neighbourhood'; SRCCSV = f'{OUTA}/supp_inflammatory_ligand_sources_cn10.csv'
    FSTEM = f'{OUTA}/supp_inflammatory_ligands_treg_suppression'; METHODS = f'{P4}/METHODS_inflammatory_axes.md'
    C1 = f'{OUTB}/axes_per_section_cn.csv'; C2 = f'{OUTB}/ligand_sources_per_section_celltype.csv'; C3 = f'{OUTB}/treg_local_inflammatory_ligands.csv.gz'
    SKIP = os.environ.get('SKIP_COMPUTE') == '1'
    UM = 0.12028; RAD = 50; MIN_NB = 5; MIN_CELLS_CN = 100; MINC = 20; MIN_TREG_SEC = 10; MIN_SEC_CN10 = 5; TARGET = 10; CNS = list(range(1, 11))
    # ================= PRE-SPECIFIED (set before any data are read)
    AXES = {'IL-6': ['IL6^IL6R', 'IL6^IL6ST'], 'IL-1': ['IL1A^IL1R1', 'IL1A^IL1RAP', 'IL1B^IL1R1', 'IL1B^IL1RAP'],
            'OSM': ['OSM^OSMR', 'OSM^LIFR', 'OSM^IL6ST'], 'IL-11': ['IL11^IL11RA', 'IL11^IL6ST'], 'TNF': ['TNF^TNFRSF1A', 'TNF^TNFRSF1B']}
    AXN = list(AXES); CONTRAST = 'TNF'
    LIG_SRC = ['IL6', 'IL1A', 'IL1B', 'OSM', 'IL11', 'TNF']             # Part 2
    LIG_T = ['IL6', 'IL1A', 'OSM', 'IL11', 'TNF']                       # Part 3 (TNF = control); 5 tests, BH over 5
    REC_T = ['IL6R', 'IL6ST', 'IL1R1', 'OSMR', 'IL11RA', 'TNFRSF1B']; REC_X = ['IL1RAP', 'LIFR', 'TNFRSF1A']   # Treg receptors (descriptive; REC_X additional)
    REC_CTX = ['IL6R', 'IL6ST', 'IL1R1', 'IL1RAP', 'OSMR', 'LIFR', 'IL11RA', 'TNFRSF1A', 'TNFRSF1B']        # Part 2 context: receptor side of the axes
    # =================
    FOLDER = {'43563_A8': 'Cos1_6KDC_strict43563_SB_A8_091725BP5_09_10_2025_13_41_29_437', '54033_B4': 'Cos1_6KDC_strict54033_B4_0715_11_08_2025_14_20_23_133',
              '43563_A11': 'Cos5_6KDC_strict43563_SB_A11_091725BP7_09_10_2025_13_42_42_775', '43563_A4': 'Cos5_6KDC_strict43563_SB_A4_091725BP8_09_10_2025_13_43_50_952',
              '11808_A4_Cos3': 'Cos3_6KDC_strict11808_A4_0716_11_08_2025_14_19_48_615', '43563_A1': 'Cos1_6KDC_strict43563_SB_A1_091725BP6_09_10_2025_13_54_55_807',
              '11808_A4_Cos4': 'Cos4_6KDC_strict11808_A4_0716_14_10_2025_12_34_38_727', '54033_B7': 'Cos4_6KDC_strict54033_B7_0716_11_08_2025_14_18_41_738',
              '54033_B5': 'Cos3_6KDC_strict54033_B5_0716_11_08_2025_14_17_56_277', '18888_A6': 'Cos2_6KDC_strict18888_A6_0715_Alan_Test_11_08_2025_14_23_09_489',
              '18888_A10': 'Cos2_6KDC_strict18888_A10_0715_Alan_data_test_11_08_2025_14_25_11_455'}
    DISP = {'OGN+RSPO3+ Fib': 'OGN⁺RSPO3⁺ Fib', 'ADAMDEC1+ Fib': 'ADAMDEC1⁺ Fib', 'PLA2G2A+ ECM Fib': 'PLA2G2A⁺ ECM Fib', 'CD9+ M2 Mac': 'CD9⁺ M2 Mac',
            'Absoprtive Epi': 'BEST4⁺ Epi', 'Gd T': 'γδ T', 'SELENOP+ Mac': 'SELENOP⁺ Mac'}
    zs = lambda x: (x - x.mean()) / x.std()
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    pd.set_option('display.width', 260); pd.set_option('display.max_rows', 400); pd.set_option('display.max_columns', 40)
    def fit(formula, d, term, groups='section'):
        for meth in ['lbfgs', 'bfgs', 'powell', 'nm', 'cg']:
            try:
                f = smf.mixedlm(formula, d, groups=d[groups]).fit(reml=True, method=meth)
                if np.isfinite(f.bse[term]) and f.bse[term] > 0: return f, meth
            except Exception: pass
        return None, 'failed'

    # ===================================================================== PART 1: axes
    memb = pd.read_csv(f'{LRC}/lr_modules_allcells_membership.csv'); MODOF = memb.set_index('lr_pair').module.to_dict()
    PAIRS = [p for ps in AXES.values() for p in ps]; AXOF = {p: a for a, ps in AXES.items() for p in ps}
    scols = set(pq.read_schema(SPOT).names)
    assert all(p in scols for p in PAIRS), [p for p in PAIRS if p not in scols]
    assert all(MODOF.get(p) == 'Inflammatory cytokine' for p in PAIRS)
    # (a) the per-pair cache of cn10_lr_pairs reproduces the cn_lr_modules_allcells module values
    L = pd.read_csv(f'{LRC}/lr_pairs_allcells_per_section_cn.csv'); V82 = pd.read_csv(f'{LRC}/lr_modules_allcells_per_section_cn.csv')
    chk = L.groupby(['section', 'CN', 'module']).agg(raw=('raw', 'mean'), n=('n_cells', 'first')).reset_index().merge(V82, on=['section', 'CN', 'module'], suffixes=('', '_82'))
    assert len(chk) == len(V82) and (chk.n == chk.n_cells).all() and np.allclose(chk.raw, chk.raw_82, rtol=1e-4, atol=1e-7), 'pair cache does not reproduce cn_lr_modules_allcells'
    print('pair cache reproduces cn_lr_modules_allcells module values:', len(chk), 'rows; max abs diff', float((chk.raw - chk.raw_82).abs().max()), flush=True)
    LA = L[L.lr_pair.isin(PAIRS)].assign(axis=lambda x: x.lr_pair.map(AXOF))
    VC = LA.groupby(['section', 'CN', 'axis']).agg(raw=('raw', 'mean'), n_cells=('n_cells', 'first'), n_pairs=('lr_pair', 'size')).reset_index()
    assert (VC.n_pairs == VC.axis.map({a: len(ps) for a, ps in AXES.items()})).all()
    # (b) direct per-cell axis scores, as cn_lr_modules_allcells
    if not SKIP:
        cn = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', 'section', '__index_level_0__']).to_pandas(); assert not cn.index.duplicated().any()
        A = pq.read_table(SPOT, columns=PAIRS + ['__index_level_0__']).to_pandas(); assert A.index.dtype == object and not A.index.duplicated().any() and len(A) > 2.9e6
        A = A[PAIRS].astype(np.float64).fillna(0.0)
        Cax = pd.DataFrame({a: A[ps].values.mean(1) for a, ps in AXES.items()}, index=A.index); del A; gc.collect()
        X = Cax.join(cn[['CN_k10', 'section']], how='inner').rename(columns={'CN_k10': 'CN'})
        print('cells with LIANA scores:', len(Cax), '| with a CN:', len(X), flush=True)
        G = X.groupby(['section', 'CN'])
        VD = G[AXN].mean().reset_index().melt(id_vars=['section', 'CN'], var_name='axis', value_name='raw').merge(G.size().rename('n_cells').reset_index(), on=['section', 'CN'])
        VD.to_csv(C1, index=False, float_format='%.8g'); del X, Cax, cn; gc.collect()
    VD = pd.read_csv(C1)
    m_ = VD.merge(VC, on=['section', 'CN', 'axis'], suffixes=('', '_cache'))
    assert len(m_) == len(VD) == len(VC) and (m_.n_cells == m_.n_cells_cache).all() and np.allclose(m_.raw, m_.raw_cache, rtol=1e-4, atol=1e-7), 'direct != cache'
    n_a = VD.drop_duplicates(['section', 'CN']).set_index(['section', 'CN']).n_cells.sort_index(); n_b = V82.drop_duplicates(['section', 'CN']).set_index(['section', 'CN']).n_cells.sort_index()
    assert n_a.index.equals(n_b.index) and (n_a.values == n_b.values).all(), 'cell counts differ from cn_lr_modules_allcells'
    print('direct per-cell axis values = mean of pair-cache values: max abs diff', float((m_.raw - m_.raw_cache).abs().max()), '| total cells', int(VD.drop_duplicates(['section', 'CN']).n_cells.sum()), flush=True)
    NCELLS = int(VD.drop_duplicates(['section', 'CN']).n_cells.sum())
    V = VD[VD.n_cells >= MIN_CELLS_CN].copy()
    V['z'] = V.groupby(['section', 'axis']).raw.transform(zs); V.loc[~np.isfinite(V.z), 'z'] = np.nan
    print('sections with defined z per axis:', V.dropna(subset=['z']).groupby('axis').section.nunique().reindex(AXN).to_dict())
    # CN names / labels (as cn_lr_modules_allcells)
    NM = pd.read_csv(f'{CND}/cn_treg_state_enrichment_all.csv'); NM = NM[(NM.definition == 'D6_tertile') & (NM.k == 10)].set_index('CN').CN_name.to_dict()
    NM = {k + 1: v.replace('OGN_RSPO3 Fib', 'OGN⁺RSPO3⁺ Fib').replace('PLA2G2A Fib', 'PLA2G2A⁺ Fib') for k, v in NM.items()}
    TRIM = [(' Fib / PLA2G2A⁺ Fib / Endothelial', '/PLA2G2A⁺ Fib'), (' / Crypt Top Fib / Cycling', ' / CryptTop / Cycl'),
            (' / ADAMDEC1 Fib', ' / ADAMDEC1'), (' / Neuron / Pericyte', ' / Neuron / Peri'), (' / Mac S_M_P / CD9 Mac', ' / Mac / CD9 Mac'),
            ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'), ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'),
            ('B / Tfh / Resting T', 'B / Tfh / Rest T')]
    def short_cn(x):
        for a, b in TRIM: x = x.replace(a, b)
        return x
    LAB = {c: short_cn(NM[c]) for c in CNS}
    H82 = pd.read_csv(f'{OUTA}/supp_lr_modules_allcells_heatmap.csv')
    assert [LAB[c] for c in CNS] == H82.drop_duplicates('CN').set_index('CN').loc[CNS, 'CN_label'].tolist()
    ST = []
    for a in AXN:
        d = V[V.axis == a].dropna(subset=['z'])
        for c in CNS:
            dd = d.assign(y=d.z, in_cn=(d.CN == c).astype(float)); zc = dd[dd.CN == c].y.values
            f, meth = fit('y ~ in_cn', dd, 'in_cn') if len(zc) >= 3 else (None, 'not fitted')
            b, se, p = (f.params['in_cn'], f.bse['in_cn'], f.pvalues['in_cn']) if f is not None else (np.nan,) * 3
            t = sst.t.ppf(0.975, len(zc) - 1) if len(zc) > 1 else np.nan
            ST.append(dict(axis=a, lr_pairs='; '.join(p_.replace('^', ' -> ') for p_ in AXES[a]), CN=c, CN_label=LAB[c], CN_name=NM[c], n_sections=len(zc),
                           mean_within_section_z=zc.mean(), mean_ci_low=zc.mean() - t * zc.std(ddof=1) / np.sqrt(len(zc)),
                           mean_ci_high=zc.mean() + t * zc.std(ddof=1) / np.sqrt(len(zc)), sections_above_0=int((zc > 0).sum()),
                           wilcoxon_p=sst.wilcoxon(zc).pvalue if len(zc) >= 3 else np.nan, lmm_estimate=b, lmm_ci_low=b - 1.96 * se, lmm_ci_high=b + 1.96 * se,
                           lmm_p=p, optimiser=meth, mean_raw_score=V[(V.axis == a) & (V.CN == c)].raw.mean()))
    ST = pd.DataFrame(ST); ST['lmm_q'] = np.nan; ST['wilcoxon_q'] = np.nan
    for a, g in ST.groupby('axis'):
        for k_ in ['lmm', 'wilcoxon']:
            ok = g[f'{k_}_p'].notna(); ST.loc[g.index[ok], f'{k_}_q'] = multipletests(g[f'{k_}_p'][ok], method='fdr_bh')[1]
    # reference: Inflammatory Cytokine module value of cn_lr_modules_allcells per CN
    ST['inflammatory_cytokine_module_z_82'] = ST.CN.map(H82[H82.module == 'Inflammatory Cytokine'].set_index('CN').mean_within_section_z)
    print('optimisers:', ST.optimiser.value_counts().to_dict())
    E_ = ST.pivot(index='CN', columns='axis', values='mean_within_section_z').loc[CNS, AXN]; Q_ = ST.pivot(index='CN', columns='axis', values='lmm_q').loc[CNS, AXN]
    assert E_.notna().all().all() and Q_.notna().all().all()
    print('==== axis matrix (Mean Within-Section z, LMM stars | Wilcoxon stars)')
    W_ = ST.pivot(index='CN', columns='axis', values='wilcoxon_q').loc[CNS, AXN]
    print(E_.round(2).astype(str).add(Q_.applymap(stars)).add(' | ').add(W_.applymap(stars)).set_axis([LAB[c] for c in CNS]).to_string())
    print(ST[ST.CN == TARGET][['axis', 'n_sections', 'mean_within_section_z', 'mean_ci_low', 'mean_ci_high', 'sections_above_0', 'lmm_estimate', 'lmm_ci_low', 'lmm_ci_high', 'lmm_q', 'wilcoxon_q']].round(3).to_string(index=False))
    # ---------------- heatmap (TNF set apart by a gap)
    VM = float(np.ceil(np.abs(E_.values).max() * 10) / 10)
    ann = np.array([[f'{E_.iloc[r, c]:+.2f}{stars(Q_.iloc[r, c])}' for c in range(len(AXN))] for r in range(len(CNS))])
    plt.rcParams.update({'font.size': 17, 'font.family': 'DejaVu Sans'})
    CW, CH, GAP, X0, Y0 = 1.06, 0.39, 0.13, 3.75, 0.62; FW, FH = 9.4, 4.75; HH = 10 * CH
    fig = plt.figure(figsize=(FW, FH))
    ax1 = fig.add_axes([X0 / FW, Y0 / FH, 4 * CW / FW, HH / FH]); ax2 = fig.add_axes([(X0 + 4 * CW + GAP) / FW, Y0 / FH, CW / FW, HH / FH])
    cax = fig.add_axes([(X0 + 5 * CW + GAP + 0.13) / FW, (Y0 + 0.075 * HH) / FH, 0.26 / FW, 0.85 * HH / FH])
    kw = dict(cmap='RdBu_r', center=0, vmin=-VM, vmax=VM, fmt='', annot_kws={'fontsize': 15}, linewidths=0.5, linecolor='white')
    sns.heatmap(E_[AXN[:4]].set_axis([LAB[c] for c in CNS]), annot=ann[:, :4], ax=ax1, cbar=False, **kw)
    sns.heatmap(E_[[CONTRAST]].set_axis([LAB[c] for c in CNS]), annot=ann[:, 4:], ax=ax2, cbar_ax=cax, cbar_kws={'label': 'Mean Within-Section z'}, **kw)
    ax1.set_xticklabels(AXN[:4], rotation=0, fontsize=14); ax2.set_xticklabels([CONTRAST], rotation=0, fontsize=14)
    ax1.set_yticklabels([LAB[c] for c in CNS], rotation=0, fontsize=13); ax2.set_yticks([])
    ax1.set_xlabel(''); ax2.set_xlabel(''); ax2.set_ylabel(''); ax1.set_ylabel('Spatial Neighborhood', fontsize=15)
    cax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−')))
    fig.savefig(f'{HSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
    OUT1 = ST.assign(cell_text=[ann[CNS.index(c), AXN.index(a)] for c, a in zip(ST.CN, ST.axis)],
                     metric='LIANA bivariate local score (raw), mean over the axis pairs per cell, mean over ALL cells of the CN per section; z across CNs within section, mean over sections',
                     model='z ~ in_CN + (1 | section) on per-section z; BH over 10 CNs within axis')
    OUT1 = pd.concat([OUT1.assign(row_type='summary'), V.rename(columns={'z': 'within_section_z', 'raw': 'raw_score'}).assign(row_type='section')], ignore_index=True)
    OUT1.apply(sig3).to_csv(f'{HSTEM}.csv', index=False)

    # ===================================================================== PARTS 2 and 3: expression per section
    if not SKIP:
        sa = pq.read_table(f'{R}/combined/all_samples.parquet', columns=['CenterX_global_px', 'CenterY_global_px', 'cell_type_general', 'patient_sec', '__index_level_0__']).to_pandas()
        sa['section'] = sa.patient_sec.astype(str); sa.index = [p + '_' + str(i) for p, i in zip(sa.section, sa.index)]; assert not sa.index.duplicated().any()
        sa['x'], sa['y'] = sa.CenterX_global_px * UM, sa.CenterY_global_px * UM; sa['ct'] = sa.cell_type_general.astype(str)
        cnt = pq.read_table(f'{CND}/cn_per_cell.parquet', columns=['CN_k10', '__index_level_0__']).to_pandas()
        assert cnt.index.isin(sa.index).all(); sa['CN'] = cnt.CN_k10.reindex(sa.index).values
        assert (sa[sa.CN.isna()].ct == 'SELENOP+ Fib').all() and not (sa[sa.CN.notna()].ct == 'SELENOP+ Fib').any()
        PC = pd.read_csv(f'{P4}/D_shared_inputs/treg_definition_comparison_per_cell.csv', index_col=0); pcr = PC[PC.cohort == 'resection']; assert len(pcr) == 6554
        T = sa.loc[pcr.index, ['section', 'x', 'y', 'CN']].copy()
        T['treg_suppression_score_z'] = -pcr.signature_score_D6.values; T['treg_suppression_score_z'] = T.groupby('section').treg_suppression_score_z.transform(zs)
        for g in LIG_T: T['lig_' + g] = np.nan
        for g in REC_T + REC_X: T['expr_' + g] = np.nan
        T['nb_n'] = np.nan; T['h5ad_log1p_total_counts'] = np.nan
        GL = list(dict.fromkeys(LIG_SRC + REC_T + REC_X)); SRC = []
        for s_, f in FOLDER.items():
            sub = sa[sa.section == s_]
            with h5py.File(f'{R}/{f}/Processed/norm_log1p.h5ad', 'r') as h:
                ob = read_elem(h['obs']); v = pd.Index(read_elem(h['var']).index.astype(str)); X = read_elem(h['X'])
            nm = pd.Index(s_ + '_' + ob.index.astype(str)); X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
            assert all(g in v for g in GL), [g for g in GL if g not in v]
            keep = np.where(nm.isin(sub.index))[0]
            E = pd.DataFrame(X[keep][:, [v.get_loc(g) for g in GL]].toarray(), index=nm[keep], columns=GL)
            lib = pd.Series(ob['log1p_total_counts'].values[keep], index=nm[keep]); del X; gc.collect()
            E = E.loc[E.index.intersection(sub.index)]; cs = sub.loc[E.index]
            # ---- Part 2: CN-assigned cells
            inc = cs.CN.notna().values; assert inc.sum() == sub.CN.notna().sum(), (s_, inc.sum(), sub.CN.notna().sum())
            ec, cc_ = E[inc], cs[inc]
            df = pd.concat([(ec[LIG_SRC + REC_CTX] > 0).astype(int).add_prefix('npos_'), ec[LIG_SRC + REC_CTX].astype(np.float64).add_prefix('sum_')], axis=1)
            df['n_cells'] = 1; df['ct'] = cc_.ct.values; df['group'] = np.where(cc_.CN.values == TARGET, 'CN10', 'All other CNs')
            SRC.append(df.groupby(['group', 'ct']).sum().reset_index().assign(section=s_))
            # ---- Part 3: local ligand availability around each Treg (as the per-Treg local table of 11_neighbourhoods_local_tables.py)
            Z = (E[LIG_T] - E[LIG_T].mean()) / E[LIG_T].std().replace(0, np.nan)
            tg = T[T.section == s_]
            tree = cKDTree(cs[['x', 'y']].values); nb = tree.query_ball_point(tg[['x', 'y']].values, r=RAD)
            T.loc[tg.index, 'nb_n'] = [len(i) for i in nb]
            for g in LIG_T:
                vals = Z[g].values; T.loc[tg.index, 'lig_' + g] = [vals[i].mean() if len(i) >= MIN_NB else np.nan for i in nb]
            common = tg.index.intersection(E.index)
            for g in REC_T + REC_X: T.loc[common, 'expr_' + g] = E.loc[common, g].values
            T.loc[common, 'h5ad_log1p_total_counts'] = lib.loc[common].values
            print(s_, 'cells', len(cs), '| CN cells', int(inc.sum()), '| Tregs', len(tg), '| detection (% cells > 0):',
                  {g: round(100 * (E[g] > 0).mean(), 2) for g in LIG_SRC}, flush=True)
            del E, Z, ec, cs, df; gc.collect()
        pd.concat(SRC, ignore_index=True).to_csv(C2, index=False, float_format='%.8g')
        rc = pd.read_csv(f'{R}/combined/raw_counts_strict_tregs.csv', index_col=0)
        assert set(rc.index) == set(T.index) and not rc.index.duplicated().any()
        rc = rc.loc[T.index]; T['total_counts'] = rc.sum(1).values
        for g in REC_T + REC_X: T['raw_' + g] = rc[g].values
        del rc; gc.collect()
        T.to_csv(C3, float_format='%.6g')
    T = pd.read_csv(C3, index_col=0)
    # ---- check against treg_local_ligand_receptor.csv.gz of 11_neighbourhoods_local_tables.py (Treg Suppression Score, local IL6 and TNF, neighbour counts)
    T64 = pd.read_csv(f'{LRS}/treg_local_ligand_receptor.csv.gz', index_col=0).loc[T.index]
    for a_, b_ in [('treg_suppression_score_z', 'act_z'), ('lig_IL6', 'lig_IL6'), ('lig_TNF', 'lig_TNF'), ('nb_n', 'nb_n')]:
        assert np.allclose(T[a_].values, T64[b_].values, rtol=1e-3, atol=2e-4, equal_nan=True), (a_, float(np.nanmax(np.abs(T[a_].values - T64[b_].values))))
    print('treg_local_ligand_receptor.csv.gz reproduced: Treg Suppression Score, lig_IL6, lig_TNF, nb_n', flush=True)
    for g in REC_T + REC_X:
        agree = ((T['raw_' + g] > 0) == (T['expr_' + g] > 0)).mean(); assert agree > 0.999, (g, agree)
    T['log_total_counts'] = np.log1p(T.total_counts); T['libsize_z'] = T.groupby('section').log_total_counts.transform(zs)
    print('raw-count library vs h5ad log1p_total_counts, Pearson r =', round(float(np.corrcoef(T.log_total_counts, T.h5ad_log1p_total_counts)[0, 1]), 4))
    T['CN'] = T.CN.astype(int); T['in_cn10'] = (T.CN == TARGET).astype(float)
    NSEC = T.section.nunique(); NTREG = len(T)

    # ===================================================================== PART 2 summary
    S2 = pd.read_csv(C2)
    secs = sorted(S2.section.unique()); cts = sorted(S2.ct.unique()); grps = ['CN10', 'All other CNs']
    full = pd.MultiIndex.from_product([secs, grps, cts], names=['section', 'group', 'ct'])
    S2 = S2.set_index(['section', 'group', 'ct']).reindex(full, fill_value=0).reset_index()
    rows = []
    for g in LIG_SRC + REC_CTX:
        s = S2[['section', 'group', 'ct', 'n_cells', f'npos_{g}', f'sum_{g}']].rename(columns={f'npos_{g}': 'npos', f'sum_{g}': 'sumx'})
        s['grp_n'] = s.groupby(['section', 'group']).n_cells.transform('sum'); s['grp_pos'] = s.groupby(['section', 'group']).npos.transform('sum')
        s['grp_sum'] = s.groupby(['section', 'group']).sumx.transform('sum')
        s['pct_expr'] = np.where(s.n_cells >= MINC, 100 * s.npos / s.n_cells.replace(0, np.nan), np.nan)
        s['share_pos'] = np.where(s.grp_pos > 0, 100 * s.npos / s.grp_pos.replace(0, np.nan), np.nan)
        s['share_sum'] = np.where(s.grp_sum > 0, 100 * s.sumx / s.grp_sum.replace(0, np.nan), np.nan)
        s['pct_of_cells'] = np.where(s.grp_n >= MIN_CELLS_CN, 100 * s.n_cells / s.grp_n.replace(0, np.nan), np.nan)
        # all-cells reference
        a = s.groupby(['section', 'group']).agg(n=('n_cells', 'sum'), p=('npos', 'sum')).reset_index(); a['pct'] = np.where(a.n >= MIN_CELLS_CN, 100 * a.p / a.n, np.nan)
        w = a.pivot(index='section', columns='group', values='pct').dropna()
        pa = a.groupby('group')[['n', 'p']].sum()
        rows.append(dict(gene=g, gene_role='ligand' if g in LIG_SRC else 'receptor (context)', cell_type_general='ALL CELLS', cell_type_label='All Cells', n_sections_pct=len(w), cn10_pct_expressing=w['CN10'].mean(),
                         other_cns_pct_expressing=w['All other CNs'].mean(), sections_cn10_higher=int((w['CN10'] > w['All other CNs']).sum()),
                         cn10_pooled_n_cells=int(pa.loc['CN10', 'n']), cn10_pooled_n_expressing=int(pa.loc['CN10', 'p']),
                         other_cns_pooled_pct_expressing=100 * pa.loc['All other CNs', 'p'] / pa.loc['All other CNs', 'n'],
                         cn10_pooled_pct_expressing=100 * pa.loc['CN10', 'p'] / pa.loc['CN10', 'n']))
        for t in cts:
            x = s[s.ct == t]; xc = x[x.group == 'CN10'].set_index('section'); xo = x[x.group == 'All other CNs'].set_index('section')
            both = xc.pct_expr.notna() & xo.pct_expr.reindex(xc.index).notna()
            rows.append(dict(gene=g, gene_role='ligand' if g in LIG_SRC else 'receptor (context)', cell_type_general=t, cell_type_label=DISP.get(t, t), n_sections_pct=int(both.sum()),
                             cn10_pct_expressing=xc.pct_expr[both].mean(), other_cns_pct_expressing=xo.pct_expr.reindex(xc.index)[both].mean(),
                             sections_cn10_higher=int((xc.pct_expr[both] > xo.pct_expr.reindex(xc.index)[both]).sum()),
                             cn10_share_of_positive_cells=xc.share_pos.mean(), other_cns_share_of_positive_cells=xo.share_pos.mean(),
                             n_sections_share_cn10=int(xc.share_pos.notna().sum()), cn10_share_of_expression=xc.share_sum.mean(),
                             other_cns_share_of_expression=xo.share_sum.mean(), cn10_pct_of_cells=xc.pct_of_cells.mean(),
                             other_cns_pct_of_cells=xo.pct_of_cells.mean(), cn10_pooled_n_cells=int(xc.n_cells.sum()), cn10_pooled_n_expressing=int(xc.npos.sum()),
                             cn10_pooled_pct_expressing=100 * xc.npos.sum() / xc.n_cells.sum() if xc.n_cells.sum() else np.nan,
                             other_cns_pooled_pct_expressing=100 * xo.npos.sum() / xo.n_cells.sum() if xo.n_cells.sum() else np.nan))
    SRCT = pd.DataFrame(rows)
    assert not SRCT.cell_type_general.isin(['SELENOP+ Fib']).any()
    SRCT['cn10_share_rank'] = SRCT[SRCT.cell_type_general != 'ALL CELLS'].groupby('gene').cn10_share_of_positive_cells.rank(ascending=False, method='min')
    SRCT['cn10_vs_other_pct_ratio'] = SRCT.cn10_pct_expressing / SRCT.other_cns_pct_expressing.replace(0, np.nan)
    SRCT['gene'] = pd.Categorical(SRCT.gene, LIG_SRC + REC_CTX); SRCT = SRCT.sort_values(['gene', 'cn10_share_rank'], na_position='first').reset_index(drop=True)
    SRCT['metric'] = ('% expressing = norm_log1p > 0 among cells of the type in the group, per-section mean over sections with >= 20 cells of the type in CN10 and in all other CNs; '
                  'share = % of the group\'s gene-positive cells (or summed norm_log1p) contributed by the type, per-section mean; receptor rows = receptor side of the axes, context only')
    SRCT.apply(sig3).to_csv(SRCCSV, index=False)
    print('==== PART 2: ligand sources in CN10 (top 5 by share of ligand+ cells; per-section means)')
    for g in LIG_SRC:
        x = SRCT[SRCT.gene == g]; al = x[x.cell_type_general == 'ALL CELLS'].iloc[0]
        print(f'-- {g}: all cells % expressing CN10 {al.cn10_pct_expressing:.2f} vs other {al.other_cns_pct_expressing:.2f} ({al.sections_cn10_higher}/{al.n_sections_pct} sections higher); '
          f'pooled CN10 n+ = {al.cn10_pooled_n_expressing}')
        print(x[x.cell_type_general != 'ALL CELLS'].nsmallest(5, 'cn10_share_rank')[['cell_type_label', 'cn10_share_of_positive_cells', 'cn10_share_of_expression',
              'cn10_pct_of_cells', 'cn10_pct_expressing', 'other_cns_pct_expressing', 'sections_cn10_higher', 'n_sections_pct', 'cn10_pooled_n_expressing']].round(2).to_string(index=False))

    print('==== receptor side (context), all cells: % expressing CN10 vs other CNs (per-section mean)')
    print(SRCT[(SRCT.gene_role != 'ligand') & (SRCT.cell_type_general == 'ALL CELLS')][['gene', 'cn10_pct_expressing', 'other_cns_pct_expressing', 'sections_cn10_higher', 'n_sections_pct']].round(2).to_string(index=False))
    for g in REC_CTX:
        x = SRCT[(SRCT.gene == g) & (SRCT.cell_type_general != 'ALL CELLS')]
        print(g, 'top CN10 holders:', '; '.join(f'{r.cell_type_label} {r.cn10_share_of_positive_cells:.1f}% of + cells ({r.cn10_pct_expressing:.1f}% expr)' for r in x.nsmallest(3, 'cn10_share_rank').itertuples()))
    # ===================================================================== PART 3
    def per_section_slopes(d, sign_ref, min_n):
        out = []
        for s_, g in d.groupby('section'):
            if len(g) < min_n or g.x.std() == 0: continue
            Xm = np.column_stack([np.ones(len(g)), g.x.values, g.libsize_z.values]); bb = np.linalg.lstsq(Xm, g.y.values, rcond=None)[0]; out.append((s_, len(g), bb[1]))
        o = pd.DataFrame(out, columns=['section', 'n', 'slope'])
        return o, int((np.sign(o.slope) == np.sign(sign_ref)).sum()), len(o)
    P3 = []; SEC3 = []
    cn10_sec = T[T.in_cn10 == 1].groupby('section').size(); NSEC_CN10_GE10 = int((cn10_sec >= MIN_TREG_SEC).sum())
    print(f'CN10 Tregs: {int(T.in_cn10.sum())} in {len(cn10_sec)} sections; sections with >= {MIN_TREG_SEC} CN10 Tregs: {NSEC_CN10_GE10}')
    for g in LIG_T:
        d = T.assign(y=T.treg_suppression_score_z, x=T['lig_' + g])[['section', 'y', 'x', 'libsize_z', 'in_cn10']].dropna()
        sdx = d.x.std()
        # primary
        f, meth = fit('y ~ x + libsize_z', d, 'x'); b, se, p = f.params['x'], f.bse['x'], f.pvalues['x']
        o, agr, nt = per_section_slopes(d, b, MIN_TREG_SEC); SEC3.append(o.assign(ligand=g, analysis='primary_all_tregs'))
        P3.append(dict(analysis='primary_all_tregs', ligand=g, role='control' if g == 'TNF' else 'test', n_tregs=len(d), n_sections=d.section.nunique(),
                       estimate=b, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se, p=p, sd_local_ligand=sdx, estimate_per_sd=b * sdx,
                       ci_low_per_sd=(b - 1.96 * se) * sdx, ci_high_per_sd=(b + 1.96 * se) * sdx, sections_agreeing=agr, sections_tested=nt,
                       libsize_estimate=f.params['libsize_z'], optimiser=meth, pct_tregs_ligand_constant_low=100 * (d.x == d.groupby('section').x.transform('min')).mean()))
        # CN10 descriptive
        dc = d[d.in_cn10 == 1]; nsec10 = int((dc.groupby('section').size() >= MIN_TREG_SEC).sum())
        r = dict(analysis='cn10_tregs_descriptive', ligand=g, role='control' if g == 'TNF' else 'test', n_tregs=len(dc), n_sections=dc.section.nunique(),
                 n_sections_ge10_tregs=nsec10, sd_local_ligand=dc.x.std())
        if nsec10 >= MIN_SEC_CN10 and dc.x.std() > 0:
            f, meth = fit('y ~ x + libsize_z', dc, 'x')
            if f is not None:
                b, se, p = f.params['x'], f.bse['x'], f.pvalues['x']; o, agr, nt = per_section_slopes(dc, b, MIN_TREG_SEC); SEC3.append(o.assign(ligand=g, analysis='cn10_tregs_descriptive'))
                r.update(estimate=b, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se, p_nominal=p, estimate_per_sd=b * sdx, ci_low_per_sd=(b - 1.96 * se) * sdx,
                         ci_high_per_sd=(b + 1.96 * se) * sdx, sections_agreeing=agr, sections_tested=nt, optimiser=meth)
        else: r.update(note=f'underpowered: fewer than {MIN_SEC_CN10} sections with >= {MIN_TREG_SEC} CN10 Tregs (or no variation); not fitted')
        P3.append(r)
        # exploratory interaction (all Tregs)
        f, meth = fit('y ~ x * in_cn10 + libsize_z', d, 'x:in_cn10')
        bi, sei, pi = f.params['x:in_cn10'], f.bse['x:in_cn10'], f.pvalues['x:in_cn10']; cv = f.cov_params()
        s10 = f.params['x'] + bi; se10 = np.sqrt(cv.loc['x', 'x'] + cv.loc['x:in_cn10', 'x:in_cn10'] + 2 * cv.loc['x', 'x:in_cn10'])
        P3.append(dict(analysis='exploratory_interaction_all_tregs', ligand=g, role='control' if g == 'TNF' else 'test', n_tregs=len(d), n_sections=d.section.nunique(),
                       estimate=bi, ci_low=bi - 1.96 * sei, ci_high=bi + 1.96 * sei, p_nominal=pi, estimate_per_sd=bi * sdx, ci_low_per_sd=(bi - 1.96 * sei) * sdx,
                       ci_high_per_sd=(bi + 1.96 * sei) * sdx, slope_outside_cn10=f.params['x'], slope_in_cn10=s10, slope_in_cn10_ci_low=s10 - 1.96 * se10,
                       slope_in_cn10_ci_high=s10 + 1.96 * se10, slope_in_cn10_p=2 * sst.norm.sf(abs(s10 / se10)), in_cn10_main_effect=f.params['in_cn10'], optimiser=meth))
        # context: local ligand availability of CN10 Tregs vs other Tregs (per-section mean difference, sections with >= 10 CN10 Tregs)
        w = d.groupby(['section', 'in_cn10']).x.mean().unstack(); w = w[w.index.isin(cn10_sec[cn10_sec >= MIN_TREG_SEC].index)]
        P3.append(dict(analysis='context_local_ligand_cn10_vs_other_tregs', ligand=g, role='control' if g == 'TNF' else 'test', n_sections=len(w),
                       cn10_tregs_mean_local_ligand=w[1.0].mean(), other_tregs_mean_local_ligand=w[0.0].mean(), estimate=(w[1.0] - w[0.0]).mean(),
                       sections_agreeing=int((w[1.0] > w[0.0]).sum()), sections_tested=len(w)))
    P3 = pd.DataFrame(P3)
    pm = P3.analysis == 'primary_all_tregs'; P3.loc[pm, 'q'] = multipletests(P3.loc[pm, 'p'], method='fdr_bh')[1]
    P3['model'] = P3.analysis.map({'primary_all_tregs': 'Treg Suppression Score ~ local ligand + libsize_z + (1 | section); BH over 5 ligands',
                                   'cn10_tregs_descriptive': 'same model, CN10 Tregs only; descriptive (no BH)',
                                   'exploratory_interaction_all_tregs': 'Treg Suppression Score ~ local ligand x in_CN10 + libsize_z + (1 | section); estimate = interaction; exploratory',
                                   'context_local_ligand_cn10_vs_other_tregs': 'mean local ligand availability, CN10 Tregs minus other Tregs, per-section mean (sections >= 10 CN10 Tregs)'})
    P3.apply(sig3).to_csv(f'{OUTB}/treg_suppression_inflammatory_ligands_stats.csv', index=False)
    pd.concat(SEC3, ignore_index=True).apply(sig3).to_csv(f'{OUTB}/treg_suppression_inflammatory_ligands_per_section_slopes.csv', index=False)
    print('==== PART 3 primary (all Tregs)')
    print(P3[pm][['ligand', 'n_tregs', 'n_sections', 'estimate', 'ci_low', 'ci_high', 'p', 'q', 'estimate_per_sd', 'ci_low_per_sd', 'ci_high_per_sd', 'sections_agreeing', 'sections_tested', 'pct_tregs_ligand_constant_low']].round(4).to_string(index=False))
    for an in ['cn10_tregs_descriptive', 'exploratory_interaction_all_tregs', 'context_local_ligand_cn10_vs_other_tregs']:
        print('====', an); x = P3[P3.analysis == an].dropna(axis=1, how='all'); print(x.drop(columns=['model', 'analysis', 'role']).round(4).to_string(index=False))
    # ---- Treg receptor expression, CN10 vs other CNs
    RC = []
    okc = cn10_sec[cn10_sec >= MIN_TREG_SEC].index
    for g in REC_T + REC_X:
        pos = (T['raw_' + g] > 0).astype(float); w = pos.groupby([T.section, T.in_cn10]).mean().unstack() * 100; w = w.loc[w.index.isin(okc)]
        RC.append(dict(receptor=g, set='pre-specified' if g in REC_T else 'additional co-receptor', pct_all_tregs=100 * pos.mean(),
                       pct_cn10_tregs_pooled=100 * pos[T.in_cn10 == 1].mean(), pct_other_tregs_pooled=100 * pos[T.in_cn10 == 0].mean(),
                       n_cn10_tregs=int((T.in_cn10 == 1).sum()), n_other_tregs=int((T.in_cn10 == 0).sum()), n_sections=len(w),
                       pct_cn10_tregs_section_mean=w[1.0].mean(), pct_other_tregs_section_mean=w[0.0].mean(), sections_cn10_higher=int((w[1.0] > w[0.0]).sum()),
                       mean_count_all_tregs=T['raw_' + g].mean()))
    RC = pd.DataFrame(RC)
    RC['note'] = 'expressing = raw count > 0 (identical to norm_log1p > 0); section means over sections with >= 10 CN10 Tregs'
    RC.apply(sig3).to_csv(f'{OUTB}/treg_receptor_expression_cn10_vs_other.csv', index=False)
    print('==== Treg receptor expression (% of Tregs with >0 counts)'); print(RC.drop(columns=['note']).round(2).to_string(index=False))
    # ---- forest plot only if a primary test is significant
    SIGP = P3[pm & (P3.q < 0.05)]
    if len(SIGP):
        x = P3[pm].iloc[::-1].reset_index(drop=True); plt.rcParams.update({'font.size': 14})
        fig, ax = plt.subplots(figsize=(4.6, 2.9)); y = np.arange(len(x))
        ax.axvline(0, color='#888888', lw=1)
        for i, r in x.iterrows():
            c = '#b2182b' if r.q < 0.05 else '#7f7f7f'
            ax.hlines(i, r.ci_low_per_sd, r.ci_high_per_sd, color=c, lw=2.4); ax.plot(r.estimate_per_sd, i, 'o', ms=9, color=c, mfc=c if r.q < 0.05 else 'white', mew=2)
            if stars(r.q): ax.text(r.ci_high_per_sd, i + 0.12, ' ' + stars(r.q), va='center', fontsize=14, color=c)
        ax.set_yticks(y); ax.set_yticklabels([f'Local {g}' + (' (control)' if g == 'TNF' else '') for g in x.ligand], fontsize=14)
        ax.set_xlabel('Treg Suppression Score per SD\nof Local Ligand (95% CI)', fontsize=14); ax.spines[['top', 'right']].set_visible(False)
        fig.savefig(f'{FSTEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.04); plt.close(fig)
        P3[pm].apply(sig3).to_csv(f'{FSTEM}.csv', index=False)
        print('forest plot written (primary q < 0.05 for', SIGP.ligand.tolist(), ')')
    else: print('no primary test significant (all q >= 0.05): no Part 3 figure')

    # ===================================================================== legends and methods
    c10 = ST[ST.CN == TARGET].set_index('axis'); pr = P3[pm].set_index('ligand'); nsec_ax = V.dropna(subset=['z']).groupby('axis').section.nunique()
    methods1 = (f'CD resections (CosMx 6K; {NSEC} sections, 4 patients; {NCELLS:,} cells with a neighbourhood '
            'assignment). Five signalling axes were defined before analysis from the NicheNet v2 prior-network ligand–receptor (LR) pairs present in '
            'the cell-type-agnostic spatial LIANA output used for the ligand-receptor modules (bivariate local scores, Gaussian kernel, 200 px bandwidth): IL-6 = IL6→IL6R, '
            'IL6→IL6ST; IL-1 = IL1A→IL1R1, IL1A→IL1RAP, IL1B→IL1R1, IL1B→IL1RAP; OSM = OSM→OSMR, OSM→LIFR, OSM→IL6ST; IL-11 = IL11→IL11RA, IL11→IL6ST; and, '
            'as a contrast, TNF = TNF→TNFRSF1A, TNF→TNFRSF1B (all members of the Inflammatory Cytokine module). Each axis was computed exactly like a '
            'ligand-receptor module: a cell\'s axis score is the mean raw local score over the axis pairs, and a neighbourhood\'s value in a section is the mean over '
            'all its cells, so the score reflects both ligand–receptor expression and the abundance of interacting cells (neighbourhood × section '
            'combinations with at least 100 cells). Values were z-standardised across the 10 neighbourhoods within each section and averaged over '
            'sections ("Mean Within-Section z"; 0 = neighbourhood average); sections in which an axis showed no variation were left out for that axis '
            f'(IL-6 {nsec_ax["IL-6"]}, IL-1 {nsec_ax["IL-1"]}, OSM {nsec_ax["OSM"]}, IL-11 {nsec_ax["IL-11"]} and TNF {nsec_ax["TNF"]} sections). '
            'Significance: linear mixed model on the per-section values, z ~ in neighbourhood + (1 | section), BH-corrected over the 10 neighbourhoods '
            'within each axis.')
    methods2 = ('Ligand sources. For IL6, IL1A, IL1B, OSM, IL11 and TNF, the percentage of cells expressing the ligand (normalised log1p expression > 0) '
            'was computed per cell type (cell_type_general) within CN10 and within all other neighbourhoods of each section and averaged over sections '
            f'with at least {MINC} cells of the type in both; the share of a group\'s ligand-positive cells (and of its summed ligand expression) '
            'contributed by each cell type was averaged over sections in the same way; the receptor genes of the axes (IL6R, IL6ST, IL1R1, IL1RAP, '
            'OSMR, LIFR, IL11RA, TNFRSF1A, TNFRSF1B) were summarised identically for context. SELENOP⁺ Fib cells carry no neighbourhood assignment and '
            'are not shown; SELENOP⁺ Mac are included.')
    methods3 = (f'Link to the Treg Suppression Score (pre-specified; {NTREG:,} strict Tregs, {NSEC} sections). For every Treg, local availability of IL6, '
            'IL1A, OSM, IL11 and TNF (TNF as a control) was the mean of the ligand\'s within-section z-standardised expression over all cells within '
            f'{RAD} µm (at least {MIN_NB} cells; as in the CXCL12 analysis). The Treg Suppression Score (scRNA-seq-derived late-versus-early pseudotime '
            'signature, sign reversed, z-scored within section) was modelled as score ~ local ligand + library size + (1 | section), where library size '
            'is the within-section z of the Treg\'s log total counts; P values were BH-corrected over the five ligands. Per-section ordinary '
            f'least-squares slopes (sections with at least {MIN_TREG_SEC} Tregs) give the sign agreement. Because CN10 contains few Tregs '
            f'({int(T.in_cn10.sum())} Tregs; {NSEC_CN10_GE10} sections with at least {MIN_TREG_SEC}), the same model restricted to CN10 Tregs is '
            'reported descriptively (estimates and 95% CIs, no multiplicity-adjusted claims), and an interaction model over all Tregs, score ~ local '
            'ligand × in CN10 + library size + (1 | section), is reported as exploratory. The Tregs\' own expression (raw count > 0) of IL6R, IL6ST, '
            'IL1R1, OSMR, IL11RA and TNFRSF1B (and the co-receptors IL1RAP, LIFR and TNFRSF1A) is given for CN10 and other neighbourhoods.')
    LIT = ('These axes have been linked to anti-TNF non-response (West et al. 2017, PMID 28368383; Friedrich et al. 2021, PMID 34675383; Smillie et al. '
       '2019, PMID 31348891; Atreya et al. 2000, PMID 10802717).')
    srcsent = ('Ligand sources in CN10 (percentage of cells expressing each ligand by cell type, CN10 versus all other neighbourhoods, and each cell type\'s '
           'share of the ligand-positive cells; per-section means; receptor genes of the axes included for context) are given in '
           'supp_inflammatory_ligand_sources_cn10.csv.')
    legH = ('Inflammatory signalling axes across cellular neighbourhoods (rows) by axis (columns), with the TNF axis set apart as a contrast. ' + methods1
            + ' * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. ' + f'Colour scale symmetric, ±{VM:.1f}. ' + LIT + ' ' + srcsent
            + ' The source data hold the per-section values, the mixed-model estimates and q values, and a Wilcoxon signed-rank test (BH within axis). '
        'Neighbourhood names are abbreviated by their dominant cell types; full names are in the source data.\n\nMethods. Inflammatory signalling axes by neighbourhood. ' + methods1 + ' ' + methods2
            + ' ' + methods3 + '\n')
    for t in (legH, methods1, methods2, methods3):
        assert 'activation score' not in t.lower() and 'quiescent' not in t.lower() and 'territor' not in t.lower() and 'to tregs' not in t.lower()
    open(f'{HSTEM}_legend.txt', 'w').write(legH)
    if len(SIGP):
        legF = ('Local inflammatory ligand availability and the Treg Suppression Score, CD resections. ' + methods3 + ' Points, mixed-model estimate per SD '
            'of local ligand availability; bars, 95% CI; filled and coloured, q < 0.05.')
        assert 'activation score' not in legF.lower() and 'territor' not in legF.lower()
        open(f'{FSTEM}_legend.txt', 'w').write(legF + '\n')
    fmt = lambda r: f'{r.estimate:+.3f} (95% CI {r.ci_low:+.3f} to {r.ci_high:+.3f}), q = {r.q:.3g}, {int(r.sections_agreeing)}/{int(r.sections_tested)} sections same sign'
    md = ('# CN10 as a TNF-independent inflammatory stroma, and the Treg Suppression Score — methods notes\n\n'
      'Script: `cosmx_cd/13_neighbourhood_signalling.py`, section `cn_inflammatory_axes`\n\n'
      'Outputs:\n'
      '- `A_fig4_panels/supp_inflammatory_axes_by_neighbourhood.png` (+ `.csv`, `_legend.txt`)\n'
      '- `A_fig4_panels/supp_inflammatory_ligand_sources_cn10.csv` (Part 2, CSV only)\n'
          + (f'- `A_fig4_panels/supp_inflammatory_ligands_treg_suppression.png` (+ `.csv`, `_legend.txt`)\n' if len(SIGP) else '')
          + '- `B_resection_analyses/inflammatory_axes/` (per section x neighbourhood axis values; per section x cell type ligand counts; per-Treg local ligand '
        'table; Part 3 statistics `treg_suppression_inflammatory_ligands_stats.csv`, per-section slopes, `treg_receptor_expression_cn10_vs_other.csv`)\n\n'
      '## Methods-ready text\n\nInflammatory signalling axes by neighbourhood. ' + methods1 + '\n\n' + methods2 + '\n\n' + methods3 + '\n\n' + LIT + '\n\n'
      '## Validation\n\n'
      '- The per-pair neighbourhood × section cache (section cn10_lr_pairs) reproduces the four module values of section cn_lr_modules_allcells (asserted).\n'
      '- The axis values computed directly from the per-cell LIANA scores equal the mean of the per-pair cache values (linearity; asserted).\n'
      '- The Treg Suppression Score, local IL6 and TNF availability and neighbour counts reproduce treg_local_ligand_receptor.csv.gz of 11_neighbourhoods_local_tables.py (asserted).\n\n'
      '## Key results\n\n'
          + '\n'.join(f'- CN10, {a}: Mean Within-Section z {c10.loc[a, "mean_within_section_z"]:+.2f} (LMM {c10.loc[a, "lmm_estimate"]:+.2f}, 95% CI '
                  f'{c10.loc[a, "lmm_ci_low"]:+.2f} to {c10.loc[a, "lmm_ci_high"]:+.2f}, q = {c10.loc[a, "lmm_q"]:.3g}; {int(c10.loc[a, "sections_above_0"])}/'
                  f'{int(c10.loc[a, "n_sections"])} sections above 0)' for a in AXN) + '\n'
          + '\n'.join(f'- Treg Suppression Score ~ local {g} (all Tregs): {fmt(pr.loc[g])}' for g in LIG_T) + '\n'
          + '\n'.join(f'- CN10 Tregs only (descriptive), local {r.ligand}: {r.estimate:+.3f} (95% CI {r.ci_low:+.3f} to {r.ci_high:+.3f}; {int(r.n_tregs)} Tregs, '
                  f'{int(r.n_sections_ge10_tregs)} sections with >= 10)' for r in P3[P3.analysis == 'cn10_tregs_descriptive'].itertuples() if np.isfinite(r.estimate)) + '\n'
          + '\n'.join(f'- Exploratory interaction local {r.ligand} x in CN10: {r.estimate:+.3f} (95% CI {r.ci_low:+.3f} to {r.ci_high:+.3f}, nominal p = {r.p_nominal:.3g})'
                       for r in P3[P3.analysis == 'exploratory_interaction_all_tregs'].itertuples()) + '\n'
          + '\n'.join(f'- {r.gene} ({r.gene_role}), all cells: {r.cn10_pct_expressing:.2f}% expressing in CN10 vs {r.other_cns_pct_expressing:.2f}% in other neighbourhoods '
                   f'({int(r.sections_cn10_higher)}/{int(r.n_sections_pct)} sections higher in CN10)' for r in SRCT[SRCT.cell_type_general == 'ALL CELLS'].itertuples()) + '\n')
    Pp = S2.groupby('ct').sum(numeric_only=True)
    rng_ = {g: (100 * Pp[f'npos_{g}'] / Pp.n_cells) for g in LIG_SRC}
    sigcn = lambda a: ', '.join(f'CN{c} ({E_.loc[c, a]:+.2f})' for c in CNS if Q_.loc[c, a] < 0.05 and E_.loc[c, a] > 0) or 'none'
    notes = ('\n## Notes for interpretation\n\n'
             + '- Ligand detection is low and broad: % of cells expressing (pooled over sections, by cell type, all neighbourhoods) ranges '
             + '; '.join(f'{g} {rng_[g].min():.1f}-{rng_[g].max():.1f}% (highest: {DISP.get(rng_[g].idxmax(), rng_[g].idxmax())})' for g in LIG_SRC)
             + '. Most cell types express each ligand at a similar low rate, so ligand-positive cells in CN10 largely follow cell abundance.\n'
             + '- The CN10 excess is larger on the receptor side (IL6ST, IL1R1, OSMR, LIFR higher in CN10 in 9-11/11 sections) than on the ligand side '
           '(only IL6 is consistently higher, 10/11 sections); the TNF receptors are also higher in CN10, so the low TNF axis reflects TNF ligand '
           'not being enriched in CN10 (it is highest around lymphoid / CD4 effector cells).\n'
             + '- Neighbourhoods significantly above average (q < 0.05) per axis: ' + '; '.join(f'{a}: {sigcn(a)}' for a in AXN)
             + '. CN10 shares the IL-6 / IL-1 / OSM elevation with CN8 (B / Tfh / Rest T), which in addition is the TNF-high neighbourhood.\n'
             + '- No pre-specified local ligand is associated with the Treg Suppression Score (all q >= ' + f'{P3.loc[pm, "q"].min():.2f}' + '); the CN10-only '
           'estimates are imprecise and centred near 0; the only nominal exploratory interaction (IL1A x CN10) is positive, i.e. opposite to '
           'IL-1 explaining the lower score of CN10 Tregs.\n')
    md = md + notes
    open(METHODS, 'w').write(md)
    print('done', flush=True)


# ======================================================================================================================
# Section S9A_cn_characterisation: CN characterisation figure: cell-type composition, LR modules and inflammatory axes per CN
# Panels: Supp 9A
# ======================================================================================================================
if _run('S9A_cn_characterisation'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/89_cn_niche_characterisation.py'   # working folder of this step (relative paths resolve from here)
    """Cellular neighbourhood characterisation figure. CD resections (CosMx 6K; 11 sections, 4 patients). Top: Cell-Type Composition (log2 enrichment table of
12_neighbourhood_treg_score.py). Bottom row, side by side and sharing the CN labels / CN colour squares on the left: 'Ligand–Receptor Modules'
(cn_lr_modules_costim: Antigen Presentation, Inflammatory Cytokine, the 7-pair APC -> T-cell Co-Stimulatory module, Matrix–Integrin) and 'Inflammatory
Signalling' (cn_inflammatory_axes: IL-6, IL-1, OSM, IL-11 | TNF). Both bottom heatmaps are Mean Within-Section z and share one diverging colourbar whose
symmetric range covers both (max |z| = CN8 Co-Stimulatory).
Sources (read only):
  A_fig4_panels/fig4_cn_composition_heatmap.csv             12_neighbourhood_treg_score.py: log2 enrichment vs all cells (unclipped; shown clipped to +/-3, own colourbar)
  A_fig4_panels/fig4_lr_modules_allcells_costim_dotplot.csv cn_lr_modules_costim: CN summary rows, mean_within_section_z, lmm_q (stars; BH over 10 CNs within module)
  A_fig4_panels/supp_inflammatory_axes_by_neighbourhood.csv cn_inflammatory_axes: summary rows, mean_within_section_z, lmm_q, cell_text (BH over 10 CNs within axis)
Cell text = mean (2 dp) + stars. The module CSV holds 3 significant figures, which is ambiguous at 2 dp for some values (e.g. 0.615), so the 2-dp text is
taken from the exact source: Antigen Presentation / Inflammatory Cytokine / Matrix–Integrin = cell_text of supp_lr_modules_allcells_clean_heatmap.csv
(cn_lr_modules_immunomodulatory; values asserted identical); Co-Stimulatory = the unrounded mean re-derived from the per-pair cache
lr_pairs_v2_per_section_cn.csv (mean over the 7 pairs per CN x section, float32, >= 100 cells, z across CNs within section, mean over sections), asserted
to reproduce the cn_lr_modules_costim CSV (summary means, n_sections and section z, 3 s.f.) and used for the 2-dp text only. Colours and the CSV use the
source CSV values; stars always from the source lmm_q.
Layout: same row height everywhere; composition cells square; bottom cells widened for the 13-pt text (LR columns further, for the 2-line 14-pt module
names); vertical gap set from the rendered geometry (rotated composition labels as exact parallelograms) with a constant clearance.
Outputs: /share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_niche_characterisation.png + .csv (all plotted values) + _legend.txt."""
    import os, sys, re
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX', 'run with OPENBLAS_CORETYPE=SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/'), sys.executable
    import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker, seaborn as sns
    from matplotlib.patches import Rectangle
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/A_fig4_panels'
    CND = f'{P4}/B_resection_analyses/cellular_neighborhoods_v2'; LRC = f'{P4}/B_resection_analyses/lr_modules_cellchat'
    C45 = f'{OUT}/fig4_cn_composition_heatmap.csv'; C86 = f'{OUT}/supp_inflammatory_axes_by_neighbourhood.csv'; L86 = f'{OUT}/supp_inflammatory_axes_by_neighbourhood_legend.txt'
    C82D = f'{OUT}/fig4_lr_modules_allcells_costim_dotplot.csv'; L82D = f'{OUT}/fig4_lr_modules_allcells_costim_dotplot_legend.txt'
    H82C = f'{OUT}/supp_lr_modules_allcells_clean_heatmap.csv'; PAIR84 = f'{LRC}/lr_pairs_v2_per_section_cn.csv'
    STEM = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_niche_characterisation'
    assert not {f'{STEM}.png', f'{STEM}.csv', f'{STEM}_legend.txt'} & {C45, C86, L86, C82D, L82D, H82C}
    CNS = list(range(1, 11))
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    MINUS = lambda t: t.replace('-', '−')

    # ============================================================== data: composition (12_neighbourhood_treg_score.py)
    E = pd.read_csv(C45, index_col=0)
    LIN = {'Epithelial': ['Epithelial'],                                   # composition table (after merging Mo-Mac, CD8 Trm, Tfh)
           'Stromal': ['OGN_RSPO3 Fib', 'PLA2G2A Fib', 'ADAMDEC1 Fib', 'Activated Fib', 'Crypt Top Fib', 'Myofibroblast', 'Pericyte', 'Endothelial', 'Neuron'],
           'Myeloid': ['DC', 'CD9 Mac', 'Mac S', 'Mac S_M_P', 'Mac S_M_S', 'Mo-Mac', 'Neutrophil', 'Granulocyte'],
           'Lymphoid': ['Treg', 'CD4 Trm', 'CD4 Effector', 'Tfh', 'Resting T', 'CD8 Trm', 'Gd T', 'NK', 'B', 'Plasma'],
           'Other': ['Cycling']}
    cols = [c for v in LIN.values() for c in v if c in E.columns]
    assert cols == list(E.columns), 'composition CSV columns differ from the LIN order'
    H = E[cols].clip(-3, 3)                                                  # display clip
    DISP = {'OGN_RSPO3 Fib': 'OGN⁺RSPO3⁺ Fib', 'PLA2G2A Fib': 'PLA2G2A⁺ Fib',          # display names
            'Mac S': 'Mac S+SG+/XS-', 'Mac S_M_P': 'Mac S+M+P+', 'Mac S_M_S': 'Mac S+M+S+',      # cell_type_short names
            'Gd T': 'γδ T'}                                                             # display map
    XL = [DISP.get(c, c) for c in cols]
    _cts = set(pq.read_table('/share/fsmresfiles/UC/AtoMx/CD_surgical_resection/combined/all_samples.parquet', columns=['cell_type_short']).column(0).to_pandas().astype(str))
    _cmp = {'Mac S': ['Mac S+SG+', 'Mac S+XS-'], 'Mac S_M_P': ['Mac S+M+P+'], 'Mac S_M_S': ['Mac S+M+S+']}
    for k_, parts in _cmp.items():
        assert all(p_ in _cts for p_ in parts), (k_, parts)
        assert DISP[k_] == (parts[0] if len(parts) == 1 else parts[0] + '/' + parts[1].replace('Mac S+', '')), (k_, DISP[k_])
    # ============================================================== data: inflammatory axes (cn_inflammatory_axes)
    A = pd.read_csv(C86, low_memory=False); S = A[A.row_type == 'summary'].copy(); S['CN'] = S.CN.astype(int)
    AXN = ['IL-6', 'IL-1', 'OSM', 'IL-11', 'TNF']; CONTRAST = 'TNF'
    assert list(dict.fromkeys(S.axis)) == AXN and len(S) == 50 and not S.duplicated(['axis', 'CN']).any()
    Z = S.pivot(index='CN', columns='axis', values='mean_within_section_z').loc[CNS, AXN]
    Q = S.pivot(index='CN', columns='axis', values='lmm_q').loc[CNS, AXN]
    TX = S.pivot(index='CN', columns='axis', values='cell_text').loc[CNS, AXN]
    assert Z.notna().all().all() and Q.notna().all().all() and TX.notna().all().all()
    for c in CNS:
        for a in AXN:
            t = TX.loc[c, a]; num = t.rstrip('*')
            assert abs(float(num) - Z.loc[c, a]) <= 0.0051 and t[len(num):] == stars(Q.loc[c, a]), (c, a, t, Z.loc[c, a], Q.loc[c, a])
    LAB = S.drop_duplicates('CN').set_index('CN').CN_label.loc[CNS].tolist(); NMF = S.drop_duplicates('CN').set_index('CN').CN_name.loc[CNS]
    assert list(E.index) == LAB, ('CN sets / order differ', list(E.index), LAB)
    assert all(l.startswith(f'CN{c}: ') for c, l in zip(CNS, LAB))
    L86T = open(L86).read().strip(); VM86 = float(np.ceil(np.abs(Z.values).max() * 10) / 10); assert f'Colour scale symmetric, ±{VM86:.1f}.' in L86T
    ANN = TX.apply(lambda s: s.map(MINUS)).values
    # ============================================================== data: LR modules (cn_lr_modules_costim)
    MODS = ['Antigen Presentation', 'Inflammatory Cytokine', 'Co-Stimulatory', 'Matrix–Integrin']
    MODL = ['Antigen\nPresentation', 'Inflammatory\nCytokine', 'Co-\nStimulatory', 'Matrix–\nIntegrin']
    LR = pd.read_csv(C82D, low_memory=False); LS = LR[(LR.row_type == 'summary') & LR.module.isin(MODS)].copy(); LS['CN'] = LS.CN.astype(int)
    assert len(LS) == 40 and not LS.duplicated(['module', 'CN']).any() and LS.shown_in_figure.astype(str).eq('True').all()
    assert (LS.drop_duplicates('CN').set_index('CN').CN_label.loc[CNS].tolist() == LAB)
    ZL = LS.pivot(index='CN', columns='module', values='mean_within_section_z').loc[CNS, MODS]
    QL = LS.pivot(index='CN', columns='module', values='lmm_q').loc[CNS, MODS]
    assert ZL.notna().all().all() and QL.notna().all().all()
    STL = LS.pivot(index='CN', columns='module', values='stars').loc[CNS, MODS].fillna('')
    assert (STL == QL.apply(lambda s: s.map(stars))).all().all(), 'CSV stars != stars(lmm_q)'
    NPAIR = LS.drop_duplicates('module').set_index('module').n_pairs.loc[MODS].astype(int).tolist(); assert NPAIR == [9, 34, 7, 227], NPAIR
    P7 = [p.replace(' → ', '^') for p in LS[LS.module == 'Co-Stimulatory'].lr_pairs.iloc[0].split('; ')]
    assert P7 == ['CD80^CD28', 'CD86^CD28', 'ICOSLG^ICOS', 'TNFSF4^TNFRSF4', 'TNFSF18^TNFRSF18', 'TNFSF9^TNFRSF9', 'CD70^CD27'], P7
    # exact 2-dp text, AP / IC / MI: cell_text of the cn_lr_modules_immunomodulatory heatmap (identical values for these modules; asserted)
    HC = pd.read_csv(H82C); TXL = pd.DataFrame(index=CNS, columns=MODS, dtype=object)
    for m in ['Antigen Presentation', 'Inflammatory Cytokine', 'Matrix–Integrin']:
        h = HC[HC.module == m].set_index('CN').loc[CNS]
        assert np.array_equal(h.mean_within_section_z.values, ZL[m].values) and np.array_equal(h.lmm_q.values, QL[m].values), m
        for c in CNS:
            t = h.loc[c, 'cell_text']; num = t.rstrip('*')
            assert abs(float(num) - ZL.loc[c, m]) <= 0.0051 and t[len(num):] == stars(QL.loc[c, m]), (m, c, t)
        TXL[m] = h.cell_text.values
    # exact 2-dp text, Co-Stimulatory: unrounded mean re-derived as in cn_lr_modules_costim, asserted against its CSV
    L84 = pd.read_csv(PAIR84)
    g = L84[L84.lr_pair.isin(P7)].groupby(['section', 'CN']).agg(raw=('raw', 'mean'), n=('lr_pair', 'size'), n_cells=('n_cells', 'first')).reset_index()
    assert (g.n == 7).all(); g['raw'] = g.raw.astype(np.float32); g = g[g.n_cells >= 100].copy()
    g['z'] = g.groupby('section').raw.transform(lambda x: (x - x.mean()) / x.std()); g.loc[~np.isfinite(g.z), 'z'] = np.nan; g = g.dropna(subset=['z'])
    mz = g.groupby('CN').z.agg(['mean', 'size']).loc[CNS]
    assert [float(f'{v:.3g}') for v in mz['mean']] == ZL['Co-Stimulatory'].tolist(), (mz, ZL['Co-Stimulatory'])
    assert mz['size'].tolist() == LS[LS.module == 'Co-Stimulatory'].set_index('CN').n_sections.loc[CNS].astype(int).tolist()
    sec = LR[(LR.row_type == 'section') & (LR.module == 'Co-Stimulatory')].merge(g, on=['section', 'CN'])
    assert len(sec) == len(g) and np.array_equal(sec.within_section_z.values, np.array([float(f'{v:.3g}') for v in sec.z])), 'section z differ from cn_lr_modules_costim'
    TXL['Co-Stimulatory'] = [f'{v:+.2f}{stars(q)}' for v, q in zip(mz['mean'], QL['Co-Stimulatory'])]
    print('Co-Stimulatory exact means (2-dp text):', dict(zip(CNS, np.round(mz['mean'].values, 4))))
    for (c, m), t in {(10, 'Antigen Presentation'): '+0.62*', (8, 'Co-Stimulatory'): '+1.92***', (9, 'Matrix–Integrin'): '+1.74***'}.items():
        assert TXL.loc[c, m] == t, (c, m, TXL.loc[c, m], t)                 # spot checks
    ANNL = TXL.apply(lambda s: s.map(MINUS)).values
    VM = float(np.ceil(max(np.abs(Z.values).max(), np.abs(ZL.values).max()) * 10) / 10)   # one shared symmetric range for both bottom heatmaps
    print(f'bottom colour range: axes alone +/-{VM86}; shared with LR modules +/-{VM} (max |z| LR {np.abs(ZL.values).max()}, axes {np.abs(Z.values).max()})')
    PAL = dict(zip(CNS, sns.color_palette('tab10', 10).as_hex()))           # CN colours

    # ============================================================== layout (inches)
    FS_TICK, FS_COL, FS_GRP, FS_NUM, FS_CB, FS_TITLE, FS_YLAB = 14, 14, 15, 13, 13, 16, 15
    plt.rcParams.update({'font.size': 14, 'font.family': 'DejaVu Sans', 'font.weight': 'normal', 'axes.titleweight': 'normal',
                         'axes.labelweight': 'normal', 'figure.titleweight': 'normal'})
    CH = 0.38; CWC = 0.38; CWI = 0.98; CWL = 1.45; GAPT = 0.13; GAPLR = 0.30; HH = 10 * CH
    NCc = len(cols); WC = NCc * CWC; WL = 4 * CWL; XI = WL + GAPLR; WI = 4 * CWI + GAPT + CWI; WB = XI + WI
    X0 = 5.2; CBW, CBP = 0.22, 0.20
    SQX, SQW = 0.05, 0.13; PADY = (SQX + SQW + 0.07) * 72
    GAPV = 0.14; ROT = 45
    FW = X0 + max(WC, WB) + CBP + CBW + 1.0; FH = 20.0
    YT = FH - 1.2 - HH
    fig = plt.figure(figsize=(FW, FH))
    ain = lambda x, y, w, h: fig.add_axes([x / FW, y / FH, w / FW, h / FH])
    fx = lambda x: x / FW; fy = lambda y: y / FH
    # ---- top: Cell-Type Composition
    axc = ain(X0, YT, WC, HH); caxc = ain(X0 + WC + CBP, YT + 0.075 * HH, CBW, 0.85 * HH)
    sns.heatmap(H, cmap='RdBu_r', center=0, vmin=-3, vmax=3, ax=axc, cbar_ax=caxc, linewidths=0.4, linecolor='white',
                xticklabels=True, yticklabels=True, cbar_kws={'label': 'log2 Enrichment vs All Cells'})
    axc.set_xticks(np.arange(NCc) + 0.5); axc.set_xticklabels(XL, rotation=ROT, ha='right', va='center', rotation_mode='anchor', fontsize=FS_COL)
    x = 0
    for lin, cs in LIN.items():
        n = len([c for c in cs if c in cols])
        if x: axc.axvline(x, color='k', lw=1.5)
        axc.text(x + n / 2, -0.25, lin, ha='center', va='bottom', fontsize=FS_GRP, fontweight='normal'); x += n
    assert x == NCc
    # ---- bottom row: Ligand–Receptor Modules | Inflammatory Signalling (shared colourbar); provisional vertical position, moved below
    YB0 = 0.5
    axl = ain(X0, YB0, WL, HH)
    ax1 = ain(X0 + XI, YB0, 4 * CWI, HH); ax2 = ain(X0 + XI + 4 * CWI + GAPT, YB0, CWI, HH); caxi = ain(X0 + WB + CBP, YB0 + 0.075 * HH, CBW, 0.85 * HH)
    kw = dict(cmap='RdBu_r', center=0, vmin=-VM, vmax=VM, fmt='', annot_kws={'fontsize': FS_NUM}, linewidths=0.5, linecolor='white')
    sns.heatmap(ZL.set_axis(LAB), annot=ANNL, ax=axl, cbar=False, xticklabels=True, yticklabels=True, **kw)
    sns.heatmap(Z[AXN[:4]].set_axis(LAB), annot=ANN[:, :4], ax=ax1, cbar=False, xticklabels=True, yticklabels=False, **kw)
    sns.heatmap(Z[[CONTRAST]].set_axis(LAB), annot=ANN[:, 4:], ax=ax2, cbar_ax=caxi, xticklabels=True, yticklabels=False,
                cbar_kws={'label': 'Mean Within-Section z'}, **kw)
    axl.set_xticks(np.arange(4) + 0.5); axl.set_xticklabels(MODL, rotation=0, fontsize=FS_COL, linespacing=1.0)
    ax1.set_xticks(np.arange(4) + 0.5); ax1.set_xticklabels(AXN[:4], rotation=0, fontsize=FS_COL)
    ax2.set_xticks([0.5]); ax2.set_xticklabels([CONTRAST], rotation=0, fontsize=FS_COL)
    for a_ in (ax1, ax2): a_.set_yticks([]); a_.set_ylabel('')
    caxi.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1.0 if VM > 1.5 else 0.5))   # integer ticks for the wider shared range
    caxi.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−')))
    for ax in (axc, axl):
        ax.set_yticks(np.arange(10) + 0.5); ax.set_yticklabels(LAB, rotation=0, fontsize=FS_TICK)
        ax.tick_params(axis='y', length=0, pad=PADY); ax.set_xlabel(''); ax.set_ylabel('Spatial Neighborhood', fontsize=FS_YLAB)
        for t in ax.get_yticklabels(): t.set_color('k')
    for ax in (axc, axl, ax1, ax2): ax.tick_params(axis='x', length=3, width=0.8, pad=2); ax.set_xlabel('')
    for cax in (caxc, caxi):
        cax.tick_params(labelsize=FS_CB); cax.yaxis.label.set_size(FS_CB)
        for sp in cax.spines.values(): sp.set_visible(False)
    assert all(t.get_fontweight() in ('normal', 400) for a_ in fig.axes for t in a_.texts + a_.get_xticklabels() + a_.get_yticklabels())

    # ============================================================== geometry: place the bottom row below the composition labels
    fig.canvas.draw(); R = fig.canvas.get_renderer(); DPI = fig.dpi
    bb_in = lambda b: (b.x0 / DPI, b.y0 / DPI, b.x1 / DPI, b.y1 / DPI)
    def label_poly(t):
        w, h, d = R.get_text_width_height_descent(t.get_text(), t.get_fontproperties(), ismath=False)
        ax_, ay_ = t.get_transform().transform(t.get_unitless_position()); th = np.deg2rad(t.get_rotation())
        uv = np.array([[-w, -h / 2], [-w, h / 2], [0, h / 2], [0, -h / 2]])
        P = np.c_[ax_ + uv[:, 0] * np.cos(th) - uv[:, 1] * np.sin(th), ay_ + uv[:, 0] * np.sin(th) + uv[:, 1] * np.cos(th)] / DPI
        b = bb_in(t.get_window_extent(R))
        assert np.allclose([P[:, 0].min(), P[:, 1].min(), P[:, 0].max(), P[:, 1].max()], b, atol=0.03), (t.get_text(), P, b)
        return P
    def poly_min_y(P, a, b):
        lo, hi = max(a, P[:, 0].min()), min(b, P[:, 0].max())
        if lo > hi: return None
        ys = []
        for xx in np.linspace(lo, hi, 200):
            for i in range(4):
                p, q = P[i], P[(i + 1) % 4]
                if min(p[0], q[0]) - 1e-9 <= xx <= max(p[0], q[0]) + 1e-9 and abs(q[0] - p[0]) > 1e-9:
                    ys.append(p[1] + (xx - p[0]) * (q[1] - p[1]) / (q[0] - p[0]))
        return min(ys) if ys else None
    def upper_shapes():
        polys = [label_poly(t) for t in axc.get_xticklabels()]
        boxes = [bb_in(axc.get_window_extent(R)), bb_in(caxc.get_tightbbox(R)), bb_in(axc.yaxis.label.get_window_extent(R))] + \
                [bb_in(t.get_window_extent(R)) for t in axc.get_yticklabels()]
        return polys, boxes
    TL = fig.text(fx(X0 + WL / 2), fy(YB0 + HH + 0.10), 'Ligand–Receptor Modules', ha='center', va='bottom', fontsize=FS_TITLE, transform=fig.transFigure)
    TI = fig.text(fx(X0 + XI + WI / 2), fy(YB0 + HH + 0.10), 'Inflammatory Signalling', ha='center', va='bottom', fontsize=FS_TITLE, transform=fig.transFigure)
    def lower_boxes():
        out = [bb_in(a_.get_window_extent(R)) for a_ in (axl, ax1, ax2)] + [bb_in(caxi.get_tightbbox(R)), bb_in(axl.yaxis.label.get_window_extent(R)),
                                                                             bb_in(TL.get_window_extent(R)), bb_in(TI.get_window_extent(R))]
        out += [bb_in(t.get_window_extent(R)) for a_ in (axl, ax1, ax2) for t in a_.get_xticklabels() + a_.get_yticklabels() if t.get_text()]
        return out
    def needed_shift():
        polys, boxes = upper_shapes(); need = -np.inf
        for (x0, y0, x1, y1) in lower_boxes():
            for P in polys:
                m = poly_min_y(P, x0 - 0.03, x1 + 0.03)
                if m is not None: need = max(need, y1 + GAPV - m)
            for (u0, v0, u1, v1) in boxes:
                if u0 <= x1 + 0.03 and x0 - 0.03 <= u1: need = max(need, y1 + GAPV - v0)
        return need
    sh = needed_shift(); YB = YB0 - sh
    for a_, (x_, w_) in ((axl, (X0, WL)), (ax1, (X0 + XI, 4 * CWI)), (ax2, (X0 + XI + 4 * CWI + GAPT, CWI))): a_.set_position([x_ / FW, YB / FH, w_ / FW, HH / FH])
    caxi.set_position([(X0 + WB + CBP) / FW, (YB + 0.075 * HH) / FH, CBW / FW, 0.85 * HH / FH])
    TL.set_position((fx(X0 + WL / 2), fy(YB + HH + 0.10))); TI.set_position((fx(X0 + XI + WI / 2), fy(YB + HH + 0.10)))
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    chk = needed_shift(); assert abs(chk) < 0.01, chk
    gtop = max(bb_in(t.get_window_extent(R))[3] for t in axc.texts)
    TC = fig.text(fx(X0 + WC / 2), fy(gtop + 0.10), 'Cell-Type Composition', ha='center', va='bottom', fontsize=FS_TITLE, transform=fig.transFigure)
    for yb in (YT, YB):
        for i, c in enumerate(CNS):
            fig.add_artist(Rectangle((fx(X0 - SQX - SQW), fy(yb + HH - (i + 1) * CH + 0.08 * CH)), fx(SQW), fy(0.84 * CH), color=PAL[c],
                                     transform=fig.transFigure, clip_on=False, lw=0))
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    SQ = [bb_in(p.get_window_extent(R)) for p in fig.artists if isinstance(p, Rectangle)]; assert len(SQ) == 20
    for k in range(10):
        a_, b_ = SQ[k], SQ[10 + k]; assert np.allclose([a_[0], a_[2], a_[3] - a_[1]], [b_[0], b_[2], b_[3] - b_[1]], atol=1e-6) and a_[2] < X0
        for blk, ax_ in ((a_, axc), (b_, axl)):
            t = bb_in(ax_.get_yticklabels()[k].get_window_extent(R)); assert abs((blk[1] + blk[3]) / 2 - (t[1] + t[3]) / 2) < 0.05 and blk[0] > t[2]
    xr_c = [bb_in(t.get_window_extent(R))[2] for t in axc.get_yticklabels()]; xr_l = [bb_in(t.get_window_extent(R))[2] for t in axl.get_yticklabels()]
    assert np.allclose(xr_c, xr_l, atol=1e-6) and np.allclose(xr_c, xr_c[0], atol=1e-6)
    # bottom row: identical row geometry in the three heatmaps
    yb_ = [bb_in(a_.get_window_extent(R)) for a_ in (axl, ax1, ax2)]; assert all(abs(b[1] - yb_[0][1]) < 1e-6 and abs(b[3] - yb_[0][3]) < 1e-6 for b in yb_)
    TXT = [TC, TL, TI, axc.yaxis.label, axl.yaxis.label, caxc.yaxis.label, caxi.yaxis.label] + list(axc.texts) + list(axc.get_yticklabels()) + \
          [t for a_ in (axl, ax1, ax2, caxc, caxi) for t in a_.get_yticklabels() + a_.get_xticklabels() if t.get_text()]
    BX = [bb_in(t.get_window_extent(R)) for t in TXT]
    for i in range(len(BX)):
        for j in range(i + 1, len(BX)):
            a, b = BX[i], BX[j]
            ov = min(a[2], b[2]) - max(a[0], b[0]) > 0.005 and min(a[3], b[3]) - max(a[1], b[1]) > 0.005
            assert not ov, ('overlap', TXT[i].get_text(), TXT[j].get_text())
    # cell text fits inside its cell (all bottom heatmaps)
    for a_ in (axl, ax1, ax2):
        cw = bb_in(a_.get_window_extent(R)); ncol = len(a_.get_xticklabels()); cwi = (cw[2] - cw[0]) / ncol
        for t in a_.texts:
            tb = bb_in(t.get_window_extent(R)); assert tb[2] - tb[0] < cwi - 0.06 and tb[3] - tb[1] < CH - 0.08, (t.get_text(), tb, cwi)
    print(f'bottom row bottom {YB:.2f} in (composition bottom {YT:.2f} in); gap between heatmap bodies {YT - (YB + HH):.2f} in')
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)

    # ============================================================== source data (one tidy CSV: every plotted cell)
    rows = []
    for c, lab in zip(CNS, LAB):
        for k, col in enumerate(cols):
            rows.append(dict(heatmap='Cell-Type Composition', CN=c, CN_label=lab, CN_name=NMF.loc[c], column=XL[k], source_column=col,
                             value=E.loc[lab, col], value_plotted=H.loc[lab, col], cell_text='',
                             metric='log2((fraction of the cell type among the CN cells + 1e-4) / (fraction among all cells + 1e-4)); colour clipped to +/-3',
                             model='descriptive', source_csv=os.path.basename(C45)))
    LSi = LS.set_index(['module', 'CN']); Si = S.set_index(['axis', 'CN'])
    KEEP = ['n_sections', 'mean_ci_low', 'mean_ci_high', 'sections_above_0', 'lmm_estimate', 'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q', 'wilcoxon_q']
    for hm, src, keys, tx, srcf, extra in (('Ligand–Receptor Modules', LSi, MODS, TXL, C82D, 'n_pairs'), ('Inflammatory Signalling', Si, AXN, TX, C86, 'lr_pairs')):
        for c, lab in zip(CNS, LAB):
            for k in keys:
                r = src.loc[(k, c)]
                rows.append(dict(heatmap=hm, CN=c, CN_label=lab, CN_name=NMF.loc[c], column=k, source_column=k, value=r.mean_within_section_z,
                                 value_plotted=r.mean_within_section_z, cell_text=tx.loc[c, k], stars=stars(r.lmm_q), **{q_: r[q_] for q_ in KEEP},
                                 n_pairs=r['n_pairs'] if extra == 'n_pairs' else len(r['lr_pairs'].split('; ')),
                                 lr_pairs=(r['lr_pairs'] if isinstance(r.get('lr_pairs'), str) else ''), metric=r.metric,
                                 model=r.model.replace('within module', 'within each module').replace('within axis', 'within each axis'),
                                 source_csv=os.path.basename(srcf)))
    OUTC = pd.DataFrame(rows)
    for c_ in ['CN', 'n_sections', 'sections_above_0', 'n_pairs']: OUTC[c_] = OUTC[c_].round().astype('Int64')
    OUTC['colour_scale'] = np.where(OUTC.heatmap == 'Cell-Type Composition', 'RdBu_r, -3 to 3 (own colourbar)', f'RdBu_r, -{VM:g} to {VM:g} (shared by both bottom heatmaps)')
    assert len(OUTC) == 10 * (NCc + 4 + 5)
    OUTC.apply(sig3).to_csv(f'{STEM}.csv', index=False)

    # ============================================================== legend (composition + LR modules + inflammatory axes)
    NCELL = pq.read_metadata(f'{CND}/cn_per_cell.parquet').num_rows
    body, meth86 = L86T.split('\n\nMethods. Inflammatory signalling axes by neighbourhood. ')
    first = 'Inflammatory signalling axes across cellular neighbourhoods (rows) by axis (columns), with the TNF axis set apart as a contrast. '
    assert body.startswith(first); body = body[len(first):]
    cells_txt = 'cells with a neighbourhood assignment)'; assert body.count(cells_txt) == 1 and meth86.count(cells_txt) == 1
    body = body.replace(cells_txt, 'cells with a neighbourhood assignment and spatial LIANA scores)'); meth86 = meth86.replace(cells_txt, 'cells with a neighbourhood assignment and spatial LIANA scores)')
    abbr = ' Neighbourhood names are abbreviated by their dominant cell types; full names are in the source data.'; assert body.endswith(abbr); body = body[:-len(abbr)]
    m_ = re.match(r'CD resections \(CosMx 6K; 11 sections, 4 patients; ([0-9,]+) cells with a neighbourhood assignment and spatial LIANA scores\)\. ', body); assert m_
    NLIANA = m_.group(1); body = body[m_.end():]
    key86 = ' * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. ' + f'Colour scale symmetric, ±{VM86:.1f}.'; assert body.count(key86) == 1
    body = body.replace(key86, '')                                          # star key and colour range are stated once for both bottom heatmaps
    src = ' The source data hold the per-section values'; assert body.count(src) == 1
    body = body.replace(src, ' supp_inflammatory_axes_by_neighbourhood.csv holds the per-section values')
    L82T = open(L82D).read().strip(); assert NLIANA in L82T                 # same cell population for modules and axes
    def cut(t, a, b):
        i = t.index(a); j = t.index(b, i) + len(b); return t[i:j]
    lr_score = cut(L82T, 'Ligand-receptor (LR) scores are cell-type-agnostic', 'computed for every cell.')
    lr_mods = cut(L82T, 'Modules (in panel order):', 'an overlapping receptor).').replace('Modules (in panel order):', 'Modules (in column order):')
    lr_drop = cut(L82T, 'Two pairs of the pre-specified 9-pair co-stimulatory set are excluded:', 'are in the source data for comparison.').replace(
        'are in the source data for comparison.', 'are in fig4_lr_modules_allcells_costim_dotplot.csv for comparison.')
    lr_comp = cut(L82T, 'No sender or receiver cell types are imposed.', '0 = neighbourhood average).')
    lr_sig = cut(L82T, 'Significance: linear mixed model on the per-section values', 'within each module;')
    assert 'Co-Stimulatory = canonical antigen-presenting cell' in lr_mods and '(7 pairs)' in lr_mods
    comp = ('log2 enrichment of each cell type in each neighbourhood relative to all cells, i.e. log2 of the cell type\'s fraction among the '
        f'neighbourhood\'s cells over its fraction among all {NCELL:,} cells with a neighbourhood assignment (pseudocount 10⁻⁴); descriptive. The '
        'colour scale is clipped at ±3 (unclipped values in fig4_cn_composition_heatmap.csv). Columns are grouped by lineage (Epithelial, Stromal, '
        'Myeloid, Lymphoid, Other). Mac S+SG+/XS-, SELENOP⁺SIGLEC1⁺ and SELENOP⁺XIST⁻ macrophages combined; Mac S+M+P+, '
        'SELENOP⁺MMP9⁺PLA2G2D⁺ macrophages; Mac S+M+S+, SELENOP⁺MMP9⁺SPP1⁺ macrophages. Mo-Mac includes inflammatory '
        'Mo-Macs, CD8 Trm combines GZMK⁺ and GZMK⁻ CD8 Trm, and Tfh includes exhausted CD4 Tfh.')
    shared = (f'Ligand–Receptor Modules and Inflammatory Signalling share one colour scale (Mean Within-Section z, symmetric ±{VM:.1f}; widened from '
          f'±{VM86:.1f} of the separate inflammatory-axes heatmap to cover the module values). Cell text: Mean Within-Section z with stars from a linear '
          'mixed model on the per-section values, z ~ in neighbourhood + (1 | section); modules, Benjamini–Hochberg (BH) correction over the 10 '
          'neighbourhoods within each module; inflammatory axes, BH over the 10 neighbourhoods within each axis; * q < 0.05, ** q < 0.01, *** q < 0.001; '
          'no star, q ≥ 0.05.')
    leg = ('Cellular neighbourhood characterisation in CD resections (CosMx 6K; 11 sections, 4 patients). All heatmaps have the same rows: the 10 '
       'neighbourhoods CN1–CN10 (k-means, k = 10, on the cell-type composition of each cell\'s 20 nearest neighbours; coloured squares as in the '
       'Fig 5 neighbourhood map), abbreviated by their dominant cell types; full names are in the source data (CN_name). '
       'Cell-Type Composition: ' + comp + ' '
       'Ligand–Receptor Modules: spatial ligand–receptor (LR) module activity by neighbourhood, based on all '
       f'{NLIANA} cells with a neighbourhood assignment and spatial LIANA scores. ' + lr_score + ' ' + lr_mods + ' ' + lr_drop + ' '
       'Inflammatory Signalling: signalling axes (columns) by neighbourhood, with the TNF axis set apart as a contrast. ' + body + ' ' + shared + ' '
       'Source data: supp_cn_niche_characterisation.csv (every plotted value, with q values and 95% CIs), drawn from fig4_cn_composition_heatmap.csv, '
       'fig4_lr_modules_allcells_costim_dotplot.csv and supp_inflammatory_axes_by_neighbourhood.csv.\n\n'
       'Methods. Cell-type composition by neighbourhood. Neighbourhoods were defined on 20-nearest-neighbour windows (k-means, k = 10; cell types '
       'with all Tregs as one type) and named by their most enriched cell types; SELENOP⁺ fibroblasts are not part of any neighbourhood and are not shown. '
       'For display, Mo-Mac and inflammatory Mo-Mac, GZMK⁺ and GZMK⁻ CD8 Trm, and Tfh and exhausted CD4 Tfh were merged and the enrichment '
       'recomputed from cell counts. Ligand–receptor modules by neighbourhood. ' + lr_score + ' ' + lr_comp + ' ' + lr_sig[:-1] + '. '
       'Inflammatory signalling axes by neighbourhood. ' + meth86 + '\n')
    for t in ('activation score', 'quiescent', 'territor', 'serial', 'blunt'): assert t not in leg.lower(), t
    assert 'panel' not in leg.lower().replace('6k panel', ''), [m.start() for m in re.finditer('panel', leg.lower())]
    legpart = leg.split('\n\nMethods. ')[0]; assert 'SELENOP⁺ fibroblast' not in legpart and 'SELENOP⁺ Fib' not in legpart
    open(f'{STEM}_legend.txt', 'w').write(leg)
    print('written', f'{STEM}.png', f'{STEM}.csv', f'{STEM}_legend.txt')


# ======================================================================================================================
# Section cn_signalling_table: Niche signalling heatmap and table (8 signalling columns x 10 CNs)
# Panels: Fig 5C, Supp 9B tables (fig4_cn_signaling_heatmap.csv)
# ======================================================================================================================
if _run('cn_signalling_table'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/89b_cn_signaling_heatmap.py'   # working folder of this step (relative paths resolve from here)
    """Compact niche signalling heatmap and its table. CD resections (CosMx 6K; 11 sections, 4 patients). Rows: all 10 niches (CN1-CN10; long compact labels
and tab10 colour squares as in S9A_cn_characterisation). Columns in three groups: T-Cell Signals (Antigen Presentation, Co-Stimulatory = 7-pair module),
Inflammatory Cytokines (IL-6, IL-1, OSM, IL-11 | TNF, set slightly apart), Stroma (Matrix–Integrin). The broad 34-pair Inflammatory Cytokine module is not
shown (the axes are its subsets). Values = Mean Within-Section z, stars from lmm_q (BH over the 10 niches within each column's family), taken from
supp_cn_niche_characterisation.csv (S9A_cn_characterisation; its cell_text holds the exact 2-dp text) and asserted against the sources (cn_lr_modules_costim
for the modules, cn_inflammatory_axes for the axes): values, q, stars and cell text. One diverging colourbar, +/-2. CN10 row outlined in red (#b2182b),
CN10 label red. Cell width = the widest cell text (measured); the two T-Cell Signals columns are wider because their 14-pt labels need it (equal within
a group).
Outputs: A_fig4_panels/fig4_cn_signaling_heatmap.png + .csv (read by S9B_cn_signalling) + _legend.txt.
Layout variants (environment variable VARIANT; the default outputs are not affected):
  VARIANT=square (or SQUARE=1) -> fig4_cn_signaling_heatmap_square.{png,csv} + _legend.txt: cells 0.60 in square (row height up to match), two-line cell
    text (value on top, stars below; largest of 13 / 12 / 11 pt that fits), column labels rotated 45 deg (their length no longer sets the column width),
    group headers wrapped when wider than their group; CN colour squares kept at the S9A_cn_characterisation size.
  VARIANT=mid -> fig4_cn_signaling_heatmap_mid.{png,csv} + _legend.txt: cells ~1.8:1 (width 0.75-0.80 in, the largest of 13 / 12 / 11 pt whose widest
    one-line text, e.g. '−1.21***', fits with 0.04 in padding per side; row height = width / 1.8 within 0.42-0.45 in), one-line cell text, column labels
    rotated 45 deg as in the square variant."""
    import os, sys
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX', 'run with OPENBLAS_CORETYPE=SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/'), sys.executable
    import numpy as np, pandas as pd, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker, seaborn as sns
    from matplotlib.patches import Rectangle
    from matplotlib.font_manager import FontProperties
    from seaborn.utils import relative_luminance
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/A_fig4_panels'
    C89 = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_niche_characterisation.csv'; C82D = f'{OUT}/fig4_lr_modules_allcells_costim_dotplot.csv'; C86 = f'{OUT}/supp_inflammatory_axes_by_neighbourhood.csv'
    VARIANT = os.environ.get('VARIANT', 'square' if os.environ.get('SQUARE') == '1' else 'default'); assert VARIANT in ('default', 'square', 'mid'), VARIANT
    SQUARE, MID = VARIANT == 'square', VARIANT == 'mid'; ROT = SQUARE or MID                   # square / mid variants; ROT = 45-deg column labels
    STEM = f'{OUT}/fig4_cn_signaling_heatmap' + ('' if VARIANT == 'default' else '_' + VARIANT)
    assert not {f'{STEM}.png', f'{STEM}.csv', f'{STEM}_legend.txt'} & {C89, C82D, C86}
    CNS = list(range(1, 11)); HI = '#b2182b'; VM = 2.0
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    GROUPS = [('T-Cell Signals', ['Antigen Presentation', 'Co-Stimulatory']),
              ('Inflammatory Cytokines', ['IL-6', 'IL-1', 'OSM', 'IL-11', 'TNF']),
              ('Stroma', ['Matrix–Integrin'])]
    COLS = [c for _, cs in GROUPS for c in cs]; GOF = {c: g for g, cs in GROUPS for c in cs}
    CLAB = {'Antigen Presentation': 'Antigen\nPresentation', 'Co-Stimulatory': 'Co-\nStimulatory', 'Matrix–Integrin': 'Matrix–\nIntegrin'}
    SRC = {c: ('Ligand–Receptor Modules' if c in ('Antigen Presentation', 'Co-Stimulatory', 'Matrix–Integrin') else 'Inflammatory Signalling') for c in COLS}

    # ============================================================== data: S9A_cn_characterisation table, asserted against cn_lr_modules_costim / cn_inflammatory_axes
    N89 = pd.read_csv(C89); N89 = N89[N89.heatmap.isin(['Ligand–Receptor Modules', 'Inflammatory Signalling'])]
    X = pd.concat([N89[(N89.heatmap == SRC[c]) & (N89.column == c)] for c in COLS], ignore_index=True)
    X['CN'] = X.CN.astype(int); assert len(X) == 10 * len(COLS) and not X.duplicated(['column', 'CN']).any()
    V = X.pivot(index='CN', columns='column', values='value').loc[CNS, COLS]; Qm = X.pivot(index='CN', columns='column', values='lmm_q').loc[CNS, COLS]
    TX = X.pivot(index='CN', columns='column', values='cell_text').loc[CNS, COLS]
    STR = X.pivot(index='CN', columns='column', values='stars').loc[CNS, COLS].fillna('')
    LR = pd.read_csv(C82D, low_memory=False); LR = LR[LR.row_type == 'summary']; AX = pd.read_csv(C86, low_memory=False); AX = AX[AX.row_type == 'summary']
    for c in COLS:
        s = (LR[LR.module == c] if SRC[c].startswith('Ligand') else AX[AX.axis == c]).assign(CN=lambda d: d.CN.astype(int)).set_index('CN').loc[CNS]
        assert np.array_equal(s.mean_within_section_z.values, V[c].values) and np.array_equal(s.lmm_q.values, Qm[c].values), c
        if 'cell_text' in s: assert (s.cell_text.values == TX[c].values).all(), c                      # axes: identical cell text
        for cn in CNS:
            t = TX.loc[cn, c]; num = t.rstrip('*')
            assert t[len(num):] == stars(Qm.loc[cn, c]) == STR.loc[cn, c] and abs(float(num) - V.loc[cn, c]) <= 0.0051, (c, cn, t)
    for (cn, c), t in {(10, 'Antigen Presentation'): '+0.62*', (8, 'Co-Stimulatory'): '+1.92***', (9, 'Matrix–Integrin'): '+1.74***',
                       (10, 'IL-6'): '+0.95**', (8, 'TNF'): '+1.30***'}.items():
        assert TX.loc[cn, c] == t, (cn, c, TX.loc[cn, c])
    assert np.abs(V.values).max() <= VM
    LAB = X.drop_duplicates('CN').set_index('CN').CN_label.loc[CNS].tolist(); NMF = X.drop_duplicates('CN').set_index('CN').CN_name.loc[CNS]
    assert all(l.startswith(f'CN{c}: ') for c, l in zip(CNS, LAB)) and LAB[9] == 'CN10: OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'
    ANN = TX.apply(lambda s: s.str.replace('-', '−'))
    PAL = dict(zip(CNS, sns.color_palette('tab10', 10).as_hex()))           # CN colours

    # ============================================================== layout (inches; x data unit = inch, y data unit = row)
    FS_NUM, FS_LAB, FS_TITLE, FS_CB = 13, 14, 16, 13
    plt.rcParams.update({'font.size': FS_LAB, 'font.family': 'DejaVu Sans', 'font.weight': 'normal', 'axes.titleweight': 'normal', 'axes.labelweight': 'normal'})
    CH = 0.60 if SQUARE else 0.38; GAPG, GAPT = 0.20, 0.13; PADC, SEPL, MINW = 0.07, 0.14, 0.70
    fig = plt.figure(figsize=(20, 10)); R = fig.canvas.get_renderer(); DPI = fig.dpi
    tw = lambda s, fs: R.get_text_width_height_descent(s, FontProperties(family='DejaVu Sans', size=fs), ismath=False)[0] / DPI
    if SQUARE:   # square cells; value / stars on two lines; largest font of 13, 12, 11 pt whose widest value fits with 0.04 in padding per side
        PADC = 0.04
        FS_NUM = next(f for f in (13, 12, 11) if max(tw(t.rstrip('*'), f) for c in COLS for t in ANN[c]) + 2 * PADC <= CH)
        W = {c: CH for c in COLS}
    elif MID:    # ~1.8:1 cells; one-line text; largest of 13 / 12 / 11 pt whose widest text fits a cell of at most 0.80 in
        PADC = 0.04
        FS_NUM = next(f for f in (13, 12, 11) if max(tw(t, f) for c in COLS for t in ANN[c]) + 2 * PADC <= 0.80)
        CWM = max(0.75, float(np.ceil((max(tw(t, FS_NUM) for c in COLS for t in ANN[c]) + 2 * PADC) * 100) / 100))
        CH = float(np.clip(round(CWM / 1.8, 2), 0.42, 0.45)); W = {c: CWM for c in COLS}
        print(f'mid: cell {CWM:.2f} x {CH:.2f} in ({CWM / CH:.2f}:1), cell text {FS_NUM} pt')
    else:
        WTXT = max(tw(t, FS_NUM) for c in COLS for t in ANN[c]) + 2 * PADC                     # widest cell text: one uniform cell width
        W = {c: max(MINW, WTXT, max(tw(l, FS_LAB) for l in CLAB.get(c, c).split('\n')) + SEPL) for c in COLS}   # wider only where a column label needs it
        for g, cs in GROUPS: W.update({c: max(W[k] for k in cs) for c in cs})                  # equal widths within a group
    XL = {}; x = 0.0
    for gi, (g, cs) in enumerate(GROUPS):
        if gi: x += GAPG
        for c in cs:
            if c == 'TNF': x += GAPT
            XL[c] = x; x += W[c]
    WT = x; HH = 10 * CH
    X0, Y0 = 4.6, 1.0; CBP, CBW = 0.22, 0.20
    ax = fig.add_axes([X0 / 20, Y0 / 10, WT / 20, HH / 10]); ax.set_xlim(0, WT); ax.set_ylim(10, 0); ax.axis('on')
    for sp in ax.spines.values(): sp.set_visible(False)
    cmap = plt.get_cmap('RdBu_r'); norm = matplotlib.colors.Normalize(-VM, VM)
    for c in COLS:
        for i, cn in enumerate(CNS):
            col = cmap(norm(V.loc[cn, c]))
            ax.add_patch(Rectangle((XL[c], i), W[c], 1, facecolor=col, edgecolor='white', lw=0.5))
            if SQUARE:   # value on top, stars below (same colour rule)
                num = ANN.loc[cn, c].rstrip('*'); st = ANN.loc[cn, c][len(num):]; tc = '.15' if relative_luminance(col) > .408 else 'w'
                ax.text(XL[c] + W[c] / 2, i + 0.40, num, ha='center', va='center', fontsize=FS_NUM, color=tc)
                if st: ax.text(XL[c] + W[c] / 2, i + 0.74, st, ha='center', va='center', fontsize=FS_NUM, color=tc)
            else:
                ax.text(XL[c] + W[c] / 2, i + 0.5, ANN.loc[cn, c], ha='center', va='center', fontsize=FS_NUM,
                        color='.15' if relative_luminance(col) > .408 else 'w')                                     # seaborn's annotation colour rule
    ax.add_patch(Rectangle((0, 9), WT, 1, fill=False, ec=HI, lw=3.0, zorder=10, clip_on=False))              # CN10 row
    SQX, SQW = 0.05, 0.13
    SQH = 0.84 if VARIANT == 'default' else 0.32 / CH                                         # square / mid: same 0.13 x 0.32 in squares as S9A_cn_characterisation
    for i, cn in enumerate(CNS): ax.add_patch(Rectangle((-SQX - SQW, i + (1 - SQH) / 2), SQW, SQH, color=PAL[cn], clip_on=False, lw=0))
    ax.set_yticks(np.arange(10) + 0.5); ax.set_yticklabels(LAB, fontsize=FS_LAB)
    ax.tick_params(axis='y', length=0, pad=(SQX + SQW + 0.07) * 72)
    ax.get_yticklabels()[9].set_color(HI)
    ax.set_xticks([XL[c] + W[c] / 2 for c in COLS])
    if ROT:
        ax.set_xticklabels(COLS, fontsize=FS_LAB, rotation=45, ha='right', va='center', rotation_mode='anchor'); ax.tick_params(axis='x', length=0, pad=7)
    else:
        ax.set_xticklabels([CLAB.get(c, c) for c in COLS], fontsize=FS_LAB, linespacing=1.0, va='top')
        ax.tick_params(axis='x', length=0, pad=4)
    ax.set_ylabel('Spatial Neighborhood', fontsize=FS_LAB)
    # group headers with a thin rule beneath, then the title
    GH = []
    for g, cs in GROUPS:
        a, b = XL[cs[0]], XL[cs[-1]] + W[cs[-1]]
        yr, yt = (-0.22, -0.32) if VARIANT == 'default' else (-0.0836 / CH, -0.1216 / CH)                # same offsets in inches as the default
        ax.plot([a + 0.04, b - 0.04], [yr, yr], color='k', lw=1.0, clip_on=False, solid_capstyle='butt')
        if ROT:
            gl = g.replace(' ', '\n', 1) if tw(g, FS_LAB) > (b - a) + 0.15 and ' ' in g else g   # wrap a header wider than its group
            GH.append(ax.text((a + b) / 2, yt, gl, ha='center', va='bottom', fontsize=FS_LAB, linespacing=1.0))
        else:
            GH.append(ax.text((a + b) / 2, -0.32, g, ha='center', va='bottom', fontsize=FS_LAB))
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    bb_in = lambda b: (b.x0 / DPI, b.y0 / DPI, b.x1 / DPI, b.y1 / DPI)
    gtop = max(bb_in(t.get_window_extent(R))[3] for t in GH)
    TT = fig.text((X0 + WT / 2) / 20, (gtop + 0.08) / 10, 'Niche Signaling', ha='center', va='bottom', fontsize=FS_TITLE)
    cax = fig.add_axes([(X0 + WT + CBP) / 20, (Y0 + 0.075 * HH) / 10, CBW / 20, 0.85 * HH / 10])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label('Mean Within-Section z', fontsize=FS_CB); cb.outline.set_visible(False); cax.tick_params(labelsize=FS_CB)
    cax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1.0))
    cax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−')))
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    # ---- checks: cell text inside cells; no overlaps among labels / headers / title / colourbar text; nothing bold
    for t in ax.texts:
        if t in GH: continue
        b = bb_in(t.get_window_extent(R)); xc, yc = t.get_position(); c = min(COLS, key=lambda k: abs(XL[k] + W[k] / 2 - xc))
        assert b[2] - b[0] <= W[c] - 2 * PADC + 0.01 and b[3] - b[1] < CH - 0.08, (t.get_text(), b, W[c])
        if ROT:   # text line box(es) inside their cell (full line boxes, not just ink)
            row = int(np.floor(yc)); cy = ax.transData.transform([(0, row + 1), (0, row)])[:, 1] / DPI; top = (None, cy.min(), None, cy.max())   # cell y-range (in)
            assert b[1] >= top[1] - 1e-6 and b[3] <= top[3] + 1e-6, (t.get_text(), row, b, top)
    if SQUARE:   # value line box above the stars line box in every cell
        vt = {}
        for t in ax.texts:
            if t in GH: continue
            vt.setdefault((round(t.get_position()[0], 4), int(np.floor(t.get_position()[1]))), []).append(t)
        for k, ts in vt.items():
            if len(ts) == 2:
                nb, sb = sorted((bb_in(t.get_window_extent(R)) for t in ts), key=lambda b: -b[3]); assert nb[1] >= sb[3] - 0.005, (k, nb, sb)
    TXT = [TT, ax.yaxis.label, cax.yaxis.label] + GH + list(ax.get_yticklabels()) + ([] if ROT else list(ax.get_xticklabels())) + [t for t in cax.get_yticklabels() if t.get_text()]
    BX = [bb_in(t.get_window_extent(R)) for t in TXT]
    for i in range(len(BX)):
        for j in range(i + 1, len(BX)):
            a, b = BX[i], BX[j]
            assert not (min(a[2], b[2]) - max(a[0], b[0]) > 0.005 and min(a[3], b[3]) - max(a[1], b[1]) > 0.005), ('overlap', TXT[i].get_text(), TXT[j].get_text())
    axb = bb_in(ax.get_window_extent(R))
    if ROT:   # rotated column labels: exact rotated boxes; no overlap with each other (checked in the label frame), below the CN10 box, clear of all other text
        def label_poly(t):
            w, h, d = R.get_text_width_height_descent(t.get_text(), t.get_fontproperties(), ismath=False)
            x_, y_ = t.get_transform().transform(t.get_unitless_position()); th = np.deg2rad(t.get_rotation())
            uv = np.array([[-w, -h / 2], [-w, h / 2], [0, h / 2], [0, -h / 2]])
            P = np.c_[x_ + uv[:, 0] * np.cos(th) - uv[:, 1] * np.sin(th), y_ + uv[:, 0] * np.sin(th) + uv[:, 1] * np.cos(th)] / DPI
            assert np.allclose([P[:, 0].min(), P[:, 1].min(), P[:, 0].max(), P[:, 1].max()], bb_in(t.get_window_extent(R)), atol=0.03), t.get_text()
            return P, th
        PL = [label_poly(t) for t in ax.get_xticklabels()]
        th = PL[0][1]; rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        RB = [(P @ rot) for P, _ in PL]; RB = [(r[:, 0].min(), r[:, 1].min(), r[:, 0].max(), r[:, 1].max()) for r in RB]
        for i in range(len(RB)):
            for j in range(i + 1, len(RB)):
                a, b = RB[i], RB[j]; assert not (min(a[2], b[2]) - max(a[0], b[0]) > 0.005 and min(a[3], b[3]) - max(a[1], b[1]) > 0.005), ('label overlap', i, j)
        for P, _ in PL:
            assert P[:, 1].max() < axb[1] - 0.03, 'column label touches the heatmap / CN10 box'
            pb = (P[:, 0].min(), P[:, 1].min(), P[:, 0].max(), P[:, 1].max())
            for t, b in zip(TXT, BX): assert not (min(pb[2], b[2]) - max(pb[0], b[0]) > 0.005 and min(pb[3], b[3]) - max(pb[1], b[1]) > 0.005), ('label vs text', t.get_text())
    assert all(bb_in(t.get_window_extent(R))[1] > axb[3] for t in GH + [TT]) and bb_in(TT.get_window_extent(R))[1] > max(bb_in(t.get_window_extent(R))[3] for t in GH)
    assert all(t.get_fontweight() in ('normal', 400) for t in TXT + list(ax.texts))
    assert ax.get_yticklabels()[9].get_color() == HI and all(t.get_color() in ('k', 'black', '#000000') or t.get_color() == plt.rcParams['ytick.labelcolor'] or t.get_color() == plt.rcParams['ytick.color'] for t in ax.get_yticklabels()[:9])
    tb = fig.get_tightbbox(R); print(f'figure size (tight): {tb.width:.2f} x {tb.height:.2f} in; column widths (in):', {c: round(W[c], 2) for c in COLS})
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)

    # ============================================================== source data
    KEEP = ['CN', 'CN_label', 'CN_name', 'column', 'value', 'cell_text', 'stars', 'n_sections', 'mean_ci_low', 'mean_ci_high', 'sections_above_0', 'lmm_estimate',
            'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q', 'wilcoxon_q', 'n_pairs', 'lr_pairs', 'metric', 'model', 'source_csv']
    O = X[KEEP].rename(columns={'value': 'mean_within_section_z'}).copy()
    O.insert(0, 'group', O.column.map(GOF)); O['column'] = pd.Categorical(O.column, COLS); O = O.sort_values(['column', 'CN']); O['column'] = O.column.astype(str)
    O['source_csv'] = O.source_csv + ' (via supp_cn_niche_characterisation.csv)'; O['highlighted'] = O.CN == 10
    for c_ in ['CN', 'n_sections', 'sections_above_0', 'n_pairs']: O[c_] = O[c_].round().astype('Int64')
    O.apply(sig3).to_csv(f'{STEM}.csv', index=False)

    # ============================================================== legend (US spelling)
    NSEC = X.groupby('column').n_sections.max().astype(int)
    NCELL = 2936007
    assert f'{NCELL:,}' in open('/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_niche_characterisation_legend.txt').read()
    leg = ('Niche signaling across the 10 spatial cellular neighborhoods (niches, CN1–CN10), CD resections (CosMx 6K; 11 sections, 4 patients; '
       f'{NCELL:,} cells assigned to a niche). Niches: k-means (k = 10) on the cell-type composition of each cell\'s 20 nearest neighbors, named by '
       'their most enriched cell types (abbreviated; full names in the source data); their cell-type composition is shown in the supplementary niche '
       'characterization heatmap (supp_cn_niche_characterisation.png). CN10, the OGN⁺RSPO3⁺/PLA2G2A⁺ fibroblast niche, is outlined in red. '
       'Columns are ligand–receptor (LR) modules and signaling axes scored from cell-type-agnostic spatial LIANA bivariate local scores (Gaussian '
       'kernel, 200 px bandwidth) for the NicheNet v2 prior-network pairs, computed for every cell; no sender or receiver cell types are imposed. '
       'T-Cell Signals: Antigen Presentation = CellChatDB MHC-II (9 pairs); Co-Stimulatory = canonical antigen-presenting cell (APC) → T-cell '
       'co-stimulation, i.e. pre-specified pairs with the ligand on the APC and the receptor on the T cell / Treg (7 pairs): CD80 → CD28; '
       'CD86 → CD28; ICOSLG → ICOS; TNFSF4 → TNFRSF4; TNFSF18 → TNFRSF18; TNFSF9 → TNFRSF9; CD70 → CD27. Inflammatory '
       'Cytokines: pre-specified axes, IL-6 = IL6 → IL6R, IL6ST; IL-1 = IL1A → IL1R1, IL1RAP and IL1B → IL1R1, IL1RAP; OSM = OSM → '
       'OSMR, LIFR, IL6ST; IL-11 = IL11 → IL11RA, IL6ST; TNF = TNF → TNFRSF1A, TNFRSF1B (set apart as a contrast). The broader 34-pair '
       'Inflammatory Cytokine LR module (CellChatDB IL1, IL6, TNF, IFN-II, IL12, IL17, OSM), of which these axes are subsets, is not shown here; it is '
       'in the supplementary niche characterization heatmap. Stroma: Matrix–Integrin = CellChatDB ECM-Receptor pathways (227 pairs). '
       'The same quantity underlies every column: a cell\'s score is the mean raw local score over the column\'s pairs; a niche\'s value in a section is '
       'the mean over all its cells (niche × section combinations with at least 100 cells), so it reflects both ligand–receptor expression and '
       'the abundance of interacting cells; values are z-standardized across the 10 niches within each section and averaged over sections (Mean '
       'Within-Section z; 0 = average of all niches). Sections in which a column shows no variation across niches are left out (number of sections: '
           + ', '.join(f'{c} {NSEC[c]}' for c in COLS) + '). Stars: linear mixed model on the per-section values, z ~ in niche + (1 | section), '
       'Benjamini–Hochberg (BH)-corrected over the 10 niches within each column; * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. '
       'Color scale symmetric, ±2. Source data (' + os.path.basename(STEM) + '.csv) hold the values, 95% CIs across sections, mixed-model estimates '
       'and q values.\n')
    for t in ('activation score', 'quiescent', 'territor', 'serial', 'blunt', 'colour', 'neighbour', 'signalling', 'characteris'):
        assert t not in leg.lower().replace('supp_cn_niche_characterisation', ''), t
    if SQUARE:
        k_ = 'Stars: linear mixed model on the per-section values'; assert leg.count(k_) == 1
        leg = leg.replace(k_, 'Each cell gives the Mean Within-Section z (top) with its stars below. ' + k_)
    open(f'{STEM}_legend.txt', 'w').write(leg)
    print('written', f'{STEM}.png', f'{STEM}.csv', f'{STEM}_legend.txt')


# ======================================================================================================================
# Section F5C_score_signalling: Treg Suppression Score forest plot and niche signalling heatmap sharing the CN rows
# Panels: Fig 5C (drawn with VM=1.2 WIDE=1 RUN_SECTIONS=F5C_score_signalling python <file>)
# ======================================================================================================================
if _run('F5C_score_signalling'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/90_cn_treg_score_and_signaling.py'   # working folder of this step (relative paths resolve from here)
    """Treg Suppression Score forest plot and niche signalling heatmap in one figure, sharing the CN rows. CD resections (CosMx 6K; 11 sections, 4 patients).
Values are read from upstream CSVs:
  forest   A_fig4_panels/fig4_cn_treg_state_forest.csv (+ _legend.txt; 12_neighbourhood_treg_score.py): estimate, 95% CI, mixed_q, tested (CN6 / CN9 not
           tested), section consistency
  heatmap  supp_cn_niche_characterisation.csv (S9A_cn_characterisation; exact 2-dp cell_text), asserted against fig4_lr_modules_allcells_costim_dotplot.csv
           (cn_lr_modules_costim, modules) and supp_inflammatory_axes_by_neighbourhood.csv (cn_inflammatory_axes, axes): values, q, stars, cell text
Layout, left to right: CN labels (long compact labels, tab10 colour squares) | forest panel ('Treg Suppression Score' / 'In vs Outside Neighborhood';
colours: q < 0.05 positive #c26a77, negative #5b7fa6, otherwise grey; CI bars, stars beside the CI, dashed 0; CN6 / CN9: 'Not tested (<100 Tregs)') |
small gap (widened only as far as needed so the rotated heatmap column labels clear the forest x label) | heatmap (cells 0.78 x 0.43 in as the
cn_signalling_table VARIANT=mid layout, one-line '+0.62*' text at the largest of 13 / 12 / 11 pt that fits, 45-deg column labels, group headers T-Cell
Signals / Inflammatory Cytokines / Stroma, TNF gap, +/-2 colourbar). Row height = heatmap cell height; forest and heatmap row centres asserted identical.
CN4, CN5 and CN10 row labels (shared by both panels) in bold (EMPH_CN); HILITE_CN10 = True adds one thin red (#b2182b) box across both panels and a red
CN10 label. Title 'Treg Suppression Score and Niche Signaling'. US spelling in all reader-facing text.
Outputs: /share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/fig4_cn_treg_score_and_signaling[_vm12|_vm10][_wide].png + .csv + _legend.txt.
Environment variables: VM=1.2 or VM=1.0 sets the heatmap colour range (default 2.0): colours saturate beyond +/-VM (colourbar extend='both'), the
cell-text colour (white on dark cells) is re-derived from the new colours by the same luminance rule and the legend states the saturation; WIDE=1 (with
VM=1.2): wider / less tall layout (rows 0.38 in, column labels at 28 deg with three long names on two lines, cells widened only as needed for label
clearance, <= 1.0 in). Fig 5C = VM=1.2 WIDE=1 (fig4_cn_treg_score_and_signaling_vm12_wide.*)."""
    import os, sys
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX', 'run with OPENBLAS_CORETYPE=SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/'), sys.executable
    import numpy as np, pandas as pd, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker, seaborn as sns
    from matplotlib.patches import Rectangle
    from matplotlib.font_manager import FontProperties
    from seaborn.utils import relative_luminance
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; P4 = f'{D}/panels/fig4_and_related'; OUT = f'{P4}/A_fig4_panels'
    C26 = f'{OUT}/fig4_cn_treg_state_forest.csv'; L26 = f'{OUT}/fig4_cn_treg_state_forest_legend.txt'
    C89 = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_niche_characterisation.csv'; L89 = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_niche_characterisation_legend.txt'
    C82D = f'{OUT}/fig4_lr_modules_allcells_costim_dotplot.csv'; C86 = f'{OUT}/supp_inflammatory_axes_by_neighbourhood.csv'
    VM = float(os.environ.get('VM', '2.0')); assert VM in (2.0, 1.2, 1.0), VM       # heatmap colour range (default +/-2; +/-1.2 and +/-1.0 variants)
    SAT = VM < 2.0                                                                    # variant: colours saturate beyond +/-VM
    WIDE = os.environ.get('WIDE') == '1'; assert not WIDE or VM == 1.2, 'the wide layout is defined for the VM=1.2 variant'
    STEM = '/share/fsmresfiles/UC/revision_jenny/Fig4A-F_treg_suppression_niches/fig4_cn_treg_score_and_signaling' + (f'_vm{int(round(VM * 10))}' if SAT else '') + ('_wide' if WIDE else '')
    assert not {f'{STEM}.png', f'{STEM}.csv', f'{STEM}_legend.txt'} & {C26, L26, C89, L89, C82D, C86}
    CNS = list(range(1, 11)); HI = '#b2182b'
    # CN4 (DC/macrophage APC niche), CN5 (inflammatory Mo-Mac niche) and CN10 row labels in bold; HILITE_CN10 = True adds a red CN10 box and label.
    EMPH_CN = (4, 5, 10); HILITE_CN10 = False
    stars = lambda q: '' if not np.isfinite(q) or q >= 0.05 else '***' if q < 0.001 else '**' if q < 0.01 else '*'
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    GROUPS = [('T-Cell Signals', ['Antigen Presentation', 'Co-Stimulatory']),
              ('Inflammatory Cytokines', ['IL-6', 'IL-1', 'OSM', 'IL-11', 'TNF']),
              ('Stroma', ['Matrix–Integrin'])]
    COLS = [c for _, cs in GROUPS for c in cs]; GOF = {c: g for g, cs in GROUPS for c in cs}
    SRC = {c: ('Ligand–Receptor Modules' if c in ('Antigen Presentation', 'Co-Stimulatory', 'Matrix–Integrin') else 'Inflammatory Signalling') for c in COLS}

    # ============================================================== data: heatmap (S9A_cn_characterisation table, asserted against cn_lr_modules_costim / cn_inflammatory_axes)
    N89 = pd.read_csv(C89); N89 = N89[N89.heatmap.isin(['Ligand–Receptor Modules', 'Inflammatory Signalling'])]
    X = pd.concat([N89[(N89.heatmap == SRC[c]) & (N89.column == c)] for c in COLS], ignore_index=True)
    X['CN'] = X.CN.astype(int); assert len(X) == 10 * len(COLS) and not X.duplicated(['column', 'CN']).any()
    V = X.pivot(index='CN', columns='column', values='value').loc[CNS, COLS]; Qm = X.pivot(index='CN', columns='column', values='lmm_q').loc[CNS, COLS]
    TX = X.pivot(index='CN', columns='column', values='cell_text').loc[CNS, COLS]
    STR = X.pivot(index='CN', columns='column', values='stars').loc[CNS, COLS].fillna('')
    LR = pd.read_csv(C82D, low_memory=False); LR = LR[LR.row_type == 'summary']; AX = pd.read_csv(C86, low_memory=False); AX = AX[AX.row_type == 'summary']
    for c in COLS:
        s = (LR[LR.module == c] if SRC[c].startswith('Ligand') else AX[AX.axis == c]).assign(CN=lambda d: d.CN.astype(int)).set_index('CN').loc[CNS]
        assert np.array_equal(s.mean_within_section_z.values, V[c].values) and np.array_equal(s.lmm_q.values, Qm[c].values), c
        if 'cell_text' in s: assert (s.cell_text.values == TX[c].values).all(), c
        for cn in CNS:
            t = TX.loc[cn, c]; num = t.rstrip('*')
            assert t[len(num):] == stars(Qm.loc[cn, c]) == STR.loc[cn, c] and abs(float(num) - V.loc[cn, c]) <= 0.0051, (c, cn, t)
    for (cn, c), t in {(10, 'Antigen Presentation'): '+0.62*', (8, 'Co-Stimulatory'): '+1.92***', (9, 'Matrix–Integrin'): '+1.74***',
                       (10, 'IL-6'): '+0.95**', (8, 'TNF'): '+1.30***'}.items():
        assert TX.loc[cn, c] == t, (cn, c, TX.loc[cn, c])
    assert np.abs(V.values).max() <= VM or SAT
    SATC = [(cn, c, TX.loc[cn, c]) for c in COLS for cn in CNS if abs(V.loc[cn, c]) > VM]   # cells beyond the colour range (saturated)
    print(f'colour range +/-{VM:g}; saturated cells ({len(SATC)}):', SATC)
    LAB = X.drop_duplicates('CN').set_index('CN').CN_label.loc[CNS].tolist(); NMF = X.drop_duplicates('CN').set_index('CN').CN_name.loc[CNS]
    assert all(l.startswith(f'CN{c}: ') for c, l in zip(CNS, LAB)) and LAB[9] == 'CN10: OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'
    ANN = TX.apply(lambda s: s.str.replace('-', '−'))
    PAL = dict(zip(CNS, sns.color_palette('tab10', 10).as_hex()))           # CN colours
    # ============================================================== data: forest (fig4_cn_treg_state_forest.csv, 12_neighbourhood_treg_score.py)
    F = pd.read_csv(C26).set_index('CN').loc[CNS]
    EST = 'treg_suppression_score_z_in_minus_out'
    assert F.tested.tolist() == [c not in (6, 9) for c in CNS] and F.loc[[6, 9], 'mixed_q'].isna().all() and F.loc[F.tested, 'mixed_q'].notna().all()
    assert ((F.ci_lo < F[EST]) & (F[EST] < F.ci_hi)).all() and (F.loc[[6, 9], 'n_tregs_in_CN'] < 100).all()
    TRIM26 = [('OGN_RSPO3 Fib / PLA2G2A Fib / Endothelial', 'OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'), ('Epithelial / Crypt Top Fib / Cycling', 'Epithelial / CryptTop / Cycl'),
              ('Plasma / ADAMDEC1 Fib', 'Plasma / ADAMDEC1'), ('Myofibroblast / Neuron / Pericyte', 'Myofibroblast / Neuron / Peri'),
              ('DC / Mac S_M_P / CD9 Mac', 'DC / Mac / CD9 Mac'), ('Crypt Top Fib / Inf Mo-Mac / Mac S_M_S', 'CryptTop / InfMo-Mac / Mac'),
              ('Activated Fib / Neutrophil / Mac S_M_S', 'Activ Fib / Neut / Mac'), ('B / Tfh / Resting T', 'B / Tfh / Rest T')]   # compact labels of the forest table
    def short26(x):
        for a, b in TRIM26: x = x.replace(a, b)
        return x
    assert [short26(x) for x in F.CN_name] == LAB, 'forest CN labels differ from the heatmap CN labels'
    FCOL = {cn: ('#9e9e9e' if not (F.loc[cn, 'mixed_q'] < 0.05) else '#c26a77' if F.loc[cn, EST] > 0 else '#5b7fa6') for cn in CNS if F.loc[cn, 'tested']}   # forest colours
    FST = {cn: stars(F.loc[cn, 'mixed_q']) for cn in CNS if F.loc[cn, 'tested']}
    assert FST[4] == '***' and FST[5] == '***' and FST[10] == '**' and FCOL[10] == '#5b7fa6' and FCOL[4] == '#c26a77'
    print('forest:', {cn: (round(F.loc[cn, EST], 3), FST[cn]) for cn in FST})

    # ============================================================== layout (inches; heatmap x data unit = inch, all y data units = rows)
    FS_LAB, FS_TITLE, FS_CB, FS_TICK, FS_STAR, FS_NT = 14, 16, 13, 13, 18, 11
    plt.rcParams.update({'font.size': FS_LAB, 'font.family': 'DejaVu Sans', 'font.weight': 'normal', 'axes.titleweight': 'normal', 'axes.labelweight': 'normal'})
    FW, FH = 24.0, 12.0
    fig = plt.figure(figsize=(FW, FH)); R = fig.canvas.get_renderer(); DPI = fig.dpi
    tw = lambda s, fs: R.get_text_width_height_descent(s, FontProperties(family='DejaVu Sans', size=fs), ismath=False)[0] / DPI
    fx = lambda x: x / FW; fy = lambda y: y / FH
    # heatmap cells: the VARIANT=mid rule of cn_signalling_table (largest of 13 / 12 / 11 pt whose widest one-line text fits a cell of at most 0.80 in; width >= 0.75; height = width / 1.8)
    PADC = 0.04; GAPG, GAPT = 0.20, 0.13
    FS_NUM = next(f for f in (13, 12, 11) if max(tw(t, f) for c in COLS for t in ANN[c]) + 2 * PADC <= 0.80)
    CW = max(0.75, float(np.ceil((max(tw(t, FS_NUM) for c in COLS for t in ANN[c]) + 2 * PADC) * 100) / 100)); CH = float(np.clip(round(CW / 1.8, 2), 0.42, 0.45))
    assert (FS_NUM, CW, CH) == (11, 0.78, 0.43), (FS_NUM, CW, CH)          # = cn_signalling_table VARIANT=mid
    if WIDE:   # wider / less tall: rows 0.38 in (as S9A_cn_characterisation); column labels at 28 deg with the three long names on two lines;
               # cells widened only as far as needed (<= 1.0 in) so neighbouring rotated labels keep >= 0.03 in clearance (exact boxes in the label frame)
        ROTA, LSP = 28, 1.0
        XLAB = {'Antigen Presentation': 'Antigen\nPresentation', 'Co-Stimulatory': 'Co-\nStimulatory', 'Matrix\u2013Integrin': 'Matrix\u2013\nIntegrin'}
        CH = 0.38
        def lab_wh(s_):
            t_ = fig.text(0, 0, s_, fontsize=FS_LAB, linespacing=LSP); b_ = t_.get_window_extent(R); t_.remove(); return b_.width / DPI, b_.height / DPI
        SZ = [lab_wh(XLAB.get(c, c)) for c in COLS]; th_ = np.deg2rad(ROTA)
        def collide(cw):
            xc = []; x_ = 0.0
            for gi, (g, cs) in enumerate(GROUPS):
                if gi: x_ += GAPG
                for c in cs:
                    if c == 'TNF': x_ += GAPT
                    xc.append(x_ + cw / 2); x_ += cw
            for i in range(len(COLS) - 1):
                d_ = xc[i + 1] - xc[i]
                if d_ * np.cos(th_) < SZ[i + 1][0] and d_ * np.sin(th_) < SZ[i][1] + 0.03: return True   # va='top': blocks hang below their anchor
            return False
        CW = next(round(float(c_), 2) for c_ in np.arange(0.78, 1.0001, 0.01) if not collide(c_))
        print(f'wide: labels at {ROTA} deg, cells {CW} x {CH} in')
    else:
        ROTA, LSP, XLAB = 45, None, {}
    XL = {}; x = 0.0
    for gi, (g, cs) in enumerate(GROUPS):
        if gi: x += GAPG
        for c in cs:
            if c == 'TNF': x += GAPT
            XL[c] = x; x += CW
    WT = x; HH = 10 * CH
    FX0, FWID, Y0 = 5.0, 2.9, 2.6                                           # forest axes left / width (in), common bottom of both panels
    GAP0 = 0.30; CBP, CBW = 0.22, 0.20
    axf = fig.add_axes([fx(FX0), fy(Y0), fx(FWID), fy(HH)])
    axh = fig.add_axes([fx(FX0 + FWID + GAP0), fy(Y0), fx(WT), fy(HH)])
    cax = fig.add_axes([fx(FX0 + FWID + GAP0 + WT + CBP), fy(Y0 + 0.075 * HH), fx(CBW), fy(0.85 * HH)])
    # ---------------- forest panel
    yc = np.arange(10) + 0.5
    axf.set_ylim(10, 0); axf.axvline(0, color='0.4', lw=1, ls='--', zorder=1)
    TST = []
    for i, cn in enumerate(CNS):
        if not F.loc[cn, 'tested']: continue
        lo, hi, e = F.loc[cn, 'ci_lo'], F.loc[cn, 'ci_hi'], F.loc[cn, EST]
        axf.hlines(yc[i], lo, hi, color='0.35', lw=2.0, zorder=2)
        axf.scatter([e], [yc[i]], c=[FCOL[cn]], s=150, zorder=3, edgecolor='k', lw=0.7)
        if FST[cn]: TST.append(axf.text(hi + 0.02 if e >= 0 else lo - 0.02, yc[i], FST[cn], ha='left' if e >= 0 else 'right', va='center', fontsize=FS_STAR))
    axf.set_yticks(yc); axf.set_yticklabels(LAB, fontsize=FS_LAB)
    SQX, SQW, SQH = 0.05, 0.13, 0.32                                        # CN colour squares (in), as cn_signalling_table
    axf.tick_params(axis='y', length=0, pad=(SQX + SQW + 0.07) * 72)
    if HILITE_CN10: axf.get_yticklabels()[9].set_color(HI)
    for i, cn in enumerate(CNS):
        if cn in EMPH_CN: axf.get_yticklabels()[i].set_fontweight('bold')   # bold EMPH_CN row labels
    axf.set_ylabel('Spatial Neighborhood', fontsize=FS_LAB)
    axf.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([-0.4, -0.2, 0.0, 0.2, 0.4]))
    axf.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−')))
    axf.tick_params(axis='x', labelsize=FS_TICK)
    axf.set_xlabel('Treg Suppression Score (z)\n← Low-suppressive   Activated →', fontsize=FS_TICK, linespacing=1.3)
    for k_ in ('top', 'right', 'left'): axf.spines[k_].set_visible(False)
    axf.spines['bottom'].set_position(('outward', 6))                      # x axis set just below the rows, clear of the CN10 box
    x0, x1 = min(-0.45, F.loc[F.tested, 'ci_lo'].min() - 0.02), max(0.45, F.loc[F.tested, 'ci_hi'].max() + 0.02)   # x range: ticks, CIs and stars
    for it in range(50):
        axf.set_xlim(x0, x1); fig.canvas.draw(); R = fig.canvas.get_renderer(); inv = axf.transData.inverted(); m = 3 * DPI / 72
        bx = [t.get_window_extent(R) for t in TST]
        n0 = min([x0] + [inv.transform((b.x0 - m, 0))[0] for b in bx]); n1 = max([x1] + [inv.transform((b.x1 + m, 0))[0] for b in bx])
        if abs(n0 - x0) < 1e-6 and abs(n1 - x1) < 1e-6: break
        x0, x1 = n0, n1
    assert it < 49; axf.set_xlim(x0, x1)
    NT = [axf.text((x0 + x1) / 2, yc[i], 'Not tested (<100 Tregs)', ha='center', va='center', fontsize=FS_NT, color='0.45', style='italic',
                   bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'), zorder=4) for i, cn in enumerate(CNS) if not F.loc[cn, 'tested']]
    # ---------------- heatmap (cn_signalling_table VARIANT=mid layout)
    axh.set_xlim(0, WT); axh.set_ylim(10, 0)
    for sp in axh.spines.values(): sp.set_visible(False)
    cmap = plt.get_cmap('RdBu_r'); norm = matplotlib.colors.Normalize(-VM, VM)
    for c in COLS:
        for i, cn in enumerate(CNS):
            col = cmap(norm(V.loc[cn, c]))
            axh.add_patch(Rectangle((XL[c], i), CW, 1, facecolor=col, edgecolor='white', lw=0.5))
            axh.text(XL[c] + CW / 2, i + 0.5, ANN.loc[cn, c], ha='center', va='center', fontsize=FS_NUM, color='.15' if relative_luminance(col) > .408 else 'w')
    axh.set_yticks([]); axh.set_xticks([XL[c] + CW / 2 for c in COLS])
    axh.set_xticklabels([XLAB.get(c, c) for c in COLS], fontsize=FS_LAB, rotation=ROTA, ha='right', va='top' if WIDE else 'center', rotation_mode='anchor', **({'linespacing': LSP} if WIDE else {})); axh.tick_params(axis='x', length=0, pad=5 if WIDE else 7)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, **({'extend': 'both'} if SAT else {}))
    cb.set_label('Mean Within-Section z', fontsize=FS_CB); cb.outline.set_visible(False); cax.tick_params(labelsize=FS_CB)
    cax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1.0 if not SAT else 0.5))
    cax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−')))
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    bb_in = lambda b: (b.x0 / DPI, b.y0 / DPI, b.x1 / DPI, b.y1 / DPI)
    def label_poly(t):
        if '\n' in t.get_text() or t.get_va() == 'top':   # multi-line / top-anchored label: unrotated size of the same laid-out text block
            t_ = fig.text(0, 0, t.get_text(), fontproperties=t.get_fontproperties(), linespacing=t._linespacing); b_ = t_.get_window_extent(R); t_.remove(); w, h = b_.width, b_.height
        else:
            w, h, d = R.get_text_width_height_descent(t.get_text(), t.get_fontproperties(), ismath=False)
        x_, y_ = t.get_transform().transform(t.get_unitless_position()); th = np.deg2rad(t.get_rotation())
        v0, v1 = (-h / 2, h / 2) if t.get_va() == 'center' else (-h, 0.0); assert t.get_va() in ('center', 'top')
        uv = np.array([[-w, v0], [-w, v1], [0, v1], [0, v0]])
        P = np.c_[x_ + uv[:, 0] * np.cos(th) - uv[:, 1] * np.sin(th), y_ + uv[:, 0] * np.sin(th) + uv[:, 1] * np.cos(th)] / DPI
        assert np.allclose([P[:, 0].min(), P[:, 1].min(), P[:, 0].max(), P[:, 1].max()], bb_in(t.get_window_extent(R)), atol=0.03), t.get_text()
        return P
    def poly_min_x(P, a, b):   # leftmost x of a convex polygon within the horizontal band y in [a, b]
        lo, hi = max(a, P[:, 1].min()), min(b, P[:, 1].max())
        if lo > hi: return None
        xs = []
        for yy in np.linspace(lo, hi, 200):
            for i in range(4):
                p, q = P[i], P[(i + 1) % 4]
                if min(p[1], q[1]) - 1e-9 <= yy <= max(p[1], q[1]) + 1e-9 and abs(q[1] - p[1]) > 1e-9:
                    xs.append(p[0] + (yy - p[1]) * (q[0] - p[0]) / (q[1] - p[1]))
        return min(xs) if xs else None
    # gap between the panels: GAP0, widened only as far as needed so the rotated heatmap labels clear the forest x label / tick labels by 0.10 in
    FB = [bb_in(axf.xaxis.label.get_window_extent(R))] + [bb_in(t.get_window_extent(R)) for t in axf.get_xticklabels() if t.get_text()]
    need = 0.0
    for t in axh.get_xticklabels():
        P = label_poly(t)
        for b in FB:
            mx = poly_min_x(P, b[1] - 0.03, b[3] + 0.03)
            if mx is not None: need = max(need, b[2] + 0.10 - mx)
    GAP = GAP0 + need
    axh.set_position([fx(FX0 + FWID + GAP), fy(Y0), fx(WT), fy(HH)]); cax.set_position([fx(FX0 + FWID + GAP + WT + CBP), fy(Y0 + 0.075 * HH), fx(CBW), fy(0.85 * HH)])
    XH = FX0 + FWID + GAP
    # headers: heatmap group headers with a rule (cn_signalling_table offsets), forest header + subheader with the same rule height; title above all
    YTOP = Y0 + HH; RULE, HBOT = 0.0836, 0.1216
    GH = []
    for g, cs in GROUPS:
        a, b = XH + XL[cs[0]], XH + XL[cs[-1]] + CW
        fig.add_artist(matplotlib.lines.Line2D([fx(a + 0.04), fx(b - 0.04)], [fy(YTOP + RULE)] * 2, color='k', lw=1.0, solid_capstyle='butt'))
        GH.append(fig.text(fx((a + b) / 2), fy(YTOP + HBOT), g, ha='center', va='bottom', fontsize=FS_LAB))
    fig.add_artist(matplotlib.lines.Line2D([fx(FX0 + 0.04), fx(FX0 + FWID - 0.04)], [fy(YTOP + RULE)] * 2, color='k', lw=1.0, solid_capstyle='butt'))
    FSUB = fig.text(fx(FX0 + FWID / 2), fy(YTOP + HBOT), 'In vs Outside Neighborhood', ha='center', va='bottom', fontsize=FS_TICK)
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    FHD = fig.text(fx(FX0 + FWID / 2), fy(bb_in(FSUB.get_window_extent(R))[3] + 0.03), 'Treg Suppression Score', ha='center', va='bottom', fontsize=FS_LAB)
    fig.canvas.draw(); R = fig.canvas.get_renderer()
    htop = max(bb_in(t.get_window_extent(R))[3] for t in GH + [FHD, FSUB])
    TT = fig.text(fx((FX0 + XH + WT) / 2), fy(htop + 0.10), 'Treg Suppression Score and Niche Signaling', ha='center', va='bottom', fontsize=FS_TITLE)
    # CN colour squares and the optional thin red box around the CN10 row across both panels
    for i, cn in enumerate(CNS):
        fig.add_artist(Rectangle((fx(FX0 - SQX - SQW), fy(YTOP - (i + 1) * CH + (CH - SQH) / 2)), fx(SQW), fy(SQH), color=PAL[cn], lw=0, transform=fig.transFigure))
    BOX = Rectangle((fx(FX0), fy(Y0)), fx(XH + WT - FX0), fy(CH), fill=False, ec=HI, lw=2.0, transform=fig.transFigure, zorder=20)
    if HILITE_CN10: fig.add_artist(BOX)   # CN10 box only when HILITE_CN10
    fig.canvas.draw(); R = fig.canvas.get_renderer()

    # ============================================================== checks
    # rows: forest and heatmap row centres / edges identical (display coordinates)
    for i in range(10):
        yf = axf.transData.transform([(0, i), (0, i + 0.5), (0, i + 1)])[:, 1]; yh = axh.transData.transform([(0, i), (0, i + 0.5), (0, i + 1)])[:, 1]
        assert np.allclose(yf, yh, atol=1e-6), ('row misaligned', i, yf, yh)
        assert abs(abs(yf[0] - yf[2]) / DPI - CH) < 1e-6
        t = bb_in(axf.get_yticklabels()[i].get_window_extent(R)); assert abs((t[1] + t[3]) / 2 - yf[1] / DPI) < 0.03, ('label not on its row', i)
    # CN10 box = the CN10 row exactly, spanning both panels
    bb = bb_in(BOX.get_window_extent(R)); r10 = axh.transData.transform([(0, 9), (0, 10)])[:, 1] / DPI
    assert abs(bb[1] - r10.min()) < 1e-3 and abs(bb[3] - r10.max()) < 1e-3 and abs(bb[0] - FX0) < 1e-3 and abs(bb[2] - (XH + WT)) < 1e-3
    assert (BOX in fig.artists) == HILITE_CN10 and (axf.get_yticklabels()[9].get_color() == HI) == HILITE_CN10
    # forest: stars inside the axes, beside (not on) their CI, clear of the zero line; dots / CIs inside the x range; not-tested text inside the axes
    AB = bb_in(axf.get_window_extent(R)); z0 = axf.transData.transform((0, 0))[0] / DPI
    for t in TST:
        b = bb_in(t.get_window_extent(R)); assert b[0] >= AB[0] and b[2] <= AB[2] and not (b[0] <= z0 <= b[2]), t.get_text()
    for t in NT:
        b = bb_in(t.get_window_extent(R)); assert b[0] >= AB[0] and b[2] <= AB[2], 'not-tested text outside the forest panel'
    assert x0 < F.loc[F.tested, 'ci_lo'].min() and x1 > F.loc[F.tested, 'ci_hi'].max()
    # heatmap cell text inside cells
    for t in axh.texts:
        b = bb_in(t.get_window_extent(R)); assert b[2] - b[0] <= CW - 2 * PADC + 0.01 and b[3] - b[1] < CH - 0.08, t.get_text()
    # no overlaps among all upright text; rotated labels checked as exact boxes
    TXT = [TT, FHD, FSUB, axf.xaxis.label, axf.yaxis.label, cax.yaxis.label] + GH + TST + NT + list(axf.get_yticklabels()) + \
          [t for t in axf.get_xticklabels() + cax.get_yticklabels() if t.get_text()]
    BX = [bb_in(t.get_window_extent(R)) for t in TXT]
    ov = lambda a, b: min(a[2], b[2]) - max(a[0], b[0]) > 0.005 and min(a[3], b[3]) - max(a[1], b[1]) > 0.005
    for i in range(len(BX)):
        for j in range(i + 1, len(BX)): assert not ov(BX[i], BX[j]), ('overlap', TXT[i].get_text(), TXT[j].get_text())
    PL = [label_poly(t) for t in axh.get_xticklabels()]; th = np.deg2rad(axh.get_xticklabels()[0].get_rotation())
    rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]); RB = [(P @ rot) for P in PL]; RB = [(r[:, 0].min(), r[:, 1].min(), r[:, 0].max(), r[:, 1].max()) for r in RB]
    for i in range(len(RB)):
        for j in range(i + 1, len(RB)): assert not ov(RB[i], RB[j]), ('column label overlap', i, j)
    HB = bb_in(axh.get_window_extent(R))
    for P in PL:
        assert P[:, 1].max() < HB[1] - 0.03, 'column label touches the heatmap / CN10 box'
        for t, b in zip(TXT, BX):
            mx = poly_min_x(P, b[1], b[3])
            if mx is not None and t in [axf.xaxis.label] + list(axf.get_xticklabels()): assert mx >= b[2] + 0.099, ('rotated label vs forest text', t.get_text())
            pb = (P[:, 0].min(), P[:, 1].min(), P[:, 0].max(), P[:, 1].max())
            if t not in [axf.xaxis.label] + list(axf.get_xticklabels()): assert not ov(pb, b), ('rotated label vs text', t.get_text())
    assert all(bb_in(t.get_window_extent(R))[1] > HB[3] for t in GH + [FHD, FSUB, TT]) and bb_in(TT.get_window_extent(R))[1] > htop
    BOLDT = [axf.get_yticklabels()[CNS.index(c)] for c in EMPH_CN]   # only the EMPH_CN (CN4 / CN5 / CN10) row labels are bold
    assert all(t.get_fontweight() == 'bold' for t in BOLDT) and all(t.get_text().startswith(f'CN{c}: ') for t, c in zip(BOLDT, EMPH_CN))
    assert all(t.get_fontweight() in ('normal', 400) for t in TXT + list(axh.texts) + list(axh.get_xticklabels()) if t not in BOLDT)
    EXTRA = fig.get_default_bbox_extra_artists() + [axf.xaxis.label]   # keep figure-level texts (title, headers) AND the long forest x label
    tb = fig.get_tightbbox(R, bbox_extra_artists=EXTRA)
    for t in [TT, FHD, FSUB] + GH: b = bb_in(t.get_window_extent(R)); assert tb.x0 <= b[0] and b[2] <= tb.x1 and tb.y0 <= b[1] and b[3] <= tb.y1, ('cropped', t.get_text())
    print(f'gap between panels {GAP:.2f} in (minimum {GAP0}); cells {CW} x {CH} in, cell text {FS_NUM} pt; figure (tight) {tb.width:.2f} x {tb.height:.2f} in')
    fig.savefig(f'{STEM}.png', dpi=300, bbox_inches='tight', pad_inches=0.05, bbox_extra_artists=EXTRA); plt.close(fig)

    # ============================================================== source data (one CSV: forest rows + heatmap cells)
    FO = F.reset_index().rename(columns={EST: 'estimate', 'ci_lo': 'ci_low', 'ci_hi': 'ci_high', 'mixed_p': 'p', 'mixed_q': 'q'})
    FO = FO.assign(panel='Treg Suppression Score (forest)', column='Treg Suppression Score, in vs outside neighborhood', CN_label=LAB,
                   stars=[FST.get(c, '') for c in CNS], color=[FCOL.get(c, 'not tested') for c in CNS],
                   model='Treg Suppression Score (z within section) ~ in_CN + (1 | section); BH over the 8 tested neighborhoods', source_csv=os.path.basename(C26))
    FO = FO[['panel', 'CN', 'CN_label', 'CN_name', 'column', 'estimate', 'ci_low', 'ci_high', 'p', 'q', 'stars', 'tested', 'n_tregs_in_CN', 'n_sections_ge10_tregs',
             'sections_agreeing', 'sections', 'color', 'model', 'source_csv']]
    KEEP = ['CN', 'CN_label', 'CN_name', 'column', 'value', 'cell_text', 'stars', 'n_sections', 'mean_ci_low', 'mean_ci_high', 'sections_above_0', 'lmm_estimate',
            'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q', 'wilcoxon_q', 'n_pairs', 'lr_pairs', 'metric', 'model', 'source_csv']
    HO = X[KEEP].rename(columns={'value': 'mean_within_section_z'}).copy(); HO.insert(0, 'panel', 'Niche Signaling (heatmap)'); HO.insert(1, 'group', HO.column.map(GOF))
    HO['column'] = pd.Categorical(HO.column, COLS); HO = HO.sort_values(['column', 'CN']); HO['column'] = HO.column.astype(str)
    HO['source_csv'] = HO.source_csv + ' (via ' + os.path.basename(C89) + ')'
    OUTC = pd.concat([FO, HO], ignore_index=True); OUTC['highlighted'] = OUTC.CN == 10
    for c_ in ['CN', 'n_tregs_in_CN', 'n_sections_ge10_tregs', 'sections_agreeing', 'sections', 'n_sections', 'sections_above_0', 'n_pairs']: OUTC[c_] = OUTC[c_].round().astype('Int64')
    assert len(OUTC) == 10 + 10 * len(COLS)
    OUTC.apply(sig3).to_csv(f'{STEM}.csv', index=False)

    # ============================================================== legend (forest legend text + niche signaling text; US spelling)
    T26 = open(L26).read().strip()
    a26 = 'Dot, difference in the continuous Treg Suppression Score'; b26 = ' Neighbourhood names are abbreviated'
    assert T26.count(a26) == 1 and T26.count(b26) == 1
    f26 = T26[T26.index(a26):T26.index(b26)]
    NTREG = T26[T26.index('all ') + 4:T26.index(' Tregs).')]; assert NTREG == '6,554', NTREG
    US = [('neighbourhoods', 'neighborhoods'), ('Neighbourhoods', 'Neighborhoods'), ('neighbourhood', 'neighborhood'), ('Neighbourhood', 'Neighborhood'), ('colour', 'color'), ('grey', 'gray')]
    for a, b in US: f26 = f26.replace(a, b)
    f26 = f26.replace('were not tested.', 'were not tested (shown as "Not tested").'); assert 'shown as "Not tested"' in f26
    NSEC = X.groupby('column').n_sections.max().astype(int)
    NCELL = 2936007; assert f'{NCELL:,}' in open(L89).read()
    leg = ('Treg Suppression Score and niche signaling across the 10 spatial cellular neighborhoods (niches, CN1–CN10), CD resections (CosMx 6K; '
       f'11 sections, 4 patients; {NTREG} Tregs; {NCELL:,} cells assigned to a niche). Both panels share the rows. Niches: k-means (k = 10) on the '
       'cell-type composition of each cell\'s 20 nearest neighbors, named by their most enriched cell types (abbreviated; full names in the source data); '
       'their cell-type composition is shown in the supplementary niche characterization heatmap (supp_cn_niche_characterisation.png). '
           + ('CN10, the OGN⁺RSPO3⁺/PLA2G2A⁺ fibroblast niche, is outlined in red across both panels. ' if HILITE_CN10 else '')
           + ((', '.join(f'CN{c}' for c in EMPH_CN[:-1]) + ' and ' if len(EMPH_CN) > 1 else '') + f'CN{EMPH_CN[-1]} labels are shown in bold. ' if EMPH_CN else '') +
           'Treg Suppression Score (left; in vs outside the neighborhood). ' + f26.strip() + ' '
       'Niche Signaling (right). Columns are ligand–receptor (LR) modules and signaling axes scored from cell-type-agnostic spatial LIANA bivariate '
       'local scores (Gaussian kernel, 200 px bandwidth) for the NicheNet v2 prior-network pairs, computed for every cell; no sender or receiver cell types '
       'are imposed. T-Cell Signals: Antigen Presentation = CellChatDB MHC-II (9 pairs); Co-Stimulatory = canonical antigen-presenting cell (APC) → '
       'T-cell co-stimulation, i.e. pre-specified pairs with the ligand on the APC and the receptor on the T cell / Treg (7 pairs): CD80 → CD28; '
       'CD86 → CD28; ICOSLG → ICOS; TNFSF4 → TNFRSF4; TNFSF18 → TNFRSF18; TNFSF9 → TNFRSF9; CD70 → CD27. Inflammatory '
       'Cytokines: pre-specified axes, IL-6 = IL6 → IL6R, IL6ST; IL-1 = IL1A → IL1R1, IL1RAP and IL1B → IL1R1, IL1RAP; OSM = OSM → '
       'OSMR, LIFR, IL6ST; IL-11 = IL11 → IL11RA, IL6ST; TNF = TNF → TNFRSF1A, TNFRSF1B (set apart as a contrast). The broader 34-pair '
       'Inflammatory Cytokine LR module (CellChatDB IL1, IL6, TNF, IFN-II, IL12, IL17, OSM), of which these axes are subsets, is not shown here; it is '
       'in the supplementary niche characterization heatmap. Stroma: Matrix–Integrin = CellChatDB ECM-Receptor pathways (227 pairs). '
       'The same quantity underlies every column: a cell\'s score is the mean raw local score over the column\'s pairs; a niche\'s value in a section is '
       'the mean over all its cells (niche × section combinations with at least 100 cells), so it reflects both ligand–receptor expression and '
       'the abundance of interacting cells; values are z-standardized across the 10 niches within each section and averaged over sections (Mean '
       'Within-Section z; 0 = average of all niches). Sections in which a column shows no variation across niches are left out (number of sections: '
           + ', '.join(f'{c} {NSEC[c]}' for c in COLS) + '). Stars: linear mixed model on the per-section values, z ~ in niche + (1 | section), '
       'Benjamini–Hochberg (BH)-corrected over the 10 niches within each column; * q < 0.05, ** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. '
           + ('Color scale symmetric, ±2. ' if not SAT else f'Color scale symmetric, ±{VM:.1f}; values beyond ±{VM:.1f} are shown in the end color (pointed '
          'colorbar ends), and the cell numbers give the exact values. ') + 'Source data (' + os.path.basename(STEM) + '.csv) hold both panels: the forest estimates, 95% CIs, q values and '
       'section counts, and the heatmap values, 95% CIs across sections, mixed-model estimates and q values.\n')
    for t in ('activation score', 'quiescent', 'territor', 'serial', 'blunt', 'colour', 'neighbour', 'signalling', 'grey'):
        assert t not in leg.lower().replace('supp_cn_niche_characterisation', ''), t
    open(f'{STEM}_legend.txt', 'w').write(leg)
    print('written', f'{STEM}.png', f'{STEM}.csv', f'{STEM}_legend.txt')


# ======================================================================================================================
# Section S9B_cn_signalling: Per-niche mixed-model estimate and 95% CI for every column of the niche signalling heatmap
# Panels: Supp 9B
# ======================================================================================================================
if _run('S9B_cn_signalling'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/82g_cn_signaling_errorbars.py'   # working folder of this step (relative paths resolve from here)
    """Per niche, the mixed-model estimate (that niche vs the other nine) with its 95% CI, i.e. the model that gives the heatmap's stars, for every column of
the niche signalling heatmap (fig4_cn_signaling_heatmap.csv, cn_signalling_table). CD resections (CosMx 6K; 11 sections, 4 patients).
Panels = the heatmap's 8 columns in its order and groups, 2 columns x 4 rows: T-Cell Signals (Antigen Presentation, Co-Stimulatory [7 pairs]),
Inflammatory Cytokines (IL-6, IL-1, OSM, IL-11, TNF), Stroma (Matrix-Integrin). Rows = the 10 niches (CN1 top, CN10 bottom; long labels + tab10 squares on
the left column only, right column shares y). Per niche: point = LMM estimate (z ~ in_CN + (1 | section); in_CN coefficient = niche minus the other
niches, within-section z), bar = its Wald 95% CI, filled if q < 0.05 else open; dashed line at 0; no section dots; stars (significant only; BH over the
10 niches within each column, as the heatmap) in one column at the right edge of each panel; shared x.
Highlight (PANEL_HI): in each panel the niche discussed in the text (CN4: Antigen Presentation, Co-Stimulatory and TNF; CN10: IL-6, IL-1, OSM, IL-11 and
Matrix-Integrin) gets a red (#b2182b), slightly larger point (ms 9.5 vs 7.5; filled / open = q; CI bar unchanged) and a red row label. Row labels are drawn
in the left column and shared by the row's two panels; both panels of a row highlight the same niche (asserted) except the last row (TNF: CN4;
Matrix-Integrin: CN10), whose red label follows the left (TNF) panel. HILITE_CN10 (light red CN10 band, red CN10 label, 'CN10' key entry) and EMPH_CN
(bold row labels) are off. The CSV keeps a 'highlighted' column (CN10).
Redraw only (nothing recomputed): estimates / CIs / P / q / stars from fig4_cn_signaling_heatmap.csv, asserted equal to the source summary rows of
cn_lr_modules_costim (modules: fig4_lr_modules_allcells_costim_dotplot.csv) and cn_inflammatory_axes (axes: supp_inflammatory_axes_by_neighbourhood.csv);
section counts from the source section rows (asserted equal to the heatmap's n_sections and to the section means).
Asserts: Wald CI symmetric and reproducing the P value; CI excludes 0 <=> P < 0.05; every starred CI excludes 0; estimate / heatmap mean within-section z
~ 10/9; drawn artists carry exactly these values; layout checks (star column right of all data, same position in all panels; no text overlaps; no bold).
Outputs: /share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_signaling_errorbars.{png,csv} + _legend.txt"""
    import os, sys, warnings; warnings.filterwarnings('ignore')
    # [env check: see README "Environment"] assert os.environ.get('OPENBLAS_CORETYPE') == 'SkylakeX'
    # [env check: see README "Environment"] assert sys.executable.startswith('/data/ydn4687/geneformer_env/')
    import numpy as np, pandas as pd, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    from matplotlib.transforms import Bbox, blended_transform_factory
    from scipy import stats as sst
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'; OUTA = f'{D}/panels/fig4_and_related/A_fig4_panels'
    C82D = f'{OUTA}/fig4_lr_modules_allcells_costim_dotplot.csv'; C86 = f'{OUTA}/supp_inflammatory_axes_by_neighbourhood.csv'
    C89B = f'{OUTA}/fig4_cn_signaling_heatmap.csv'; P89B = f'{OUTA}/fig4_cn_signaling_heatmap.png'
    STEM = '/share/fsmresfiles/UC/revision_jenny/SuppFig9A-F_niche_cxcl12_apc/supp_cn_signaling_errorbars'; OUTFILES = [f'{STEM}.png', f'{STEM}.csv', f'{STEM}_legend.txt']
    assert not set(OUTFILES) & {C82D, C86, C89B, P89B} and os.path.exists(P89B)
    CNS = list(range(1, 11)); MIN_CELLS = 100; HI = '#b2182b'; BAND = '#f9dede'; GREY = '#6e6e6e'; PTC = '#252525'
    # per-panel highlight of the niche discussed in the text (PANEL_HI); CN10 band, bold labels and CN10 key entry off
    EMPH_CN = (); HILITE_CN10 = False
    PANEL_HI = {'Antigen Presentation': 4, 'Co-Stimulatory': 4, 'IL-6': 10, 'IL-1': 10, 'OSM': 10, 'IL-11': 10, 'TNF': 4, 'Matrix–Integrin': 10}
    stars = lambda q: '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''
    sig3 = lambda s: s.map(lambda v: float(f'{v:.3g}') if isinstance(v, (float, np.floating)) and np.isfinite(v) else v)
    GROUPS = [('T-Cell Signals', ['Antigen Presentation', 'Co-Stimulatory']), ('Inflammatory Cytokines', ['IL-6', 'IL-1', 'OSM', 'IL-11', 'TNF']),
              ('Stroma', ['Matrix–Integrin'])]
    COLS = [c for _, cs in GROUPS for c in cs]; GOF = {c: g for g, cs in GROUPS for c in cs}
    MODS = ('Antigen Presentation', 'Co-Stimulatory', 'Matrix–Integrin')
    LMMC = ['lmm_estimate', 'lmm_ci_low', 'lmm_ci_high', 'lmm_p', 'lmm_q']
    # ================= data
    H = pd.read_csv(C89B); H['CN'] = H.CN.astype(int)
    assert len(H) == 80 and list(dict.fromkeys(H.column)) == COLS and not H.duplicated(['column', 'CN']).any()
    assert all(H[H.column == c].group.iloc[0] == GOF[c] for c in COLS)
    M = pd.read_csv(C82D, low_memory=False); X = pd.read_csv(C86, low_memory=False)
    LAB = H.drop_duplicates('CN').set_index('CN').CN_label.loc[CNS].tolist()
    assert all(l.startswith(f'CN{c}: ') for c, l in zip(CNS, LAB)) and LAB[9] == 'CN10: OGN⁺RSPO3⁺/PLA2G2A⁺ Fib'
    SUMR = []
    for c in COLS:
        src, col = (M, 'module') if c in MODS else (X, 'axis')
        srcname = 'fig4_lr_modules_allcells_costim_dotplot.csv (cn_lr_modules_costim)' if c in MODS else 'supp_inflammatory_axes_by_neighbourhood.csv (cn_inflammatory_axes)'
        sm = src[(src.row_type == 'summary') & (src[col] == c)].assign(CN=lambda d: d.CN.astype(int)).set_index('CN').loc[CNS]
        hc = H[H.column == c].set_index('CN').loc[CNS]
        for k_ in ['mean_within_section_z', 'n_sections'] + LMMC:
            assert np.array_equal(sm[k_].values.astype(float), hc[k_].values.astype(float)), (c, k_)      # heatmap == source summary rows
        if c == 'Co-Stimulatory': assert (sm.n_pairs == 7).all() and (hc.n_pairs == 7).all()
        ss = src[(src.row_type == 'section') & (src[col] == c)].assign(CN=lambda d: d.CN.astype(int))
        assert not ss.duplicated(['section', 'CN']).any() and (ss.n_cells >= MIN_CELLS).all(), c
        for cn in CNS:
            h = hc.loc[cn]; v = ss[ss.CN == cn].within_section_z.dropna().values.astype(float)
            st = h.stars if isinstance(h.stars, str) else ''
            assert len(v) == int(h.n_sections) >= 3 and abs(v.mean() - h.mean_within_section_z) < 0.01, (c, cn)        # section counts (legend)
            assert st == stars(h.lmm_q) == h.cell_text[len(h.cell_text.rstrip('*')):], (c, cn)
            # --- the plotted estimate and CI: Wald, symmetric, reproduce P; excludes 0 <=> P < 0.05; every starred CI excludes 0
            e, lo, hi, p, q = (float(h[k_]) for k_ in LMMC)
            assert lo < e < hi and abs((hi - e) - (e - lo)) < 0.012 * max(1.0, abs(e)), (c, cn, e, lo, hi)
            zc_, zp_ = abs(e / ((hi - lo) / 3.92)), sst.norm.isf(p / 2)
            assert abs(zc_ - zp_) <= 0.015 * zp_ + 0.02, (c, cn, zc_, zp_)          # Wald z from the CI == z from P (3-s.f. inputs)
            excl = bool(lo > 0 or hi < 0); assert excl == (p < 0.05), (c, cn, excl, p)
            if st: assert excl, ('STARRED CI INCLUDES 0', c, cn, lo, hi)
            ratio = e / h.mean_within_section_z if abs(h.mean_within_section_z) > 0.01 else np.nan
            if np.isfinite(ratio): assert 1.09 < ratio < 1.125, (c, cn, ratio)
            SUMR.append(dict(group=GOF[c], column=c, CN=cn, CN_label=LAB[cn - 1], CN_name=h.CN_name, n_sections=int(h.n_sections),
                             lmm_estimate=e, lmm_ci_low=lo, lmm_ci_high=hi, lmm_p=p, lmm_q=q, stars=st, filled_point=bool(q < 0.05), ci_excludes_0=excl,
                             heatmap_mean_within_section_z=h.mean_within_section_z, estimate_to_heatmap_ratio=ratio, heatmap_cell_text=h.cell_text,
                             source_csv='fig4_cn_signaling_heatmap.csv (cn_signalling_table) = ' + srcname))
    SUMR = pd.DataFrame(SUMR)
    NSEC = SUMR.groupby('column').n_sections.agg(['min', 'max']).loc[COLS]
    RAT = SUMR.estimate_to_heatmap_ratio.dropna()
    NOSTAR_EXCL = SUMR[(SUMR.stars == '') & SUMR.ci_excludes_0]
    assert not len(SUMR[(SUMR.stars != '') & ~SUMR.ci_excludes_0])
    print('all 80 cells: estimate / CI / q == heatmap == source CSVs; every starred CI excludes 0 (', int((SUMR.stars != '').sum()), 'starred )')
    print('estimate / heatmap mean ratio: %.3f-%.3f' % (RAT.min(), RAT.max()))
    print('CI excludes 0 without a star (P < 0.05, q >= 0.05):', [(r.column, r.CN) for r in NOSTAR_EXCL.itertuples()], flush=True)
    # ================= layout (inches)
    FS_LAB, FS_TICK, FS_TITLE, FS_GROUP, FS_STAR, FS_KEY = 13, 13, 14, 13, 14, 13
    plt.rcParams.update({'font.size': FS_TICK, 'font.family': 'DejaVu Sans', 'font.weight': 'normal', 'axes.titleweight': 'normal', 'axes.labelweight': 'normal'})
    PAL = dict(zip(CNS, sns.color_palette('tab10', 10).as_hex()))
    fig0 = plt.figure(); R0 = fig0.canvas.get_renderer()
    tw = lambda s, fs: R0.get_text_width_height_descent(s, FontProperties(family='DejaVu Sans', size=fs), ismath=False)[0] / fig0.dpi; plt.close(fig0)
    SQX, SQW, PADL = 0.05, 0.13, 0.07
    L0 = max(tw(l, FS_LAB) for l in LAB) + SQX + SQW + PADL + 0.04
    # L0 (panel positions / widths) is set from the regular-weight labels; a bold EMPH_CN label can reach left of the canvas edge, and the tight
    # bounding box at save then extends the PNG to the left to include it (asserted at save).
    FW = 11.0; RM = 0.06; GX = 0.30; PW = (FW - L0 - RM - GX) / 2
    PH, TH, BOT, TOP = 2.43, 0.40, 0.72, 0.42; FH = BOT + 4 * (PH + TH) + TOP
    fig = plt.figure(figsize=(FW, FH))
    AX = {}
    for k, c in enumerate(COLS):
        r_, c_ = divmod(k, 2)
        AX[c] = fig.add_axes([(L0 + c_ * (PW + GX)) / FW, (BOT + (3 - r_) * (PH + TH)) / FH, PW / FW, PH / FH])
    MS, MEW, ELW = 7.5, 1.6, 1.9
    MS_HI = 9.5                                                          # highlighted point, slightly larger
    R_PT = max(MS, MS_HI if PANEL_HI else MS) / 2 + MEW / 2; R_CI = ELW / 2; GAP, MARG = 6.0, 3.0      # pt
    XLO = SUMR.lmm_ci_low.min() - 0.15
    FP = FontProperties(family='DejaVu Sans', size=FS_STAR)
    GLY = {s: TextPath((0, 0), s, size=FS_STAR, prop=FP).get_extents() for s in ['*', '**', '***']}
    SW = max(GLY[s].x1 - GLY[s].x0 for s in SUMR.stars if s); W_pt = PW * 72
    def column(x0, x1):
        ppx = W_pt / (x1 - x0)
        right = max((SUMR.lmm_ci_high + R_CI / ppx).max(), (SUMR.lmm_estimate + R_PT / ppx).max())
        xc = right + (GAP + SW / 2) / ppx
        return xc, xc + (SW / 2 + MARG) / ppx
    XHI = SUMR.lmm_ci_high.max() + 0.15
    for it in range(100):
        xc, n1 = column(XLO, XHI)
        if abs(n1 - XHI) < 1e-9: break
        XHI = n1
    assert it < 99
    YLIM = (9.55, -0.55)
    STXT = {c: {} for c in COLS}; SQUARES = {c: [] for c in COLS}; PTS = {c: [] for c in COLS}; BARS = {c: [] for c in COLS}
    for k, c in enumerate(COLS):
        ax = AX[c]; left = k % 2 == 0; bottom = k >= 6
        ax.set_xlim(XLO, XHI); ax.set_ylim(*YLIM)
        if HILITE_CN10: ax.axhspan(8.5, 9.5, color=BAND, lw=0, zorder=0)   # CN10 band only when HILITE_CN10
        ax.axvline(0, ls=(0, (4, 3)), color='#555555', lw=1.0, zorder=1)
        h_ = SUMR[SUMR.column == c].set_index('CN')
        for i, cn in enumerate(CNS):
            r = h_.loc[cn]
            BARS[c].append(ax.hlines(i, r.lmm_ci_low, r.lmm_ci_high, color=PTC, lw=ELW, zorder=3))
            hc_ = PTC if PANEL_HI.get(c) != cn else HI                               # per-panel highlight (red, larger point)
            PTS[c].append(ax.plot(r.lmm_estimate, i, 'o', ms=MS if hc_ == PTC else MS_HI, color=hc_, mfc=hc_ if r.filled_point else 'white', mew=MEW,
                                  zorder=4 if hc_ == PTC else 5)[0])
            if r.stars:
                gb = GLY[r.stars]
                STXT[c][cn] = ax.annotate(r.stars, xy=(xc, i), xytext=(-(gb.x0 + gb.x1) / 2, -(gb.y0 + gb.y1) / 2), textcoords='offset points',
                                          ha='left', va='baseline', fontsize=FS_STAR, color='k', annotation_clip=False)
        ax.set_title(c, loc='left', fontsize=FS_TITLE, pad=4)
        ax.set_title(GOF[c], loc='right', fontsize=FS_GROUP, color=GREY, pad=4)
        ax.set_yticks(range(10)); ax.tick_params(axis='y', length=0)
        if left:
            ax.set_yticklabels(LAB, fontsize=FS_LAB); ax.tick_params(axis='y', pad=(SQX + SQW + PADL) * 72)
            if HILITE_CN10: ax.get_yticklabels()[9].set_color(HI)
            if PANEL_HI:   # shared row labels: red for the niche highlighted in the row's panels
                # only the TNF / Matrix-Integrin row differs (TNF CN4, Matrix-Integrin CN10); the label follows the LEFT panel
                assert PANEL_HI[c] == PANEL_HI[COLS[k + 1]] or (c, COLS[k + 1]) == ('TNF', 'Matrix–Integrin'), ('row panels highlight different niches', c, COLS[k + 1])
                ax.get_yticklabels()[CNS.index(PANEL_HI[c])].set_color(HI)
            for i, cn in enumerate(CNS):
                if cn in EMPH_CN: ax.get_yticklabels()[i].set_fontweight('bold')   # bold EMPH_CN row labels
            tr = blended_transform_factory(ax.transAxes, ax.transData)
            for i, cn in enumerate(CNS):
                SQUARES[c].append(ax.add_patch(Rectangle((-(SQX + SQW) / PW, i - 0.40), SQW / PW, 0.80, transform=tr, color=PAL[cn], clip_on=False, lw=0)))
        else:
            ax.set_yticklabels([])
        ax.set_xticks([-2, -1, 0, 1, 2]); ax.tick_params(axis='x', labelsize=FS_TICK, length=3)
        if not bottom: ax.set_xticklabels([])
        ax.spines[['top', 'right', 'left']].set_visible(False)
    fig.canvas.draw(); rd = fig.canvas.get_renderer()
    xt_bottom = min(t.get_window_extent(rd).y0 for c in COLS[6:] for t in AX[c].get_xticklabels() if t.get_text())
    XLAB = fig.text((L0 + PW + GX / 2) / FW, (xt_bottom - 6) / fig.dpi / FH, 'Niche vs Other Niches (Mixed-Model Estimate, z)', ha='center', va='top', fontsize=FS_TITLE)
    key = [Line2D([], [], marker='o', ls='-', color=PTC, mfc=PTC, ms=MS, mew=MEW, lw=ELW, label='Estimate ± 95% CI, q < 0.05'),
           Line2D([], [], marker='o', ls='-', color=PTC, mfc='white', ms=MS, mew=MEW, lw=ELW, label='q ≥ 0.05'),
           Line2D([], [], ls=(0, (4, 3)), color='#555555', lw=1.0, label='0 = No Difference'),
           Rectangle((0, 0), 1, 1, color=BAND, lw=0, label='CN10')][:4 if HILITE_CN10 else 3]   # no CN10 key entry without the band
    KEY = fig.legend(handles=key, loc='lower right', bbox_to_anchor=((L0 + 2 * PW + GX) / FW, (FH - TOP + 0.02) / FH), ncol=4, frameon=False,
                     fontsize=FS_KEY, handletextpad=0.4, columnspacing=1.2, borderaxespad=0, handlelength=1.8)
    # ================= checks
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); px = fig.dpi / 72
    ov = lambda a, b: a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1
    TXT = {'key': KEY.get_window_extent(rd), 'xlabel': XLAB.get_window_extent(rd)}
    xcs = []
    for k, c in enumerate(COLS):
        ax = AX[c]; AB = ax.get_window_extent(rd); T = ax.transData.transform; h_ = SUMR[SUMR.column == c].set_index('CN')
        for i, cn in enumerate(CNS):                                            # drawn values == heatmap / source values
            r = h_.loc[cn]; ln = PTS[c][i]; sg = BARS[c][i].get_segments()[0]
            assert ln.get_xdata()[0] == r.lmm_estimate and ln.get_ydata()[0] == i and (ln.get_markerfacecolor() != 'white') == r.filled_point, (c, cn)
            assert (ln.get_markeredgecolor() == HI and ln.get_markersize() == MS_HI) == (PANEL_HI.get(c) == cn), (c, cn)   # highlight check
            assert sg[0][0] == r.lmm_ci_low and sg[1][0] == r.lmm_ci_high and sg[0][1] == sg[1][1] == i, (c, cn)
        assert {cn: t.get_text() for cn, t in STXT[c].items()} == {cn: h_.loc[cn, 'stars'] for cn in CNS if h_.loc[cn, 'stars']}, c
        data_right = max(max(T((h_.loc[cn, 'lmm_ci_high'], 0))[0] + R_CI * px, T((h_.loc[cn, 'lmm_estimate'], 0))[0] + R_PT * px) for cn in CNS)
        for cn, t in STXT[c].items():
            gb = GLY[t.get_text()]; cx, cy = T((xc, cn - 1)); hw = (gb.x1 - gb.x0) / 2 * px; hh = (gb.y1 - gb.y0) / 2 * px
            G_ = Bbox([[cx - hw, cy - hh], [cx + hw, cy + hh]]); xcs.append(cx - AB.x0)
            tb = t.get_window_extent(rd); assert abs(tb.x0 - (cx - (gb.x0 + gb.x1) / 2 * px)) < 1.0, (c, cn)
            assert G_.x0 > data_right and G_.x1 <= AB.x1 - 1 and G_.y0 >= AB.y0 and G_.y1 <= AB.y1, ('star placement', c, cn)
            TXT[f'star {c} CN{cn}'] = G_
        assert ax.title.get_text() == '' and ax._left_title.get_text() == c
        TXT[f'title {c}'] = ax._left_title.get_window_extent(rd); TXT[f'group {c}'] = ax._right_title.get_window_extent(rd)
        assert TXT[f'title {c}'].x1 + 6 < TXT[f'group {c}'].x0, ('title and group tag overlap', c)
        for i, t in enumerate(ax.get_yticklabels()):
            if t.get_text(): TXT[f'ylab {c} {i}'] = t.get_window_extent(rd)
        for i, t in enumerate(ax.get_xticklabels()):
            if t.get_text(): TXT[f'xtick {c} {i}'] = t.get_window_extent(rd)
        if k % 2 == 0:
            assert [t.get_text() for t in ax.get_yticklabels()] == LAB and len(SQUARES[c]) == 10
            assert [cn for cn, t in zip(CNS, ax.get_yticklabels()) if t.get_color() == HI] == ([10] if HILITE_CN10 else []) + ([PANEL_HI[c]] if PANEL_HI and PANEL_HI[c] != 10 or (PANEL_HI and not HILITE_CN10) else [])
            assert len(ax.patches) == 10 + int(HILITE_CN10)   # 10 colour squares (+ the CN10 band only if highlighted)
            for i, p in enumerate(SQUARES[c]): TXT[f'square {c} {i}'] = p.get_window_extent(rd)
    assert max(xcs) - min(xcs) < 0.5, 'star column differs between panels'
    ks = list(TXT)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            assert not ov(TXT[ks[i]], TXT[ks[j]]), ('overlap', ks[i], ks[j])
    AXB = {c: AX[c].get_window_extent(rd) for c in COLS}
    for k, c in enumerate(COLS):
        for t in (f'title {c}', f'group {c}'):
            assert TXT[t].y0 >= AXB[c].y1 - 0.5 and (k < 2 or TXT[t].y1 < AXB[COLS[k - 2]].y0), (t,)
        for key_ in [kk for kk in ks if kk.split(' ')[0] in ('key', 'xlabel', 'title', 'group')]:
            if key_.endswith(' ' + c): continue
            assert not ov(TXT[key_], AXB[c]), ('text on axes', key_, c)
    assert TXT['key'].x0 >= min(TXT[k_].x0 for k_ in ks if k_.startswith('ylab')) and TXT['key'].x1 <= AXB[COLS[1]].x1 + 1, 'key wider than the figure body'
    assert TXT['key'].y0 > max(TXT[f'title {c}'].y1 for c in COLS[:2]) and TXT['key'].y0 > max(TXT[f'group {c}'].y1 for c in COLS[:2])
    assert TXT['xlabel'].y1 < min(TXT[k_].y0 for k_ in ks if k_.startswith('xtick'))
    BOLDT = [AX[c].get_yticklabels()[CNS.index(cn)] for c in COLS[::2] for cn in EMPH_CN]   # only the EMPH_CN row labels are bold
    assert all(t.get_fontweight() == 'bold' for t in BOLDT) and [t.get_text()[:4] for t in BOLDT] == [f'CN{cn}:' for c in COLS[::2] for cn in EMPH_CN]
    assert all(t.get_fontweight() in ('normal', 400) for a in AX.values() for t in a.texts + [a._left_title, a._right_title] + a.get_yticklabels() + a.get_xticklabels() if t not in BOLDT)
    assert len(KEY.get_texts()) == 3 + int(HILITE_CN10)
    print('checks passed: points / CI bars / stars == heatmap values; star column right of all data at the same position in all panels; no overlaps', flush=True)
    if BOLDT: assert fig.get_tightbbox(rd).x0 * fig.dpi <= min(t.get_window_extent(rd).x0 for t in BOLDT) + 0.5   # bold labels inside the saved area
    fig.savefig(OUTFILES[0], dpi=300, bbox_inches='tight', pad_inches=0.05)
    tb = fig.get_tightbbox(rd); print(f'figure size (tight): {tb.width:.2f} x {tb.height:.2f} in', flush=True); plt.close(fig)
    # ================= source data
    OUTC = SUMR.copy(); OUTC['highlighted'] = OUTC.CN == 10
    OUTC = OUTC.assign(metric='per section: mean LIANA bivariate local score over the column\'s pairs and over ALL cells of the niche (>= 100 cells), z across the '
                          '10 niches within the section; PLOTTED: lmm_estimate (niche minus the other niches, within-section z) with its Wald 95% CI',
                       model='z ~ in_CN + (1 | section) on per-section z; BH over the 10 niches within each column (values from the heatmap; not recomputed)')
    for c_ in ['CN', 'n_sections']: OUTC[c_] = OUTC[c_].round().astype('Int64')
    OUTC.apply(sig3).to_csv(OUTFILES[1], index=False)
    # ================= legend (US spelling)
    nmax = NSEC['max'].astype(int); nmin = NSEC['min'].astype(int)
    fewer = [f'{c} {nmin[c]}–{nmax[c]}' for c in COLS if nmin[c] < nmax[c]]
    ne = '; '.join(f'{r.column} CN{r.CN}' for r in NOSTAR_EXCL.itertuples())
    leg = ('Niche-versus-other-niche differences behind the niche signaling heatmap (Fig 5C; fig4_cn_signaling_heatmap.png), CD resections (CosMx 6K; '
       '11 sections, 4 patients; 2,936,007 cells assigned to a niche). Panels are the heatmap\'s columns, in its order and groups (group named at the '
       'top right of each panel): T-Cell Signals, Antigen Presentation (CellChatDB MHC-II, 9 pairs) and Co-Stimulatory (canonical antigen-presenting '
       'cell → T-cell co-stimulation, 7 pairs: CD80 → CD28; CD86 → CD28; ICOSLG → ICOS; TNFSF4 → TNFRSF4; TNFSF18 → TNFRSF18; TNFSF9 → TNFRSF9; '
       'CD70 → CD27); Inflammatory Cytokines, the pre-specified axes IL-6 (IL6 → IL6R, IL6ST), IL-1 (IL1A, IL1B → IL1R1, IL1RAP), OSM (OSM → OSMR, '
       'LIFR, IL6ST), IL-11 (IL11 → IL11RA, IL6ST) and TNF (TNF → TNFRSF1A, TNFRSF1B; contrast); Stroma, Matrix–Integrin (CellChatDB ECM-Receptor '
       'pathways, 227 pairs). Scores are cell-type-agnostic spatial LIANA bivariate local scores (Gaussian kernel, 200 px bandwidth) for the NicheNet '
       'v2 prior-network pairs. Rows are the 10 niches in the heatmap\'s order (CN1 at the top; colored squares as in the heatmap)'
           + ('; CN10, the OGN⁺RSPO3⁺/PLA2G2A⁺ fibroblast niche, is on a light red band with its label in red' if HILITE_CN10 else '')
           + ('; ' + ' and '.join(f'CN{c}' for c in EMPH_CN) + ' row labels are shown in bold' if EMPH_CN else '')
           + ('. Red, CN4 (T-cell signals and TNF) or CN10 (other cytokine and stroma panels)' if PANEL_HI else '')
           + '. For each section, a niche\'s value is the mean score over '
       'all its cells (niche × section combinations with at least 100 cells), z-standardized across the 10 niches within the section (within-section z). '
       'Points: for each niche, the estimate of a linear mixed model fitted to the per-section values of all niches, z ~ in niche + (1 | section), '
       'i.e. the difference in within-section z between that niche and the other nine niches (dashed line: 0, no difference); filled points, q < 0.05; '
       'open points, q ≥ 0.05. Because the within-section z values of the 10 niches sum to zero in each section, the other nine niches average −z/9 '
       'when a niche has value z, so the estimate is about 10/9 of the heatmap\'s mean within-section z '
       f'(observed ratio {RAT.min():.2f}–{RAT.max():.2f}; exactly 10/9 only when every niche has a value in every section). Bars: 95% CI of the '
       'estimate (Wald, ±1.96 SE), from the same model that gives the stars. Number of sections per panel: '
           + ', '.join(f'{c} {nmax[c]}' for c in COLS) + ' (sections in which a column shows no variation across niches are left out'
           + (f'; a niche with fewer than 100 cells in a section has no value there, so some niches have fewer sections: {", ".join(fewer)}' if fewer else '')
           + '). Stars, as in the heatmap: q values of this model, Benjamini–Hochberg-corrected over the 10 niches within each column; * q < 0.05, '
       '** q < 0.01, *** q < 0.001; no star, q ≥ 0.05. Every starred CI excludes 0; the CIs are not adjusted for multiplicity, so a CI can exclude 0 '
       f'without a star ({ne}; nominal P < 0.05, q ≥ 0.05). Source data (supp_cn_signaling_errorbars.csv) hold the estimates, CIs, P and q values '
       'and the heatmap\'s mean within-section z for every niche and column.')
    for t in ('activation score', 'quiescent', 'territor', 'serial', 'blunt', 'colour', 'neighbour', 'signalling', 'characteris', 'standardis'):
        assert t not in leg.lower(), t
    f_ = open(OUTFILES[2], 'w'); f_.write(leg + '\n'); f_.close()
    print('wrote', *OUTFILES, sep='\n  '); print('done', flush=True)
