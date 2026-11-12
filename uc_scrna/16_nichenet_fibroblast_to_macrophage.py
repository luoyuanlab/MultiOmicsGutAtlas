#!/usr/bin/env python
"""UC scRNA-seq atlas: NicheNet ligands from OGN+RSPO3+ fibroblasts to SELENOP+MMP9+SPP1+ macrophages (patient pseudobulk).

Panels made
  Supp 5C  NicheNet ligand activity of OGN+RSPO3+ fibroblast ligands for genes up in Mac S+M+S+ vs resident SELENOP+
           macrophages (patient pseudobulk DESeq2), fibroblast enrichment and a patient-paired target-response test
  Supp tables: NicheNet ligands, ligand rank stability

Inputs
  /share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad (UC atlas)
  /share/fsmresfiles/ydn4687/nichenet/ (NicheNet ligand_target_matrix_nsga2r_final.rds, lr_network_human_21122021.rds)

Outputs
  /share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands/ (intermediate/, CSVs, panel)

Run order of sections
  S5C_pseudobulk_de -> S5C_geneset_size_sensitivity -> S5C_nichenet_ranking -> S5C_patient_bootstrap
  -> S5C_ligand_response_test -> S5C_ligand_activity_plot

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist.
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
# Section S5C_pseudobulk_de: patient pseudobulk DESeq2 of Mac S+M+S+ (SMS) vs pooled resident SELENOP+ macrophages
# Panels: Supp Fig 5C (receiver program)
# ======================================================================================================================
if _run('S5C_pseudobulk_de'):
    """Receiver program for NicheNet: patient-level pseudobulk DE (pydeseq2) of SELENOP+MMP9+SPP1+ (SMS) macrophages vs the
pooled homeostatic resident SELENOP+ macrophages (Mac S+XS- + Mac S+SG+). Raw counts are summed per patient (obs['Patient'])
and group; pseudobulks with >= 10 cells and genes with >= 10 counts are kept. Outputs: pseudobulk_DE_homeoSp_pt.csv,
intermediate/geneset_up_in_SMS_vs_homeoSp_pt.csv (padj < 0.1, log2FC > 0.5)."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, pyreadr
    from sklearn.metrics import average_precision_score, roc_auc_score
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands'
    NN='/share/fsmresfiles/ydn4687/nichenet'
    REF=['Mac S+XS-','Mac S+SG+']; TAG='homeoSp_pt'   # pooled homeostatic SELENOP+ macs
    scrna=sc.read_h5ad('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad')
    mye=scrna[scrna.obs['cell_category']=='Myeloid'].copy()

    # ---- patient-level pseudobulk DE: SMS vs pooled homeostatic SELENOP+ macs ----
    from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
    def pseudobulk_de(target, ref, tag):
        refs=[ref] if isinstance(ref,str) else list(ref)
        sub=mye[mye.obs['cell_type_short'].isin([target]+refs)].copy()
        raw=sub.raw.to_adata()[sub.obs_names]; genes=sub.raw.var['features'].astype(str).values
        M=raw.X; M=M.toarray() if sp.issparse(M) else np.asarray(M)
        cnt=pd.DataFrame(M, index=sub.obs_names, columns=genes); cnt=cnt.loc[:,~cnt.columns.duplicated()]
        meta=sub.obs[['Patient','cell_type_short']].copy(); meta['grp']=np.where(meta['cell_type_short']==target,'target','ref')
        cnt['pb']=meta['Patient'].astype(str)+'__'+meta['grp']; pb=cnt.groupby('pb').sum()
        pbmeta=pd.DataFrame({'grp':[i.split('__')[1] for i in pb.index],'patient':[i.split('__')[0] for i in pb.index]}, index=pb.index)
        ncell=meta.groupby(meta['Patient'].astype(str)+'__'+meta['grp']).size()
        keep=ncell[ncell>=10].index; pb=pb.loc[pb.index.isin(keep)]; pbmeta=pbmeta.loc[pb.index]
        print(f'  [{tag}] pseudobulk samples: {pbmeta.grp.value_counts().to_dict()}')
        pb=pb.loc[:, pb.sum(0)>=10].round().astype(int)
        dds=DeseqDataSet(counts=pb, metadata=pbmeta, design_factors='grp', ref_level=['grp','ref'], quiet=True)
        dds.deseq2(); st=DeseqStats(dds, contrast=['grp','target','ref'], quiet=True); st.summary()
        res=st.results_df.sort_values('padj'); res.to_csv(f'{O}/pseudobulk_DE_{tag}.csv')
        up=res[(res.padj<0.1)&(res.log2FoldChange>0.5)].index.tolist()
        pd.Series(up).to_csv(f'{O}/intermediate/geneset_up_in_SMS_vs_{tag}.csv', index=False, header=['gene'])
        print(f'  [{tag}] genes UP in SMS (padj<0.1, LFC>0.5): {len(up)}'); return up

    up=pseudobulk_de('Mac S+M+S+', REF, TAG)
    print('SAVED')


# ======================================================================================================================
# Section S5C_geneset_size_sensitivity: NicheNet ligand ranks for four definitions of the SMS target gene set
# Panels: ligand rank stability supplementary table (input)
# ======================================================================================================================
if _run('S5C_geneset_size_sensitivity'):
    """NicheNet ligand ranking sensitivity to the target gene-set size, using the patient-level pseudobulk DE
(pseudobulk_DE_homeoSp_pt.csv; SMS mac vs pooled resident SELENOP+ macs; 20 vs 54 patient pseudobulks). The default threshold
(padj<0.1, LFC>0.5) gives 1456 genes, which flattens the NicheNet AUPR, so stricter sets are also scored: padj<0.05 & LFC>1,
and the top 300 / top 100 genes by LFC (padj<0.05). Ligands: expressed (> 10% of cells) by OGN+RSPO3+ fibroblasts with a
receptor expressed by SMS macrophages. Output: intermediate/intermediate_ligand_rank_by_geneset_size.csv."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, pyreadr
    from sklearn.metrics import average_precision_score, roc_auc_score
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands'; NN='/share/fsmresfiles/ydn4687/nichenet'
    de=pd.read_csv(f'{O}/pseudobulk_DE_homeoSp_pt.csv',index_col=0).dropna(subset=['padj'])
    scrna=sc.read_h5ad('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad',backed='r')
    def expressed(ct,thr=0.10):
        idx=np.flatnonzero((scrna.obs['cell_type_short']==ct).values); X=scrna[idx].to_memory().X
        return set(scrna.var_names[np.asarray((X>0).mean(0)).flatten()>thr])
    fib_expr=expressed('OGN+RSPO3+ Fib'); mac_expr=expressed('Mac S+M+S+')
    ltm=list(pyreadr.read_r(f'{NN}/ligand_target_matrix_nsga2r_final.rds').values())[0]
    lr=list(pyreadr.read_r(f'{NN}/lr_network_human_21122021.rds').values())[0]; lr.columns=[c.lower() for c in lr.columns]
    lr_ok=lr[lr['from'].isin(fib_expr)&lr['to'].isin(mac_expr)]; ligs=sorted(set(lr_ok['from'])&set(ltm.columns))
    bg=sorted(set(ltm.index)&mac_expr)
    up=de[de.log2FoldChange>0]
    SETS={'padj<0.1 & LFC>0.5 (default)':up[(up.padj<0.1)&(up.log2FoldChange>0.5)].index,
          'padj<0.05 & LFC>1':up[(up.padj<0.05)&(up.log2FoldChange>1)].index,
          'top 300 by LFC (padj<0.05)':up[up.padj<0.05].sort_values('log2FoldChange',ascending=False).index[:300],
          'top 100 by LFC (padj<0.05)':up[up.padj<0.05].sort_values('log2FoldChange',ascending=False).index[:100]}
    out=[]
    for name,gs in SETS.items():
        gs=set(gs)&set(bg); resp=np.array([g in gs for g in bg],int)
        rows=[]
        for L in ligs:
            tv=ltm.loc[bg,L].values.astype(float)
            if len(np.unique(tv))<2: continue
            rows.append([L,average_precision_score(resp,tv),roc_auc_score(resp,tv)])
        a=pd.DataFrame(rows,columns=['ligand','aupr','auroc']).sort_values('aupr',ascending=False).reset_index(drop=True)
        a['rank']=a.index+1; a['set']=name; a['n_genes']=len(gs); a['aupr_baseline']=resp.mean(); out.append(a)
        spread=a['aupr'].iloc[0]/a['aupr'].iloc[19]
        print(f"\n== {name}: {len(gs)} genes in background (chance AUPR {resp.mean():.3f}); top1/top20 AUPR ratio {spread:.2f}")
        print(a.head(12)[['rank','ligand','aupr','auroc']].round(3).to_string(index=False))
    A=pd.concat(out); A.to_csv(f'{O}/intermediate/intermediate_ligand_rank_by_geneset_size.csv',index=False)
    key=['CSF1','TGFB1','IL33','TNFSF12','FGF7','FN1','ADAM17','CCL2','TNC','COL15A1']
    print('\nRank of key ligands per gene-set definition:')
    print(A[A.ligand.isin(key)].pivot(index='ligand',columns='set',values='rank').to_string())
    print('SAVED')


# ======================================================================================================================
# Section S5C_nichenet_ranking: NicheNet ligand activity (AUPR, AUROC) of OGN+RSPO3+ fibroblast ligands for the SMS program
# Panels: Supp Fig 5C (ranking), ligand rank stability supplementary table
# ======================================================================================================================
if _run('S5C_nichenet_ranking'):
    """NicheNet ranking with the patient as replicate. Receiver program = genes up in SELENOP+MMP9+SPP1+ macrophages vs pooled
resident SELENOP+ macrophages (Mac S+XS- + Mac S+SG+) in the patient pseudobulk DE (pseudobulk_DE_homeoSp_pt.csv, 20 vs 54
patient pseudobulks), padj<0.05 & log2FC>1 (494 genes in the NicheNet background of SMS-expressed genes). Ligand activity =
AUPR (and AUROC) of the ligand's regulatory potential for the program; also records the mean expression in OGN+RSPO3+ vs
other fibroblasts and the matching SMS receptors.
Outputs: nichenet_ligand_activity_patient.csv (full ranking), supp_table_ligand_rank_stability.csv (top 15 ligands, ranks
across the gene-set definitions of S5C_geneset_size_sensitivity)."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, pyreadr, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from sklearn.metrics import average_precision_score, roc_auc_score
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands'; NN='/share/fsmresfiles/ydn4687/nichenet'
    de=pd.read_csv(f'{O}/pseudobulk_DE_homeoSp_pt.csv',index_col=0).dropna(subset=['padj'])
    geneset=set(de[(de.padj<0.05)&(de.log2FoldChange>1)].index)
    scrna=sc.read_h5ad('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad',backed='r')
    ct=scrna.obs['cell_type_short']
    def sub(cts): return scrna[np.flatnonzero(ct.isin(cts).values)].to_memory()
    def expressed(name,thr=0.10):
        X=sub([name]).X; return set(scrna.var_names[np.asarray((X>0).mean(0)).flatten()>thr])
    fib_expr=expressed('OGN+RSPO3+ Fib'); mac_expr=expressed('Mac S+M+S+')
    ltm=list(pyreadr.read_r(f'{NN}/ligand_target_matrix_nsga2r_final.rds').values())[0]
    lr=list(pyreadr.read_r(f'{NN}/lr_network_human_21122021.rds').values())[0]; lr.columns=[c.lower() for c in lr.columns]
    lr_ok=lr[lr['from'].isin(fib_expr)&lr['to'].isin(mac_expr)]; ligs=sorted(set(lr_ok['from'])&set(ltm.columns))
    bg=sorted(set(ltm.index)&mac_expr); resp=np.array([g in geneset for g in bg],int)
    print(f'program genes in background: {resp.sum()} / {len(bg)}; potential ligands: {len(ligs)}')
    rows=[]
    for L in ligs:
        tv=ltm.loc[bg,L].values.astype(float)
        if len(np.unique(tv))<2: continue
        rows.append([L,average_precision_score(resp,tv),roc_auc_score(resp,tv)])
    act=pd.DataFrame(rows,columns=['ligand','aupr','auroc']).sort_values('aupr',ascending=False).reset_index(drop=True)
    FIB=['ADAMDEC1+ Fib','VSTM2A+ Crypt Top Fib','OGN+RSPO3+ Fib','NRG1+ Crypt Top Fib','Activ Fib','CXCL5+ Activ Fib','T-interact Fib','SELENOP+ Fib']
    fibs=sub(FIB); isogn=(fibs.obs['cell_type_short']=='OGN+RSPO3+ Fib').values
    def mexpr(mask,g): return float(np.asarray(fibs[mask,g].X.mean()).flatten()[0])
    act['ogn_mean']=[mexpr(isogn,g) for g in act.ligand]; act['otherfib_mean']=[mexpr(~isogn,g) for g in act.ligand]
    act['fib_enriched']=act['ogn_mean']>act['otherfib_mean']
    act['mac_receptors']=[','.join(sorted(set(lr_ok[lr_ok['from']==g]['to'])&mac_expr)) for g in act.ligand]
    act['rank']=act.index+1
    act.to_csv(f'{O}/nichenet_ligand_activity_patient.csv',index=False)
    s3=lambda v: f'{v:#.3g}'.rstrip('.')
    top=act.head(15)
    st=pd.read_csv(f'{O}/intermediate/intermediate_ligand_rank_by_geneset_size.csv')
    keep=top['ligand'].tolist()
    piv=st[st.ligand.isin(keep)].pivot(index='ligand',columns='set',values='rank').reindex(keep)
    order=['padj<0.05 & LFC>1','padj<0.1 & LFC>0.5 (default)','top 300 by LFC (padj<0.05)','top 100 by LFC (padj<0.05)']
    nm={'padj<0.05 & LFC>1':'Rank, padj<0.05 and LFC>1 (main)','padj<0.1 & LFC>0.5 (default)':'Rank, padj<0.1 and LFC>0.5',
        'top 300 by LFC (padj<0.05)':'Rank, top 300 genes','top 100 by LFC (padj<0.05)':'Rank, top 100 genes'}
    piv=piv[order].rename(columns=nm).reset_index().rename(columns={'ligand':'Ligand'})
    piv.to_csv(f'{O}/supp_table_ligand_rank_stability.csv',index=False)
    print(top[['rank','ligand','aupr','auroc','fib_enriched']].round(3).to_string(index=False)); print(piv.to_string(index=False))
    print('SAVED')


# ======================================================================================================================
# Section S5C_patient_bootstrap: patient-bootstrap stability of the NicheNet ligand ranks (100 resamples)
# Panels: Supp Fig 5C (rank stability)
# ======================================================================================================================
if _run('S5C_patient_bootstrap'):
    """Patient-bootstrap stability of the NicheNet ligand ranking. Per iteration: resample the 54 patients with replacement
(duplicates = separate pseudobulk units), redo the pseudobulk DE (pydeseq2, SMS mac vs pooled resident SELENOP+ macs, >= 10
cells per unit), take padj<0.05 & log2FC>1 as the program, recompute the NicheNet AUPR for the same 175 ligands and record the
ranks. NBOOT=100. Outputs: intermediate/intermediate_patient_pseudobulk_counts.csv.gz (pseudobulk cache),
intermediate/intermediate_bootstrap_ranks.csv (ligand x iteration ranks); adds % of resamples in top 5 / top 10 and the
median rank to nichenet_ligand_activity_patient.csv."""
    import os, warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, pyreadr
    from sklearn.metrics import average_precision_score
    from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands'; NN='/share/fsmresfiles/ydn4687/nichenet'
    NBOOT=100; rng=np.random.default_rng(0); CACHE=f'{O}/intermediate/intermediate_patient_pseudobulk_counts.csv.gz'
    TARGET='Mac S+M+S+'; REF=['Mac S+XS-','Mac S+SG+']
    act=pd.read_csv(f'{O}/nichenet_ligand_activity_patient.csv'); LIGS=act['ligand'].tolist()
    scrna=sc.read_h5ad('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad')
    ob=scrna.obs
    if os.path.exists(CACHE):
        PB=pd.read_csv(CACHE,index_col=0)
    else:
        sub=scrna[ob['cell_type_short'].isin([TARGET]+REF).values].copy()
        raw=sub.raw.to_adata()[sub.obs_names]; genes=sub.raw.var['features'].astype(str).values
        grp=np.where(sub.obs['cell_type_short']==TARGET,'target','ref'); key=sub.obs['Patient'].astype(str).values+'__'+grp
        keys=pd.Series(key); n=keys.value_counts(); keep=n[n>=10].index
        rows={}
        M=raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
        for k in keep:
            rows[k]=np.asarray(M[np.flatnonzero(keys.values==k)].sum(0)).flatten()
        PB=pd.DataFrame(rows,index=genes).T; PB=PB.loc[:,~PB.columns.duplicated()]
        PB.to_csv(CACHE)
    print('pseudobulk units:',PB.shape)
    units=pd.Series(PB.index); pat=units.str.split('__').str[0]; grp=units.str.split('__').str[1]
    patients=sorted(pat.unique()); print('patients:',len(patients),'| target units:',(grp=='target').sum(),'ref units:',(grp=='ref').sum())
    # NicheNet background (SMS-mac-expressed genes in ligand-target matrix), as in S5C_nichenet_ranking
    X=scrna[np.flatnonzero((ob['cell_type_short']==TARGET).values)].X
    mac_expr=set(scrna.var_names[np.asarray((X>0).mean(0)).flatten()>0.10])
    ltm=list(pyreadr.read_r(f'{NN}/ligand_target_matrix_nsga2r_final.rds').values())[0]
    bg=sorted(set(ltm.index)&mac_expr); Mlt=ltm.loc[bg,LIGS].values.astype(float)
    def rank_once(units_idx,tag):
        C=PB.iloc[units_idx].copy(); C.index=[f'{tag}{i}' for i in range(len(C))]
        meta=pd.DataFrame({'grp':grp.iloc[units_idx].values},index=C.index)
        if (meta.grp=='target').sum()<3 or (meta.grp=='ref').sum()<3: return None
        C=C.loc[:,C.sum(0)>=10].round().astype(int)
        dds=DeseqDataSet(counts=C,metadata=meta,design_factors='grp',ref_level=['grp','ref'],quiet=True); dds.deseq2()
        st=DeseqStats(dds,contrast=['grp','target','ref'],quiet=True); st.summary(); r=st.results_df.dropna(subset=['padj'])
        gs=set(r[(r.padj<0.05)&(r.log2FoldChange>1)].index); resp=np.array([g in gs for g in bg],int)
        if resp.sum()<10: return None
        au=np.array([average_precision_score(resp,Mlt[:,j]) for j in range(len(LIGS))])
        return pd.Series(au,index=LIGS).rank(ascending=False,method='min')
    ranks={}
    for b in range(NBOOT):
        pick=rng.choice(patients,size=len(patients),replace=True)
        idx=np.concatenate([np.flatnonzero(pat.values==p) for p in pick])
        r=rank_once(idx,f'b{b}_')
        if r is not None: ranks[b]=r
        if b%10==0: print(f'iter {b}: done {len(ranks)}',flush=True)
    R=pd.DataFrame(ranks); R.to_csv(f'{O}/intermediate/intermediate_bootstrap_ranks.csv')
    n=R.shape[1]
    stab=pd.DataFrame({'top5':(R<=5).mean(1)*100,'top10':(R<=10).mean(1)*100,'median_rank':R.median(1)})
    act=act.set_index('ligand').join(stab).reset_index(); act.to_csv(f'{O}/nichenet_ligand_activity_patient.csv',index=False)
    top=act.sort_values('rank').head(15)
    print(f'successful resamples: {n}'); print(top[['rank','ligand','aupr','top10','top5','median_rank','fib_enriched']].to_string(index=False)); print('SAVED')


# ======================================================================================================================
# Section S5C_ligand_response_test: patient-level ligand-response scores, SMS vs resident macrophages (paired Wilcoxon and mixed model, BH)
# Panels: Supp Fig 5C (significance stars)
# ======================================================================================================================
if _run('S5C_ligand_response_test'):
    """Patient-level significance per ligand. For each of the 175 OGN+RSPO3+ fibroblast ligands, score a ligand-response program
in every patient pseudobulk and compare SELENOP+MMP9+SPP1+ macrophages vs pooled resident SELENOP+ macrophages
(Mac S+XS- + Mac S+SG+).
Pseudobulks: intermediate/intermediate_patient_pseudobulk_counts.csv.gz (patient x group, >= 10 cells; 20 SMS / 54 ref).
Expression: log1p CPM. Ligand net: top 250 NicheNet targets per ligand (ligand_target_matrix), weight = regulatory
potential; score = decoupler ULM per pseudobulk.
Tests per ligand: (a) paired two-sided Wilcoxon, SMS vs resident within the 20 patients with both; (b) mixed model
score ~ group + (1|patient) on all 74 pseudobulks. BH-FDR across ligands for each.
Output: ligand_response_patient_test.csv."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, anndata as ad, decoupler as dc, pyreadr
    import statsmodels.formula.api as smf
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands'; NN='/share/fsmresfiles/ydn4687/nichenet'; NTOP=250
    PB=pd.read_csv(f'{O}/intermediate/intermediate_patient_pseudobulk_counts.csv.gz',index_col=0)
    units=pd.Series(PB.index,index=PB.index); pat=units.str.split('__').str[0]; grp=units.str.split('__').str[1]
    cpm=PB.div(PB.sum(1),axis=0)*1e6; L=np.log1p(cpm)
    L=L.loc[:,(PB>0).mean(0)>=0.2]
    act=pd.read_csv(f'{O}/nichenet_ligand_activity_patient.csv'); LIGS=act['ligand'].tolist()
    ltm=list(pyreadr.read_r(f'{NN}/ligand_target_matrix_nsga2r_final.rds').values())[0]
    genes=[g for g in ltm.index if g in L.columns]; rows=[]
    for lg in LIGS:
        w=ltm.loc[genes,lg].sort_values(ascending=False).head(NTOP)
        rows+= [(lg,g,float(v)) for g,v in w.items() if v>0]
    net=pd.DataFrame(rows,columns=['source','target','weight'])
    A=ad.AnnData(L.values.astype(float),obs=pd.DataFrame({'patient':pat.values,'group':grp.values},index=L.index),var=pd.DataFrame(index=L.columns))
    dc.mt.ulm(data=A,net=net,tmin=10)
    S=dc.pp.get_obsm(A,key='score_ulm'); S=pd.DataFrame(S.X,index=A.obs_names,columns=S.var_names)
    S['patient']=pat.values; S['group']=grp.values
    both=[p for p,g in S.groupby('patient') if set(g.group)=={'target','ref'}]
    print('pseudobulks:',S.group.value_counts().to_dict(),'| patients with both:',len(both),'| unique patients:',S.patient.nunique())
    out=[]
    for lg in [c for c in S.columns if c not in ('patient','group')]:
        w=S[S.patient.isin(both)].pivot(index='patient',columns='group',values=lg); d=w['target']-w['ref']
        pw=stats.wilcoxon(w['target'],w['ref']).pvalue
        est,pm=np.nan,np.nan; t='C(group, Treatment("ref"))[T.target]'
        for meth in ['lbfgs','powell','nm','bfgs']:
            try:
                m=smf.mixedlm(f'Q("{lg}") ~ C(group, Treatment("ref"))',S,groups=S['patient']).fit(reml=True,method=meth)
                if np.isfinite(m.bse[t]): est,pm=m.params[t],m.pvalues[t]; break
            except Exception: pass
        out.append([lg,len(both),d.mean(),int((d>0).sum()),pw,est,pm])
    R=pd.DataFrame(out,columns=['ligand','patients_paired','mean_diff_sms_minus_resident','patients_higher_in_sms','p_paired_wilcoxon','mixed_estimate','p_mixed'])
    R['q_paired']=multipletests(R['p_paired_wilcoxon'],method='fdr_bh')[1]
    R['q_mixed']=multipletests(R['p_mixed'].fillna(1),method='fdr_bh')[1]
    R=act[['ligand','rank','aupr']].merge(R,on='ligand').sort_values('rank')
    R.to_csv(f'{O}/ligand_response_patient_test.csv',index=False)
    print(f"ligands with q_paired < 0.05 and higher in SMS: {int(((R.q_paired<0.05)&(R.mean_diff_sms_minus_resident>0)).sum())} of {len(R)}")
    pd.set_option('display.width',200); print(R.head(15).round(4).to_string(index=False)); print('SAVED')


# ======================================================================================================================
# Section S5C_ligand_activity_plot: fibroblast specificity of each ligand, NicheNet ligand supplementary table and the Supp 5C bar plot
# Panels: Supp Fig 5C, NicheNet ligand supplementary table
# ======================================================================================================================
if _run('S5C_ligand_activity_plot'):
    """NicheNet panel and supplementary table (patient as replicate throughout).
Inputs: nichenet_ligand_activity_patient.csv (S5C_nichenet_ranking + S5C_patient_bootstrap, 100 resamples),
ligand_response_patient_test.csv (S5C_ligand_response_test).
OGN+RSPO3+ fibroblast specificity of each ligand: UC scRNA atlas, fibroblast subtypes (8 fine types). Per patient
(obs['Patient']) with >= 10 OGN+RSPO3+ and >= 10 other fibroblasts: mean normalized expression (expm1 of log-normalized X)
in each group; log2FC = log2((OGN + 0.1) / (other + 0.1)). Paired two-sided Wilcoxon across patients, BH-FDR over all
ligands. Enriched = median per-patient log2FC >= 0.5 AND FDR q < 0.1.
Figure: top 25 ligands by AUPR; stars = FDR q of the patient-level ligand-response test (paired Wilcoxon SMS vs resident
macs, 20 patients; * <0.1 ** <0.05 *** <0.01); AUROC column. Bootstrap top-10 stability is reported in the CSV.
Outputs: ligand_fib_enrichment.csv, supp_table_nichenet_ligands.csv, intermediate/intermediate_nichenet_full_with_enrichment.csv,
fig_nichenet_ligand_activity_patient.png."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    O='/share/fsmresfiles/UC/revision_jenny/SuppFig5C_nichenet_fib_to_mac_ligands'
    LFC_MIN=0.5; Q_MAX=0.10; NTOP=25; MIN=10; PC=0.1
    act=pd.read_csv(f'{O}/nichenet_ligand_activity_patient.csv').sort_values('rank')
    FIB=['ADAMDEC1+ Fib','VSTM2A+ Crypt Top Fib','OGN+RSPO3+ Fib','NRG1+ Crypt Top Fib','Activ Fib','CXCL5+ Activ Fib','T-interact Fib','SELENOP+ Fib']
    scrna=sc.read_h5ad('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad',backed='r')
    idx=np.flatnonzero(scrna.obs['cell_type_short'].isin(FIB).values)
    ligs=[g for g in act['ligand'] if g in scrna.var_names]
    f=scrna[idx,ligs].to_memory()
    X=f.X.toarray() if sp.issparse(f.X) else np.asarray(f.X)
    print('X max (log-normalized expected < ~10):',X.max().round(2))
    E=np.expm1(X); isogn=(f.obs['cell_type_short']=='OGN+RSPO3+ Fib').values; pat=f.obs['Patient'].astype(str).values
    rows=[]
    for p in np.unique(pat):
        m=pat==p; a=m&isogn; b=m&~isogn
        if a.sum()>=MIN and b.sum()>=MIN: rows.append((p,E[a].mean(0),E[b].mean(0)))
    print('patients with >=10 OGN+RSPO3+ and >=10 other fibroblasts:',len(rows))
    OG=np.vstack([r[1] for r in rows]); OT=np.vstack([r[2] for r in rows]); LFC=np.log2((OG+PC)/(OT+PC))
    res=[]
    for j,g in enumerate(ligs):
        d=OG[:,j]-OT[:,j]
        p=stats.wilcoxon(OG[:,j],OT[:,j]).pvalue if np.any(d!=0) else 1.0
        res.append([g,np.median(LFC[:,j]),int((LFC[:,j]>0).sum()),len(rows),p])
    EN=pd.DataFrame(res,columns=['ligand','median_log2FC_ogn_vs_other_fib','n_patients_higher_in_ogn','n_patients','p'])
    EN['q']=multipletests(EN['p'],method='fdr_bh')[1]; EN['enriched']=(EN.median_log2FC_ogn_vs_other_fib>=LFC_MIN)&(EN.q<Q_MAX)
    EN.to_csv(f'{O}/ligand_fib_enrichment.csv',index=False)
    act=act.drop(columns=[c for c in ['fib_enriched','ogn_mean','otherfib_mean'] if c in act]).merge(EN,on='ligand',how='left')
    LR=pd.read_csv(f'{O}/ligand_response_patient_test.csv')[['ligand','patients_paired','patients_higher_in_sms','q_paired']]
    act=act.merge(LR,on='ligand',how='left')
    print(f"enriched ligands (log2FC >= {LFC_MIN}, q < {Q_MAX}): {int(act.enriched.sum())} of {len(act)}")
    s3=lambda v: f'{v:#.3g}'.rstrip('.')
    SUP=str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    def pf(v):
        if v>=0.001: return s3(v)
        m,e=f'{v:.2e}'.split('e'); return f'{m} × 10{str(int(e)).translate(SUP)}'
    top=act.head(NTOP).copy(); n=100
    pd.DataFrame({'Rank':top['rank'],'Ligand':top['ligand'],'AUPR':top['aupr'].map(s3),'AUROC':top['auroc'].map(s3),
      f'Top 10 in patient resamples (% of {n})':top['top10'].round(0).astype(int),f'Top 5 in patient resamples (% of {n})':top['top5'].round(0).astype(int),
      'Median rank in resamples':top['median_rank'].map(lambda v: f'{v:g}'),
      'Patients with higher ligand response in SMS mac':top['patients_higher_in_sms'].astype(int).astype(str)+' of '+top['patients_paired'].astype(int).astype(str),
      'Ligand response FDR q (paired Wilcoxon)':top['q_paired'].map(pf),
      'log2FC, OGN⁺RSPO3⁺ vs other fibroblasts (median over patients)':top['median_log2FC_ogn_vs_other_fib'].map(s3),
      'Enrichment FDR q':top['q'].map(pf),'Enriched in OGN⁺RSPO3⁺ fibroblasts':np.where(top['enriched'],'Yes','No'),
      'Receptors on SELENOP⁺MMP9⁺SPP1⁺ macrophages':top['mac_receptors']}).to_csv(f'{O}/supp_table_nichenet_ligands.csv',index=False)
    act.to_csv(f'{O}/intermediate/intermediate_nichenet_full_with_enrichment.csv',index=False)
    print(top[['rank','ligand','aupr','top10','median_log2FC_ogn_vs_other_fib','q','enriched']].round(3).to_string(index=False))
    # figure
    BLACK='#1b1b1f'; TEAL='#1C7C74'; GRAY='#B9C2C7'
    plt.rcParams.update({'font.size':14,'axes.titlesize':16,'axes.labelsize':15,'axes.titleweight':'normal',
        'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':300,'text.color':BLACK,
        'axes.labelcolor':BLACK,'xtick.color':BLACK,'ytick.color':BLACK,'axes.edgecolor':BLACK})
    t=top.iloc[::-1].reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(5.6,0.24*NTOP+2.0))
    ax.barh(range(len(t)),t['aupr'],color=[TEAL if e else GRAY for e in t['enriched']],height=0.5)
    ax.set_yticks(range(len(t))); ax.set_yticklabels(t['ligand'],fontsize=13)
    star=lambda q: '***' if q<0.01 else '**' if q<0.05 else '*' if q<0.1 else ''
    xa=t['aupr'].max()*1.04
    for i,(a,au,q) in enumerate(zip(t['aupr'],t['auroc'],t['q_paired'])):
        ax.text(a+0.002,i,star(q),va='center',ha='left',fontsize=13); ax.text(xa+0.012,i,f'{au:.2f}',va='center',ha='left',fontsize=12)
    ax.text(xa+0.012,len(t)-0.55,'AUROC',va='bottom',ha='left',fontsize=12)
    ax.set_xlabel('Ligand Activity (AUPR)'); ax.set_ylabel('Ligand')
    ax.set_title('OGN+RSPO3+Fib Ligands Predicting SELENOP+MMP9+SPP1+ Mac Programs',fontsize=12)
    ax.set_xlim(0,t['aupr'].max()*1.32); ax.set_xticks([0,0.1,0.2]); ax.set_ylim(-0.6,len(t)+0.3)
    leg=ax.legend(handles=[Patch(color='none',label='Fibroblast specificity:'),Patch(color=TEAL,label='OGN RSPO3 enriched'),Patch(color=GRAY,label='Not enriched')],
                  loc='upper center',bbox_to_anchor=(0.42,-0.1),ncol=3,frameon=False,fontsize=11,handlelength=1.0,handletextpad=0.4,columnspacing=0.9)
    plt.tight_layout(); fig.savefig(f'{O}/fig_nichenet_ligand_activity_patient.png',bbox_inches='tight'); plt.close(fig)
    print('SAVED')
