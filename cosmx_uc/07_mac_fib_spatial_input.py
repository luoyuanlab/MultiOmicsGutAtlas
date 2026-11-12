#!/usr/bin/env python
"""UC CosMx 6K: shared spatial input and the per-macrophage distance table.

One h5ad with the UC and control CosMx 6K pieces (middle section z = 2 of each serial block; raw counts, AUCell labels,
coordinates in um) and, for every SELENOP+ macrophage, its distance to the nearest ligand-expressing OGN+RSPO3+
fibroblast ('union-both': label OR fibroblast detecting both OGN and RSPO3; capped at 100 um), receptor detection and
NicheNet response-program counts.

Inputs:  UC 6K+custom h5ad (Processed_merged_upd; its processing is not in this repository), CosMx 6K ProSeg *_proseg_anno.h5ad, NicheNet ligand_target_matrix_nsga2r_final.rds
Outputs: UC_shared_inputs/uc_cosmx_colocal_input.h5ad (read by 09, 10);
         Fig2H_cosmx_fib_mac_distance/intermediate/ligand_sender_cells.csv.gz (read by 11, 12), ligand_programs.csv
Run:     python <file>, or RUN_SECTIONS=<id>[,<id>] python <file>
"""
import os as _os
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# Section F3H_distance_input: shared input h5ad (8 UC + 3 control patients) and per-macrophage distance table
if _run('F3H_distance_input'):
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
    from scipy.spatial import cKDTree
    INPUT='/share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad'
    UM=0.12028; ZSEL=2
    T='/share/fsmresfiles/UC/AtoMx/07212026_transfer'
    SERIAL={'i900':('/share/fsmresfiles/UC/AtoMx/03312026/stim_alignment_ar2_250hvg/proseg_viewer/i900_proseg_anno.h5ad','i900','Control',np.nan),
     'i0548':('/share/fsmresfiles/UC/AtoMx/Cos16KAddOnDI0548RSDF00AMayo2Tissues1841726_01_05_2026_22_52_39_21/processed/proseg_viewer/i0548_proseg_anno.h5ad','i0548','UC',2),
     'i0606':('/share/fsmresfiles/UC/AtoMx/Cos16KADDOnDI0606RSDF00AMayo3Tissues1641726_2_07_05_2026_15_44_18_459/processed/proseg_viewer/i0606_proseg_anno.h5ad','i0606','UC',3)}
    for pc,run,cond,mayo in [('i0463','Cos56KAddOnI0463RSDF00AMayo2Tissues165142026_29_05_2026_23_15_18_655','UC',2),
                             ('i0350','Cos56KAddOnI0350RSDF00AMayo1Tissues1651426_29_05_2026_23_14_18_419','UC',1),
                             ('i0599','Cos46KAddOnDI0599RSCF00AControlTissues1651426_29_05_2026_23_18_03_995','Control',np.nan),
                             ('i0465','Cos46KAddOnI0465RSCF00AControlTissues1651426_29_05_2026_23_16_52_323','Control',np.nan)]:
        for k in [1,2]: SERIAL[f'{pc}_piece{k}']=(f'{T}/{run}/processed/proseg_viewer/{pc}_piece{k}_proseg_anno.h5ad',pc,cond,mayo)
    parts=[]
    # UC 6K+custom
    a=sc.read_h5ad('/share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/6k/Processed_merged_upd/h5ad_obj/qc_6k_merged_log1p_anno_6375genes.h5ad')
    xy=np.asarray(a.obsm['spatial_fov'],float)*UM
    b=ad.AnnData(X=sp.csr_matrix(a.layers['raw']),var=pd.DataFrame(index=a.var_names))
    b.obs=pd.DataFrame({'cohort':'UC 6K+custom','patient':a.obs['patient'].astype(str).values,'block':a.obs['slide'].astype(str).values,
        'section_id':a.obs['slide'].astype(str).values,'z_index':np.nan,'condition':np.where(a.obs['patient'].astype(str).values=='I0346','Control','UC'),'mayo':np.nan,
        'cell_type_short':a.obs['aucell_cell_type_short'].astype(str).values,'cell_type_general':a.obs['aucell_cell_type_general'].astype(str).values,'x_um':xy[:,0],'y_um':xy[:,1],
        'fov':a.obs['fov'].astype(str).values},index=('6K_'+a.obs_names.astype(str)).values)
    parts.append(b)
    for piece,(path,pat,cond,mayo) in SERIAL.items():
        s=sc.read_h5ad(path,backed='r'); z=s.obs['z_index'].astype(int).values; idx=np.flatnonzero(z==ZSEL)
        s=s[idx].to_memory()
        c=ad.AnnData(X=sp.csr_matrix(s.layers['raw']),var=pd.DataFrame(index=s.var_names))
        c.obs=pd.DataFrame({'cohort':'Serial-section ProSeg 6K','patient':pat,'block':piece,'section_id':f'{piece}_z{ZSEL}','z_index':ZSEL,
            'condition':cond,'mayo':mayo,'cell_type_short':s.obs['cell_type_short'].astype(str).values,'cell_type_general':s.obs['cell_type_general'].astype(str).values,
            'x_um':s.obs['centroid_x'].astype(float).values,'y_um':s.obs['centroid_y'].astype(float).values,
            'fov':s.obs['fov'].astype(str).values},index=(f'{piece}_'+s.obs_names.astype(str)).values)
        parts.append(c)
    genes=sorted(set.intersection(*[set(p.var_names) for p in parts]))
    A=ad.concat([p[:,genes] for p in parts],join='inner',merge='same'); A.obs_names_make_unique()
    for c in ['cohort','patient','block','section_id','condition','cell_type_short','cell_type_general','fov']: A.obs[c]=A.obs[c].astype('category')
    A.uns['notes']='Built by cosmx_uc/07_mac_fib_spatial_input.py. X = raw counts. Labels = AUCell. Coordinates in um, per section_id.'
    A.write_h5ad(INPUT,compression='gzip')
    summ=(A.obs.groupby(['cohort','patient','condition','mayo'],observed=True,dropna=False)
          .agg(sections=('section_id','nunique'),cells=('cell_type_short','size'),
               sms_mac=('cell_type_short',lambda x:(x=='Mac S+M+S+').sum()),ogn_fib=('cell_type_short',lambda x:(x=='OGN+RSPO3+ Fib').sum()))
          .reset_index())
    summ.to_csv('/share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance/cohort_summary.csv',index=False)

    # per-macrophage distance, receptor and response-program table
    import pyreadr
    O='/share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance'; NN='/share/fsmresfiles/ydn4687/nichenet'
    OGN='OGN+RSPO3+ Fib'; CAP=100; NTOP=50
    REC={'IL33':['IL1RL1','IL1RAP'],'CXCL12':['CXCR4','ACKR3'],'FGF7':['FGFR1','FGFR2','NRP1'],'IGF1':['IGF1R','IGF2R']}
    A=sc.read_h5ad(INPUT); o=A.obs.copy()
    X=A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    col=lambda g: np.asarray(X[:,A.var_names.get_loc(g)].todense()).flatten()
    ct=o['cell_type_short'].astype(str).values; isfib=pd.Series(ct).str.contains('Fib').values; isogn=ct==OGN
    bothogn=isogn|(isfib&(col('OGN')>0)&(col('RSPO3')>0))   # union-both: label OR fibroblast detecting BOTH OGN and RSPO3
    mac=(col('SELENOP')>0)&(o['cell_type_general']=='SELENOP+ Mac').values
    # response programs
    ltm=list(pyreadr.read_r(f'{NN}/ligand_target_matrix_nsga2r_final.rds').values())[0]
    excl=set(['SPP1','MMP9','SELENOP'])|{g for v in REC.values() for g in v}
    det=pd.Series(np.asarray((X[mac]>0).mean(0)).flatten(),index=A.var_names)
    progs={}
    for L in REC:
        top=ltm[L].sort_values(ascending=False).head(250); keep=[g for g in top.index if g in A.var_names and g not in excl and det.get(g,0)>=0.05]
        progs[L]=keep[:NTOP]
    pd.DataFrame([(L,i+1,g,float(ltm.loc[g,L])) for L,gs in progs.items() for i,g in enumerate(gs)],columns=['ligand','rank','gene','regulatory_potential']).to_csv(f'{O}/ligand_programs.csv',index=False)
    # per-cell table
    o['unit']=o['section_id'].astype(str)+'|'+o['patient'].astype(str); xy=o[['x_um','y_um']].values
    lpos={L:col(L)>0 for L in REC}; tot=np.asarray(X.sum(1)).flatten()
    B=(X>0).tocsr()
    rows=[]
    for u,d in o.groupby('unit'):
        ii=np.flatnonzero((o['unit']==u).values); mi=ii[mac[ii]]
        if len(mi)==0: continue
        def nd(mask):
            pts=xy[ii][mask[ii]]
            return np.minimum(cKDTree(pts).query(xy[mi])[0],CAP) if len(pts) else np.full(len(mi),float(CAP))
        df=pd.DataFrame({'sample':u,'patient':d['patient'].iat[0],'condition':d['condition'].iat[0],'log_total':np.log1p(tot[mi])})
        df['dub_ogn']=nd(bothogn); df['dub_otherfib']=nd(isfib&~bothogn)
        for g in ['SPP1','MMP9']: df[g]=(col(g)[mi]>0).astype(int)
        for L,rs in REC.items():
            ubs=bothogn&lpos[L]; df[f'dub_{L}']=nd(ubs); df[f'dub_other_{L}']=nd(isfib&~ubs)   # sender = label OR (OGN AND RSPO3), expressing L
            for r in rs: df[r]=(col(r)[mi]>0).astype(int)
            df[f'any_{L}']=df[rs].max(1)
            gi=[A.var_names.get_loc(g) for g in progs[L]]
            df[f'prog_{L}_k']=np.asarray(B[mi][:,gi].sum(1)).flatten(); df[f'prog_{L}_n']=len(gi)
        rows.append(df)
    C=pd.concat(rows,ignore_index=True); C.to_csv(f'{O}/intermediate/ligand_sender_cells.csv.gz',index=False)
