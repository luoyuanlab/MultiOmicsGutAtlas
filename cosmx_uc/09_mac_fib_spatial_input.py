#!/usr/bin/env python
"""UC CosMx 6K: shared spatial input and the per-macrophage distance table.

Builds one h5ad with the UC and control CosMx 6K pieces (middle section of each block) and, for every SELENOP+
macrophage, its distance to the nearest (ligand-expressing) OGN+RSPO3+ fibroblast and its receptor detection.

Panels made
  none (input of Fig 3E-F, 3H-J and Supp 6E)

Inputs
  /share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/6k/Processed_merged_upd/h5ad_obj/qc_6k_merged_log1p_anno_6375genes.h5ad (02_6k_processing.ipynb)
  /share/fsmresfiles/UC/AtoMx/<run folders>/processed/proseg_viewer/*_proseg_anno.h5ad (CosMx 6K ProSeg sections)
  /share/fsmresfiles/ydn4687/nichenet/ligand_target_matrix_nsga2r_final.rds

Outputs
  /share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad (1.8 GB)
  /share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance/ (intermediate/ligand_sender_cells.csv.gz, cohort_summary.csv, ligand_programs.csv)

Run order of sections
  F3H_distance_input

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
# Section F3H_distance_input: shared CosMx 6K input h5ad (8 UC + 3 control patients, middle section per block) and per-macrophage distance/receptor table
# Panels: Fig 3H inputs; shared UC CosMx input read by the Fig 3E/3F/3I/3J and Supp 6E sections
# ======================================================================================================================
if _run('F3H_distance_input'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance/01_cosmx_build_input.py'   # working folder of this step (relative paths resolve from here)
    """Builds the shared CosMx 6K input h5ad (8 UC + 3 control patients, middle section of each block) and the
per-macrophage distance/receptor table. Two parts, run in order:
  1. combined input h5ad (uc_cosmx_colocal_input.h5ad) and cohort summary
  2. per-macrophage ligand-sender distance table and ligand response programs"""

    # ====================================================================================================
    # PART 1: combined CosMx input h5ad
    # ====================================================================================================
    """CosMx input for the OGN+RSPO3+ fibroblast <-> SELENOP+MMP9+SPP1+ macrophage co-localization analyses.
One combined h5ad, written to UC_shared_inputs/.

Cohorts:
  * UC CosMx, 6K+custom panel, 4 unique patients (I0261, I0331, I0346, I0351), all slides.
    File qc_6k_merged_log1p_anno_6375genes.h5ad; labels = aucell_cell_type_short (AUCell);
    coords = obsm['spatial_fov'] (global px) x 0.12028 um/px; section unit = obs['slide'].
  * CosMx 6K ProSeg, 11 tissue blocks / 7 patients, ONE section per block: z_index == 2 (middle section).
    Files <piece>_proseg_anno.h5ad; labels = cell_type_short (AUCell); coords = obs centroid_x/centroid_y (um);
    section unit = piece. Condition and Mayo score per patient:
    i900, i0599, i0465 = Control; i0350 = Mayo 1; i0463, i0548 = Mayo 2; i0606 = Mayo 3.
X = raw counts (layers['raw']) on genes shared by both panels.
obs: cohort, patient, block, section_id, z_index, condition (UC/Control; 6K+custom I0346 = healthy control, I0351 = non-inflamed UC, I0261/I0331 = inflamed UC), mayo, cell_type_short, cell_type_general
(both AUCell), x_um, y_um, fov."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
    O='/share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance'; UM=0.12028; ZSEL=2
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
    parts.append(b); print('6K+custom:',b.shape,b.obs.patient.value_counts().to_dict())
    for piece,(path,pat,cond,mayo) in SERIAL.items():
        s=sc.read_h5ad(path,backed='r'); z=s.obs['z_index'].astype(int).values; idx=np.flatnonzero(z==ZSEL)
        s=s[idx].to_memory()
        c=ad.AnnData(X=sp.csr_matrix(s.layers['raw']),var=pd.DataFrame(index=s.var_names))
        c.obs=pd.DataFrame({'cohort':'Serial-section ProSeg 6K','patient':pat,'block':piece,'section_id':f'{piece}_z{ZSEL}','z_index':ZSEL,
            'condition':cond,'mayo':mayo,'cell_type_short':s.obs['cell_type_short'].astype(str).values,'cell_type_general':s.obs['cell_type_general'].astype(str).values,
            'x_um':s.obs['centroid_x'].astype(float).values,'y_um':s.obs['centroid_y'].astype(float).values,
            'fov':s.obs['fov'].astype(str).values},index=(f'{piece}_'+s.obs_names.astype(str)).values)
        parts.append(c); print(piece,c.shape)
    genes=sorted(set.intersection(*[set(p.var_names) for p in parts])); print('shared genes:',len(genes))
    A=ad.concat([p[:,genes] for p in parts],join='inner',merge='same'); A.obs_names_make_unique()
    for c in ['cohort','patient','block','section_id','condition','cell_type_short','cell_type_general','fov']: A.obs[c]=A.obs[c].astype('category')
    A.uns['notes']='Built by cosmx_uc/09_mac_fib_spatial_input.py. X = raw counts. Labels = AUCell. Coordinates in um, per section_id.'
    A.write_h5ad('/share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad',compression='gzip')
    o=A.obs
    summ=(o.groupby(['cohort','patient','condition','mayo'],observed=True,dropna=False)
          .agg(sections=('section_id','nunique'),cells=('cell_type_short','size'),
               sms_mac=('cell_type_short',lambda x:(x=='Mac S+M+S+').sum()),ogn_fib=('cell_type_short',lambda x:(x=='OGN+RSPO3+ Fib').sum()))
          .reset_index())
    summ.to_csv(f'{O}/cohort_summary.csv',index=False); print(summ.to_string(index=False)); print('SAVED',A.shape)


    # ====================================================================================================
    # PART 2: per-macrophage ligand-sender distance table
    # ====================================================================================================
    """Per-macrophage table for the ligand-specific sender analysis (UC 6K+custom + CosMx 6K ProSeg middle sections;
8 UC + 3 control patients; AUCell labels). Cells: SELENOP+ macrophages ('SELENOP+ Mac', raw SELENOP > 0).
For each ligand L in IL33, CXCL12, FGF7, IGF1 (within section x patient):
  d_<L>      = distance (um) to nearest OGN+RSPO3+ fibroblast with raw L > 0 (the ligand-expressing sender), cap 100
  d_other_<L> = distance to nearest fibroblast that is NOT an L+ OGN+RSPO3+ fibroblast (incl. L- OGN+RSPO3+), cap 100
  dmk_<L> = sender is a fibroblast detecting OGN or RSPO3 (marker-defined) and L; dun_<L> = label OR marker OGN+RSPO3+ and L
  (CosMx fine labels are fuzzy: only 24% of labelled OGN+RSPO3+ cells detect OGN/RSPO3)
  daf_<L> / daf_other_<L> = same with the sender broadened to ANY fibroblast expressing L (the OGN+RSPO3+ source
  is established in scRNA-seq)
Also d_ogn / d_otherfib (any OGN+RSPO3+ fibroblast) for SPP1 / MMP9.
Receptors (raw > 0): IL33 IL1RL1, IL1RAP; CXCL12 CXCR4, ACKR3; FGF7 FGFR1, FGFR2, NRP1; IGF1 IGF1R, IGF2R; any_<L>.
Response program per ligand (standard NicheNet ligand-target links): the ligand's top 250 targets by regulatory potential
(ligand_target_matrix); keep those on the panel, detected (raw > 0) in >= 5% of SELENOP+ macrophages and not SPP1, MMP9,
SELENOP or any receptor above; up to 50 per ligand, in order of potential. Per cell: prog_<L>_k = number of program genes detected, prog_<L>_n = genes in set.
Outputs: intermediate/ligand_sender_cells.csv.gz, ligand_programs.csv."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, pyreadr
    from scipy.spatial import cKDTree
    O='/share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance'; NN='/share/fsmresfiles/ydn4687/nichenet'
    OGN='OGN+RSPO3+ Fib'; CAP=100; NTOP=50
    REC={'IL33':['IL1RL1','IL1RAP'],'CXCL12':['CXCR4','ACKR3'],'FGF7':['FGFR1','FGFR2','NRP1'],'IGF1':['IGF1R','IGF2R']}
    A=sc.read_h5ad('/share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad'); o=A.obs.copy()
    X=A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    col=lambda g: np.asarray(X[:,A.var_names.get_loc(g)].todense()).flatten()
    ct=o['cell_type_short'].astype(str).values; isfib=pd.Series(ct).str.contains('Fib').values; isogn=ct==OGN
    mkogn=isfib&((col('OGN')>0)|(col('RSPO3')>0))   # marker-defined OGN/RSPO3+ fibroblast
    unogn=isogn|mkogn                               # union: label OR marker
    bothogn=isogn|(isfib&(col('OGN')>0)&(col('RSPO3')>0))   # union-both: label OR fibroblast detecting BOTH OGN and RSPO3
    mac=(col('SELENOP')>0)&(o['cell_type_general']=='SELENOP+ Mac').values
    print('SELENOP+ macrophages:',int(mac.sum()))
    # response programs
    ltm=list(pyreadr.read_r(f'{NN}/ligand_target_matrix_nsga2r_final.rds').values())[0]
    excl=set(['SPP1','MMP9','SELENOP'])|{g for v in REC.values() for g in v}
    det=pd.Series(np.asarray((X[mac]>0).mean(0)).flatten(),index=A.var_names)
    progs={}
    for L in REC:
        top=ltm[L].sort_values(ascending=False).head(250); keep=[g for g in top.index if g in A.var_names and g not in excl and det.get(g,0)>=0.05]
        progs[L]=keep[:NTOP]
    pd.DataFrame([(L,i+1,g,float(ltm.loc[g,L])) for L,gs in progs.items() for i,g in enumerate(gs)],columns=['ligand','rank','gene','regulatory_potential']).to_csv(f'{O}/ligand_programs.csv',index=False)
    Ls=list(progs); print('program overlap (shared genes):',{f'{a}-{b}':len(set(progs[a])&set(progs[b])) for i,a in enumerate(Ls) for b in Ls[i+1:]})
    print({L:len(v) for L,v in progs.items()}); print({L:v[:10] for L,v in progs.items()})
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
        df['d_ogn']=nd(isogn); df['d_otherfib']=nd(isfib&~isogn)
        df['dub_ogn']=nd(bothogn); df['dub_otherfib']=nd(isfib&~bothogn)
        df['dmk_ogn']=nd(mkogn); df['dmk_otherfib']=nd(isfib&~mkogn); df['dun_ogn']=nd(unogn); df['dun_otherfib']=nd(isfib&~unogn)
        for g in ['SPP1','MMP9']: df[g]=(col(g)[mi]>0).astype(int)
        for L,rs in REC.items():
            sender=isogn&lpos[L]; df[f'd_{L}']=nd(sender); df[f'd_other_{L}']=nd(isfib&~sender)
            afs=isfib&lpos[L]; df[f'daf_{L}']=nd(afs); df[f'daf_other_{L}']=nd(isfib&~afs)   # sender = ANY fibroblast expressing L
            mks=mkogn&lpos[L]; df[f'dmk_{L}']=nd(mks); df[f'dmk_other_{L}']=nd(isfib&~mks)   # sender = OGN/RSPO3-marker fibroblast expressing L
            uns=unogn&lpos[L]; df[f'dun_{L}']=nd(uns); df[f'dun_other_{L}']=nd(isfib&~uns)   # sender = label OR marker OGN+RSPO3+, expressing L
            ubs=bothogn&lpos[L]; df[f'dub_{L}']=nd(ubs); df[f'dub_other_{L}']=nd(isfib&~ubs)   # sender = label OR (OGN AND RSPO3), expressing L
            for r in rs: df[r]=(col(r)[mi]>0).astype(int)
            df[f'any_{L}']=df[rs].max(1)
            gi=[A.var_names.get_loc(g) for g in progs[L]]
            df[f'prog_{L}_k']=np.asarray(B[mi][:,gi].sum(1)).flatten(); df[f'prog_{L}_n']=len(gi)
        rows.append(df)
    C=pd.concat(rows,ignore_index=True); C.to_csv(f'{O}/intermediate/ligand_sender_cells.csv.gz',index=False)
    print(C.groupby('condition').agg(cells=('SPP1','size'),patients=('patient','nunique'),samples=('sample','nunique')))
    print('UC macrophages with a ligand+ OGN+RSPO3+ sender within 100 / 50 um (%):',{L:(round((C[C.condition=="UC"][f'd_{L}']<CAP).mean()*100,1),round((C[C.condition=="UC"][f'd_{L}']<=50).mean()*100,1)) for L in REC})
    print('SAVED')
