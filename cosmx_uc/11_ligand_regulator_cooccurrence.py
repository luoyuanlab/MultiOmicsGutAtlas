#!/usr/bin/env python
"""UC CosMx (1K imputed, 6K, WTX): macrophage ligand vs fibroblast regulator co-expression, co-occurrence and neighbourhoods.

Panels made
  Fig 3D  NicheNet ligand activity of Mac S+M+S+ ligands (scRNA-seq), with % expressing in WTX
  Fig 3E  ligand in the nearest Mac S+M+S+ vs regulator in each OGN+RSPO3+ fibroblast, Spearman r pooled by random-effects meta-analysis
  Fig 3I  co-occurrence of activated fibroblasts with myeloid and fibroblast populations
  Fig 3J  odds of fibroblast activation per neighbouring myeloid cell (20 nearest cells, SEs clustered by piece)
  Supp tables for Fig 3E and 3J

Inputs
  /share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad (09_mac_fib_spatial_input.py)
  /share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/I0286RI0277R1119/ (1K validation object and imputed AUCell labels, 04-05)
  /share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/whole_trans/Processed_merged/v7/h5ad_obj/qc_wtx_v7_merged_anno.h5ad (01_wtx_processing.ipynb, 05_aucell_annotation.R)
  /share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/ (cell_type_name_mapper.csv; ogn_spp1_mac_nichenet_ligand_activities.csv from
    uc_scrna/14_nichenet_macrophage_to_fibroblast.R; ogn_spp1_mac_nichenet_tf_spatial_validation_combined_v2.csv, the
    ligand-regulator pair list from 10_wtx_ligand_regulator_pairs.ipynb)

Outputs
  /share/fsmresfiles/UC/revision_jenny/Fig2E_I_J_ligand_regulator_cooccurrence_neighbourhood/ (panels, supp tables, intermediate/fig2c_per_piece.csv)
  /share/fsmresfiles/UC/revision_jenny/Fig2D_nichenet_ligand_activity/fig1i_ligand_activity.png

Run order of sections
  F3EIJ_ligand_regulator_cooccurrence -> F3D_ligand_activity

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
# Section F3EIJ_ligand_regulator_cooccurrence: ligand-regulator spatial correlation, activated-fibroblast co-occurrence and neighbour models (inflamed UC)
# Panels: Fig 3E, 3I, 3J
# ======================================================================================================================
if _run('F3EIJ_ligand_regulator_cooccurrence'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig2E_I_J_ligand_regulator_cooccurrence_neighbourhood/01_fig2cmn_allsamples.py'   # working folder of this step (relative paths resolve from here)
    """Ligand-regulator spatial correlation, activated-fibroblast co-occurrence and neighbour models in all inflamed UC samples:
1K validation (I0277, I0286; 4 pieces; imputed-1K AUCell labels) + CosMx 6K ProSeg (middle section; i0350, i0463, i0548, i0606;
6 pieces; AUCell labels) + WTX (I0261, I0331; inflamed; aucell_cell_type_short). Piece = tissue piece within a patient.

Fig 3E: top 50 NicheNet ligands (Mac S+M+S+ -> OGN+RSPO3+ Fib, AUROC) present in the panel; regulators per ligand from
ogn_spp1_mac_nichenet_tf_spatial_validation_combined_v2.csv (exact ligand match); per piece (>= 10 OGN+RSPO3+ Fib,
>= 5 Mac S+M+S+), Spearman r of the ligand in the nearest Mac S+M+S+ vs the regulator in each OGN+RSPO3+ Fib. WTX and
6K pieces combined; per-piece r pooled by random-effects meta-analysis (95% CI); pairs with |r| >= 0.03 shown (top 19),
grouped by P bin (<0.05, <0.1, >=0.1) and ranked by |r| within bin; bars colored by P; BH-FDR across all pairs in the CSV.
Expression: log-normalized.
Fig 3J: fibroblast activation ~ counts of each myeloid type among the 20 nearest cells (per piece, global coordinates),
logistic regression, SEs clustered by piece; BH-FDR over the 8 predictors.
Fig 3I: co-occurrence score (squidpy definition) of activated fibroblasts with fibroblast/myeloid groups, pair counts within
pieces pooled over pieces, 20-um shells to 500 um; DC, ADAMDEC1+ Fib, Trans cDC2/Mac, T-interact Fib and Crypt Top Fib not
plotted.
Outputs: fig2c_bar_combined.png, fig2n_forest_allsamples.png, fig2m_cooccurrence_allsamples.png,
supp_table_fig2c_allsamples.csv, supp_table_fig2n_allsamples.csv."""
    import warnings; warnings.filterwarnings('ignore')
    import os, numpy as np, pandas as pd, scanpy as sc, squidpy as sq, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from scipy.spatial import cKDTree
    from scipy import stats
    from sklearn.neighbors import NearestNeighbors
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
    O='/share/fsmresfiles/UC/revision_jenny/Fig2E_I_J_ligand_regulator_cooccurrence_neighbourhood'; PX=0.12028; NN='/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip'
    MAP=pd.read_csv(f'{NN}/cell_type_name_mapper.csv'); CAT=dict(zip(MAP['cell_type_short'],MAP['cell_category']))
    def lognorm(a):
        a=a.copy(); sc.pp.normalize_total(a); sc.pp.log1p(a); return a
    # ---------------- data (inflamed UC only) ----------------
    S=sc.read_h5ad('/share/fsmresfiles/UC/revision_jenny/UC_shared_inputs/uc_cosmx_colocal_input.h5ad')
    S=lognorm(S[(S.obs['cohort'].astype(str).str.contains('Serial'))&(S.obs['condition'].astype(str)=='UC')].copy())
    S.obs['ct']=S.obs['cell_type_short'].astype(str); S.obs['piece']=S.obs['block'].astype(str); S.obsm['xy']=S.obs[['x_um','y_um']].values
    V=sc.read_h5ad('/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/I0286RI0277R1119/h5ad_obj/qc_norm_log1p.h5ad')
    V.obs['ct']=pd.read_csv('/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/I0286RI0277R1119/anno/imp_log1p_counts_aucell.csv')['label'].astype(str).values
    V.obs['piece']=V.obs['slide'].astype(str); V.obsm['xy']=V.obsm['spatial_fov']*PX
    W=sc.read_h5ad('/share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/whole_trans/Processed_merged/v7/h5ad_obj/qc_wtx_v7_merged_anno.h5ad')
    W=W[W.obs['patient'].astype(str).isin(['I0261','I0331'])].copy()
    W.obs['ct']=W.obs['aucell_cell_type_short'].astype(str); W.obs['piece']=W.obs['slide'].astype(str)+'|'+W.obs['patient'].astype(str); W.obsm['xy']=W.obsm['spatial_fov']*PX
    for a,name in [(V,'1K validation'),(S,'Serial 6K'),(W,'WTX')]: a.obs['patient']=a.obs['patient'].astype(str); a.obs['cohort2']=name
    COH=[V,S,W]

    # ---------------- Fig 3E: ligand-regulator spatial correlation ----------------
    la=pd.read_csv(f'{NN}/ogn_spp1_mac_nichenet_ligand_activities.csv'); combined=pd.read_csv(f'{NN}/ogn_spp1_mac_nichenet_tf_spatial_validation_combined_v2.csv')
    TOP=la.nlargest(50,'auroc')['test_ligand'].tolist()
    def ligand_tf_corr(a,slide_col):
        out=[]
        for ligand in [l for l in TOP if l in a.var_names]:
            ligand_tfs=combined[combined['upstream_ligands'].fillna('').map(lambda x: ligand in [t.strip() for t in x.split(',')])]['TF'].tolist()   # exact match (substring matching would also match e.g. MMP14 for MMP1)
            for sid in a.obs[slide_col].unique():
                sd=a[a.obs[slide_col]==sid]; om=(sd.obs['ct']=='OGN+RSPO3+ Fib').values; mm=(sd.obs['ct']=='Mac S+M+S+').values
                if om.sum()<10 or mm.sum()<5: continue
                lx=sd[mm,ligand].X; lx=lx.toarray().ravel() if hasattr(lx,'toarray') else np.asarray(lx).ravel()
                dist,nn=cKDTree(sd.obsm['xy'][mm]).query(sd.obsm['xy'][om],k=1); lnear=lx[nn]
                for tf in ligand_tfs:
                    if tf not in sd.var_names: continue
                    tx=sd[om,tf].X; tx=tx.toarray().ravel() if hasattr(tx,'toarray') else np.asarray(tx).ravel()
                    if len(np.unique(lnear))>1:
                        r,p=stats.spearmanr(lnear,tx); out.append({'slide':sid,'ligand':ligand,'TF':tf,'correlation':r,'pvalue':p,'n_fibs':len(tx)})
        return pd.DataFrame(out)
    # per-piece Spearman r pooled across pieces by random-effects meta-analysis (DerSimonian-Laird) of Fisher-z, var = 1.06/(n-3)
    # (Spearman); gives pooled r, 95% CI and P; bars colored by nominal P; BH-FDR across all pairs in the CSV.
    def meta(d):
        z=np.arctanh(np.clip(d['correlation'].values,-0.999,0.999)); v=1.06/(d['n_fibs'].values-3); w=1/v
        zf=(w*z).sum()/w.sum(); Q=(w*(z-zf)**2).sum(); k=len(z)
        tau2=max(0.0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())) if k>1 else 0.0
        ws=1/(v+tau2); zr=(ws*z).sum()/ws.sum(); se=np.sqrt(1/ws.sum()); p=2*stats.norm.sf(abs(zr/se))
        return pd.Series({'r':np.tanh(zr),'r_lo':np.tanh(zr-1.96*se),'r_hi':np.tanh(zr+1.96*se),'p':p,'pieces':k,'n_fibs':int(d['n_fibs'].sum())})
    RW=ligand_tf_corr(W,'slide').assign(cohort='WTX'); RS=ligand_tf_corr(S,'piece').assign(cohort='Serial 6K')   # WTX per slide
    RAW=pd.concat([RW,RS],ignore_index=True); RAW['piece']=RAW['cohort']+'|'+RAW['slide'].astype(str); RAW.to_csv(f'{O}/intermediate/fig2c_per_piece.csv',index=False)
    res=[]
    for name,d0 in [('Combined',RAW),('WTX',RW),('Serial 6K',RS)]:
        m=d0.groupby(['ligand','TF']).apply(meta).reset_index(); m['q']=stats.false_discovery_control(m['p']); m['analysis']=name; res.append(m)   # BH-FDR across all pairs (reported in CSV/text)
    M=pd.concat(res,ignore_index=True)
    CMB=M[(M.analysis=='Combined')&(M.r.abs()>=0.03)].assign(nr=lambda x: -x.r.abs(),pbin=lambda x: np.where(x.p<0.05,0,np.where(x.p<0.1,1,2))).sort_values(['pbin','nr']).head(19)   # grouped by P bin, |r| descending within bin; display |r| >= 0.03; FDR over all pairs
    npat_c=len(set(W.obs.patient)|set(S.obs.patient))
    fig,ax=plt.subplots(figsize=(8.2,max(3.6,0.24*len(CMB)+1.4))); y=np.arange(len(CMB))
    cols=['#D62828' if p<0.05 else '#F77F00' if p<0.1 else '#FCBF49' for p in CMB['p']]
    ax.barh(y,CMB['r'],color=cols,edgecolor='black',linewidth=1.2)
    ax.errorbar(CMB['r'],y,xerr=np.vstack([CMB['r']-CMB['r_lo'],CMB['r_hi']-CMB['r']]),fmt='none',ecolor='black',elinewidth=1,capsize=2.5)
    ax.set_yticks(y); ax.set_yticklabels([f"{l} → {t}" for l,t in zip(CMB.ligand,CMB.TF)],fontsize=14)
    ax.set_ylabel('Ligand (Mac S+M+S+) →\nTF (OGN+RSPO3+ Fib)',fontsize=12); ax.set_xlabel('Spatial Correlation (r) of Ligand and TF (95% CI)',fontsize=13)
    ax.set_title('Spatial Validation of Ligand → TF Expression: WTX + 6K',fontsize=13,loc='left'); ax.axvline(0,color='black',linewidth=1)
    lim=float(max(CMB['r_hi'].abs().max(),CMB['r_lo'].abs().max()))*1.25; ax.set_xlim(-lim,lim); ax.invert_yaxis()
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    for i,(_,r) in enumerate(CMB.iterrows()):
        xe=r['r_hi'] if r['r']>0 else r['r_lo']
        ax.text(xe+(0.006 if r['r']>0 else -0.006),i,f"r:{r['r']:.3f}",va='center',ha='left' if r['r']>0 else 'right',fontsize=9.5)
    ax.legend(handles=[Patch(facecolor='#D62828',edgecolor='black',label='P < 0.05'),Patch(facecolor='#F77F00',edgecolor='black',label='P < 0.1'),
                       Patch(facecolor='#FCBF49',edgecolor='black',label='P ≥ 0.1')],loc='lower right',frameon=False,fontsize=11)
    plt.tight_layout(); fig.savefig(f'{O}/fig2c_bar_combined.png',bbox_inches='tight',dpi=300); plt.close(fig)
    g3=lambda x: f'{x:#.3g}'.rstrip('.')
    M.assign(nr=lambda x: -x.r.abs(),analysis=lambda x: x.analysis.replace({'Serial 6K':'6K','Combined':'WTX + 6K'})).sort_values(['analysis','p','nr']).assign(**{'Spearman r (95% CI)':lambda x: [f'{g3(a)} ({g3(b)}–{g3(c)})' for a,b,c in zip(x.r,x.r_lo,x.r_hi)],
        'P':lambda x: x.p.map(g3),'FDR q (all pairs)':lambda x: x.q.map(g3)})[['analysis','ligand','TF','Spearman r (95% CI)','P','FDR q (all pairs)','pieces','n_fibs']].rename(columns={
        'analysis':'Cohort','ligand':'Ligand (Mac S+M+S+)','TF':'Regulator (OGN+RSPO3+ Fib)','pieces':'Tissue pieces','n_fibs':'Fibroblasts'}).to_csv(f'{O}/supp_table_fig2c_allsamples.csv',index=False)
    print(M[(M.ligand=='MMP1')&(M.TF=='TNC')].round(4).to_string(index=False)); print(CMB.round(4).to_string(index=False))

    # ---------------- Fig 3J: myeloid neighbours and fibroblast activation ----------------
    EXC={'Cycl Fib','CD36+ Peri','RERGL+ Contr Peri','SELENOP+ Fib'}
    GRP={'density_Mac_SEL':['Mac S+XS-','Mac S+SG+'],'density_DC':['cDC1','cDC2','CCR7+ DC'],'density_Mac_SMS':['Mac S+M+S+'],
         'density_Neutrophil':['Neutrophil'],'density_Mac_SMP':['Mac S+M+P+'],'density_Trans_cDC2_Mac':['Trans cDC2/Mac'],
         'density_Inf_Mo_Mac':['Inf Mo-Mac'],'density_Mo_Mac':['Mo-Mac']}
    label_map={"density_Mo_Mac":"Mo-Mac","density_Inf_Mo_Mac":"Inf Mo-Mac","density_Neutrophil":"Neutrophil","density_DC":"DC",
               "density_Mac_SMS":"Mac S+M+S+","density_Mac_SMP":"Mac S+M+P+","density_Mac_SEL":"Mac S+SG+","density_Trans_cDC2_Mac":"Trans cDC2/Mac"}
    rows=[]
    for a in COH:
        for pc in a.obs['piece'].unique():
            pm=(a.obs['piece']==pc).values; lab=a.obs['ct'].values[pm]; xy=a.obsm['xy'][pm]
            _,idx=NearestNeighbors(n_neighbors=20).fit(xy).kneighbors(xy)
            isfib=np.array([CAT.get(c)=='Connective' and c not in EXC for c in lab])
            df=pd.DataFrame({'act_fib':(lab=='Activ Fib').astype(int)[isfib],'sample':pc,'patient':a.obs['patient'].values[pm][isfib],'cohort':a.obs['cohort2'].iat[0]})
            for col,types in GRP.items(): df[col]=np.isin(lab,types)[idx].sum(1)[isfib]
            rows.append(df)
    df=pd.concat(rows,ignore_index=True); F='act_fib ~ '+' + '.join(GRP)
    res=[]
    for cl in ['sample','patient']:
        r=smf.logit(F,data=df).fit(disp=False,cov_type='cluster',cov_kwds={'groups':pd.factorize(df[cl])[0]}); ci=np.exp(r.conf_int())
        t=pd.DataFrame({'clustered_by':cl,'Predictor':[label_map[c] for c in GRP],'OR':np.exp(r.params[list(GRP)]).values,'CI_low':ci.loc[list(GRP),0].values,
                        'CI_high':ci.loc[list(GRP),1].values,'Pval':r.pvalues[list(GRP)].values})
        t['adj_p']=multipletests(t['Pval'],method='fdr_bh')[1]; res.append(t)
    R=pd.concat(res,ignore_index=True)
    npat,npc=df.patient.nunique(),df['sample'].nunique()
    R.assign(**{'OR (95% CI)':[f'{g3(a)} ({g3(b)}–{g3(c)})' for a,b,c in zip(R.OR,R.CI_low,R.CI_high)],'FDR q':R.adj_p.map(g3),
                'SE clustered by':R.clustered_by.map({'sample':'tissue piece','patient':'patient'}),'Patients':npat,'Tissue pieces':npc,'Fibroblasts':len(df)})[
        ['SE clustered by','Predictor','OR (95% CI)','FDR q','Patients','Tissue pieces','Fibroblasts']].rename(columns={'Predictor':'Myeloid neighbour'}).to_csv(f'{O}/supp_table_fig2n_allsamples.csv',index=False)
    print(R.round(4).to_string(index=False))
    # forest plot; stars = BH-FDR q < 0.1 / 0.05 / 0.01
    custom_palette=["#FF9999","#66B3FF","#99FF99","#FFD700","#FF69B4","#8A2BE2","#00BFFF","#ADFF2F","#FF4500","#DDA0DD","#20B2AA","#778899",
                    "#FFA500","#7B68EE","#D2691E","#00FA9A","#B22222","#32CD32"]
    colormap=ListedColormap(custom_palette[1:9])
    cell=R[R.clustered_by=='sample'].sort_values('OR',ascending=False).reset_index(drop=True)
    color_dict=dict(zip(sorted(cell['Predictor'].unique()),colormap.colors[:len(cell)])); colors=[color_dict[p] for p in cell['Predictor']]
    xerr=np.vstack([cell['OR']-cell['CI_low'],cell['CI_high']-cell['OR']])
    fig,ax=plt.subplots(figsize=(4,5))
    ax.errorbar(cell['OR'],cell['Predictor'],xerr=xerr,fmt='o',color='black',ecolor='gray',capsize=3,linestyle='none',markerfacecolor='none',markeredgecolor='none')
    ax.scatter(cell['OR'],cell['Predictor'],c=colors,s=70,edgecolor='black')
    for i,row in cell.iterrows():
        s='***' if row['adj_p']<0.01 else '**' if row['adj_p']<0.05 else '*' if row['adj_p']<0.1 else ''
        if s: ax.text(row['CI_high']+0.03,i,s,va='center',ha='left',fontsize=13,fontweight='bold')
    ax.axvline(1,color='red',linestyle='--')
    plt.xlim(0.6,float(cell['CI_high'].max())+0.35)
    ax.set_xlabel('OR for Fib Activation',fontsize=13); ax.set_ylabel('Cell Type (Neighborhood Density)',fontsize=14)
    ax.tick_params(axis='x',labelsize=12); ax.tick_params(axis='y',labelsize=12); ax.invert_yaxis()
    plt.tight_layout(); fig.savefig(f'{O}/fig2n_forest_allsamples.png',bbox_inches='tight',dpi=300); plt.close(fig)

    # ---------------- Fig 3I: co-occurrence with activated fibroblasts ----------------
    # Co-occurrence score as defined in squidpy (P(group | partner of an activated fibroblast at distance d) / P(group | any pair
    # at distance d)), computed directly: pairs counted within each piece in 10-um distance shells, counts pooled over pieces.
    # (computed directly because squidpy co_occurrence on one object holding spatially offset pieces does not give scores centred on 1)
    def fm_group(i):
        if i in ['Mac S+XS-','Mac S+SG+']: return 'Mac S+SG+'
        if i in ['cDC1','cDC2','CCR7+ DC']: return 'DC'
        if 'Peri' in i: return 'Pericyte'
        if 'Myofib' in i: return 'Myofib'
        if 'Crypt Top' in i: return 'Crypt Top Fib'
        return i
    EDGES=np.arange(0,520,20.0); cats=None; N=None; ncell=0
    pieces=[]
    for a in COH:
        for pc in a.obs['piece'].unique():
            pm=(a.obs['piece']==pc).values; lab=a.obs['ct'].values[pm]
            keep=np.array([CAT.get(c) in ('Myeloid','Connective') for c in lab]); g=np.array([fm_group(c) for c in lab[keep]])
            ok=~np.isin(g,['Cycl Myeloid','Cycl Fib','Pericyte','SELENOP+ Fib']); pieces.append((g[ok],a.obsm['xy'][pm][keep][ok])); ncell+=ok.sum()
    cats=sorted(set(np.concatenate([p[0] for p in pieces]))); ci={c:i for i,c in enumerate(cats)}
    N=np.zeros((len(EDGES)-1,len(cats),len(cats)))
    for g,xy in pieces:
        pr=cKDTree(xy).query_pairs(r=EDGES[-1],output_type='ndarray')
        if len(pr)==0: continue
        d=np.linalg.norm(xy[pr[:,0]]-xy[pr[:,1]],axis=1); k=np.digitize(d,EDGES)-1; gi=np.array([ci[x] for x in g]); nc=len(cats)
        ok=(k>=0)&(k<len(EDGES)-1); k,a_,b_=k[ok],gi[pr[ok,0]],gi[pr[ok,1]]
        N+=(np.bincount(k*nc*nc+a_*nc+b_,minlength=N.size)+np.bincount(k*nc*nc+b_*nc+a_,minlength=N.size)).reshape(N.shape)
    i=ci['Activ Fib']
    score=np.array([(N[k][i]/N[k][i].sum())/(N[k].sum(0)/N[k].sum()) for k in range(len(EDGES)-1)])
    SC=pd.DataFrame(score,columns=cats,index=EDGES[1:]).rename_axis('distance_um'); SC.to_csv(f'{O}/intermediate/fig2m_allsamples_scores.csv')
    plt.rcParams.update({'font.size':14,'axes.titlesize':14,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11,'legend.title_fontsize':12})
    fig,ax=plt.subplots(figsize=(4,5))
    DROP={'DC','ADAMDEC1+ Fib','Trans cDC2/Mac','T-interact Fib','Crypt Top Fib'}   # groups not shown in the panel
    for j,c in enumerate(cats):
        if c not in DROP: ax.plot(SC.index,SC[c],color=custom_palette[j%len(custom_palette)],label=c)
    ax.axhline(1,color='gray',linestyle='--',linewidth=1)
    ax.set_title(''); ax.set_ylabel('Co-occurrence Score w/ Activ Fib'); ax.set_xlabel('Distance (µm)'); ax.set_xlim(20,500)
    ax.legend(title='Cell Type',frameon=False,loc='center left',bbox_to_anchor=(1.02,0.5),fontsize=9,title_fontsize=9,handlelength=1.2,handletextpad=0.4)
    plt.tight_layout(); fig.savefig(f'{O}/fig2m_cooccurrence_allsamples.png',bbox_inches='tight',dpi=300); plt.close(fig)
    print('patients',npat,'pieces',npc,'cells in 2M',ncell); print('SAVED')


# ======================================================================================================================
# Section F3D_ligand_activity: NicheNet ligand activity of Mac S+M+S+ ligands, highlighting ligands with spatial support in Fig 3E
# Panels: Fig 3D
# ======================================================================================================================
if _run('F3D_ligand_activity'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/Fig2D_nichenet_ligand_activity/01_fig1i_ligand_activity.py'   # working folder of this step (relative paths resolve from here)
    """NicheNet ligand activity plot. Highlighted: ligands with at least one ligand-regulator pair at P < 0.1 in the pooled
WTX + 6K Fig 3E analysis (supp_table_fig2c_allsamples.csv, cohort 'Combined') -> MMP1, CTSD, CALM1.
Bars: NicheNet AUROC (scRNA-seq; Mac S+M+S+ ligands vs genes up in OGN+RSPO3+ vs ADAMDEC1+ fibroblasts), top 20 ligands
present on the WTX panel; % = Mac S+M+S+ cells expressing the ligand in inflamed WTX (I0261, I0331).
Output: fig1i_ligand_activity.png."""
    import warnings; warnings.filterwarnings('ignore')
    import numpy as np, pandas as pd, scanpy as sc, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    O='/share/fsmresfiles/UC/revision_jenny/Fig2D_nichenet_ligand_activity'; F2='/share/fsmresfiles/UC/revision_jenny/Fig2E_I_J_ligand_regulator_cooccurrence_neighbourhood'; NN='/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip'
    ligand_activities=pd.read_csv(f'{NN}/ogn_spp1_mac_nichenet_ligand_activities.csv')
    T=pd.read_csv(f'{F2}/supp_table_fig2c_allsamples.csv'); T=T[T.Cohort=='Combined']
    validated=sorted(set(T.loc[T['P'].astype(float)<0.1,'Ligand (Mac S+M+S+)'])); print('highlighted:',validated)
    adata_wtx=sc.read_h5ad('/share/fsmresfiles/UC/AtoMx/UC_1k_6k_wholetrans/whole_trans/Processed_merged/v7/h5ad_obj/qc_wtx_v7_merged_anno.h5ad')
    adata_wtx=adata_wtx[adata_wtx.obs['patient'].isin(['I0261','I0331'])]
    top_ligands=ligand_activities.nlargest(50,'auroc')
    spp1_mac_mask=adata_wtx.obs['aucell_cell_type_short']=='Mac S+M+S+'
    ligand_in_macs=[]
    for ligand in top_ligands['test_ligand'].tolist():
        if ligand in adata_wtx.var_names:
            e=adata_wtx[spp1_mac_mask,ligand].X; e=e.toarray().flatten() if hasattr(e,'toarray') else np.asarray(e).flatten()
            ligand_in_macs.append({'ligand':ligand,'auroc':ligand_activities[ligand_activities['test_ligand']==ligand]['auroc'].values[0],
                                   'mean_expr_in_macs':e.mean(),'pct_macs_expressing':(e>0).sum()/len(e)*100})
    validated_ligands=pd.DataFrame(ligand_in_macs).sort_values('auroc',ascending=False)
    ligand_df=ligand_activities.sort_values('auroc',ascending=False)[:20]
    ligand_df=ligand_df.merge(validated_ligands,how='left',left_on='test_ligand',right_on='ligand')
    ligand_df=ligand_df.dropna(subset=['mean_expr_in_macs','pct_macs_expressing'])
    # ---- lollipop plot ----
    fig,ax=plt.subplots(figsize=(8.7,4.5)); x=list(range(len(ligand_df)))
    stem_colors=['#E63946' if lig in validated else '#A8DADC' for lig in ligand_df['test_ligand']]
    baseline=0.5
    for i,(idx,row) in enumerate(ligand_df.iterrows()):
        ax.plot([i,i],[baseline,row['auroc_x']],color=stem_colors[i],linewidth=4,alpha=0.85,solid_capstyle='round',zorder=2)
    ax.scatter(x,ligand_df['auroc_x'],s=150,c=stem_colors,edgecolors='black',linewidths=1,zorder=3)
    ax.axhline(baseline,color='gray',linestyle='--',alpha=0.5,linewidth=1.5,zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(ligand_df['test_ligand'],rotation=45,ha='right',fontsize=10)
    ax.set_xlabel('Mac S+M+S+ Ligands',fontsize=12); ax.set_ylabel('AUROC for OGN+RSPO3+ Fib DEGs',fontsize=12); ax.set_title('')
    ymax=float(ligand_df['auroc_x'].max()); ax.set_ylim([baseline-0.02,ymax+0.05])
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    for i,(idx,row) in enumerate(ligand_df.iterrows()):
        ax.text(i,row['auroc_x']+0.012,f"{row['pct_macs_expressing']:.1f}%",ha='center',va='bottom',fontsize=9,alpha=0.75,zorder=4)
    ax.text(-0.6,ymax+0.045,"% Mac S+M+S+ Expressing",ha='left',va='top',fontsize=10,alpha=0.85)
    legend_elements=[Line2D([0],[0],marker='o',color='w',markerfacecolor='#E63946',markersize=10,markeredgecolor='black',markeredgewidth=1,label='Spatially Supported (P < 0.1)'),
                     Line2D([0],[0],marker='o',color='w',markerfacecolor='#A8DADC',markersize=10,markeredgecolor='black',markeredgewidth=1,label='Predicted')]
    ax.legend(handles=legend_elements,loc='upper left',bbox_to_anchor=(0.72,1.0),fontsize=9.5,frameon=True,borderaxespad=0)
    plt.tight_layout(); fig.savefig(f'{O}/fig1i_ligand_activity.png',bbox_inches='tight',dpi=300); plt.close(fig)
    print(ligand_df[['test_ligand','auroc_x','pct_macs_expressing']].round(3).to_string(index=False)); print('SAVED')
