#!/usr/bin/env python
"""UC CosMx 6K: forest plot and supplementary table of the macrophage distance models (Fig 3H).

Input:   Fig2H_cosmx_fib_mac_distance/intermediate/glmm_final.csv (11) and ligand_sender_cells.csv.gz (07)
Outputs: fig_distance_forest.png (Fig 3H), supp_table_distance_models_cosmx.csv (all outcomes and models)
Run:     python <file>, or RUN_SECTIONS=<id>[,<id>] python <file>
"""
import os as _os
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# Section F3H_distance_forest: forest plot (Fig 3H) and supplementary table
if _run('F3H_distance_forest'):
    import warnings; warnings.filterwarnings('ignore')
    import pandas as pd, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    O='/share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance'
    C=pd.read_csv(f'{O}/intermediate/ligand_sender_cells.csv.gz'); G=pd.read_csv(f'{O}/intermediate/glmm_final.csv')
    GROUPS=[('SPP1 and MMP9','ogn',['SPP1','MMP9'],'#C1436D'),('IL33 Receptors','IL33',['IL1RL1','IL1RAP'],'#1C7C74'),
            ('CXCL12 Receptors','CXCL12',['CXCR4','ACKR3'],'#4F5FA0'),('FGF7 Receptors','FGF7',['FGFR1','FGFR2','NRP1'],'#D4952A')]   # IGF1 not plotted; reported in the CSV
    star=lambda q: 'ns' if pd.isna(q) else '***' if q<0.01 else '**' if q<0.05 else '*' if q<0.1 else 'ns'
    BLACK='#1b1b1f'
    plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':12,'axes.titleweight':'normal','figure.dpi':300,
        'axes.spines.top':False,'axes.spines.right':False,'text.color':BLACK,'axes.labelcolor':BLACK,'xtick.color':BLACK,'ytick.color':BLACK,'axes.edgecolor':BLACK})
    # (1) forest
    rows=[]; y=0; ylab=[]
    for grp,key,gs,col in GROUPS:
        for g in gs: rows.append((y,g,col)); ylab.append((y,g)); y+=1
        y+=0.7
    NS=C.groupby('condition')['sample'].nunique()
    fig,axes=plt.subplots(1,2,figsize=(8.6,4.6),sharey=True)
    for ax,(model,title) in zip(axes,[('UC',f"Ulcerative Colitis ({NS['UC']} Samples)"),('HC',f"Healthy Controls ({NS['Control']} Samples)")]):
        Gm=G[G.model==model].set_index('outcome')
        for yy,g,col in rows:
            r=Gm.loc[g]; ax.plot([r.OR_lo,r.OR_hi],[yy,yy],color=col,lw=2); ax.scatter(r.OR,yy,s=45,color=col,edgecolor=BLACK,lw=0.6,zorder=3)
            ax.text(r.OR_hi*1.07,yy,star(r.q),va='center',ha='left',fontsize=11 if star(r.q)!='ns' else 9)
            if g=='ACKR3': ax.text(0.32,yy,'decoy',va='center',ha='left',fontsize=9,style='italic',color=col)
        ax.axvline(1,color='#8A9399',ls='--',lw=1); ax.set_xscale('log'); ax.set_xticks([0.5,1,2]); ax.set_xticklabels(['0.5','1','2'])
        ax.set_xlim(0.3,4.6); ax.set_title(title,fontsize=12); ax.minorticks_off()
    axes[0].set_yticks([t for t,_ in ylab]); axes[0].set_yticklabels([l for _,l in ylab]); axes[0].invert_yaxis()
    for t,(yy,_g) in zip(axes[0].get_yticklabels(),ylab): t.set_color(next(c for y_,g_,c in rows if y_==yy))
    axes[0].set_ylabel('Gene in SELENOP⁺\nMacrophages',fontsize=11.5,linespacing=1.1)   # model uses all SELENOP+ macrophages (Methods)
    fig.supxlabel('Odds Ratio per 50 µm Closer to OGN+RSPO3+ Fibroblast Sender',fontsize=11.5,y=0.265,va='top')
    fig.suptitle('SELENOP+ Macrophage Expression vs Distance to Ligand-Expressing OGN+RSPO3+ Fibroblasts',fontsize=12.5,y=0.985)
    fig.legend(handles=[Line2D([],[],color=c,lw=3,label=g) for g,_,_,c in GROUPS[1:]],loc='upper center',bbox_to_anchor=(0.5,0.225),ncol=3,frameon=False,fontsize=10,handlelength=1.6,columnspacing=1.2)
    fig.subplots_adjust(left=0.17,right=0.97,top=0.86,bottom=0.34,wspace=0.12); fig.savefig(f'{O}/fig_distance_forest.png',bbox_inches='tight'); plt.close(fig)
    # (2) supplementary table
    lab=lambda g: 'Any receptor' if g.startswith('any_') else 'Response' if g.startswith('prog_') else 'ACKR3 (decoy)' if g=='ACKR3' else g
    SUP=str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    def qf(q):
        if pd.isna(q): return ''
        if q>=0.001: return f'{q:#.3g}'.rstrip('.')
        m,e=f'{q:.2e}'.split('e'); return f'{m} × 10{str(int(e)).translate(SUP)}'
    NP=G.groupby('model')['patients'].first()
    MN={'UC':f"Ulcerative colitis ({NP['UC']} patients)",'HC':f"Healthy controls ({NP['HC']} patients)",'ALL':f"UC + healthy combined ({NP['ALL']} patients; condition-adjusted)",'INT':f"Difference UC vs healthy (interaction, {NP['INT']} patients)"}
    rows=[{'Model':MN[r.model],'Family':r.family,'Gene or program':lab(r.outcome),
           ('OR per 50 µm closer (95% CI)' if r.model!='INT' else 'Ratio of ORs, UC vs healthy (95% CI)'):f'{r.OR:#.3g} ({r.OR_lo:#.3g}–{r.OR_hi:#.3g})','FDR q':qf(r.q)} for _,r in G.iterrows()]
    pd.DataFrame(rows).to_csv(f'{O}/supp_table_distance_models_cosmx.csv',index=False)
