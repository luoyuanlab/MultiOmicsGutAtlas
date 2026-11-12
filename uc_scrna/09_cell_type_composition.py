#!/usr/bin/env python
"""UC scRNA-seq atlas: cell-type composition across healthy, non-inflamed and inflamed UC.

Panels: Fig 2E (macrophage subtype composition per patient; Kruskal-Wallis, BH over 6 subtypes);
        Supp 3D (scCODA differential abundance, est. FDR 0.1) and its supp table.
Inputs: UC atlas h5ad (06_atlas_assembly.R).
Outputs: SuppFig3D_sccoda/ (effects CSV, panel, supp table); Fig1E_macrophage_composition/ (panel, per-patient CSV).
Run: python <file>, or RUN_SECTIONS=<id>[,<id>...] python <file>. S3D_sccoda runs in the pertpy 1.0.2 environment.
"""
import os as _os
import numpy as np
import pandas as pd
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# Section S3D_sccoda: scCODA of UC atlas cell types, studies P1/P4/P5/P8, reference Cycl Epi, est. FDR 0.1 (Supp 3D; pertpy env)
if _run('S3D_sccoda'):
    import sys, types, os
    if 'platformdirs' not in sys.modules:   # stub platformdirs if absent (only used by a plotting dependency)
        try: import platformdirs  # noqa: F401
        except ImportError:
            import pathlib; m = types.ModuleType('platformdirs'); m.user_cache_dir = lambda *a, **k: os.path.expanduser('~/.cache'); m.user_cache_path = lambda *a, **k: pathlib.Path(os.path.expanduser('~/.cache')); sys.modules['platformdirs'] = m
    import scanpy as sc, pertpy as pt
    O = '/share/fsmresfiles/UC/revision_jenny/SuppFig3D_sccoda'
    scrna = sc.read_h5ad('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad')
    scrna = scrna[scrna.obs['orig.ident'].isin(['P1', 'P4', 'P5', 'P8'])]
    scc = pt.tl.Sccoda()
    mdata = scc.load(scrna, type='cell_level', generate_sample_level=True, cell_type_identifier='cell_type_short',
                     sample_identifier='sample', covariate_obs=['condition'])
    mdata = scc.prepare(mdata, formula='C(condition)', reference_cell_type='Cycl Epi')
    scc = pt.tl.Sccoda()
    scc.run_nuts(mdata, rng_key=1234)
    scc.set_fdr(mdata, est_fdr=0.1)
    scc.get_effect_df(mdata).to_csv(f'{O}/sccoda_effects_rerun_fdr01.csv')


# Section S3D_differential_abundance_plot: scCODA effects per cell type, inflamed and non-inflamed vs healthy (HDI 3-97%) (Supp 3D)
if _run('S3D_differential_abundance_plot'):
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    O = '/share/fsmresfiles/UC/revision_jenny/SuppFig3D_sccoda'
    HIGHLIGHT = ['Activ Fib', 'CXCL5+ Activ Fib', 'Inf Mo-Mac', 'Mo-Mac', 'Mac S+M+S+']          # boxed cell types
    df = pd.read_csv(f'{O}/sccoda_effects_rerun_fdr01.csv')
    cov_inflamed, cov_noninflamed = 'C(condition)T.Inflamed', 'C(condition)T.Non-inflamed'
    order = df[df.Covariate == cov_inflamed].sort_values('HDI 97%', ascending=True)['Cell Type'].tolist()
    piv = df.pivot(index='Cell Type', columns='Covariate', values=['Final Parameter', 'HDI 3%', 'HDI 97%']).reindex(order)
    y = np.arange(len(piv.index)); offset = 0.15
    arr = lambda m, c: piv[(m, c)].to_numpy()
    xerr = lambda e, lo, hi: np.vstack([e - lo, hi - e])
    fig, ax = plt.subplots(figsize=(6, 0.27 * len(piv.index)))
    for cov, sgn, col, lab in [(cov_inflamed, -1, '#ab210c', 'Inflamed'), (cov_noninflamed, 1, '#a88c0f', 'Non-inflamed')]:
        e = arr('Final Parameter', cov)
        ax.errorbar(e, y + sgn * offset, xerr=xerr(e, arr('HDI 3%', cov), arr('HDI 97%', cov)), fmt='o', capsize=2, elinewidth=1,
                    markersize=4, color=col, ecolor=col, alpha=0.9, label=lab)
    ax.axvline(0, linestyle='--', linewidth=1)
    ax.set_yticks(y); ax.set_yticklabels(piv.index)
    ax.set_ylim(-0.5 - offset, (len(piv.index) - 0.5) + offset); ax.invert_yaxis(); ax.margins(y=0)
    ax.set_xlabel('Final parameter (3–97% HDI)'); ax.set_title('Differential Abundance'); ax.legend(frameon=False, loc='upper right')
    fig.tight_layout(); fig.canvas.draw()
    inv = fig.transFigure.inverted()
    for lbl in ax.get_yticklabels():
        if lbl.get_text() in HIGHLIGHT:
            bb = lbl.get_window_extent().expanded(1.06, 1.25); (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
            fig.patches.append(Rectangle((x0, y0), x1 - x0, y1 - y0, transform=fig.transFigure, fill=False, edgecolor='#9b1b1b', lw=1.4))
    fig.savefig(f'{O}/ed3d_differential_abundance_fdr01.png', dpi=300, bbox_inches='tight')


# Section S3D_sccoda_table: supp table of the scCODA effects at est. FDR 0.1 (3 significant figures)
if _run('S3D_sccoda_table'):
    O = '/share/fsmresfiles/UC/revision_jenny/SuppFig3D_sccoda'
    b = pd.read_csv(f'{O}/sccoda_effects_rerun_fdr01.csv')
    b['Covariate'] = b['Covariate'].str.replace('C(condition)T.', '', regex=False)
    out = b[['Covariate', 'Cell Type', 'Final Parameter', 'HDI 3%', 'HDI 97%', 'Inclusion probability', 'log2-fold change']].copy()
    for c in out.columns[2:]: out[c] = out[c].map(lambda x: float(f'{x:.3g}'))
    out['Credible (est. FDR 0.1)'] = out['Final Parameter'] != 0
    out.rename(columns={'Covariate': 'Condition vs healthy'}).to_csv(f'{O}/supp_table_sccoda_uc_fdr01.csv', index=False)


# Section F2E_macrophage_subtype_composition: macrophage subtype composition (% of macrophages) per patient x condition (>= 20 macs); Kruskal-Wallis per subtype, BH over 6 (Fig 2E)
if _run('F2E_macrophage_subtype_composition'):
    import warnings; warnings.filterwarnings('ignore')
    import anndata as ad, matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from scipy.stats import kruskal
    from statsmodels.stats.multitest import multipletests
    O='/share/fsmresfiles/UC/revision_jenny/Fig1E_macrophage_composition'
    obs=ad.read_h5ad('/share/fsmresfiles/UC/scRNA-seq/merged_9p_wip/pdated_annoed_all_cell_types.h5ad',backed='r').obs[['label','Patient','condition']].astype(str)
    MAC=['SELENOP+SIGLEC1+ Mac','SELENOP+XIST- Mac','SELENOP+MMP9+PLA2G2D+ Mac','SELENOP+MMP9+SPP1+ Mac','Mo-Mac','Inflammatory Mo-Mac']
    SHORT={'SELENOP+SIGLEC1+ Mac':'Mac S+SG+','SELENOP+XIST- Mac':'Mac S+XS-','SELENOP+MMP9+PLA2G2D+ Mac':'Mac S+M+P+',
           'SELENOP+MMP9+SPP1+ Mac':'Mac S+M+S+','Mo-Mac':'Mo-Mac','Inflammatory Mo-Mac':'Inf Mo-Mac'}
    G=['Healthy','Non-inflamed','Inflamed']
    m=obs[obs.label.isin(MAC)&obs.condition.isin(G)]
    ct=pd.crosstab([m.Patient,m.condition],m.label); ct=ct[ct.sum(1)>=20]; fr=ct.div(ct.sum(1),axis=0)*100
    fr.reset_index().to_csv(f'{O}/macrophage_composition_per_patient.csv',index=False)
    p={l:kruskal(*[fr.xs(g,level='condition')[l] for g in G]).pvalue for l in MAC}; q=dict(zip(p,multipletests(list(p.values()),method='fdr_bh')[1]))
    star=lambda v: '***' if v<0.01 else '**' if v<0.05 else '*' if v<0.1 else 'ns'
    COL={'Healthy':'#4F9D69','Non-inflamed':'#E0A93B','Inflamed':'#C1436D'}; BLACK='#1b1b1f'
    plt.rcParams.update({'font.size':11,'figure.dpi':300,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(9,4.2)); x0=np.arange(len(MAC))*1.0; w=0.26
    for j,g in enumerate(G):
        for i,l in enumerate(MAC):
            v=fr.xs(g,level='condition')[l].values; xp=x0[i]+(j-1)*w
            bp=ax.boxplot(v,positions=[xp],widths=w*0.85,patch_artist=True,showfliers=False,medianprops=dict(color=BLACK),
                          boxprops=dict(edgecolor=BLACK),whiskerprops=dict(color=BLACK),capprops=dict(color=BLACK))
            bp['boxes'][0].set_facecolor(COL[g]); bp['boxes'][0].set_alpha(0.55)
            ax.scatter(np.random.default_rng(10*i+j).normal(xp,0.035,len(v)),v,s=6,color=BLACK,alpha=0.6,lw=0,zorder=3)
    for i,l in enumerate(MAC):
        top=fr[l].max(); ax.text(x0[i],top+3,star(q[l]),ha='center',va='bottom',fontsize=11 if star(q[l])!='ns' else 9)
    ax.set_xticks(x0); ax.set_xticklabels([SHORT[l] for l in MAC]); ax.set_ylabel('% of Macrophages (per patient)'); ax.set_ylim(0,108)
    ax.set_title('Macrophage Composition by Condition (UC scRNA-seq Atlas)',fontsize=12)
    n={g:len(fr.xs(g,level='condition')) for g in G}
    ax.legend(handles=[plt.Rectangle((0,0),1,1,facecolor=COL[g],alpha=0.55,edgecolor=BLACK,label=f'{g} (n = {n[g]})') for g in G],
              loc='upper center',bbox_to_anchor=(0.5,-0.12),ncol=3,frameon=False,fontsize=10)
    ax.text(0.99,0.98,'Kruskal–Wallis, BH-FDR: * q<0.1, ** q<0.05, *** q<0.01',transform=ax.transAxes,ha='right',va='top',fontsize=8.5,color='#777777')
    plt.tight_layout(); fig.savefig(f'{O}/fig_macrophage_composition.png',bbox_inches='tight'); plt.close(fig)
