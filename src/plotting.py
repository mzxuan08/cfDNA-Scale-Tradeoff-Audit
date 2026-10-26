"""Vector PDF and 300 DPI PNG. Functions only plot observed analysis outputs."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve


def save(fig, directory, name, dpi=300):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    plt.rcParams['pdf.fonttype'] = 42
    fig.savefig(directory / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(directory / f'{name}.png', dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def fig1(embedding, cfg):
    colors = cfg['plotting']['colors']
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for i, (center, g) in enumerate(embedding.groupby('center')):
        axes[0].scatter(g.PC1, g.PC2, s=16, alpha=.7, color=colors[i % len(colors)], label=center)
    for label, g in embedding.groupby('label'):
        axes[1].scatter(g.PC1, g.PC2, s=16, alpha=.7, label='Cancer' if label else 'Control')
    for ax, title in zip(axes, ['A  Center / declared study proxy', 'B  Clinical outcome']):
        ax.set(xlabel='PC1', ylabel='PC2', title=title)
        ax.legend(fontsize=7)
    save(fig, Path(cfg['paths']['results']) / 'figures', 'Fig1_batch_structure', cfg['plotting']['dpi'])


def fig2(summary, cfg):
    g = summary[(summary.method == 'raw') & summary.status.eq('complete')]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, metric, title in zip(axes, ['auroc', 'macro_center_auroc'], ['A  Pooled OOF AUROC', 'B  Macro center AUROC']):
        for mode, m in g.groupby('mode'):
            ax.plot(m.model, m[metric], 'o-', label=mode)
        ax.set(ylabel='AUROC', ylim=(0, 1), title=title)
        ax.tick_params(axis='x', rotation=35)
        ax.legend()
    save(fig, Path(cfg['paths']['results']) / 'figures', 'Fig2_generalization', cfg['plotting']['dpi'])


def fig3(summary, predictions, cfg):
    complete = summary[(summary['mode'] == 'loco') & summary.status.eq('complete')]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for method, g in complete.groupby('method'):
        axes[0].plot(g.model, g.auroc, 'o-', label=method)
    axes[0].set(ylabel='LOCO OOF AUROC', ylim=(0, 1), title='A  Prespecified ablations')
    axes[0].tick_params(axis='x', rotation=35)
    if len(complete):
        axes[0].legend(fontsize=7)
    axes[1].plot([0, 1], [0, 1], '--', color='gray')
    # Model chosen a priori, never chosen by best test performance.
    for method in sorted(complete.loc[complete.model.eq('elasticnet'), 'method'].unique()):
        p = predictions[(predictions['mode'] == 'loco') & (predictions['model'] == 'elasticnet') & (predictions['method'] == method)]
        p = p.groupby('patient_id', as_index=False).agg(label=('label', 'first'), probability=('probability', 'mean'))
        observed, predicted = calibration_curve(p.label, p.probability, n_bins=10, strategy='uniform')
        axes[1].plot(predicted, observed, 'o-', label=method)
    axes[1].set(xlabel='Predicted probability', ylabel='Observed cancer frequency', title='B  ElasticNet calibration', xlim=(0, 1), ylim=(0, 1))
    if axes[1].get_legend_handles_labels()[0]:
        axes[1].legend(fontsize=7)
    save(fig, Path(cfg['paths']['results']) / 'figures', 'Fig3_mitigation', cfg['plotting']['dpi'])


def fig4(table, cfg):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    g = table[table.feature.str.startswith('ratio__')].sort_values('feature')
    axes[0].scatter(np.arange(len(g)), g.disease_beta, c=cfg['plotting']['colors'][0], s=8)
    axes[0].axhline(0, color='gray', linewidth=.5)
    axes[0].set(xlabel='Genomic bin (lexical order)', ylabel='Adjusted disease coefficient', title='A  Exploratory disease association')
    if 'center_fraction' in table:
        axes[1].hist(table.center_fraction.dropna(), bins=20, color=cfg['plotting']['colors'][1])
    axes[1].set(xlabel='Conditional center variance fraction', ylabel='Features', title='B  Center variance distribution')
    save(fig, Path(cfg['paths']['results']) / 'figures', 'Fig4_exploratory_biology', cfg['plotting']['dpi'])


def short_center(value):
    return value.replace(' et al, ',' ')


def major_revision_figures(analysis):
    """Read revision outputs and save new figures; frozen v1 panels stay intact."""
    from .revision_plotting import revision_figures
    return revision_figures(analysis)


def clinical_style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,
        'axes.labelsize':8,'legend.fontsize':7,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})


def clinical_figure1(meta,ratio,lengths,bins,embedding,perma,partition,drift,cfg):
    clinical_style()
    colors=cfg['plotting']['colors'];centers=sorted(meta.center.unique())
    fig=plt.figure(figsize=(13,10),layout='constrained');outer=fig.add_gridspec(2,2)
    a=outer[0,0].subgridspec(2,1,height_ratios=[1,2]);composition=fig.add_subplot(a[0]);curve=fig.add_subplot(a[1])
    counts=pd.crosstab(meta.center,meta.label).reindex(centers)
    composition.barh(range(3),counts[0],color=colors[1],label='Healthy')
    composition.barh(range(3),counts[1],left=counts[0],color=colors[0],label='Cancer')
    composition.set(yticks=range(3),yticklabels=[short_center(c) for c in centers],xlabel='Independent participants',title='a  Cohort and healthy fragment profiles')
    composition.legend(loc='lower right',ncol=2)
    for i,c in enumerate(centers):
        mask=meta.center.eq(c)&meta.label.eq(0)
        curve.plot(np.arange(80,401),lengths[mask].mean(axis=0),color=colors[i],label=f'{short_center(c)} (n={mask.sum()})')
    for length in [143,167]:curve.axvline(length,color='gray',ls=':',lw=.8)
    curve.set(xlim=(100,220),xlabel='Fragment length (bp)',ylabel='Global frequency');curve.legend()
    b=outer[0,1].subgridspec(1,2);axes=[fig.add_subplot(b[0]),fig.add_subplot(b[1])]
    e=embedding[embedding.family.eq('raw_ratio')]
    for i,c in enumerate(centers):
        g=e[e.center.eq(c)];axes[0].scatter(g.UMAP1,g.UMAP2,s=9,color=colors[i],alpha=.65,label=short_center(c))
    for y,color in [(0,colors[1]),(1,colors[0])]:
        g=e[e.label.eq(y)];axes[1].scatter(g.UMAP1,g.UMAP2,s=9,color=color,alpha=.65,label='Cancer' if y else 'Healthy')
    for ax,term,title in zip(axes,['center','disease'],['Study','Disease']):
        row=perma[perma.family.eq('raw_ratio')&perma.term.eq(term)].iloc[0]
        ax.set(title=f'{title}: unique Type III R²={row.R2:.3f}; p={row.p:.3f}',xlabel='UMAP1',ylabel='UMAP2');ax.legend(loc='best')
    cax=fig.add_subplot(outer[1,0]);g=partition[partition.family.eq('raw_ratio')]
    offset=0;ticks=[];labels=[]
    for k,(chrom,group) in enumerate(g.groupby('chrom',sort=False)):
        size=int(group.chrom_size_bp.iloc[0]);pos=(group.midpoint_bp+offset)/1e9
        cax.scatter(pos,group.eta2_batch*100,s=6,color=colors[0],label='Center' if k==0 else None)
        cax.scatter(pos,group.eta2_bio*100,s=6,color=colors[1],label='Disease' if k==0 else None)
        ticks.append((offset+size/2)/1e9);labels.append(chrom[3:]);offset+=size
    cax.set(title='c  Unique variance by genomic window',ylabel='Type III eta² (% of total SS)',xticks=ticks,xticklabels=labels,xlabel='Autosome');cax.tick_params(axis='x',labelsize=6);cax.legend()
    d=outer[1,1].subgridspec(1,2);box=fig.add_subplot(d[0]);density=fig.add_subplot(d[1])
    values=[np.median(ratio[meta.center.eq(c)&meta.label.eq(0)],axis=1) for c in centers]
    box.boxplot(values,tick_labels=[short_center(c) for c in centers],showfliers=False)
    box.tick_params(axis='x',rotation=35);box.set(title='d  Healthy baseline drift',ylabel='Within-sample median raw ratio')
    density.hist(drift.center_icc1,bins=25,color=colors[2],alpha=.8,density=True)
    density.set(title=f'Healthy center ICC(1)\nFDR<0.01: {drift.q.lt(.01).mean():.1%} bins',xlabel='Center clustering ICC(1)',ylabel='Density')
    save(fig,Path(cfg['paths']['results'])/'figures','Figure1_Batch_Vulnerability_Audit',300)


def clinical_figures23(meta,R,bins,pred,metrics,significance,ablation,center_diagnostic,cfg):
    from sklearn.metrics import roc_curve,precision_recall_curve
    from scipy.stats import beta
    from .debiasing import ArmLocalLogRatio
    from .study import PIPELINES,MODELS
    clinical_style();colors=cfg['plotting']['colors'];destination=Path(cfg['paths']['results'])/'figures'
    allm=metrics[metrics.center.eq('AllCenters')]
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    raw=allm[allm.pipeline.eq('Raw_Ratio')].pivot(index='model',columns='mode',values='macro_center_auroc').reindex(MODELS)
    for i,(model,row) in enumerate(raw.iterrows()):
        axes[0,0].plot([row['random'],row['loco']],[i,i],color='gray',lw=1)
        axes[0,0].scatter(row['random'],i,color=colors[0],label='Random-CV' if i==0 else None)
        axes[0,0].scatter(row['loco'],i,color=colors[1],label='LOCO' if i==0 else None)
    axes[0,0].set(yticks=range(6),yticklabels=MODELS,xlabel='Macro study AUROC',xlim=(0,1),title='a  Random versus unseen-study performance');axes[0,0].legend()
    groups=[g.balanced_accuracy.to_numpy() for _,g in center_diagnostic.groupby('scope',sort=False)]
    axes[0,1].boxplot(groups,tick_labels=['All subjects (3 studies)','Healthy (2 studies)'])
    axes[0,1].set(ylim=(0,1),ylabel='RF study-prediction balanced accuracy',title='b  Study information in fragmentation features')
    for i,pipeline in enumerate(PIPELINES):
        g=metrics[metrics['mode'].eq('loco')&metrics.center.ne('AllCenters')&metrics.pipeline.eq(pipeline)&metrics.model.eq('elasticnet')]
        locations=np.arange(len(g))+(i-1.5)*.06
        lows=[];highs=[]
        for row in g.itertuples():
            controls=pred[pred['mode'].eq('loco')&pred.pipeline.eq(pipeline)&pred.model.eq('elasticnet')&pred.center.eq(row.center)&pred.label.eq(0)]
            n=len(controls);success=int((controls.probability<=controls.threshold).sum())
            lows.append(float(beta.ppf(.025,success,n-success+1)) if success else 0.)
            highs.append(float(beta.ppf(.975,success+1,n-success)) if success<n else 1.)
        axes[1,0].errorbar(locations,g.specificity95,yerr=[g.specificity95.to_numpy()-lows,np.array(highs)-g.specificity95.to_numpy()],fmt='o',color=colors[i],capsize=2,label=pipeline,alpha=.8)
    axes[1,0].set_xticks(range(len(g)),[short_center(c) for c in g.center])
    axes[1,0].axhline(.95,ls='--',color='gray');axes[1,0].set(ylim=(0,1.03),ylabel='Achieved test specificity',title='c  Source-locked 95% specificity threshold');axes[1,0].legend(fontsize=6)
    h=metrics[metrics['mode'].eq('loco')&metrics.pipeline.eq('Raw_Ratio')&metrics.center.ne('AllCenters')].pivot(index='model',columns='center',values='auroc').reindex(MODELS)
    image=axes[1,1].imshow(h,aspect='auto',vmin=0,vmax=1,cmap='viridis')
    axes[1,1].set(yticks=range(len(h)),yticklabels=h.index,xticks=range(len(h.columns)),xticklabels=[short_center(c) for c in h.columns],title='d  Held-out study AUROC (Raw Ratio)')
    for i in range(len(h)):
        for j in range(len(h.columns)):axes[1,1].text(j,i,f'{h.iloc[i,j]:.3f}',ha='center',va='center',color='white' if h.iloc[i,j]<.8 else 'black',fontsize=7)
    fig.colorbar(image,ax=axes[1,1],shrink=.7,label='AUROC')
    save(fig,destination,'Figure2_Random_vs_LOCO_Collapse',300)
    fig=plt.figure(figsize=(14,10),layout='constrained');outer=fig.add_gridspec(2,2)
    a=outer[0,0].subgridspec(2,1,height_ratios=[1,3]);diagram=fig.add_subplot(a[0]);diagram.axis('off')
    diagram.text(.5,.6,'S/(L+1)  →  same-arm log₂ flank contrast  →  train-healthy center filter  →  train RobustScaler',ha='center',va='center',wrap=True,fontsize=9)
    track=fig.add_subplot(a[1]);lsn=ArmLocalLogRatio(bins,3).fit(R).transform(R)
    for i,c in enumerate(sorted(meta.center.unique())):
        track.plot(np.nanmedian(lsn[meta.center.eq(c)&meta.label.eq(0)],axis=0),color=colors[i],lw=.8,label=short_center(c))
    track.set(title='a  LSN pipeline and observed healthy baselines',xlabel='Genomic window order',ylabel='Median healthy LSN');track.legend()
    b=outer[0,1].subgridspec(1,2);rocax=fig.add_subplot(b[0]);prax=fig.add_subplot(b[1])
    for i,pipeline in enumerate(PIPELINES):
        g=pred[pred['mode'].eq('loco')&pred.model.eq('elasticnet')&pred.pipeline.eq(pipeline)]
        fpr,tpr,_=roc_curve(g.label,g.probability);precision,recall,_=precision_recall_curve(g.label,g.probability)
        row=allm[allm['mode'].eq('loco')&allm.model.eq('elasticnet')&allm.pipeline.eq(pipeline)].iloc[0]
        rocax.plot(fpr,tpr,color=colors[i],label=f'{pipeline}: {row.auroc:.3f}')
        prax.plot(recall,precision,color=colors[i],label=f'{pipeline}: {row.auprc:.3f}')
    rocax.plot([0,1],[0,1],'--',color='gray',lw=.5);rocax.set(title='b  Pooled LOCO ROC',xlabel='FPR',ylabel='TPR');rocax.legend(fontsize=5)
    prax.set(title='LOCO precision–recall',xlabel='Recall',ylabel='Precision');prax.legend(fontsize=5)
    sig=significance[(significance.test=='paired_DeLong')&(significance.model=='elasticnet')&(significance.comparator=='Raw_Ratio')]
    label='Held-out DeLong p: '+', '.join(f'{short_center(r.center)} {r.p:.3g}' if np.isfinite(r.p) else f'{short_center(r.center)} NA (1 HC)' for r in sig.itertuples())
    rocax.text(0,-.25,label,transform=rocax.transAxes,fontsize=6,wrap=True)
    cal=fig.add_subplot(outer[1,0]);cal.plot([0,1],[0,1],'--',color='gray')
    for i,pipeline in enumerate(PIPELINES):
        g=pred[pred['mode'].eq('loco')&pred.model.eq('elasticnet')&pred.pipeline.eq(pipeline)]
        observed,expected=calibration_curve(g.label,g.probability,n_bins=10)
        row=allm[allm['mode'].eq('loco')&allm.model.eq('elasticnet')&allm.pipeline.eq(pipeline)].iloc[0]
        cal.plot(expected,observed,'o-',color=colors[i],label=f'{pipeline}; Brier={row.brier:.3f}')
    cal.set(title='c  ElasticNet calibration',xlabel='Predicted cancer probability',ylabel='Observed proportion',xlim=(0,1),ylim=(0,1));cal.legend(fontsize=6)
    surface=fig.add_subplot(outer[1,1])
    grid=ablation[(ablation.center=='MacroCenters')&(ablation.pipeline=='Proposed_LSN_ICC')].pivot(index='alpha',columns='radius',values='auroc')
    im=surface.imshow(grid,origin='lower',aspect='auto',vmin=.5,vmax=1,cmap='viridis')
    surface.set(title='d  Prespecified LOCO sensitivity (ElasticNet)',xticks=range(len(grid.columns)),xticklabels=grid.columns,
                yticks=range(len(grid.index)),yticklabels=[f'{v:.0%}' for v in grid.index],xlabel='Flank radius k',ylabel='Removed healthy-batch fraction α')
    for i in range(len(grid)):
        for j in range(len(grid.columns)):surface.text(j,i,f'{grid.iloc[i,j]:.3f}',ha='center',va='center',color='white',fontsize=7)
    fig.colorbar(im,ax=surface,shrink=.7,label='Macro study AUROC')
    save(fig,destination,'Figure3_LSN_Debiasing_Performance',300)


def clinical_figure4(data,correlation,enrichment,meta,R,L,bins,importance,cfg):
    from statsmodels.nonparametric.smoothers_lowess import lowess
    clinical_style();colors=cfg['plotting']['colors']
    fig=plt.figure(figsize=(14,10),layout='constrained');outer=fig.add_gridspec(2,2)
    a=outer[0,0].subgridspec(1,2);axgc=fig.add_subplot(a[0]);axt=fig.add_subplot(a[1])
    for ax,variable,label in [(axgc,'gc_fraction','GC fraction'),(axt,'log_telomere_distance','log(1 + same-arm telomere distance)')]:
        x=data[variable].to_numpy();y=data.eta2_batch.to_numpy()*100
        ax.scatter(x,y,s=7,color=colors[1],alpha=.6)
        smooth=lowess(y,x,frac=.4,return_sorted=True);ax.plot(smooth[:,0],smooth[:,1],color=colors[0],lw=1.3)
        row=correlation[correlation.annotation.eq(variable)].iloc[0]
        ax.set(xlabel=label,ylabel='Unique center eta² (%)',title=f'{"a  " if variable=="gc_fraction" else ""}rho={row.rho:.3f}; p={row.p:.2g}')
    forest=fig.add_subplot(outer[0,1]);ys=np.arange(len(enrichment))
    forest.errorbar(enrichment.odds_ratio,ys,xerr=[enrichment.odds_ratio-enrichment.ci_low,enrichment.ci_high-enrichment.odds_ratio],fmt='o',color=colors[0],capsize=2)
    forest.axvline(1,color='gray',ls='--');forest.set(xscale='log',yticks=ys,
        yticklabels=[s.replace('_',' ') for s in enrichment.region],xlabel='Odds ratio (95% interval)',title='b  Majority-fold removed versus retained windows')
    c=outer[1,0].subgridspec(2,1);rawax=fig.add_subplot(c[0]);lsnax=fig.add_subplot(c[1])
    chrom='chr8';mask=bins.chrom.eq(chrom).to_numpy();positions=bins.loc[mask,'midpoint_bp'].to_numpy()/1e6
    breaks=np.flatnonzero(np.diff(positions)>5)+1
    track_positions=np.insert(positions,breaks,np.nan)
    for i,center in enumerate(sorted(meta.center.unique())):
        healthy=meta.center.eq(center)&meta.label.eq(0);cancer=meta.center.eq(center)&meta.label.eq(1)
        for ax,x in [(rawax,R),(lsnax,L)]:
            ax.plot(track_positions,np.insert(np.median(x[healthy][:,mask],axis=0),breaks,np.nan),color=colors[i],lw=1,label=short_center(center)+' healthy')
            ax.plot(track_positions,np.insert(np.median(x[cancer][:,mask],axis=0),breaks,np.nan),color=colors[i],ls='--',alpha=.65,lw=.9)
    rawax.set(title='c  chr8: observed healthy (solid) and cancer (dashed)',ylabel='Raw S/(L+1)');rawax.legend(fontsize=6,ncol=2)
    lsnax.set(xlabel='chr8 coordinate (Mb)',ylabel='Local log₂ ratio')
    heat=fig.add_subplot(outer[1,1])
    top=importance.groupby('bin_id').importance.sum().nlargest(20).index
    ranks=importance.pivot(index='bin_id',columns='fold',values='rank').reindex(top)
    cmap=plt.get_cmap('viridis_r').copy();cmap.set_bad('lightgray')
    im=heat.imshow(ranks,aspect='auto',cmap=cmap,vmin=1,vmax=len(bins))
    heat.set(yticks=range(len(ranks)),yticklabels=ranks.index,xticks=[0,1,2],xticklabels=['Held out Cristiano','Held out Jiang','Held out Snyder'],title='d  Top 20 ElasticNet coefficient rank by LOCO fold')
    heat.tick_params(axis='y',labelsize=6);heat.tick_params(axis='x',labelsize=6,rotation=15)
    fig.colorbar(im,ax=heat,shrink=.7,label='Rank (gray: filtered out)')
    save(fig,Path(cfg['paths']['results'])/'figures','Figure4_Genomic_Interpretability',300)
