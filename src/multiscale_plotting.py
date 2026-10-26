"""Step 7 figures: distinguish bilateral inference from sensitivity-only follow-up."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve,precision_recall_curve
from sklearn.calibration import calibration_curve
from .plotting import clinical_style,save,short_center
from .multiscale_benchmark import PRIMARY,COMBINED,EXTERNAL,TUNED,read_protocol


def step7_figures(cfg,pred,metrics,significance,counts,pairs,stress):
    clinical_style();protocol=read_protocol(cfg);colors=cfg['plotting']['colors']
    pipelines=protocol['pipelines']+[TUNED]
    color={p:colors[i] for i,p in enumerate(pipelines)}
    labels={'Raw_Ratio':'Raw Ratio','Global_Zscore':'Global Z-score','InFold_ComBat':'In-fold ComBat',
            'Proposed_LSN_ICC':'Fixed LSN (ICC bypassed*)','Hybrid_MultiScale':'Hybrid MultiScale',TUNED:'Nested-tuned LSN*'}
    folder=Path(cfg['paths']['results'])/'figures'
    primary=metrics[metrics.cohort_scope.eq(PRIMARY)]
    fig,axes=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    primary_n=int(counts[counts.center.isin(protocol['primary_studies'])].n.sum())
    fig.suptitle(f'Primary Bilateral LOCO Benchmark (n={primary_n}); separate external sensitivity cohort',fontsize=12)
    for offset,pipeline in [(-.13,'Raw_Ratio'),(.13,'Hybrid_MultiScale')]:
        table=primary[primary.pipeline.eq(pipeline)&primary.center.eq('AllCenters')].pivot(index='model',columns='mode',values='macro_bilateral_auroc')
        for i,model in enumerate(protocol['models']):
            if model not in table.index or not np.isfinite(table.loc[model,'loco']):continue
            a,b=table.loc[model,['random','loco']]
            axes[0,0].plot([a,b],[i+offset,i+offset],color=color[pipeline],alpha=.6)
            axes[0,0].scatter(a,i+offset,marker='o',facecolors='none',edgecolors=color[pipeline],s=25)
            axes[0,0].scatter(b,i+offset,marker='o',color=color[pipeline],s=25,label=labels[pipeline] if i==0 else None)
    axes[0,0].set(yticks=range(len(protocol['models'])),yticklabels=protocol['models'],xlim=(0,1),
        xlabel='Mean AUROC of the two primary target studies',title='a  Primary cohort: random versus bilateral LOCO')
    axes[0,0].text(.02,.02,'Hollow: random CV; solid: bilateral LOCO\nDANN unavailable with one source study',transform=axes[0,0].transAxes,fontsize=7)
    axes[0,0].legend(fontsize=7)
    table=primary[primary['mode'].eq('loco')&primary.center.eq('AllCenters')].pivot(index='pipeline',columns='model',values='macro_bilateral_auroc').reindex(index=pipelines,columns=protocol['models'])
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#EEEEEE')
    im=axes[0,1].imshow(table,aspect='auto',vmin=0,vmax=1,cmap=cmap)
    axes[0,1].set(yticks=range(len(table)),yticklabels=[labels[p] for p in table.index],xticks=range(len(table.columns)),xticklabels=table.columns,
        title='b  Bilateral LOCO: macro-study AUROC')
    for i in range(len(table)):
        for j in range(len(table.columns)):
            v=table.iloc[i,j];axes[0,1].text(j,i,f'{v:.3f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=7,
                                                color='black' if not np.isfinite(v) or v>.72 else 'white')
    axes[0,1].tick_params(axis='x',labelsize=7);axes[0,1].tick_params(axis='y',labelsize=7)
    fig.colorbar(im,ax=axes[0,1],shrink=.7,label='AUROC')
    for i,pipeline in enumerate(pipelines):
        group=primary[primary['mode'].eq('loco')&primary.model.eq('elasticnet')&primary.pipeline.eq(pipeline)&primary.center.ne('AllCenters')]
        group=group[group.status.eq('complete')].sort_values('center')
        if group.empty:continue
        x=np.arange(len(group))+(i-2.5)*.055
        axes[1,0].errorbar(x,group.specificity_deviation_pp,
            yerr=[(group.specificity95-group.specificity_exact_low)*100,(group.specificity_exact_high-group.specificity95)*100],
            fmt='o',capsize=2,color=color[pipeline],label=labels[pipeline])
    axes[1,0].axhline(0,color='gray',ls='--')
    axes[1,0].set(xticks=[0,1],xticklabels=[short_center(c) for c in protocol['primary_studies']],
        ylabel='Specificity minus source target 95% (percentage points)',title='c  Primary ElasticNet threshold transfer (exact 95% CI)')
    axes[1,0].legend(fontsize=6,loc='lower right')
    ext=metrics[metrics.cohort_scope.eq(EXTERNAL)&metrics.model.eq('elasticnet')].set_index('pipeline').reindex(pipelines)
    ys=np.arange(len(ext));axes[1,1].errorbar(ext.sensitivity95,ys,
        xerr=[ext.sensitivity95-ext.sensitivity_exact_low,ext.sensitivity_exact_high-ext.sensitivity95],fmt='o',color=colors[4],capsize=3)
    axes[1,1].set(yticks=ys,yticklabels=[labels[p].replace(' (ICC bypassed*)','').replace('*','') for p in ext.index],
        xlim=(0,1.02),xlabel='Cancer sensitivity at primary-source-locked threshold',title='d  External high-burden cohort: sensitivity focus')
    axes[1,1].text(.02,.04,'Snyder: 36 cancer / 1 healthy\nSpecificity has no stable population estimate\n*ICC bypass applies only to one-source bilateral folds',transform=axes[1,1].transAxes,fontsize=7)
    save(fig,folder,'Figure2_Random_vs_LOCO_Collapse',300)

    fig=plt.figure(figsize=(14,11),layout='constrained');outer=fig.add_gridspec(2,2)
    fig.suptitle(f'Multiscale validation in the primary bilateral cohort (n={primary_n})',fontsize=12)
    a=outer[0,0].subgridspec(2,1,height_ratios=[1,2]);diagram=fig.add_subplot(a[0]);diagram.axis('off')
    diagram.text(.5,.8,'a  Sample-wise two-scale representation',ha='center',fontsize=11)
    diagram.text(.5,.44,'Global window Z-scores  +  same-arm local log contrasts\nTrain-only column scaling → concatenation → ElasticNet',ha='center',va='center',fontsize=8)
    diagram.text(.5,.04,'Local radius is fixed for Hybrid; nested tuning is a separate LSN analysis',ha='center',fontsize=7)
    bars=fig.add_subplot(a[1]);g=primary[primary['mode'].eq('loco')&primary.model.eq('elasticnet')&primary.center.ne('AllCenters')&primary.status.eq('complete')]
    available=[p for p in pipelines if p in set(g.pipeline)]
    for i,center in enumerate(protocol['primary_studies']):
        values=g[g.center.eq(center)].set_index('pipeline').reindex(available)
        bars.bar(np.arange(len(available))+(i-.5)*.35,values.auroc,width=.35,color=colors[i],label='Test '+short_center(center))
    bars.set(xticks=range(len(available)),xticklabels=[labels[p] for p in available],ylim=(0,1),ylabel='Held-out study AUROC')
    bars.tick_params(axis='x',rotation=17,labelsize=6);bars.legend(fontsize=7)
    b=outer[0,1].subgridspec(1,2);roc=fig.add_subplot(b[0]);pr=fig.add_subplot(b[1])
    for pipeline in available:
        g=pred[pred.cohort_scope.eq(PRIMARY)&pred['mode'].eq('loco')&pred.model.eq('elasticnet')&pred.pipeline.eq(pipeline)]
        row=primary[primary['mode'].eq('loco')&primary.model.eq('elasticnet')&primary.pipeline.eq(pipeline)&primary.center.eq('AllCenters')].iloc[0]
        fpr,tpr,_=roc_curve(g.label,g.probability);prec,recall,_=precision_recall_curve(g.label,g.probability)
        roc.plot(fpr,tpr,color=color[pipeline],label=f'{labels[pipeline]} {row.auroc:.3f}')
        pr.plot(recall,prec,color=color[pipeline],label=f'{labels[pipeline]} {row.auprc:.3f}')
    roc.plot([0,1],[0,1],'--',color='gray',lw=.6);roc.set(title='b  Primary pooled LOCO ROC',xlabel='FPR',ylabel='TPR');roc.legend(fontsize=5)
    pr.set(title='Primary precision–recall',xlabel='Recall',ylabel='Precision');pr.legend(fontsize=5)
    test=significance[significance.cohort_scope.eq(PRIMARY)&significance.model.eq('elasticnet')&significance.method.eq('Hybrid_MultiScale')&significance.comparator.eq('Global_Zscore')]
    annotation='Hybrid vs Global, held-out DeLong p: '+', '.join(f'{short_center(r.center)} {r.p:.3g}' for r in test.itertuples())
    roc.text(0,-.22,annotation,transform=roc.transAxes,fontsize=6)
    cal=fig.add_subplot(outer[1,0]);cal.plot([0,1],[0,1],'--',color='gray')
    for pipeline in available:
        g=pred[pred.cohort_scope.eq(PRIMARY)&pred['mode'].eq('loco')&pred.model.eq('elasticnet')&pred.pipeline.eq(pipeline)]
        observed,expected=calibration_curve(g.label,g.probability,n_bins=10)
        row=primary[primary['mode'].eq('loco')&primary.model.eq('elasticnet')&primary.pipeline.eq(pipeline)&primary.center.eq('AllCenters')].iloc[0]
        cal.plot(expected,observed,'o-',color=color[pipeline],ms=3,label=f'{labels[pipeline]}: Brier {row.brier:.3f}')
    cal.set(xlim=(0,1),ylim=(0,1),title='c  Primary bilateral calibration',xlabel='Predicted cancer probability',ylabel='Observed cancer proportion');cal.legend(fontsize=6)
    drift=fig.add_subplot(outer[1,1]);summary=stress.groupby('pipeline').max_anchor_probability_range.max().reindex(protocol['pipelines'])
    drift.barh(range(len(summary)),summary,color=[color[p] for p in summary.index])
    drift.set(yticks=range(len(summary)),yticklabels=[labels[p].replace(' (ICC bypassed*)','') for p in summary.index],
        xlabel='Maximum fixed-anchor probability change',title='d  Target-composition stress using real held-out people')
    xmax=max(.01,float(summary.max())*1.2);drift.set_xlim(0,xmax)
    for i,value in enumerate(summary):drift.text(value+xmax*.01,i,f'{value:.3g}',va='center',fontsize=8)
    drift.text(.02,.02,'Identical anchors; target cancer fractions 0.25, 0.50, 0.75\nPost-evaluation diagnostic; never used for model selection',transform=drift.transAxes,fontsize=7)
    save(fig,folder,'Figure3_LSN_Debiasing_Performance',300)

    # Complete within/cross-study rank contribution, not a causal destruction claim.
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    table=pairs.groupby(['pipeline','pair_type']).pooled_auc_contribution.sum().unstack().reindex(pipelines)
    table.plot.bar(stacked=True,ax=axes[0],color=colors[:2]);axes[0].set(ylabel='Contribution to pooled three-study LOCO AUROC',title='Within- versus cross-study cancer–control ranking')
    axes[0].tick_params(axis='x',rotation=30,labelsize=7);axes[0].legend(fontsize=6)
    g=metrics[metrics.cohort_scope.eq(COMBINED)&metrics['mode'].eq('loco')&metrics.model.eq('elasticnet')&metrics.center.eq('AllCenters')].set_index('pipeline').reindex(pipelines)
    axes[1].plot(g.auroc,np.arange(len(g)),'o',label='Pooled three-study AUROC')
    axes[1].plot(g.macro_bilateral_auroc,np.arange(len(g)),'s',label='Mean of two primary held-out-study AUROCs')
    axes[1].set(yticks=range(len(g)),yticklabels=[labels[p].replace(' (ICC bypassed*)','').replace('*','') for p in g.index],xlim=(0,1),title='Discrimination and cross-fold score comparability');axes[1].legend(fontsize=6)
    save(fig,folder,'FigureS5_ComBat_Score_Comparability',300)
