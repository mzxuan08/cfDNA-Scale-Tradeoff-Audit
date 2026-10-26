"""Compact revision figures with a physical minimum text size of 8.5 points."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
from .major_revision import PRIMARY,COMBINED,EXTERNAL,PRIMARY_STUDIES,METHODS,display_method

COLORS=['#3C5488','#00A087','#E64B35','#8491B4','#F39B7F']
SHORT={'Raw_Ratio':'Raw','Global_Zscore':'Global','Local_LSN':'Local','Local_LSN_ICC':'Local + ICC','Hybrid_MultiScale':'Hybrid','EB-TargetMap':'EB-TargetMap'}

def revision_figures(analysis, figure1_only=False):
    a=analysis;root=a.root;folder=a.out/'figures';folder.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':9.5,'axes.labelsize':9,
        'xtick.labelsize':8.8,'ytick.labelsize':8.8,'legend.fontsize':8.8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    table=lambda n:pd.read_csv(a.tables/(n+'.csv'))
    validations=[]
    def save(fig,name):
        for text in fig.findobj(match=Text):
            if text.get_text():text.set_fontsize(max(8.8,text.get_fontsize()))
        fig.savefig(folder/(name+'.pdf'),bbox_inches='tight',pad_inches=.06)
        fig.savefig(folder/(name+'.png'),dpi=300,bbox_inches='tight',pad_inches=.06)
        validations.append({'figure':name,'minimum_text_pt':min(t.get_fontsize() for t in fig.findobj(match=Text) if t.get_text()),'width_inches':fig.get_figwidth(),'height_inches':fig.get_figheight(),'png_dpi':300})
        plt.close(fig)
    def short(s):return s.split(' et al')[0]
    pred=a.pred[a.pred.cohort_scope.eq(PRIMARY)&a.pred['mode'].eq('loco')]
    methods=['Raw_Ratio','Global_Zscore','Proposed_LSN_ICC','Hybrid_MultiScale']
    pmetrics=[]
    for (study,method),g in pred[pred.pipeline.isin(methods)].groupby(['center','pipeline']):
        from sklearn.metrics import roc_auc_score,average_precision_score
        h=g.label.eq(0)
        pmetrics.append(dict(study=study,method=display_method(method,PRIMARY),auroc=roc_auc_score(g.label,g.probability),auprc=average_precision_score(g.label,g.probability),specificity=(g.loc[h,'probability']<=g.loc[h,'threshold']).mean()))
    metrics=pd.DataFrame(pmetrics)
    # Figure 1 keeps descriptive study geometry, then pairs it with dispersion.
    fig=plt.figure(figsize=(6.5,7.5),layout='constrained');gs=fig.add_gridspec(2,2)
    sub=gs[0,0].subgridspec(2,1,height_ratios=[1,2]);ax=fig.add_subplot(sub[0]);curve=fig.add_subplot(sub[1])
    counts=pd.crosstab(a.meta.center,a.meta.label)
    ax.barh(range(3),counts[0],color=COLORS[0],label='Healthy');ax.barh(range(3),counts[1],left=counts[0],color=COLORS[2],label='Cancer')
    ax.set(yticks=range(3),yticklabels=[short(c) for c in counts.index],title='a  Public study composition',xlabel='Participants');ax.legend(ncol=2,loc='upper right')
    for i,c in enumerate(counts.index):curve.plot(range(80,401),a.lengths[a.meta.center.eq(c)&a.meta.label.eq(0)].mean(axis=0),color=COLORS[i],label=short(c))
    curve.set(xlim=(100,220),xlabel='Fragment length (bp)',ylabel='Healthy frequency');curve.legend()
    emb=pd.read_parquet(root/'results/tables/clinical_embeddings.parquet');emb=emb[emb.family.eq('raw_ratio')]
    sub=gs[0,1].subgridspec(2,1)
    for j,col in enumerate(['center','label']):
        ax=fig.add_subplot(sub[j]);values=sorted(emb[col].unique())
        for i,v in enumerate(values):
            g=emb[emb[col].eq(v)];ax.scatter(g.UMAP1,g.UMAP2,s=5,alpha=.7,color=COLORS[i],label=short(v) if isinstance(v,str) else ('Cancer' if v else 'Healthy'))
        perma=pd.read_csv(root/'results/tables/Audit_PERMANOVA.csv');test=perma[perma.family.eq('raw_ratio')&perma.term.eq('center' if j==0 else 'disease')].iloc[0]
        title=('Study' if j==0 else 'Disease')+f': unique Type III R²={test.R2:.3f}; p={test.p:.3f}'
        ax.set(xlabel='UMAP 1',ylabel='UMAP 2',title=title.replace('; p=', ';\np='));ax.legend(loc='best',framealpha=.8)
    ax=fig.add_subplot(gs[1,0]);v=pd.read_csv(root/'results/tables/Bin_Variance_Partition.csv');v=v[v.family.eq('raw_ratio')]
    ax.plot(np.arange(len(v)),v.eta2_batch*100,'.',ms=2,color=COLORS[2],label='Study');ax.plot(np.arange(len(v)),v.eta2_bio*100,'.',ms=2,color=COLORS[0],label='Disease')
    ticks=[];labels=[]
    for chrom,g in v.reset_index(drop=True).groupby('chrom',sort=False):ticks.append(np.mean(g.index));labels.append(chrom[3:])
    ax.set(title='c  Unique window variance',xlabel='Autosomal genomic order',ylabel='Type III eta² (%)',xticks=ticks[::4],xticklabels=labels[::4]);ax.legend()
    ax=fig.add_subplot(gs[1,1]);d=table('PERMDISP_Patient_Distances');d=d[d.family.eq('raw_ratio')&d.term.eq('study')]
    ax.boxplot([g.distance for _,g in d.groupby('group')],tick_labels=[short(c) for c in sorted(d.group.unique())],showfliers=False)
    ax.tick_params(axis='x',rotation=30);ax.set(title='d  Within-study dispersion',ylabel='Distance to centroid')
    save(fig,'Figure1_Batch_Vulnerability_Audit')
    if figure1_only:
        from .common import write_json
        previous=a.out/'reports/figure_typography.json'
        old=pd.read_json(previous).to_dict('records') if previous.exists() else []
        write_json(validations+[r for r in old if r['figure']!='Figure1_Batch_Vulnerability_Audit'],previous)
        return validations
    # Primary per-target and equal-weight summaries are the main performance views.
    fig,axes=plt.subplots(2,2,figsize=(6.5,6.7),layout='constrained')
    display=[display_method(m,PRIMARY) for m in methods]
    for i,c in enumerate(PRIMARY_STUDIES):
        g=metrics[metrics.study.eq(c)].set_index('method').loc[display]
        axes[0,0].plot(g.auroc,np.arange(4)+(i-.5)*.15,'o-',color=COLORS[i],label='Target '+short(c))
    axes[0,0].set(yticks=range(4),yticklabels=[SHORT[m] for m in display],xlim=(.5,1),xlabel='Held-out AUROC',title='a  Direction-specific ranking');axes[0,0].legend()
    mat=metrics.pivot(index='method',columns='study',values='auroc').loc[display,PRIMARY_STUDIES]
    ax=axes[0,1];im=ax.imshow(mat,cmap='Blues',vmin=.5,vmax=1,aspect='auto')
    ax.set(yticks=range(4),yticklabels=[SHORT[m] for m in display],xticks=[0,1],xticklabels=['Cristiano','Jiang'],title='b  Primary AUROC')
    for i in range(4):
        for j in range(2):ax.text(j,i,f'{mat.iloc[i,j]:.3f}',ha='center',va='center',color='white' if mat.iloc[i,j]>.82 else 'black')
    ax.set_xticks(np.arange(-.5,2,1),minor=True);ax.set_yticks(np.arange(-.5,4,1),minor=True);ax.grid(which='minor',color='white',lw=1);ax.tick_params(which='minor',length=0)
    thresholds=table('Threshold_Calibration_Comparison');g=thresholds[thresholds.scope.eq(PRIMARY)&thresholds.scheme.eq('PooledSource3Fold')&thresholds.status.eq('complete')]
    for i,m in enumerate(display):
        h=g[g.method.eq(m)].set_index('target').loc[PRIMARY_STUDIES]
        axes[1,0].plot(np.arange(2)+(i-1.5)*.06,h.target_specificity_shift_pp,'o-',color=COLORS[i],label=SHORT[m])
    axes[1,0].axhline(0,ls='--',color='gray');axes[1,0].set(xticks=[0,1],xticklabels=['Cristiano','Jiang'],ylabel='Specificity shift (pp)',title='c  Locked source threshold');axes[1,0].legend(ncol=2)
    means=metrics.groupby('method').auroc.mean().reindex(display)
    axes[1,1].barh(range(4),means,color=COLORS[:4]);axes[1,1].set(yticks=range(4),yticklabels=[SHORT[m] for m in display],xlim=(0,1),xlabel='Equal-weight mean AUROC',title='d  Bilateral mean')
    for i,value in enumerate(means):axes[1,1].text(value+.02,i,f'{value:.3f}',va='center')
    save(fig,'Figure2_Primary_Ranking_Threshold_Transport')
    fig,axes=plt.subplots(2,2,figsize=(6.5,7),layout='constrained')
    axes[0,0].axis('off');axes[0,0].text(.02,.95,'a  Patient-wise invariance',va='top',weight='bold')
    axes[0,0].text(.03,.75,'Frozen source parameters θ\n\nzᵢ = fθ(xᵢ)\npᵢ = gθ(zᵢ)\n\nChanging companion set B\nleaves pᵢ unchanged.\n\nMathematical property\nverified by unit tests.',va='top',linespacing=1.4)
    d=table('Anchor_All_Predictions');eb=d[d.method.eq('EB-TargetMap')]
    ax=axes[0,1]
    groups=[];positions=[];tick=[]
    for j,c in enumerate(PRIMARY_STUDIES):
        for i,pi in enumerate([.25,.5,.75]):groups.append(eb[eb.target.eq(c)&eb.target_case_fraction.eq(pi)].probability_shift);positions.append(j*4+i);tick.append(f'{pi:.2f}')
    ax.boxplot(groups,positions=positions,tick_labels=tick,showfliers=False);ax.axhline(0,color='gray',ls='--')
    ax.set(title='b  EB-TargetMap anchor shifts',xlabel='Case fraction: Cristiano | Jiang',ylabel='Probability minus reference')
    ax=axes[1,0];fixed=eb[eb.target_case_fraction.eq(.5)]
    for i,(c,h) in enumerate(fixed.groupby('target')):ax.hist(h.probability_shift,bins=20,alpha=.5,color=COLORS[i],label=short(c),density=True)
    ax.set(title='c  Fixed case fraction 0.50',xlabel='Companion-identity shift',ylabel='Density');ax.legend()
    ax=axes[1,1];summ=table('Anchor_Distribution_Summary');summ=summ[summ.method.eq('EB-TargetMap')]
    for i,(c,h) in enumerate(summ.groupby('target')):ax.plot(h.case_fraction,h.threshold_crossing_flip_rate,'o-',color=COLORS[i],label=short(c))
    ax.set(title='d  Decision flips',xlabel='Target case fraction',ylabel='Threshold-crossing rate',ylim=(-.01,1.01));ax.legend()
    save(fig,'Figure3_Companion_Dependence_and_Invariance')
    fig,axes=plt.subplots(2,2,figsize=(6.5,7.6),layout='constrained')
    anno=pd.read_parquet(Path(a.cfg['paths']['reference'])/'hg19_window_annotations.parquet').set_index('bin_id').loc[a.bins.bin_id]
    variance=v.set_index('bin_id').loc[a.bins.bin_id].eta2_batch.to_numpy();block=table('Chromosome_Block_Associations')
    for ax,col,label in [(axes[0,0],'gc_fraction','GC fraction'),(axes[0,1],'distance_to_telomere_bp','Telomere distance (Mb)')]:
        x=anno[col].to_numpy()/(1e6 if col!='gc_fraction' else 1);ax.scatter(x,variance*100,s=7,alpha=.5,color=COLORS[0])
        row=block[block.representation.eq('Raw_Ratio')&block.annotation.eq(col)].iloc[0]
        ax.set(xlabel=label,ylabel='Unique study eta² (%)',title=f'{"a" if col=="gc_fraction" else "b"}  ρ={row.rho:.3f}\nBlock CI [{row.block_ci_low:.3f}, {row.block_ci_high:.3f}]')
    from .debiasing import ArmLocalLogRatio
    local=ArmLocalLogRatio(a.bins,3).fit(a.R).transform(a.R);mask=a.bins.chrom.eq('chr8').to_numpy();pos=a.bins.loc[mask,'midpoint_bp'].to_numpy()/1e6
    spec=axes[1,0].get_subplotspec();axes[1,0].remove();sub=spec.subgridspec(2,1);rawax=fig.add_subplot(sub[0]);ax=fig.add_subplot(sub[1])
    breaks=np.flatnonzero(np.diff(pos)>5)+1;positions=np.insert(pos,breaks,np.nan)
    for i,c in enumerate(PRIMARY_STUDIES):
        h=a.meta.center.eq(c)&a.meta.label.eq(0)
        rawax.plot(positions,np.insert(np.median(a.R[h][:,mask],axis=0),breaks,np.nan),color=COLORS[i],label=short(c))
        ax.plot(positions,np.insert(np.median(local[h][:,mask],axis=0),breaks,np.nan),color=COLORS[i],label=short(c))
    rawax.set(title='c  Healthy chr8 trajectories',ylabel='Raw ratio');rawax.legend()
    ax.set(xlabel='chr8 coordinate (Mb)',ylabel='Local log₂ ratio')
    imp=pd.read_parquet(root/'results/tables/Clinical_Feature_Importance.parquet')
    imp=imp[imp.pipeline.eq('Proposed_LSN_ICC')&imp.model.eq('elasticnet')&imp['mode'].eq('loco')]
    top=imp.groupby('bin_id').importance.sum().nlargest(20).index;rank=imp.pivot(index='bin_id',columns='fold',values='rank').reindex(top)
    # Compact bin labels preserve all coordinates in the accompanying rank table.
    compact=[s.split(':')[0]+':'+str(int(s.split(':')[1].split('-')[0])//1000000)+'–'+str(int(s.split('-')[1])//1000000) for s in rank.index]
    ax=axes[1,1];cmap=plt.get_cmap('Blues_r').copy();cmap.set_bad('#EEEEEE');im=ax.imshow(rank,cmap=cmap,aspect='auto',vmin=1,vmax=477)
    ax.set(yticks=range(20),yticklabels=compact,xticks=range(3),xticklabels=['C','J','S'],title='d  Feature ranks (dark = high)',ylabel='Window coordinates (Mb)')
    ax.set_xticks(np.arange(-.5,3,1),minor=True);ax.set_yticks(np.arange(-.5,20,1),minor=True);ax.grid(which='minor',color='white',lw=.5);ax.tick_params(which='minor',length=0)
    fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.045,pad=.04,ticks=[1,120,240,360,477],label='Feature rank (1 = highest)')
    save(fig,'Figure4_Genomic_Context_and_Ranks');a.save(rank.reset_index(),'Feature_Rank_Heatmap')
    # Secondary pooled scores explicitly serve score-comparability diagnostics.
    fig,axes=plt.subplots(2,1,figsize=(6.5,6.5),layout='constrained')
    pairs=pd.read_csv(root/'results/tables/Step7_AUC_Pair_Decomposition.csv');pairs=pairs[pairs.pipeline.isin(METHODS)]
    pivot=pairs.groupby(['pipeline','pair_type']).pooled_auc_contribution.sum().unstack().reindex(METHODS)
    pivot.index=[SHORT.get(display_method(x),display_method(x)) for x in pivot.index]
    pivot.plot.bar(stacked=True,ax=axes[0],color=COLORS[:2],rot=0);axes[0].set(title='a  Secondary pooled score comparability',ylabel='Pooled AUROC contribution');axes[0].legend(['Cross-study pairs','Within-study pairs'])
    threshold=table('Threshold_Calibration_Comparison');g=threshold[threshold.scope.eq(EXTERNAL)&threshold.status.eq('complete')&threshold.method.ne('EB-TargetMap')]
    for i,(scheme,h) in enumerate(g.groupby('scheme')):
        axes[1].plot(h.threshold,np.arange(len(h))+(i-.5)*.1,'o',color=COLORS[i],label=scheme.replace('LeaveOneSourceStudyOut','Source-study held out').replace('PooledSource3Fold','Mixed source OOF'))
        axes[1].set(yticks=range(len(h)),yticklabels=[SHORT.get(m,m) for m in h.method])
    axes[1].set(title='b  Source calibration before Snyder testing',xlabel='Source-only 95% specificity threshold');axes[1].legend()
    save(fig,'FigureS5_Score_Comparability_and_Calibration')
    # Secondary classifier comparisons belong in the supplement.
    allmetrics=pd.read_csv(root/'results/tables/Table3_Benchmark_Metrics_Full.csv')
    g=allmetrics[allmetrics.cohort_scope.eq(PRIMARY)&allmetrics['mode'].eq('loco')&allmetrics.center.eq('AllCenters')]
    mat=g[g.pipeline.isin(METHODS)].pivot(index='pipeline',columns='model',values='macro_bilateral_auroc').reindex(METHODS)
    fig,ax=plt.subplots(figsize=(6.5,3.3),layout='constrained');cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#EEEEEE');im=ax.imshow(mat,aspect='auto',vmin=.5,vmax=1,cmap=cmap)
    ax.set(yticks=range(len(mat)),yticklabels=[SHORT.get(display_method(m,PRIMARY),display_method(m,PRIMARY)) for m in mat.index],xticks=range(len(mat.columns)),xticklabels=mat.columns,title='Secondary classifier bilateral mean AUROC')
    for i in range(len(mat)):
        for j in range(len(mat.columns)):
            value=mat.iloc[i,j];ax.text(j,i,f'{value:.3f}' if np.isfinite(value) else 'NA',ha='center',va='center',color='white' if np.isfinite(value) and value>.83 else 'black')
    ax.set_xticks(np.arange(-.5,len(mat.columns),1),minor=True);ax.set_yticks(np.arange(-.5,len(mat),1),minor=True);ax.grid(which='minor',color='white',lw=1);ax.tick_params(which='minor',length=0)
    save(fig,'FigureS6_Secondary_Classifiers')
    from .common import write_json
    write_json(validations,a.out/'reports/figure_typography.json')
