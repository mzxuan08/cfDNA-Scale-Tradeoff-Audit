"""Test observed effect retention instead of inferring signal loss from pooled AUC."""
import numpy as np
import pandas as pd
from .study import load_clinical,export_table
from .multiscale_benchmark import read_protocol,fit_preprocessor
from .debiasing import ArmLocalLogRatio


def cohen_d_by_feature(x,y):
    x=np.asarray(x,dtype='float64');y=np.asarray(y)
    a=x[y==1];b=x[y==0]
    sd=np.sqrt(((len(a)-1)*a.var(axis=0,ddof=1)+(len(b)-1)*b.var(axis=0,ddof=1))/(len(a)+len(b)-2))
    return (a.mean(axis=0)-b.mean(axis=0))/np.maximum(sd,1e-12)


def heldout_effect_diagnostic(cfg):
    meta,R,_,bins=load_clinical(cfg,True);rows=[]
    for center in read_protocol(cfg)['primary_studies']:
        train=meta.center.ne(center).to_numpy();test=~train;y=meta.loc[test,'label'].to_numpy()
        raw=cohen_d_by_feature(R[test],y)
        for pipeline in ['Raw_Ratio','Global_Zscore','InFold_ComBat','Proposed_LSN_ICC','Hybrid_MultiScale']:
            pre,_=fit_preprocessor(pipeline,bins,R[train],meta.loc[train,'label'],meta.loc[train,'center'])
            transformed=pre.transform(R[test],meta.loc[test,'center'])
            effects=cohen_d_by_feature(transformed,y)
            if pipeline=='Hybrid_MultiScale':
                blocks=[('global',effects[:len(bins)]),('local',effects[len(bins):])]
            else:blocks=[('all',effects)]
            for block,d in blocks:
                rows.append(dict(target_study=center,pipeline=pipeline,block=block,windows=len(d),
                    median_abs_raw_d=float(np.median(np.abs(raw))),median_abs_transformed_d=float(np.median(np.abs(d))),
                    max_abs_d_change=float(np.max(np.abs(d-raw))) if len(d)==len(raw) else np.nan,
                    interpretation='positive_per_feature_affine_maps_preserve_within_study_standardized_mean_differences' if pipeline=='InFold_ComBat' else 'samplewise_nonlinear_or_filtered_representation'))
    result=pd.DataFrame(rows);export_table(result,cfg,'Step7_Heldout_Effect_Preservation');return result


def multiscale_signal_diagnostics(cfg):
    from scipy.stats import spearmanr
    from threadpoolctl import threadpool_limits
    from .variance_audit import marginal_variance
    from .evaluation import bh_fdr
    from pathlib import Path
    meta,R,_,bins=load_clinical(cfg,True)
    annotation=pd.read_parquet(Path(cfg['paths']['reference'])/'hg19_window_annotations.parquet').set_index('bin_id').loc[bins.bin_id]
    representations={'Raw_Ratio':R,'Global_Zscore':((R-R.mean(axis=1,keepdims=True))/R.std(axis=1,keepdims=True)).astype('float32'),
                     'LSN':ArmLocalLogRatio(bins,3).fit(R).transform(R)}
    correlations=[];partitions=[];effects=[]
    with threadpool_limits(limits=1):
        for name,x in representations.items():
            variance=marginal_variance(x,meta)
            rho,p=spearmanr(annotation.gc_fraction,variance.eta2_batch)
            correlations.append(dict(representation=name,rho=float(rho),p=float(p),windows=len(bins),
                median_eta2_batch=float(variance.eta2_batch.median()),median_eta2_bio=float(variance.eta2_bio.median()),
                scope='descriptive_coverage_and_diagnosis_adjusted_unique_variance_not_CV_preprocessing'))
            partitions.append(pd.concat([bins[['bin_id']].reset_index(drop=True),variance],axis=1).assign(representation=name))
            for center in read_protocol(cfg)['primary_studies']:
                subset=meta.center.eq(center).to_numpy();y=meta.loc[subset,'label'].to_numpy()
                for arm,group in bins.groupby(['chrom','arm'],sort=False):
                    values=np.mean(x[subset][:,group.index],axis=1,keepdims=True)
                    d=cohen_d_by_feature(values,y)[0]
                    effects.append(dict(representation=name,center=center,arm=''.join(arm),windows=len(group),cohen_d=float(d),
                        interpretation='arm_averaged_fragmentation_contrast_not_measured_copy_number'))
    frame=pd.DataFrame(correlations);frame['q']=bh_fdr(frame.p)
    export_table(frame,cfg,'Step7_Residual_Batch_GC_Associations')
    export_table(pd.concat(partitions,ignore_index=True),cfg,'Step7_Representation_Variance')
    export_table(pd.DataFrame(effects),cfg,'Step7_Arm_Level_Disease_Contrasts')
    return frame


def qualify_paired_inference(cfg,frame):
    """Ordinary fixed-rule DeLong does not cover a target-batch fitted mapping."""
    from .evaluation import bh_fdr
    result=frame.copy()
    adapted=result['method'].eq('InFold_ComBat')|result.comparator.eq('InFold_ComBat')
    computed=adapted&result['p'].notna()
    result.loc[computed,'status']='DeLong_not_applicable_transductive_target_batch_dependency'
    result.loc[adapted,'p']=np.nan
    result['q']=np.nan
    for _,idx in result.groupby('cohort_scope').groups.items():result.loc[idx,'q']=bh_fdr(result.loc[idx,'p'])
    export_table(result,cfg,'Step7_Paired_Significance')
    return result


def source_threshold_alignment_diagnostic(cfg,pred):
    from scipy.special import logit
    from sklearn.metrics import roc_auc_score
    rows=[]
    subset=pred[pred.model.eq('elasticnet')&pred['mode'].eq('loco')]
    for (scope,pipeline),group in subset.groupby(['cohort_scope','pipeline']):
        p=np.clip(group.probability.to_numpy(dtype='float64'),np.finfo(float).tiny,np.nextafter(1.,0.))
        t=np.clip(group.threshold.to_numpy(dtype='float64'),np.finfo(float).tiny,np.nextafter(1.,0.))
        margin=logit(p)-logit(t)
        raw_auc=roc_auc_score(group.label,p);aligned_auc=roc_auc_score(group.label,margin)
        rows.append(dict(cohort_scope=scope,pipeline=pipeline,n=len(group),original_pooled_auroc=raw_auc,
            source_threshold_margin_auroc=aligned_auc,delta=aligned_auc-raw_auc,
            scope='post_evaluation_score_comparability_diagnostic_source_thresholds_only_not_a_replacement_primary_result'))
    result=pd.DataFrame(rows);export_table(result,cfg,'Step7_Source_Threshold_Alignment');return result
