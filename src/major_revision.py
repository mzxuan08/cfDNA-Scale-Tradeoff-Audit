"""Real-patient methodological extensions, isolated from submitted v1 outputs."""
from __future__ import annotations
import copy,hashlib,json,pickle,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr,f_oneway
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
from threadpoolctl import threadpool_limits
from .common import sha256,write_json
from .study import load_clinical
from .models import make_model
from .multiscale_benchmark import fit_preprocessor,PRIMARY,COMBINED,EXTERNAL,exact_interval
from .evaluation import (source_calibration_splits,threshold_bootstrap,
                         paired_auc_difference_ci,bh_fdr,ece)

METHODS=['Raw_Ratio','Global_Zscore','Proposed_LSN_ICC','Hybrid_MultiScale','InFold_ComBat']
PRIMARY_STUDIES=['Cristiano et al, 2019','Jiang et al, 2015']

def display_method(method,scope=COMBINED):
    if method=='InFold_ComBat':return 'EB-TargetMap'
    if method=='Proposed_LSN_ICC':return 'Local_LSN' if scope==PRIMARY else 'Local_LSN_ICC'
    if method=='Tuned_LSN_ICC':return 'Tuned_Local_LSN' if scope==PRIMARY else 'Tuned_Local_LSN_ICC'
    return method

def block_rank_inference(x,y,blocks,repeats=1000,permutations=999,seed=2026):
    """Chromosome cluster bootstrap and block-wild sign randomization.

    Whole chromosomes are resampled for the CI. For the null test, sum centered
    rank cross-products within each chromosome, then independently flip their
    signs. It assumes sign-symmetric independent chromosome contributions and
    retains arbitrary within-chromosome dependence; it is not an IID-bin test.
    """
    x=np.asarray(x,float);y=np.asarray(y,float);blocks=np.asarray(blocks)
    ok=np.isfinite(x)&np.isfinite(y);x,y,blocks=x[ok],y[ok],blocks[ok]
    ids=[np.flatnonzero(blocks==b) for b in np.unique(blocks)]
    rng=np.random.default_rng(seed);rho=float(spearmanr(x,y).statistic);values=[]
    for _ in range(repeats):
        draw=np.concatenate([ids[j] for j in rng.integers(len(ids),size=len(ids))])
        values.append(spearmanr(x[draw],y[draw]).statistic)
    rx=rankdata(x)-np.mean(rankdata(x));ry=rankdata(y)-np.mean(rankdata(y))
    contributions=np.array([np.sum(rx[g]*ry[g]) for g in ids]);observed=abs(contributions.sum())
    signs=rng.choice([-1,1],size=(permutations,len(ids)))
    p=(1+np.sum(np.abs(signs@contributions)>=observed-1e-12))/(permutations+1)
    return dict(rho=rho,block_ci_low=float(np.quantile(values,.025)),block_ci_high=float(np.quantile(values,.975)),
                block_permutation_p=float(p),blocks=len(ids),windows=len(x),bootstrap_replicates=repeats,permutations=permutations,
                null='chromosome_block_wild_sign_of_centered_rank_crossproducts',assumption='independent_sign_symmetric_chromosome_contributions')

def permdisp_euclidean(x,groups,permutations=999,seed=2026):
    """Bias-adjusted centroid distances, with intercept-only residual permutation.

    Distances are fixed relative to the original group centroids. Permuting the
    intercept-only distance residuals is equivalent to permuting these distances
    across group labels. It avoids mixing locations by relabelling raw vectors.
    This is an unrestricted, unadjusted dispersion diagnostic, distinct from the
    coverage-adjusted restricted Freedman-Lane location/association analysis.
    """
    x=np.asarray(x,float);groups=np.asarray(groups);labels=np.unique(groups)
    if min(np.sum(groups==c) for c in labels)<2:raise ValueError('PERMDISP requires two members per group')
    norms=np.sum(x*x,axis=1);D=np.maximum(norms[:,None]+norms[None,:]-2*x@x.T,0)
    def distances(g):
        d=np.empty(len(x))
        for c in labels:
            ids=np.flatnonzero(g==c);n=len(ids)
            d[ids]=np.sqrt(np.maximum(D[np.ix_(ids,ids)].mean(axis=1)-D[np.ix_(ids,ids)].mean()/2,0))*np.sqrt(n/(n-1))
        return d
    d=distances(groups);stat=float(f_oneway(*[d[groups==c] for c in labels]).statistic)
    rng=np.random.default_rng(seed);exceed=0
    for _ in range(permutations):
        v=rng.permutation(d-d.mean())+d.mean()
        exceed+=f_oneway(*[v[groups==c] for c in labels]).statistic>=stat-1e-12
    return dict(F=stat,p=(exceed+1)/(permutations+1),permutations=permutations,center='centroid',
                bias_adjusted=True,permutation='intercept_only_residuals_of_fixed_centroid_distances'),d

class RevisionAnalysis:
    def __init__(self,cfg):
        self.cfg=cfg;self.root=Path(cfg['_root']);self.out=self.root/'results/step09_major'
        self.tables=self.out/'tables';self.cache=self.out/'checkpoints'
        for p in [self.tables,self.cache,self.out/'figures',self.out/'reports']:p.mkdir(parents=True,exist_ok=True)
        # load_clinical may export identity exclusions: redirect those bookkeeping
        # writes into the new revision directory to preserve frozen v1 tables.
        local=copy.deepcopy(cfg);local['paths']['results']=str(self.out)
        self.meta,self.R,self.lengths,self.bins=load_clinical(local,True)
        assert len(self.meta)==585 and self.meta.patient_id.nunique()==585
        self.pred=pd.read_parquet(self.root/'results/tables/Step7_OOF_Predictions.parquet')
        self.pred=self.pred[self.pred.model.eq('elasticnet')].copy()
        self.fit_diagnostics=[];self.calibration_records=[];self.split_records=[]
        self.data_signature=sha256(Path(cfg['paths']['processed'])/cfg['extraction']['final_file'])

    def save(self,frame,name):
        frame.to_csv(self.tables/(name+'.csv'),index=False)
        frame.to_parquet(self.tables/(name+'.parquet'),index=False)
        return frame

    def fit(self,train,test,method):
        train=np.asarray(train);test=np.asarray(test);a=self.meta.iloc[train];b=self.meta.iloc[test]
        assert not set(a.patient_id)&set(b.patient_id)
        pre,info=fit_preprocessor(method,self.bins,self.R[train],a.label,a.center,3,.25)
        model=make_model('elasticnet',self.cfg)
        model.fit(pre.transform(self.R[train],a.center),a.label)
        probability=model.predict_proba(pre.transform(self.R[test],b.center))[:,1].astype('float32')
        self.fit_diagnostics.append(dict(method=display_method(method,PRIMARY if a.center.nunique()==1 else COMBINED),
            train_n=len(train),test_n=len(test),train_studies=';'.join(sorted(a.center.unique())),
            test_studies=';'.join(sorted(b.center.unique())),iterations=int(model.n_iter_.max()),
            converged=bool(model.n_iter_.max()<model.max_iter),**info))
        return probability,pre,model

    def calibration(self,train,method,scheme,outer):
        a=self.meta.iloc[train];scores=np.full(len(train),np.nan,dtype='float32');notes=[]
        if scheme=='PooledSource3Fold':splits,rare=source_calibration_splits(a,3,2026)
        else:
            if a.center.nunique()<2:return None,'single_source_LOSO_unavailable'
            splits=[(np.flatnonzero(a.center.ne(c)),np.flatnonzero(a.center.eq(c))) for c in sorted(a.center.unique())]
            rare=[]
        for k,(it,iv) in enumerate(splits):
            ti,vi=train[it],train[iv]
            if method=='InFold_ComBat' and a.iloc[it].center.nunique()<2:
                return None,'EB-TargetMap_source_harmonization_unidentifiable_in_one_source_LOSO'
            scores[iv],_,_=self.fit(ti,vi,method)
            for role,ids in [('train',ti),('validation',vi)]:
                self.split_records.extend(dict(outer=outer,scheme=scheme,method=display_method(method,PRIMARY if a.center.nunique()==1 else COMBINED),fold=k,role=role,patient_id=self.meta.iloc[j].patient_id,study=self.meta.iloc[j].center) for j in ids)
        controls=np.isfinite(scores)&a.label.eq(0).to_numpy()
        for j in np.flatnonzero(np.isfinite(scores)):
            self.calibration_records.append(dict(outer=outer,scheme=scheme,method=display_method(method,PRIMARY if a.center.nunique()==1 else COMBINED),patient_id=a.iloc[j].patient_id,study=a.iloc[j].center,label=int(a.iloc[j].label),probability=float(scores[j])))
        if not controls.any():return None,'no_heldout_source_controls'
        return scores[controls],'complete'

    def threshold_analysis(self):
        rows=[];drawrows=[]
        for scope in [PRIMARY,COMBINED,EXTERNAL]:
            chosen=self.pred[self.pred.cohort_scope.eq(scope)&self.pred['mode'].isin(['loco','external'])]
            for (target,method),g in chosen.groupby(['center','pipeline'],sort=True):
                if method not in METHODS:continue
                eligible=self.meta.center.isin(PRIMARY_STUDIES) if scope in [PRIMARY,EXTERNAL] else np.ones(len(self.meta),bool)
                train=np.flatnonzero(eligible&self.meta.center.ne(target));test=np.flatnonzero(self.meta.patient_id.isin(g.patient_id))
                assert not set(self.meta.iloc[train].center)&{target}
                outer=scope+'|'+target
                for scheme in ['PooledSource3Fold','LeaveOneSourceStudyOut']:
                    controls,status=self.calibration(train,method,scheme,outer)
                    row=dict(scope=scope,target=target,method=display_method(method,scope),scheme=scheme,source_n=len(train),target_n=len(g),status=status)
                    if controls is None:rows.append(row);continue
                    threshold=float(np.quantile(controls,.95,method='higher'))
                    draws,spec,sens=threshold_bootstrap(controls,g.probability,g.label,1000,2026)
                    healthy=g.label.eq(0);cancer=g.label.eq(1);n=int(healthy.sum());success=int((g.loc[healthy,'probability']<=threshold).sum())
                    lo,hi=exact_interval(success,n)
                    if scheme=='PooledSource3Fold':
                        row['frozen_threshold_difference']=threshold-float(g.threshold.iloc[0])
                        if abs(row['frozen_threshold_difference'])>2e-5:raise ValueError(('source calibration differs from frozen',row))
                    row.update(source_control_n=len(controls),threshold=threshold,threshold_low=float(np.quantile(draws,.025)),threshold_high=float(np.quantile(draws,.975)),
                        target_specificity=success/n if n else np.nan,target_specificity_low_exact=lo,target_specificity_high_exact=hi,
                        target_healthy_success=success,target_healthy_n=n,target_specificity_description=f'{success}/{n}'+(' (descriptive only)' if n==1 else ''),
                        target_specificity_shift_pp=(success/n-.95)*100 if n else np.nan,
                        specificity_from_threshold_low=float(np.quantile(spec,.025)),specificity_from_threshold_high=float(np.quantile(spec,.975)),
                        target_sensitivity=float((g.loc[cancer,'probability']>threshold).mean()),
                        sensitivity_from_threshold_low=float(np.quantile(sens,.025)),sensitivity_from_threshold_high=float(np.quantile(sens,.975)),
                        threshold_interval_scope='source_control_score_resampling_conditional_on_fixed_calibration_models',
                        local_calibration_branch='Local_LSN when one source; Local_LSN_ICC when multiple healthy sources')
                    rows.append(row)
                    drawrows.extend(dict(scope=scope,target=target,method=row['method'],scheme=scheme,replicate=j,threshold=float(v),target_specificity=float(spec[j])) for j,v in enumerate(draws))
        self.save(pd.DataFrame(self.calibration_records),'Source_Calibration_Predictions')
        self.save(pd.DataFrame(self.split_records),'Source_Calibration_Splits')
        self.save(pd.DataFrame(drawrows),'Threshold_Bootstrap_Distribution')
        return self.save(pd.DataFrame(rows),'Threshold_Calibration_Comparison')

    def pair_deltas(self):
        old=pd.read_csv(self.root/'results/tables/Step7_Paired_Significance.csv')
        family=old[old.cohort_scope.eq(PRIMARY)&old.p.notna()].copy();assert len(family)==38
        family['hypothesis_id']=['Primary_H'+str(i+1).zfill(2) for i in range(len(family))]
        family['family_m']=38;family['q_recomputed']=bh_fdr(family.p)
        np.testing.assert_allclose(family.q,family.q_recomputed)
        family['method']=[display_method(m,PRIMARY) for m in family.method]
        family['comparator']=[display_method(m,PRIMARY) for m in family.comparator]
        self.save(family,'Primary_BH_Hypothesis_Family')
        rows=[]
        for target,g in self.pred[self.pred.cohort_scope.eq(PRIMARY)&self.pred['mode'].eq('loco')].groupby('center'):
            wide=g.pivot(index=['patient_id','label'],columns='pipeline',values='probability').reset_index()
            match=family[family.center.eq(target)&family.model.eq('elasticnet')&family.method.eq('Hybrid_MultiScale')&family.comparator.eq('Global_Zscore')].iloc[0]
            rows.append(dict(target=target,method='Hybrid_MultiScale',comparator='Global_Zscore',p=match.p,q=match.q,family_m=38,
                **paired_auc_difference_ci(wide.label,wide.Hybrid_MultiScale,wide.Global_Zscore,1000,2026)))
        return self.save(pd.DataFrame(rows),'Primary_Hybrid_Global_Delta')

    def cohort_table(self):
        rows=[]
        for study,g in self.meta.groupby('center'):
            stages=g.loc[g.label.eq(1),'stage']
            types=g.loc[g.label.eq(1),'cancer_type'].value_counts()
            row=dict(study=study,n=len(g),cancer=int(g.label.sum()),healthy=int(g.label.eq(0).sum()),
                cancer_types='; '.join(f'{k}: {v}' for k,v in types.items()),
                stage_recorded_n=int(stages.notna().sum()),stage_missing_n=int(stages.isna().sum()),
                stage_distribution='; '.join(f'{k}: {v}' for k,v in stages.dropna().astype(str).value_counts().items()) or 'Not recorded')
            for name,col in [('sequencing_depth','seq_depth'),('accepted_fragments','accepted_fragments')]:
                for stat,number in zip(['min','q1','median','q3','max'],np.quantile(g[col],[0,.25,.5,.75,1])):row[name+'_'+stat]=float(number)
            rows.append(row)
        self.save(self.meta.groupby(['center','label','cancer_type'],dropna=False).size().reset_index(name='n'),'Cancer_Type_Composition')
        return self.save(pd.DataFrame(rows),'Table2_Study_Depth_Fragments_Composition')

    def sensitivities(self):
        rows=[];predictions=[]
        primary=self.meta.center.isin(PRIMARY_STUDIES)
        gi=['Bile Duct','Colorectal','Duodenal','Gastric','Pancreatic']
        scenarios={'Depth_2to4x':primary&self.meta.seq_depth.between(2.,4.),
                   'GI_vs_HCC':primary&(self.meta.label.eq(0)|self.meta.center.eq(PRIMARY_STUDIES[1])|self.meta.cancer_type.isin(gi))}
        for scenario,eligible in scenarios.items():
            for target in PRIMARY_STUDIES:
                train=np.flatnonzero(eligible&self.meta.center.ne(target));test=np.flatnonzero(eligible&self.meta.center.eq(target))
                for method in METHODS[:-1]:
                    controls,status=self.calibration(train,method,'PooledSource3Fold',scenario+'|'+target)
                    if controls is None:raise ValueError(status)
                    threshold=float(np.quantile(controls,.95,method='higher'));p,_,_=self.fit(train,test,method)
                    b=self.meta.iloc[test];y=b.label.to_numpy();h=y==0
                    point=dict(scenario=scenario,target=target,method=display_method(method,PRIMARY),train_n=len(train),target_n=len(test),target_cancer=int(y.sum()),target_healthy=int(h.sum()),
                        auroc=roc_auc_score(y,p),auprc=average_precision_score(y,p),threshold=threshold,
                        specificity=float((p[h]<=threshold).mean()),sensitivity=float((p[~h]>threshold).mean()),brier=brier_score_loss(y,p),ece=ece(y,p),
                        analysis_scope='retrained_source_only_sensitivity_no_target_label_selection',
                        cancer_subgroup='Cristiano GI incl Bile Duct/Colorectal/Duodenal/Gastric/Pancreatic vs Jiang HCC; no shared-HCC claim' if scenario=='GI_vs_HCC' else 'fixed_2_to_4x_original_depth_both_studies')
                    rows.append(point)
                    predictions.extend(dict(scenario=scenario,target=target,method=point['method'],patient_id=r.patient_id,label=int(r.label),probability=float(prob),threshold=threshold) for r,prob in zip(b.itertuples(),p))
        self.save(pd.DataFrame(predictions),'Sensitivity_Predictions')
        # Save the additional sensitivity calibration folds as well.
        self.save(pd.DataFrame(self.calibration_records),'Source_Calibration_Predictions')
        self.save(pd.DataFrame(self.split_records),'Source_Calibration_Splits')
        return self.save(pd.DataFrame(rows),'Depth_CancerType_Sensitivity')

    def anchors(self,repeats=50):
        records=[];members=[]
        for target in PRIMARY_STUDIES:
            train=np.flatnonzero(self.meta.center.ne(target));ids=np.flatnonzero(self.meta.center.eq(target))
            rng=np.random.default_rng(2026);groups={y:rng.permutation(ids[self.meta.iloc[ids].label.to_numpy()==y]) for y in [0,1]}
            anchors=np.r_[groups[0][:5],groups[1][:5]];pool={y:groups[y][5:] for y in [0,1]};batches=[]
            for rep in range(repeats):
                # Matched random priorities couple the three prevalence conditions.
                priorities={y:rng.permutation(pool[y]) for y in [0,1]}
                for pi in [.25,.5,.75]:
                    nc=int(32*pi);batch=np.r_[anchors,priorities[0][:32-nc-5],priorities[1][:nc-5]]
                    assert len(batch)==len(np.unique(batch))==32
                    batches.append((rep,pi,batch))
                    members.extend(dict(target=target,replicate=rep,target_case_fraction=pi,patient_id=self.meta.iloc[i].patient_id,label=int(self.meta.iloc[i].label),anchor=bool(i in anchors)) for i in batch)
            for method in METHODS:
                _,pre,model=self.fit(train,ids,method);frozen=pickle.dumps(pre)
                original=self.pred[self.pred.cohort_scope.eq(COMBINED)&self.pred.pipeline.eq(method)&self.pred.center.eq(target)]
                threshold=float(original.threshold.iloc[0])
                baseline=model.predict_proba(pre.transform(self.R[batches[1][2]],self.meta.iloc[batches[1][2]].center))[:10,1]
                max_invariance=0.
                if method!='InFold_ComBat':
                    independent=model.predict_proba(pre.transform(self.R[anchors],self.meta.iloc[anchors].center))[:,1]
                for rep,pi,batch in batches:
                    p=model.predict_proba(pre.transform(self.R[batch],self.meta.iloc[batch].center))[:10,1]
                    if method!='InFold_ComBat':
                        max_invariance=max(max_invariance,float(np.max(np.abs(p-independent))))
                    for j,idx in enumerate(anchors):
                        records.append(dict(target=target,method=display_method(method),replicate=rep,target_case_fraction=pi,
                            patient_id=self.meta.iloc[idx].patient_id,label=int(self.meta.iloc[idx].label),probability=float(p[j]),threshold=threshold,
                            reference_probability=float(baseline[j]),probability_shift=float(p[j]-baseline[j]),
                            threshold_crossing_flip=bool((p[j]>threshold)!=(baseline[j]>threshold)),
                            reference='replicate_0_case_fraction_0.50',experiment='fixed_real_anchors_50_companion_resamples'))
                assert pickle.dumps(pre)==frozen
                if method!='InFold_ComBat':assert max_invariance<1e-6
        data=pd.DataFrame(records);self.save(data,'Anchor_All_Predictions');self.save(pd.DataFrame(members),'Anchor_Real_Batch_Members')
        rows=[];effects=[]
        for (target,method),g in data.groupby(['target','method']):
            # Ranges over resamples per fixed prevalence isolate companion identity.
            for pi,h in g.groupby('target_case_fraction'):
                ranges=h.groupby('patient_id').probability.agg(lambda s:s.max()-s.min())
                shifts=h.probability_shift.abs();base=abs(h.probability_shift)
                rows.append(dict(target=target,method=method,case_fraction=pi,resamples=repeats,anchors=h.patient_id.nunique(),
                    median_abs_probability_shift=float(base.median()),q1_abs_probability_shift=float(base.quantile(.25)),q3_abs_probability_shift=float(base.quantile(.75)),
                    iqr_abs_probability_shift=float(base.quantile(.75)-base.quantile(.25)),max_abs_probability_shift=float(base.max()),
                    median_identity_range=float(ranges.median()),max_identity_range=float(ranges.max()),
                    threshold_crossing_flip_rate=float(h.threshold_crossing_flip.mean()),
                    fraction_anchors_ever_flip=float(h.groupby('patient_id').threshold_crossing_flip.any().mean())))
            means=g.groupby(['patient_id','target_case_fraction']).probability.mean().unstack()
            effects.append(dict(target=target,method=method,
                median_abs_mean_075_minus_025=float(abs(means[.75]-means[.25]).median()),
                max_abs_mean_075_minus_025=float(abs(means[.75]-means[.25]).max()),
                max_probability_range=float(g.groupby('patient_id').probability.agg(lambda s:s.max()-s.min()).max()),
                distinction='mean_prevalence_contrast_vs_within_fixed_prevalence_identity_spread; finite_real_patient_diagnostic_not_causal'))
        self.save(pd.DataFrame(effects),'Anchor_Prevalence_Contrasts')
        return self.save(pd.DataFrame(rows),'Anchor_Distribution_Summary')

    def spatial(self):
        anno=pd.read_parquet(Path(self.cfg['paths']['reference'])/'hg19_window_annotations.parquet').set_index('bin_id').loc[self.bins.bin_id]
        partition=pd.read_csv(self.root/'results/tables/Step7_Representation_Variance.csv');rows=[]
        for method,g in partition.groupby('representation'):
            v=g.set_index('bin_id').loc[self.bins.bin_id].eta2_batch.to_numpy()
            for variable in ['gc_fraction','distance_to_telomere_bp','distance_to_centromere_bp']:
                rows.append(dict(representation=method,annotation=variable,**block_rank_inference(anno[variable],v,self.bins.chrom)))
        result=pd.DataFrame(rows);result['q_block']=bh_fdr(result.block_permutation_p)
        return self.save(result,'Chromosome_Block_Associations')

    def dispersion(self):
        rows=[];individual=[]
        for family,x in [('raw_ratio',self.R),('length_frequency',self.lengths)]:
            z=((x-x.mean(axis=0))/np.maximum(x.std(axis=0),1e-8)).astype('float32')
            for term,col in [('study','center'),('disease','label')]:
                result,d=permdisp_euclidean(z,self.meta[col].to_numpy())
                rows.append(dict(family=family,term=term,**result))
                individual.extend(dict(family=family,term=term,patient_id=r.patient_id,group=str(getattr(r,col)),distance=float(v)) for r,v in zip(self.meta.itertuples(),d))
        result=pd.DataFrame(rows);result['q']=bh_fdr(result.p);self.save(pd.DataFrame(individual),'PERMDISP_Patient_Distances')
        return self.save(result,'PERMDISP_Dispersion')

    def run(self):
        start=time.monotonic()
        with threadpool_limits(limits=1):
            outputs={}
            for name,function in [('cohort',self.cohort_table),('deltas',self.pair_deltas),('thresholds',self.threshold_analysis),('sensitivity',self.sensitivities),('anchors',self.anchors),('spatial',self.spatial),('dispersion',self.dispersion)]:
                print('Running '+name,flush=True);outputs[name]=function();print('Completed '+name,flush=True)
        self.save(pd.DataFrame(self.fit_diagnostics),'Revision_Fit_Diagnostics')
        quality=dict(participants=len(self.meta),windows=len(self.bins),ratio_missing=int(np.isnan(self.R).sum()),ratio_total=int(self.R.size),
            length_missing=int(np.isnan(self.lengths).sum()),seconds=time.monotonic()-start,seed=2026,all_fits_converged=all(r['converged'] for r in self.fit_diagnostics),
            stage_missing=int(self.meta.stage.isna().sum()),matrix_sha256=self.data_signature,
            clinical_refitting='new source-only ElasticNet calibration and sensitivity analyses; old frozen predictions retained',
            imputation='per-feature source-training median; all-missing source columns use zero; post-transform medians; no target fit',
            elasticnet={'C':1.,'l1_ratio':.5,'solver':'saga','max_iter':20000,'tol':.001})
        write_json(quality,self.out/'reports/analysis_summary.json');return outputs
