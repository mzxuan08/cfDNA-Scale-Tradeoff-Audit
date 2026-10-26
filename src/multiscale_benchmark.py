"""Step 7: prespecified bilateral testing and source-only multiscale selection.

No downloads or generated subjects. An unidentifiable one-source method is never
given an artificial center; ICC bypass and unavailable ComBat/DANN are explicit.
"""
import json
import pickle
import time
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from scipy.stats import beta,chi2_contingency
from sklearn.metrics import roc_auc_score
from threadpoolctl import threadpool_limits
from .common import fingerprint,sha256,write_json,write_parquet
from .study import load_clinical,export_table
from .debiasing import ClinicalPreprocessor
from .evaluation import (make_splits,source_calibration_splits,fit_clinical_model,
                         clinical_metric_ci,delong_test,bh_fdr)
from .models import make_model

PRIMARY='Primary_Bilateral'
COMBINED='Combined_LOCO'
EXTERNAL='External_HighBurden'
TUNED='Tuned_LSN_ICC'


def read_protocol(cfg):
    path=Path(cfg['_config_path']).parent/'step7_protocol.yaml'
    return yaml.safe_load(path.read_text(encoding='utf-8'))


def fit_preprocessor(pipeline,bins,x,y,centers,radius=3,alpha=.25):
    """Identity of absent ICC is reported, never estimated from the test study."""
    healthy_centers=np.unique(np.asarray(centers)[np.asarray(y)==0])
    effective=pipeline;icc_estimable=len(healthy_centers)>=2
    if pipeline in ('Proposed_LSN_ICC',TUNED):
        effective='Proposed_LSN_ICC' if icc_estimable else 'LSN_only'
        if not icc_estimable:alpha=0.
    if pipeline=='InFold_ComBat' and len(np.unique(centers))<2:
        raise ValueError('not_estimable_ComBat_requires_multiple_source_studies')
    pre=ClinicalPreprocessor(effective,bins,radius,alpha).fit(x,y,center_train=centers)
    info=dict(effective_pipeline=effective,radius=int(radius),alpha=float(alpha),
              icc_estimable=bool(icc_estimable),icc_bypassed=bool(pipeline in ('Proposed_LSN_ICC',TUNED) and not icc_estimable))
    return pre,info


def parameter_candidates(protocol,centers,y):
    alphas=protocol['alpha_grid'] if len(np.unique(np.asarray(centers)[np.asarray(y)==0]))>=2 else [0.]
    return [(k,a) for k in protocol['radius_grid'] for a in alphas]


def macro_eligible_auc(labels,scores,centers):
    values=[]
    for center in np.unique(centers):
        m=np.asarray(centers)==center;y=np.asarray(labels)[m]
        if min(np.sum(y==0),np.sum(y==1))>=2:values.append(roc_auc_score(y,np.asarray(scores)[m]))
    if not values:raise ValueError('No source validation stratum has two subjects per class')
    return float(np.mean(values))


def bilateral_splits(meta,mode,seed=2026):
    if mode=='random':return make_splits(meta,mode,5,seed)
    if meta.center.nunique()!=2 or meta.patient_id.duplicated().any():
        raise ValueError('Bilateral design requires two studies and unique patients')
    return [(np.flatnonzero(meta.center.ne(c)),np.flatnonzero(meta.center.eq(c))) for c in sorted(meta.center.unique())]


class Step7Engine:
    def __init__(self,cfg):
        self.cfg=cfg;self.protocol=read_protocol(cfg)
        self.meta,self.R,_,self.bins=load_clinical(cfg,True)
        if self.meta.patient_id.nunique()!=len(self.meta):raise ValueError('Patient isolation contract')
        self.root=Path(cfg['paths']['results']);self.cache=self.root/'step7_checkpoints';self.cache.mkdir(exist_ok=True)
        learning=['multiscale_benchmark.py','debiasing.py','evaluation.py','models.py','study.py']
        self.signature=fingerprint({'protocol':self.protocol,'matrix':sha256(Path(cfg['paths']['processed'])/cfg['extraction']['final_file']),
            'metadata':sha256(cfg['paths']['metadata']),'windows':self.bins.to_dict('list'),
            'code':{name:sha256(Path(__file__).parent/name) for name in learning},'evaluation':cfg['evaluation']})
        self.predictions=[];self.diagnostics=[];self.unavailable=[];self.splits=[]
        self.tuning_records=[]

    def tune(self,ids):
        key=fingerprint([self.signature,'tuning',self.meta.iloc[ids].patient_id.to_list()]);path=self.cache/(key+'.tuning.json')
        if path.exists():
            result=json.loads(path.read_text());self.tuning_records.append(result);return result['radius'],result['alpha']
        m=self.meta.iloc[ids];splits,rare=source_calibration_splits(m,3,self.cfg['seed'])
        records=[]
        with threadpool_limits(limits=1):
            for k,alpha in parameter_candidates(self.protocol,m.center,m.label):
                scores=np.full(len(ids),np.nan,dtype='float32')
                for tr,va in splits:
                    pre,_=fit_preprocessor(TUNED,self.bins,self.R[ids[tr]],m.label.to_numpy()[tr],m.center.to_numpy()[tr],k,alpha)
                    scores[va],_,_=fit_clinical_model('elasticnet',self.cfg,pre.transform(self.R[ids[tr]],m.center.to_numpy()[tr]),
                        m.label.to_numpy()[tr],m.center.to_numpy()[tr],pre.transform(self.R[ids[va]],m.center.to_numpy()[va]))
                valid=np.isfinite(scores)
                score=macro_eligible_auc(m.label.to_numpy()[valid],scores[valid],m.center.to_numpy()[valid])
                records.append(dict(radius=int(k),alpha=float(alpha),inner_macro_auroc=score))
        best=sorted(records,key=lambda r:(-r['inner_macro_auroc'],r['alpha'],abs(r['radius']-3),r['radius']))[0]
        result=dict(key=key,signature=self.signature,source_patients=m.patient_id.to_list(),source_studies=sorted(m.center.unique()),
                    radius=best['radius'],alpha=best['alpha'],candidates=records,
                    objective='source_only_3fold_macro_study_AUROC_min2_per_class',rare_train_only=m.iloc[rare].patient_id.to_list())
        write_json(result,path);self.tuning_records.append(result)
        print(f"Nested selection n={len(ids)} k={best['radius']} alpha={best['alpha']}",flush=True)
        return best['radius'],best['alpha']

    def fit_fold(self,scope,mode,fold,train,test,pipeline,model):
        a=self.meta.iloc[train];b=self.meta.iloc[test]
        if set(a.patient_id)&set(b.patient_id):raise ValueError('Patient leakage')
        if mode!='random' and set(a.center)&set(b.center):raise ValueError('Study leakage')
        if len(a.center.unique())<2 and (pipeline=='InFold_ComBat' or model=='dann'):
            self.unavailable.append(dict(cohort_scope=scope,mode=mode,fold=int(fold),pipeline=pipeline,model=model,
                center=b.center.iloc[0],status='not_estimable_single_source_'+('ComBat' if pipeline=='InFold_ComBat' else 'DANN')))
            return
        key=fingerprint([self.signature,scope,mode,fold,pipeline,model]);path=self.cache/(key+'.parquet')
        if path.exists() and path.with_suffix('.json').exists():
            stamp=json.loads(path.with_suffix('.json').read_text())
            if sha256(path)!=stamp['sha256']:raise ValueError('Step 7 checkpoint checksum')
            self.predictions.append(pd.read_parquet(path));self.diagnostics.append(stamp['diagnostics']);return
        k,alpha=self.tune(train) if pipeline==TUNED else (self.protocol['hybrid_radius'],.25)
        pre,info=fit_preprocessor(pipeline,self.bins,self.R[train],a.label,a.center,k,alpha)
        calibration=np.full(len(train),np.nan,dtype='float32');inner,rare=source_calibration_splits(a,3,self.cfg['seed'])
        inner_notes=[]
        with threadpool_limits(limits=1):
            for ti,vi in inner:
                ii,jj=train[ti],train[vi]
                # Selection is repeated WITHOUT inner calibration subjects.
                ik,ia=self.tune(ii) if pipeline==TUNED else (k,alpha)
                ip,note=fit_preprocessor(pipeline,self.bins,self.R[ii],self.meta.iloc[ii].label,self.meta.iloc[ii].center,ik,ia)
                calibration[vi],_,_=fit_clinical_model(model,self.cfg,ip.transform(self.R[ii],self.meta.iloc[ii].center),
                    self.meta.iloc[ii].label.to_numpy(),self.meta.iloc[ii].center.to_numpy(),ip.transform(self.R[jj],self.meta.iloc[jj].center))
                inner_notes.append(dict(train_ids=self.meta.iloc[ii].patient_id.to_list(),validation_ids=self.meta.iloc[jj].patient_id.to_list(),**note))
            controls=np.isfinite(calibration)&a.label.eq(0).to_numpy()
            threshold=float(np.quantile(calibration[controls],.95,method='higher'))
            probability,importance,model_info=fit_clinical_model(model,self.cfg,pre.transform(self.R[train],a.center),
                a.label.to_numpy(),a.center.to_numpy(),pre.transform(self.R[test],b.center))
        prediction=b[['patient_id','sample_id','center','label']].copy()
        prediction['cohort_scope']=scope;prediction['mode']=mode;prediction['fold']=fold
        prediction['pipeline']=pipeline;prediction['model']=model;prediction['probability']=probability;prediction['threshold']=threshold
        prediction['protocol']='transductive_unlabeled_target_batch' if pipeline=='InFold_ComBat' else 'inductive_no_target_population_fit'
        for name,value in info.items():prediction[name]=value
        write_parquet(prediction,path)
        diagnostic=dict(cohort_scope=scope,mode=mode,fold=fold,pipeline=pipeline,model=model,
            train_patients=a.patient_id.to_list(),test_patients=b.patient_id.to_list(),train_studies=sorted(a.center.unique()),test_studies=sorted(b.center.unique()),
            threshold=threshold,inner_control_n=int(controls.sum()),inner_source_specificity=float((calibration[controls]<=threshold).mean()),
            rare_train_only=a.iloc[rare].patient_id.to_list(),inner_folds=inner_notes,**info,**model_info)
        write_json(dict(sha256=sha256(path),signature=self.signature,diagnostics=diagnostic),path.with_suffix('.json'))
        self.predictions.append(prediction);self.diagnostics.append(diagnostic)
        print(f'Finished {scope} {mode} fold={fold} {pipeline} {model}',flush=True)

    def add_archived_combined(self):
        folder=self.root/'archive_step06';manifest={k.replace('\\','/'):v for k,v in json.loads((folder/'snapshot_manifest.json').read_text()).items()}
        name='tables/Clinical_OOF_Predictions.parquet';path=folder/name
        if sha256(path)!=manifest[name]:raise ValueError('Step 6 archive checksum')
        old=pd.read_parquet(path);old=old[old['mode'].eq('loco')].copy()
        # Algorithm changes only add a new sample-wise branch; source subjects,
        # all old model settings, thresholds and preprocessing remain identical.
        if sha256(Path(__file__).parent/'models.py')!=sha256(folder/'source/src/models.py'):raise ValueError('Cannot reuse old fitted-model predictions after a model change')
        expected=set(self.meta.patient_id)
        for _,g in old.groupby(['pipeline','model']):
            if set(g.patient_id)!=expected:raise ValueError('Archived participant mismatch')
        old['cohort_scope']=COMBINED;old['effective_pipeline']=old.pipeline
        old['radius']=3;old['alpha']=np.where(old.pipeline.eq('Proposed_LSN_ICC'),.25,0.)
        old['icc_estimable']=True;old['icc_bypassed']=False
        self.predictions.append(old)
        external=old[old.center.eq(self.protocol['external_study'])].copy()
        external['cohort_scope']=EXTERNAL;external['mode']='external';external['fold']=0
        self.predictions.append(external)
        write_json({'source':str(path),'sha256':sha256(path),'scope':'unchanged four-pipeline three-study LOCO; external fold source patients equal primary cohort',
                    'model_file_sha256':sha256(Path(__file__).parent/'models.py')},self.root/'reports/step7_archived_prediction_reuse.json')

    def run(self):
        primary=np.flatnonzero(self.meta.center.isin(self.protocol['primary_studies']))
        external=np.flatnonzero(self.meta.center.eq(self.protocol['external_study']))
        if len(primary)!=548 or len(external)!=37:raise ValueError('Locked cohort size differs from the existing verified cohort')
        tasks=[]
        for mode in ['random','loco']:
            for fold,(a,b) in enumerate(bilateral_splits(self.meta.iloc[primary].reset_index(drop=True),mode,self.cfg['seed'])):
                train,test=primary[a],primary[b]
                for role,ids in [('train',train),('test',test)]:
                    self.splits.extend(dict(cohort_scope=PRIMARY,mode=mode,fold=fold,role=role,patient_id=self.meta.iloc[i].patient_id,center=self.meta.iloc[i].center) for i in ids)
                for pipeline in self.protocol['pipelines']:
                    for model in self.protocol['models']:tasks.append((PRIMARY,mode,fold,train,test,pipeline,model))
                tasks.append((PRIMARY,mode,fold,train,test,TUNED,'elasticnet'))
        self.add_archived_combined()
        # The external fit uses ALL primary subjects and no external labels.
        for pipeline in ['Hybrid_MultiScale',TUNED]:
            for model in (self.protocol['models'] if pipeline!='Tuned_LSN_ICC' else ['elasticnet']):
                tasks.append((EXTERNAL,'external',0,primary,external,pipeline,model))
        for fold,(train,test) in enumerate(make_splits(self.meta,'loco',5,self.cfg['seed'])):
            for pipeline in ['Hybrid_MultiScale',TUNED]:
                for model in (self.protocol['models'] if pipeline!=TUNED else ['elasticnet']):
                    if self.meta.iloc[test].center.eq(self.protocol['external_study']).all():continue
                    tasks.append((COMBINED,'loco',fold,train,test,pipeline,model))
        for task in tasks:self.fit_fold(*task)
        pred=pd.concat(self.predictions,ignore_index=True)
        add=pred[pred.cohort_scope.eq(EXTERNAL)&pred.pipeline.isin(['Hybrid_MultiScale',TUNED])].copy()
        add['cohort_scope']=COMBINED;add['mode']='loco';add['fold']=2;pred=pd.concat([pred,add],ignore_index=True)
        keys=['cohort_scope','mode','pipeline','model','patient_id']
        if pred.duplicated(keys).any():raise ValueError('Duplicate outer predictions')
        write_parquet(pred,self.root/'tables/Step7_OOF_Predictions.parquet')
        export_table(pd.DataFrame(self.splits),self.cfg,'Step7_Primary_Splits')
        write_json(self.diagnostics,self.root/'reports/step7_fit_diagnostics.json')
        export_table(pd.DataFrame(self.unavailable),self.cfg,'Step7_Unavailable_Methods')
        tune=[]
        for file in self.cache.glob('*.tuning.json'):
            item=json.loads(file.read_text())
            if item['signature']!=self.signature:continue
            tune.extend(dict(key=item['key'],training_n=len(item['source_patients']),training_studies='; '.join(item['source_studies']),
                             selected_radius=item['radius'],selected_alpha=item['alpha'],**row) for row in item['candidates'])
        export_table(pd.DataFrame(tune),self.cfg,'Step7_Nested_Parameter_Search')
        write_json({'signature':self.signature,'protocol':self.protocol,'independent_patients':len(self.meta),
                    'input_feature_sha256':sha256(Path(self.cfg['paths']['processed'])/self.cfg['extraction']['final_file'])},self.root/'reports/step7_run_provenance.json')
        return pred


def exact_interval(success,n):
    return (float(beta.ppf(.025,success,n-success+1)) if success else 0.,
            float(beta.ppf(.975,success+1,n-success)) if success<n else 1.)


def summarize_step7(cfg,pred):
    protocol=read_protocol(cfg);rows=[]
    for (scope,mode,pipeline,model),g in pred.groupby(['cohort_scope','mode','pipeline','model']):
        groups=[('AllCenters',g),*list(g.groupby('center'))] if scope!=EXTERNAL else [(protocol['external_study'],g)]
        for center,sub in groups:
            result=clinical_metric_ci(sub,protocol['bootstrap_replicates'],cfg['seed'])
            h=sub[sub.label.eq(0)];c=sub[sub.label.eq(1)]
            result['specificity_exact_low'],result['specificity_exact_high']=exact_interval(int((h.probability<=h.threshold).sum()),len(h))
            result['sensitivity_exact_low'],result['sensitivity_exact_high']=exact_interval(int((c.probability>c.threshold).sum()),len(c))
            result['specificity_deviation_pp']=(result['specificity95']-.95)*100
            if len(h)<2:
                result['exploratory_single_control_auroc']=result['auroc']
                for key in ['auroc','auroc_low','auroc_high','auprc_low','auprc_high']:result[key]=np.nan
                result['auc_scope']='not_reported_one_control_sensitivity_focus'
            else:result['auc_scope']='two_class_benchmark'
            rows.append(dict(cohort_scope=scope,mode=mode,pipeline=pipeline,model=model,center=center,
                protocol=sub.protocol.iloc[0],effective_pipeline=sub.effective_pipeline.iloc[0],icc_bypassed=bool(sub.icc_bypassed.any()),
                radius=';'.join(map(str,sorted(sub.radius.unique()))),alpha=';'.join(map(str,sorted(sub.alpha.unique()))),**result))
    # Explicit unavailable rows keep the original requested method identity intact.
    unavailable=pd.read_csv(Path(cfg['paths']['results'])/'tables/Step7_Unavailable_Methods.csv')
    for (scope,mode,pipeline,model),g in unavailable.groupby(['cohort_scope','mode','pipeline','model']):
        for center in ['AllCenters',*protocol['primary_studies']]:
            rows.append(dict(cohort_scope=scope,mode=mode,pipeline=pipeline,model=model,center=center,
                status='not_estimable_single_source',reason=';'.join(sorted(g.status.unique())),n=548 if center=='AllCenters' else int((pred[(pred.cohort_scope==PRIMARY)&(pred['mode']==mode)&pred.pipeline.eq('Raw_Ratio')&pred.model.eq('elasticnet')].center==center).sum())))
    metrics=pd.DataFrame(rows)
    macros=metrics[metrics.center.isin(protocol['primary_studies'])].groupby(['cohort_scope','mode','pipeline','model']).auroc.mean().to_dict()
    metrics['macro_bilateral_auroc']=[macros.get((r.cohort_scope,r.mode,r.pipeline,r.model),np.nan) for r in metrics.itertuples()]
    export_table(metrics,cfg,'Table3_Benchmark_Metrics_Full')
    tests=[]
    for scope in [PRIMARY,COMBINED]:
        for (model,center),g in pred[pred.cohort_scope.eq(scope)&pred['mode'].eq('loco')].groupby(['model','center']):
            wide=g.pivot(index=['patient_id','label'],columns='pipeline',values='probability').reset_index()
            comparisons=[('Hybrid_MultiScale',c) for c in ['Raw_Ratio','Global_Zscore','Proposed_LSN_ICC','InFold_ComBat']]
            if model=='elasticnet':comparisons += [(TUNED,c) for c in ['Raw_Ratio','Global_Zscore','Proposed_LSN_ICC','Hybrid_MultiScale']]
            for proposed,comparator in comparisons:
                result={'p':np.nan,'delta':np.nan};status='method_not_estimable_in_this_design'
                if proposed in wide and comparator in wide:
                    try:result=delong_test(wide.label,wide[proposed],wide[comparator]);status='fixed_models_independent_heldout_study'
                    except ValueError as exc:status=str(exc)
                tests.append(dict(cohort_scope=scope,model=model,center=center,method=proposed,comparator=comparator,status=status,**result))
    significance=pd.DataFrame(tests);significance['q']=np.nan
    for _,idx in significance.groupby('cohort_scope').groups.items():significance.loc[idx,'q']=bh_fdr(significance.loc[idx,'p'])
    export_table(significance,cfg,'Step7_Paired_Significance')
    return metrics,significance


def prevalence_diagnostics(cfg,pred):
    meta,R,_,bins=load_clinical(cfg,True);root=Path(cfg['paths']['results']);protocol=read_protocol(cfg)
    counts=meta.groupby('center',as_index=False).agg(n=('patient_id','size'),cancer=('label','sum'))
    counts['healthy']=counts.n-counts.cancer;counts['prevalence']=counts.cancer/counts.n
    export_table(counts,cfg,'Step7_Prevalence_By_Study')
    chi,p,df,expected=chi2_contingency(counts[['healthy','cancer']].to_numpy(),correction=False)
    export_table(pd.DataFrame([dict(chi2=chi,df=df,p=p,cramers_v=np.sqrt(chi/counts.n.sum()),n=int(counts.n.sum()))]),cfg,'Step7_Prevalence_Association')
    pairs=[];combined=pred[pred.cohort_scope.eq(COMBINED)&pred.model.eq('elasticnet')]
    for pipeline,g in combined.groupby('pipeline'):
        cases=g[g.label.eq(1)];controls=g[g.label.eq(0)];den=len(cases)*len(controls)
        for a,x in cases.groupby('center'):
            for b,y in controls.groupby('center'):
                diff=x.probability.to_numpy()[:,None]-y.probability.to_numpy()[None,:]
                wins=float(np.sum(diff>0)+.5*np.sum(diff==0));total=diff.size
                pairs.append(dict(pipeline=pipeline,case_study=a,control_study=b,pair_type='within_study' if a==b else 'cross_study_different_fitted_model',
                                  pair_n=total,correct_order_fraction=wins/total,pooled_auc_contribution=wins/den,pooled_pair_weight=total/den))
    export_table(pd.DataFrame(pairs),cfg,'Step7_AUC_Pair_Decomposition')
    scores=[];members=[];summary=[]
    # Real held-out participants compose diagnostic target batches. Labels are
    # used ONLY for this post-evaluation stress experiment, never for fitting.
    for target in protocol['primary_studies']:
        train=np.flatnonzero(meta.center.ne(target));target_ids=np.flatnonzero(meta.center.eq(target))
        rng=np.random.default_rng(cfg['seed']);groups={label:rng.permutation(target_ids[meta.iloc[target_ids].label.to_numpy()==label]) for label in [0,1]}
        anchors=np.r_[groups[0][:5],groups[1][:5]];background={label:groups[label][5:] for label in [0,1]}
        for pipeline in protocol['pipelines']:
            pre,_=fit_preprocessor(pipeline,bins,R[train],meta.iloc[train].label,meta.iloc[train].center,3,.25)
            model=make_model('elasticnet',cfg)
            with threadpool_limits(limits=1):model.fit(pre.transform(R[train],meta.iloc[train].center),meta.iloc[train].label)
            frozen=pickle.dumps(pre);anchor_values=[]
            for prevalence in protocol['prevalence_stress_fractions']:
                ncase=int(32*prevalence);nhealthy=32-ncase
                batch=np.r_[anchors,background[0][:nhealthy-5],background[1][:ncase-5]]
                if len(batch)!=32 or len(np.unique(batch))!=32:raise ValueError('Real-patient mixture cannot be constructed')
                prob=model.predict_proba(pre.transform(R[batch],meta.iloc[batch].center))[:len(anchors),1]
                if frozen!=pickle.dumps(pre):raise ValueError('Transform modified source-fitted state')
                anchor_values.append(prob)
                for idx,probability in zip(anchors,prob):scores.append(dict(pipeline=pipeline,target_study=target,patient_id=meta.iloc[idx].patient_id,
                    label=int(meta.iloc[idx].label),target_prevalence=prevalence,probability=float(probability),diagnostic_only=True))
                if pipeline==protocol['pipelines'][0]:
                    members.extend(dict(target_study=target,target_prevalence=prevalence,patient_id=meta.iloc[idx].patient_id,label=int(meta.iloc[idx].label),anchor=bool(idx in anchors)) for idx in batch)
            anchor_values=np.stack(anchor_values)
            summary.append(dict(pipeline=pipeline,target_study=target,anchors=len(anchors),batch_n=32,
                max_anchor_probability_range=float(np.ptp(anchor_values,axis=0).max()),mean_anchor_probability_range=float(np.ptp(anchor_values,axis=0).mean()),
                source_training_studies=';'.join(sorted(meta.iloc[train].center.unique())),
                scope='diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models'))
    export_table(pd.DataFrame(scores),cfg,'Step7_Prevalence_Stress_Anchor_Predictions')
    export_table(pd.DataFrame(members),cfg,'Step7_Prevalence_Stress_Members')
    export_table(pd.DataFrame(summary),cfg,'Step7_Prevalence_Stress_Summary')
    return counts,pd.DataFrame(pairs),pd.DataFrame(summary)


def run_step7(cfg,render_only=False):
    start=time.perf_counter();root=Path(cfg['paths']['results'])
    if render_only:pred=pd.read_parquet(root/'tables/Step7_OOF_Predictions.parquet')
    else:pred=Step7Engine(cfg).run()
    metrics,significance=summarize_step7(cfg,pred)
    counts,pairs,stress=prevalence_diagnostics(cfg,pred)
    from .multiscale_plotting import step7_figures
    from .multiscale_manuscript import compile_step7_manuscript
    step7_figures(cfg,pred,metrics,significance,counts,pairs,stress)
    compile_step7_manuscript(cfg,pred,metrics,significance,counts,pairs,stress)
    report=dict(status='complete',runtime_seconds=round(time.perf_counter()-start,3),render_only=render_only,
                independent_patients=int(counts.n.sum()),counts_per_study=counts.to_dict('records'),
                primary_patients=int(counts[counts.center.isin(read_protocol(cfg)['primary_studies'])].n.sum()),
                total_prediction_rows=len(pred),metric_rows=len(metrics),
                cohort_patient_ids_sha256=sha256(root/'tables/Analysis_Participants.parquet'),
                no_additional_clinical_data_downloaded=True,no_synthetic_subjects=True)
    write_json(report,root/'reports/step_07_report.json')
    from .manuscript import refresh_clinical_manifest
    refresh_clinical_manifest(cfg)
    return report
