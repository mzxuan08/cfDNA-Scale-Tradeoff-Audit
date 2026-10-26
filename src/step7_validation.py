"""Audit the real Step 7 deliverables and every nested source-patient boundary."""
import json
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
from .common import sha256,write_json
from .study import load_clinical
from .debiasing import ClinicalPreprocessor
from .evaluation import bh_fdr
from .multiscale_benchmark import PRIMARY,COMBINED,EXTERNAL,TUNED,read_protocol


def validate_step7_outputs(cfg):
    root=Path(cfg['paths']['results']);tables=root/'tables';protocol=read_protocol(cfg)
    meta,R,_,bins=load_clinical(cfg,True);known=set(meta.patient_id)
    primary=set(meta.loc[meta.center.isin(protocol['primary_studies']),'patient_id'])
    external=known-primary
    archive=json.loads((root/'archive_step06/asset_manifest.json').read_text())
    for path,digest in archive['inputs'].items():assert sha256(path)==digest,('Input changed',path)
    pred=pd.read_parquet(tables/'Step7_OOF_Predictions.parquet')
    assert pred.probability.between(0,1).all()
    assert not pred.duplicated(['cohort_scope','mode','pipeline','model','patient_id']).any()
    for (scope,mode,pipeline,model),g in pred.groupby(['cohort_scope','mode','pipeline','model']):
        expected=primary if scope==PRIMARY else (external if scope==EXTERNAL else known)
        assert set(g.patient_id)==expected,(scope,mode,pipeline,model)
    primary_loco=pred[pred.cohort_scope.eq(PRIMARY)&pred['mode'].eq('loco')]
    assert not primary_loco.pipeline.eq('InFold_ComBat').any() and not primary_loco.model.eq('dann').any()
    fallback=primary_loco[primary_loco.pipeline.isin(['Proposed_LSN_ICC',TUNED])]
    assert fallback.icc_bypassed.all() and fallback.alpha.eq(0).all()
    splits=pd.read_csv(tables/'Step7_Primary_Splits.csv');split_count=0
    for (mode,fold),g in splits.groupby(['mode','fold']):
        a=g[g.role.eq('train')];b=g[g.role.eq('test')]
        assert set(a.patient_id).isdisjoint(b.patient_id)
        assert set(g.patient_id)<=primary
        if mode=='loco':assert set(a.center).isdisjoint(b.center)
        split_count+=1
    assert split_count==7
    provenance=json.loads((root/'reports/step7_run_provenance.json').read_text())
    tuning={}
    for path in (root/'step7_checkpoints').glob('*.tuning.json'):
        item=json.loads(path.read_text())
        if item['signature']!=provenance['signature']:continue
        assert set(item['source_patients'])<=known
        tuning[tuple(item['source_patients'])]=item
        assert item['radius'] in protocol['radius_grid'] and item['alpha'] in protocol['alpha_grid']
    diag=json.loads((root/'reports/step7_fit_diagnostics.json').read_text());inner_folds=0;max_vram=0
    for item in diag:
        train=set(item['train_patients']);test=set(item['test_patients'])
        assert train.isdisjoint(test) and train|test<=known
        if item['mode']!='random':assert set(item['train_studies']).isdisjoint(item['test_studies'])
        if item['cohort_scope']==PRIMARY:assert train|test<=primary
        if item['cohort_scope']==EXTERNAL:assert train==primary and test==external
        assert item['inner_source_specificity']>=.95
        assert not item['max_iterations_reached']
        max_vram=max(max_vram,item['peak_vram_bytes'] or 0)
        if item['pipeline']==TUNED:
            selection=tuning[tuple(item['train_patients'])]
            assert selection['radius']==item['radius'] and selection['alpha']==item['alpha']
        for inner in item['inner_folds']:
            a=set(inner['train_ids']);b=set(inner['validation_ids'])
            assert a.isdisjoint(b) and a|b<=train and (a|b).isdisjoint(test)
            if item['pipeline']==TUNED:
                chosen=tuning[tuple(inner['train_ids'])]
                assert chosen['radius']==inner['radius'] and chosen['alpha']==inner['alpha']
            inner_folds+=1
    assert max_vram<4*1024**3
    # The reused external predictions have exactly the primary source population.
    oldsplit=pd.read_csv(root/'archive_step06/tables/Clinical_CV_Splits.csv')
    source=oldsplit[oldsplit['mode'].eq('loco')&oldsplit.fold.eq(2)&oldsplit.role.eq('train')]
    assert set(source.patient_id)==primary
    metrics=pd.read_csv(tables/'Table3_Benchmark_Metrics_Full.csv')
    assert len(metrics)==341 and metrics.status.eq('complete').sum()==311
    assert metrics[metrics.status.ne('complete')].auroc.isna().all()
    assert metrics[metrics.cohort_scope.eq(EXTERNAL)].auroc.isna().all()
    sig=pd.read_csv(tables/'Step7_Paired_Significance.csv')
    assert sig.loc[sig.comparator.eq('InFold_ComBat'),'p'].isna().all()
    for _,g in sig.groupby('cohort_scope'):np.testing.assert_allclose(g.q,bh_fdr(g.p),equal_nan=True)
    pair=pd.read_csv(tables/'Step7_AUC_Pair_Decomposition.csv')
    for pipeline,g in pair.groupby('pipeline'):
        row=metrics[metrics.cohort_scope.eq(COMBINED)&metrics.pipeline.eq(pipeline)&metrics.model.eq('elasticnet')&metrics.center.eq('AllCenters')].iloc[0]
        assert abs(g.pooled_pair_weight.sum()-1)<1e-12
        assert abs(g.pooled_auc_contribution.sum()-row.auroc)<1e-10
    stress=pd.read_csv(tables/'Step7_Prevalence_Stress_Summary.csv')
    independent=stress[stress.pipeline.ne('InFold_ComBat')]
    assert independent.max_anchor_probability_range.max()<1e-7
    members=pd.read_csv(tables/'Step7_Prevalence_Stress_Members.csv')
    assert set(members.patient_id)<=primary
    for (study,fraction),g in members.groupby(['target_study','target_prevalence']):
        assert g.patient_id.nunique()==32 and g.anchor.sum()==10
        assert abs(g.label.mean()-fraction)<1e-12
        assert set(g.patient_id)<=set(meta.loc[meta.center.eq(study),'patient_id'])
    f=ClinicalPreprocessor('Hybrid_MultiScale',bins).fit(R[meta.patient_id.isin(primary)])
    frozen=pickle.dumps(f);test=R[meta.patient_id.isin(external)]
    value=f.transform(test)
    assert value.dtype==np.float32 and value.shape==(len(external),2*len(bins))
    np.testing.assert_allclose(value[0],f.transform(test[:1])[0],atol=1e-7)
    assert frozen==pickle.dumps(f)
    text=(root/'MANUSCRIPT_FULL_TEXT.md').read_text(encoding='utf-8')
    for heading in ['## Abstract','## 1. Introduction','## 2. Materials and methods','## 3. Results','## 4. Discussion']:
        assert heading in text
    assert '### 4.1. Limitations' in text or '## 5. Limitations' in text
    assert '## 5. Conclusions' in text or '## 6. Conclusions' in text
    assert text in (root/'reports/MANUSCRIPT_MASTER_DRAFT.md').read_text(encoding='utf-8')
    assert chr(92)*2 not in text
    result=dict(status='passed',independent_patients=len(meta),primary_patients=len(primary),external_patients=len(external),
        outer_primary_split_count=split_count,new_outer_model_fits=len(diag),inner_calibration_folds_checked=inner_folds,
        source_only_tuning_jobs=len(tuning),oof_prediction_rows=len(pred),metric_rows=len(metrics),estimable_metric_rows=int(metrics.status.eq('complete').sum()),
        unavailable_metric_rows=int(metrics.status.ne('complete').sum()),max_torch_reserved_bytes=max_vram,
        predictor_state_frozen=True,raw_data_and_identity_hashes_unchanged=True,
        ordinary_DeLong_excludes_target_batch_adaptation=True,full_manuscript_word_count=len(text.split()))
    write_json(result,root/'reports/step7_validation.json');return result
