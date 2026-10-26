import json
import pickle
import numpy as np
import pandas as pd
import pytest
from src.debiasing import ClinicalPreprocessor,ArmLocalLogRatio
from src.multiscale_benchmark import (Step7Engine,fit_preprocessor,parameter_candidates,
    read_protocol,PRIMARY,TUNED,exact_interval,bilateral_splits)
from src.common import load_config


def schema():
    return pd.DataFrame({'chrom':['chr1']*10,'arm':['p']*5+['q']*5,'start':np.arange(10)*5000000})


def test_hybrid_preserves_both_observable_scales_and_is_batch_independent():
    rng=np.random.default_rng(2026);x=rng.lognormal(size=(40,10)).astype('float32')
    f=ClinicalPreprocessor('Hybrid_MultiScale',schema()).fit(x)
    frozen=pickle.dumps(f);result=f.transform(x[:7])
    recovered=f.scaler_.inverse_transform(result)
    global_z=(x[:7]-x[:7].mean(axis=1,keepdims=True))/x[:7].std(axis=1,keepdims=True)
    local=ArmLocalLogRatio(schema()).fit(x).transform(x[:7])
    assert result.shape==(7,20) and result.dtype==np.float32
    np.testing.assert_allclose(recovered[:,:10],global_z,atol=2e-6)
    np.testing.assert_allclose(recovered[:,10:],local,atol=2e-6)
    np.testing.assert_allclose(result[0],f.transform(x[:1])[0])
    assert pickle.dumps(f)==frozen


def test_single_source_does_not_invent_icc_or_combat_parameters():
    rng=np.random.default_rng(2026);x=rng.lognormal(size=(40,10)).astype('float32')
    y=np.tile([0,1],20);centers=np.repeat('source_A',40)
    pre,note=fit_preprocessor('Proposed_LSN_ICC',schema(),x,y,centers)
    assert note['icc_bypassed'] and not note['icc_estimable'] and note['alpha']==0
    assert pre.filter_ is None and pre.transform(x).shape==x.shape
    with pytest.raises(ValueError,match='multiple_source_studies'):
        fit_preprocessor('InFold_ComBat',schema(),x,y,centers)
    assert {a for k,a in parameter_candidates(read_protocol(load_config()),centers,y)}=={0.}


def test_nested_selection_excludes_outer_test_and_inner_calibration_patients(tmp_path,monkeypatch):
    import src.multiscale_benchmark as module
    e=Step7Engine.__new__(Step7Engine);e.cfg={'seed':2026};e.protocol={'hybrid_radius':3}
    e.signature='test_only';e.cache=tmp_path;e.bins=schema()
    e.meta=pd.DataFrame({'patient_id':[f'patient{i}' for i in range(60)],'sample_id':[f's{i}' for i in range(60)],
                         'center':['source_A']*40+['test_B']*20,'label':np.tile([0,1],30)})
    e.R=np.random.default_rng(2026).lognormal(size=(60,10)).astype('float32')
    e.predictions=[];e.diagnostics=[];e.unavailable=[];selected=[]
    def tune(ids):selected.append(set(e.meta.iloc[ids].patient_id));return 3,0.
    e.tune=tune
    def fake_model(name,cfg,x,y,centers,test):
        return np.linspace(.1,.9,len(test),dtype='float32'),None,{'peak_vram_bytes':None,'max_iterations_reached':False}
    monkeypatch.setattr(module,'fit_clinical_model',fake_model)
    e.fit_fold(PRIMARY,'loco',0,np.arange(40),np.arange(40,60),TUNED,'elasticnet')
    assert len(selected)==4
    heldout=set(e.meta.iloc[40:].patient_id)
    assert all(s.isdisjoint(heldout) for s in selected)
    for chosen,fold in zip(selected[1:],e.diagnostics[0]['inner_folds']):
        assert chosen==set(fold['train_ids']) and chosen.isdisjoint(fold['validation_ids'])
    assert e.diagnostics[0]['icc_bypassed']


def test_singleton_control_uncertainty_is_not_a_degenerate_bootstrap_interval():
    lo,hi=exact_interval(1,1)
    assert lo==pytest.approx(.025) and hi==1


def test_bilateral_split_uses_only_two_actual_studies():
    m=pd.DataFrame({'patient_id':['a','b','c','d'],'center':['A','A','B','B'],'label':[0,1,0,1]})
    splits=bilateral_splits(m,'loco')
    assert len(splits)==2
    for train,test in splits:
        assert set(m.iloc[train].center).isdisjoint(m.iloc[test].center)
        assert set(m.iloc[train].patient_id).isdisjoint(m.iloc[test].patient_id)


def test_affine_harmonization_preserves_within_study_standardized_effects():
    from src.prevalence_analysis import cohen_d_by_feature
    x=np.random.default_rng(2026).normal(size=(40,4));y=np.tile([0,1],20)
    x[y==1]+=np.array([1,2,3,4])
    transformed=x*np.array([.2,2,4,10])+np.array([100,-20,.3,4])
    np.testing.assert_allclose(cohen_d_by_feature(x,y),cohen_d_by_feature(transformed,y),atol=1e-12)


def test_target_batch_adaptation_is_excluded_from_fixed_rule_delong_family(tmp_path):
    from src.prevalence_analysis import qualify_paired_inference
    f=pd.DataFrame({'cohort_scope':[PRIMARY]*3,'method':['Hybrid_MultiScale']*3,
                    'comparator':['InFold_ComBat','Global_Zscore','Raw_Ratio'],
                    'p':[.001,.05,.02],'status':['computed']*3})
    result=qualify_paired_inference({'paths':{'results':str(tmp_path)}},f)
    assert np.isnan(result.loc[0,'p']) and np.isnan(result.loc[0,'q'])
    assert 'transductive' in result.loc[0,'status']
    np.testing.assert_allclose(result.loc[[1,2],'q'],[.05,.04])
