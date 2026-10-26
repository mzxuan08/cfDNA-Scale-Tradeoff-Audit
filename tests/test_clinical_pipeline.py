import pickle
import numpy as np
import pandas as pd
import pytest
from src.debiasing import ArmLocalLogRatio,ClinicalPreprocessor,healthy_center_components


def bins_fixture():
    return pd.DataFrame([dict(chrom='chr1',arm=arm,start=i*5000000,end=(i+1)*5000000)
                         for arm in ['p','q'] for i in range(5)])


def test_arm_lsn_excludes_focal_and_never_crosses_centromere():
    bins=bins_fixture();x=np.ones((1,10),dtype='float32');x[0,4]=100;x[0,5:]=1000
    f=ArmLocalLogRatio(bins,3).fit(x)
    assert 4 not in f.neighbors_[4] and all(j<5 for j in f.neighbors_[4])
    assert all(j>=5 for j in f.neighbors_[5])
    assert f.transform(x)[0,4]>6


def test_lsn_global_multiplicative_invariance_with_relative_epsilon():
    x=np.arange(1,31,dtype='float32').reshape(3,10)
    f=ArmLocalLogRatio(bins_fixture(),3).fit(x)
    np.testing.assert_allclose(f.transform(x),f.transform(x*np.array([.01,2,100],dtype='float32')[:,None]),atol=2e-6)


def test_arbitrary_additive_or_nonlinear_invariance_is_not_falsely_claimed():
    x=np.arange(1,31,dtype='float32').reshape(3,10)
    f=ArmLocalLogRatio(bins_fixture(),3).fit(x)
    assert not np.allclose(f.transform(x),f.transform(x+5))
    assert not np.allclose(f.transform(x),f.transform(x**2))


@pytest.mark.parametrize('pipeline',['Raw_Ratio','Global_Zscore','InFold_ComBat','Proposed_LSN_ICC','LSN_only','ICC_only','Hybrid_MultiScale'])
def test_clinical_transform_has_frozen_training_state(pipeline):
    rng=np.random.default_rng(2026);x=rng.lognormal(size=(48,10)).astype('float32')
    y=np.tile([0,1],24);centers=np.repeat(['A','B','C'],16)
    f=ClinicalPreprocessor(pipeline,bins_fixture()).fit(x,y,center_train=centers)
    snapshot=pickle.dumps(f)
    result=f.transform(x[:8],center_test=np.repeat('unseen',8))
    assert np.isfinite(result).all() and result.dtype==np.float32
    assert snapshot==pickle.dumps(f)
    if pipeline!='InFold_ComBat':
        np.testing.assert_allclose(result[0],f.transform(x[:1])[0])


def test_center_clustering_icc_is_not_a_repeatability_filter():
    rng=np.random.default_rng(2026)
    x=rng.normal(size=(60,3));x[20:40,0]+=20;x[40:,0]+=40
    icc,between,within,n=healthy_center_components(x,np.zeros(60),np.repeat(['A','B','C'],20))
    assert icc[0]>.99 and np.mean(icc[1:])<.2
    assert n.tolist()==[20,20,20]
