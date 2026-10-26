"""Inference and companion-invariance contracts for the major revision."""
import pickle
import numpy as np
import pandas as pd
import pytest
from src.evaluation import threshold_bootstrap,paired_auc_difference_ci
from src.debiasing import ClinicalPreprocessor

def test_calibration_bootstrap_uses_source_scores_and_does_not_change_state():
    source=np.array([.1,.2,.3,.4]);before=source.copy()
    t,s,c=threshold_bootstrap(source,[.15,.35,.8,.9],[0,0,1,1],100,2026)
    np.testing.assert_array_equal(source,before)
    assert set(t)<=set(source) and np.all((s>=0)&(s<=1)) and np.all(c==1)
    t2,_,_=threshold_bootstrap(source,[.15,.35,.8,.9],[1,1,0,0],100,2026)
    np.testing.assert_array_equal(t,t2)

def test_paired_delta_is_zero_for_identical_scores():
    r=paired_auc_difference_ci([0,0,1,1],[.1,.3,.5,.7],[.1,.3,.5,.7],100)
    assert r['delta_auroc']==r['delta_low']==r['delta_high']==0

@pytest.mark.parametrize('method',['Raw_Ratio','Global_Zscore','LSN_only','Proposed_LSN_ICC','Hybrid_MultiScale'])
def test_patientwise_invariance_to_companion_identity_and_frozen_state(method):
    rng=np.random.default_rng(2026)
    bins=pd.DataFrame({'chrom':['chr1']*12,'arm':['p']*6+['q']*6,'start':np.arange(12)*5000000})
    train=rng.uniform(.1,2,(40,12)).astype('float32');y=np.tile([0,1],20);centers=np.repeat(['A','B'],20)
    pre=ClinicalPreprocessor(method,bins).fit(train,y,centers)
    state=pickle.dumps(pre);anchor=train[:1]
    alone=pre.transform(anchor,['C'])
    together=pre.transform(np.concatenate([anchor,train[10:30]*3]),['C']*21)[:1]
    np.testing.assert_allclose(alone,together,rtol=0,atol=0)
    assert state==pickle.dumps(pre)

    # A frozen fitted classifier must retain the same patient probability too.
    from sklearn.linear_model import LogisticRegression
    clf=LogisticRegression(random_state=2026).fit(pre.transform(train,centers),y)
    np.testing.assert_array_equal(clf.predict_proba(alone),clf.predict_proba(together))


def test_permdisp_distance_and_F_agree_with_direct_euclidean_centroids():
    from src.major_revision import permdisp_euclidean
    from scipy.stats import f_oneway
    x=np.array([[0,0],[1,0],[3,2],[5,1],[10,0],[11,0],[18,5],[15,1]],float)
    groups=np.repeat(['a','b'],4)
    expected=np.concatenate([np.linalg.norm(v-v.mean(axis=0),axis=1)*np.sqrt(4/3) for v in [x[:4],x[4:]]])
    record,d=permdisp_euclidean(x,groups,permutations=19)
    np.testing.assert_allclose(d,expected,atol=1e-12)
    assert record['F']==pytest.approx(f_oneway(expected[:4],expected[4:]).statistic)
    assert 1/20<=record['p']<=1
    shifted=x.copy();shifted[4:]+=1000
    other,_=permdisp_euclidean(shifted,groups,permutations=19)
    assert other['F']==pytest.approx(record['F'],rel=1e-7)
    assert other['p']==record['p']


def test_chromosome_block_inference_preserves_correlated_cluster_size():
    from src.major_revision import block_rank_inference
    x=np.arange(30);blocks=np.repeat(np.arange(6),5)
    result=block_rank_inference(x,x,blocks,repeats=30,permutations=31)
    assert result['blocks']==6 and result['windows']==30
    assert result['rho']==pytest.approx(1)
    assert result['block_ci_low']==pytest.approx(1)
    assert result['block_ci_high']==pytest.approx(1)
    assert result['block_permutation_p']>=1/32


def test_real_revision_calibration_is_source_only_and_study_heldout():
    from pathlib import Path
    folder=Path(__file__).resolve().parents[1]/'results/step09_major/tables'
    if not folder.exists():pytest.skip('real-patient revision outputs are not part of minimal code distribution')
    splits=pd.read_parquet(folder/'Source_Calibration_Splits.parquet')
    assert len(splits)>0
    for (outer,scheme,method,fold),g in splits.groupby(['outer','scheme','method','fold']):
        train=g[g.role.eq('train')];valid=g[g.role.eq('validation')]
        assert not set(train.patient_id)&set(valid.patient_id)
        target=outer.split('|')[-1]
        assert target not in set(g.study),outer
        if scheme=='LeaveOneSourceStudyOut':assert not set(train.study)&set(valid.study)


def test_real_anchors_have_exact_composition_and_invariant_inductive_probabilities():
    from pathlib import Path
    folder=Path(__file__).resolve().parents[1]/'results/step09_major/tables'
    if not folder.exists():pytest.skip('clinical diagnostic results not bundled')
    members=pd.read_parquet(folder/'Anchor_Real_Batch_Members.parquet')
    groups=members.groupby(['target','replicate','target_case_fraction'])
    assert len(groups)==300
    for (_,_,fraction),g in groups:
        assert len(g)==g.patient_id.nunique()==32
        assert g.label.mean()==fraction
        assert g.anchor.sum()==10
    p=pd.read_parquet(folder/'Anchor_All_Predictions.parquet')
    assert p[p.method.ne('EB-TargetMap')].probability_shift.abs().max()<1e-6
    assert p[p.method.eq('EB-TargetMap')&p.target_case_fraction.eq(.5)].threshold_crossing_flip.any()


def test_real_primary_multiple_comparison_family_is_complete():
    from pathlib import Path
    from src.evaluation import bh_fdr
    path=Path(__file__).resolve().parents[1]/'results/step09_major/tables/Primary_BH_Hypothesis_Family.csv'
    if not path.exists():pytest.skip('clinical inference tables not bundled')
    table=pd.read_csv(path)
    assert len(table)==table.hypothesis_id.nunique()==38
    np.testing.assert_allclose(table.q,bh_fdr(table.p),atol=1e-14)


def test_revision_manuscript_math_has_no_escaped_control_characters():
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    path=root/'results/MANUSCRIPT_FULL_TEXT.md'
    if not path.exists() or not (root/'results/step09_major').exists():pytest.skip('publication artifacts not bundled')
    text=path.read_text(encoding='utf-8')
    assert all(ord(c)>=32 or c=='\n' for c in text)
    assert r'\frac{S_{si}}{L_{si}+1}' in text
    assert r'g_\theta(f_\theta(x_i))' in text
    assert text.count('Declaration of generative AI')==1
    assert 'Proposed' not in text and 'InFold_ComBat' not in text
    assert '### 2.8.' not in text
    assert '0.384' not in text.split('## 1. Introduction')[1].split('## 2.')[0]
