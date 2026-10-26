import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
from src.evaluation import weighted_clinical_metrics,clinical_metric_ci,source_calibration_splits


def test_weighted_bootstrap_statistics_equal_explicit_resampled_subjects():
    y=np.array([0,1,0,1,0,1]);p=np.array([.2,.8,.5,.5,.1,.9]);threshold=np.full(6,.6)
    weights=np.array([[1,2,3,4,2,1],[3,1,1,2,4,3]])
    result=weighted_clinical_metrics(y,p,threshold,weights)
    for j,w in enumerate(weights):
        ids=np.repeat(np.arange(6),w)
        np.testing.assert_allclose(result['auroc'][j],roc_auc_score(y[ids],p[ids]))
        np.testing.assert_allclose(result['auprc'][j],average_precision_score(y[ids],p[ids]))
        np.testing.assert_allclose(result['brier'][j],brier_score_loss(y[ids],p[ids]))


def test_bootstrap_metric_scope_and_fixed_training_threshold():
    f=pd.DataFrame(dict(patient_id=[f'p{i}' for i in range(20)],center=['A']*10+['B']*10,
                        label=[0,1]*10,probability=[.2,.8]*10,threshold=[.5]*20))
    result=clinical_metric_ci(f,repeats=100)
    assert result['auroc']==result['auroc_low']==result['auroc_high']==1
    assert result['specificity95']==1 and result['sensitivity95']==1


def test_singleton_healthy_source_is_train_only_not_leaked_into_calibration():
    meta=pd.DataFrame(dict(patient_id=[f'p{i}' for i in range(31)],center=['A']*20+['B']*11,
                           label=[0,1]*10+[0]+[1]*10))
    splits,rare=source_calibration_splits(meta)
    assert rare.tolist()==[20]
    tested=[]
    for train,test in splits:
        assert 20 in train and 20 not in test
        assert meta.iloc[train].query('label==0').center.nunique()==2
        tested.extend(test)
    assert sorted(tested)==[i for i in range(31) if i!=20]
