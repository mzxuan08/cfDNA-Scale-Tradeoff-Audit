import pandas as pd
from src.metadata_scope_audit import study_counts


def test_healthy_threshold_counts_unique_known_people_not_runs():
    f=pd.DataFrame({'study':['A']*12,'patient_id':[f'h{i}' for i in range(9)]+['h0',None,'c1'],
                    'type':['Healthy']*11+['Breast']})
    result=study_counts(f).iloc[0]
    assert result.healthy_runs==11 and result.healthy_patients==9
    assert not result.meets_cancer_and_10_healthy
    f.loc[10,'patient_id']='h9'
    assert study_counts(f).iloc[0].meets_cancer_and_10_healthy


def test_disease_controls_are_not_relabelled_as_healthy():
    f=pd.DataFrame({'study':['A']*21,'patient_id':[f'p{i}' for i in range(21)],
                    'type':['Cirrhosis']*10+['Liver Transplant']*10+['Colorectal']})
    result=study_counts(f).iloc[0]
    assert result.healthy_patients==0 and result.other_condition_runs==20
    assert not result.meets_cancer_and_10_healthy


def test_controls_are_not_borrowed_from_another_source_and_blank_ids_are_unknown():
    f=pd.DataFrame({'study':['A']+['B']*11,'patient_id':['c1']+[f'h{i}' for i in range(10)]+['  '],
                    'type':['Breast']+['Healthy']*11})
    result=study_counts(f).set_index('study')
    assert not result.meets_cancer_and_10_healthy.any()
    assert result.loc['B','healthy_patients']==10 and result.loc['B','unresolved_patient_id_runs']==1
