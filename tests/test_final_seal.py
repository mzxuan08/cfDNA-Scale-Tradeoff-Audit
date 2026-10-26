"""Read-only publication checks; no estimator or classifier is fitted."""
from pathlib import Path
import ast,hashlib,json,re
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1];FINAL=ROOT/'results/step10_final'
pytestmark=pytest.mark.skipif(not FINAL.exists(),reason='Local final publication workspace required')

def test_stage_counts_reconcile_with_raw_selected_metadata():
    source=pd.read_parquet(ROOT/'data/metadata/cohort.parquet')
    scores=pd.read_parquet(ROOT/'results/tables/Step7_OOF_Predictions.parquet')
    raw=source[source.sample_id.isin(scores.sample_id.unique())&source.center.eq('Cristiano et al, 2019')&source.label.eq(1)]
    assert len(raw)==raw.patient_id.nunique()==218
    actual=raw.stage.astype('string').fillna('Unknown/Missing').value_counts().to_dict()
    audit=pd.read_csv(FINAL/'tables/Cristiano_Stage_Count_Audit.csv')
    for row in audit.itertuples():assert row.cancer_patients==actual.get(row.source_stage,0)
    assert audit.cancer_patients.sum()==218
    assert actual['III']==34 and actual['0']==2 and actual['X']==1

def test_permanova_qualifier_and_exploratory_scope_are_in_the_manuscript():
    text=(ROOT/'results/MANUSCRIPT_FULL_TEXT.md').read_text(encoding='utf-8')
    assert 'marginal' not in text.lower()
    first=text.split('### 3.1.')[1].split('### 3.2.')[0]
    quote='These fractions describe study-associated multivariate variation and should not be interpreted as pure centroid shifts because dispersion heterogeneity was also detected.'
    assert quote in first and re.search(r'0\.008.*?\)\. '+re.escape(quote),first)
    for section in ['## 2. Materials and methods','## 3. Results','## 4. Discussion']:
        value=text.split(section)[1].split('\n## ')[0]
        assert 'post-hoc exploratory sensitivity analyses' in value.lower()
        assert 'confirmatory evidence' in value
    assert 'v0.2.0' in text
    state=json.loads((ROOT/'config/release.json').read_text(encoding='utf-8'))
    if state.get('status')=='published':
        assert 'https://doi.org/'+state['zenodo_doi'] in text
        assert 'No Zenodo DOI has yet been assigned to v0.2.0' not in text
    else:assert 'No Zenodo DOI has yet been assigned to v0.2.0' in text

def test_all_previously_frozen_scientific_files_are_unchanged():
    before=json.loads((FINAL/'reports/frozen_inputs_sha256.json').read_text(encoding='utf-8'))
    assert len(before)>1000
    for rel,digest in before.items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel

def test_polish_entrypoint_has_no_training_calls():
    tree=ast.parse((ROOT/'scripts/14_final_polish.py').read_text(encoding='utf-8'))
    calls=[x for x in ast.walk(tree) if isinstance(x,ast.Call)]
    assert not any(isinstance(x.func,ast.Attribute) and x.func.attr in ['fit','fit_transform','partial_fit'] for x in calls)
    state=json.loads((FINAL/'reports/publication_manifest.json').read_text(encoding='utf-8'))
    assert state['classifier_fits']==state['model_selection_runs']==0

def test_only_figure1_has_changed_and_all_numeric_csvs_are_frozen():
    previous=ROOT/'results/step09_major'
    for p in (previous/'figures').glob('*'):
        if p.stem=='Figure1_Batch_Vulnerability_Audit':continue
        assert (FINAL/'figures'/p.name).read_bytes()==p.read_bytes()
    for p in (previous/'tables').glob('*.csv'):
        if p.name=='Table2_Study_Depth_Fragments_Composition.csv':
            a=pd.read_csv(p);b=pd.read_csv(FINAL/'tables'/p.name)
            pd.testing.assert_frame_equal(a.drop(columns='stage_distribution'),b.drop(columns='stage_distribution'))
        else:assert (FINAL/'tables'/p.name).read_bytes()==p.read_bytes(),p

def test_versioned_release_contains_all_revision_outputs():
    release=ROOT/'submission_package/github_zenodo_release'
    if not (release/'results/step10_final').exists():pytest.skip('Release synchronization is the final packaging phase')
    citation=json.loads((release/'CITATION.cff').read_text(encoding='utf-8'))
    assert citation['version']=='0.2.0'
    state=json.loads((ROOT/'config/release.json').read_text(encoding='utf-8'))
    if state.get('status')=='published':assert citation['doi']==state['zenodo_doi']
    else:assert 'doi' not in citation
    for p in (FINAL/'tables').glob('*.csv'):assert (release/'results/step10_final/tables'/p.name).read_bytes()==p.read_bytes()
    for p in (FINAL/'figures').glob('*'):
        if p.suffix in ['.png','.pdf']:assert (release/'results/step10_final/figures'/p.name).read_bytes()==p.read_bytes()
    for name in ['11_major_revision_analysis.py','12_build_major_revision_package.py','13_finalize_major_revision_package.py','14_final_polish.py','15_seal_v020_package.py']:
        assert (release/'scripts'/name).read_bytes()==(ROOT/'scripts'/name).read_bytes()
