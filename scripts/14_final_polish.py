"""Prepare v0.2.0 from frozen outputs and audit source stage labels.

This command reads existing patient metadata, scores and plotted coordinates.
It does not instantiate, fit or select a classifier.
"""
from pathlib import Path
from types import SimpleNamespace
import copy,json,shutil
import pandas as pd
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.common import load_config,read_metadata,sha256,write_json

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'submission_package'
BASE=ROOT/'results/step09_major';NEW=ROOT/'results/step10_final'
ARCHIVE=ROOT/'results/archive_before_step10'

def snapshot():
    if (ARCHIVE/'manifest.json').exists():return
    paths=list(OUT.glob('*.docx'))+list(OUT.glob('*.pdf'))+list(OUT.glob('*.md'))
    paths+=[ROOT/'results/MANUSCRIPT_FULL_TEXT.md',ROOT/'results/reports/MANUSCRIPT_MASTER_DRAFT.md']
    manifest={}
    for p in paths:
        target=ARCHIVE/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,target);manifest[str(p.relative_to(ROOT))]=sha256(target)
    # The existing distribution retains the exact pre-polish source versions.
    release=OUT/'github_zenodo_release'
    for folder in ['src','scripts','config']:
        for p in (release/folder).rglob('*'):
            if p.is_file() and p.suffix in ['.py','.json','.yaml','.yml']:
                target=ARCHIVE/p.relative_to(release);target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(p,target);manifest[str(p.relative_to(release))]=sha256(target)
    write_json(manifest,ARCHIVE/'manifest.json')

def prepare():
    snapshot()
    for sub in ['tables','figures','reports']:(NEW/sub).mkdir(parents=True,exist_ok=True)
    protected_file=NEW/'reports/frozen_inputs_sha256.json'
    if not protected_file.exists():
        paths=list((ROOT/'data').rglob('*'))+list((ROOT/'results/tables').glob('*'))
        paths+=list((BASE/'tables').glob('*'))+list((BASE/'figures').glob('*'))
        paths+=[ROOT/'src/models.py',ROOT/'src/debiasing.py',ROOT/'src/evaluation.py',ROOT/'config/config.yaml']
        write_json({str(p.relative_to(ROOT)):sha256(p) for p in paths if p.is_file()},protected_file)
    for p,h in json.loads(protected_file.read_text(encoding='utf-8')).items():assert sha256(ROOT/p)==h,p
    for sub in ['tables','figures']:
        for p in (BASE/sub).iterdir():
            if p.is_file():shutil.copy2(p,NEW/sub/p.name)
    shutil.copy2(BASE/'reports/figure_typography.json',NEW/'reports/figure_typography.json')
    cfg=load_config();metadata=read_metadata(cfg,allow_unresolved_patients=True)
    scores=pd.read_parquet(ROOT/'results/tables/Step7_OOF_Predictions.parquet')
    ids=scores[['sample_id','patient_id']].drop_duplicates()
    assert len(ids)==ids.patient_id.nunique()==585
    selected=metadata[metadata.sample_id.isin(ids.sample_id)].copy()
    assert len(selected)==585 and selected.patient_id.nunique()==585
    cancer=selected[selected.center.eq('Cristiano et al, 2019')&selected.label.eq(1)].copy()
    assert len(cancer)==218
    cancer['raw_stage']=cancer.stage.astype('string').fillna('Unknown/Missing').str.strip().replace('', 'Unknown/Missing')
    counts=cancer.raw_stage.value_counts().to_dict()
    assert set(counts)<=set(['0','I','II','III','IV','X','Unknown/Missing'])
    stages=['0','I','II','III','IV','X','Unknown/Missing']
    audit=pd.DataFrame([{'study':'Cristiano et al, 2019','source_stage':x,'cancer_patients':int(counts.get(x,0))} for x in stages])
    assert audit.cancer_patients.sum()==218
    audit.to_csv(NEW/'tables/Cristiano_Stage_Count_Audit.csv',index=False)
    cancer[['sample_id','patient_id','raw_stage']].sort_values(['raw_stage','patient_id']).to_csv(NEW/'tables/Cristiano_Stage_Patient_Trace.csv',index=False)
    cohort=pd.read_csv(BASE/'tables/Table2_Study_Depth_Fragments_Composition.csv')
    wording='; '.join(('Unknown/Missing' if x=='Unknown/Missing' else 'Stage '+x)+': '+str(counts.get(x,0)) for x in stages)+'; Total: 218'
    cohort.loc[cohort.study.eq('Cristiano et al, 2019'),'stage_distribution']=wording
    cohort.to_csv(NEW/'tables/Table2_Study_Depth_Fragments_Composition.csv',index=False)
    cohort.to_parquet(NEW/'tables/Table2_Study_Depth_Fragments_Composition.parquet',index=False)
    write_json({'version':'v0.2.0','status':'prepared_locally_not_published','classifier_fits':0,
        'model_selection_runs':0,'clinical_statistics_recomputed':False,'patients':585,
        'stage_counts':counts,'stage_total':218,'source_metadata_sha256':sha256(cfg['paths']['metadata']),
        'selected_scores_sha256':sha256(ROOT/'results/tables/Step7_OOF_Predictions.parquet')},NEW/'reports/publication_manifest.json')
    from src.revision_manuscript import compile_revision
    compile_revision(ROOT)
    # Read-only feature loading for the observed healthy length curves.
    from src.study import load_clinical
    local=copy.deepcopy(cfg);local['paths']['results']=str(NEW)
    meta,ratios,lengths,bins=load_clinical(local,True)
    a=SimpleNamespace(root=ROOT,out=NEW,tables=NEW/'tables',meta=meta,R=ratios,lengths=lengths,bins=bins,
        pred=scores[scores.model.eq('elasticnet')])
    from src.revision_plotting import revision_figures
    revision_figures(a,figure1_only=True)
    for p,h in json.loads(protected_file.read_text(encoding='utf-8')).items():assert sha256(ROOT/p)==h,p
    report='# Step 10 final polish\n\n'
    report+='No model or classifier was trained. All discrimination, calibration, anchor, spatial and dispersion statistics were retained. Only Figure 1 terminology was redrawn from frozen coordinates/tables.\n\n'
    report+='Cristiano stage counts from the 218 selected cancer participants: '+wording+'. Stage III was 34, not 3; no unobserved cases were imputed.\n\n'
    report+='The shared-depth and GI-versus-HCC results are post-hoc exploratory sensitivity analyses from the same public cohorts with altered target size/composition, not independent external validation or confirmatory evidence.\n\n'
    report+='Software v0.2.0 is a complete local release package. No new Zenodo DOI or public release is asserted. The v0.1.0 DOI remains historical.\n'
    (ROOT/'results/reports/step_10_report.md').write_text(report,encoding='utf-8')
    print(audit.to_string(index=False));print('Classifier fits: 0; frozen scientific inputs unchanged.')

if __name__=='__main__':prepare()
