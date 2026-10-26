"""Offline census of source eligibility; no live-catalog or object-access claims."""
import io
import json
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
from .cohort import CANCERS,KNOWN_QC_FAILURES
from .common import sha256,write_json
from .study import export_table


def study_counts(frame):
    """Count known unique people, never equate runs with healthy participants."""
    frame=frame.copy()
    frame['patient_id']=frame.patient_id.astype('string').str.strip().replace('',pd.NA)
    rows=[]
    for study,group in frame.groupby('study',sort=True):
        healthy=group.type.eq('Healthy');cancer=group.type.isin(CANCERS)
        h=int(group.loc[healthy,'patient_id'].nunique());c=int(group.loc[cancer,'patient_id'].nunique())
        rows.append(dict(study=study,runs=len(group),known_patients=int(group.patient_id.nunique()),
            healthy_runs=int(healthy.sum()),healthy_patients=h,cancer_runs=int(cancer.sum()),cancer_patients=c,
            other_condition_runs=int((~healthy&~cancer).sum()),unresolved_patient_id_runs=int(group.patient_id.isna().sum()),
            meets_cancer_and_10_healthy=bool(h>=10 and c>=1)))
    return pd.DataFrame(rows)


def publication_inventory(path):
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    with zipfile.ZipFile(path) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read('btaa999_supplementary_data.zip'))) as inner:
            with zipfile.ZipFile(io.BytesIO(inner.read('FinaleDB.supplementary.final.docx'))) as doc:
                xml=ET.fromstring(doc.read('word/document.xml'))
    table=xml.findall('.//w:tbl',ns)[0];rows=[]
    for row in table.findall('w:tr',ns)[1:]:
        cells=[' '.join(t.text or '' for t in cell.findall('.//w:t',ns)) for cell in row.findall('w:tc',ns)]
        if len(cells)!=2 or cells[0].strip()=='Total':continue
        publication,count=cells
        surname=publication.split(',')[0].strip()
        year=next((str(year) for year in range(2015,2020) if str(year) in publication),None)
        if year is None:raise ValueError('Unexpected publication year in archived source table')
        rows.append(dict(study=f'{surname} et al, {year}',publication=publication,publication_records=int(count.strip())))
    return pd.DataFrame(rows)


def audit_local_finaledb(cfg):
    root=Path(cfg['_root']);folder=root/'data/metadata'
    path=folder/'finaledb_published_metadata.parquet';frame=pd.read_parquet(path)
    stamp=json.loads((folder/'finaledb_published_metadata.json').read_text())
    if len(frame)!=stamp['rows'] or not frame.data_source.eq('FINALEDB').all():raise ValueError('Unexpected metadata snapshot')
    if sha256(folder/'PMC13353424_supplement.zip')!=stamp['archive_sha256']:raise ValueError('Clinical source archive checksum')
    plasma=frame[frame.sample_type.str.casefold().eq('plasma')].copy()
    broad=study_counts(plasma).assign(screen='plasma_snapshot_no_depth_restriction')
    depth=plasma[plasma.seq_depth.between(cfg['cohort']['min_depth'],cfg['cohort']['max_depth'])&~plasma.bam_id.isin(KNOWN_QC_FAILURES)]
    locked=study_counts(depth).assign(screen='locked_depth_and_published_QC')
    census=pd.concat([broad,locked],ignore_index=True)
    current=set(pd.read_parquet(cfg['paths']['metadata']).center)
    census['already_in_analysis']=census.study.isin(current)
    census['additional_eligible']=census.meets_cancer_and_10_healthy&~census.already_in_analysis
    census['assay_scope']='published human plasma cfDNA WGS clinical snapshot; no new object availability check'
    export_table(census,cfg,'Final_Alignment_FinaleDB_Source_Census')
    publication=publication_inventory(folder/'PMC8388032_supplement.zip')
    publication=publication.merge(broad.drop(columns='screen'),on='study',how='left',validate='one_to_one')
    reasons=[]
    for row in publication.itertuples():
        if pd.notna(row.runs):
            if row.meets_cancer_and_10_healthy:reason='metadata supports cancer plus at least 10 unique healthy participants'
            elif row.healthy_patients==0:reason='no participants explicitly labelled Healthy; other conditions are not substituted'
            else:reason='fewer than 10 unique healthy participants in the full snapshot'
        elif 'Urinary' in row.publication:reason='reported urinary study; fails plasma criterion; no run-level healthy/cancer counts asserted'
        elif 'prenatal' in row.publication:reason='reported prenatal single-gene study; no qualifying cancer/control cohort documented locally; run-level counts unavailable'
        else:reason='run-level metadata unavailable; eligibility not verified'
        reasons.append(reason)
    publication['eligibility_evidence']=reasons
    export_table(publication,cfg,'Final_Alignment_Original_Publication_Inventory')
    summary=dict(scope='local_historical_snapshot_not_complete_current_FinaleDB_or_all_public_datasets',
        snapshot_runs=len(frame),plasma_runs=len(plasma),run_level_studies=int(broad.study.nunique()),
        archived_publication_studies=len(publication),archived_publication_records=int(publication.publication_records.sum()),
        qualifying_studies=broad.loc[broad.meets_cancer_and_10_healthy,'study'].tolist(),
        qualifying_study_count=int(broad.meets_cancer_and_10_healthy.sum()),
        additional_qualifying_studies=broad.loc[broad.meets_cancer_and_10_healthy&~broad.study.isin(current),'study'].tolist(),
        additional_qualifying_count=int(census[census.screen.eq('plasma_snapshot_no_depth_restriction')].additional_eligible.sum()),
        depth_interval=[cfg['cohort']['min_depth'],cfg['cohort']['max_depth']],
        healthy_definition='type == Healthy; unique nonmissing source patient_id within study; no cross-study controls borrowed',
        unknown_counts_policy='missing run-level records remain unavailable, not zero',
        object_availability='not rechecked; no new clinical downloads or feature extraction',
        source_files={str(p):sha256(p) for p in [path,folder/'finaledb_published_metadata.json',folder/'PMC13353424_supplement.zip',folder/'PMC8388032_supplement.zip']})
    write_json(summary,Path(cfg['paths']['results'])/'reports/final_alignment_metadata_audit.json')
    return census,publication,summary
