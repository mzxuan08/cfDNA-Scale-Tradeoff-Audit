"""Public-source cohort recovery, source-specific labels, and deterministic screening."""
from __future__ import annotations

import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd

from .common import sha256, write_json, write_parquet
from .stream_extractor import small_download, RemoteResourceError, remote_operation

CANCERS = {'Bile Duct', 'Breast', 'Colorectal', 'Duodenal', 'Gastric', 'Lung', 'Ovarian',
           'Pancreatic', 'Hepatocellular', 'Bladder', 'Esophageal', 'Head And Neck', 'Melanoma',
           'Prostate', 'Renal', 'Testicular', 'Uterine'}
KNOWN_QC_FAILURES = {'EE87920', 'EE87921', 'EE86749', 'EE86532'}


def evidence_download(url, path, network, limit=75 * 1024 ** 2):
    path = Path(path)
    stamp = path.with_suffix(path.suffix + '.json')
    if path.exists() and stamp.exists():
        old = json.loads(stamp.read_text(encoding='utf-8'))
        if old.get('url') == url and old.get('sha256') == sha256(path):
            return path
    small_download(url, path, network, limit=limit)
    return path


def published_metadata(cfg):
    """Recover a published FinaleDB snapshot (2362 runs), not a live full catalog.

    UNITE data file S1 sheet s2 is a source-curated clinical metadata table. Its
    FINALEDB subset contains actual EE run IDs, disease, study, patient/specimen,
    plasma status and measured sequencing depth. No values are inferred from IDs.
    """
    folder = Path(cfg['_root']) / 'data/metadata'
    archive = evidence_download(cfg['metadata']['published_supplement_url'],
                               folder / 'PMC13353424_supplement.zip', cfg['network'])
    with zipfile.ZipFile(archive) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read('sciadv.ady9432_data_files_s1_and_s2.zip'))) as inner:
            workbook = inner.read('ady9432_data_file_s1.xlsx')
    frame = pd.read_excel(io.BytesIO(workbook), sheet_name='s2', engine='openpyxl')
    for column in ('type', 'study', 'sample_type', 'data_source', 'bam_id'):
        frame[column] = frame[column].astype(str).str.strip()
    frame = frame[frame.data_source.eq('FINALEDB')].copy()
    if frame.empty or not frame.bam_id.str.fullmatch(r'EE\d+').all() or frame.bam_id.duplicated().any():
        raise ValueError('Published metadata has unexpected or duplicate FinaleDB identifiers')
    write_parquet(frame, folder / 'finaledb_published_metadata.parquet')
    export_published_catalog(frame, cfg)
    write_json({'url': cfg['metadata']['published_supplement_url'], 'archive_sha256': sha256(archive),
                'nested_member': 'sciadv.ady9432_data_files_s1_and_s2.zip/ady9432_data_file_s1.xlsx',
                'sheet': 's2', 'rows': len(frame), 'scope': 'published_snapshot_not_complete_live_catalog'},
               folder / 'finaledb_published_metadata.json')
    return frame


def export_published_catalog(raw, cfg):
    """Enrich the entire published snapshot without inventing availability/labels."""
    catalog = raw.copy()
    catalog['sample_id'] = catalog.bam_id
    catalog['disease_state'] = 'other_or_ineligible_condition'
    catalog.loc[catalog.type.isin(CANCERS), 'disease_state'] = 'cancer'
    catalog.loc[catalog.type.eq('Healthy'), 'disease_state'] = 'healthy_control'
    catalog['cancer_type'] = catalog.type
    catalog['publication'] = catalog.study
    catalog['GSE_id'] = pd.Series(pd.NA, index=catalog.index, dtype='string')
    catalog['genome'] = 'hg19'
    base = cfg['metadata']['public_s3_base'].rstrip('/')
    catalog['fragment_url'] = [f'{base}/entries/{i}/hg19/{i}.hg19.frag.tsv.bgz' for i in catalog.sample_id]
    catalog['object_availability'] = 'candidate_url_not_yet_stream_verified'
    catalog['source'] = cfg['metadata']['published_metadata_citation'] + '; data file S1/s2'
    catalog['catalog_scope'] = '2362_run_published_clinical_snapshot_not_live_complete_database'
    write_parquet(catalog, Path(cfg['_root']) / 'data/metadata/finaledb_catalog_recovered.parquet')
    return catalog


def inspect_live_sources(cfg):
    """Probe with required backoff; preserve exact HTTP/TLS evidence on failure."""
    attempts = []
    for url in cfg['metadata']['api_urls']:
        try:
            def consume(response):
                payload = response.json()
                return payload
            payload = remote_operation(url, cfg['network'], consume)
            # API shape must be interpreted explicitly rather than silently guessing.
            attempts.append({'url': url, 'status': 'reachable', 'response_type': type(payload).__name__})
            write_json(payload, Path(cfg['_root']) / 'data/metadata/live_api_first_page.json')
        except (RemoteResourceError, ValueError) as exc:
            attempts.append({'url': url, 'status': 'unavailable', 'error': str(exc)})
    write_json(attempts, Path(cfg['_root']) / 'data/metadata/finaledb_live_access.json')
    return attempts


def stratified_cap(frame, limit, seed):
    """Proportional label/cancer strata; keep at least one member of each stratum."""
    if len(frame) <= limit:
        return frame.sort_values('sample_id').reset_index(drop=True)
    groups = list(frame.groupby(['center', 'label', 'cancer_type'], sort=True, dropna=False))
    if len(groups) > limit:
        raise ValueError('Sample cap cannot represent every clinical stratum')
    weights = np.array([len(g) for _, g in groups], dtype=float)
    ideal = weights / weights.sum() * limit
    counts = np.maximum(1, np.floor(ideal).astype(int))
    while counts.sum() > limit:
        eligible = np.flatnonzero(counts > 1)
        j = eligible[np.argmax(counts[eligible] - ideal[eligible])]
        counts[j] -= 1
    while counts.sum() < limit:
        eligible = np.flatnonzero(counts < weights)
        j = eligible[np.argmax(ideal[eligible] - counts[eligible])]
        counts[j] += 1
    return pd.concat([g.sample(n=int(n), random_state=seed) for (_, g), n in zip(groups, counts)],
                     ignore_index=True).sort_values('sample_id').reset_index(drop=True)


def select_cohort(raw, cfg):
    c = cfg['cohort']
    raw = raw.copy()
    clinical = raw.type.isin(CANCERS | {c['healthy_label']})
    checks = pd.DataFrame({'not_plasma': ~raw.sample_type.str.casefold().eq('plasma'),
                           'not_cancer_or_healthy': ~clinical,
                           'outside_depth_range': ~raw.seq_depth.between(c['min_depth'], c['max_depth']),
                           'published_qc_failure': raw.bam_id.isin(KNOWN_QC_FAILURES)})
    excluded = raw[['bam_id', 'patient_id', 'study', 'type']].copy()
    excluded['reason'] = checks.apply(lambda row: ';'.join(row.index[row.to_numpy()]), axis=1)
    kept = raw.loc[~checks.any(axis=1)].copy()
    kept['sample_id'] = kept.bam_id
    kept['center'] = kept.study
    kept['label'] = kept.type.ne(c['healthy_label']).astype('int8')
    kept['cancer_type'] = kept.type
    if c['one_measurement_per_specimen']:
        kept = kept.sort_values(['seq_depth', 'sample_id']).drop_duplicates(['study', 'specimen_id'], keep='first')
    eligible = []
    for center, group in kept.groupby('center'):
        if len(group) < c['min_per_center'] or (c['require_both_labels_per_study'] and group.label.nunique() < 2):
            excluded.loc[excluded.bam_id.isin(group.sample_id), 'reason'] = 'study_insufficient_samples_or_no_healthy_controls'
        else:
            eligible.append(stratified_cap(group, c['max_per_center'], cfg['seed']))
    if not eligible:
        raise ValueError('No source studies meet the cohort criteria')
    cohort = stratified_cap(pd.concat(eligible, ignore_index=True), c['max_samples'], cfg['seed'])
    validate_cohort_size(cohort, cfg)
    base = cfg['metadata']['public_s3_base'].rstrip('/')
    cohort['fragment_url'] = [f'{base}/entries/{i}/hg19/{i}.hg19.frag.tsv.bgz' for i in cohort.sample_id]
    citation = cfg['metadata']['published_metadata_citation']
    cohort['center_source'] = citation + '; data file S1/s2 study field; STUDY PROXY, not verified physical center'
    cohort['label_source'] = citation + '; data file S1/s2 type field, Healthy exact match'
    cohort['study_id'], cohort['publication'] = cohort.study, cohort.study
    cohort['genome'], cohort['genome_source'] = 'hg19', 'Explicit FinaleDB hg19 object path; coordinates checked against hg19 reference'
    cohort['verified'] = True
    cohort['sex'] = cohort.gender
    cohort['sequencing_depth'] = cohort.seq_depth.astype('float32')
    excluded.loc[excluded.bam_id.isin(cohort.sample_id), 'reason'] = 'selected'
    excluded.loc[excluded.reason.eq(''), 'reason'] = 'deduplication_or_seeded_stratified_cap'
    return cohort.reset_index(drop=True), excluded


def validate_cohort_size(cohort, cfg):
    c = cfg['cohort']
    if not c['min_samples'] <= len(cohort) <= c['max_samples']:
        raise ValueError(f'Cohort size {len(cohort)} outside {c["min_samples"]}--{c["max_samples"]}')
    counts = cohort.groupby('center').size()
    if len(counts) < c['min_centers'] or counts.lt(c['min_per_center']).any():
        raise ValueError('Not enough independent studies or fewer than 20 samples in one study')


def import_published_fallback(cfg, cause):
    """Acquire a real author-published DELFI table; NEVER fabricate absent fields.

    It is a separate alternative cohort, not replacement measurements of failed EE
    patients. LUCAS is a diagnostic risk cohort, not an interchangeable healthy cohort.
    Without 321-bin histograms, IQR and globally filtered depth, canonical admission
    is refused; the downloaded counts still remain reusable secondary evidence.
    """
    import pyreadr
    folder = Path(cfg['_root']) / 'data/metadata/published_fallback'
    folder.mkdir(parents=True, exist_ok=True)
    path = evidence_download(cfg['extraction']['fallback_lucas_url'], folder / 'lucas_5mb.rda', cfg['network'])
    frame = pyreadr.read_r(path)['bins5mb']
    write_parquet(frame, Path(cfg['paths']['processed']) / 'published_lucas_original_bins.parquet')
    record = {'url': cfg['extraction']['fallback_lucas_url'], 'sha256': sha256(path),
              'source': 'Mathios-related LUCAS 5 Mb counts republished in authors artemis2024/rlucas',
              'rows': len(frame), 'samples': int(frame.id.nunique()), 'trigger': cause,
              'canonical_eligible': False,
              'reason': 'Different cohort and interval/feature definitions; missing global MAPQ30 80--400 counts, IQR and 321-length vectors; controls are not all healthy',
              'available_columns': list(frame.columns)}
    write_json(record, Path(cfg['paths']['results']) / 'reports/step_02_fallback.json')
    return record
