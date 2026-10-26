"""Shared contracts for the fully extracted clinical study (no simulated inputs)."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .common import sha256, fingerprint, write_json, write_parquet, read_metadata
from .reference import load_reference_windows

PIPELINES = ['Raw_Ratio', 'Global_Zscore', 'InFold_ComBat', 'Proposed_LSN_ICC']
MODELS = ['elasticnet', 'rf', 'xgboost', 'lightgbm', 'mlp', 'dann']


def clinical_available(cfg):
    path = Path(cfg['paths']['processed']) / cfg.get('extraction', {}).get('final_file', 'cfdna_multicenter_features.parquet')
    return path.exists()


def load_clinical(cfg, patient_level=True):
    path = Path(cfg['paths']['processed']) / cfg['extraction']['final_file']
    stamp = json.loads(path.with_suffix('.json').read_text())
    if sha256(path) != stamp['features_sha256']:
        raise ValueError('Clinical feature matrix checksum mismatch')
    if sha256(cfg['paths']['metadata']) != stamp['metadata_sha256']:
        raise ValueError('Clinical metadata changed after extraction')
    bins = load_reference_windows(cfg)
    if stamp['window_sha256'] != sha256(Path(cfg['paths']['reference']) / cfg['reference']['window_table']):
        raise ValueError('Clinical windows checksum mismatch')
    meta = read_metadata(cfg,allow_unresolved_patients=True)
    unresolved=meta.patient_id.isna()
    if unresolved.any():
        excluded=meta.loc[unresolved,['sample_id','specimen_id','center']].assign(reason='source_patient_identity_missing_no_alias_inferred')
        export_table(excluded,cfg,'Analysis_Identity_Exclusions')
        meta=meta.loc[~unresolved].copy()
    features = pd.read_parquet(path)
    if not set(meta.sample_id)<=set(features.sample_id):
        raise ValueError('Clinical sample identities differ')
    joined = meta.merge(features, on='sample_id', validate='one_to_one')
    # Prespecified earliest specimen, then shallowest depth, then run ID;
    # no feature/outcome performance influences this independent-subject cohort.
    if patient_level:
        joined = joined.sort_values(['patient_id', 'timepoint', 'seq_depth', 'sample_id'],
                                    kind='stable', na_position='last').drop_duplicates('patient_id')
    joined = joined.reset_index(drop=True)
    short = joined[[f'short_count__{b}' for b in bins.bin_id]].to_numpy(dtype=np.float32)
    long = joined[[f'long_count__{b}' for b in bins.bin_id]].to_numpy(dtype=np.float32)
    ratio = short / (long + np.float32(1))
    lengths = joined[[f'length__{n}' for n in range(80, 401)]].to_numpy(dtype=np.float32)
    selected_meta = joined[list(meta.columns) + ['accepted_fragments']].copy()
    selected_meta['center_id'] = selected_meta.center
    selected_meta['disease_status'] = np.where(selected_meta.label.eq(1), 'Cancer', 'Healthy')
    return selected_meta, ratio, lengths, bins


def study_signature(cfg):
    return fingerprint({'matrix': sha256(Path(cfg['paths']['processed']) / cfg['extraction']['final_file']),
        'metadata': sha256(cfg['paths']['metadata']), 'bins': sha256(Path(cfg['paths']['reference']) / cfg['reference']['window_table']),
        'study_settings': cfg.get('study', {}), 'evaluation': cfg['evaluation'], 'version': 'clinical_v3'})


def export_table(frame, cfg, name):
    folder = Path(cfg['paths']['results']) / 'tables'
    folder.mkdir(parents=True, exist_ok=True)
    frame.to_csv(folder / f'{name}.csv', index=False)
    write_parquet(frame, folder / f'{name}.parquet')


def record_outputs(cfg, step, paths):
    write_json({'signature': study_signature(cfg), 'files': {str(Path(p).resolve()): sha256(p) for p in paths}},
               Path(cfg['paths']['results']) / 'reports' / f'{step}_provenance.json')
