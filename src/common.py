from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import sys
import time
import traceback
from contextlib import contextmanager

import numpy as np
import pandas as pd
import yaml

SEED = 2026
ROOT = Path(__file__).resolve().parents[1]


def seed_everything(seed=SEED, include_torch=True):
    # PYTHONHASHSEED affects a new interpreter only; launcher sets it before import.
    os.environ['PYTHONHASHSEED'] = str(seed)
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    random.seed(seed)
    np.random.seed(seed)
    if not include_torch:
        return
    try:
        import torch
    except ImportError:
        return
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)


def load_config(path=None):
    path = Path(path or ROOT / 'config/config.yaml').resolve()
    cfg = yaml.safe_load(path.read_text(encoding='utf-8'))
    root = path.parent.parent
    cfg['_root'] = str(root)
    cfg['_config_path'] = str(path)
    for key, value in cfg['paths'].items():
        cfg['paths'][key] = str((root / value).resolve())
    for key in ('reference', 'processed', 'results'):
        Path(cfg['paths'][key]).mkdir(parents=True, exist_ok=True)
    # Neural models seed torch immediately before fit; reference/audit scripts
    # must also work when an optional CUDA runtime cannot be imported.
    seed_everything(cfg['seed'], include_torch=False)
    return cfg


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def fingerprint(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


def analysis_signature(cfg):
    root = ROOT
    code = {str(p.relative_to(root)): sha256(p) for folder in ('src', 'scripts')
            for p in sorted((root / folder).glob('*.py'))}
    return fingerprint({'config': cfg, 'code': code,
                        'metadata': sha256(cfg['paths']['metadata']),
                        'features': sha256(Path(cfg['paths']['processed']) / 'features.parquet')})


def write_json(obj, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, default=str, allow_nan=False), encoding='utf-8')
    os.replace(tmp, path)


def write_parquet(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.parquet')
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)


def read_metadata(cfg, allow_unresolved_patients=False):
    path = Path(cfg['paths']['metadata'])
    if not path.exists():
        raise FileNotFoundError(f'Verified cohort manifest required: {path}. See data/metadata/SCHEMA.md')
    df = pd.read_parquet(path)
    required = ['sample_id', 'patient_id', 'center', 'center_source', 'study_id',
                'label', 'label_source', 'genome', 'fragment_url', 'verified']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f'Missing metadata columns: {sorted(missing)}')
    nonmissing = [c for c in required if c!='patient_id' or not allow_unresolved_patients]
    if df.empty or df[nonmissing].isna().any().any():
        raise ValueError('Manifest must contain nonmissing verified records')
    if not df['verified'].eq(True).all() or not df.label.isin([0, 1]).all():
        raise ValueError('Require verified=True and explicit binary clinical labels 0/1')
    for col in required:
        if col not in ('label', 'verified') and df[col].astype(str).str.strip().eq('').any():
            raise ValueError(f'Empty {col}')
    if df.sample_id.duplicated().any():
        raise ValueError('sample_id must identify one selected sequencing run uniquely')
    if df.groupby('patient_id').label.nunique().gt(1).any():
        raise ValueError('Inconsistent patient labels')
    if not df.genome.eq(cfg['reference']['genome']).all():
        raise ValueError('Mixed/mismatched genome builds are forbidden')
    if not df.fragment_url.str.match(r'^https?://').all():
        raise ValueError('fragment_url must be an explicit HTTP(S) URL')
    return df.reset_index(drop=True)


def load_dataset(cfg):
    meta = read_metadata(cfg)
    path = Path(cfg['paths']['processed']) / 'features.parquet'
    features = pd.read_parquet(path)
    provenance = json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
    if provenance['metadata_sha256'] != sha256(cfg['paths']['metadata']):
        raise ValueError('Manifest changed after extraction; rerun step 02')
    if provenance['feature_config'] != fingerprint(cfg['features']):
        raise ValueError('Feature settings changed; rerun step 02')
    ref, genome = Path(cfg['paths']['reference']), cfg['reference']['genome']
    for key, suffix in [('sizes', '.chrom.sizes'), ('blacklist', '.blacklist.bed.gz')]:
        if provenance['reference_sha256'][key] != sha256(ref / (genome + suffix)):
            raise ValueError('Reference changed after extraction; rerun step 02')
    if provenance['features_sha256'] != sha256(path):
        raise ValueError('Feature matrix checksum mismatch')
    if features.sample_id.duplicated().any() or set(features.sample_id) != set(meta.sample_id):
        raise ValueError('Feature/sample identities do not match manifest exactly')
    df = meta.merge(features, on='sample_id', validate='one_to_one', sort=False)
    cols = [c for c in features if c.startswith(('ratio__', 'length__'))]
    if not cols:
        raise ValueError('No fragmentomic features')
    X = df[cols].astype('float32')
    if np.isinf(X.to_numpy()).any():
        raise ValueError('Infinite features')
    return df, X


@contextmanager
def diagnostic_report(cfg, step):
    start = time.perf_counter()
    report = {'step': step, 'status': 'running', 'seed': cfg['seed'],
              'config_sha256': sha256(cfg['_config_path']),
              'metadata_missingness': 'Unavailable until a verified cohort manifest is supplied',
              'feature_missingness': 'Not computed for this step',
              'statistics': 'Not computed for this step; no inferred p values or CIs.'}
    try:
        p = Path(cfg['paths']['metadata'])
        if p.exists():
            m = pd.read_parquet(p)
            report['metadata_missingness'] = m.isna().mean().to_dict()
            if {'center', 'sample_id', 'patient_id', 'label'} <= set(m):
                report['counts_per_center'] = m.groupby('center').agg(
                    samples=('sample_id', 'size'), patients=('patient_id', 'nunique'),
                    cases=('label', 'sum')).to_dict('index')
        else:
            report['counts_per_center'] = 'Unavailable: verified cohort not yet supplied (0 samples analyzed).'
        yield report
        if report['status'] == 'running':
            report['status'] = 'complete'
    except BaseException as exc:
        report['status'] = 'failed'
        report['error'] = f'{type(exc).__name__}: {exc}'
        report['traceback'] = traceback.format_exc()
        raise
    finally:
        report['runtime_seconds'] = round(time.perf_counter() - start, 3)
        report['python'] = sys.version
        report['packages'] = {p: importlib.metadata.version(p) for p in
                              ('numpy', 'pandas', 'scikit-learn', 'scipy', 'requests')}
        dest = Path(cfg['paths']['results']) / 'reports' / f'step_{step:02d}_report.md'
        dest.parent.mkdir(parents=True, exist_ok=True)
        markdown = report.pop('_markdown', '')
        dest.write_text(f'# Step {step:02d} diagnostic report\n\n' + markdown + '\n```json\n' +
                        json.dumps(report, indent=2, default=str, ensure_ascii=False) + '\n```\n', encoding='utf-8')

