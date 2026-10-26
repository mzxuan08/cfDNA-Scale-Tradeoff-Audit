"""Synthetic inputs exist ONLY in tests and temporary test directories."""
import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.common import ROOT, load_config, fingerprint, sha256, write_parquet, write_json, analysis_signature, diagnostic_report
from src.evaluation import run_benchmarks, summarize_predictions
from src.variance_audit import audit
from src.plotting import fig1, fig2, fig3, fig4


def test_small_full_statistical_pipeline_and_assets(tmp_path, cohort):
    cfg = yaml.safe_load((ROOT / 'config/config.yaml').read_text())
    (tmp_path / 'config').mkdir()
    cfg['evaluation'].update(models=['elasticnet'], methods=['raw', 'lsn', 'icc', 'lsn_icc', 'combat'],
                             random_folds=3, bootstrap_replicates=30)
    cfg['audit'].update(permutations=19, pca_components=3, umap=False, covariates=[])
    (tmp_path / 'config/config.yaml').write_text(yaml.safe_dump(cfg))
    cfg = load_config(tmp_path / 'config/config.yaml')
    rng = np.random.default_rng(2026)
    names = ['ratio__chr1:0-5', 'ratio__chr1:5-10', 'ratio__chr2:0-5', 'length__80']
    X = pd.DataFrame((rng.lognormal(size=(len(cohort), 4)) + np.repeat([0, 1, 2], 12)[:, None] * 4).astype('float32'), columns=names)
    predictions, failures = run_benchmarks(cohort, X, cfg, 'test-only-signature')
    assert any('ICC requires' in f['reason'] for f in failures)
    assert any('ComBat' in f['reason'] for f in failures)
    summary, centers, pairs, mitigation = summarize_predictions(predictions, cohort, cfg)
    assert len(summary[summary.status.eq('complete')]) == 5
    assert set(pairs.method) == {'raw', 'lsn'}
    assert set(mitigation.method) == {'lsn', 'combat'}
    assert len(centers) == 15
    embedding, lmm, details = audit(cohort, X, cfg)
    assert 0 <= details['permanova']['p'] <= 1
    assert len(lmm) == 4
    fig1(embedding, cfg)
    fig2(summary, cfg)
    fig3(summary, predictions, cfg)
    valid = lmm[lmm.status.eq('converged')]
    if len(valid) and 'disease_beta' in valid:
        fig4(valid, cfg)
    for name in ('Fig1_batch_structure', 'Fig2_generalization', 'Fig3_mitigation'):
        pdf = Path(cfg['paths']['results']) / 'figures' / (name + '.pdf')
        png = pdf.with_suffix('.png')
        assert pdf.read_bytes().startswith(b'%PDF') and png.stat().st_size > 1000
    # Exercise the real saved-input contract and manuscript CLI stages. All
    # generated subjects and outputs stay in pytest's temporary directory.
    manifest = cohort.copy()
    manifest['center_source'] = manifest['label_source'] = 'test fixture only'
    manifest['study_id'] = manifest.center
    manifest['genome'], manifest['verified'] = 'hg19', True
    manifest['fragment_url'] = 'http://test.invalid/test-only.frag.gz'
    write_parquet(manifest, cfg['paths']['metadata'])
    ref = Path(cfg['paths']['reference'])
    (ref / 'hg19.chrom.sizes').write_text('chr1\t10000000\nchr2\t5000000\n')
    (ref / 'hg19.blacklist.bed.gz').write_bytes(b'test reference identity only')
    features = X.copy()
    features.insert(0, 'sample_id', cohort.sample_id)
    fpath = Path(cfg['paths']['processed']) / 'features.parquet'
    write_parquet(features, fpath)
    write_json({'metadata_sha256': sha256(cfg['paths']['metadata']),
                'feature_config': fingerprint(cfg['features']), 'features_sha256': sha256(fpath),
                'reference_sha256': {'sizes': sha256(ref / 'hg19.chrom.sizes'),
                                     'blacklist': sha256(ref / 'hg19.blacklist.bed.gz')}}, fpath.with_suffix('.json'))
    import sys
    sys.path.insert(0, str(ROOT / 'scripts'))
    old_argv = sys.argv[:]
    try:
        for step in range(3, 7):
            script = next((ROOT / 'scripts').glob(f'{step:02d}_*.py'))
            spec = importlib.util.spec_from_file_location(f'test_cli_step_{step}', script)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            sys.argv = [str(script), '--config', cfg['_config_path']]
            module.main()
            report = Path(cfg['paths']['results']) / f'reports/step_{step:02d}_report.md'
            assert report.exists()
        asset = Path(cfg['paths']['results']) / 'asset_manifest.json'
        assert len(json.loads(asset.read_text())['files']) > 10
    finally:
        sys.argv = old_argv


def test_failure_report_is_written_even_without_cohort(tmp_path):
    path = tmp_path / 'config/config.yaml'
    path.parent.mkdir()
    cfg = yaml.safe_load((ROOT / 'config/config.yaml').read_text())
    path.write_text(yaml.safe_dump(cfg))
    cfg = load_config(path)
    with pytest.raises(RuntimeError):
        with diagnostic_report(cfg, 2):
            raise RuntimeError('test diagnostic')
    output = (tmp_path / 'results/reports/step_02_report.md').read_text()
    assert 'test diagnostic' in output and '0 samples analyzed' in output and 'runtime_seconds' in output
