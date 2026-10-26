"""Opt-in validation against the authors' public chr22 tutorial fragment file.

No raw intervals or clinical labels are saved. This validates IO/aggregation only.
"""
import gzip
import os
from pathlib import Path
import time

import numpy as np
import pytest

from src.common import ROOT, load_config, write_json, fingerprint
from src.stream_extractor import remote_operation, HashReader, bounded_lines, aggregate_fragments, Blacklist, make_bins


@pytest.mark.skipif(os.environ.get('RUN_PUBLIC_STREAM') != '1', reason='opt-in real public network transfer (~52 MB)')
def test_public_chr22_stream_without_raw_storage():
    cfg = load_config()
    ref = Path(cfg['paths']['reference'])
    bins = make_bins(ref / 'hg19.chrom.sizes', cfg['features']['bin_size'], ['chr22'])
    url = 'https://zenodo.org/records/20724659/files/BH01.hg19.mdups.22.frag.gz?download=1'
    report = {'url': url, 'kind': 'real_public_tutorial_data_IO_validation_no_clinical_labels',
              'raw_files_stored': 0, 'feature_parameters': cfg['features']}
    start = time.perf_counter()
    try:
        def consume(response):
            reader = HashReader(response.raw)
            with gzip.GzipFile(fileobj=reader) as handle:
                values, qc = aggregate_fragments(bounded_lines(handle), bins, cfg['features'], Blacklist(ref / 'hg19.blacklist.bed.gz'))
            return values, {**qc, 'compressed_bytes': reader.count, 'compressed_sha256': reader.digest.hexdigest()}
        values, qc = remote_operation(url, cfg['network'], consume)
        report.update(status='complete', qc=qc, feature_count=len(values), finite_features=int(sum(np.isfinite(v) for v in values.values())))
        assert qc['accepted'] >= 1000
        assert sum(v for k, v in values.items() if k.startswith('length__')) == pytest.approx(1, abs=1e-6)
    except Exception as exc:
        report.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        report['runtime_seconds'] = time.perf_counter() - start
        write_json(report, ROOT / 'results/reports/real_stream_validation.json')
