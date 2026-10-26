"""Mathematical fixtures remain in tests; production windows use UCSC/ENCODE."""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.common import sha256, write_json
from src.reference import (build_quality_windows, overlap_bp, read_chrominfo, read_cytobands,
                           arm_distribution, fetch_reference_resource)
from src.stream_extractor import Blacklist, RemoteResourceError


def test_overlap_union_and_half_open_coordinates():
    assert overlap_bp(0, 100, [(20, 60), (40, 80)]) == 60
    assert overlap_bp(0, 20, [(20, 30)]) == 0


def test_fixed_5mb_quality_filters_exact_blacklist_threshold_and_arm_distances():
    million = 1000000
    sizes = {'chr1': 62_500_000}
    cyto = pd.DataFrame([('chr1', 10 * million, 15 * million, 'p21', 'stalk')],
                        columns=['chrom', 'start', 'end', 'band', 'stain'])
    masks = Blacklist()
    masks.intervals['chr1'] = np.array([[15_000_000, 15_600_000], [15_400_000, 16_000_000],
                                       [35_000_000, 36_100_000]])
    filters = dict(full_windows_only=True, telomere_buffer_bp=million,
                   max_blacklist_fraction=.2, excluded_cytoband_stains=['gvar', 'stalk'])
    bins, candidates = build_quality_windows(sizes, cyto, {'chr1': (25 * million, 30 * million)}, masks, 5 * million, filters)
    assert len(candidates) == 13 and len(bins) == 8
    assert (bins.end - bins.start).eq(5 * million).all()
    assert bins.start.min() >= million
    assert not bins.start.isin([0, 10 * million, 25 * million, 35 * million, 60 * million]).any()
    # Overlapping blacklist intervals cover UNION=1 Mb, exactly 20%: retain it.
    boundary = bins[bins.start.eq(15 * million)].iloc[0]
    assert boundary.blacklist_overlap_bp == million
    assert boundary.blacklist_fraction == pytest.approx(.2)
    p = bins[bins.start.eq(20 * million)].iloc[0]
    q = bins[bins.start.eq(30 * million)].iloc[0]
    assert p.arm == 'p' and q.arm == 'q'
    assert p.distance_to_centromere_bp == q.distance_to_centromere_bp == 0
    assert p.distance_to_telomere_bp == 20 * million
    assert q.distance_to_telomere_bp == 27_500_000
    dist = arm_distribution(bins, ['chr1', 'chr2'])
    assert len(dist) == 4 and dist[dist.chrom.eq('chr2')].windows.eq(0).all()
    assert dist.windows.sum() == len(bins)


def test_ucsc_chrominfo_cytobands_parse_and_centromere_union(tmp_path):
    ci, cb = tmp_path / 'chromInfo.txt.gz', tmp_path / 'cytoBand.txt.gz'
    ci.write_bytes(gzip.compress(b'chr1\t100\t/path\nchrX\t90\t/path\n'))
    cb.write_bytes(gzip.compress(b'chr1\t0\t40\tp1\tgneg\nchr1\t40\t45\tp11\tacen\n'
                                 b'chr1\t45\t50\tq11\tacen\nchr1\t50\t100\tq1\tgneg\n'))
    sizes = read_chrominfo(ci, ['chr1'])
    table, centers = read_cytobands(cb, sizes)
    assert sizes == {'chr1': 100} and centers == {'chr1': (40, 50)}
    assert len(table) == 4


def test_same_provider_mirror_and_checksum_checkpoint_are_explicit(tmp_path, monkeypatch):
    import src.reference as reference
    calls = []
    primary, mirror = 'https://primary.test/file.gz', 'https://mirror.test/file.gz'
    def download(url, path, network):
        calls.append(url)
        if url == primary:
            raise RemoteResourceError(f'URL={url} HTTP=NO_HTTP_RESPONSE TLS mismatch')
        path.write_bytes(b'test reference only')
        write_json({'url': url, 'sha256': sha256(path), 'http_status': 200}, path.with_suffix('.gz.json'))
    monkeypatch.setattr(reference, 'small_download', download)
    spec = dict(filename='test.gz', url=primary, mirrors=[mirror])
    path, stamp = fetch_reference_resource(spec, tmp_path, {})
    assert stamp['mirror_used'] and stamp['requested_url'] == primary and len(stamp['failures']) == 1
    _, cached = fetch_reference_resource(spec, tmp_path, {})
    assert cached['reused_checkpoint'] and calls == [primary, mirror]
    path.write_bytes(b'corrupted')
    _, repaired = fetch_reference_resource(spec, tmp_path, {})
    assert not repaired['reused_checkpoint'] and calls == [primary, mirror, primary, mirror]
