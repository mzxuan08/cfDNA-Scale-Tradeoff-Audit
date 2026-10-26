"""Source-traceable hg19 5 Mb windows; all exclusions precede cohort analysis."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .common import fingerprint, sha256, write_json, write_parquet
from .stream_extractor import Blacklist, RemoteResourceError, small_download, make_bins


def initialize_directories(root):
    for relative in ('config', 'data/ref', 'data/metadata', 'data/processed/checkpoints',
                     'src', 'scripts', 'tests', 'results/tables', 'results/figures',
                     'results/reports', 'results/checkpoints'):
        path = Path(root) / relative
        path.mkdir(parents=True, exist_ok=True)
        if not any(path.iterdir()):
            (path / '.gitkeep').touch()


def fetch_reference_resource(spec, directory, network):
    """Retry each configured source, then an explicit same-provider mirror.

    Persist the exact primary URL, the successful URL, and earlier failures.
    Reuse only an intact checkpoint belonging to the current source list.
    """
    path = Path(directory) / spec['filename']
    stamp = path.with_suffix(path.suffix + '.json')
    sources = [spec['url'], *spec.get('mirrors', [])]
    if path.exists() and stamp.exists():
        previous = json.loads(stamp.read_text(encoding='utf-8'))
        if (previous.get('requested_url') == spec['url'] and previous.get('url') in sources
                and previous.get('sha256') == sha256(path)):
            return path, {**previous, 'reused_checkpoint': True}
    failures = []
    for url in sources:
        try:
            small_download(url, path, network)
        except RemoteResourceError as exc:
            failures.append({'url': url, 'error': str(exc)})
            if url == sources[-1]:
                write_json({'requested_url': spec['url'], 'failures': failures},
                           path.with_suffix(path.suffix + '.failure.json'))
                raise
            continue
        record = json.loads(stamp.read_text(encoding='utf-8'))
        record.update(requested_url=spec['url'], failures=failures,
                      mirror_used=url != spec['url'], reused_checkpoint=False)
        write_json(record, stamp)
        return path, record
    raise AssertionError('At least one reference URL is required')


def read_chrominfo(path, chromosomes):
    table = pd.read_csv(path, sep='\t', header=None, names=['chrom', 'size', 'sequence_file'],
                        usecols=[0, 1, 2], dtype={'chrom': str, 'size': 'int64'})
    table = table[table.chrom.isin(chromosomes)].copy()
    if table.chrom.duplicated().any() or set(table.chrom) != set(chromosomes) or table['size'].le(0).any():
        raise ValueError('chromInfo must contain exactly one positive size for each requested chromosome')
    return table.set_index('chrom').loc[chromosomes, 'size'].to_dict()


def read_cytobands(path, chrom_sizes):
    table = pd.read_csv(path, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'band', 'stain'],
                        dtype={'chrom': str, 'start': 'int64', 'end': 'int64', 'band': str, 'stain': str})
    table = table[table.chrom.isin(chrom_sizes)].copy()
    centromeres = {}
    for chrom, size in chrom_sizes.items():
        group = table[table.chrom.eq(chrom)].sort_values('start')
        if (group.empty or int(group.start.iloc[0]) != 0 or int(group.end.iloc[-1]) != size
                or (group.end <= group.start).any()
                or not np.array_equal(group.start.to_numpy()[1:], group.end.to_numpy()[:-1])):
            raise ValueError(f'Invalid/noncontiguous cytobands for {chrom}')
        acen = group[group.stain.eq('acen')]
        if len(acen) != 2 or set(acen.band.str[0]) != {'p', 'q'} or acen.end.iloc[0] != acen.start.iloc[1]:
            raise ValueError(f'Expected adjacent p/q acen bands for {chrom}')
        centromeres[chrom] = (int(acen.start.min()), int(acen.end.max()))
    return table, centromeres


def union_intervals(intervals):
    merged = []
    for start, end in sorted((int(a), int(b)) for a, b in intervals):
        if start < 0 or end <= start:
            raise ValueError('Invalid genomic interval')
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def overlap_bp(start, end, intervals):
    """Count union overlap of zero-based, half-open intervals without double counting."""
    return sum(max(0, min(end, b) - max(start, a)) for a, b in union_intervals(intervals))


def build_quality_windows(chrom_sizes, cytobands, centromeres, blacklist, bin_size, filters):
    if bin_size <= 0 or filters['telomere_buffer_bp'] < 0 or not 0 <= filters['max_blacklist_fraction'] <= 1:
        raise ValueError('Invalid bin/filter settings')
    records = []
    proxy_stains = set(filters.get('excluded_cytoband_stains', []))
    for chrom, size in chrom_sizes.items():
        cstart, cend = centromeres[chrom]
        if not 0 < cstart < cend < size:
            raise ValueError(f'Invalid centromere coordinates for {chrom}')
        cb = cytobands[cytobands.chrom.eq(chrom) & cytobands.stain.isin(proxy_stains)]
        proxies = union_intervals(cb[['start', 'end']].itertuples(index=False, name=None))
        masks = blacklist.intervals.get(chrom, [])
        for start in range(0, size, bin_size):
            end = min(start + bin_size, size)
            center_overlap = max(0, min(end, cend) - max(start, cstart))
            black_bp = overlap_bp(start, end, masks)
            proxy_bp = overlap_bp(start, end, proxies)
            # Every retained interval fits entirely in one arm.
            arm = 'p' if end <= cstart else ('q' if start >= cend else None)
            terminal_buffer = filters['telomere_buffer_bp']
            flags = {
                'partial': filters.get('full_windows_only', True) and end - start != bin_size,
                'centromere': center_overlap > 0,
                'terminal': start < terminal_buffer or end > size - terminal_buffer,
                'low_mappability_proxy': proxy_bp > 0,
                # Strictly GREATER than 20%; exactly 20% is allowed.
                'blacklist': black_bp > filters['max_blacklist_fraction'] * (end - start),
            }
            distance_c = cstart - end if arm == 'p' else (start - cend if arm == 'q' else None)
            distance_t = start if arm == 'p' else (size - end if arm == 'q' else None)
            mid = (start + end) // 2
            records.append(dict(chrom=chrom, start=start, end=end, bin_id=f'{chrom}:{start}-{end}',
                 length_bp=end - start, arm=arm, chromosome_arm=f'{chrom}{arm}' if arm else None,
                 chrom_size_bp=size, centromere_start_bp=cstart, centromere_end_bp=cend,
                 midpoint_bp=mid, distance_to_telomere_bp=distance_t,
                 distance_to_centromere_bp=distance_c,
                 nearest_telomere_distance_bp=min(start, size - end),
                 midpoint_to_telomere_bp=(mid if arm == 'p' else (size - mid if arm == 'q' else None)),
                 midpoint_to_centromere_bp=(cstart - mid if arm == 'p' else (mid - cend if arm == 'q' else None)),
                 blacklist_overlap_bp=black_bp, blacklist_fraction=black_bp / (end - start),
                 centromere_overlap_bp=center_overlap, low_mappability_proxy_overlap_bp=proxy_bp,
                 retained=not any(flags.values()), exclusion_reasons=';'.join(k for k, v in flags.items() if v),
                 **{f'excluded_{key}': value for key, value in flags.items()}))
    candidates = pd.DataFrame(records)
    for column in ('distance_to_telomere_bp', 'distance_to_centromere_bp',
                   'midpoint_to_telomere_bp', 'midpoint_to_centromere_bp'):
        candidates[column] = candidates[column].astype('Int64')
    candidates['blacklist_fraction'] = candidates.blacklist_fraction.astype('float32')
    retained = candidates[candidates.retained].reset_index(drop=True)
    if retained.empty or retained.arm.isna().any():
        raise ValueError('No retained windows or ambiguous chromosome arm')
    validate_windows(retained, chrom_sizes, bin_size, filters)
    return retained, candidates


def validate_windows(bins, chrom_sizes, bin_size, filters):
    if bins.bin_id.duplicated().any() or not bins.chrom.isin(chrom_sizes).all():
        raise ValueError('Duplicate windows or unexpected chromosome')
    if filters.get('full_windows_only', True) and not (bins.end - bins.start).eq(bin_size).all():
        raise ValueError('Retained windows must have exactly the requested size')
    for chrom, group in bins.groupby('chrom', sort=False):
        group = group.sort_values('start')
        if ((group.start < 0).any() or (group.end > chrom_sizes[chrom]).any()
                or (group.end.to_numpy()[:-1] > group.start.to_numpy()[1:]).any()):
            raise ValueError('Overlapping/out-of-bounds windows')
    if not bins.arm.isin(['p', 'q']).all() or bins.centromere_overlap_bp.gt(0).any():
        raise ValueError('Retained window crosses a centromere')
    # Compare integer overlap to avoid a float32 rounding decision at exactly 20%.
    if (bins.blacklist_overlap_bp > filters['max_blacklist_fraction'] * bins.length_bp).any():
        raise ValueError('Retained window exceeds blacklist fraction')


def arm_distribution(bins, chromosomes):
    index = pd.MultiIndex.from_product([chromosomes, ['p', 'q']], names=['chrom', 'arm'])
    counts = bins.groupby(['chrom', 'arm']).size().reindex(index, fill_value=0).rename('windows').reset_index()
    totals = bins.assign(_bp=bins.end - bins.start).groupby(['chrom', 'arm'])._bp.sum().reindex(index, fill_value=0)
    counts['base_pairs'] = totals.to_numpy(dtype='int64')
    return counts


def window_contract(cfg):
    return fingerprint({'genome': cfg['reference']['genome'], 'assembly': cfg['reference']['assembly'],
                        'chromosomes': cfg['features']['chromosomes'], 'bin_size': cfg['features']['bin_size'],
                        'filters': cfg['reference']['window_filter'], 'algorithm': 'aligned_full_5mb_v1'})


def load_reference_windows(cfg):
    ref = Path(cfg['paths']['reference'])
    path = ref / cfg['reference']['window_table']
    stamp = json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
    if stamp['contract'] != window_contract(cfg) or stamp['sha256'] != sha256(path):
        raise ValueError('Reference windows/config changed; rerun step 01')
    for name, digest in stamp['source_sha256'].items():
        if sha256(ref / name) != digest:
            raise ValueError('Reference source changed; rerun step 01')
    return pd.read_parquet(path)
