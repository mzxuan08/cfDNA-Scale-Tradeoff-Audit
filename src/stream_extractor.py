"""Bounded-memory HTTP gzip/BGZF extraction. No raw interval file touches disk."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import logging
from pathlib import Path
import time
import zlib

import numpy as np
import pandas as pd
import requests
from urllib3.exceptions import HTTPError as UrllibHTTPError

from .common import fingerprint, sha256, write_json, write_parquet

LOG = logging.getLogger(__name__)


class RemoteResourceError(RuntimeError):
    pass


def remote_operation(url, network, consume, session=None):
    """Restart the entire aggregation after a stream failure, never double-count.

    max_retries=5 means one initial attempt plus at most five retries. The callback
    must create fresh aggregation buffers on every call; partial vectors never commit.
    """
    own = session is None
    session = session or requests.Session()
    try:
        for attempt in range(network['max_retries'] + 1):
            status = 'NO_HTTP_RESPONSE'
            try:
                with session.get(url, stream=True, timeout=(network['connect_timeout'],
                                                          network['read_timeout'])) as response:
                    status = response.status_code
                    response.raise_for_status()
                    response.raw.decode_content = False
                    return consume(response)
            except (requests.RequestException, UrllibHTTPError, OSError, EOFError, zlib.error) as exc:
                message = f'URL={url} HTTP={status} attempt={attempt + 1}: {type(exc).__name__}: {exc}'
                LOG.error(message)
                if attempt == network['max_retries']:
                    raise RemoteResourceError(message) from exc
                time.sleep(network['backoff_seconds'] * 2 ** attempt)
    finally:
        if own:
            session.close()


def small_download(url, path, network, limit=10 * 1024 * 1024):
    """Only small reference/metadata files may be materialized."""
    def consume(response):
        result = bytearray()
        for block in response.iter_content(65536):
            result.extend(block)
            if len(result) > limit:
                raise ValueError(f'Small reference exceeds {limit} bytes: {url}')
        return bytes(result), {'http_status': response.status_code, 'response_url': response.url}
    payload, transport = remote_operation(url, network, consume)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_bytes(payload)
    tmp.replace(path)
    write_json({'url': url, 'sha256': sha256(path), 'bytes': len(payload), **transport},
               path.with_suffix(path.suffix + '.json'))


def make_bins(chrom_sizes, size, chromosomes):
    sizes = {}
    for line in Path(chrom_sizes).read_text().splitlines():
        chrom, length = line.split()[:2]
        sizes[chrom] = int(length)
    rows = []
    for chrom in chromosomes:
        if chrom not in sizes:
            raise ValueError(f'Missing chromosome size: {chrom}')
        for start in range(0, sizes[chrom], size):
            end = min(start + size, sizes[chrom])
            rows.append((chrom, start, end, f'{chrom}:{start}-{end}'))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'bin_id'])


class Blacklist:
    def __init__(self, path=None):
        self.intervals = {}
        if path is None:
            return
        opener = gzip.open if str(path).endswith('.gz') else open
        intervals = {}
        with opener(path, 'rt') as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                chrom, start, end, *_ = line.split()
                intervals.setdefault(chrom, []).append((int(start), int(end)))
        for chrom, items in intervals.items():
            merged = []
            for start, end in sorted(items):
                if end <= start or start < 0:
                    raise ValueError('Invalid blacklist interval')
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.intervals[chrom] = np.asarray(merged, dtype=np.int64)

    def overlaps(self, chrom, start, end):
        intervals = self.intervals.get(chrom)
        if intervals is None:
            return False
        i = np.searchsorted(intervals[:, 0], end, side='left') - 1
        return i >= 0 and intervals[i, 1] > start


class HashReader:
    def __init__(self, raw):
        self.raw, self.digest, self.count = raw, hashlib.sha256(), 0

    def read(self, size=-1):
        data = self.raw.read(size)
        self.digest.update(data)
        self.count += len(data)
        return data


def histogram_quantiles(hist, lower=80):
    """Exact NumPy-linear quartiles from integer histograms; no fragment lists."""
    hist = np.asarray(hist, dtype=np.int64)
    output = np.full((len(hist), 3), np.nan, dtype=np.float32)
    for i, counts in enumerate(hist):
        n = int(counts.sum())
        if n == 0:
            continue
        cumulative = np.cumsum(counts)
        for j, quantile in enumerate((.25, .5, .75)):
            position = (n - 1) * quantile
            a, b = int(np.floor(position)), int(np.ceil(position))
            av, bv = np.searchsorted(cumulative, [a, b], side='right') + lower
            output[i, j] = av + (bv - av) * (position - a)
    return output


def features_from_histograms(hist, per_bin, bins, params):
    lower, upper = params['histogram']
    short = per_bin[:, params['short'][0] - lower:params['short'][1] - lower + 1].sum(axis=1)
    long = per_bin[:, params['long'][0] - lower:params['long'][1] - lower + 1].sum(axis=1)
    ratios = np.full(len(bins), np.nan, dtype=np.float32)
    valid = long >= params['min_long_count']
    a = params['pseudocount']
    ratios[valid] = (short[valid] + a) / (long[valid] + a)
    quartiles = histogram_quantiles(per_bin, lower)
    values = {}
    for j, name in enumerate(bins.bin_id):
        values.update({f'short_count__{name}': np.float32(short[j]),
                       f'long_count__{name}': np.float32(long[j]),
                       f'total_count__{name}': np.float32(short[j] + long[j]),
                       f'ratio__{name}': ratios[j],
                       f'iqr__{name}': np.float32(quartiles[j, 2] - quartiles[j, 0])})
    total = int(hist.sum())
    for length in range(lower, upper + 1):
        values[f'length__{length}'] = np.float32(hist[length - lower] / total) if total else np.float32(np.nan)
    # Disjoint 10 bp bins; the inclusive 400 bp endpoint has its own singleton bin.
    for start in range(lower, upper + 1, 10):
        end = min(start + 9, upper)
        values[f'length10__{start}_{end}'] = np.float32(hist[start-lower:end-lower+1].sum() / total) if total else np.float32(np.nan)
    return values


def aggregate_fragment_batches(handle, bins, params, blacklist, chunk_rows=250000):
    """C-parser + vectorized aggregation lets six IO threads avoid per-row Python.

    Global histograms use all filtered autosomal fragments. Window features use
    only the retained Step 1 intervals; holes never shift fragments to neighbors.
    IQR uses every 80--400 bp high-quality fragment whose midpoint is in that bin.
    """
    bins = bins.reset_index(drop=True)
    chromosomes = params['chromosomes']
    mapping, chrom_sizes = [], []
    for chrom in chromosomes:
        group = bins[bins.chrom.eq(chrom)]
        if group.empty:
            raise ValueError(f'No retained windows for {chrom}')
        size = int(group.chrom_size_bp.iloc[0]) if 'chrom_size_bp' in group else int(group.end.max())
        chrom_sizes.append(size)
        slots = np.full((size + params['bin_size'] - 1) // params['bin_size'], -1, dtype=np.int32)
        slots[(group.start // params['bin_size']).to_numpy(dtype=int)] = group.index.to_numpy()
        mapping.append(slots)
    lower, upper = params['histogram']
    width = upper - lower + 1
    hist = np.zeros(width, dtype=np.int64)
    per_bin = np.zeros((len(bins), width), dtype=np.int64)
    qc = dict(rows=0, accepted=0, excluded_chrom=0, low_mapq=0, blacklisted=0,
              length_outside=0, renamed_chromosome=0, outside_retained_windows=0, assigned_to_windows=0)
    reader = pd.read_csv(handle, sep='\t', header=None, comment='#', usecols=[0, 1, 2, 3],
                         names=['chrom', 'start', 'end', 'mapq'], chunksize=chunk_rows,
                         dtype={'chrom': str, 'start': 'int64', 'end': 'int64', 'mapq': 'int32'},
                         compression=None, engine='c')
    for frame in reader:
        qc['rows'] += len(frame)
        names = frame.chrom
        if params.get('normalize_chr_prefix', False):
            numeric = names.isin([c.removeprefix('chr') for c in chromosomes])
            qc['renamed_chromosome'] += int(numeric.sum())
            names = names.where(~numeric, 'chr' + names)
        codes = pd.Categorical(names, categories=chromosomes).codes
        start, end, mapq = (frame[c].to_numpy() for c in ('start', 'end', 'mapq'))
        if (start < 0).any() or (end <= start).any():
            raise ValueError('Malformed fragment coordinates')
        valid = codes >= 0
        qc['excluded_chrom'] += int((~valid).sum())
        if np.any(valid & (end > np.asarray(chrom_sizes)[np.maximum(codes, 0)])):
            raise ValueError('Fragment exceeds chromosome size: possible genome mismatch')
        qc['low_mapq'] += int(np.sum(valid & (mapq < params['mapq_min'])))
        valid &= mapq >= params['mapq_min']
        lengths = end - start
        qc['length_outside'] += int(np.sum(valid & ((lengths < lower) | (lengths > upper))))
        valid &= (lengths >= lower) & (lengths <= upper)
        for code in np.unique(codes[valid]):
            selected = np.flatnonzero(valid & (codes == code))
            masks = blacklist.intervals.get(chromosomes[code])
            if masks is not None:
                index = np.searchsorted(masks[:, 0], end[selected], side='left') - 1
                overlap = (index >= 0) & (masks[np.maximum(index, 0), 1] > start[selected])
                qc['blacklisted'] += int(overlap.sum())
                selected = selected[~overlap]
            if not len(selected):
                continue
            qc['accepted'] += len(selected)
            length_code = lengths[selected] - lower
            hist += np.bincount(length_code, minlength=width)
            slot = ((start[selected] + end[selected]) // 2) // params['bin_size']
            target = mapping[code][slot]
            keep = target >= 0
            qc['outside_retained_windows'] += int((~keep).sum())
            qc['assigned_to_windows'] += int(keep.sum())
            combined = target[keep] * width + length_code[keep]
            per_bin += np.bincount(combined, minlength=len(bins) * width).reshape(len(bins), width)
    return features_from_histograms(hist, per_bin, bins, params), qc


def aggregate_fragments(lines, bins, params, blacklist):
    """0-based half-open fragments; midpoint assigns exactly one genomic bin.

    Length bounds are inclusive. Any fragment overlapping a blacklist is excluded.
    Undefined ratios and empty histograms stay NaN until train-only imputation.
    """
    n = len(bins)
    short = np.zeros(n, dtype=np.int64)
    long = np.zeros(n, dtype=np.int64)
    lo, hi = params['histogram']
    hist = np.zeros(hi - lo + 1, dtype=np.int64)
    chrom_bins = {c: (g.index.to_numpy(), g.end.max()) for c, g in bins.groupby('chrom', sort=False)}
    qc = dict(rows=0, accepted=0, malformed=0, excluded_chrom=0, low_mapq=0,
              blacklisted=0, out_of_bounds=0, length_outside=0, renamed_chromosome=0)
    for line in lines:
        if isinstance(line, bytes):
            line = line.decode('ascii')
        if not line.strip() or line.startswith('#'):
            continue
        qc['rows'] += 1
        parts = line.split()
        if len(parts) < 5:
            raise ValueError(f'Fragment row {qc["rows"]}: expected chrom,start,end,mapq,strand')
        chrom = parts[0]
        if params.get('normalize_chr_prefix', False) and chrom not in chrom_bins and ('chr' + chrom) in chrom_bins:
            chrom = 'chr' + chrom
            qc['renamed_chromosome'] += 1
        try:
            start, end, mapq = map(int, parts[1:4])
        except ValueError as exc:
            raise ValueError(f'Invalid numeric fragment fields at row {qc["rows"]}') from exc
        if parts[4] not in ('+', '-') or start < 0 or end <= start:
            raise ValueError(f'Invalid BED3+2 interval at row {qc["rows"]}')
        if chrom not in chrom_bins:
            qc['excluded_chrom'] += 1
            continue
        indices, chrom_end = chrom_bins[chrom]
        if end > chrom_end:
            raise ValueError(f'Coordinate exceeds {chrom} reference length; check genome build')
        if mapq < params['mapq_min']:
            qc['low_mapq'] += 1
            continue
        if blacklist.overlaps(chrom, start, end):
            qc['blacklisted'] += 1
            continue
        length = end - start
        if not lo <= length <= hi:
            qc['length_outside'] += 1
            continue
        qc['accepted'] += 1
        hist[length - lo] += 1
        idx = indices[((start + end) // 2) // params['bin_size']]
        if params['short'][0] <= length <= params['short'][1]:
            short[idx] += 1
        if params['long'][0] <= length <= params['long'][1]:
            long[idx] += 1
    if qc['accepted'] < params['min_accepted_fragments']:
        raise ValueError(f'Insufficient accepted fragments: {qc["accepted"]}; QC={qc}')
    ratio = np.full(n, np.nan, dtype=np.float32)
    valid = long >= params['min_long_count']
    a = params['pseudocount']
    ratio[valid] = ((short[valid] + a) / (long[valid] + a)).astype(np.float32)
    result = {f'ratio__{b}': ratio[i] for i, b in enumerate(bins.bin_id)}
    result.update({f'short_count__{b}': short[i] for i, b in enumerate(bins.bin_id)})
    result.update({f'long_count__{b}': long[i] for i, b in enumerate(bins.bin_id)})
    result.update({f'length__{length}': np.float32(hist[length - lo] / hist.sum())
                   for length in range(lo, hi + 1)})
    return result, qc


def bounded_lines(handle):
    while True:
        line = handle.readline(4097)
        if not line:
            break
        if len(line) > 4096:
            raise ValueError('Unexpectedly long fragment line; possible wrong remote resource')
        yield line


def extract_sample(row, bins, blacklist_path, cfg):
    dest = Path(cfg['paths']['processed']) / 'checkpoints' / (fingerprint(row.sample_id)[:24] + '.parquet')
    provenance_path = dest.with_suffix('.json')
    contract = fingerprint({'sample': row.sample_id, 'url': row.fragment_url, 'genome': row.genome,
                            'params': cfg['features'], 'bins': bins.to_dict('list'),
                            'blacklist_sha256': sha256(blacklist_path), 'extractor_version': 2,
                            'chunk_rows': cfg.get('extraction', {}).get('chunk_rows', 250000)})
    if dest.exists() and provenance_path.exists():
        saved = json.loads(provenance_path.read_text())
        if saved['contract'] == contract and saved['sha256'] == sha256(dest):
            return pd.read_parquet(dest), saved
    blacklist = Blacklist(blacklist_path)
    def consume(response):
        reader = HashReader(response.raw)
        with gzip.GzipFile(fileobj=reader, mode='rb') as gz:
            values, qc = aggregate_fragment_batches(gz, bins, cfg['features'], blacklist,
                          cfg.get('extraction', {}).get('chunk_rows', 250000))
        return values, {**qc, 'compressed_bytes': reader.count,
                        'compressed_sha256': reader.digest.hexdigest(),
                        'etag': response.headers.get('ETag'),
                        'last_modified': response.headers.get('Last-Modified')}
    values, qc = remote_operation(row.fragment_url, cfg['network'], consume)
    values['sample_id'] = row.sample_id
    values['accepted_fragments'] = qc['accepted']
    df = pd.DataFrame([values])
    write_parquet(df, dest)
    provenance = {'contract': contract, 'url': row.fragment_url, 'qc': qc, 'sha256': sha256(dest)}
    write_json(provenance, provenance_path)
    return df, provenance


def stream_gc_reference(url, bins, cfg):
    """Exact called-base GC fractions from reference FASTA, with no FASTA on disk."""
    def consume(response):
        gc_count = np.zeros(len(bins), dtype=np.int64)
        called = np.zeros(len(bins), dtype=np.int64)
        grouped = {c: g.index.to_numpy() for c, g in bins.groupby('chrom', sort=False)}
        sizes = bins.groupby('chrom').end.max().to_dict()
        seen, chrom, pos = set(), None, 0
        reader = HashReader(response.raw)
        with gzip.GzipFile(fileobj=reader) as handle:
            for raw in handle:
                if raw.startswith(b'>'):
                    if chrom in sizes and pos != sizes[chrom]:
                        raise ValueError('Reference FASTA length mismatch')
                    chrom, pos = raw[1:].split()[0].decode(), 0
                    if chrom in grouped:
                        if chrom in seen:
                            raise ValueError('Duplicate reference contig')
                        seen.add(chrom)
                    continue
                seq = raw.strip().upper()
                if chrom not in grouped:
                    continue
                offset = 0
                while offset < len(seq):
                    k = pos // cfg['features']['bin_size']
                    length = min(len(seq) - offset, cfg['features']['bin_size'] - pos % cfg['features']['bin_size'])
                    fragment = seq[offset:offset + length]
                    idx = grouped[chrom][k]
                    gc_count[idx] += fragment.count(b'G') + fragment.count(b'C')
                    called[idx] += sum(fragment.count(b) for b in (b'A', b'C', b'G', b'T'))
                    offset += length
                    pos += length
        if chrom in sizes and pos != sizes[chrom] or seen != set(grouped):
            raise ValueError('Reference FASTA is incomplete')
        result = bins.copy()
        result['gc_fraction'] = np.divide(gc_count, called, out=np.full(len(bins), np.nan), where=called > 0).astype('float32')
        result['called_fraction'] = (called / (bins.end - bins.start)).astype('float32')
        return result, {'url': url, 'compressed_sha256': reader.digest.hexdigest(), 'compressed_bytes': reader.count}
    return remote_operation(url, cfg['network'], consume)
