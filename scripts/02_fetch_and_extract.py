from _bootstrap import parser
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import logging
from pathlib import Path
import time
import pandas as pd

from src.common import load_config, diagnostic_report, write_json, write_parquet, sha256, fingerprint
from src.stream_extractor import extract_sample
from src.reference import load_reference_windows
from src.cohort import (published_metadata, inspect_live_sources, select_cohort,
                         validate_cohort_size, import_published_fallback)


def cohort_table(meta, features, threshold):
    if features.empty:
        ready = meta.iloc[:0].copy()
        ready['accepted_fragments'] = pd.Series(dtype='int64')
    else:
        ready = meta.merge(features[['sample_id', 'accepted_fragments']], on='sample_id', validate='one_to_one')
    records = []
    for center, group in meta.groupby('center', sort=True):
        measured = ready[ready.center.eq(center)]
        good = measured[measured.accepted_fragments.ge(threshold)]
        fail = int(measured.accepted_fragments.lt(threshold).sum())
        records.append(dict(center=center, selected_samples=len(group), selected_patients=group.patient_id.nunique(),
                  selected_cancer=int(group.label.eq(1).sum()), selected_healthy=int(group.label.eq(0).sum()),
                  extracted_samples=len(measured), pending_samples=len(group)-len(measured), retained_samples=len(good),
                  retained_cancer=int(good.label.eq(1).sum()), retained_healthy=int(good.label.eq(0).sum()),
                  cancer_fraction=float(good.label.mean()) if len(good) else None,
                  median_effective_fragments=float(good.accepted_fragments.median()) if len(good) else None,
                  qc_exclusions=fail, qc_exclusion_rate=fail/len(measured) if len(measured) else None,
                  median_original_sequencing_depth=float(group.sequencing_depth.median())))
    return pd.DataFrame(records)


def main():
    p = parser('Step 2: recover public clinical metadata and stream six samples concurrently')
    p.add_argument('--prepare-only', action='store_true', help='Recover and screen real metadata without fragment transfers')
    p.add_argument('--skip-live-probe', action='store_true', help='Reuse published source directly after an earlier documented API failure')
    p.add_argument('--max-new-samples', type=int, help='Explicit diagnostic subset; never publish it as a complete cohort')
    args = p.parse_args()
    cfg = load_config(args.config)
    results, processed = Path(cfg['paths']['results']), Path(cfg['paths']['processed'])
    metadata = Path(cfg['_root']) / 'data/metadata'
    for folder in (results/'reports', results/'tables', processed, metadata):
        folder.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(results/'reports/step_02_download.log', encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    logger = logging.getLogger('src.stream_extractor')
    logger.addHandler(handler)
    try:
        with diagnostic_report(cfg, 2) as report:
            bins = load_reference_windows(cfg)
            if not args.skip_live_probe:
                report['live_source_access'] = inspect_live_sources(cfg)
            raw = published_metadata(cfg)
            meta, screening = select_cohort(raw, cfg)
            write_parquet(screening, metadata/'cohort_screening.parquet')
            write_parquet(meta, metadata/'cohort_selected.parquet')
            write_parquet(meta, cfg['paths']['metadata'])
            meta.to_csv(metadata/'cohort_selected.csv', index=False)
            report.update(metadata_snapshot_rows=len(raw), selected_samples=len(meta),
                          selected_patients=meta.patient_id.nunique(),
                          counts_per_center=meta.groupby('center').label.agg(['count','sum']).to_dict('index'),
                          center_definition='Published study proxy; not a confirmed physical collection center',
                          depth_policy=f'Original plasma WGS depth {cfg["cohort"]["min_depth"]}--{cfg["cohort"]["max_depth"]}x; no random depth values or downsampling',
                          repeated_patient_specimens=int(meta.patient_id.duplicated().sum()),
                          raw_files_stored=0, metadata_missingness=meta.isna().mean().to_dict())
            print('Selected studies: '+str(meta.groupby('center').size().to_dict()), flush=True)
            if args.prepare_only:
                report['status'] = 'metadata_prepared_extraction_not_run'
                return
            blacklist = Path(cfg['paths']['reference'])/'hg19.blacklist.bed.gz'
            checkpoint = processed/cfg['extraction']['checkpoint_file']
            contract = fingerprint({'selected': meta[['sample_id','patient_id','center','label','fragment_url']].to_dict('list'),
                        'windows_sha256': sha256(Path(cfg['paths']['reference'])/cfg['reference']['window_table']),
                        'features': cfg['features'], 'extractor': 2, 'blacklist': sha256(blacklist)})
            chunks, failures, qcs = [], [], []
            if checkpoint.exists() and checkpoint.with_suffix('.json').exists():
                saved = json.loads(checkpoint.with_suffix('.json').read_text())
                if saved['contract'] == contract and saved['sha256'] == sha256(checkpoint):
                    chunks = [pd.read_parquet(checkpoint)]
            completed = set(chunks[0].sample_id) if chunks else set()
            pending = meta[~meta.sample_id.isin(completed)]
            if args.max_new_samples is not None:
                pending = pending.iloc[:args.max_new_samples]
            print(f'Resume: {len(completed)} complete; scheduling {len(pending)} with 6 workers', flush=True)
            if cfg['extraction']['workers'] != 6:
                raise ValueError('Step 2 requires exactly six concurrent extraction threads')
            def snapshot():
                matrix = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
                if not matrix.empty:
                    matrix = matrix.drop_duplicates('sample_id').sort_values('sample_id').reset_index(drop=True)
                    write_parquet(matrix, checkpoint)
                    write_json({'contract': contract, 'sha256': sha256(checkpoint), 'completed': len(matrix)}, checkpoint.with_suffix('.json'))
                table = cohort_table(meta, matrix, cfg['extraction']['min_effective_fragments'])
                table.to_csv(results/'tables/Table1_Cohort_Characteristics.csv', index=False)
                records = json.loads(table.to_json(orient='records'))
                write_json({'selected': len(meta), 'completed': len(matrix), 'failed_attempts': len(failures),
                            'pending': len(meta)-len(matrix), 'state': 'running' if len(matrix)<len(meta) else 'extracted',
                            'study_counts': records}, results/'reports/step_02_progress.json')
                (results/'reports/step_02_report.md').write_text(
                    '# Step 02 — extraction in progress\n\n'
                    f'Real selected samples: {len(meta)}; completed: {len(matrix)}; pending: {len(meta)-len(matrix)}.\n\n'
                    'Canonical cleaned cohort is published only after all selected samples complete and the cohort-size/QC checks pass.\n\n'
                    + table.to_markdown(index=False) + '\n', encoding='utf-8')
                return matrix, table
            start = time.perf_counter()
            try:
                with ThreadPoolExecutor(max_workers=6, thread_name_prefix='cfdna-stream') as executor:
                    future_rows = {executor.submit(extract_sample, row, bins, blacklist, cfg): row.sample_id
                                   for row in pending.itertuples(index=False)}
                    for future in as_completed(future_rows):
                        sample = future_rows[future]
                        try:
                            frame, provenance = future.result()
                            chunks.append(frame)
                            qcs.append({'sample_id': sample, **provenance['qc']})
                            completed.add(sample)
                            print(f'[{len(completed)}/{len(meta)}] {sample}: {int(frame.accepted_fragments.iloc[0]):,} fragments; elapsed={time.perf_counter()-start:.1f}s', flush=True)
                        except Exception as exc:
                            failures.append({'sample_id': sample, 'error': f'{type(exc).__name__}: {exc}'})
                            print(f'FAILED {sample}: {exc}', flush=True)
                        if (len(completed)+len(failures)) % cfg['extraction']['checkpoint_every'] == 0:
                            snapshot()
            finally:
                matrix, table = snapshot()
                write_json(failures, results/'reports/step_02_sample_failures.json')
                if qcs:
                    write_parquet(pd.DataFrame(qcs), processed/'extraction_qc.parquet')
            if failures:
                try:
                    report['published_fallback'] = import_published_fallback(cfg, f'{len(failures)} remote extraction failures')
                except Exception as exc:
                    report['fallback_error'] = f'{type(exc).__name__}: {exc}'
            report['cohort_characteristics'] = json.loads(table.to_json(orient='records'))
            report['extracted_samples'] = len(matrix)
            report['_markdown'] = ('## 真实来源与队列\n\n'
                '元数据来自公开论文 UNITE 的 data file S1 / s2，保留原研究、患者和标本编号。'
                '研究作为中心代理；Healthy 标签严格匹配，排除肝硬化、乙肝、移植、炎症及自身免疫病例。'
                '同一患者的不同标本保留身份以供患者分组交叉验证；未将其当作独立患者。\n\n'
                + table.to_markdown(index=False) + '\n\n'
                '全局长度分布使用全常染色体 MAPQ≥30、80–400 bp 且不与黑名单重叠的片段；'
                '窗口特征仅使用 Step 1 保留窗口的中点归属。T=S+L；IQR 为窗口内全部 80–400 bp 长度的线性插值四分位距。'
                '10 bp 频率采用 80–89、90–99、…、390–399、400；最后一格为端点。'
                '单样本与每 25 例全局检查点均原子保存。\n\n')
            if len(matrix) != len(meta):
                report['status'] = 'partial_extraction_checkpoint_saved'
                report['incomplete_reason'] = 'Some selected samples have not completed; canonical final matrix is NOT published'
                return
            good_ids = matrix.loc[matrix.accepted_fragments.ge(cfg['extraction']['min_effective_fragments']), 'sample_id']
            clean = meta[meta.sample_id.isin(good_ids)].copy()
            validate_cohort_size(clean, cfg)
            matrix = matrix[matrix.sample_id.isin(good_ids)].reset_index(drop=True)
            final = processed/cfg['extraction']['final_file']
            write_parquet(matrix, final)
            write_parquet(matrix, processed/'features.parquet')
            write_parquet(clean, cfg['paths']['metadata'])
            clean = clean.merge(matrix[['sample_id','accepted_fragments']], on='sample_id', validate='one_to_one')
            clean.to_csv(metadata/'cohort_summary.csv', index=False)
            for path in (final, processed/'features.parquet'):
                write_json({'metadata_sha256': sha256(cfg['paths']['metadata']), 'feature_config': fingerprint(cfg['features']),
                        'features_sha256': sha256(path), 'reference_sha256': {'sizes': sha256(blacklist.parent/'hg19.chrom.sizes'),
                        'blacklist': sha256(blacklist)}, 'window_sha256': sha256(blacklist.parent/cfg['reference']['window_table']),
                        'contract': contract, 'n_samples': len(matrix)}, path.with_suffix('.json'))
            report['final_samples'] = len(matrix)
            report['feature_missingness'] = matrix.isna().mean().to_dict()
            report['statistics'] = 'Descriptive cohort/QC statistics only; no case-control hypothesis test or model performance generated'
    finally:
        logger.removeHandler(handler)
        handler.close()


if __name__ == '__main__':
    main()
