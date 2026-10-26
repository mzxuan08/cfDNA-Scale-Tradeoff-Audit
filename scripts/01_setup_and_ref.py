from _bootstrap import parser
import json
import logging
from pathlib import Path
import subprocess

from src.common import load_config, diagnostic_report, write_json, write_parquet, sha256, seed_everything
from src.stream_extractor import small_download, make_bins, stream_gc_reference, Blacklist
from src.reference import (initialize_directories, fetch_reference_resource, read_chrominfo,
                           read_cytobands, build_quality_windows, arm_distribution, window_contract)


def main():
    p = parser('Step 1: initialize project and build source-verified hg19 autosomal 5 Mb windows')
    p.add_argument('--discover', action='store_true', help='Optional FinaleDB query; not required for Step 1')
    p.add_argument('--with-gc', action='store_true', help='Optional streamed reference FASTA GC calculation')
    args = p.parse_args()
    cfg = load_config(args.config)
    initialize_directories(cfg['_root'])
    reports = Path(cfg['paths']['results']) / 'reports'
    handler = logging.FileHandler(reports / 'step_01_download.log', encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    logger = logging.getLogger('src.stream_extractor')
    logger.addHandler(handler)
    try:
        with diagnostic_report(cfg, 1) as report:
            if cfg['reference']['genome'] != 'hg19' or cfg['reference']['assembly'] != 'GRCh37':
                raise ValueError('This Step 1 protocol requires hg19 / GRCh37')
            seed_everything(cfg['seed'])
            ref = Path(cfg['paths']['reference'])
            resources, files = {}, {}
            report['references'] = resources
            for key in ('cytoband', 'chrominfo', 'blacklist'):
                files[key], resources[key] = fetch_reference_resource(
                    cfg['reference']['resources'][key], ref, cfg['network'])
                print(f'Reference ready: {files[key].name}', flush=True)
            chromosomes = cfg['features']['chromosomes']
            sizes = read_chrominfo(files['chrominfo'], chromosomes)
            cytobands, centromeres = read_cytobands(files['cytoband'], sizes)
            sizes_path = ref / 'hg19.chrom.sizes'
            sizes_path.write_text(''.join(f'{chrom}\t{size}\n' for chrom, size in sizes.items()), encoding='ascii')
            write_json({'derived_from': files['chrominfo'].name, 'source_sha256': sha256(files['chrominfo']),
                        'sha256': sha256(sizes_path)}, sizes_path.with_suffix('.sizes.json'))
            bins, candidates = build_quality_windows(sizes, cytobands, centromeres, Blacklist(files['blacklist']),
                                      cfg['features']['bin_size'], cfg['reference']['window_filter'])
            output = ref / cfg['reference']['window_table']
            write_parquet(bins, output)
            write_parquet(candidates, ref / 'hg19_5mb_candidates.parquet')
            write_parquet(cytobands, ref / 'hg19_cytobands.parquet')
            # This legacy table is an explicitly unfiltered candidate grid.
            full_grid = make_bins(sizes_path, cfg['features']['bin_size'], chromosomes)
            write_parquet(full_grid, ref / 'hg19.bins.parquet')
            source_hashes = {path.name: sha256(path) for path in files.values()}
            write_json({'genome': 'hg19', 'assembly': 'GRCh37', 'contract': window_contract(cfg),
                        'sha256': sha256(output), 'source_sha256': source_hashes,
                        'n_candidates': len(candidates), 'n_retained': len(bins),
                        'filter_parameters': cfg['reference']['window_filter'],
                        'distance_definition': 'Window edge to same-arm telomere or centromere boundary, in bp; midpoint distances stored separately'},
                       output.with_suffix('.json'))
            arms = arm_distribution(bins, chromosomes)
            write_parquet(arms, ref / 'hg19_5mb_arm_distribution.parquet')
            exclusions = {column.removeprefix('excluded_'): int(candidates[column].sum())
                          for column in candidates if column.startswith('excluded_')}
            report.update(genome='hg19', assembly='GRCh37', n_bins=len(bins),
                          n_candidate_windows=len(candidates), n_retained_windows=len(bins),
                          n_rejected_windows=int((~candidates.retained).sum()),
                          excluded_reason_counts_nonexclusive=exclusions,
                          arm_distribution=arms.to_dict('records'),
                          feature_missingness=bins.isna().mean().to_dict(),
                          statistics='Reference-only step: clinical samples=0; p values, effect sizes and CIs not applicable',
                          raw_fragment_files_stored=0, source_sha256=source_hashes,
                          retained_table_sha256=sha256(output), window_filter=cfg['reference']['window_filter'])
            report['_markdown'] = (
                '## 工程初始化与参考窗口构建\n\n'
                f'基因组：hg19 / GRCh37；随机种子：{cfg["seed"]}；配色：NPG nrc。\n\n'
                f'常染色体原始网格 **{len(candidates)}** 个；剔除 **{len(candidates) - len(bins)}** 个；'
                f'最终保留 **{len(bins)}** 个完整、互不重叠的 5,000,000 bp 窗口。\n\n'
                '筛选固定在标签与队列分析之前，不按目标数量调阈值。黑名单按区间并集计算，严格大于 20% 才剔除。'
                '着丝粒以相邻 p/q 两条 acen 带的并集定义；与其重叠的窗口全部剔除。'
                '染色体两端各 1 Mb 为边界风险代理；gvar/stalk 带为难比对区域代理，均按重叠筛除。'
                '这些代理不是基于实测 mappability 分数的过滤。\n\n'
                '距离以 bp 计：窗口边缘到同一臂端粒或着丝粒边缘的距离；另附窗口中点距离。'
                'chr13/14/15/21/22 的重复性短臂可能没有保留窗口，表中保留零计数。\n\n'
                '## 各染色体臂的窗口分布\n\n' + arms.to_markdown(index=False) + '\n\n'
                '## 可核验记录\n\n'
                '下载来源、备用地址使用、HTTP 状态和哈希见下方 JSON 及 data/ref/ 来源文件。'
                '每个地址最多重试 5 次；配置的同一提供方备用地址也失败后中止。'
                '本步骤未读取临床病例，统计检验不适用。\n\n')
            if args.with_gc:
                gc, provenance = stream_gc_reference(cfg['reference']['gc_urls']['hg19'], full_grid, cfg)
                gc = bins[['bin_id']].merge(gc, on='bin_id', validate='one_to_one')
                write_parquet(gc, ref / 'hg19.bins_gc.parquet')
                write_json(provenance, ref / 'hg19.bins_gc.json')
                report['gc'] = provenance
            else:
                report['gc'] = 'Not requested; GC is outside this Step 1 acceptance criteria'
            if args.discover:
                path = Path(cfg['_root']) / 'data/metadata/finaledb_response.json'
                report['discovery_url'] = cfg['metadata']['discovery_url']
                small_download(cfg['metadata']['discovery_url'], path, cfg['network'])
                report['discovery_response_type'] = type(json.loads(path.read_text(encoding='utf-8'))).__name__
            try:
                gpu = subprocess.run(['nvidia-smi', '--query-gpu=name,memory.total,driver_version', '--format=csv'],
                                     capture_output=True, text=True, check=False)
                report['gpu_inventory'] = gpu.stdout.strip() or gpu.stderr.strip()
            except FileNotFoundError:
                report['gpu_inventory'] = 'nvidia-smi unavailable'
            print(f'Completed: {len(bins)} retained 5 Mb windows; report={reports / "step_01_report.md"}', flush=True)
    finally:
        logger.removeHandler(handler)
        handler.close()


if __name__ == '__main__':
    main()
