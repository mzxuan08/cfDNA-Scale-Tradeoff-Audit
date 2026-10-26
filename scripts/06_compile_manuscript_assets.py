from _bootstrap import parser
import json
from pathlib import Path
import pandas as pd

from src.common import load_config, diagnostic_report, load_dataset, analysis_signature, sha256, write_json


def main():
    p=parser('Compile source-verified tables and figure inventory without inventing findings')
    p.add_argument('--text-only',action='store_true',help='Align manuscript scope using existing results; no fitting or extraction')
    args=p.parse_args()
    cfg = load_config(args.config)
    if args.text_only:
        from src.manuscript_alignment import finalize_scope_alignment
        finalize_scope_alignment(cfg)
        return
    if (Path(cfg['paths']['results'])/'tables/Step7_OOF_Predictions.parquet').exists():
        first=Path(cfg['paths']['results'])/'reports/step7_first_complete_execution.json'
        current=Path(cfg['paths']['results'])/'reports/step_07_report.json'
        if current.exists() and not first.exists():
            previous=json.loads(current.read_text())
            if previous.get('status')=='complete' and not previous.get('render_only'):
                write_json(previous,first)
        from src.multiscale_benchmark import run_step7
        run_step7(cfg,render_only=True)
        return
    from src.study import clinical_available
    if clinical_available(cfg):
        from src.manuscript import compile_clinical_assets,refresh_clinical_manifest
        with diagnostic_report(cfg, 6) as report:
            report['statistics']=compile_clinical_assets(cfg)
        # Include the final report written by the context manager, not its prior version.
        refresh_clinical_manifest(cfg)
        return
    with diagnostic_report(cfg, 6) as report:
        meta, X = load_dataset(cfg)
        root = Path(cfg['paths']['results'])
        table = root / 'tables'
        signature = analysis_signature(cfg)
        for name in ('audit', 'cv', 'biology'):
            record = json.loads((table / f'{name}_provenance.json').read_text())
            if record['signature'] != signature:
                raise ValueError(f'{name} output is stale; rerun corresponding step')
            for file, digest in record['files'].items():
                if sha256(table / file) != digest:
                    raise ValueError(f'Output checksum mismatch: {file}')
        cohort = meta.groupby('center', as_index=False).agg(samples=('sample_id', 'size'),
                    patients=('patient_id', 'nunique'), cases=('label', 'sum'))
        cohort['controls'] = cohort.samples - cohort.cases
        cohort.to_csv(table / 'Table1_cohort.csv', index=False)
        cohort.to_latex(table / 'Table1_cohort.tex', index=False, escape=True)
        (table / 'Table1_cohort.md').write_text(cohort.to_markdown(index=False), encoding='utf-8')
        summary = pd.read_parquet(table / 'cv_summary.parquet')
        complete = summary[summary.status.eq('complete')]
        complete.to_latex(table / 'Table2_performance.tex', index=False, escape=True, float_format='%.4f')
        (table / 'Table2_performance.md').write_text(complete.to_markdown(index=False), encoding='utf-8')
        unavailable = json.loads((table / 'unavailable_comparisons.json').read_text())
        content = '# Manuscript asset inventory\n\n'
        content += f'Analyzed {len(meta)} verified samples, {meta.patient_id.nunique()} patients, {meta.center.nunique()} declared centers/proxies.\n\n'
        content += 'Scientific hypotheses are not predetermined conclusions. Review source identity and center confounding before interpretation.\n\n'
        content += complete.to_markdown(index=False) + '\n\n'
        content += 'OOF score bootstrap intervals are descriptive and conditional on fitted models. DeLong is reserved for independent external subjects and fixed models.\n\n'
        content += f'Unavailable fold/model/method entries: {len(unavailable)}. See tables/unavailable_comparisons.json.\n\n'
        content += 'LSN radius and ICC cutoff are frozen design parameters; the pipeline is not literally hyperparameter-free. No restoration or clinical screening claim is generated automatically.\n'
        (root / 'reports/manuscript_assets.md').write_text(content, encoding='utf-8')
        files = sorted(p for folder in ('tables', 'figures') for p in (root / folder).glob('*') if p.is_file())
        inventory = {str(p.relative_to(root)): {'bytes': p.stat().st_size, 'sha256': sha256(p)} for p in files}
        write_json({'analysis_signature': signature, 'files': inventory}, root / 'asset_manifest.json')
        report['statistics'] = {'complete_comparisons': len(complete), 'unavailable_entries': len(unavailable),
                                 'source': 'tables/cv_summary.parquet'}
        report['assets'] = inventory
        report['feature_missingness'] = X.isna().mean().to_dict()


if __name__ == '__main__':
    main()
