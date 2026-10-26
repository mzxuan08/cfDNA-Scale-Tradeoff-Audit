from _bootstrap import parser
from pathlib import Path

from src.common import load_config, diagnostic_report, load_dataset, write_parquet, write_json, analysis_signature, sha256
from src.evaluation import run_benchmarks, summarize_predictions
from src.plotting import fig2, fig3


def main():
    p=parser('Patient-purged Random-CV vs LOCO-CV with fold-specific preprocessing')
    p.add_argument('--step7',action='store_true',help='Run the locked multiscale and bilateral-cohort revision')
    p.add_argument('--step7-render-only',action='store_true',help='Refresh Step 7 tables/figures/reports from completed predictions')
    args=p.parse_args()
    cfg = load_config(args.config)
    if args.step7 or args.step7_render_only or Path(cfg['_config_path']).with_name('step7_protocol.yaml').exists():
        if args.step7_render_only:
            import json
            previous=Path(cfg['paths']['results'])/'reports/step_07_report.json'
            first=previous.with_name('step7_first_complete_execution.json')
            if previous.exists() and not first.exists():
                record=json.loads(previous.read_text())
                if record.get('status')=='complete' and not record.get('render_only'):write_json(record,first)
        from src.multiscale_benchmark import run_step7
        run_step7(cfg,render_only=args.step7_render_only)
        return
    from src.study import clinical_available
    if clinical_available(cfg):
        from src.evaluation import clinical_benchmarks,clinical_ablation,center_predictability
        from src.plotting import clinical_figures23
        with diagnostic_report(cfg,5) as report:
            meta,R,bins,pred,metrics,significance,signature=clinical_benchmarks(cfg)
            ablation=clinical_ablation(cfg,meta,R,bins,pred,signature)
            center_diagnostic=center_predictability(cfg,meta,R)
            clinical_figures23(meta,R,bins,pred,metrics,significance,ablation,center_diagnostic,cfg)
            report['statistics']=metrics[metrics.center.eq('AllCenters')].to_dict('records')
            report['combinations']=24
            report['analysed_patients']=len(meta)
            report['signature']=signature
            report['inference']='DeLong only within held-out centers with >=2 per class; Wilcoxon n=3 centers; 1000 fixed-score patient bootstraps'
            report['_markdown']='## Clinical benchmark\n\n24 combinations, 5-fold patient-isolated Random-CV and 3-fold LOCO; source-only DANN. Thresholds use nested 3-fold source-only OOF controls at 95% specificity. InFold_ComBat uses explicitly transductive unlabeled test-batch alignment.\n\n' + metrics[metrics.center.eq('AllCenters')][['mode','pipeline','model','auroc','macro_center_auroc','sensitivity95','specificity95','brier','ece']].to_markdown(index=False)+'\n\nSensitivity grid is descriptive across held-out studies and never selects the locked k=3/alpha=.25 configuration.\n\n'
        return
    with diagnostic_report(cfg, 4) as report:
        meta, X = load_dataset(cfg)
        signature = analysis_signature(cfg)
        predictions, unavailable = run_benchmarks(meta, X, cfg, signature)
        summary, centers, paired, mitigation = summarize_predictions(predictions, meta, cfg)
        dest = Path(cfg['paths']['results']) / 'tables'
        for frame, name in [(summary, 'cv_summary'), (centers, 'cv_per_center'), (paired, 'generalization_delta'),
                            (mitigation, 'mitigation_delta')]:
            write_parquet(frame, dest / f'{name}.parquet')
            frame.to_csv(dest / f'{name}.csv', index=False)
        fig2(summary, cfg)
        fig3(summary, predictions, cfg)
        write_json({'signature': signature, 'files': {name: sha256(dest / name) for name in
                    ('oof_predictions.parquet', 'cv_summary.parquet', 'cv_per_center.parquet',
                     'generalization_delta.parquet', 'mitigation_delta.parquet', 'selected_features.parquet', 'unavailable_comparisons.json')}},
                   dest / 'cv_provenance.json')
        report['statistics'] = summary.to_dict('records')
        report['status'] = 'complete_with_unavailable_comparisons' if unavailable else 'complete'
        report['unavailable_comparisons'] = unavailable
        report['requested_models'] = cfg['evaluation']['models']
        report['complete_comparisons'] = int(summary.status.eq('complete').sum())
        report['feature_missingness'] = X.isna().mean().to_dict()
        report['inferential_limit'] = 'OOF confidence intervals conditional on trained models; no model-retraining or external validation uncertainty'
        report['operating_threshold'] = '0.5 fixed; no target-center threshold fitting'


if __name__ == '__main__':
    main()
