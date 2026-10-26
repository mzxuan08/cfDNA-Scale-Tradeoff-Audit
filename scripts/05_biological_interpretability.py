from _bootstrap import parser
import json
from pathlib import Path
import pandas as pd

from src.common import load_config, diagnostic_report, load_dataset, analysis_signature, sha256, write_json, write_parquet
from src.plotting import fig4


def main():
    p = parser('Exploratory disease associations and train-fold selection stability')
    p.add_argument('--annotations-only',action='store_true')
    args=p.parse_args()
    cfg = load_config(args.config)
    from src.study import clinical_available
    if clinical_available(cfg):
        from src.interpretability import build_annotations,biology_analysis
        if args.annotations_only:build_annotations(cfg)
        else:biology_analysis(cfg)
        return
    with diagnostic_report(cfg, 5) as report:
        meta, X = load_dataset(cfg)
        table = Path(cfg['paths']['results']) / 'tables'
        provenance = json.loads((table / 'audit_provenance.json').read_text())
        if provenance['signature'] != analysis_signature(cfg):
            raise ValueError('Audit is stale; rerun step 03')
        for name, digest in provenance['files'].items():
            if sha256(table / name) != digest:
                raise ValueError(f'Audit table checksum mismatch: {name}')
        lmm = pd.read_parquet(table / 'variance_partition.parquet')
        if 'disease_q' not in lmm:
            raise ValueError('No estimable LMM disease association; cannot create biological claims')
        valid = lmm[lmm.status.eq('converged')].copy()
        if valid.empty:
            raise ValueError('No converged LMM associations')
        write_parquet(valid, table / 'exploratory_associations.parquet')
        fig4(valid, cfg)
        selected = table / 'selected_features.parquet'
        if selected.exists():
            cp = json.loads((table / 'cv_provenance.json').read_text())
            if cp['signature'] != analysis_signature(cfg) or cp['files']['selected_features.parquet'] != sha256(selected):
                raise ValueError('CV selection table is stale; rerun step 04')
            s = pd.read_parquet(selected)
            selection = s.groupby(['mode', 'method', 'feature'], as_index=False).agg(selected_folds=('fold', 'nunique'))
            denominator = s.groupby(['mode', 'method'])['fold'].nunique().to_dict()
            selection['eligible_folds'] = [denominator[(m, t)] for m, t in zip(selection['mode'], selection['method'])]
            selection['selection_fraction'] = selection.selected_folds / selection.eligible_folds
            write_parquet(selection, table / 'feature_selection_stability.parquet')
        report['statistics'] = {'converged_associations': len(valid), 'BH_q_below_0.05': int(valid.disease_q.lt(.05).sum()),
                                 'multiplicity': 'BH within all converged LMM disease-feature tests',
                                 'effect_and_CI_table': str(table / 'exploratory_associations.parquet')}
        report['interpretation'] = 'Exploratory adjusted associations, not causal/pathway validation; external genomic annotation not inferred'
        report['feature_missingness'] = X.isna().mean().to_dict()
        write_json({'signature': analysis_signature(cfg), 'files': {'exploratory_associations.parquet': sha256(table / 'exploratory_associations.parquet')}}, table / 'biology_provenance.json')


if __name__ == '__main__':
    main()
