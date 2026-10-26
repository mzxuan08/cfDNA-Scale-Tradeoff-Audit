"""Patient-purged splits, paired OOF evaluation, and explicitly qualified inference."""
from __future__ import annotations

import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.model_selection import StratifiedGroupKFold

from .common import fingerprint, write_json, write_parquet
from .debiasing import FoldPreprocessor, ICCDesignError
from .models import make_model, TorchClassifier


def threshold_bootstrap(source_controls, target_scores, target_labels, repeats=1000, seed=2026):
    """Conditional calibration uncertainty with source-only control resampling.

    The outer model, calibration scores and target patients stay fixed. This is
    not a model-refitting interval or a target-calibrated operating threshold.
    """
    controls=np.asarray(source_controls,dtype=float)
    p=np.asarray(target_scores,dtype=float);y=np.asarray(target_labels)
    if not len(controls) or not np.isfinite(controls).all():
        raise ValueError('Finite held-out source-control scores required')
    rng=np.random.default_rng(seed)
    draws=rng.choice(controls,size=(repeats,len(controls)),replace=True)
    thresholds=np.quantile(draws,.95,axis=1,method='higher')
    h=p[y==0];c=p[y==1]
    specificity=(h[None,:]<=thresholds[:,None]).mean(axis=1) if len(h) else np.full(repeats,np.nan)
    sensitivity=(c[None,:]>thresholds[:,None]).mean(axis=1) if len(c) else np.full(repeats,np.nan)
    return thresholds,specificity,sensitivity


def paired_auc_difference_ci(labels, scores_a, scores_b, repeats=1000, seed=2026):
    """Paired diagnosis-stratified patient bootstrap for two fixed target rules."""
    y=np.asarray(labels);a=np.asarray(scores_a);b=np.asarray(scores_b)
    if min(np.sum(y==0),np.sum(y==1))<2:raise ValueError('Two patients per class required')
    rng=np.random.default_rng(seed);weights=np.zeros((repeats,len(y)),dtype=np.int32)
    for label in [0,1]:
        ids=np.flatnonzero(y==label)
        weights[:,ids]=rng.multinomial(len(ids),np.full(len(ids),1/len(ids)),size=repeats)
    delta=weighted_clinical_metrics(y,a,.5,weights)['auroc']-weighted_clinical_metrics(y,b,.5,weights)['auroc']
    return {'delta_auroc':float(roc_auc_score(y,a)-roc_auc_score(y,b)),
            'delta_low':float(np.quantile(delta,.025)),'delta_high':float(np.quantile(delta,.975)),
            'bootstrap_replicates':repeats,'bootstrap_scope':'paired_patient_diagnosis_stratified_fixed_scores'}


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=np.float64)
    out = np.full(p.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    if ((p[valid] < 0) | (p[valid] > 1)).any():
        raise ValueError('p values must be in [0,1]')
    order = valid[np.argsort(p[valid])]
    m = len(order)
    if m:
        ranked = p[order] * m / np.arange(1, m + 1)
        out[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return out


def delong_test(y, prediction_a, prediction_b):
    """Paired DeLong U-statistic placement covariance (ties receive half credit).

    Only valid for independent subjects evaluated by two FIXED scoring rules.
    Shared CV training makes OOF DeLong p values anti-conservative; the benchmark
    uses cluster bootstrap descriptively instead, and never reports this as its test.
    """
    y = np.asarray(y)
    predictions = np.asarray([prediction_a, prediction_b], dtype=np.float64)
    if not np.isin(y, [0, 1]).all() or predictions.shape != (2, len(y)) or not np.isfinite(predictions).all():
        raise ValueError('Invalid DeLong inputs')
    positives, negatives = predictions[:, y == 1], predictions[:, y == 0]
    m, n = positives.shape[1], negatives.shape[1]
    if min(m, n) < 2:
        raise ValueError('DeLong needs >=2 positives and negatives')
    v_pos, v_neg, aucs = [], [], []
    for pos, neg in zip(positives, negatives):
        compare = (pos[:, None] > neg).astype(float) + 0.5 * (pos[:, None] == neg)
        v_pos.append(compare.mean(axis=1))
        v_neg.append(compare.mean(axis=0))
        aucs.append(compare.mean())
    covariance = np.cov(v_pos, ddof=1) / m + np.cov(v_neg, ddof=1) / n
    variance = float(np.array([1, -1]) @ covariance @ np.array([1, -1]))
    delta = float(aucs[0] - aucs[1])
    p = float(2 * norm.sf(abs(delta) / np.sqrt(variance))) if variance > 1e-15 else (1.0 if abs(delta) < 1e-15 else 0.0)
    return dict(auc_a=float(aucs[0]), auc_b=float(aucs[1]), delta=delta, variance=max(variance, 0), p=p)


def ece(y, p, bins=10):
    y, p = np.asarray(y), np.asarray(p)
    if len(y) == 0 or not np.isfinite(p).all() or (p < 0).any() or (p > 1).any():
        raise ValueError('Invalid calibration inputs')
    membership = np.minimum((p * bins).astype(int), bins - 1)
    return float(sum(np.mean(membership == b) * abs(y[membership == b].mean() - p[membership == b].mean())
                     for b in range(bins) if np.any(membership == b)))


def metric_values(y, p, bins=10):
    y, p = np.asarray(y), np.asarray(p)
    if not np.isin(y, [0, 1]).all():
        raise ValueError('Binary outcomes required')
    both = len(np.unique(y)) == 2
    return {'auroc': float(roc_auc_score(y, p)) if both else np.nan,
            'auprc': float(average_precision_score(y, p)) if both else np.nan,
            'brier': float(brier_score_loss(y, p)), 'ece': ece(y, p, bins),
            'sensitivity_0.5': float(np.mean(p[y == 1] >= 0.5)) if np.any(y == 1) else np.nan,
            'specificity_0.5': float(np.mean(p[y == 0] < 0.5)) if np.any(y == 0) else np.nan}


def assert_split(meta, train, test, loco=False):
    if set(train) & set(test) or not len(train) or not len(test):
        raise ValueError('Empty/overlapping sample split')
    a, b = meta.iloc[train], meta.iloc[test]
    if set(a.patient_id) & set(b.patient_id):
        raise ValueError('Patient leakage')
    if set(a.sample_id) & set(b.sample_id):
        raise ValueError('Sample leakage')
    if loco and set(a.center) & set(b.center):
        raise ValueError('LOCO center leakage')
    if a.label.nunique() < 2:
        raise ValueError('Training fold lacks one outcome class')


def make_splits(meta, mode, n_splits=5, seed=2026):
    if mode == 'random':
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        folds = splitter.split(meta, meta.label, groups=meta.patient_id)
    elif mode == 'loco':
        folds = []
        for center in sorted(meta.center.unique()):
            test = np.flatnonzero(meta.center.eq(center))
            patients = set(meta.iloc[test].patient_id)
            train = np.flatnonzero(~meta.center.eq(center) & ~meta.patient_id.isin(patients))
            folds.append((train, test))
    else:
        raise ValueError('mode must be random or loco')
    result = []
    for train, test in folds:
        assert_split(meta, train, test, loco=mode == 'loco')
        result.append((np.asarray(train), np.asarray(test)))
    if mode == 'loco' and len(result) < 3:
        raise ValueError('At least 3 centers required for the multi-center LOCO benchmark')
    test_ids = np.concatenate([test for _, test in result])
    if len(test_ids) != len(meta) or len(np.unique(test_ids)) != len(meta):
        raise ValueError('Each sample must appear in exactly one test fold')
    return result


def bootstrap_ci(frame, metric='auroc', repeats=2000, seed=2026):
    """Hierarchical center then patient bootstrap, conditional on fixed OOF scores.

    This does not incorporate model-retraining uncertainty. Patient-center rows
    with shared patients across centers are rejected to prevent false independence.
    """
    if frame.patient_id.duplicated().any():
        raise ValueError('Bootstrap expects one row per patient, globally reconciled')
    rng = np.random.default_rng(seed)
    groups = [g for _, g in frame.groupby('center', sort=True)]
    if not groups:
        raise ValueError('Empty bootstrap input')
    values = []
    for _ in range(repeats):
        pieces = []
        for j in rng.integers(len(groups), size=len(groups)):
            group = groups[j]
            pieces.append(group.iloc[rng.integers(len(group), size=len(group))])
        draw = pd.concat(pieces, ignore_index=True)
        value = metric_values(draw.label, draw.probability)[metric]
        if np.isfinite(value):
            values.append(value)
    if len(values) < max(20, repeats // 2):
        raise ValueError('Insufficient class-informative bootstrap draws')
    return dict(low=float(np.quantile(values, .025)), high=float(np.quantile(values, .975)),
                valid_replicates=len(values), resampling='center_then_patient_conditional_on_fixed_OOF')


def paired_bootstrap_delta(frame, repeats=2000, seed=2026):
    """Paired Random-minus-LOCO AUROC difference; fixed-score descriptive CI."""
    if frame.patient_id.duplicated().any():
        raise ValueError('Paired comparison requires one independent row per patient')
    rng = np.random.default_rng(seed)
    groups = [g for _, g in frame.groupby('center', sort=True)]
    values = []
    for _ in range(repeats):
        pieces = []
        for j in rng.integers(len(groups), size=len(groups)):
            g = groups[j]
            pieces.append(g.iloc[rng.integers(len(g), size=len(g))])
        draw = pd.concat(pieces)
        if draw.label.nunique() == 2:
            values.append(roc_auc_score(draw.label, draw.random_probability) - roc_auc_score(draw.label, draw.loco_probability))
    if len(values) < max(20, repeats // 2):
        raise ValueError('Insufficient paired bootstrap draws')
    return {'delta_random_minus_loco': float(roc_auc_score(frame.label, frame.random_probability) - roc_auc_score(frame.label, frame.loco_probability)),
            'ci_low': float(np.quantile(values, .025)), 'ci_high': float(np.quantile(values, .975)),
            'valid_replicates': len(values), 'inference': 'descriptive_fixed_OOF_no_retraining'}


def fit_predict(name, preprocessor, cfg, Xtrain, ytrain, centers_train, Xtest, centers_test):
    train = preprocessor.transform(Xtrain, centers_train)
    test = preprocessor.transform(Xtest, centers_test)
    model = make_model(name, cfg)
    try:
        if isinstance(model, TorchClassifier):
            model.fit(train, ytrain, center_train=centers_train)
        else:
            model.fit(train, ytrain)
        if name == 'xgboost':
            device = json.loads(model.get_booster().save_config())['learner']['generic_param']['device']
            if not device.startswith('cuda'):
                raise RuntimeError(f'XGBoost silently fell back to {device}; CUDA result rejected')
        return model.predict_proba(test)[:, 1], getattr(model, 'peak_vram_bytes_', None)
    finally:
        if isinstance(model, TorchClassifier):
            model.close()


def run_benchmarks(meta, X, cfg, checkpoint_signature):
    e, d = cfg['evaluation'], cfg['debiasing']
    results = Path(cfg['paths']['results'])
    failures, predictions, selections = [], [], []
    all_splits = {mode: make_splits(meta, mode, e['random_folds'], cfg['seed']) for mode in ('random', 'loco')}
    split_rows = []
    for mode, folds in all_splits.items():
        for fold, (train, test) in enumerate(folds):
            for role, indices in (('train', train), ('test', test)):
                split_rows.extend({'mode': mode, 'fold': fold, 'role': role,
                                   'sample_id': meta.iloc[i].sample_id,
                                   'patient_id': meta.iloc[i].patient_id,
                                   'center': meta.iloc[i].center} for i in indices)
    write_parquet(pd.DataFrame(split_rows), results / 'tables/split_manifest.parquet')
    for mode, folds in all_splits.items():
        for fold, (train, test) in enumerate(folds):
            a, b = meta.iloc[train], meta.iloc[test]
            for method in e['methods']:
                if method == 'combat' and (set(b.center) - set(a.center)):
                    failures.append(dict(mode=mode, fold=fold, method=method, model='all',
                                         reason='UNAVAILABLE: ComBat cannot fit unseen test-center parameters'))
                    continue
                p = FoldPreprocessor(method, list(X.columns), d['local_radius'], d['icc_threshold'], d['icc_min_subjects'])
                try:
                    p.fit(X.iloc[train], a.label, a.center, a.patient_id)
                except (ICCDesignError, ValueError, RuntimeError) as exc:
                    failures.append(dict(mode=mode, fold=fold, method=method, model='all', reason=str(exc)))
                    continue
                names = np.asarray(X.columns)[p.support_]
                for name in names:
                    selections.append(dict(mode=mode, fold=fold, method=method, feature=name))
                for model_name in e['models']:
                    key = fingerprint([checkpoint_signature, mode, fold, method, model_name])
                    path = results / 'checkpoints' / (key + '.parquet')
                    try:
                        if path.exists() and path.with_suffix('.json').exists():
                            from .common import sha256
                            stamp = json.loads(path.with_suffix('.json').read_text())
                            if stamp['sha256'] != sha256(path):
                                raise ValueError('CV checkpoint checksum mismatch')
                            pred = pd.read_parquet(path)
                        else:
                            prob, peak = fit_predict(model_name, p, cfg, X.iloc[train], a.label.to_numpy(),
                                                     a.center.to_numpy(), X.iloc[test], b.center.to_numpy())
                            if not np.isfinite(prob).all() or (prob < 0).any() or (prob > 1).any():
                                raise ValueError('Invalid classifier probabilities')
                            pred = b[['sample_id', 'patient_id', 'center', 'label']].copy()
                            pred['probability'] = prob.astype('float32')
                            pred['mode'], pred['fold'], pred['method'], pred['model'] = mode, fold, method, model_name
                            pred['peak_vram_bytes'] = peak
                            write_parquet(pred, path)
                            from .common import sha256
                            write_json({'sha256': sha256(path), 'signature': key}, path.with_suffix('.json'))
                        predictions.append(pred)
                    except (ImportError, RuntimeError, ValueError, MemoryError) as exc:
                        failures.append(dict(mode=mode, fold=fold, method=method, model=model_name,
                                             reason=f'{type(exc).__name__}: {exc}'))
    write_json(failures, results / 'tables/unavailable_comparisons.json')
    if selections:
        write_parquet(pd.DataFrame(selections), results / 'tables/selected_features.parquet')
    if not predictions:
        raise RuntimeError('No benchmark completed; see unavailable_comparisons.json')
    predictions = pd.concat(predictions, ignore_index=True)
    write_parquet(predictions, results / 'tables/oof_predictions.parquet')
    return predictions, failures


def summarize_predictions(predictions, meta, cfg):
    records, per_center, patients = [], [], []
    e = cfg['evaluation']
    for (mode, method, model), g in predictions.groupby(['mode', 'method', 'model']):
        complete = len(g) == len(meta) and set(g.sample_id) == set(meta.sample_id) and not g.sample_id.duplicated().any()
        if not complete:
            records.append(dict(mode=mode, method=method, model=model, status='incomplete_excluded_from_comparison', n=len(g)))
            continue
        # Cross-center repeated subjects cannot be assigned a unique bootstrap center.
        n_centers = g.groupby('patient_id').center.nunique()
        patient = g.groupby('patient_id', as_index=False).agg(center=('center', 'first'), label=('label', 'first'), probability=('probability', 'mean'))
        patient['mode'], patient['method'], patient['model'] = mode, method, model
        patients.append(patient)
        row = dict(mode=mode, method=method, model=model, status='complete', n=len(patient),
                   **metric_values(patient.label, patient.probability, e['ece_bins']))
        if n_centers.gt(1).any():
            row['ci_status'] = 'unavailable: repeated subjects cross bootstrap center blocks'
        elif patient.label.nunique() == 2:
            ci = bootstrap_ci(patient, repeats=e['bootstrap_replicates'], seed=cfg['seed'])
            row.update(auroc_ci_low=ci['low'], auroc_ci_high=ci['high'], ci_status=ci['resampling'])
        center_aucs = []
        for center, cg in g.groupby('center'):
            cg = cg.groupby('patient_id', as_index=False).agg(label=('label', 'first'), probability=('probability', 'mean'))
            metrics = metric_values(cg.label, cg.probability, e['ece_bins'])
            per_center.append(dict(mode=mode, method=method, model=model, center=center, n=len(cg), **metrics))
            if np.isfinite(metrics['auroc']):
                center_aucs.append(metrics['auroc'])
        row['macro_center_auroc'] = np.mean(center_aucs) if center_aucs else np.nan
        row['centers_with_estimable_auc'] = len(center_aucs)
        records.append(row)
    summary, center_table = pd.DataFrame(records), pd.DataFrame(per_center)
    pairs, mitigation = [], []
    patient_table = pd.concat(patients, ignore_index=True) if patients else pd.DataFrame()
    if not patient_table.empty and not meta.groupby('patient_id').center.nunique().gt(1).any():
        for (method, model), g in patient_table.groupby(['method', 'model']):
            wide = g.pivot(index=['patient_id', 'center', 'label'], columns='mode', values='probability').reset_index()
            if {'random', 'loco'} <= set(wide) and not wide[['random', 'loco']].isna().any().any():
                wide = wide.rename(columns={'random': 'random_probability', 'loco': 'loco_probability'})
                pairs.append(dict(method=method, model=model, **paired_bootstrap_delta(wide, e['bootstrap_replicates'], cfg['seed'])))
        for (mode, model), g in patient_table.groupby(['mode', 'model']):
            wide = g.pivot(index=['patient_id', 'center', 'label'], columns='method', values='probability').reset_index()
            if 'raw' not in wide:
                continue
            methods = [m for m in e['methods'] if m != 'raw' and m in wide]
            for method in methods:
                if wide[['raw', method]].isna().any().any():
                    continue
                paired = wide[['patient_id', 'center', 'label', method, 'raw']].rename(
                    columns={method: 'random_probability', 'raw': 'loco_probability'})
                record = paired_bootstrap_delta(paired, e['bootstrap_replicates'], cfg['seed'])
                record['delta_method_minus_raw'] = record.pop('delta_random_minus_loco')
                mitigation.append(dict(mode=mode, method=method, model=model, **record))
    return summary, center_table, pd.DataFrame(pairs), pd.DataFrame(mitigation)


def weighted_clinical_metrics(y,p,threshold,weights):
    """Exact weighted AUROC/AP (including ties), and vectorized bootstrap metrics."""
    y=np.asarray(y,dtype=float);p=np.asarray(p,dtype=float);threshold=np.asarray(threshold)
    w=np.atleast_2d(np.asarray(weights,dtype=float))
    order=np.argsort(-p,kind='stable');scores=p[order]
    starts=np.r_[0,np.flatnonzero(np.diff(scores)!=0)+1]
    pos=np.add.reduceat(w[:,order]*y[order],starts,axis=1)
    neg=np.add.reduceat(w[:,order]*(1-y[order]),starts,axis=1)
    tp,fp=np.cumsum(pos,axis=1),np.cumsum(neg,axis=1)
    npos,nneg=tp[:,-1],fp[:,-1]
    auc=np.sum((tp+np.column_stack([np.zeros(len(w)),tp[:,:-1]]))*.5*neg,axis=1)/(npos*nneg)
    ap=np.sum(np.divide(tp,tp+fp,out=np.zeros_like(tp),where=tp+fp>0)*pos,axis=1)/npos
    total=w.sum(axis=1);positive=p>threshold
    sensitivity=(w@(y*positive))/npos
    specificity=(w@((1-y)*(~positive)))/nneg
    brier=(w@((p-y)**2))/total
    calibration=np.zeros(len(w))
    membership=np.minimum((p*10).astype(int),9)
    for b in range(10):
        indicator=membership==b
        calibration+=np.abs(w@((y-p)*indicator))/total
    return dict(auroc=auc,auprc=ap,sensitivity95=sensitivity,specificity95=specificity,
                fpr95=1-specificity,brier=brier,ece=calibration)


def clinical_metric_ci(frame,repeats=1000,seed=2026):
    if frame.patient_id.duplicated().any():raise ValueError('Independent patients required for clinical bootstrap')
    y=frame.label.to_numpy();p=frame.probability.to_numpy();t=frame.threshold.to_numpy()
    if len(np.unique(y))<2:return {'n':len(frame),'status':'single_class_not_estimable'}
    rng=np.random.default_rng(seed)
    weights=np.zeros((repeats,len(frame)),dtype=np.int32)
    # Patient bootstrap stratified by study AND diagnosis; no simulated subjects.
    for _,group in frame.reset_index(drop=True).groupby(['center','label']):
        ids=group.index.to_numpy();weights[:,ids]=rng.multinomial(len(ids),np.full(len(ids),1/len(ids)),size=repeats)
    draws=weighted_clinical_metrics(y,p,t,weights)
    point=weighted_clinical_metrics(y,p,t,np.ones(len(frame)))
    out={'n':len(frame),'cancer':int(y.sum()),'healthy':int((y==0).sum()),'status':'complete',
         'ci_scope':'patient_bootstrap_stratified_study_diagnosis_fixed_fitted_scores',
         'bootstrap_replicates':repeats}
    for metric,value in point.items():
        out[metric]=float(value[0]);out[metric+'_low']=float(np.quantile(draws[metric],.025));out[metric+'_high']=float(np.quantile(draws[metric],.975))
    return out


def fit_clinical_model(name,cfg,train,y,centers,test):
    model=make_model(name,cfg)
    try:
        if isinstance(model,TorchClassifier):model.fit(train,y,center_train=centers)
        else:model.fit(train,y)
        if name=='xgboost':
            import xgboost as xgb
            device=json.loads(model.get_booster().save_config())['learner']['generic_param']['device']
            if not device.startswith('cuda'):raise RuntimeError('XGBoost CUDA fallback rejected')
            probability=model.get_booster().predict(xgb.DMatrix(test))
        else:probability=model.predict_proba(test)[:,1]
        importance=(np.abs(model.coef_[0]) if hasattr(model,'coef_') else
                    model.feature_importances_ if hasattr(model,'feature_importances_') else None)
        diagnostics={'peak_vram_bytes':getattr(model,'peak_vram_bytes_',None),
                     'max_iterations_reached':bool(hasattr(model,'n_iter_') and hasattr(model,'max_iter') and np.max(model.n_iter_)>=model.max_iter)}
        return probability.astype('float32'),importance,diagnostics
    finally:
        if isinstance(model,TorchClassifier):model.close()


def source_calibration_splits(meta,n_splits=3,seed=2026):
    """Rare study/diagnosis strata stay training-only, never scored in calibration.

    A singleton healthy source subject cannot simultaneously supply two-center
    ICC training and a held-out calibration observation. We do not invent another
    control or substitute in-sample predictions for its missing OOF score.
    """
    from sklearn.model_selection import StratifiedKFold
    strata=meta.center.astype(str)+'__'+meta.label.astype(str)
    counts=strata.map(strata.value_counts()).to_numpy()
    eligible=np.flatnonzero(counts>=n_splits);always=np.flatnonzero(counts<n_splits)
    splits=[]
    for train,test in StratifiedKFold(n_splits,shuffle=True,random_state=seed).split(eligible,strata.iloc[eligible]):
        a=np.r_[eligible[train],always];b=eligible[test]
        assert not set(meta.iloc[a].patient_id)&set(meta.iloc[b].patient_id)
        splits.append((a,b))
    return splits,always


def clinical_benchmarks(cfg):
    from .study import load_clinical,study_signature,export_table,PIPELINES,MODELS
    from .debiasing import ClinicalPreprocessor
    from .common import sha256
    from sklearn.model_selection import StratifiedKFold
    from threadpoolctl import threadpool_limits
    from scipy.stats import wilcoxon
    meta,R,_,bins=load_clinical(cfg,True)
    destination=Path(cfg['paths']['results']);tables=destination/'tables'
    cache=destination/'clinical_cv_checkpoints';cache.mkdir(parents=True,exist_ok=True)
    # Only learning code contributes to fit-cache identity; plot edits do not force retraining.
    code={name:sha256(Path(__file__).parent/name) for name in ['evaluation.py','debiasing.py','models.py','study.py']}
    signature=fingerprint([study_signature(cfg),code])
    predictions=[];importances=[];diagnostics=[];split_rows=[]
    def preprocess(pipeline,indices,radius=None,alpha=None):
        return ClinicalPreprocessor(pipeline,bins,radius or cfg['study']['radius'],
                  cfg['study']['exclusion_fraction'] if alpha is None else alpha,cfg['study']['relative_epsilon']).fit(
                      R[indices],meta.label.to_numpy()[indices],center_train=meta.center.to_numpy()[indices])
    for mode in ['random','loco']:
        splits=make_splits(meta,mode,5,cfg['seed'])
        for fold,(train,test) in enumerate(splits):
            for role,ids in [('train',train),('test',test)]:
                split_rows.extend(dict(mode=mode,fold=fold,role=role,patient_id=meta.iloc[i].patient_id,
                                       sample_id=meta.iloc[i].sample_id,center=meta.iloc[i].center) for i in ids)
            a,b=meta.iloc[train],meta.iloc[test]
            inner_splits,rare_calibration=source_calibration_splits(a,cfg['study']['inner_calibration_folds'],cfg['seed'])
            for pipeline in PIPELINES:
                missing=[]
                for name in MODELS:
                    key=fingerprint([signature,mode,fold,pipeline,name]);path=cache/(key+'.parquet')
                    if path.exists() and path.with_suffix('.json').exists():
                        stamp=json.loads(path.with_suffix('.json').read_text())
                        if stamp['sha256']!=sha256(path):raise ValueError('CV checkpoint corrupted')
                        predictions.append(pd.read_parquet(path));diagnostics.append(stamp['diagnostics'])
                        if path.with_suffix('.importance.parquet').exists():importances.append(pd.read_parquet(path.with_suffix('.importance.parquet')))
                    else:missing.append((name,path))
                if not missing:continue
                pre=preprocess(pipeline,train)
                xt=pre.transform(R[train],a.center.to_numpy());xv=pre.transform(R[test],b.center.to_numpy())
                inner_data=[]
                for it,iv in inner_splits:
                    ii,jj=train[it],train[iv];assert_split(meta,ii,jj)
                    ip=preprocess(pipeline,ii)
                    inner_data.append((it,iv,ip.transform(R[ii],meta.iloc[ii].center.to_numpy()),
                                        ip.transform(R[jj],meta.iloc[jj].center.to_numpy())))
                with threadpool_limits(limits=1):
                    for name,path in missing:
                        calibration=np.full(len(a),np.nan,dtype='float32')
                        for it,iv,ix,iy in inner_data:
                            calibration[iv],_,_=fit_clinical_model(name,cfg,ix,a.label.to_numpy()[it],a.center.to_numpy()[it],iy)
                        eligible=np.ones(len(a),dtype=bool);eligible[rare_calibration]=False
                        if not np.isfinite(calibration[eligible]).all():raise ValueError('Incomplete eligible inner OOF calibration')
                        control=eligible & (a.label.to_numpy()==0)
                        threshold=float(np.quantile(calibration[control],.95,method='higher'))
                        prob,importance,diag=fit_clinical_model(name,cfg,xt,a.label.to_numpy(),a.center.to_numpy(),xv)
                        pred=b[['sample_id','patient_id','center','label']].copy()
                        pred['probability']=prob;pred['threshold']=threshold
                        pred['mode'],pred['fold'],pred['pipeline'],pred['model']=mode,fold,pipeline,name
                        pred['protocol']='transductive_unlabeled_target_batch' if pipeline=='InFold_ComBat' else 'inductive_no_target_population_fit'
                        write_parquet(pred,path)
                        info=dict(mode=mode,fold=fold,pipeline=pipeline,model=name,threshold=threshold,
                            inner_control_n=int(control.sum()),inner_training_specificity=float(np.mean(calibration[control]<=threshold)),
                            training_only_rare_stratum_patients=a.iloc[rare_calibration].patient_id.to_list(),
                            training_healthy_per_center=a[a.label.eq(0)].groupby('center').size().to_dict(),**diag)
                        write_json({'sha256':sha256(path),'diagnostics':info},path.with_suffix('.json'))
                        if importance is not None:
                            imp=pd.DataFrame({'bin_id':bins.loc[pre.support_,'bin_id'].to_numpy(),'importance':importance})
                            imp['mode'],imp['fold'],imp['pipeline'],imp['model']=mode,fold,pipeline,name
                            imp['rank']=imp.importance.rank(ascending=False,method='min')
                            write_parquet(imp,path.with_suffix('.importance.parquet'));importances.append(imp)
                        predictions.append(pred);diagnostics.append(info)
                        print(f'CV completed {mode} fold={fold} {pipeline} {name}; threshold={threshold:.5f}',flush=True)
                if pre.filter_ is not None:
                    selection=bins[['bin_id']].copy();selection['retained']=pre.support_;selection['healthy_center_icc']=pre.filter_.icc_
                    selection['mode'],selection['fold'],selection['pipeline']=mode,fold,pipeline
                    write_parquet(selection,cache/f'selection_{mode}_{fold}_{pipeline}.parquet')
    pred=pd.concat(predictions,ignore_index=True)
    export_table(pd.DataFrame(split_rows),cfg,'Clinical_CV_Splits')
    write_parquet(pred,tables/'Clinical_OOF_Predictions.parquet')
    write_json(diagnostics,tables/'Clinical_CV_Diagnostics.json')
    if importances:export_table(pd.concat(importances,ignore_index=True),cfg,'Clinical_Feature_Importance')
    rows=[]
    for (mode,pipeline,model),g in pred.groupby(['mode','pipeline','model']):
        if len(g)!=len(meta) or g.sample_id.duplicated().any():raise ValueError('Incomplete or duplicate clinical OOF predictions')
        for center,sub in [('AllCenters',g),*list(g.groupby('center'))]:
            rows.append(dict(mode=mode,pipeline=pipeline,model=model,center=center,
                        protocol=sub.protocol.iloc[0],**clinical_metric_ci(sub,cfg['study']['bootstrap_replicates'],cfg['seed'])))
    metrics=pd.DataFrame(rows)
    macro=metrics[metrics.center.ne('AllCenters')].groupby(['mode','pipeline','model']).auroc.mean().to_dict()
    metrics['macro_center_auroc']=[macro[(mode,pipeline,model)] for mode,pipeline,model in zip(metrics['mode'],metrics.pipeline,metrics.model)]
    export_table(metrics,cfg,'Table3_Benchmark_Metrics_Full')
    significance=[]
    loco=pred[pred['mode'].eq('loco')]
    for (model,center),g in loco.groupby(['model','center']):
        wide=g.pivot(index=['patient_id','label'],columns='pipeline',values='probability').reset_index()
        for comparator in ['Raw_Ratio','Global_Zscore']:
            try:r=delong_test(wide.label,wide.Proposed_LSN_ICC,wide[comparator]);status='computed_heldout_center_fixed_models'
            except ValueError as exc:r={'p':np.nan,'delta':np.nan};status=str(exc)
            significance.append(dict(test='paired_DeLong',model=model,center=center,comparator=comparator,status=status,**r))
    for model,g in metrics[metrics['mode'].eq('loco')&metrics.center.ne('AllCenters')].groupby('model'):
        wide=g.pivot(index='center',columns='pipeline',values='auroc')
        for comparator in ['Raw_Ratio','Global_Zscore']:
            delta=(wide.Proposed_LSN_ICC-wide[comparator]).dropna().to_numpy()
            p=1. if np.allclose(delta,0) else float(wilcoxon(delta,zero_method='pratt',method='exact').pvalue)
            significance.append(dict(test='paired_Wilcoxon_centers',model=model,center='AllCenters',comparator=comparator,
                                      p=p,delta=float(np.median(delta)),status=f'{len(delta)}_paired_centers'))
    significance=pd.DataFrame(significance);significance['q']=bh_fdr(significance.p)
    export_table(significance,cfg,'Clinical_Paired_Significance')
    return meta,R,bins,pred,metrics,significance,signature


def clinical_ablation(cfg,meta,R,bins,main_pred,signature):
    from .debiasing import ClinicalPreprocessor
    from .study import export_table
    from sklearn.metrics import roc_auc_score
    from .common import sha256
    cache=Path(cfg['paths']['results'])/'clinical_cv_checkpoints'
    settings=[('LSN_only',3,0.),('ICC_only',3,.25)]+[
        ('Proposed_LSN_ICC',k,alpha) for k in cfg['study']['radius_grid'] for alpha in cfg['study']['exclusion_grid']]
    rows=[];prediction_frames=[]
    for pipeline,k,alpha in settings:
        scores=[]
        for fold,(train,test) in enumerate(make_splits(meta,'loco')):
            key=fingerprint([signature,'ablation',pipeline,k,alpha,fold]);path=cache/(key+'.parquet')
            if path.exists():g=pd.read_parquet(path)
            elif pipeline=='Proposed_LSN_ICC' and k==cfg['study']['radius'] and alpha==cfg['study']['exclusion_fraction']:
                g=main_pred[(main_pred['mode']=='loco')&(main_pred.pipeline==pipeline)&(main_pred.model=='elasticnet')&(main_pred.fold==fold)].copy()
                write_parquet(g,path)
            else:
                pre=ClinicalPreprocessor(pipeline,bins,k,alpha,cfg['study']['relative_epsilon']).fit(
                    R[train],meta.iloc[train].label.to_numpy(),center_train=meta.iloc[train].center.to_numpy())
                a=pre.transform(R[train],meta.iloc[train].center.to_numpy());b=pre.transform(R[test],meta.iloc[test].center.to_numpy())
                prob,_,_=fit_clinical_model(cfg['study']['sensitivity_model'],cfg,a,meta.iloc[train].label.to_numpy(),meta.iloc[train].center.to_numpy(),b)
                g=meta.iloc[test][['sample_id','patient_id','center','label']].copy();g['probability']=prob
                write_parquet(g,path)
            score=float(roc_auc_score(g.label,g.probability));scores.append(score)
            g=g.copy();g['ablation'],g['radius'],g['alpha'],g['fold']=pipeline,k,alpha,fold;prediction_frames.append(g)
            rows.append(dict(pipeline=pipeline,radius=k,alpha=alpha,fold=fold,center=g.center.iloc[0],
                             model=cfg['study']['sensitivity_model'],auroc=score,n=len(g),selection='prespecified_no_outer_score_parameter_selection'))
        rows.append(dict(pipeline=pipeline,radius=k,alpha=alpha,fold=-1,center='MacroCenters',model=cfg['study']['sensitivity_model'],
                         auroc=float(np.mean(scores)),n=len(meta),selection='prespecified_no_outer_score_parameter_selection'))
        print(f'Ablation completed {pipeline} k={k} alpha={alpha}',flush=True)
    table=pd.DataFrame(rows);export_table(table,cfg,'Table4_Ablation_Study')
    write_parquet(pd.concat(prediction_frames,ignore_index=True),Path(cfg['paths']['results'])/'tables/Ablation_Predictions.parquet')
    return table


def center_predictability(cfg,meta,R):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import balanced_accuracy_score
    from .study import export_table
    rows=[]
    for scope,mask in [('all_participants',np.ones(len(meta),dtype=bool)),
                       ('healthy_centers_with_at_least_5_controls',meta.label.eq(0).to_numpy())]:
        ids=np.flatnonzero(mask)
        counts=meta.iloc[ids].groupby('center').size()
        ids=ids[meta.iloc[ids].center.isin(counts[counts>=5].index)]
        target=meta.iloc[ids].center.to_numpy()
        for fold,(train,test) in enumerate(StratifiedKFold(5,shuffle=True,random_state=cfg['seed']).split(R[ids],target)):
            model=RandomForestClassifier(n_estimators=300,max_depth=8,n_jobs=cfg['evaluation']['n_jobs'],random_state=cfg['seed'])
            model.fit(R[ids[train]],target[train])
            score=balanced_accuracy_score(target[test],model.predict(R[ids[test]]))
            rows.append(dict(scope=scope,fold=fold,balanced_accuracy=score,n_centers=len(counts[counts>=5]),
                             n=len(ids),interpretation='center_predictability_not_causal_proof_of_overfitting'))
    frame=pd.DataFrame(rows);export_table(frame,cfg,'Center_Predictability_Diagnostics')
    return frame
