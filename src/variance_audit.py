"""Descriptive embeddings and outcome-adjusted center variance analyses."""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from .debiasing import TrainMedianScaler
from .evaluation import bh_fdr


def design_matrix(meta, covariates):
    base = pd.DataFrame({'intercept': np.ones(len(meta)), 'label': meta.label.to_numpy()})
    missing = []
    for name in covariates:
        if name not in meta:
            missing.append(name)
            continue
        value = meta[name].reset_index(drop=True)
        if value.isna().any():
            raise ValueError(f'Audit covariate {name} has missing values; declare a complete-case policy')
        if pd.api.types.is_numeric_dtype(value):
            base[name] = value.astype(float)
        else:
            encoded = pd.get_dummies(value.astype(str), prefix=name, drop_first=True, dtype=float)
            base = pd.concat([base, encoded], axis=1)
    return base.to_numpy(dtype=float), missing, list(base)


def projection(Z):
    # SVD drops linearly dependent terms and avoids explicit n x n hat matrices.
    u, s, _ = np.linalg.svd(Z, full_matrices=False)
    rank = int(np.sum(s > (s[0] * 1e-10))) if len(s) else 0
    return u[:, :rank], rank


def partial_permanova(X, reduced, center_design, strata, permutations=999, seed=2026):
    """Partial Euclidean PERMANOVA, Freedman-Lane residual permutations within label.

    Tests the center contribution conditional on diagnosis and declared covariates.
    Exact randomization interpretation requires exchangeable residuals within strata.
    One row per independent patient is required by the caller. This is association,
    not proof that pre-analytic handling caused a center effect.
    """
    X = np.asarray(X, dtype=np.float64)
    q0, rank0 = projection(reduced)
    q1, rank1 = projection(np.column_stack([reduced, center_design]))
    df1, df2 = rank1 - rank0, len(X) - rank1
    if df1 <= 0 or df2 <= 0:
        raise ValueError('Center effect not identifiable after covariate adjustment')
    fitted = q0 @ (q0.T @ X)
    residual = X - fitted
    def statistic(values):
        ss0 = np.sum((values - q0 @ (q0.T @ values)) ** 2)
        ss1 = np.sum((values - q1 @ (q1.T @ values)) ** 2)
        effect = max(0., ss0 - ss1)
        return (effect / df1) / (ss1 / df2) if ss1 > 0 else np.inf, effect, ss1
    observed, effect, ss1 = statistic(X)
    rng = np.random.default_rng(seed)
    strata = np.asarray(strata)
    blocks = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    exceed = 0
    for _ in range(permutations):
        indices = np.arange(len(X))
        for block in blocks:
            indices[block] = rng.permutation(block)
        permuted, _, _ = statistic(fitted + residual[indices])
        exceed += permuted >= observed - 1e-12
    total = np.sum((X - X.mean(axis=0)) ** 2)
    return {'pseudo_F': float(observed), 'partial_R2': float(effect / (effect + ss1)),
            'total_R2_center_increment': float(effect / total) if total else np.nan,
            'p': (exceed + 1) / (permutations + 1), 'permutations': permutations,
            'df_center': df1, 'df_residual': df2,
            'inference': 'Freedman_Lane_within_diagnosis_exchangeable_patient_residuals',
            'ci': 'Not estimated for PERMANOVA effect; point estimate and permutation p only'}


def mixed_variance_partition(X, meta, reduced):
    """Center random-intercept LMM; disease/covariates fixed effects.

    center_fraction is variance(center)/(variance(center)+variance(residual)),
    conditional on fixed effects, NOT the repeated-measurement reliability ICC.
    Wald intervals for fixed disease effects; no spurious boundary variance p test.
    """
    import statsmodels.api as sm
    rows = []
    for name in X.columns:
        values = X[name].to_numpy(dtype=float)
        valid = np.isfinite(values)
        if np.sum(valid) <= reduced.shape[1] + 2 or np.std(values[valid]) < 1e-8:
            rows.append(dict(feature=name, status='insufficient_or_constant'))
            continue
        try:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter('always')
                fit = sm.MixedLM(values[valid], reduced[valid], groups=meta.center.to_numpy()[valid]).fit(
                    reml=True, method='lbfgs', maxiter=300, disp=False)
            between, residual = float(fit.cov_re[0, 0]), float(fit.scale)
            conf = np.asarray(fit.conf_int())
            rows.append(dict(feature=name, status='converged' if fit.converged else 'not_converged',
                             center_variance=between, residual_variance=residual,
                             center_fraction=between / (between + residual),
                             disease_beta=float(fit.fe_params[1]), disease_p=float(fit.pvalues[1]),
                             disease_ci_low=float(conf[1, 0]), disease_ci_high=float(conf[1, 1]),
                             warnings=' | '.join(str(w.message) for w in captured)))
        except (ValueError, np.linalg.LinAlgError) as exc:
            rows.append(dict(feature=name, status=f'not_estimable: {exc}'))
    table = pd.DataFrame(rows)
    if 'disease_p' in table:
        valid = table.status.eq('converged')
        table.loc[~valid, 'disease_p'] = np.nan
        table['disease_q'] = bh_fdr(table.disease_p)
    return table


def audit(meta, X, cfg):
    if meta.patient_id.duplicated().any():
        raise ValueError('Audit requires one selected sample per patient for independent permutations/LMM. '
                         'Use a prespecified audit cohort; do not silently drop repeat measurements.')
    if meta.center.nunique() < 3:
        raise ValueError('At least three independent centers/study proxies required')
    contingency = pd.crosstab(meta.center, meta.label).reindex(columns=[0, 1], fill_value=0)
    if not contingency.gt(0).all(axis=None):
        raise ValueError('Primary audit requires cases AND controls in each center; diagnosis-center confounding detected')
    meta = meta.copy()
    if 'accepted_fragments' in meta:
        meta['log_fragments'] = np.log1p(meta.accepted_fragments)
    reduced, missing, names = design_matrix(meta, cfg['audit']['covariates'])
    center = pd.get_dummies(meta.center, drop_first=True, dtype=float).to_numpy()
    scaler = TrainMedianScaler().fit(X)
    clean = scaler.transform(X)
    components = min(cfg['audit']['pca_components'], min(clean.shape) - 1)
    if components < 2:
        raise ValueError('Insufficient rank for PCA')
    pca = PCA(n_components=components, svd_solver='full', random_state=cfg['seed'])
    scores = pca.fit_transform(clean).astype('float32')
    embedding = meta[['sample_id', 'center', 'label']].copy()
    for i in range(components):
        embedding[f'PC{i + 1}'] = scores[:, i]
    umap_status = 'not_requested'
    if cfg['audit']['umap']:
        try:
            from umap import UMAP
            xy = UMAP(n_neighbors=min(15, len(meta) - 1), random_state=cfg['seed'], n_jobs=1).fit_transform(clean)
            embedding['UMAP1'], embedding['UMAP2'] = xy[:, 0].astype('float32'), xy[:, 1].astype('float32')
            umap_status = 'computed_descriptive_only'
        except ImportError:
            umap_status = 'unavailable_missing_umap_learn'
    permanence = partial_permanova(clean, reduced, center, meta.label.to_numpy(),
                                   cfg['audit']['permutations'], cfg['seed'])
    lmm = mixed_variance_partition(X, meta, reduced)
    details = dict(permanova=permanence, unmeasured_covariates=missing, fixed_effects=names,
                   pca_variance_explained=pca.explained_variance_ratio_.tolist(), umap_status=umap_status,
                   normalization='full-cohort scaling is descriptive audit only; never reused in CV',
                   causal_interpretation='center association may reflect protocol, composition, ancestry, or other confounding')
    return embedding, lmm, details


def clinical_design(meta):
    center=pd.get_dummies(meta.center,drop_first=True,dtype=float).to_numpy()
    disease=meta.label.to_numpy(dtype=float)[:,None]
    depth=meta.seq_depth.to_numpy(dtype=float)[:,None]
    return np.column_stack([np.ones(len(meta)),center,disease,depth]), center, disease, depth


def marginal_variance(X,meta):
    """Drop-one/type-III SS divided by total SS; unique fractions need not sum to 1."""
    from scipy.stats import f
    x=np.asarray(X,dtype=float)
    full,center,disease,depth=clinical_design(meta)
    q,rank=projection(full)
    total=np.sum((x-x.mean(axis=0))**2,axis=0)
    residual=np.sum((x-q@(q.T@x))**2,axis=0)
    values={}
    for name,reduced in [('batch',np.column_stack([np.ones(len(meta)),disease,depth])),
                          ('bio',np.column_stack([np.ones(len(meta)),center,depth]))]:
        qr,rr=projection(reduced)
        ss0=np.sum((x-qr@(qr.T@x))**2,axis=0)
        increment=np.maximum(ss0-residual,0)
        df=rank-rr
        effect=np.divide(increment,total,out=np.zeros_like(total),where=total>0)
        F=np.divide(increment/df,residual/(len(meta)-rank),out=np.zeros_like(total),where=residual>0)
        values['eta2_'+name]=effect
        values['F_'+name]=F
        values['p_'+name]=f.sf(F,df,len(meta)-rank)
    values['vulnerability_ratio']=values['eta2_batch']/(values['eta2_bio']+1e-8)
    values['vulnerable']=values['vulnerability_ratio']>5
    return pd.DataFrame(values)


def marginal_permanova(X,meta,term,permutations=999,seed=2026):
    full,center,disease,depth=clinical_design(meta)
    reduced=(np.column_stack([np.ones(len(meta)),disease,depth]) if term=='center'
             else np.column_stack([np.ones(len(meta)),center,depth]))
    x=np.asarray(X,dtype=float)
    q0,r0=projection(reduced);q1,r1=projection(full)
    fitted=q0@(q0.T@x);residual=x-fitted
    total=float(np.sum((x-x.mean(axis=0))**2))
    def statistic(z):
        ss0=np.sum(z*z)-np.sum((q0.T@z)**2)
        ss1=np.sum(z*z)-np.sum((q1.T@z)**2)
        increment=max(0.,ss0-ss1)
        return increment/(r1-r0)/(ss1/(len(x)-r1)),increment,ss1
    observed,increment,ss1=statistic(x)
    rng=np.random.default_rng(seed)
    strata=meta.label.to_numpy() if term=='center' else meta.center.to_numpy()
    blocks=[np.flatnonzero(strata==value) for value in np.unique(strata)]
    exceed=0
    for _ in range(permutations):
        order=np.arange(len(x))
        for block in blocks:order[block]=rng.permutation(block)
        permuted,_,_=statistic(fitted+residual[order])
        exceed+=permuted>=observed-1e-10
    return dict(term=term,pseudo_F=observed,R2=increment/total,partial_R2=increment/(increment+ss1),
                p=(exceed+1)/(permutations+1),permutations=permutations,
                adjustment='other_term_and_original_coverage_depth',
                permutation='Freedman_Lane_within_disease' if term=='center' else 'Freedman_Lane_within_center')


def clinical_audit(cfg):
    from pathlib import Path
    from scipy.stats import kruskal
    from threadpoolctl import threadpool_limits
    from .study import load_clinical,export_table,record_outputs
    from .debiasing import healthy_center_components
    from .common import write_parquet,write_json
    from .plotting import clinical_figure1
    from umap import UMAP
    meta,ratio,lengths,bins=load_clinical(cfg,patient_level=True)
    families={'raw_ratio':ratio,'median_centered_ratio':ratio-np.median(ratio,axis=1,keepdims=True),
              'length_frequency':lengths}
    embeddings=[];tests=[];partitions=[]
    with threadpool_limits(limits=1):
        for family,x in families.items():
            clean=TrainMedianScaler().fit(x).transform(x)
            pc=PCA(n_components=min(10,min(clean.shape)-1),svd_solver='full').fit_transform(clean)
            xy=UMAP(n_neighbors=15,min_dist=.1,random_state=cfg['seed'],n_jobs=1).fit_transform(clean)
            frame=meta[['sample_id','patient_id','center','label']].copy()
            frame['family']=family;frame['PC1']=pc[:,0];frame['PC2']=pc[:,1]
            frame['UMAP1']=xy[:,0];frame['UMAP2']=xy[:,1];embeddings.append(frame)
            for term in ('center','disease'):
                tests.append(dict(family=family,**marginal_permanova(clean,meta,term,cfg['study']['permutations'],cfg['seed'])))
            if family!='length_frequency':
                table=pd.concat([bins.reset_index(drop=True),marginal_variance(x,meta)],axis=1)
                table['family']=family;partitions.append(table)
    perma=pd.DataFrame(tests);perma['q']=bh_fdr(perma.p)
    partition=pd.concat(partitions,ignore_index=True)
    for term in ['batch','bio']:partition['q_'+term]=bh_fdr(partition['p_'+term])
    healthy=meta.label.eq(0).to_numpy()
    counts=meta[healthy].groupby('center').size().to_dict()
    drift=[]
    for j,bin_id in enumerate(bins.bin_id):
        groups=[ratio[healthy & meta.center.eq(c).to_numpy(),j] for c in sorted(counts)]
        H,p=kruskal(*groups)
        stable=[g for g in groups if len(g)>=5]
        sh,sp=kruskal(*stable) if len(stable)>=2 else (np.nan,np.nan)
        drift.append(dict(bin_id=bin_id,H=H,p=p,H_min5=sh,p_min5=sp))
    drift=pd.DataFrame(drift)
    drift['q']=bh_fdr(drift.p);drift['q_min5']=bh_fdr(drift.p_min5)
    icc,between,within,_=healthy_center_components(ratio,meta.label,meta.center)
    drift['center_icc1']=icc;drift['between_center_var']=between;drift['within_center_var']=within
    for frame,name in [(perma,'Audit_PERMANOVA'),(partition,'Bin_Variance_Partition'),(drift,'Healthy_Drift')]:
        export_table(frame,cfg,name)
    combined=pd.concat([partition.assign(record_type='bin_variance'),perma.assign(record_type='permanova'),
                        drift.assign(record_type='healthy_drift')],ignore_index=True)
    export_table(combined,cfg,'Table2_Variance_Partition_Summary')
    embedding=pd.concat(embeddings,ignore_index=True)
    write_parquet(embedding,Path(cfg['paths']['results'])/'tables/clinical_embeddings.parquet')
    write_parquet(meta,Path(cfg['paths']['results'])/'tables/analysis_participants.parquet')
    clinical_figure1(meta,ratio,lengths,bins,embedding,perma,partition,drift,cfg)
    details={'independent_patients':len(meta),'original_specimens':len(pd.read_parquet(cfg['paths']['metadata'])),
             'participant_rule':cfg['study']['participant_policy'],'healthy_counts':counts,
             'healthy_significant_windows':int(drift.q.lt(.01).sum()),
             'healthy_significant_fraction':float(drift.q.lt(.01).mean()),
             'healthy_min5_significant_fraction':float(drift.q_min5.lt(.01).mean()),
             'vulnerable_raw_windows':int(partition.query('family=="raw_ratio"').vulnerable.sum()),
             'permanova':perma.to_dict('records'),
             'cautions':['Study proxies, not physical center IDs','One healthy Snyder participant; KW asymptotic inference is fragile; min5 sensitivity supplied',
                         'Center ICC(1) is healthy-person clustering, not measurement repeatability','Type III unique eta-squared components need not sum to one',
                         'Descriptive full-cohort embeddings are never used to fit CV models']}
    write_json(details,Path(cfg['paths']['results'])/'tables/clinical_audit_summary.json')
    record_outputs(cfg,'clinical_audit',[Path(cfg['paths']['results'])/'tables/clinical_audit_summary.json',
                   Path(cfg['paths']['results'])/'tables/Table2_Variance_Partition_Summary.csv'])
    return details
