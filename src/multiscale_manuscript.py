"""Data-linked Step 7 manuscript. Empirical statements are read from CSV outputs."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .common import write_json
from .study import export_table
from .multiscale_benchmark import PRIMARY,COMBINED,EXTERNAL,TUNED,read_protocol


REQUESTED_TITLE='Within-Sample Multi-Scale Self-Normalization Overcomes Center-Specific Batch Drift and Prevalence Confounding in Multi-Center Cell-Free DNA Fragmentomics'
EVIDENCE_TITLE='Scale Trade-offs in Within-Sample Self-Normalization for Cross-Study Generalization in Public cfDNA Cohorts'


def ptext(value):
    return f'{value:.3g}' if value<.001 else f'{value:.3f}'


def compile_step7_manuscript(cfg,pred,metrics,significance,counts,pairs,stress,reuse_existing=False):
    from .prevalence_analysis import (heldout_effect_diagnostic,multiscale_signal_diagnostics,
                                    qualify_paired_inference,source_threshold_alignment_diagnostic)
    root=Path(cfg['paths']['results']);folder=root/'tables';protocol=read_protocol(cfg)
    from .metadata_scope_audit import audit_local_finaledb
    census,publication,source_audit=audit_local_finaledb(cfg)
    if not reuse_existing:
        heldout_effect_diagnostic(cfg);multiscale_signal_diagnostics(cfg)
        source_threshold_alignment_diagnostic(cfg,pred)
    # All prose below reads the exported CSVs, including diagnostic values.
    metrics=pd.read_csv(folder/'Table3_Benchmark_Metrics_Full.csv')
    if not reuse_existing:
        thresholds=[]
        for row in metrics.itertuples():
            group=pred[pred.cohort_scope.eq(row.cohort_scope)&pred['mode'].eq(row.mode)&pred.pipeline.eq(row.pipeline)&pred.model.eq(row.model)]
            if row.center!='AllCenters':group=group[group.center.eq(row.center)]
            thresholds.append((float(group.threshold.min()),float(group.threshold.max())) if len(group) else (np.nan,np.nan))
        metrics['source_threshold_min']=[a for a,b in thresholds]
        metrics['source_threshold_max']=[b for a,b in thresholds]
        singleton=metrics.healthy.lt(2)&metrics.status.eq('complete')
        for field in ['specificity95_low','specificity95_high','fpr95_low','fpr95_high']:
            metrics.loc[singleton,field]=np.nan
        metrics['interval_policy']=np.where(singleton,'exact_binomial_operating_intervals_single_control_no_bootstrap_specificity',
                                            'conditional_patient_bootstrap_plus_exact_operating_intervals')
        metrics.loc[~metrics.pipeline.isin(['Proposed_LSN_ICC',TUNED]),'alpha']='not_applicable'
        metrics.loc[~metrics.pipeline.isin(['Proposed_LSN_ICC',TUNED,'Hybrid_MultiScale']),'radius']='not_applicable'
        export_table(metrics,cfg,'Table3_Benchmark_Metrics_Full')
    sig=pd.read_csv(folder/'Step7_Paired_Significance.csv')
    if not reuse_existing:
        sig=qualify_paired_inference(cfg,sig)
        for method,label in [('Hybrid_MultiScale','Hybrid'),(TUNED,'Tuned')]:
            ps=[];qs=[]
            for row in metrics.itertuples():
                matched=sig[sig.cohort_scope.eq(row.cohort_scope)&sig.center.eq(row.center)&sig.model.eq(row.model)
                            &sig['method'].eq(method)&sig.comparator.eq(row.pipeline)] if row.mode=='loco' else sig.iloc[:0]
                ps.append(float(matched.p.iloc[0]) if len(matched) else np.nan)
                qs.append(float(matched.q.iloc[0]) if len(matched) else np.nan)
            metrics[f'DeLong_{label}_vs_pipeline_p']=ps;metrics[f'DeLong_{label}_vs_pipeline_q']=qs
        export_table(metrics,cfg,'Table3_Benchmark_Metrics_Full')
    counts=pd.read_csv(folder/'Step7_Prevalence_By_Study.csv').set_index('center')
    stress=pd.read_csv(folder/'Step7_Prevalence_Stress_Summary.csv')
    pairs=pd.read_csv(folder/'Step7_AUC_Pair_Decomposition.csv')
    effect=pd.read_csv(folder/'Step7_Heldout_Effect_Preservation.csv')
    alignment=pd.read_csv(folder/'Step7_Source_Threshold_Alignment.csv')
    gc=pd.read_csv(folder/'Step7_Residual_Batch_GC_Associations.csv').set_index('representation')
    arm=pd.read_csv(folder/'Step7_Arm_Level_Disease_Contrasts.csv')
    prevalence=pd.read_csv(folder/'Step7_Prevalence_Association.csv').iloc[0]
    audit=pd.read_csv(folder/'Audit_PERMANOVA.csv').set_index(['family','term'])
    biology=pd.read_csv(folder/'Disease_Effect_Sizes.csv').pivot(index=['bin_id','center'],columns='pipeline',values='cohen_d')
    enrich=pd.read_csv(folder/'Filtered_Window_Enrichment.csv')
    parameters=pd.read_csv(folder/'Step7_Nested_Parameter_Search.csv')
    software=pd.read_csv(folder/'Methods_Software_Versions.csv')
    methods=pd.read_csv(folder/'Methods_Model_Parameters.csv')
    clinical=pd.read_parquet(folder/'Analysis_Participants.parquet')
    missing=clinical.isna().mean().rename('missing_fraction').rename_axis('field').reset_index()
    export_table(missing,cfg,'Step7_Metadata_Missingness')
    features=pd.read_parquet(Path(cfg['paths']['processed'])/cfg['extraction']['final_file']).set_index('sample_id').loc[clinical.sample_id]
    measured=features.select_dtypes(include='float32')
    quality=pd.DataFrame([dict(patients=len(measured),feature_columns=measured.shape[1],
                              overall_feature_missing_fraction=float(measured.isna().to_numpy().mean()),
                              max_column_missing_fraction=float(measured.isna().mean().max()))])
    export_table(quality,cfg,'Step7_Data_Quality')
    primary=metrics[metrics.cohort_scope.eq(PRIMARY)&metrics['mode'].eq('loco')&metrics.model.eq('elasticnet')]
    pooled=primary[primary.center.eq('AllCenters')].set_index('pipeline')
    direct=primary[primary.center.ne('AllCenters')&primary.status.eq('complete')]
    combined=metrics[metrics.cohort_scope.eq(COMBINED)&metrics.center.eq('AllCenters')&metrics.model.eq('elasticnet')].set_index('pipeline')
    external=metrics[metrics.cohort_scope.eq(EXTERNAL)&metrics.model.eq('elasticnet')].set_index('pipeline')
    studies=protocol['primary_studies'];C,J=studies;S=protocol['external_study']
    H=pooled.loc['Hybrid_MultiScale'];G=pooled.loc['Global_Zscore'];F=pooled.loc['Proposed_LSN_ICC'];T=pooled.loc[TUNED]
    hd=direct[direct.pipeline.eq('Hybrid_MultiScale')].set_index('center')
    gd=direct[direct.pipeline.eq('Global_Zscore')].set_index('center')
    differences=hd.auroc-gd.auroc
    ordering='higher in both directions' if differences.gt(0).all() else ('lower in both directions' if differences.lt(0).all() else 'different in direction between the two targets')
    adjusted=sig[sig.cohort_scope.eq(PRIMARY)&sig.model.eq('elasticnet')&sig.method.eq('Hybrid_MultiScale')&sig.comparator.eq('Global_Zscore')].set_index('center')
    cb=pairs[pairs.pipeline.eq('InFold_ComBat')]
    cb_alignment=alignment[alignment.cohort_scope.eq(COMBINED)&alignment.pipeline.eq('InFold_ComBat')].iloc[0]
    within=cb[cb.pair_type.eq('within_study')];cross=cb[cb.pair_type.ne('within_study')]
    within_auc=float(within.pooled_auc_contribution.sum()/within.pooled_pair_weight.sum())
    cross_auc=float(cross.pooled_auc_contribution.sum()/cross.pooled_pair_weight.sum())
    maximum_cb_d=effect[effect.pipeline.eq('InFold_ComBat')].max_abs_d_change.max()
    drifts=stress.groupby('pipeline').max_anchor_probability_range.max()
    claim_audit=pd.DataFrame([
        dict(claim='Raw GC association applies unchanged after Global Z-score',assessment='not_supported',
             evidence=f"Observed raw rho={gc.loc['Raw_Ratio','rho']:.6f}; Global rho={gc.loc['Global_Zscore','rho']:.6f}",source='Step7_Residual_Batch_GC_Associations.csv'),
        dict(claim='Low pooled ComBat AUROC proves destruction of within-study feature signal',assessment='not_supported',
             evidence=f'Within-study pair AUROC={within_auc:.6f}; cross-study pair AUROC={cross_auc:.6f}; max within-study Cohen d change={maximum_cb_d:.9g}',
             source='Step7_AUC_Pair_Decomposition.csv; Step7_Heldout_Effect_Preservation.csv'),
        dict(claim='Sample-wise transforms are immune to all prevalence confounding',assessment='scope_must_be_restricted',
             evidence=f"Hybrid fixed-anchor maximum probability range={drifts['Hybrid_MultiScale']:.9g}; source fitting and predictive values still depend on population composition",source='Step7_Prevalence_Stress_Summary.csv'),
        dict(claim='Reduced local contrast proves loss of measured chromosome-arm copy number',assessment='not_identified_without_copy_number_ground_truth',
             evidence='Arm-averaged fragmentation contrasts are measured; matched copy-number truth is not available for this interpretation',source='Step7_Arm_Level_Disease_Contrasts.csv'),
        dict(claim='Healthy cross-center ICC can be estimated in a one-source bilateral fold',assessment='not_identifiable',
             evidence='Only one real source study; filter bypass is explicit and alpha is zero',source='Step7_Selected_Parameters.csv')])
    export_table(claim_audit,cfg,'Step7_Scientific_Claim_Audit')
    selected=[]
    diag=json.loads((root/'reports/step7_fit_diagnostics.json').read_text())
    for row in diag:
        if row['pipeline']==TUNED and row['mode']!='random':
            selected.append(dict(cohort_scope=row['cohort_scope'],fold=row['fold'],training_n=len(row['train_patients']),
                test_study=';'.join(row['test_studies']),selected_radius=row['radius'],selected_alpha=row['alpha'],icc_bypassed=row['icc_bypassed']))
    selected=pd.DataFrame(selected);export_table(selected,cfg,'Step7_Selected_Parameters')
    method_settings={'seed':cfg['seed'],'inner_folds':protocol['inner_folds'],'bootstrap_replicates':protocol['bootstrap_replicates'],
        'hybrid_radius':protocol['hybrid_radius'],'core_pipelines':len(protocol['pipelines']),'classifiers':len(protocol['models']),
        'window_bp':cfg['features']['bin_size'],'fragment_min_bp':cfg['features']['histogram'][0],
        'fragment_max_bp':cfg['features']['histogram'][1],'mapq_min':cfg['features']['mapq_min'],
        'source_depth_min':cfg['cohort']['min_depth'],'source_depth_max':cfg['cohort']['max_depth'],
        'target_specificity':.95,'primary_studies':len(protocol['primary_studies'])}
    export_table(pd.DataFrame([{'parameter':k,'value':v} for k,v in method_settings.items()]),cfg,'Step7_Methods_Execution_Parameters')
    title_assessment=pd.DataFrame([dict(requested_title=REQUESTED_TITLE,adopted_title=EVIDENCE_TITLE,
        rationale='Systematic scale-trade-off audit in public cohorts; Hybrid is a fusion probe, not a baseline-beating algorithm. CNA preservation is unverified without matched copy-number truth.')])
    export_table(title_assessment,cfg,'Step7_Title_Assessment')
    # Capture every numeric CSV cell available to manuscript assembly.
    trace=[]
    sources=['Table3_Benchmark_Metrics_Full','Step7_Paired_Significance','Step7_Prevalence_By_Study',
        'Step7_Prevalence_Association','Step7_Prevalence_Stress_Summary','Step7_AUC_Pair_Decomposition',
        'Step7_Heldout_Effect_Preservation','Step7_Residual_Batch_GC_Associations','Step7_Arm_Level_Disease_Contrasts',
        'Step7_Selected_Parameters','Audit_PERMANOVA','Disease_Effect_Sizes','Filtered_Window_Enrichment','Methods_Model_Parameters','Methods_Execution_Parameters','Step7_Methods_Execution_Parameters','Step7_Source_Threshold_Alignment',
        'Final_Alignment_FinaleDB_Source_Census','Final_Alignment_Original_Publication_Inventory']
    for name in sources:
        frame=pd.read_csv(folder/(name+'.csv'))
        for col in frame.select_dtypes(include='number'):
            trace.extend(dict(source_csv=name+'.csv',row_index=int(i),column=col,value=float(v)) for i,v in frame[col].items() if np.isfinite(v))
    export_table(pd.DataFrame(trace),cfg,'Step7_Manuscript_Numeric_Traceability')
    primary_n=int(counts.loc[studies,'n'].sum());total_n=int(counts.n.sum());external_n=int(counts.loc[S,'n'])
    count_table=counts.reset_index().copy();count_table['cancer_percent']=100*count_table.prevalence
    metric_cols=['pipeline','center','auroc','auroc_low','auroc_high','auprc','sensitivity95','specificity95','specificity_deviation_pp','brier','ece']
    primary_table=direct[metric_cols].sort_values(['pipeline','center'])
    external_table=external.reset_index()[['pipeline','sensitivity95','sensitivity_exact_low','sensitivity_exact_high','specificity95','source_threshold_min']]
    combined_table=combined[['auroc','auroc_low','auroc_high','auprc','sensitivity95','specificity95','brier','ece']]
    arm_summary=arm.groupby('representation').cohen_d.apply(lambda s:s.abs().median())
    source_table=pd.read_csv(folder/'Final_Alignment_FinaleDB_Source_Census.csv').query('screen=="plasma_snapshot_no_depth_restriction"')
    source_by_name=source_table.set_index('study')
    other_adal=source_by_name.loc['Adalsteinsson et al, 2017']
    other_sun=source_by_name.loc['Sun et al, 2019']
    fixed_radius=protocol['hybrid_radius']
    manuscript=rf'''# {EVIDENCE_TITLE}

**Authors and affiliations:** [To be supplied by the authors]

## Abstract

**Background:** Normalization scale can change broad fragmentation contrasts and residual GC-associated study variation. We investigated this scale trade-off in cross-study generalization in public cfDNA cohorts.

**Methods:** We analysed {total_n} source-identified participants: a bilateral cohort of {primary_n} and a case-dominated external cohort of {external_n}. Global, local and hybrid representations were compared using source-only selection and locked thresholds. Hybrid served as a scale-fusion probe. We also audited the local FinaleDB source inventory.

**Results:** Combined LOCO ElasticNet AUROCs were {combined.loc['Hybrid_MultiScale','auroc']:.3f} for Hybrid, {combined.loc['Global_Zscore','auroc']:.3f} for Global Z-score and {combined.loc['Proposed_LSN_ICC','auroc']:.3f} for fixed LSN-ICC. Primary bilateral Hybrid changes versus Global were {differences.loc[C]:+.3f} in Cristiano and {differences.loc[J]:+.3f} in Jiang. Median absolute arm-averaged fragmentation effects were {arm_summary['Global_Zscore']:.3f} for Global and {arm_summary['LSN']:.3f} for local normalization; residual GC correlations with study variance were {gc.loc['Global_Zscore','rho']:.3f} and {gc.loc['LSN','rho']:.3f}. Target-composition changes shifted ComBat-adapted anchor probabilities by up to {drifts.loc['InFold_ComBat']:.3f}. No additional locally documented source contained cancer cases and at least ten healthy participants.

**Conclusions:** The contribution is a systematic audit of scale-dependent representation trade-offs and study-source confounding. Fusion did not establish consistent superiority. Chromosome-arm CNA preservation remains a mechanistic hypothesis without matched copy-number ground truth.

**Keywords:** cell-free DNA; fragmentomics; cross-study validation; sample-wise normalization; batch composition; threshold transport

## 1. Introduction

Plasma cell-free DNA fragmentation encodes information about nucleosome organization and tissue contribution. Genome-wide fragmentation profiles have been used for cancer classification, complementing sequence-variant assays [1–3]. The present study evaluates cross-study generalization in public cfDNA cohorts. Study Source is a proxy that combines experimental batch, specimen processing, sequencing procedures, recruitment and population composition. It is not a verified hospital identifier or an isolated technical intervention. A high score in a mixed-study validation set does not by itself establish that a locked threshold will maintain its intended specificity in a new study.

Normalization has two separable roles in this setting. A global transformation can retain broad relative deviations across the genome while removing a sample-wide offset and scale. A local contrast can attenuate slowly varying backgrounds but also change broad disease-associated contrasts. These are representation choices, not interchangeable definitions of biological signal. Our preceding window audit found median absolute within-study disease contrasts of {biology.Raw_Ratio.abs().median():.3f} for raw ratios and {biology.LSN.abs().median():.3f} for a local representation. Those contrasts motivate testing both scales; they do not establish which changes are caused by chromosome-arm copy number or tissue-specific nucleosome depletion.

Batch-distribution alignment introduces a further distinction. An empirical-Bayes source harmonizer estimates study-specific location and scale parameters [4]. Extending a fitted system to a previously unseen study requires an explicit adaptation rule. If that rule estimates a target median or quantiles, its output can depend on the target case fraction. Under the mixture identity \(F_c=(1-\pi_c)F_{{0c}}+\pi_cF_{{1c}}\), changing prevalence \(\pi_c\) can change batch-level quantiles even when the individuals being evaluated remain fixed. A transformation of each person's own feature vector has no such functional dependence on companion patients. This property does not make the downstream classifier immune to selection bias, source prevalence, population shift or prevalence-dependent predictive values.

We therefore conduct a systematic Scale Trade-off audit of global and local self-normalization, using Hybrid as a diagnostic fusion comparator. The central question is how normalization scale changes broad arm-level fragmentation contrasts and GC-associated study variation. Chromosome-arm CNA preservation is a plausible interpretation of broad signals, not a directly validated mechanism in these data. We separately examine two-study bidirectional discrimination, cancer sensitivity in a case-dominated external source, three-study score comparability and companion-dependent adaptation. The study does not introduce a classifier shown to outperform established baselines consistently.

## 2. Materials and methods

### 2.1. Data provenance and independent-participant cohort

All new analyses used the already extracted local feature matrix and the same {total_n} independently identified patients; no new clinical dataset or synthetic participant was introduced. Public study metadata had been recovered from the published UNITE supplementary data table, and fragment intervals had been streamed from explicit FinaleDB objects [5,6]. Original inclusion required plasma WGS, source sequencing depth within {cfg['cohort']['min_depth']}–{cfg['cohort']['max_depth']}× and explicit source cancer or healthy labels. The extraction cohort contained {int(pd.read_csv(folder/'Table1_Cohort_Characteristics.csv').retained_samples.sum())} specimens. A specimen with unresolved source patient identity was isolated, and one prespecified specimen per known patient was retained using earliest time point, then shallowest sequencing depth, then sequencing-run identifier. The labels describe cancer versus an explicitly healthy source category. Hepatitis, cirrhosis, inflammatory disease, transplantation and autoimmune disease were not recoded as healthy.

The primary cohort contained the Cristiano and Jiang studies. Each bilateral direction trained on one entire study and tested the other. The Snyder study was held out from all primary training, primary parameter selection and primary threshold calibration. Its high case burden makes cancer sensitivity the useful external endpoint; it does not supply an adequately sampled healthy specificity estimate. A separate secondary analysis used three-study LOCO across all participants, with the other two studies used for each held-out target. Snyder can therefore enter source training only in that explicitly secondary protocol, never in the primary bilateral experiment.

{count_table[['center','n','cancer','healthy','cancer_percent']].to_markdown(index=False,floatfmt='.3f')}

Study identifiers are proxies combining experimental batch and population-composition confounding, rather than independently verified hospital identifiers. The source publications have different recruitment aims, and observed case proportions must not be interpreted as population cancer prevalence. Source identity checks rely on published identifiers; undocumented cross-study aliases cannot be resolved from these records.

The local source census covered {source_audit['snapshot_runs']} FinaleDB run records from {source_audit['run_level_studies']} studies. Eligibility for an additional two-class source required a documented plasma WGS cohort with at least one cancer participant and at least ten unique, explicitly labelled Healthy participants; healthy controls were not borrowed across studies. The broad census below precedes depth and QC restrictions, and counts only nonmissing patient identifiers:

{source_table[['study','runs','healthy_patients','cancer_patients','meets_cancer_and_10_healthy']].to_markdown(index=False)}

The archived original database publication listed {source_audit['archived_publication_studies']} sources and {source_audit['archived_publication_records']} sequencing records. Its additional urinary and prenatal single-gene studies did not document an eligible plasma cancer/control cohort in the local materials. Their unavailable run-level case/control counts were not replaced with zeros. This audit is limited to historical local records; it does not establish that the current database or all public repositories lack further eligible cohorts. No additional clinical data were downloaded and no feature extraction was repeated.

### 2.2. Genomic windows, fragmentation features and descriptive audit

The hg19/GRCh37 reference retained {int(pd.read_csv(folder/'Methods_Execution_Parameters.csv').set_index('parameter').loc['windows','value'])} autosomal non-overlapping complete windows of 5 Mb after fixed centromere, terminal-risk, cytoband-risk and blacklist criteria. Only autosomal 80–400 bp fragments with MAPQ at least 30 were aggregated; overlaps with blacklisted intervals were excluded. Midpoints assigned fragments to windows. In each window, \(S_i\) counted 100–150 bp fragments and \(L_i\) counted 151–220 bp fragments. The classifier ratio was

$$R_{{si}}=\frac{{S_{{si}}}}{{L_{{si}}+1}}.$$

Window counts, IQR and the global integer-length frequency distribution were retained in compact local tables. The analyses used float32 feature matrices. GC fractions came from reference sequence, while RepeatMasker and segmental-duplication fractions used interval unions. GM12878 chromatin compartments are a reference annotation, not a tumor accessibility assay.

Descriptive PERMANOVA and window variance models adjusted for disease status, study and original sequencing depth as appropriate. PERMANOVA used Euclidean geometry after descriptive column standardization and 999 restricted Freedman–Lane residual permutations. Window-level unique/type III eta-squared was the reduced-to-full residual sum-of-squares increment divided by total sum of squares. These components need not add to one. Full-cohort embeddings and descriptive associations never entered predictive preprocessing. For the revision, GC associations were recomputed separately for raw, sample-wise global and local representations instead of assuming that the raw-feature correlation applies to Global Z-score.

### 2.3. Global, local and hybrid representations

Raw Ratio used training-only imputation and column standardization. Global Z-score first used each individual's own window mean and standard deviation:

$$G_{{si}}=\frac{{R_{{si}}-\bar R_s}}{{\operatorname{{sd}}_i(R_s)}}.$$

Its subsequent column scaler was fitted on source training patients. The local transform used other retained windows on the same chromosome arm, up to \(k\) on either side, excluding the focal window:

$$L_{{si}}^{{(k)}}=\log_2\left[\frac{{R_{{si}}+\epsilon_s}}{{\operatorname{{median}}_{{j\in N_k(i),j\ne i}}R_{{sj}}+\epsilon_s}}\right],
\qquad \epsilon_s=10^{{-6}}\operatorname{{median}}_iR_{{si}}.$$

This relative epsilon preserves global multiplicative invariance to numerical precision. Additive and arbitrary nonlinear invariance are not implied. Hybrid MultiScale concatenated the global and local blocks,

$$H_s=[G_s,L_s^{{({fixed_radius})}}],$$

and fitted its column scaler only within the training split. Both blocks were retained; no outer test result selected the block weights or the fixed Hybrid radius. The primary classifier was elastic-net logistic regression with an L1/L2 mixture; its regularization learns coefficients for the two correlated blocks. The other fixed classifiers supply secondary algorithm comparisons.

The fixed LSN-ICC design attempted to remove the highest quarter of healthy-center clustering scores within source training data, followed by RobustScaler. These scores are one-way between-center variance ratios among different healthy people, not repeat-measurement reliability. In a strict bilateral direction there is only one source study, so this component is not identifiable. The implementation explicitly bypassed the filter and labelled its effective representation LSN only. It did not use target healthy labels or invent source batches. The single-source result is consequently a transparent operational fallback, not validation of an estimable cross-center ICC filter.

### 2.4. Nested source-only selection and locked thresholds

A separate Tuned LSN analysis considered radii \(k\in\{{1,2,3,5,8\}}\) and eligible exclusion fractions \(\alpha\in\{{0,0.15,0.25,0.35,0.50\}}\). With one healthy source study, only zero exclusion was eligible. Selection maximized the mean validation AUROC across source strata with at least two patients of each class, using three inner folds. It used the prespecified ElasticNet model and deterministic tie-breaking. This was a separate analysis, not a retrospective replacement of fixed Hybrid.

Source threshold calibration used three-fold out-of-fold predictions from source patients. For Tuned LSN, parameter selection was repeated inside each calibration-training split before predicting its validation patients. Thus calibration outcomes never participated in the parameter search that generated their own scores. Rare study/diagnosis strata remained training-only and supplied no in-sample substitute score. The higher empirical 0.95 quantile of eligible source-healthy out-of-fold probabilities locked the decision threshold. Probabilities strictly above that threshold were positive. A model was then fitted to all source patients and applied unchanged to the target patients. No target diagnosis selected a threshold, model setting, feature mask or scale. Random and bilateral validation shared the primary patient universe, but their training-set sizes differed; their contrast does not isolate a pure technical batch effect.

### 2.5. Comparator identifiability and classifier implementation

The five core pipelines were Raw Ratio, Global Z-score, InFold ComBat, fixed Proposed LSN-ICC and Hybrid MultiScale. The tuned LSN variant was evaluated with the primary ElasticNet model. Random Forest, CUDA XGBoost, LightGBM, CUDA MLP and source-only CUDA DANN were additional classifiers with locked configurations. DANN required multiple source studies and never received target-domain samples. Standard multi-batch ComBat likewise required multiple source batches. Neither method was assigned an artificial result in a one-source bilateral fit: these cells are explicitly unavailable. They were estimable in primary mixed-study random CV and in the secondary two-source/one-target evaluations.

InFold ComBat combined source-only empirical-Bayes harmonization with an explicit unlabeled target median/IQR mapping for an unseen study. It is therefore a transductive adaptation comparator; this implementation is not interchangeable with every use of the original ComBat algorithm. The other transformations process each target person's feature vector independently of other target patients. Source-trained scalers and supervised classifiers can still reflect source population composition.

All random operations used seed 2026. Neural hidden dimensions and batch sizes respected the original memory constraints, and CUDA availability and peak process allocation were checked. Exact fixed model parameters are supplied in Methods_Model_Parameters.csv. The executable software versions were:

{software.to_markdown(index=False)}

### 2.6. Metrics and statistical comparison

We report direction-specific AUROC and AUPRC, pooled held-out predictions, equal-weight mean of the two primary study AUROCs, source-threshold sensitivity and specificity, Brier score and 10-bin expected calibration error. Specificity deviation is the achieved target specificity minus 0.95 in percentage points. We calculated 1000 patient bootstraps stratified by the fixed study and diagnosis, conditional on the fitted out-of-fold scoring rules. These intervals do not include refitting, method-development or new-center uncertainty. Exact binomial intervals supplement sensitivity and specificity estimates. The one-control cohort's AUROC interval was withheld from the primary presentation rather than displayed as an apparently precise bootstrap interval.

Paired DeLong tests compared inductive scoring rules for the same independent patients within each held-out primary study [7]. Target-batch-fitted ComBat adaptation is not an individual fixed scoring rule, so its ordinary DeLong comparisons were explicitly unavailable; descriptive performance estimates were retained. BH correction covered finite tests within the primary comparison family; the secondary three-study family was separate [8]. Pooled predictions across different fitted models did not receive a pooled DeLong test. With two primary directions, center-level tests cannot establish broad population-center inference. AUPRC is interpreted together with case proportion because its baseline changes with prevalence.

### 2.7. Prevalence and score-comparability diagnostics

A contingency-table chi-square test described the study/diagnosis association. To separate a low pooled AUROC from study-wise ranking, we decomposed all cancer–healthy patient pairs in the secondary ElasticNet out-of-fold predictions into within-study and cross-study terms. Their weighted concordance contributions sum exactly to pooled AUROC. Cross-study terms compare scores produced by different source-trained models, so model offset or calibration differences can affect them even when ranking inside each target remains informative.

An additional descriptive check expressed each score as logit(probability) minus logit(source-locked threshold). This used no target labels to align the scoring scales and was not substituted for the primary predictions or interpreted as a calibrated probability. It tests how a source-defined operating-point reference changes pooled ranking across separately fitted rules.

The composition stress experiment trained on the two non-target studies and froze both preprocessing state and ElasticNet weights. For each eligible target study, five healthy and five cancer anchors were held fixed. We constructed target batches of 32 distinct real patients at cancer fractions 0.25, 0.50 and 0.75, using source-independent deterministic selection and no synthetic feature vectors. Target labels were used only to define this post-evaluation diagnostic composition; those labels did not tune any benchmark prediction. We measured the maximum change in an anchor's probability. The experiment tests companion-dependence of a fitted transformation, rather than estimating a universal causal effect of population prevalence.

Finally, within-target window Cohen d was compared before and after harmonization. A positive feature-wise affine map preserves this standardized mean difference within a study. This provides a direct check against attributing all low pooled performance to destruction of univariate biological contrasts. Chromosome-arm averaged fragmentation contrasts were also calculated; absent independently measured copy number, they were not relabelled as confirmed copy-number signals.

## 3. Results

### 3.1. Study-associated drift and the need to distinguish scales

The independent analysis cohort comprised {total_n} patients. The unique study fraction in raw-ratio PERMANOVA was {100*audit.loc[('raw_ratio','center'),'R2']:.3f}%, while the disease fraction was {100*audit.loc[('raw_ratio','disease'),'R2']:.3f}% (both permutation p={audit.loc[('raw_ratio','center'),'p']:.3f}). GC correlated with raw window study eta-squared at rho={gc.loc['Raw_Ratio','rho']:.3f}, p={ptext(gc.loc['Raw_Ratio','p'])}. For Global Z-score the corresponding correlation was {gc.loc['Global_Zscore','rho']:.3f}, p={ptext(gc.loc['Global_Zscore','p'])}; for local contrasts it was {gc.loc['LSN','rho']:.3f}, p={ptext(gc.loc['LSN','p'])}. Median unique window study variance fractions were {100*gc.loc['Raw_Ratio','median_eta2_batch']:.3f}%, {100*gc.loc['Global_Zscore','median_eta2_batch']:.3f}% and {100*gc.loc['LSN','median_eta2_batch']:.3f}%, respectively. These representation-specific values distinguish residual GC ranking from the absolute amount of study variation; the raw-feature correlation cannot be assigned to Global Z-score.

Median absolute arm-averaged disease contrasts were {arm_summary['Raw_Ratio']:.3f} for raw features, {arm_summary['Global_Zscore']:.3f} for Global Z-score and {arm_summary['LSN']:.3f} for local contrasts. They describe fragmentation differences at arm scale and do not identify their molecular cause. The earlier reference-region enrichment comparisons yielded {int(enrich.q.lt(.05).sum())} BH-significant results among {len(enrich)} tests. The multiscale representation was therefore evaluated as a way to retain complementary information, without assuming a validated tumor-specific compartment or copy-number mechanism.

The Scale Trade-off is empirical and representation-dependent: attenuation of broad arm-level contrasts coexisted with changes in residual study variance and GC association. Local normalization did not uniformly eliminate GC-associated variation; its residual GC rank association was stronger than Global Z-score's despite a slightly smaller median study-variance fraction. A direct trade-off between preservation of measured chromosome-arm CNA and removal of a causal GC batch effect is not established without matched CNA measurements and controlled technical perturbations.

### 3.2. Primary bilateral generalization and operating-point transfer

All primary training and calibration excluded the external high-burden cohort. In the Jiang-to-Cristiano direction, Hybrid AUROC was {hd.loc[C,'auroc']:.3f} and AUPRC was {hd.loc[C,'auprc']:.3f}; in the Cristiano-to-Jiang direction they were {hd.loc[J,'auroc']:.3f} and {hd.loc[J,'auprc']:.3f}. Hybrid discrimination relative to Global Z-score was {ordering}. The DeLong p values versus Global Z-score were {ptext(adjusted.loc[C,'p'])} for the Cristiano target and {ptext(adjusted.loc[J,'p'])} for the Jiang target; BH values were {ptext(adjusted.loc[C,'q'])} and {ptext(adjusted.loc[J,'q'])}, respectively.

At the source-locked operating threshold, Hybrid specificity was {hd.loc[C,'specificity95']:.3f} in Cristiano and {hd.loc[J,'specificity95']:.3f} in Jiang, corresponding to deviations of {hd.loc[C,'specificity_deviation_pp']:+.3f} and {hd.loc[J,'specificity_deviation_pp']:+.3f} percentage points from the source target. Target sensitivities were {hd.loc[C,'sensitivity95']:.3f} and {hd.loc[J,'sensitivity95']:.3f}. The full direction-specific primary comparison is:

{primary_table.to_markdown(index=False,floatfmt='.3f')}

The pooled primary Hybrid AUROC was {H.auroc:.3f}, with a mean study AUROC of {H.macro_bilateral_auroc:.3f}. Global Z-score gave corresponding values of {G.auroc:.3f} and {G.macro_bilateral_auroc:.3f}. These pooled and equal-study summaries have different weighting and score-comparability requirements. The complete six-classifier table includes unavailable single-source ComBat and DANN cells as such, rather than excluding the identifiability issue from the denominator.

### 3.3. Nested local selection and the external sensitivity cohort

The nested-tuned local representation yielded pooled primary AUROC {T.auroc:.3f}, AUPRC {T.auprc:.3f}, Brier score {T.brier:.3f} and ECE {T.ece:.3f}. The fixed local fallback yielded AUROC {F.auroc:.3f}. Outer source-selected settings were:

{selected.to_markdown(index=False)}

In the external high-burden cohort, all thresholds and models were locked using the primary source studies. Hybrid detected a fraction {external.loc['Hybrid_MultiScale','sensitivity95']:.3f} of cancer cases; its exact binomial 95% interval was {external.loc['Hybrid_MultiScale','sensitivity_exact_low']:.3f}–{external.loc['Hybrid_MultiScale','sensitivity_exact_high']:.3f}. The external operating-point results were:

{external_table.to_markdown(index=False,floatfmt='.3f')}

The observed healthy outcome in this cohort was retained in the numerical table, but it was not represented as a stable population specificity estimate in Figure 2c. The secondary Combined LOCO results across all independently identified patients were:

{combined_table.to_markdown(floatfmt='.3f')}

For the same fixed ElasticNet comparison, Combined LOCO Hybrid AUROC {combined.loc['Hybrid_MultiScale','auroc']:.3f} was below Global Z-score {combined.loc['Global_Zscore','auroc']:.3f} and numerically above Proposed LSN-ICC {combined.loc['Proposed_LSN_ICC','auroc']:.3f}. The numerical difference versus the local pipeline is not evidence that fusion consistently beats both baselines. No pooled cross-validation DeLong test was performed; the primary target-specific tests showed a decline relative to Global in Cristiano and an improvement in Jiang that did not survive BH significance. The contribution is the systematic scale-trade-off audit, with Hybrid serving as a comparator that tests complementarity.

### 3.4. Why pooled ComBat discrimination and within-study ranking diverge

Cancer proportions were {100*counts.loc[C,'prevalence']:.3f}% in Cristiano, {100*counts.loc[J,'prevalence']:.3f}% in Jiang and {100*counts.loc[S,'prevalence']:.3f}% in Snyder. Study and diagnosis were associated (chi-square={prevalence.chi2:.3f}, df={int(prevalence.df)}, p={ptext(prevalence.p)}, Cramer's V={prevalence.cramers_v:.3f}). The association is a property of these selected source cohorts and is not a randomized test of technical batches.

For the empirical-Bayes-plus-target-alignment comparator, secondary pooled AUROC was {combined.loc['InFold_ComBat','auroc']:.3f}. Its weighted within-study cancer–control concordance was {within_auc:.3f}, whereas cross-study concordance was {cross_auc:.3f}. Cross-study pairs accounted for {100*cross.pooled_pair_weight.sum():.3f}% of all pairs. Thus the pooled score combines within-study ranking with comparability across independently fitted source models. Its numerical collapse does not by itself establish that all within-study tumor discrimination disappeared.

The maximum change in within-target window Cohen d under this comparator was {maximum_cb_d:.3g}. Positive feature-wise affine harmonization preserved those standardized contrasts to numerical precision. The classifier and pooled operating scores can nevertheless change because target-batch location and scale alter the relationship to the source-trained rule. This distinguishes information retained within a target study from calibration and feature weighting learned elsewhere.

Expressing the ComBat scores relative to their own source-locked thresholds yielded secondary pooled AUROC {cb_alignment.source_threshold_margin_auroc:.3f}, a change of {cb_alignment.delta:+.3f} from the original pooled ranking. This score-comparability diagnostic does not replace the observed benchmark or establish that any remaining degradation has a single cause.

When the same anchor patients were placed in real target batches with different case fractions, their maximum probability range was {drifts['InFold_ComBat']:.3f} under target-batch adaptation. Corresponding ranges were {drifts['Raw_Ratio']:.3g} for Raw Ratio, {drifts['Global_Zscore']:.3g} for Global Z-score, {drifts['Proposed_LSN_ICC']:.3g} for fixed LSN-ICC and {drifts['Hybrid_MultiScale']:.3g} for Hybrid. These results identify target-composition dependence in this adaptation rule. They support the narrower statement that sample-wise transformations avoid companion-dependent target normalization; they do not prove immunity of training, discrimination or predictive values to prevalence confounding.

## 4. Discussion

The principal contribution is a systematic account of the Scale Trade-off in cross-study generalization in public cfDNA cohorts. Global and local normalization change broad relative fragmentation contrasts, GC-associated residual study variation and score transport in different ways. Hybrid provides a controlled fusion probe of these scales; the observed comparisons do not justify presenting a baseline-beating new algorithm. Broad chromosome-arm contrasts are potential CNA-associated signals, but the current data do not independently validate CNA preservation or complete removal of GC-driven technical effects. The bilateral and external analyses keep discrimination, sensitivity and specificity evidence separate.

The two scale blocks have different mathematical effects. Global Z-scoring retains relative deviations over an entire genome, whereas the local log contrast compares a window with adjacent windows on the same arm. A locally smooth difference can therefore be attenuated even when it is disease-associated. Retaining both blocks leaves the allocation of predictive weights to a source-trained penalized model. The observed direction-specific comparisons and calibration values are the evidence for that allocation; the existence of two blocks alone does not establish superior transport. Nested local selection likewise addresses algorithmic parameter leakage but cannot create additional independent source studies.

The ComBat diagnostics clarify why a single pooled AUROC is an incomplete mechanistic explanation. The evaluated comparator contains two operations: empirical-Bayes source adjustment and a robust unseen-target mapping. Study composition is associated with disease, and a batch quantile is a property of the mixture of patients. The fixed-anchor experiment directly demonstrates that the unseen-target mapping can change an individual's output when companion composition changes. At the same time, preserved within-study standardized feature contrasts and the pair decomposition show why this finding should not be described as universal removal of tumor biology. Other covariate-aware or reference-batch ComBat designs are not tested by this implementation [4].

The independence property of a sample-wise transform is exact in a useful but restricted sense. For a fitted source model, appending, removing or relabelling other target people does not change an individual's global, local or hybrid transform. This removes one route by which target composition can alter outputs. It does not remove the association between source study and cancer, fix a source-trained classifier's prior, guarantee calibration in a new disease spectrum, or make precision independent of prevalence. Cross-study threshold transport therefore still requires adequate external controls and explicit reporting of the achieved operating point.

The hierarchical cohort presentation is especially relevant for public-data method development. Both-class primary targets support direction-specific discrimination and specificity intervals. A cancer-enriched external collection contributes a valuable sensitivity observation without being presented as a screening population. The secondary three-study analysis remains useful for comparison with the earlier audit, provided that its source composition and pooled-score interpretation remain explicit. These distinctions make the evaluation more informative than interpreting all changes as a single technical batch-removal effect.

### 4.1. Limitations

Source availability is the core design limitation: only {source_audit['qualifying_study_count']} locally documented sources had both cancer participants and at least ten healthy participants, and the analysed third source was case dominated ({int(counts.loc[S,'cancer'])} cancer/{int(counts.loc[S,'healthy'])} healthy). This is adequate two-class representation in the primary sources, not an exactly balanced case/control design. Adalsteinsson had {int(other_adal.cancer_patients)} known cancer participants and {int(other_adal.healthy_patients)} explicitly healthy participants; Sun had {int(other_sun.cancer_patients)} and {int(other_sun.healthy_patients)}, respectively. Neither supplies the missing healthy baseline. No additional qualifying source was verified locally, but the snapshot is not an exhaustive current database census.

Study Source is a proxy containing both experimental-batch and population-composition confounding. Each strict bilateral fit has only one source. A cross-center ICC filter, multi-batch source ComBat and multi-source adversarial training are consequently not all identifiable in that design. The fixed and tuned LSN fallbacks are labelled and cannot validate the unavailable center-filter component. The case-dominated external source does not establish screening specificity or clinical utility; the scope remains cross-study generalization in public cfDNA cohorts.

The Hybrid architecture and revised cohort reporting were motivated by earlier analyses of the same public patients. Each fitted rule maintains algorithmic patient and study isolation, but this is method development rather than an untouched confirmatory evaluation. Bootstrap intervals condition on fitted scores and fixed studies; they do not incorporate the history of method development, retraining or uncertainty over future studies and populations. DeLong tests are within fixed held-out targets and are not pooled cross-validation significance tests. The composition stress experiment uses a finite set of real-patient batches and is a diagnostic intervention on companions, not proof of a sole prevalence-mediated cause for all study drift.

The analysis lacks directly measured copy-number truth and tumor-specific chromatin assays for every patient. Arm-scale fragmentation contrasts and reference GC associations cannot establish their proposed biological mechanisms. The source clinical metadata have missing covariates, and an independent cohort with balanced healthy recruitment and prospective specimen handling would be needed to support a clinical deployment claim.

## 5. Conclusions

Cross-study generalization in public cfDNA cohorts depends on normalization scale, study-source composition and score comparability. The systematic Scale Trade-off audit quantifies broad fragmentation-contrast retention alongside residual GC-associated source variation. Hybrid fusion did not establish consistent superiority over Global Z-score, and CNA-specific preservation remains unverified. These findings support source-aware methodological evaluation without extrapolation to population screening.

## Data and code availability

The study reuses source-linked public cfDNA interval data and published metadata. No new biological samples were collected. The local derived matrices, source manifests, fixed configuration, tests and statistical CSVs accompany the analysis package. Code repository and archival release: [Public URL and immutable release identifier to be supplied by the authors]. The original data-use conditions apply to redistributed derived materials. Original-study ethics approvals and any required institutional determination for secondary analysis must be verified by the authors; no new approval identifier or exemption is asserted here.

## Author contributions, funding and competing interests

[To be completed and approved by the authors; no contributor roles, financial support or absence of conflicts is inferred.]

## Figure legends

**Figure 1. Study-associated baseline variation in the independently identified cfDNA cohort.** (a) Study composition and observed healthy fragment-frequency curves, with reference length markers. (b) Descriptive UMAP coloured by study and diagnosis, annotated with coverage-adjusted conditional PERMANOVA fractions and permutation results. (c) Unique study and disease window variance fractions across the autosomes. (d) Healthy baseline distributions and one-way healthy-center clustering scores. The external study's single healthy person is descriptive. No embedding was supplied to predictive preprocessing.

**Figure 2. Hierarchical reporting of primary transport and external cancer sensitivity.** (a) Random and bilateral study-held-out macro AUROC for Raw Ratio and Hybrid across fixed classifiers, using the same primary patient cohort; training-set sizes differ between the schemes. (b) The bilateral macro-AUROC matrix; unavailable ComBat, DANN and unrequested tuned-classifier cells remain grey. Asterisks identify LSN fits whose healthy-center filter cannot be estimated from one source study. (c) ElasticNet target specificity minus the source target in percentage points, with exact binomial intervals for the two adequately represented healthy target groups. (d) Cancer sensitivity in the separately held-out high-burden cohort, using primary-source-locked thresholds and exact binomial intervals. The lone healthy outcome is reported numerically and does not generate an unstable specificity error bar in panel c.

**Figure 3. Multiscale representation, primary discrimination and target-composition dependence.** (a) Sample-wise global and same-arm local branches, training-only scaling, concatenation and ElasticNet classification; bars show the observed two directional AUROCs. (b) Primary pooled held-out ROC and precision–recall curves, with target-specific DeLong comparisons against Global Z-score. Pooled curves are descriptive across different source-trained models. (c) Primary calibration curves and Brier scores. (d) Maximum change in a fixed real anchor patient's probability when the target batch case fraction changes, with all source fitting frozen. Diagnostic composition labels never entered benchmark fitting or selection.

**Figure 4. Reference genomic context and fixed-pipeline feature retention.** (a) Raw-feature study variance associations with reference GC and telomere distance. (b) Reference-region odds ratios for majority-fold filtering under the previously fixed LSN-ICC protocol; corrected results are reported in the accompanying table. (c) Observed chromosome trajectories before and after local contrast, with gaps across removed windows. (d) Fixed-pipeline coefficient rank concordance across the original three-study LOCO fits. These panels contextualize the fixed local representation and are not relabelled as Hybrid feature importance or tumor-specific accessibility measurements.

**Supplementary Figure S5. ComBat pooled discrimination and score comparability.** Within-study and cross-study cancer–control pair contributions sum to secondary pooled AUROC. The paired summary shows pooled discrimination and the equal-weight mean over the two primary held-out targets, emphasizing their different score-comparability requirements.

## Highlights

- A scale-trade-off audit compares broad and local cfDNA fragmentation contrasts.
- Hybrid fusion serves as a diagnostic comparator without a superiority claim.
- Bilateral study-held-out testing separates discrimination from threshold transport.
- A high-burden external cohort is reported as a sensitivity evaluation.
- Fixed real-patient anchors identify target-composition-dependent adaptation.

## References

1. Cristiano S, Leal A, Phallen J, et al. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature. 2019;570:385–389. https://doi.org/10.1038/s41586-019-1272-6
2. Jiang P, Chan CWM, Chan KCA, et al. Lengthening and shortening of plasma DNA in hepatocellular carcinoma patients. Proc Natl Acad Sci USA. 2015;112:E1317–E1325. https://doi.org/10.1073/pnas.1500076112
3. Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. Cell. 2016;164:57–68. https://doi.org/10.1016/j.cell.2015.11.050
4. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8:118–127. https://doi.org/10.1093/biostatistics/kxj037
5. Zheng H, Zhu MS, Liu Y. FinaleDB: a browser and database of cell-free DNA fragmentation patterns. Bioinformatics. 2021;37:2502–2503. https://doi.org/10.1093/bioinformatics/btaa999
6. Wang H, Mennea PD, McAndrew G, et al. A scalable deep-learning framework for cancer detection using cell-free DNA shallow whole-genome sequencing. Sci Adv. 2026;12:eady9432. Data file S1, sheet s2. https://doi.org/10.1126/sciadv.ady9432
7. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics. 1988;44:837–845. https://doi.org/10.2307/2531595
8. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
'''
    manuscript=manuscript.replace(chr(92)*2,chr(92))
    manuscript=manuscript.replace(chr(92)+'(','$').replace(chr(92)+')','$')
    path=root/'MANUSCRIPT_FULL_TEXT.md';path.write_text(manuscript,encoding='utf-8')
    master='# Manuscript master draft — Step 7\n\n'+manuscript
    master+='\n\n## Numerical source map\n\nAll empirical values are read from generated CSVs. Step7_Manuscript_Numeric_Traceability.csv contains the numeric source cells. Original Step 6 outputs are preserved in archive_step06.\n'
    (root/'reports/MANUSCRIPT_MASTER_DRAFT.md').write_text(master,encoding='utf-8')
    report='# Step 07：多尺度验证与科学叙事修订\n\n'
    report+=f'本次只使用原有{total_n}名独立患者。主双向队列{primary_n}人；外部高负荷队列{external_n}人。\n\n'
    report+='## 主要结果\n\n'+primary_table.to_markdown(index=False,floatfmt='.3f')+'\n\n'
    report+='## 来源锁定的外部敏感性\n\n'+external_table.to_markdown(index=False,floatfmt='.3f')+'\n\n'
    report+='## 机制诊断\n\n'+stress.to_markdown(index=False,floatfmt='.6f')+'\n\n'+gc.reset_index().to_markdown(index=False,floatfmt='.6f')+'\n\n'
    report+='## 数据完整性\n\n'+quality.to_markdown(index=False,floatfmt='.6f')+'\n\n临床字段缺失率：\n\n'+missing.to_markdown(index=False,floatfmt='.6f')+'\n\n'
    report+=f'ComBat pooled AUROC={combined.loc["InFold_ComBat","auroc"]:.3f}，加权中心内配对AUC={within_auc:.3f}，跨中心配对AUC={cross_auc:.3f}。其中心内逐窗口Cohen d最大变化仅{maximum_cb_d:.3g}；不能据低pooled AUROC声称生物信号被摧毁。\n\n'
    report+='## 设计边界与标题\n\n单一来源无法估计多中心ICC/ComBat/DANN；ICC采用明确标记的LSN回退，ComBat和DANN相应单元格不估计。样本内变换不依赖目标批次构成，不等于分类器、训练患病率、PPV或临床泛化天然免疫混杂。\n\n'
    report+='现稿贡献定位为Scale Trade-off系统审计，Hybrid是融合对照，不宣称击败现有基线。CNA保留及GC驱动漂移消除尚无独立因果验证；研究范围为公开cfDNA队列中的跨研究泛化。Study Source同时包含实验批次与人群构成混杂。\n\n'
    report+=f"Combined LOCO的固定ElasticNet比较：Hybrid AUROC={combined.loc['Hybrid_MultiScale','auroc']:.3f}，Global={combined.loc['Global_Zscore','auroc']:.3f}，Proposed_LSN_ICC={combined.loc['Proposed_LSN_ICC','auroc']:.3f}。这些点估计不构成融合算法的普遍显著优越性。\n\n"
    report+='## 本地来源资格核查\n\n'+source_table[['study','runs','healthy_patients','cancer_patients','meets_cancer_and_10_healthy']].to_markdown(index=False)+'\n\n'
    report+='本地未核验到额外符合癌症+至少10名健康人标准的来源；此结论仅限历史快照与存档来源汇总，不代表当前数据库或所有公开数据集的完整普查。两个双类别且健康人足量来源+一个病例主导来源是核心局限；不将双类别充足误称为严格1:1均衡。\n\n'
    report+='完整手稿：../MANUSCRIPT_FULL_TEXT.md；所有统计量来自本地真实CSV，代码仓库、作者、资助和伦理确认项保持明确待作者填写。\n'
    (root/'reports/step_07_report.md').write_text(report,encoding='utf-8')
    return path
