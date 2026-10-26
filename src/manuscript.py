"""Create a manuscript-ready data draft exclusively from verified study outputs."""
import importlib.metadata
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import beta
from .common import sha256,write_json
from .study import load_clinical,export_table
from .models import make_model


def compile_clinical_assets(cfg):
    root=Path(cfg['paths']['results']);folder=root/'tables'
    meta,R,lengths,bins=load_clinical(cfg,True)
    export_table(meta,cfg,'Analysis_Participants')
    peaks=[]
    for center,group in meta[meta.label.eq(0)].groupby('center'):
        for bp in [143,167]:
            peaks.append(dict(center=center,length_bp=bp,healthy_n=len(group),mean_frequency=float(lengths[group.index,bp-80].mean())))
    export_table(pd.DataFrame(peaks),cfg,'Healthy_Reference_Lengths')
    metrics=pd.read_csv(folder/'Table3_Benchmark_Metrics_Full.csv')
    main=metrics[metrics.center.eq('AllCenters')]
    if len(main)!=48 or not main.status.eq('complete').all():raise ValueError('All 24 combinations x 2 schemes must be complete')
    predictions=pd.read_parquet(folder/'Clinical_OOF_Predictions.parquet')
    exact=[]
    for row in metrics.itertuples():
        g=predictions[(predictions['mode']==row.mode)&(predictions.pipeline==row.pipeline)&(predictions.model==row.model)]
        if row.center!='AllCenters':g=g[g.center.eq(row.center)]
        healthy=g[g.label.eq(0)];n=len(healthy);success=int((healthy.probability<=healthy.threshold).sum())
        lo=float(beta.ppf(.025,success,n-success+1)) if success else 0.
        hi=float(beta.ppf(.975,success+1,n-success)) if success<n else 1.
        exact.append((lo,hi))
    metrics['specificity_exact_low']=[a for a,b in exact];metrics['specificity_exact_high']=[b for a,b in exact]
    export_table(metrics,cfg,'Table3_Benchmark_Metrics_Full')
    parameters={'seed':cfg['seed'],'participants':len(meta),'cancer_participants':int(meta.label.sum()),
        'healthy_participants':int(meta.label.eq(0).sum()),'windows':len(bins),'studies':meta.center.nunique(),
        'radius':cfg['study']['radius'],'alpha':cfg['study']['exclusion_fraction'],'relative_epsilon':cfg['study']['relative_epsilon'],
        'permutations':cfg['study']['permutations'],'bootstrap':cfg['study']['bootstrap_replicates'],
        'random_outer_folds':5,'inner_calibration_folds':cfg['study']['inner_calibration_folds'],
        'target_specificity':.95,'healthy_FDR_threshold':.01,'vulnerability_ratio_threshold':5,
        'excluded_unresolved_patient_identity':len(pd.read_csv(folder/'Analysis_Identity_Exclusions.csv'))}
    export_table(pd.DataFrame([{'parameter':k,'value':v} for k,v in parameters.items()]),cfg,'Methods_Execution_Parameters')
    software=[]
    for package in ['numpy','pandas','pyarrow','scipy','scikit-learn','statsmodels','torch','xgboost','lightgbm','umap-learn','matplotlib','seaborn','requests']:
        software.append(dict(package=package,version=importlib.metadata.version(package)))
    software.append(dict(package='Python',version=sys.version.split()[0]));export_table(pd.DataFrame(software),cfg,'Methods_Software_Versions')
    model_parameters=[]
    for model in ['elasticnet','rf','xgboost','lightgbm','mlp','dann']:
        for parameter,value in make_model(model,cfg).get_params().items():
            model_parameters.append(dict(model=model,parameter=parameter,value=str(value)))
    export_table(pd.DataFrame(model_parameters),cfg,'Methods_Model_Parameters')
    perma=pd.read_csv(folder/'Audit_PERMANOVA.csv');part=pd.read_csv(folder/'Bin_Variance_Partition.csv')
    drift=pd.read_csv(folder/'Healthy_Drift.csv');ablation=pd.read_csv(folder/'Table4_Ablation_Study.csv')
    correlations=pd.read_csv(folder/'Genomic_Spearman_Associations.csv');effects=pd.read_csv(folder/'Disease_Effect_Sizes.csv')
    enrich=pd.read_csv(folder/'Filtered_Window_Enrichment.csv');concord=pd.read_csv(folder/'Feature_Rank_Concordance.csv')
    original=pd.read_csv(folder/'Table1_Cohort_Characteristics.csv');significance=pd.read_csv(folder/'Clinical_Paired_Significance.csv')
    audit=perma.set_index(['family','term'])
    raw=audit.loc[('raw_ratio','center')];disease=audit.loc[('raw_ratio','disease')]
    centered=audit.loc[('median_centered_ratio','center')]
    core=main[(main['mode']=='loco')&(main.model=='elasticnet')].set_index('pipeline')
    rawmean=main[main.pipeline.eq('Raw_Ratio')].groupby('mode').macro_center_auroc.mean()
    grid=ablation[(ablation.center=='MacroCenters')&(ablation.pipeline=='Proposed_LSN_ICC')]
    effect=effects.pivot(index=['bin_id','center'],columns='pipeline',values='cohen_d').dropna()
    gc=correlations[correlations.annotation.eq('gc_fraction')].iloc[0]
    gext=correlations[correlations.annotation.eq('gc_extremeness')].iloc[0]
    tel=correlations[correlations.annotation.eq('log_telomere_distance')].iloc[0]
    centromere=correlations[correlations.annotation.eq('log_centromere_distance')].iloc[0]
    proposed=core.loc['Proposed_LSN_ICC'];base=core.loc['Raw_Ratio']
    clinical_counts=meta.groupby('center').agg(n=('patient_id','size'),cancer=('label','sum')).reset_index()
    clinical_counts['healthy']=clinical_counts.n-clinical_counts.cancer
    export_table(clinical_counts,cfg,'Analysis_Cohort_Counts')
    # Every inserted statistical value below is read from generated CSV tables.
    n=int(pd.read_csv(folder/'Methods_Execution_Parameters.csv').set_index('parameter').loc['participants','value'])
    text=f'''# Manuscript master draft

Data-derived working title: Study-specific batch vulnerability and threshold transfer in plasma cfDNA fragmentomics

## Abstract

**Background:** Fragmentomic classifiers require evaluation across source studies with different baseline fragmentation profiles.
**Methods:** We audited {int(original.selected_samples.sum())} publicly available plasma specimens and analysed one prespecified specimen for each of {n} identified participants across {len(original)} studies. {len(bins)} retained autosomal windows were evaluated using conditional Euclidean PERMANOVA, coverage-adjusted type III variance fractions, and healthy-only rank tests. Four preprocessing pipelines and six classifiers underwent patient-isolated random and leave-one-study-out validation. Source thresholds were locked using eligible inner out-of-fold healthy predictions; uncertainty was estimated by 1000 stratified patient bootstraps.
**Results:** Study provenance explained {raw.R2*100:.3f}% of raw-ratio variation, compared with {disease.R2*100:.3f}% attributed uniquely to disease status (permutation p={raw.p:.3f} and p={disease.p:.3f}). Sample-wise median centering reduced the unique study fraction to {centered.R2*100:.3f}%. Healthy baseline differences met FDR<0.01 in {drift.q.lt(.01).sum()} of {len(drift)} windows. Across the six classifiers, mean macro-study Raw Ratio AUROC was {rawmean['random']:.3f} in random validation and {rawmean['loco']:.3f} in held-out-study validation. For the prespecified ElasticNet comparison, pooled LOCO AUROC was {base.auroc:.3f} (95% CI {base.auroc_low:.3f}–{base.auroc_high:.3f}) for Raw Ratio and {proposed.auroc:.3f} ({proposed.auroc_low:.3f}–{proposed.auroc_high:.3f}) for LSN–ICC. Achieved external specificity was {base.specificity95:.3f} and {proposed.specificity95:.3f}, respectively.
**Conclusions:** Source-dependent baseline variation and threshold transfer are measurable components of cfDNA classifier evaluation. The observed pipeline comparisons identify where local contrasts change predictive discrimination and calibration across studies.

## Methods

### Clinical source and independent-participant cohort

Clinical labels, original patient/specimen identities, study provenance and depth were recovered from the published UNITE data file S1/s2. Only plasma WGS specimens with original depth in the locked 0.1–5× interval, explicit cancer or Healthy labels, and an eligible study were selected. Cirrhosis, hepatitis, transplantation, inflammatory and autoimmune labels were not pooled into Healthy. Interval objects were streamed from explicit FinaleDB hg19 URLs. No raw BAM, FASTQ or compressed fragment files were stored. The completed extraction cohort contained {int(original.selected_samples.sum())} specimens. {parameters['excluded_unresolved_patient_identity']} specimen with an unresolved source patient ID was excluded from patient-independent analysis without inferring an alias. Earliest timepoint, then shallowest depth, then run ID determined the representative specimen. The analysis cohort was:

{clinical_counts.to_markdown(index=False)}

These are study proxies, not independently verified physical hospital centers. The Snyder healthy stratum contains one participant; its uncertainty and unavailable two-class covariance estimates are reported explicitly.

### Fragmentomic features and batch audit

Autosomal MAPQ≥30, 80–400 bp intervals were assigned by midpoint to the fixed filtered hg19 grid. For each window, S counted 100–150 bp and L counted 151–220 bp fragments. The audit and classifiers used R=S/(L+1), recomputed from counts; the extraction-stage 0.5-pseudocount ratio was not silently substituted. The centered audit ratio was R minus the sample-wide median. Global fragment frequencies retained each integer length from 80 to 400 bp. Euclidean analyses used descriptive column standardization, independently of classifier preprocessing. PERMANOVA tested the unique study or disease term conditional on the other term and original coverage depth with 999 Freedman–Lane permutations, restricted within disease for the study test and within study for the disease test. The per-window model was Feature ~ C(study) + disease + original_coverage_depth. Unique/type III eta² was the reduced-to-full residual SS increment divided by total SS; components need not sum to one. Vulnerable windows met eta²_batch/(eta²_bio+1e-8)>5. Healthy-only Kruskal–Wallis p values were BH corrected across windows; a sensitivity analysis removed studies with fewer than five healthy participants.

### Preprocessing and classifier models

Raw Ratio used train-only imputation and StandardScaler. Global Zscore first standardized each sample across its own windows and then applied train-only StandardScaler. InFold ComBat estimated source location/scale parameters from training data; its unseen-center robust median/IQR mapping used the unlabeled target batch and was labelled **transductive**, separately from the inductive pipelines. LSN used other retained windows on the same chromosome arm, up to k neighbors on each side:

$$LSN_{{si}}=\\log_2\\frac{{R_{{si}}+\\epsilon_s}}{{\\mathrm{{median}}_{{j\\in N_k(i),j\\ne i}}R_{{sj}}+\\epsilon_s}},\\qquad \\epsilon_s=10^{{-6}}\\mathrm{{median}}_jR_{{sj}}.$$

The relative epsilon preserves global multiplicative scale invariance to numerical precision. Additive and arbitrary nonlinear invariance do not follow from this formula; counterexample tests are included. The locked radius was k={cfg['study']['radius']}. The source-healthy filter estimated one-way center clustering ICC(1), not repeated-subject reliability, using n0=(N−sum(n_g²)/N)/(G−1), sigma²_center=max((MS_between−MS_within)/n0,0), and ICC_center=sigma²_center/(sigma²_center+MS_within). The highest {cfg['study']['exclusion_fraction']*100:.3f}% of windows were removed; training signal/noise broke equal-ICC ties. RobustScaler was fitted on the retained training features. All stateful steps were refitted inside each outer and inner training split.

ElasticNet, Random Forest, CUDA XGBoost, LightGBM, two-layer CUDA MLP and source-only CUDA DANN were evaluated with fixed settings. Neural models used dropout and weight decay and never consumed the held-out study during adversarial training. Exact software and model settings are supplied in Methods_Software_Versions.csv and Methods_Model_Parameters.csv:

{pd.DataFrame(software).to_markdown(index=False)}

### Validation, operating thresholds and inference

Random validation used five patient-isolated stratified folds; LOCO held out each full study. Three-fold inner calibration used study/diagnosis strata with enough participants; rarer strata remained training-only and supplied no in-sample substitute calibration scores. This preserves two-source healthy-filter estimation when a source has one healthy participant. A higher empirical 0.95 quantile of eligible source-healthy OOF scores locked the threshold, with p>threshold defining a positive result. Test labels never selected thresholds, filters or hyperparameters. AUROC, AUPRC, transferred sensitivity and specificity, Brier score and 10-bin ECE were calculated. Bootstrap intervals use 1000 patient resamples stratified by fixed study and diagnosis, conditional on fitted OOF scores; they do not include model-retraining or population-center uncertainty. Exact binomial specificity intervals supplement bootstrap intervals for the singleton-control study. DeLong compared fixed model pairs within the same held-out study when both classes contained at least two people. Wilcoxon comparisons across three overlapping-training LOCO folds are exploratory; BH correction was applied to finite paired-test p values. The fixed k/alpha configuration was not chosen from the outer sensitivity surface.

### Reference annotations and interpretation

GC fractions were calculated from streamed UCSC hg19 sequence, conditional on called A/C/G/T bases. RepeatMasker and segmental-duplication fractions used interval unions. Distances came from the fixed cytoband-derived arm and centromere boundaries. A/B reference coverage used GM12878 subcompartments; it is not a tumor-specific chromatin assay. Genomic correlations used Spearman tests and BH correction, and regression used chromosome-cluster-robust standard errors. The removed-versus-retained comparison used majority agreement across the three LOCO training masks. Fisher enrichment p values used the original contingency counts; zero-cell odds-ratio intervals used a 0.5 correction. Cancer-versus-healthy Cohen d was calculated within studies with at least two members of each class. Sequence and reference-chromatin associations are exploratory mechanisms, not direct proof of tumor nucleosome depletion.

## Results

### Source-study variation and healthy baseline displacement

The independently identified cohort comprised {n} participants. Conditional study and disease fractions in raw ratios were {raw.R2*100:.3f}% and {disease.R2*100:.3f}%, respectively. In median-centered ratios, the study fraction was {centered.R2*100:.3f}%. The global-length study fraction was {audit.loc[('length_frequency','center'),'R2']*100:.3f}%. Healthy-only baseline testing identified {drift.q.lt(.01).sum()}/{len(drift)} significant windows at FDR<0.01; the minimum-five-controls sensitivity analysis identified {drift.q_min5.lt(.01).sum()}/{len(drift)}. {int(part.query('family=="raw_ratio"').vulnerable.sum())} raw-ratio windows exceeded the locked batch-to-disease variance ratio criterion. These source associations include collection, processing, population and study-composition differences.

### Random mixing and transfer to held-out studies

Across the six Raw Ratio classifiers, mean macro-study AUROC was {rawmean['random']:.3f} under random validation and {rawmean['loco']:.3f} under LOCO. Individual classifier/study results, uncertainty, achieved sensitivity and external specificity are reported in Table3_Benchmark_Metrics_Full.csv. The following ElasticNet summary uses the same prespecified model for all pipeline comparisons:

{core[['auroc','auroc_low','auroc_high','auprc','sensitivity95','specificity95','brier','ece']].to_markdown(floatfmt='.3f')}

### Local normalization, filtering and operating-point behavior

For ElasticNet, LSN–ICC yielded pooled LOCO AUROC {proposed.auroc:.3f}, AUPRC {proposed.auprc:.3f}, sensitivity {proposed.sensitivity95:.3f}, achieved specificity {proposed.specificity95:.3f}, Brier score {proposed.brier:.3f}, and ECE {proposed.ece:.3f}. Raw Ratio yielded corresponding AUROC {base.auroc:.3f}, AUPRC {base.auprc:.3f}, sensitivity {base.sensitivity95:.3f}, specificity {base.specificity95:.3f}, Brier {base.brier:.3f}, and ECE {base.ece:.3f}. Global Zscore had pooled LOCO AUROC {core.loc['Global_Zscore','auroc']:.3f}; the local pipeline did not exceed this comparator for the prespecified model. Over the prespecified radius/filter grid, macro-study AUROC ranged from {grid.auroc.min():.3f} to {grid.auroc.max():.3f}, indicating material parameter dependence rather than demonstrated insensitivity. Table4 provides all fold-specific values and the isolated LSN/filter ablations. The sensitivity grid was not used to reselect the primary configuration.

### Genomic associations and distribution of disease contrasts

GC fraction correlated with raw-ratio center eta² at rho={gc.rho:.3f}, p={gc.p:.3g}; GC extremeness gave rho={gext.rho:.3f}, p={gext.p:.3g}. Log telomere distance gave rho={tel.rho:.3f}, p={tel.p:.3g}. Log centromere distance gave rho={centromere.rho:.3f}, p={centromere.p:.3f}; this association was not significant. Enrichment results for the actual training-fold filter masks were:

{enrich[['region','odds_ratio','ci_low','ci_high','p','q']].to_markdown(index=False,floatfmt='.3f')}

Of the {len(enrich)} reference-defined enrichment tests, {enrich.q.lt(.05).sum()} met BH FDR<0.05. Nominal odds-ratio patterns therefore did not establish corrected region enrichment in this cohort.

Across estimable within-study window contrasts, median absolute Cohen d was {effect.Raw_Ratio.abs().median():.3f} for raw ratios and {effect.LSN.abs().median():.3f} for LSN; {100*(effect.LSN.abs()>effect.Raw_Ratio.abs()).mean():.3f}% of contrasts had a larger absolute local-normalization effect. Effects therefore varied by window and study. Shared-feature rank concordance across the three LOCO fits ranged from {concord.rho.min():.3f} to {concord.rho.max():.3f}. Reference associations and ranking patterns quantify genomic context without asserting universal preservation or amplification of tumor-specific signals.

## Figure legends

**Figure 1. Study-associated batch vulnerability in plasma cfDNA.** (a) Independent-participant composition and observed healthy mean 80–400 bp frequency curves; vertical markers identify 143 and 167 bp reference lengths. (b) Raw-ratio UMAP coordinates colored by study or diagnosis, with coverage-adjusted unique Type III R² and 999-permutation p values. (c) Unique study and disease eta² across retained autosomal windows in genomic order; these are type III components, not an additive causal decomposition. (d) Healthy sample-wide median ratio distributions and center-clustering ICC(1) density. The single Snyder healthy participant is shown as observed; the minimum-five-controls sensitivity table is supplied separately. Descriptive embedding transformations were not passed into classifier fitting.

**Figure 2. Random mixing, study predictability and transferred thresholds.** (a) Raw Ratio macro-study AUROC for each of six classifiers under random and held-out-study validation. (b) Random Forest study-prediction balanced accuracy in all participants and in healthy cohorts with adequate size; study predictability alone does not establish a causal technical effect. (c) Achieved ElasticNet specificity when source-locked 95%-specificity thresholds are transferred unchanged, with exact binomial 95% intervals. The singleton-control study has correspondingly broad uncertainty. (d) Raw Ratio held-out-study AUROC heatmap on the full zero-to-one scale. Numerical scores are observed evaluations, without requiring a drop for every classifier.

**Figure 3. Local contrast preprocessing and observed cross-study performance.** (a) Same-arm local log-ratio, source-healthy center-clustering filter and train-only robust scaling, with observed healthy baseline traces. (b) Pooled LOCO ROC and precision–recall curves for the prespecified ElasticNet model; DeLong comparisons are held-out-study-specific, and the one-control study is marked non-estimable. Pooling scores from different trained models is descriptive. InFold ComBat uses unlabeled target-batch adaptation and is explicitly transductive. (c) Probability calibration and Brier scores. (d) The complete prespecified radius/filter surface for ElasticNet macro-study AUROC; it is a sensitivity analysis, not outer-test-driven model selection.

**Figure 4. Reference genomic context and fold-dependent feature retention.** (a) Window-level GC and telomere-distance associations with source eta²; smoothed curves are descriptive LOWESS fits. (b) Odds ratios and intervals comparing majority-fold removed windows with retained windows for reference-defined regions; Fisher tests and BH values are supplied. (c) Observed chr8 healthy and cancer median trajectories before and after LSN, with solid and dashed lines. (d) Top ElasticNet window coefficient ranks across LOCO training folds; gray cells indicate filtering, not observed zero rank. A/B annotations refer to GM12878 reference chromatin and do not directly measure tumor-specific accessibility.

## Highlights — Computers in Biology and Medicine

- Study provenance explains {raw.R2*100:.3f}% of raw cfDNA ratio variation.
- Patient-isolated validation quantifies cross-study operating-point transfer.
- Local contrast and healthy-baseline filtering are tested without target labels.
- All comparisons retain source provenance and explicit uncertainty scopes.

## Highlights — Computer Methods and Programs in Biomedicine

- Streaming extraction avoids storing raw sequencing files.
- Four preprocessing pipelines share locked patient-isolated splits.
- Threshold calibration and healthy-baseline filtering remain inside training.
- Checkpointed model evaluation supports reproducible study-held-out audits.

## Code availability — template requiring author completion

The analysis code, configuration, dependency specification and test suite are available at **[PUBLIC REPOSITORY URL AND RELEASE TAG TO BE CONFIRMED BY AUTHORS]**. Feature tables are derived from source-linked public interval objects; the metadata and annotation source manifests identify the corresponding publications and files. No public repository URL or release status is asserted by this template.

## Publication-route notes

The intended route is subscription publication in either target journal. Elsevier's general subscription policy states that no publication fee applies to that route; target-journal options and any optional charges must be checked against the current author instructions before submission. The specific target options pages returned access restrictions during this run, so no target-specific fee guarantee is asserted. Official general policy: https://www-prod.elsevier.com/publishing/publish-in-a-journal/journal-selection/special-issues ; target pages: https://www.sciencedirect.com/journal/computers-in-biology-and-medicine/publish/open-access-options and https://www.sciencedirect.com/journal/computer-methods-and-programs-in-biomedicine/publish/open-access-options .

## Statistical traceability

All quantitative Results and Abstract values are programmatically read from generated CSV outputs. Methods constants and software versions are exported to dedicated CSV tables. The paired significance table is reproduced in full below; unavailable entries remain missing.

{significance.to_markdown(index=False,floatfmt='.3f')}
'''
    path=root/'reports/MANUSCRIPT_MASTER_DRAFT.md';path.write_text(text,encoding='utf-8')
    inventory=refresh_clinical_manifest(cfg)
    return {'manuscript':str(path),'independent_patients':n,'complete_method_scheme_rows':len(main),'assets':len(inventory)}


def refresh_clinical_manifest(cfg):
    root=Path(cfg['paths']['results'])
    paths=sorted(p for sub in ['tables','figures','reports'] for p in (root/sub).glob('*') if p.is_file() and p.name!='asset_manifest.json')
    if (root/'MANUSCRIPT_FULL_TEXT.md').exists():paths.append(root/'MANUSCRIPT_FULL_TEXT.md')
    inventory={str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in paths}
    project=root.parent
    sources=sorted(p for sub in ['src','scripts','tests','config'] for p in (project/sub).glob('*') if p.is_file())
    sources += [project/name for name in ['AGENTS.md','README.md','requirements.txt','environment.yml','requirements-tested.txt']]
    code={str(p.relative_to(project)):sha256(p) for p in sources if p.exists()}
    inputs=[Path(cfg['paths']['metadata']),Path(cfg['paths']['processed'])/cfg['extraction']['final_file'],
            Path(cfg['paths']['reference'])/cfg['reference']['window_table'],Path(cfg['paths']['reference'])/'hg19_window_annotations.parquet']
    write_json({'assets':inventory,'source_files':code,'inputs':{str(p):sha256(p) for p in inputs},'raw_intervals_stored':False,
                'manuscript':str(root/'reports/MANUSCRIPT_MASTER_DRAFT.md')},root/'asset_manifest.json')
    return inventory
