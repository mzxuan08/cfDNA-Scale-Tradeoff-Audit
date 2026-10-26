# Scale Trade-offs in Within-Sample Self-Normalization for Cross-Study Generalization in Public cfDNA Cohorts

**Author:** Zhenxuan Ma

**Affiliation:** Hubei University of Medicine, Shiyan 442000, Hubei, China

**Corresponding author:** Zhenxuan Ma; mzxmc08@gmail.com

**ORCID:** Not registered.

## Abstract

**Background:** Normalization determines which cfDNA fragmentation contrasts a model retains and how its scores transfer across studies. We audited the relationship between normalization scale, study-source variation and target-cohort composition in public cfDNA data.

**Methods:** We analysed 585 source-identified participants using a primary bilateral benchmark of 548 and a separately held-out, case-dominated cohort of 37. Global, local and hybrid representations were evaluated with source-only selection and locked thresholds. Cancer–control pair concordance and fixed-patient anchor experiments separated within-study ranking, cross-study score comparability and target-composition dependence.

**Results:** Combined LOCO ElasticNet AUROCs were 0.759 for Global Z-score, 0.731 for Hybrid and 0.724 for fixed LSN-ICC. Median absolute arm-level fragmentation effects were 0.270 for Global and 0.244 for local normalization; residual GC correlations with study variance were 0.158 and 0.405. For ComBat with unseen-target alignment, within-study and cross-study concordances were 0.669 and 0.366. Changing the companions of fixed real patients shifted their probabilities by up to 0.214 under target-batch adaptation and by 0 under within-sample transformations.

**Conclusions:** Normalization scale, study-wise discrimination and score comparability are distinct targets for cross-study cfDNA evaluation. The audit makes their trade-offs measurable and identifies patient-wise self-normalization as a direct way to remove companion-dependent target preprocessing.

**Keywords:** cell-free DNA; fragmentomics; cross-study validation; sample-wise normalization; batch composition; threshold transport

## 1. Introduction

Cross-study cfDNA fragmentomics requires representations that retain disease-associated fragmentation contrasts while controlling study-source variation and maintaining interpretable scores in a new cohort. Plasma DNA fragmentation reflects nucleosome organization and tissue contribution, and genome-wide profiles support cancer classification alongside sequence-variant assays [1–3]. Study Source combines specimen processing, sequencing procedures, recruitment and population composition. Evaluating transport therefore requires both patient-level discrimination and the achieved operating point of a threshold locked in the source study.

Normalization scale offers a way to examine these requirements directly. Global normalization retains relative deviations across the genome while removing a sample-wide offset and scale. A local contrast compares adjacent windows and attenuates slowly varying backgrounds, including broad disease-associated differences. In the window-level audit, median absolute within-study disease contrasts were 0.384 for raw ratios and 0.263 for a local representation. This motivates a comparison of global, local and fused representations, with fragmentation contrast retention as the measured endpoint.

Cross-sample alignment adds a separate dependence on cohort composition. An empirical-Bayes harmonizer estimates study-specific location and scale parameters [4]. Adapting it to an unseen study with target medians or quantiles makes those estimates properties of the target patient mixture. Under $F_c=(1-\pi_c)F_{0c}+\pi_cF_{1c}$, changing the case fraction $\pi_c$ can change batch-level quantiles while the evaluated individuals remain fixed. Within-sample transformations offer a precise alternative: each patient supplies the reference for their own transformation. This functional distinction can be tested by holding real patients fixed and changing their companions.

We present a Scale Trade-off audit that connects representation, discrimination and score transport in public cfDNA cohorts. Global and same-arm local normalization quantify how scale changes broad fragmentation contrasts and residual GC-associated study variation; Hybrid tests their complementarity in a source-trained classifier. Bidirectional study-held-out testing and a separate external sensitivity cohort establish the evaluation scope. Pair concordance decomposition and fixed-patient composition stress tests then identify whether pooled score differences arise from within-study ranking, comparability across fitted models or target-batch adaptation.

## 2. Materials and methods

### 2.1. Data provenance and independent-participant cohort

The analysis used the existing extracted feature matrix for 585 source-identified participants. Clinical metadata came from UNITE supplementary data file S1, sheet s2, and fragment intervals were streamed from explicit FinaleDB objects [5,6]. Inclusion required plasma WGS, original sequencing depth within 0.1–5.0× and explicit cancer or Healthy source labels. The extraction cohort contained 628 specimens. One prespecified specimen per known patient was retained by earliest time point, then shallowest sequencing depth, then sequencing-run identifier; a specimen with unresolved patient identity was isolated. Healthy controls required an explicit Healthy source category. Hepatitis, cirrhosis, inflammatory disease, transplantation and autoimmune disease retained their original categories.

The primary cohort contained the Cristiano and Jiang studies. Each bilateral direction trained on one entire study and tested the other. The Snyder study was held out from all primary training, primary parameter selection and primary threshold calibration. Its high case burden makes cancer sensitivity the useful external endpoint; it does not supply an adequately sampled healthy specificity estimate. A separate secondary analysis used three-study LOCO across all participants, with the other two studies used for each held-out target. Snyder can therefore enter source training only in that explicitly secondary protocol, never in the primary bilateral experiment.

| center                |   n |   cancer |   healthy |   cancer_percent |
|:----------------------|----:|---------:|----------:|-----------------:|
| Cristiano et al, 2019 | 427 |      218 |       209 |           51.054 |
| Jiang et al, 2015     | 121 |       89 |        32 |           73.554 |
| Snyder et al, 2016    |  37 |       36 |         1 |           97.297 |

Study Source is an observational proxy for experimental and population-composition differences. The source publications have different recruitment aims, so their case proportions describe the selected cohorts. Participant identity follows published source identifiers; undocumented aliases across studies remain unresolved.

The local source census covered 2362 FinaleDB run records from 5 studies. Eligibility for an additional two-class source required a documented plasma WGS cohort with at least one cancer participant and at least ten unique, explicitly labelled Healthy participants; healthy controls were not borrowed across studies. The broad census below precedes depth and QC restrictions, and counts only nonmissing patient identifiers:

| study                     |   runs |   healthy_patients |   cancer_patients | meets_cancer_and_10_healthy   |
|:--------------------------|-------:|-------------------:|------------------:|:------------------------------|
| Adalsteinsson et al, 2017 |   1508 |                  0 |               528 | False                         |
| Cristiano et al, 2019     |    538 |                262 |               231 | True                          |
| Jiang et al, 2015         |    225 |                 32 |                90 | True                          |
| Snyder et al, 2016        |     60 |                  4 |                48 | False                         |
| Sun et al, 2019           |     31 |                  0 |                12 | False                         |

The archived database publication listed 7 sources and 2579 sequencing records. Its additional urinary and prenatal single-gene studies did not document an eligible plasma cancer/control cohort in the local materials. Unavailable run-level case/control counts remain unknown. The census describes the archived local inventory at the stated eligibility threshold.

### 2.2. Genomic windows, fragmentation features and descriptive audit

The hg19/GRCh37 reference retained 477 autosomal non-overlapping complete windows of 5 Mb after fixed centromere, terminal-risk, cytoband-risk and blacklist criteria. Only autosomal 80–400 bp fragments with MAPQ at least 30 were aggregated; overlaps with blacklisted intervals were excluded. Midpoints assigned fragments to windows. In each window, $S_i$ counted 100–150 bp fragments and $L_i$ counted 151–220 bp fragments. The classifier ratio was

$$R_{si}=\frac{S_{si}}{L_{si}+1}.$$

Window counts, IQR and the global integer-length frequency distribution were retained in compact local tables. The analyses used float32 feature matrices. GC fractions came from reference sequence, while RepeatMasker and segmental-duplication fractions used interval unions. GM12878 chromatin compartments are a reference annotation, not a tumor accessibility assay.

Descriptive PERMANOVA and window variance models adjusted for disease status, study and original sequencing depth as appropriate. PERMANOVA used Euclidean geometry after descriptive column standardization and 999 restricted Freedman–Lane residual permutations. Window-level unique/type III eta-squared was the reduced-to-full residual sum-of-squares increment divided by total sum of squares; overlapping explanatory components can yield a sum different from one. Descriptive embeddings were kept separate from predictive preprocessing. GC associations were calculated independently for raw, global and local representations.

### 2.3. Global, local and hybrid representations

Raw Ratio used training-only imputation and column standardization. Global Z-score first used each individual's own window mean and standard deviation:

$$G_{si}=\frac{R_{si}-\bar R_s}{\operatorname{sd}_i(R_s)}.$$

Its subsequent column scaler was fitted on source training patients. The local transform used other retained windows on the same chromosome arm, up to $k$ on either side, excluding the focal window:

$$L_{si}^{(k)}=\log_2\left[\frac{R_{si}+\epsilon_s}{\operatorname{median}_{j\in N_k(i),j\ne i}R_{sj}+\epsilon_s}\right],
\qquad \epsilon_s=10^{-6}\operatorname{median}_iR_{si}.$$

The relative epsilon preserves global multiplicative invariance to numerical precision. The local operation measures a relative contrast; additive shifts and arbitrary nonlinear distortions are outside that invariance. Hybrid MultiScale concatenated the global and local blocks,

$$H_s=[G_s,L_s^{(3)}],$$

and fitted its column scaler within the source training split. Both blocks and the fixed Hybrid radius were specified before outer evaluation. Elastic-net logistic regression supplied the primary classifier, learning regularized coefficients for the two correlated blocks. Additional fixed classifiers provided secondary algorithm comparisons.

The fixed LSN-ICC pipeline ranked healthy-center clustering scores within source training data, removed the highest quarter and applied RobustScaler. The scores are one-way between-center variance ratios among different healthy people. In a bilateral direction, the single source study makes cross-center filtering unidentifiable; the recorded operational representation is therefore LSN only. This fallback is labelled in the results and uses source data exclusively.

### 2.4. Nested source-only selection and locked thresholds

Tuned LSN considered radii $k\in\{1,2,3,5,8\}$ and eligible exclusion fractions $\alpha\in\{0,0.15,0.25,0.35,0.50\}$. With one healthy source study, zero exclusion was the eligible setting. Three inner folds selected the highest mean validation AUROC across source strata containing at least two participants of each class, with deterministic tie-breaking and the prespecified ElasticNet model. This analysis evaluated source-only scale selection separately from fixed Hybrid.

Source threshold calibration used three-fold out-of-fold predictions from source patients. For Tuned LSN, parameter selection was repeated inside each calibration-training split before predicting its validation patients. Thus calibration outcomes never participated in the parameter search that generated their own scores. Rare study/diagnosis strata remained training-only and supplied no in-sample substitute score. The higher empirical 0.95 quantile of eligible source-healthy out-of-fold probabilities locked the decision threshold. Probabilities strictly above that threshold were positive. A model was then fitted to all source patients and applied unchanged to the target patients. No target diagnosis selected a threshold, model setting, feature mask or scale. Random and bilateral validation shared the primary patient universe, but their training-set sizes differed; their contrast does not isolate a pure technical batch effect.

### 2.5. Comparator identifiability and classifier implementation

The five core pipelines were Raw Ratio, Global Z-score, InFold ComBat, fixed Proposed LSN-ICC and Hybrid MultiScale. Tuned LSN used the primary ElasticNet classifier. Random Forest, CUDA XGBoost, LightGBM, CUDA MLP and source-only CUDA DANN supplied locked secondary models. DANN used training-source domains exclusively. Multi-source DANN and standard multi-batch source ComBat were estimable in mixed-study random CV and two-source/one-target Combined LOCO. Their one-source bilateral cells were marked unavailable.

InFold ComBat combined source-only empirical-Bayes harmonization with an unlabeled target median/IQR mapping for an unseen study. This explicitly transductive comparator evaluates a particular adaptation rule. The other fitted transformations process each target patient independently of companion patients; their source-trained scalers and classifiers still reflect the source training population.

All random operations used seed 2026. Neural hidden dimensions and batch sizes respected the original memory constraints, and CUDA availability and peak process allocation were checked. Exact fixed model parameters are supplied in Methods_Model_Parameters.csv. The executable software versions were:

| package      | version     |
|:-------------|:------------|
| numpy        | 2.2.6       |
| pandas       | 2.3.3       |
| pyarrow      | 24.0.0      |
| scipy        | 1.15.3      |
| scikit-learn | 1.7.2       |
| statsmodels  | 0.15.0      |
| torch        | 2.5.1+cu121 |
| xgboost      | 2.1.4       |
| lightgbm     | 4.6.0       |
| umap-learn   | 0.5.12      |
| matplotlib   | 3.10.9      |
| seaborn      | 0.13.2      |
| requests     | 2.34.2      |
| Python       | 3.10.10     |

### 2.6. Metrics and statistical comparison

We report direction-specific AUROC and AUPRC, pooled held-out predictions, equal-weight mean of the two primary study AUROCs, source-threshold sensitivity and specificity, Brier score and 10-bin expected calibration error. Specificity deviation is the achieved target specificity minus 0.95 in percentage points. We calculated 1000 patient bootstraps stratified by the fixed study and diagnosis, conditional on the fitted out-of-fold scoring rules. These intervals do not include refitting, method-development or new-center uncertainty. Exact binomial intervals supplement sensitivity and specificity estimates. The one-control cohort's AUROC interval was withheld from the primary presentation rather than displayed as an apparently precise bootstrap interval.

Paired DeLong tests compared inductive scoring rules for the same independent patients within each held-out primary study [7]. Target-batch-fitted ComBat adaptation is not an individual fixed scoring rule, so its ordinary DeLong comparisons were explicitly unavailable; descriptive performance estimates were retained. BH correction covered finite tests within the primary comparison family; the secondary three-study family was separate [8]. Pooled predictions across different fitted models did not receive a pooled DeLong test. With two primary directions, center-level tests cannot establish broad population-center inference. AUPRC is interpreted together with case proportion because its baseline changes with prevalence.

### 2.7. Prevalence and score-comparability diagnostics

A contingency-table chi-square test described the study/diagnosis association. To separate a low pooled AUROC from study-wise ranking, we decomposed all cancer–healthy patient pairs in the secondary ElasticNet out-of-fold predictions into within-study and cross-study terms. Their weighted concordance contributions sum exactly to pooled AUROC. Cross-study terms compare scores produced by different source-trained models, so model offset or calibration differences can affect them even when ranking inside each target remains informative.

An additional descriptive check expressed each score as logit(probability) minus logit(source-locked threshold). This used no target labels to align the scoring scales and was not substituted for the primary predictions or interpreted as a calibrated probability. It tests how a source-defined operating-point reference changes pooled ranking across separately fitted rules.

The composition stress experiment trained on the two non-target studies and froze both preprocessing state and ElasticNet weights. For each eligible target study, five healthy and five cancer anchors were held fixed. We constructed target batches of 32 distinct real patients at cancer fractions 0.25, 0.50 and 0.75, using source-independent deterministic selection and no synthetic feature vectors. Target labels were used only to define this post-evaluation diagnostic composition; those labels did not tune any benchmark prediction. We measured the maximum change in an anchor's probability. The experiment tests companion-dependence of a fitted transformation, rather than estimating a universal causal effect of population prevalence.

Finally, within-target window Cohen d was compared before and after harmonization. A positive feature-wise affine map preserves this standardized mean difference within a study. This provides a direct check against attributing all low pooled performance to destruction of univariate biological contrasts. Chromosome-arm averaged fragmentation contrasts were also calculated; absent independently measured copy number, they were not relabelled as confirmed copy-number signals.

### 2.8. Code assistance and figure provenance

OpenAI Codex (OpenAI; desktop application used in October 2026; exact application build and underlying model version not recorded) assisted implementation and refactoring, testing, statistical-workflow code and figure-generation code. The author checked the computational logic, statistical outputs and interpretation. All quantitative panels were generated from the recorded real-data outputs using the supplied Python plotting scripts. The explanatory workflow in Figure 3a was drawn programmatically with the same code assistance and reviewed by the author.

## 3. Results

### 3.1. Study-associated drift and the need to distinguish scales

The analysis cohort comprised 585 participants. In raw-ratio PERMANOVA, the unique study fraction was 65.628% and the disease fraction was 0.833% (both permutation p=0.001). GC correlated with raw window study eta-squared at rho=0.606, p=3.21e-49. For Global Z-score, the correlation was 0.158, p=0.000525; for local contrasts, it was 0.405, p=3.26e-20. Median unique window study variance fractions were 65.545%, 2.553% and 2.456%, respectively. Normalization therefore changed both the amount of source-associated variation and its genomic ordering (Figure 1; Figure 4).

Median absolute arm-averaged disease contrasts were 0.381 for raw ratios, 0.270 for Global Z-score and 0.244 for local contrasts. The reference-region comparisons yielded 0 BH-significant enrichments among 7 tests. These results quantify how normalization redistributes broad fragmentation contrasts and provide the measured context for the multiscale comparison.

The Scale Trade-off has two observable dimensions: broad fragmentation-contrast retention and the magnitude and genomic pattern of residual study variation. Local normalization produced a slightly smaller median study-variance fraction than Global Z-score, alongside a stronger residual GC rank association. Thus the preferred representation depends on the property being evaluated. Copy-number and tissue-accessibility interpretations are considered as mechanistic hypotheses in the Discussion.

### 3.2. Primary bilateral generalization and operating-point transfer

All primary training and calibration excluded the external high-burden cohort. In the Jiang-to-Cristiano direction, Hybrid AUROC was 0.608 and AUPRC was 0.671; in the Cristiano-to-Jiang direction they were 0.877 and 0.952. Hybrid discrimination relative to Global Z-score was different in direction between the two targets. The DeLong p values versus Global Z-score were 0.001 for the Cristiano target and 0.083 for the Jiang target; BH values were 0.003 and 0.144, respectively.

At the source-locked operating threshold, Hybrid specificity was 0.632 in Cristiano and 0.531 in Jiang, corresponding to deviations of -31.842 and -41.875 percentage points from the source target. Target sensitivities were 0.514 and 0.910. The full direction-specific primary comparison is:

| pipeline          | center                |   auroc |   auroc_low |   auroc_high |   auprc |   sensitivity95 |   specificity95 |   specificity_deviation_pp |   brier |   ece |
|:------------------|:----------------------|--------:|------------:|-------------:|--------:|----------------:|----------------:|---------------------------:|--------:|------:|
| Global_Zscore     | Cristiano et al, 2019 |   0.632 |       0.576 |        0.682 |   0.694 |           0.555 |           0.660 |                    -28.971 |   0.297 | 0.241 |
| Global_Zscore     | Jiang et al, 2015     |   0.859 |       0.783 |        0.925 |   0.946 |           0.944 |           0.406 |                    -54.375 |   0.171 | 0.166 |
| Hybrid_MultiScale | Cristiano et al, 2019 |   0.608 |       0.553 |        0.658 |   0.671 |           0.514 |           0.632 |                    -31.842 |   0.316 | 0.255 |
| Hybrid_MultiScale | Jiang et al, 2015     |   0.877 |       0.807 |        0.939 |   0.952 |           0.910 |           0.531 |                    -41.875 |   0.151 | 0.155 |
| Proposed_LSN_ICC  | Cristiano et al, 2019 |   0.607 |       0.550 |        0.659 |   0.675 |           0.619 |           0.512 |                    -43.804 |   0.327 | 0.281 |
| Proposed_LSN_ICC  | Jiang et al, 2015     |   0.816 |       0.730 |        0.890 |   0.925 |           0.753 |           0.812 |                    -13.750 |   0.159 | 0.134 |
| Raw_Ratio         | Cristiano et al, 2019 |   0.725 |       0.678 |        0.771 |   0.782 |           0.794 |           0.450 |                    -50.024 |   0.279 | 0.208 |
| Raw_Ratio         | Jiang et al, 2015     |   0.759 |       0.671 |        0.837 |   0.917 |           0.326 |           1.000 |                      5.000 |   0.231 | 0.248 |
| Tuned_LSN_ICC     | Cristiano et al, 2019 |   0.622 |       0.566 |        0.674 |   0.681 |           0.619 |           0.574 |                    -37.584 |   0.306 | 0.247 |
| Tuned_LSN_ICC     | Jiang et al, 2015     |   0.798 |       0.717 |        0.875 |   0.925 |           0.708 |           0.812 |                    -13.750 |   0.181 | 0.155 |

The pooled primary Hybrid AUROC was 0.689 and its equal-weight mean study AUROC was 0.743; corresponding Global Z-score values were 0.702 and 0.745. These summaries evaluate different weighting and score-comparability properties. Figure 2 and Figure 3 separate direction-specific discrimination, calibration and locked-threshold transport. The complete classifier table reports estimable and unavailable combinations explicitly.

### 3.3. Nested local selection and the external sensitivity cohort

The nested-tuned local representation yielded pooled primary AUROC 0.676, AUPRC 0.770, Brier score 0.278 and ECE 0.225. The fixed local fallback yielded AUROC 0.670. Outer source-selected settings were:

| cohort_scope        |   fold |   training_n | test_study            |   selected_radius |   selected_alpha | icc_bypassed   |
|:--------------------|-------:|-------------:|:----------------------|------------------:|-----------------:|:---------------|
| Primary_Bilateral   |      0 |          121 | Cristiano et al, 2019 |                 8 |             0    | True           |
| Primary_Bilateral   |      1 |          427 | Jiang et al, 2015     |                 2 |             0    | True           |
| External_HighBurden |      0 |          548 | Snyder et al, 2016    |                 5 |             0    | False          |
| Combined_LOCO       |      0 |          158 | Cristiano et al, 2019 |                 8 |             0.15 | False          |
| Combined_LOCO       |      1 |          464 | Jiang et al, 2015     |                 2 |             0    | False          |

In the external high-burden cohort, all thresholds and models were locked using the primary source studies. Hybrid detected a fraction 0.944 of cancer cases; its exact binomial 95% interval was 0.813–0.993. The external operating-point results were:

| pipeline          |   sensitivity95 |   sensitivity_exact_low |   sensitivity_exact_high |   specificity95 |   source_threshold_min |
|:------------------|----------------:|------------------------:|-------------------------:|----------------:|-----------------------:|
| Global_Zscore     |           0.944 |                   0.813 |                    0.993 |           0.000 |                  0.771 |
| Hybrid_MultiScale |           0.944 |                   0.813 |                    0.993 |           0.000 |                  0.777 |
| InFold_ComBat     |           0.028 |                   0.001 |                    0.145 |           1.000 |                  0.658 |
| Proposed_LSN_ICC  |           0.917 |                   0.775 |                    0.982 |           0.000 |                  0.914 |
| Raw_Ratio         |           0.972 |                   0.855 |                    0.999 |           0.000 |                  0.567 |
| Tuned_LSN_ICC     |           0.972 |                   0.855 |                    0.999 |           0.000 |                  0.885 |

The observed healthy outcome in this cohort was retained in the numerical table, but it was not represented as a stable population specificity estimate in Figure 2c. The secondary Combined LOCO results across all independently identified patients were:

| pipeline          |   auroc |   auroc_low |   auroc_high |   auprc |   sensitivity95 |   specificity95 |   brier |   ece |
|:------------------|--------:|------------:|-------------:|--------:|----------------:|----------------:|--------:|------:|
| Global_Zscore     |   0.759 |       0.721 |        0.791 |   0.846 |           0.636 |           0.802 |   0.243 | 0.203 |
| Hybrid_MultiScale |   0.731 |       0.694 |        0.767 |   0.832 |           0.644 |           0.674 |   0.255 | 0.207 |
| InFold_ComBat     |   0.543 |       0.508 |        0.577 |   0.707 |           0.397 |           0.855 |   0.266 | 0.225 |
| Proposed_LSN_ICC  |   0.724 |       0.685 |        0.761 |   0.819 |           0.592 |           0.719 |   0.261 | 0.215 |
| Raw_Ratio         |   0.555 |       0.522 |        0.591 |   0.732 |           0.714 |           0.174 |   0.308 | 0.284 |
| Tuned_LSN_ICC     |   0.738 |       0.697 |        0.775 |   0.835 |           0.510 |           0.913 |   0.246 | 0.189 |

For fixed ElasticNet, Combined LOCO AUROC was 0.759 for Global Z-score, 0.731 for Hybrid and 0.724 for Proposed LSN-ICC. Hybrid therefore occupied an intermediate position between the global and fixed local representations. In the primary paired comparison, its AUROC was lower than Global in Cristiano and numerically higher in Jiang; the latter comparison did not reach BH significance. These direction-specific results support the scale-complementarity audit rather than a claim of consistent Hybrid superiority. DeLong inference remains target-specific.

### 3.4. Why pooled ComBat discrimination and within-study ranking diverge

Cancer proportions were 51.054% in Cristiano, 73.554% in Jiang and 97.297% in Snyder. Study and diagnosis were associated (chi-square=44.024, df=2, p=2.76e-10, Cramer's V=0.274), establishing the composition contrast examined in the adaptation diagnostics.

The empirical-Bayes-plus-target-alignment comparator had secondary pooled AUROC 0.543. Pair decomposition separated weighted within-study concordance of 0.669 from cross-study concordance of 0.366; cross-study pairs contributed 41.636% of all pairs. The pooled AUROC therefore combines study-wise ranking with comparability between scores produced by independently source-trained models (Supplementary Figure S5).

The maximum change in within-target window Cohen d under this comparator was 3.04e-08. Positive feature-wise affine harmonization preserved those standardized contrasts to numerical precision. The classifier and pooled operating scores can nevertheless change because target-batch location and scale alter the relationship to the source-trained rule. This distinguishes information retained within a target study from calibration and feature weighting learned elsewhere.

Referencing each ComBat score to its source-locked threshold yielded secondary pooled AUROC 0.616, a change of +0.073. This descriptive alignment preserved the benchmark predictions and used source-defined operating points to examine comparability across fitted score scales.

The fixed-anchor experiment directly identified target-composition dependence. Across batches of real target patients with different case fractions, the maximum anchor probability range was 0.214 under target-batch adaptation. Corresponding ranges were 0 for Raw Ratio, 0 for Global Z-score, 0 for fixed LSN-ICC and 0 for Hybrid. Within-sample transformations thus removed companion-dependent target normalization for the frozen scoring rules (Figure 3d).

## 4. Discussion

This audit separates three properties that are often combined in cross-study cfDNA evaluation: normalization scale, study-wise discrimination and comparability across fitted score scales. Global and local transformations change broad fragmentation contrasts and residual source variation in different ways. Pair decomposition explains how pooled AUROC combines within-study ranking with cross-study score comparability, while fixed real-patient anchors isolate companion dependence in unseen-target adaptation. Together, these analyses turn a single benchmark score into an interpretable account of representation and transport.

Global and local normalization impose different reference scales. Global Z-scoring retains deviations relative to the full genome; the local log contrast references adjacent windows on the same arm. Attenuation of a smooth background can therefore accompany attenuation of a broad disease-associated contrast. Hybrid exposes both scales to a source-trained penalized classifier, making their complementarity testable. The observed target-specific comparisons identify the conditions in which the fused representation changes ranking and calibration; the Combined LOCO estimates favour Global Z-score over Hybrid for the fixed ElasticNet comparison.

The ComBat diagnostics locate a specific source of composition dependence. The evaluated system combines empirical-Bayes source adjustment and robust unseen-target mapping. A target quantile follows the mixture of patients, so changing that mixture can move the adapted features of a fixed individual. The anchor experiment measures the resulting probability change with all source fitting held constant. Preserved within-study Cohen d and the pair decomposition distinguish this adaptation effect from removal of univariate disease contrasts. The conclusion applies to the evaluated adaptation rule; other covariate-aware or reference-batch ComBat designs remain separate methods [4].

Within-sample transformation makes a patient's representation independent of other target patients. For a frozen source model, adding or removing companions leaves the global, local and hybrid feature vectors and predictions unchanged. This property removes a defined route of target-composition dependence. Training-population composition, disease spectrum and prevalence-dependent predictive values remain separate determinants of performance, making calibration and source-locked threshold transport complementary evaluation targets.

The cohort hierarchy aligns each public source with the endpoint it can support. The two primary targets permit bidirectional discrimination and specificity estimation. The separately held-out, cancer-enriched source contributes sensitivity at primary-source-locked thresholds. Combined LOCO then evaluates the three-source setting with its score-comparability requirements made explicit. This structure provides a reusable evaluation design for public datasets with unequal case fractions and control representation.

### 4.1. Limitations

The available source-level replication comprises 2 cohorts with substantial healthy-control representation and a case-dominated external source (36 cancer/1 healthy). The primary sources contain both classes but are not exactly balanced. Adalsteinsson had 528 known cancer participants and 0 explicitly healthy participants; Sun had 12 and 0. The archived local inventory yielded no additional source meeting the eligibility criterion. This inventory bounds the present analysis and is distinct from a current exhaustive database census.

Study Source combines experimental-batch and population-composition differences. Each bilateral fit has one training study, making multi-source ICC filtering, multi-batch source ComBat and multi-source adversarial training unidentifiable in that protocol. The LSN-only fallback and unavailable method cells are recorded explicitly. The external cohort's single healthy participant supports a descriptive outcome, while generalization inference remains confined to the public-cohort evaluation.

The analysis is methodological development on the same public participants that motivated the representation comparison. Fold-level patient isolation addresses algorithmic leakage; independent confirmatory evaluation remains a separate requirement. Bootstrap intervals condition on the fixed studies and fitted scores and exclude refitting, method-development history and future-study uncertainty. DeLong tests are target-specific. The composition stress experiment uses a finite set of real-patient batches to measure adaptation dependence under the tested compositions.

Matched copy-number measurements and tumor-specific chromatin assays were unavailable for the complete cohort. Broad arm-level fragmentation contrasts and reference GC associations therefore remain descriptive endpoints rather than direct CNA or accessibility measurements. Source metadata also have missing clinical covariates. Further source-level replication with balanced control recruitment and prospective specimen handling would extend the present cross-study conclusions.

## 5. Conclusions

Normalization scale, source composition and score comparability define distinct axes of cross-study cfDNA evaluation. The Scale Trade-off audit quantifies broad fragmentation contrasts and residual genomic source variation, while pair decomposition connects study-wise ranking to pooled discrimination. Fixed-patient anchors identify companion-dependent target adaptation and demonstrate its absence for within-sample transformations. This framework makes the behaviour of normalization and locked thresholds directly testable across public studies.

## Data and code availability

Source data are described in Cristiano et al. (2019), Jiang et al. (2015) and Snyder et al. (2016). Fragment intervals were accessed through FinaleDB (https://finaledb.research.cchmc.org/); source-linked metadata are documented in UNITE supplementary data file S1, sheet s2 (https://doi.org/10.1126/sciadv.ady9432). Source URLs, eligibility rules and provenance are recorded in the extraction configuration and manifests.

Code and compact statistical outputs are distributed at https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit. Zenodo registration is in progress for the v0.1.0 software release. The release contains the MIT-licensed original software, frozen aggregate statistical tables, tests and compact benchmark provenance. Large participant feature matrices, raw fragment streams and fitted-model caches are excluded from the lightweight distribution. Streaming scripts reconstruct the derived inputs from the published sources subject to availability and their original data-use terms. Third-party data and annotations retain those terms.

## CRediT authorship contribution statement

Zhenxuan Ma: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Resources, Data curation, Writing – original draft, Writing – review & editing, Visualization, Supervision, Project administration. All applicable roles were undertaken by the sole author, who approved the final manuscript.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Declaration of competing interests

The author declares that he has no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Ethics statement

This study is a secondary analysis of previously published, de-identified public human plasma cfDNA fragment-interval data and associated metadata obtained through FinaleDB and the published UNITE supplementary data. No new participants were recruited, no new specimens were collected, and no attempt was made to re-identify participants. No additional institutional ethics approval or formal exemption determination was obtained for this secondary analysis. Original-study consent and ethics protections remain those of the source investigations; no new consent was sought. Public-data exemption provisions may apply subject to institutional requirements; this statement does not assert an institutional exemption decision.

## Acknowledgements

The author thanks the investigators who generated the original public datasets and the developers of FinaleDB, UNITE and the open-source scientific software used in this analysis.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work, the author used OpenAI Codex (OpenAI; desktop application used in October 2026; exact application build and underlying model version not recorded) to assist with code implementation and refactoring, unit tests, statistical-workflow and figure-generation code, reference verification, manuscript organization and drafting, and language revision. The author reviewed the code logic, statistical outputs and text, edited the material as needed, approved the manuscript, and takes full responsibility for its accuracy, integrity and interpretation. The AI tool is not an author. Figures were generated by reproducible code from the recorded data; no clinical observations were fabricated.

## Figure legends

**Figure 1. Study-associated baseline variation in the independently identified cfDNA cohort.** (a) Study composition and observed healthy fragment-frequency curves, with reference length markers. (b) Descriptive UMAP coloured by study and diagnosis, annotated with coverage-adjusted conditional PERMANOVA fractions and permutation results. (c) Unique study and disease window variance fractions across the autosomes. (d) Healthy baseline distributions and one-way healthy-center clustering scores. The external study's single healthy person is descriptive. No embedding was supplied to predictive preprocessing.

**Figure 2. Hierarchical reporting of primary transport and external cancer sensitivity.** (a) Random and bilateral study-held-out macro AUROC for Raw Ratio and Hybrid across fixed classifiers, using the same primary patient cohort; training-set sizes differ between the schemes. (b) The bilateral macro-AUROC matrix; unavailable ComBat, DANN and unrequested tuned-classifier cells remain grey. Asterisks identify LSN fits whose healthy-center filter cannot be estimated from one source study. (c) ElasticNet target specificity minus the source target in percentage points, with exact binomial intervals for the two adequately represented healthy target groups. (d) Cancer sensitivity in the separately held-out high-burden cohort, using primary-source-locked thresholds and exact binomial intervals. The lone healthy outcome is reported numerically and does not generate an unstable specificity error bar in panel c.

**Figure 3. Multiscale representation, primary discrimination and target-composition dependence.** (a) Sample-wise global and same-arm local branches, training-only scaling, concatenation and ElasticNet classification; bars show the observed two directional AUROCs. (b) Primary pooled held-out ROC and precision–recall curves, with target-specific DeLong comparisons against Global Z-score. Pooled curves are descriptive across different source-trained models. (c) Primary calibration curves and Brier scores. (d) Maximum change in a fixed real anchor patient's probability when the target batch case fraction changes, with all source fitting frozen. Diagnostic composition labels never entered benchmark fitting or selection. The explanatory workflow in panel a was drawn with code assisted by OpenAI Codex (OpenAI; October 2026 desktop session; exact build and model version not recorded) and checked by the author.

**Figure 4. Reference genomic context and fixed-pipeline feature retention.** (a) Raw-feature study variance associations with reference GC and telomere distance. (b) Reference-region odds ratios for majority-fold filtering under the previously fixed LSN-ICC protocol; corrected results are reported in the accompanying table. (c) Observed chromosome trajectories before and after local contrast, with gaps across removed windows. (d) Fixed-pipeline coefficient rank concordance across the original three-study LOCO fits. These panels contextualize the fixed local representation and are not relabelled as Hybrid feature importance or tumor-specific accessibility measurements.

**Supplementary Figure S5. ComBat pooled discrimination and score comparability.** Within-study and cross-study cancer–control pair contributions sum to secondary pooled AUROC. The paired summary shows pooled discrimination and the equal-weight mean over the two primary held-out targets, emphasizing their different score-comparability requirements.

## Highlights

- Normalization scale shapes cross-study cfDNA discrimination and calibration.
- Pair concordance separates study-wise ranking from pooled score comparability.
- Fixed patient anchors expose composition-dependent batch adaptation.
- Within-sample transforms preserve anchor scores across target compositions.
- Source-locked thresholds reveal operating-point drift across public studies.

## References

1. Cristiano S, Leal A, Phallen J, et al. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature. 2019;570(7761):385-389. https://doi.org/10.1038/s41586-019-1272-6

2. Jiang P, Chan CWM, Chan KCA, et al. Lengthening and shortening of plasma DNA in hepatocellular carcinoma patients. Proceedings of the National Academy of Sciences. 2015;112(11):E1317-E1325. https://doi.org/10.1073/pnas.1500076112

3. Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. Cell. 2016;164(1-2):57-68. https://doi.org/10.1016/j.cell.2015.11.050

4. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8(1):118-127. https://doi.org/10.1093/biostatistics/kxj037

5. Zheng H, Zhu MS, Liu Y. FinaleDB: a browser and database of cell-free DNA fragmentation patterns. Bioinformatics. 2021;37(16):2502-2503. https://doi.org/10.1093/bioinformatics/btaa999

6. Wang H, Mennea PD, McAndrew G, et al. A scalable deep-learning framework for cancer detection using cell-free DNA shallow whole-genome sequencing. Science Advances. 2026;12(28):eady9432. https://doi.org/10.1126/sciadv.ady9432

7. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the Areas under Two or More Correlated Receiver Operating Characteristic Curves: A Nonparametric Approach. Biometrics. 1988;44(3):837-845. https://doi.org/10.2307/2531595

8. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society: Series B (Methodological). 1995;57(1):289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
