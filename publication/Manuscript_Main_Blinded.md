# Scale Trade-offs in Within-Sample Self-Normalization for Cross-Study Generalization in Public cfDNA Cohorts

## Abstract

**Background:** Within-sample cfDNA normalization makes target preprocessing independent of companion patients. Its implications for genomic scale, discrimination and threshold transport require separate evaluation.

**Methods:** We analysed 585 source-identified participants, prioritizing a bidirectional study-held-out benchmark of 548. Source-only calibration, paired patient bootstraps and repeated real-patient companion resampling evaluated operating points and score dependence; a separate case-dominated cohort contributed sensitivity evidence.

**Results:** In Cristiano, Global and Hybrid AUROC/AUPRC were 0.632/0.694 and 0.608/0.671; Hybrid−Global ΔAUROC was -0.024 (95% CI [-0.039, -0.008]; BH q=0.003). In Jiang, corresponding AUROC/AUPRC were 0.859/0.946 and 0.877/0.952, with ΔAUROC +0.018 ([-0.003, 0.040]; q=0.144). Equal-weight AUROCs were 0.742 for Raw, 0.745 for Global and 0.743 for Hybrid. Source thresholds calibrated to 95% specificity yielded target specificities 0.632–0.660 in Cristiano and 0.406–0.531 in Jiang. At fixed target case fraction 0.50, EB-TargetMap companion resampling produced decision-flip rates 0.140 and 0.126.

**Conclusions:** Companion independence is a mathematical property of within-sample transformations. Normalization-scale discrimination is target-dependent, and companion independence alone does not ensure source-threshold transport.

**Keywords:** cell-free DNA; fragmentomics; normalization scale; study-held-out validation; threshold transport; companion dependence

## 1. Introduction

Cross-study cfDNA fragmentomics requires representations that retain disease-associated fragmentation contrasts while controlling study-source variation and maintaining interpretable scores in a new cohort. Plasma DNA fragmentation reflects nucleosome organization and tissue contribution, and genome-wide profiles support cancer classification alongside sequence-variant assays [1,2,3,4]. Fragment-size selection, transcription-associated coverage and lung-cancer fragmentome studies establish complementary uses of these signals [5,6,7]. Doebley et al. further link nucleosome protection profiles with cancer phenotypes while explicitly addressing GC bias [8]. Study Source combines specimen processing, sequencing procedures, recruitment and population composition. This joint structure is consistent with the broader challenge of batch effects in high-throughput molecular data [9,10]. Evaluating transport therefore requires both patient-level discrimination and the achieved operating point of a threshold locked in the source study.

Global and local normalization define different genomic references. A genome-wide reference preserves broad relative deviations, while a same-arm local contrast emphasizes departures from nearby windows. Combining the scales makes their complementary information testable. A useful evaluation must distinguish the disease contrasts preserved by a representation from the study on which a classifier is tested and the threshold at which it is operated.

Cross-sample alignment adds dependence on the target patient mixture. Empirical-Bayes source harmonization estimates study-level location and scale [11]; an additional unseen-target median/IQR map uses the distribution of the patients currently available. Under $F_c=(1-\pi_c)F_{0c}+\pi_cF_{1c}$, both case fraction and individual companion identities can change these summaries. Within-sample transformations provide an exact comparison because the individual supplies their own target reference.

We evaluate normalization scale through bidirectional discrimination, source-locked threshold transfer and companion dependence. A paired analysis identifies target-dependent differences between Global and Hybrid. Calibration experiments distinguish uncertainty in the source threshold from the specificity actually achieved in the target. Repeated fixed-anchor batches separate changes in case proportion from companion-identity variability. Genomic block inference and multivariate dispersion diagnostics place study-associated variation in its spatial and distributional context.

## 2. Materials and methods

### 2.1. Cohort and source provenance

Published UNITE supplementary metadata and FinaleDB fragment intervals supplied the clinical and technical inputs [12,13]. Plasma WGS specimens required an explicit cancer or Healthy source label, original sequencing depth 0.1–5.0× and the recorded extraction QC. From 628 extracted specimens, one prespecified specimen per known participant was selected by earliest time point, shallowest depth and run identifier. The unresolved-identity specimen was excluded, leaving 585 source-identified participants. Undocumented cross-study aliases cannot be reconciled from the source inventory.

Cristiano and Jiang constituted the Primary Bilateral Benchmark (548 participants). Each direction trained on one entire study and tested the other. Snyder was excluded from primary training, selection and calibration, providing a separate external sensitivity cohort (37 participants, 36 cancer/1 healthy). Secondary leave-one-center-out cross-validation (LOCO-CV), with study source as the grouping variable, used all three studies. Study Source is a proxy combining experimental handling, recruitment and population composition.

**Table 1. Independent participant cohort.**

| Study                 | Participants   | Cancer   | Healthy   |
|:----------------------|:---------------|:---------|:----------|
| Cristiano et al, 2019 | 427            | 218      | 209       |
| Jiang et al, 2015     | 121            | 89       | 32        |
| Snyder et al, 2016    | 37             | 36       | 1         |

**Table 2. Study-specific sequencing depth, fragment counts, and cancer-type composition.**

| Characteristic              | Cristiano                                                                                                       | Jiang                  | Snyder                                                                                                                                                                                      |
|:----------------------------|:----------------------------------------------------------------------------------------------------------------|:-----------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| n (cancer/healthy)          | 427 (218/209)                                                                                                   | 121 (89/32)            | 37 (36/1)                                                                                                                                                                                   |
| Depth median [Q1,Q3] ×      | 3.591 [3.175,4.017]                                                                                             | 3.332 [2.787,3.839]    | 4.089 [3.385,4.712]                                                                                                                                                                         |
| Fragments median [Q1,Q3], M | 33.304 [29.328,37.250]                                                                                          | 30.927 [25.800,35.451] | 28.910 [24.698,32.550]                                                                                                                                                                      |
| Cancer stage                | Stage 0: 2; Stage I: 44; Stage II: 102; Stage III: 34; Stage IV: 35; Stage X: 1; Unknown/Missing: 0; Total: 218 | Not recorded           | IV: 36                                                                                                                                                                                      |
| Cancer types (n)            | Breast: 51; Lung: 35; Pancreatic: 31; Ovarian: 27; Gastric: 26; Colorectal: 24; Bile Duct: 23; Duodenal: 1      | Hepatocellular: 89     | Lung: 7; Breast: 5; Pancreatic: 3; Melanoma: 2; Uterine: 2; Ovarian: 2; Colorectal: 2; Head And Neck: 2; Renal: 2; Bladder: 2; Testicular: 2; Hepatocellular: 2; Esophageal: 2; Prostate: 1 |

Fragment counts are accepted autosomal 80–400 bp intervals after QC. Stage categories are preserved as supplied, including 0 and X; Jiang stage is unrecorded.

### 2.2. Features and source-fitted preprocessing

The hg19/GRCh37 reference contains 477 complete autosomal 5 Mb windows after fixed centromere, terminal-risk, cytoband-risk and blacklist exclusions. Autosomal fragments with MAPQ≥30 and length 80–400 bp contribute to the global integer-length distribution; midpoint assignment supplies window counts. Short and long counts cover 100–150 and 151–220 bp:

$$R_{si}=\frac{S_{si}}{L_{si}+1}.$$

Global Z-score references the individual's window mean and standard deviation. Local self-normalization (LSN) references other retained windows on the same chromosome arm, up to three on each side:

$$G_{si}=\frac{R_{si}-\overline{R}_{s}}{\operatorname{sd}_{i}(R_{si})}.$$

$$\mathrm{LSN}_{si}=\log_2\frac{R_{si}+\epsilon_s}{\operatorname{median}_{j\in N_3(i),j\ne i}R_{sj}+\epsilon_s}.$$

The relative offset is $\epsilon_s=10^{-6}\operatorname{median}_iR_{si}$.

Hybrid concatenates $[G_s,\mathrm{LSN}_s]$ before source-fitted column scaling. Raw, Global and Hybrid use StandardScaler; Local uses RobustScaler. All imputers and scalers fit source-training patients exclusively, including inside each calibration fold. Per-column source medians replace post-transform NaNs; an entirely missing source column has zero as its documented fallback. The ratio and length matrices each had 0 missing cells; the ratio matrix contained 279045 cells. Feature matrices use float32. In multiple-source fits, a one-way healthy-center intraclass correlation (ICC) clustering score ranks windows, and Local_LSN_ICC removes the highest quarter. In one-source bilateral fits, this score is unidentifiable and the representation is precisely Local_LSN, without an ICC estimate.

Elastic-net logistic regression used $C=1.0$, `l1_ratio=0.5`, solver saga, tolerance 0.001 and maximum 20000 iterations, with seed 2026 [14,15]. The fixed Hybrid radius was specified before outer evaluation. Source-only nested local selection evaluated radii 1, 2, 3, 5 and 8 and exclusion fractions 0, 0.15, 0.25, 0.35 and 0.50 as a supplementary analysis. A single healthy source restricts exclusion to zero. Three-fold source-only mean eligible-study AUROC selected the settings; ties favored lower exclusion, radius nearer three, then smaller radius. Random Forest, CUDA XGBoost, LightGBM, MLP and domain-adversarial neural network comparisons are supplied in the supplement [16,17,18]. NumPy and UMAP supported numerical computation and descriptive embeddings [19,20]. Exact package versions are listed in Supplementary Table S2.

### 2.3. Mathematical property and the EB-TargetMap comparator

For frozen source parameters $\theta$, an inductive sample-wise rule has $p_i=g_\theta(f_\theta(x_i))$. For any companion batches $B$ and $B'$ containing $i$,

$$p_i(B)-p_i(B')=0.$$

This follows because neither function uses companion feature vectors. Automated tests vary companion identity and scale while checking predictions and serialized fitted state. The property concerns target preprocessing with a fixed source model; it does not imply invariant training populations, discrimination or operating thresholds. The local ratio's relative epsilon also preserves global multiplicative scaling to numerical precision, without extending this guarantee to arbitrary additive or nonlinear distortions.

EB-TargetMap denotes the specific empirical-Bayes source harmonization plus unlabeled target median/interquartile-range (IQR) map used here. Source harmonization contains study indicators and no disease covariate. With pooled source mean $\mu_j$, residual scale $\sigma_j$ and empirical-Bayes source parameters $(\gamma_{cj},\delta_{cj})$, source features are

$$x^{EB}_{ij}=\mu_j+\sigma_j\frac{(x_{ij}-\mu_j)/\sigma_j-\gamma_{cj}}{\sqrt{\delta_{cj}}}.$$

For an unseen target batch $B$, the rule maps raw imputed features to the harmonized-source reference:

$$x^{map}_{ij}=\frac{x_{ij}-\operatorname{median}_{k\in B}x_{kj}}{\max(\operatorname{IQR}_{k\in B}x_{kj},10^{-6})}\operatorname{IQR}_{source,j}^{EB}+\operatorname{median}_{source,j}^{EB}.$$

The source IQR is also floored at $10^{-6}$. The final StandardScaler and classifier use harmonized source patients only. Target labels do not enter the map, but companion features do. EB-TargetMap is a particular transductive adaptation rule, rather than a label for all ComBat variants.

### 2.4. Calibration and statistical comparisons

Three-fold mixed-source out-of-fold (OOF) predictions supplied source control scores; rare study/diagnosis strata remained training-only. The higher empirical 0.95 quantile fixed the threshold, and scores strictly above it were positive. Tuned local selection was repeated inside calibration-training folds. The source model was then refitted on all source patients and applied unchanged to the target. In multiple-source experiments, a second calibration scheme left one entire source study out at a time. Each calibration model fit the remaining source data, including its imputer, feature filter and scaler. EB-TargetMap cannot estimate source harmonization in these one-source calibration-training folds, so its source-study-held-out calibration is unavailable. Local_LSN_ICC correspondingly uses its documented Local_LSN branch inside such folds.

We resampled held-out source control scores 1000 times to estimate threshold uncertainty and propagated these thresholds through the same fixed target scores. These intervals condition on fitted calibration models and target patients. They isolate calibration-score sampling uncertainty rather than refitting or future-study uncertainty. Source-study-held-out calibration changes the validation domain and training size; its comparison with mixed-source OOF does not isolate either component causally.

Primary outcomes were target-specific AUROC, area under the precision–recall curve (AUPRC), equal-weight mean target AUROC and achieved specificity/sensitivity at the locked source threshold. Secondary measures included Brier score and 10-bin expected calibration error (ECE) [21,22]. Patient bootstrap intervals and paired ΔAUROC intervals used 1000 diagnosis-stratified resamples conditional on fixed scoring rules. Exact binomial intervals accompany operating-point counts. Snyder specificity is reported as a fraction with denominator one, solely descriptively.

DeLong tests compare fixed individual scoring rules within each independent held-out primary study [23]. The retained primary Benjamini–Hochberg (BH) family has 38 finite tests: Hybrid versus Raw, Global and Local for each of five estimable classifiers in both targets (30), plus Tuned Local versus Raw, Global, Local and Hybrid for ElasticNet in both targets (8) [24]. Supplementary Table S1 enumerates every hypothesis. Unavailable EB-TargetMap and DANN cells contribute no p value. Secondary families are separate; no pooled cross-model DeLong test is asserted. Source isolation addresses algorithmic leakage [25].

### 2.5. Companion resampling and cohort sensitivity analyses

For each adequately represented target, the other two studies trained a frozen ElasticNet scoring rule. Five healthy and five cancer anchors were held fixed. At target case fractions 0.25, 0.50 and 0.75, each diagnostic batch contained 32 distinct real patients; 50 independent companion priority draws per fraction generated the batch series. Matched priorities coupled the fractions within a replicate. Outcomes defined diagnostic compositions after evaluation and selected neither predictive settings nor benchmark thresholds. At fixed fraction 0.50, identity resampling measures sampling variability independently of changing the fraction. Probability shifts reference each anchor's prediction in replicate 0 at fraction 0.50. We report absolute-shift median, quartiles, IQR, full patient-level distributions, within-fraction ranges and decision flips relative to this same reference. Repeated anchor observations are dependent and do not provide independent patient replication.

Post-hoc exploratory sensitivity analyses restricted both primary studies to a fixed 2–4× original-depth interval, then repeated source-only fitting and calibration. This is a common-range restriction, not exact distribution matching. The post-hoc exploratory GI-versus-HCC analysis retained all controls and restricted Cristiano cancer to Bile Duct, Colorectal, Duodenal, Gastric or Pancreatic labels; Jiang remained HCC. No HCC-labelled Cristiano participant was available, so this gastrointestinal restriction evaluates a narrower disease spectrum without claiming a matched HCC comparison. These post-hoc exploratory sensitivity analyses reuse the same public cohorts and change target sample size and composition. They explore the potential direction of confounding and do not constitute independent external validation or confirmatory evidence.

### 2.6. Spatial and multivariate diagnostics

Coverage-adjusted PERMANOVA used Euclidean geometry and 999 restricted Freedman–Lane residual permutations [26,27]. Throughout, unique Type III $R^2=(SS_{reduced}-SS_{full})/SS_{total}$ denotes the drop-one contribution conditional on the other covariates, not partial $R^2$. Window eta-squared uses the analogous unique increment. For the descriptive dispersion diagnostic, columns were standardized across all audit participants, without contributing to prediction fitting. Group-centroid distances, bias-adjusted by $\sqrt{n_g/(n_g-1)}$, supplied an additional PERMDISP test with 999 unrestricted permutations of intercept-only residuals of the fixed centroid distances [28]. This dispersion diagnostic is unadjusted and does not decompose adjusted PERMANOVA into causal location and dispersion components.

Spearman associations with GC and telomere/centromere distance used a 1000-replicate whole-chromosome bootstrap over 22 autosomes. For a block-wild null, centered-rank cross-products were summed within chromosomes and their signs independently randomized in 999 replicates. This retains within-chromosome dependence under independent, sign-symmetric chromosome contributions. BH correction covered the nine representation-by-annotation tests; the four PERMDISP tests formed a separate BH family. GM12878 compartments, repeats and segmental duplications provide reference context; arm-level fragmentation effects are not independently measured CNA or tumor accessibility.

## 3. Results

### 3.1. Primary discrimination depends on target study and reference scale

Across the recorded raw-ratio geometry, unique Type III R² was 0.656 for study source (p=0.001) and 0.008 for disease status (p=0.001). These fractions describe study-associated multivariate variation and should not be interpreted as pure centroid shifts because dispersion heterogeneity was also detected. Figure 1 presents these descriptive associations alongside dispersion.

Global and Hybrid changed rank across targets (Figure 2). In Cristiano, Hybrid−Global ΔAUROC was -0.024, 95% CI [-0.039, -0.008], DeLong p=0.001042, BH q=0.003301. In Jiang, ΔAUROC was +0.018, 95% CI [-0.003, 0.040], p=0.083436, q=0.144117. Equal-weight AUROC was 0.742 for Raw, 0.745 for Global and 0.743 for Hybrid. These target-specific results identify a scale trade-off without a consistent discrimination winner.

**Table 3A. Primary bilateral fixed ElasticNet results.**

| Study     | Method            | AUROC   | AUROC 95% CI   | AUPRC   | Specificity shift (pp)   | Hybrid−Global Δ [95% CI]   | BH q   |
|:----------|:------------------|:--------|:---------------|:--------|:-------------------------|:---------------------------|:-------|
| Cristiano | Global_Zscore     | 0.632   | [0.576, 0.682] | 0.694   | -28.971                  | —                          | —      |
| Jiang     | Global_Zscore     | 0.859   | [0.783, 0.925] | 0.946   | -54.375                  | —                          | —      |
| Cristiano | Hybrid_MultiScale | 0.608   | [0.553, 0.658] | 0.671   | -31.842                  | -0.024 [-0.039, -0.008]    | 0.003  |
| Jiang     | Hybrid_MultiScale | 0.877   | [0.807, 0.939] | 0.952   | -41.875                  | +0.018 [-0.003, 0.040]     | 0.144  |
| Cristiano | Local_LSN         | 0.607   | [0.550, 0.659] | 0.675   | -43.804                  | —                          | —      |
| Jiang     | Local_LSN         | 0.816   | [0.730, 0.890] | 0.925   | -13.750                  | —                          | —      |
| Cristiano | Raw_Ratio         | 0.725   | [0.678, 0.771] | 0.782   | -50.024                  | —                          | —      |
| Jiang     | Raw_Ratio         | 0.759   | [0.671, 0.837] | 0.917   | 5.000                    | —                          | —      |

Paired difference intervals condition on the fixed target scoring rules; BH q values belong to the 38-test primary family.

### 3.2. Companion independence leaves threshold transport unresolved

At the source-locked threshold, Global/Hybrid specificity was 0.660/0.632 in Cristiano and 0.406/0.531 in Jiang. These achieved values differ substantially from the 0.950 source target despite patient-wise preprocessing independence.

For Global Zscore targeting Cristiano, the source threshold was 0.629 (conditional bootstrap 95% CI [0.407, 0.744]). Propagating source-control score uncertainty gave target specificity [0.512, 0.737] with target patients fixed.

For Hybrid MultiScale targeting Cristiano, the source threshold was 0.559 (conditional bootstrap 95% CI [0.334, 0.617]). Propagating source-control score uncertainty gave target specificity [0.517, 0.675] with target patients fixed.

For Global Zscore targeting Jiang, the source threshold was 0.729 (conditional bootstrap 95% CI [0.576, 0.881]). Propagating source-control score uncertainty gave target specificity [0.281, 0.594] with target patients fixed.

For Hybrid MultiScale targeting Jiang, the source threshold was 0.764 (conditional bootstrap 95% CI [0.624, 0.883]). Propagating source-control score uncertainty gave target specificity [0.406, 0.719] with target patients fixed.

**Table 3B. External Snyder sensitivity and descriptive single-control outcomes.**

| Method            | Source calibration    | Threshold [95% CI]   | Cancer sensitivity   | Healthy result (descriptive only)   |
|:------------------|:----------------------|:---------------------|:---------------------|:------------------------------------|
| Global_Zscore     | Mixed OOF             | 0.771 [0.641, 0.904] | 0.944                | 0/1; exact 95% CI [0.000, 0.975]    |
| Global_Zscore     | Source-study held out | 0.969 [0.940, 0.991] | 0.806                | 0/1; exact 95% CI [0.000, 0.975]    |
| Hybrid_MultiScale | Mixed OOF             | 0.777 [0.687, 0.877] | 0.944                | 0/1; exact 95% CI [0.000, 0.975]    |
| Hybrid_MultiScale | Source-study held out | 0.960 [0.926, 0.984] | 0.861                | 0/1; exact 95% CI [0.000, 0.975]    |
| EB-TargetMap      | Mixed OOF             | 0.658 [0.630, 0.678] | 0.028                | 1/1; exact 95% CI [0.025, 1.000]    |
| Local_LSN_ICC     | Mixed OOF             | 0.914 [0.837, 0.932] | 0.917                | 0/1; exact 95% CI [0.000, 0.975]    |
| Local_LSN_ICC     | Source-study held out | 0.978 [0.960, 0.986] | 0.806                | 0/1; exact 95% CI [0.000, 0.975]    |
| Raw_Ratio         | Mixed OOF             | 0.567 [0.551, 0.585] | 0.972                | 0/1; exact 95% CI [0.000, 0.975]    |
| Raw_Ratio         | Source-study held out | 0.721 [0.717, 0.728] | 0.944                | 0/1; exact 95% CI [0.000, 0.975]    |

Source-study-held-out EB-TargetMap calibration is unidentifiable because each calibration-training fold contains one study. A denominator of one supplies no stable population specificity estimate.

For external Global Zscore, source-study-held-out calibration changed the threshold from 0.771 to 0.969 and Snyder cancer sensitivity from 0.944 to 0.806. This demonstrates calibration-design dependence; single-control outcomes remain descriptive.

For external Hybrid MultiScale, source-study-held-out calibration changed the threshold from 0.777 to 0.960 and Snyder cancer sensitivity from 0.944 to 0.861. This demonstrates calibration-design dependence; single-control outcomes remain descriptive.

### 3.3. Repeated anchors distinguish case fraction from companion identity

At fixed fraction 0.50, EB-TargetMap median absolute anchor shifts were 0.016 (Q1–Q3 0.007–0.027; IQR 0.020) in Cristiano and 0.055 (0.019–0.091; IQR 0.073) in Jiang. Decision-flip rates relative to the reference batch were 0.140 and 0.126; the fractions of anchors ever changing decisions were 0.300 and 0.600. Thus companion identity alone changes adapted scores at a fixed case fraction (Figure 3).

For Cristiano, the median absolute difference between mean anchor probabilities at fractions 0.75 and 0.25 was 0.018; the maximum per-anchor range across all tested batches was 0.262. The distribution table separates this mean case-fraction contrast from within-fraction identity ranges.

For Jiang, the median absolute difference between mean anchor probabilities at fractions 0.75 and 0.25 was 0.066; the maximum per-anchor range across all tested batches was 0.692. The distribution table separates this mean case-fraction contrast from within-fraction identity ranges.

Raw, Global, Local and Hybrid preserved anchor probabilities under companion replacement to numerical tolerance, as required by the mathematical property and checked in unit tests. This invariant is a functional guarantee for frozen rules, not an empirical discrimination advantage. Full per-anchor predictions and all real batch memberships are supplied.

### 3.4. Post-hoc exploratory sensitivity analyses

These post-hoc exploratory sensitivity analyses use the same public cohorts and alter target sample size and composition. They explore potential confounding directions and provide neither independent external validation nor confirmatory evidence.

In the post-hoc exploratory sensitivity analyses within the 2–4× common-depth interval, Cristiano had 306 target participants. Global/Hybrid AUROCs were 0.663/0.632, with achieved specificities 0.795/0.627.

In the post-hoc exploratory sensitivity analyses within the 2–4× common-depth interval, Jiang had 93 target participants. Global/Hybrid AUROCs were 0.927/0.931, with achieved specificities 0.500/0.625.

In the post-hoc exploratory sensitivity analyses under the gastrointestinal-versus-HCC restriction, Cristiano had 314 target participants. Global/Hybrid AUROCs were 0.628/0.611, with achieved specificities 0.660/0.632.

In the post-hoc exploratory sensitivity analyses under the gastrointestinal-versus-HCC restriction, Jiang had 121 target participants. Global/Hybrid AUROCs were 0.879/0.879, with achieved specificities 0.562/0.719.

These post-hoc exploratory sensitivity analyses retained different source/target sizes and do not equalize stage or cancer composition. The gastrointestinal analysis contains bile duct and other digestive tumors rather than matched HCC cases. Both analyses retain the distinction between ranking and source-threshold transport (Supplementary Table S4).

### 3.5. Spatial inference and dispersion qualify study-associated variation

For Raw Ratio, the GC–unique-study-variance Spearman correlation was 0.606, chromosome-bootstrap 95% CI [0.402, 0.775], block-wild p=0.001, BH q=0.004.

For Global Zscore, the GC–unique-study-variance Spearman correlation was 0.158, chromosome-bootstrap 95% CI [-0.003, 0.352], block-wild p=0.052, BH q=0.094.

For LSN, the GC–unique-study-variance Spearman correlation was 0.405, chromosome-bootstrap 95% CI [0.307, 0.494], block-wild p=0.001, BH q=0.004.

Raw-ratio PERMDISP detected study dispersion heterogeneity (F=78.364, p=0.001) and disease dispersion heterogeneity (F=94.874, p=0.001). The adjusted unique Type III PERMANOVA study fraction remains the recorded association measure, but its significance cannot be attributed exclusively to a centroid shift. Figure 1 and Supplementary Table S5 present dispersion alongside the genomic variance audit.

The window-level audit measured median absolute disease contrasts of 0.384 for Raw and 0.263 for local normalization. Broad arm-averaged effects were 0.381, 0.270 and 0.244 for Raw, Global and Local, respectively. These fragmentation measurements describe the reference-scale trade-off; no independently measured CNA or tumor chromatin assay is inferred from them. The 7 reference-region enrichment tests yielded 0 BH-significant associations.

### 3.6. Secondary pooled scores quantify cross-model comparability

Cancer–control pair decomposition gives within-study concordance 0.669 and cross-study concordance 0.366 for EB-TargetMap, with cross-study pairs accounting for 41.636% of all pairs. Pooled three-study AUROC compares probabilities from separately fitted rules and is consequently a score-comparability diagnostic. The table retains those descriptive results separately from the primary direction-specific performance evidence (Supplementary Figure S5).

**Table 3C. Secondary Combined LOCO pooled-score comparability diagnostics.**

| Method              | Pooled AUROC   | 95% CI         | AUPRC   | Brier   | ECE   |
|:--------------------|:---------------|:---------------|:--------|:--------|:------|
| Global_Zscore       | 0.759          | [0.721, 0.791] | 0.846   | 0.243   | 0.203 |
| Hybrid_MultiScale   | 0.731          | [0.694, 0.767] | 0.832   | 0.255   | 0.207 |
| EB-TargetMap        | 0.543          | [0.508, 0.577] | 0.707   | 0.266   | 0.225 |
| Local_LSN_ICC       | 0.724          | [0.685, 0.761] | 0.819   | 0.261   | 0.215 |
| Raw_Ratio           | 0.555          | [0.522, 0.591] | 0.732   | 0.308   | 0.284 |
| Tuned_Local_LSN_ICC | 0.738          | [0.697, 0.775] | 0.835   | 0.246   | 0.189 |

## 4. Discussion

The study establishes a precise distinction between a mathematical property of preprocessing and the outcomes required for cross-study prediction. Within-sample normalization guarantees independence from companion patients for a frozen rule. Primary bidirectional results show that neither scale fusion nor that invariance establishes a consistent discrimination winner. Source thresholds calibrated at a nominal operating point can still produce substantially different target specificity. These endpoints answer different methodological questions and must be evaluated separately.

Global and local references redistribute the spatial contrasts available to a classifier. The local transform attenuates broad arm-level fragmentation effects along with slowly varying backgrounds. The retained block associations also show that a smaller median source-variance fraction need not remove the genomic ordering of residual variation. Hybrid exposes both scales to a regularized classifier; its paired target-specific comparisons quantify the dependence on evaluation context. Matched copy-number measurements, such as those enabled by ichorCNA [29], would directly test the genomic interpretation of the observed arm-level fragmentation contrasts.

Source calibration introduces a second transport problem. The threshold bootstrap measures finite-control uncertainty within a fixed calibration design, while source-study-held-out calibration changes the domain on which validation scores are produced. The resulting threshold and sensitivity changes demonstrate why calibration must be evaluated at the study level as well as the patient level. A higher transported threshold trades sensitivity for specificity; its utility cannot be selected using the target diagnosis.

Repeated anchor batches make the evaluated target-map dependence explicit. Holding case fraction fixed still permits substantial probability variability and changes in threshold decisions. Comparing mean anchor probabilities across fractions then separates a composition contrast from identity resampling. These are diagnostics of the EB-TargetMap rule under the tested real-patient mixtures; covariate-aware and reference-batch harmonization designs require their own evaluation. The retained within-study standardized disease effects and pair concordance distinguish target mapping and score comparability from a universal claim of biological signal removal.

### 4.1. Scope and limitations

The evidence comprises two primary sources with substantial control representation and one external source with 36 cancer participants and one healthy participant. The archived local FinaleDB/UNITE census yielded no additional human plasma WGS source with cancer and at least ten explicitly healthy participants. Cristiano lacked HCC labels. The post-hoc exploratory sensitivity analyses involving the shared-depth interval and the GI-versus-HCC restriction reuse the same public cohorts and change target sample size and composition. They explore potential confounding directions rather than providing independent external validation or confirmatory evidence. Stage was unrecorded for Jiang and other covariates remain incomplete. Study effects therefore combine technical and population differences.

Method development used the same public cohorts that motivated the representation comparison. Fold isolation addresses algorithmic leakage; independent confirmatory studies are needed to evaluate development-history and future-study uncertainty. Patient and threshold bootstraps condition on fixed fitted scores. Single-source calibration folds cannot identify multi-source harmonization or ICC filtering. The post-hoc exploratory shared-depth restriction is not exact distribution matching, and repeated anchors are dependent observations rather than additional independent patients.

Genomic block inference assumes independent chromosome contributions and sign symmetry for its null randomization. The unadjusted centroid-dispersion diagnostic and adjusted PERMANOVA address different hypotheses; detected heterogeneity precludes a pure-location explanation. Arm-level fragmentation and reference genomic annotations remain descriptive without complete matched CNA or tumor accessibility assays.

## 5. Conclusions

Within-sample normalization guarantees companion-independent target preprocessing for a frozen scoring rule. Its genomic reference scale produces target-dependent discrimination trade-offs, and the resulting classifier still requires explicit assessment of source-threshold transport. Bidirectional evaluation, source-study calibration and repeated real-patient anchors make these distinct properties testable in public cfDNA cohorts.

## Declarations

Author-identifying availability, authorship and compliance statements are supplied on the separate title page.

## Declaration of generative AI and AI-assisted technologies in the writing process


During the preparation of this work, the author used OpenAI Codex (OpenAI; desktop application used in October 2026; exact application build and underlying model version not recorded) to assist with code implementation and refactoring, unit tests, statistical-workflow and figure-generation code, reference verification, manuscript organization and drafting, and language revision. The author takes responsibility for reviewing the code logic, statistical outputs and text, for the decision to submit, and for the accuracy, integrity and interpretation of the work. The AI tool is not an author. Figures were generated by reproducible code from the recorded data.

## Figure legends

**Figure 1. Study-associated geometry and multivariate dispersion.** (a) Participant composition and observed healthy fragment-length profiles. The Snyder healthy profile is one descriptive observation. (b) Descriptive UMAP by study and disease; unique Type III R² is the drop-one, covariate-adjusted increment divided by total SS. (c) Unique study and disease window variance fractions. (d) Bias-adjusted Euclidean distances to each study centroid; study PERMDISP accompanies the adjusted association test and does not isolate a pure location effect.

**Figure 2. Primary bilateral discrimination and threshold transport.** (a) Target-specific ElasticNet AUROCs. (b) The same primary target-by-representation estimates with a sequential blue palette. (c) Achieved target specificity minus the nominal 95% source operating point. Exact operating-point intervals and conditional source-threshold bootstrap intervals are supplied in the accompanying tables. (d) Equal-weight mean AUROC across the two primary targets. Local denotes Local_LSN; source ICC filtering is unidentifiable in this one-source design. Secondary classifiers and random-CV comparisons appear in Supplementary Figure S6 and its table.

**Figure 3. Mathematical companion independence and empirical target-map dependence.** (a) The sample-wise functional property for frozen source parameters. (b) EB-TargetMap probability shifts of fixed real anchors across 50 companion draws at each case fraction, relative to replicate 0 at fraction 0.50; boxes show median and IQR, whiskers extend to 1.5 IQR. (c) The complete signed shift distribution at fraction 0.50 isolates companion-identity variability. (d) Anchor-evaluation decision-flip rates relative to the same reference. Observations repeat ten anchors per target and are not independent patients. Individual predictions, distributions and batch memberships are supplied.

**Figure 4. Spatial context and local feature ranks.** (a,b) Raw unique study-window variance against GC and telomere distance, with whole-chromosome bootstrap intervals on Spearman rho; block-wild null assumptions are specified in Methods. (c) Healthy median chr8 raw and local trajectories; gaps across excluded windows remain gaps. (d) Top 20 fixed Local_LSN_ICC coefficient ranks across the original three-study fits, shown with a blue palette and visible cell borders. Darker values indicate higher importance (lower rank), gray indicates unavailable/filtered; C/J/S denote held-out Cristiano/Jiang/Snyder. These coefficients are the fixed local model's descriptive feature context.

**Supplementary Figure S5. Secondary score comparability and source calibration.** (a) Within-study and cross-study pair contributions to pooled Combined LOCO AUROC. (b) Mixed-source OOF versus source-study-held-out thresholds fitted before testing Snyder. Pooled discrimination is a cross-model score diagnostic; the primary performance evidence is direction-specific.

**Supplementary Figure S6. Secondary classifier comparisons.** Bilateral mean target AUROC by representation and fixed classifier. Gray cells denote unidentifiable EB-TargetMap/DANN configurations. Local uses its Local_LSN branch in the primary one-source design. Full random and study-held-out metrics are in Supplementary Table S6.

## References

1. Cristiano S, Leal A, Phallen J, et al. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature. 2019;570(7761):385-389. https://doi.org/10.1038/s41586-019-1272-6

2. Jiang P, Chan CWM, Chan KCA, et al. Lengthening and shortening of plasma DNA in hepatocellular carcinoma patients. Proceedings of the National Academy of Sciences. 2015;112(11):E1317-E1325. https://doi.org/10.1073/pnas.1500076112

3. Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. Cell. 2016;164(1-2):57-68. https://doi.org/10.1016/j.cell.2015.11.050

4. Lo YMD, Han DSC, Jiang P, Chiu RWK. Epigenetics, fragmentomics, and topology of cell-free DNA in liquid biopsies. Science. 2021;372(6538):eaaw3616. https://doi.org/10.1126/science.aaw3616

5. Mouliere F, Chandrananda D, Piskorz AM, et al. Enhanced detection of circulating tumor DNA by fragment size analysis. Science Translational Medicine. 2018;10(466):eaat4921. https://doi.org/10.1126/scitranslmed.aat4921

6. Ulz P, Thallinger GG, Auer M, et al. Inferring expressed genes by whole-genome sequencing of plasma DNA. Nature Genetics. 2016;48(10):1273-1278. https://doi.org/10.1038/ng.3648

7. Mathios D, Johansen JS, Cristiano S, et al. Detection and characterization of lung cancer using cell-free DNA fragmentomes. Nature Communications. 2021;12(1):5060. https://doi.org/10.1038/s41467-021-24994-w

8. Doebley AL, Ko M, Liao H, et al. A framework for clinical cancer subtyping from nucleosome profiling of cell-free DNA. Nature Communications. 2022;13(1):7475. https://doi.org/10.1038/s41467-022-35076-w

9. Leek JT, Scharpf RB, Bravo HC, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nature Reviews Genetics. 2010;11(10):733-739. https://doi.org/10.1038/nrg2825

10. Goh WWB, Wang W, Wong L. Why Batch Effects Matter in Omics Data, and How to Avoid Them. Trends in Biotechnology. 2017;35(6):498-507. https://doi.org/10.1016/j.tibtech.2017.02.012

11. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8(1):118-127. https://doi.org/10.1093/biostatistics/kxj037

12. Zheng H, Zhu MS, Liu Y. FinaleDB: a browser and database of cell-free DNA fragmentation patterns. Bioinformatics. 2021;37(16):2502-2503. https://doi.org/10.1093/bioinformatics/btaa999

13. Wang H, Mennea PD, McAndrew G, et al. A scalable deep-learning framework for cancer detection using cell-free DNA shallow whole-genome sequencing. Science Advances. 2026;12(28):eady9432. https://doi.org/10.1126/sciadv.ady9432

14. Zou H, Hastie T. Regularization and Variable Selection Via the Elastic Net. Journal of the Royal Statistical Society Series B: Statistical Methodology. 2005;67(2):301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x

15. Friedman J, Hastie T, Tibshirani R. Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software. 2010;33(1):1-22. https://doi.org/10.18637/jss.v033.i01

16. Breiman L. Random Forests. Machine Learning. 2001;45(1):5-32. https://doi.org/10.1023/A:1010933404324

17. Chen T, Guestrin C. XGBoost: A Scalable Tree Boosting System. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 2016. p. 785-794. https://doi.org/10.1145/2939672.2939785

18. Ganin Y, Ustinova E, Ajakan H, et al. Domain-Adversarial Training of Neural Networks. In: Csurka G, editor. Domain Adaptation in Computer Vision Applications. Cham: Springer International Publishing; 2017. p. 189-209. https://doi.org/10.1007/978-3-319-58347-1_10

19. McInnes L, Healy J, Saul N, Großberger L. UMAP: Uniform Manifold Approximation and Projection. Journal of Open Source Software. 2018;3(29):861. https://doi.org/10.21105/joss.00861

20. Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature. 2020;585(7825):357-362. https://doi.org/10.1038/s41586-020-2649-2

21. Saito T, Rehmsmeier M. The Precision-Recall Plot Is More Informative than the ROC Plot When Evaluating Binary Classifiers on Imbalanced Datasets. PLOS ONE. 2015;10(3):e0118432. https://doi.org/10.1371/journal.pone.0118432

22. Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW; on behalf of the STRATOS Topic Group on diagnostic tests and prediction models. Calibration: the Achilles heel of predictive analytics. BMC Medicine. 2019;17(1):230. https://doi.org/10.1186/s12916-019-1466-7

23. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the Areas under Two or More Correlated Receiver Operating Characteristic Curves: A Nonparametric Approach. Biometrics. 1988;44(3):837-845. https://doi.org/10.2307/2531595

24. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society: Series B (Methodological). 1995;57(1):289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

25. Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. Patterns. 2023;4(9):100804. https://doi.org/10.1016/j.patter.2023.100804

26. Anderson MJ. A new method for non‐parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32-46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x

27. Freedman D, Lane D. A Nonstochastic Interpretation of Reported Significance Levels. Journal of Business & Economic Statistics. 1983;1(4):292-298. https://doi.org/10.1080/07350015.1983.10509354

28. Anderson MJ. Distance-Based Tests for Homogeneity of Multivariate Dispersions. Biometrics. 2006;62(1):245-253. https://doi.org/10.1111/j.1541-0420.2005.00440.x

29. Adalsteinsson VA, Ha G, Freeman SS, et al. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. Nature Communications. 2017;8(1):1324. https://doi.org/10.1038/s41467-017-00965-y
