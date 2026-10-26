# Supplementary material

Scale Trade-offs in Within-Sample Self-Normalization for Cross-Study Generalization in Public cfDNA Cohorts

**Supplementary Figure S5. Secondary score comparability and source calibration.** (a) Within-study and cross-study pair contributions to pooled Combined LOCO AUROC. (b) Mixed-source OOF versus source-study-held-out thresholds fitted before testing Snyder. Pooled discrimination is a cross-model score diagnostic; the primary performance evidence is direction-specific.



**Supplementary Figure S6. Secondary classifier comparisons.** Bilateral mean target AUROC by representation and fixed classifier. Gray cells denote unidentifiable EB-TargetMap/DANN configurations. Local uses its Local_LSN branch in the primary one-source design. Full random and study-held-out metrics are in Supplementary Table S6.



## Statistical table guide

**Supplementary Table S1: Primary_BH_Hypothesis_Family.csv.** All 38 primary hypotheses, models, target studies, DeLong p and BH q values. The complete family remains the historical prespecified family.

**Supplementary Table S9: Threshold_Calibration_Comparison.csv.** Mixed-source OOF and source-study-held-out calibration, source threshold intervals and propagated operating points. Missing one-source EB-TargetMap cells are explicit.

**Supplementary Table S3: Primary_Hybrid_Global_Delta.csv.** Paired Hybrid−Global ΔAUROC and 1000 diagnosis-stratified patient-bootstrap intervals.

**Supplementary Table S4: Depth_CancerType_Sensitivity.csv.** Post-hoc exploratory sensitivity analyses of the shared-depth interval and GI-versus-HCC restriction, using the same public cohorts with altered target size/composition; not independent external validation or confirmatory evidence.

**Supplementary Table S5: Chromosome_Block_Associations.csv.** Nine genomic block tests with full block bootstrap intervals and block-wild BH correction. PERMDISP_Dispersion.csv gives the four separate dispersion tests.

**Supplementary Table S6: Table3_Benchmark_Metrics_Full_Revision.csv.** Complete fixed-classifier random and study-held-out metrics, with current method names and explicit unidentifiable cells.

**Supplementary Table S7: Anchor_Distribution_Summary.csv.** Median, Q1, Q3, IQR, full range and threshold-crossing rates. Anchor_All_Predictions.csv gives all repeated observations and Anchor_Real_Batch_Members.csv supplies patient identities for reproducibility.

**Supplementary Table S8: Table2_Study_Depth_Fragments_Composition.csv.** Exact depth/fragment ranges and cancer/stage categories for all sources. Source stage codes 0 and X remain as supplied.

**Supplementary Table S2. Software versions.**

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

**Supplementary Table S9. Source threshold calibration (primary).**

| Target    | Method            | Source threshold [95% CI]   | Target specificity   | Propagated specificity [95% CI]   |
|:----------|:------------------|:----------------------------|:---------------------|:----------------------------------|
| Cristiano | Global_Zscore     | 0.629 [0.407, 0.744]        | 0.660                | [0.512, 0.737]                    |
| Cristiano | Hybrid_MultiScale | 0.559 [0.334, 0.617]        | 0.632                | [0.517, 0.675]                    |
| Cristiano | Local_LSN         | 0.671 [0.517, 0.805]        | 0.512                | [0.445, 0.632]                    |
| Cristiano | Raw_Ratio         | 0.695 [0.679, 0.709]        | 0.450                | [0.105, 0.770]                    |
| Jiang     | Global_Zscore     | 0.729 [0.576, 0.881]        | 0.406                | [0.281, 0.594]                    |
| Jiang     | Hybrid_MultiScale | 0.764 [0.624, 0.883]        | 0.531                | [0.406, 0.719]                    |
| Jiang     | Local_LSN         | 0.866 [0.820, 0.915]        | 0.812                | [0.750, 0.812]                    |
| Jiang     | Raw_Ratio         | 0.537 [0.506, 0.570]        | 1.000                | [1.000, 1.000]                    |

