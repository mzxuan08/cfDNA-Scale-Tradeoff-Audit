"""Compose the major revision directly from frozen and newly computed CSVs."""
import json,re
from pathlib import Path
import numpy as np
import pandas as pd
from .major_revision import PRIMARY,COMBINED,EXTERNAL,display_method,PRIMARY_STUDIES

def compile_revision(root):
    root=Path(root);base=root/'results/step09_major';out=root/'submission_package';tables=base/'tables'
    get=lambda n:pd.read_csv(tables/(n+'.csv'))
    old=(root/'results/archive_before_step09_major/results/MANUSCRIPT_FULL_TEXT.md').read_text(encoding='utf-8')
    author=old.split('## Abstract')[0]
    title=author.splitlines()[0][2:]
    metrics=pd.read_csv(root/'results/tables/Table3_Benchmark_Metrics_Full.csv')
    primary=metrics[metrics.cohort_scope.eq(PRIMARY)&metrics['mode'].eq('loco')&metrics.model.eq('elasticnet')&metrics.center.isin(PRIMARY_STUDIES)&metrics.status.eq('complete')].copy()
    row=lambda study,method:primary[primary.center.eq(study)&primary.pipeline.eq(method)].iloc[0]
    c,j=PRIMARY_STUDIES;cg,ch=row(c,'Global_Zscore'),row(c,'Hybrid_MultiScale');jg,jh=row(j,'Global_Zscore'),row(j,'Hybrid_MultiScale')
    delta=get('Primary_Hybrid_Global_Delta').set_index('target');dc,dj=delta.loc[c],delta.loc[j]
    means=primary.groupby('pipeline').auroc.mean()
    thresholds=get('Threshold_Calibration_Comparison');anchors=get('Anchor_Distribution_Summary');contrasts=get('Anchor_Prevalence_Contrasts');blocks=get('Chromosome_Block_Associations');disp=get('PERMDISP_Dispersion');sens=get('Depth_CancerType_Sensitivity');cohort=pd.read_csv(root/'results/step10_final/tables/Table2_Study_Depth_Fragments_Composition.csv') if (root/'results/step10_final/tables/Table2_Study_Depth_Fragments_Composition.csv').exists() else get('Table2_Study_Depth_Fragments_Composition')
    fixed=anchors[anchors.method.eq('EB-TargetMap')&anchors.case_fraction.eq(.5)].set_index('target');fc,fj=fixed.loc[c],fixed.loc[j]
    ci=lambda lo,hi:f'[{lo:.3f}, {hi:.3f}]'
    abstract=f'''**Background:** Within-sample cfDNA normalization makes target preprocessing independent of companion patients. Its implications for genomic scale, discrimination and threshold transport require separate evaluation.

**Methods:** We analysed 585 source-identified participants, prioritizing a bidirectional study-held-out benchmark of 548. Source-only calibration, paired patient bootstraps and repeated real-patient companion resampling evaluated operating points and score dependence; a separate case-dominated cohort contributed sensitivity evidence.

**Results:** In Cristiano, Global and Hybrid AUROC/AUPRC were {cg.auroc:.3f}/{cg.auprc:.3f} and {ch.auroc:.3f}/{ch.auprc:.3f}; Hybrid−Global ΔAUROC was {dc.delta_auroc:.3f} (95% CI {ci(dc.delta_low,dc.delta_high)}; BH q={dc.q:.3f}). In Jiang, corresponding AUROC/AUPRC were {jg.auroc:.3f}/{jg.auprc:.3f} and {jh.auroc:.3f}/{jh.auprc:.3f}, with ΔAUROC {dj.delta_auroc:+.3f} ({ci(dj.delta_low,dj.delta_high)}; q={dj.q:.3f}). Equal-weight AUROCs were {means['Raw_Ratio']:.3f} for Raw, {means['Global_Zscore']:.3f} for Global and {means['Hybrid_MultiScale']:.3f} for Hybrid. Source thresholds calibrated to 95% specificity yielded target specificities {ch.specificity95:.3f}–{cg.specificity95:.3f} in Cristiano and {jg.specificity95:.3f}–{jh.specificity95:.3f} in Jiang. At fixed target case fraction 0.50, EB-TargetMap companion resampling produced decision-flip rates {fc.threshold_crossing_flip_rate:.3f} and {fj.threshold_crossing_flip_rate:.3f}.

**Conclusions:** Companion independence is a mathematical property of within-sample transformations. Normalization-scale discrimination is target-dependent, and companion independence alone does not ensure source-threshold transport.

**Keywords:** cell-free DNA; fragmentomics; normalization scale; study-held-out validation; threshold transport; companion dependence
'''
    intro=old.split('## 1. Introduction\n\n')[1].split('\n\nNormalization scale offers')[0]
    intro=intro.replace('Griffin further links nucleosome protection profiles with cancer phenotypes while explicitly addressing GC bias [8].', 'Doebley et al. further link nucleosome protection profiles with cancer phenotypes while explicitly addressing GC bias [8].')
    body=author+'## Abstract\n\n'+abstract+'\n## 1. Introduction\n\n'+intro+r'''

Global and local normalization define different genomic references. A genome-wide reference preserves broad relative deviations, while a same-arm local contrast emphasizes departures from nearby windows. Combining the scales makes their complementary information testable. A useful evaluation must distinguish the disease contrasts preserved by a representation from the study on which a classifier is tested and the threshold at which it is operated.

Cross-sample alignment adds dependence on the target patient mixture. Empirical-Bayes source harmonization estimates study-level location and scale [11]; an additional unseen-target median/IQR map uses the distribution of the patients currently available. Under $F_c=(1-\pi_c)F_{0c}+\pi_cF_{1c}$, both case fraction and individual companion identities can change these summaries. Within-sample transformations provide an exact comparison because the individual supplies their own target reference.

We evaluate normalization scale through bidirectional discrimination, source-locked threshold transfer and companion dependence. A paired analysis identifies target-dependent differences between Global and Hybrid. Calibration experiments distinguish uncertainty in the source threshold from the specificity actually achieved in the target. Repeated fixed-anchor batches separate changes in case proportion from companion-identity variability. Genomic block inference and multivariate dispersion diagnostics place study-associated variation in its spatial and distributional context.

## 2. Materials and methods

### 2.1. Cohort and source provenance

Published UNITE supplementary metadata and FinaleDB fragment intervals supplied the clinical and technical inputs [12,13]. Plasma WGS specimens required an explicit cancer or Healthy source label, original sequencing depth 0.1–5.0× and the recorded extraction QC. From 628 extracted specimens, one prespecified specimen per known participant was selected by earliest time point, shallowest depth and run identifier. The unresolved-identity specimen was excluded, leaving 585 source-identified participants. Undocumented cross-study aliases cannot be reconciled from the source inventory.

Cristiano and Jiang constituted the Primary Bilateral Benchmark (548 participants). Each direction trained on one entire study and tested the other. Snyder was excluded from primary training, selection and calibration, providing a separate external sensitivity cohort (37 participants, 36 cancer/1 healthy). Secondary leave-one-center-out cross-validation (LOCO-CV), with study source as the grouping variable, used all three studies. Study Source is a proxy combining experimental handling, recruitment and population composition.

'''
    counts=pd.read_csv(root/'results/tables/Analysis_Cohort_Counts.csv');counts=counts.rename(columns={'center':'Study','n':'Participants','cancer':'Cancer','healthy':'Healthy'})
    body+='**Table 1. Independent participant cohort.**\n\n'+counts.to_markdown(index=False,disable_numparse=True)+'\n\n'
    body+=r'''### 2.2. Features and source-fitted preprocessing

The hg19/GRCh37 reference contains 477 complete autosomal 5 Mb windows after fixed centromere, terminal-risk, cytoband-risk and blacklist exclusions. Autosomal fragments with MAPQ≥30 and length 80–400 bp contribute to the global integer-length distribution; midpoint assignment supplies window counts. Short and long counts cover 100–150 and 151–220 bp:

$$R_{si}=\frac{S_{si}}{L_{si}+1}.$$

Global Z-score references the individual's window mean and standard deviation. Local self-normalization (LSN) references other retained windows on the same chromosome arm, up to three on each side:

$$G_{si}=\frac{R_{si}-\overline{R}_{s}}{\operatorname{sd}_{i}(R_{si})}.$$

$$\mathrm{LSN}_{si}=\log_2\frac{R_{si}+\epsilon_s}{\operatorname{median}_{j\in N_3(i),j\ne i}R_{sj}+\epsilon_s}.$$

The relative offset is $\epsilon_s=10^{-6}\operatorname{median}_iR_{si}$.

Hybrid concatenates $[G_s,\mathrm{LSN}_s]$ before source-fitted column scaling. Raw, Global and Hybrid use StandardScaler; Local uses RobustScaler. All imputers and scalers fit source-training patients exclusively, including inside each calibration fold. Per-column source medians replace post-transform NaNs; an entirely missing source column has zero as its documented fallback. The ratio and length matrices each had 0 missing cells; the ratio matrix contained 279045 cells. Feature matrices use float32. In multiple-source fits, a one-way healthy-center intraclass correlation (ICC) clustering score ranks windows, and Local_LSN_ICC removes the highest quarter. In one-source bilateral fits, this score is unidentifiable and the representation is precisely Local_LSN, without an ICC estimate.

Elastic-net logistic regression used $C=1.0$, `l1_ratio=0.5`, solver saga, tolerance 0.001 and maximum 20000 iterations, with seed 2026 [18,19]. The fixed Hybrid radius was specified before outer evaluation. Source-only nested local selection evaluated radii 1, 2, 3, 5 and 8 and exclusion fractions 0, 0.15, 0.25, 0.35 and 0.50 as a supplementary analysis. A single healthy source restricts exclusion to zero. Three-fold source-only mean eligible-study AUROC selected the settings; ties favored lower exclusion, radius nearer three, then smaller radius. Random Forest, CUDA XGBoost, LightGBM, MLP and domain-adversarial neural network comparisons are supplied in the supplement [21,22,23]. NumPy and UMAP supported numerical computation and descriptive embeddings [14,15]. Exact package versions are listed in Supplementary Table S2.

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

Primary outcomes were target-specific AUROC, area under the precision–recall curve (AUPRC), equal-weight mean target AUROC and achieved specificity/sensitivity at the locked source threshold. Secondary measures included Brier score and 10-bin expected calibration error (ECE) [26,27]. Patient bootstrap intervals and paired ΔAUROC intervals used 1000 diagnosis-stratified resamples conditional on fixed scoring rules. Exact binomial intervals accompany operating-point counts. Snyder specificity is reported as a fraction with denominator one, solely descriptively.

DeLong tests compare fixed individual scoring rules within each independent held-out primary study [24]. The retained primary Benjamini–Hochberg (BH) family has 38 finite tests: Hybrid versus Raw, Global and Local for each of five estimable classifiers in both targets (30), plus Tuned Local versus Raw, Global, Local and Hybrid for ElasticNet in both targets (8) [25]. Supplementary Table S1 enumerates every hypothesis. Unavailable EB-TargetMap and DANN cells contribute no p value. Secondary families are separate; no pooled cross-model DeLong test is asserted. Source isolation addresses algorithmic leakage [20].

### 2.5. Companion resampling and cohort sensitivity analyses

For each adequately represented target, the other two studies trained a frozen ElasticNet scoring rule. Five healthy and five cancer anchors were held fixed. At target case fractions 0.25, 0.50 and 0.75, each diagnostic batch contained 32 distinct real patients; 50 independent companion priority draws per fraction generated the batch series. Matched priorities coupled the fractions within a replicate. Outcomes defined diagnostic compositions after evaluation and selected neither predictive settings nor benchmark thresholds. At fixed fraction 0.50, identity resampling measures sampling variability independently of changing the fraction. Probability shifts reference each anchor's prediction in replicate 0 at fraction 0.50. We report absolute-shift median, quartiles, IQR, full patient-level distributions, within-fraction ranges and decision flips relative to this same reference. Repeated anchor observations are dependent and do not provide independent patient replication.

Post-hoc exploratory sensitivity analyses restricted both primary studies to a fixed 2–4× original-depth interval, then repeated source-only fitting and calibration. This is a common-range restriction, not exact distribution matching. The post-hoc exploratory GI-versus-HCC analysis retained all controls and restricted Cristiano cancer to Bile Duct, Colorectal, Duodenal, Gastric or Pancreatic labels; Jiang remained HCC. No HCC-labelled Cristiano participant was available, so this gastrointestinal restriction evaluates a narrower disease spectrum without claiming a matched HCC comparison. These post-hoc exploratory sensitivity analyses reuse the same public cohorts and change target sample size and composition. They explore the potential direction of confounding and do not constitute independent external validation or confirmatory evidence.

### 2.6. Spatial and multivariate diagnostics

Coverage-adjusted PERMANOVA used Euclidean geometry and 999 restricted Freedman–Lane residual permutations [16,17]. Throughout, unique Type III $R^2=(SS_{reduced}-SS_{full})/SS_{total}$ denotes the drop-one contribution conditional on the other covariates, not partial $R^2$. Window eta-squared uses the analogous unique increment. For the descriptive dispersion diagnostic, columns were standardized across all audit participants, without contributing to prediction fitting. Group-centroid distances, bias-adjusted by $\sqrt{n_g/(n_g-1)}$, supplied an additional PERMDISP test with 999 unrestricted permutations of intercept-only residuals of the fixed centroid distances [29]. This dispersion diagnostic is unadjusted and does not decompose adjusted PERMANOVA into causal location and dispersion components.

Spearman associations with GC and telomere/centromere distance used a 1000-replicate whole-chromosome bootstrap over 22 autosomes. For a block-wild null, centered-rank cross-products were summed within chromosomes and their signs independently randomized in 999 replicates. This retains within-chromosome dependence under independent, sign-symmetric chromosome contributions. BH correction covered the nine representation-by-annotation tests; the four PERMDISP tests formed a separate BH family. GM12878 compartments, repeats and segmental duplications provide reference context; arm-level fragmentation effects are not independently measured CNA or tumor accessibility.

## 3. Results

### 3.1. Primary discrimination depends on target study and reference scale

'''
    perma=pd.read_csv(root/'results/tables/Audit_PERMANOVA.csv');perma=perma[perma.family.eq('raw_ratio')].set_index('term')
    body+=f"Across the recorded raw-ratio geometry, unique Type III R² was {perma.loc['center','R2']:.3f} for study source (p={perma.loc['center','p']:.3f}) and {perma.loc['disease','R2']:.3f} for disease status (p={perma.loc['disease','p']:.3f}). These fractions describe study-associated multivariate variation and should not be interpreted as pure centroid shifts because dispersion heterogeneity was also detected. Figure 1 presents these descriptive associations alongside dispersion.\n\n"
    body+=f'''Global and Hybrid changed rank across targets (Figure 2). In Cristiano, Hybrid−Global ΔAUROC was {dc.delta_auroc:.3f}, 95% CI {ci(dc.delta_low,dc.delta_high)}, DeLong p={dc.p:.6f}, BH q={dc.q:.6f}. In Jiang, ΔAUROC was {dj.delta_auroc:+.3f}, 95% CI {ci(dj.delta_low,dj.delta_high)}, p={dj.p:.6f}, q={dj.q:.6f}. Equal-weight AUROC was {means['Raw_Ratio']:.3f} for Raw, {means['Global_Zscore']:.3f} for Global and {means['Hybrid_MultiScale']:.3f} for Hybrid. These target-specific results identify a scale trade-off without a consistent discrimination winner.\n\n'''
    t=[]
    for r in primary[primary.pipeline.isin(['Raw_Ratio','Global_Zscore','Proposed_LSN_ICC','Hybrid_MultiScale'])].itertuples():
        d=delta.loc[r.center] if r.pipeline=='Hybrid_MultiScale' else None
        t.append({'Study':r.center.split(' et al')[0],'Method':display_method(r.pipeline,PRIMARY),'AUROC':f'{r.auroc:.3f}','AUROC 95% CI':ci(r.auroc_low,r.auroc_high),'AUPRC':f'{r.auprc:.3f}','Specificity shift (pp)':f'{r.specificity_deviation_pp:.3f}','Hybrid−Global Δ [95% CI]':f'{d.delta_auroc:+.3f} {ci(d.delta_low,d.delta_high)}' if d is not None else '—','BH q':f'{d.q:.3f}' if d is not None else '—'})
    body+='**Table 3A. Primary bilateral fixed ElasticNet results.**\n\n'+pd.DataFrame(t).to_markdown(index=False,disable_numparse=True)+'\n\nPaired difference intervals condition on the fixed target scoring rules; BH q values belong to the 38-test primary family.\n\n'
    # Table 2 is positioned in the cohort sensitivity subsection but numbered logically.
    t2=[]
    for r in cohort.itertuples():
        t2.append({'Study':r.study.split(' et al')[0],'n (cancer/healthy)':f'{r.n} ({r.cancer}/{r.healthy})','Depth median [Q1,Q3] ×':f'{r.sequencing_depth_median:.3f} [{r.sequencing_depth_q1:.3f},{r.sequencing_depth_q3:.3f}]','Fragments median [Q1,Q3], M':f'{r.accepted_fragments_median/1e6:.3f} [{r.accepted_fragments_q1/1e6:.3f},{r.accepted_fragments_q3/1e6:.3f}]','Cancer stage':r.stage_distribution,'Cancer types (n)':r.cancer_types})
    t2=pd.DataFrame(t2).set_index('Study').T.reset_index().rename(columns={'index':'Characteristic'})
    table2='**Table 2. Study-specific sequencing depth, fragment counts, and cancer-type composition.**\n\n'+t2.to_markdown(index=False,disable_numparse=True)+'\n\nFragment counts are accepted autosomal 80–400 bp intervals after QC. Stage categories are preserved as supplied, including 0 and X; Jiang stage is unrecorded.\n\n'
    body=body.replace('### 2.2. Features and source-fitted preprocessing',table2+'### 2.2. Features and source-fitted preprocessing')
    body+='### 3.2. Companion independence leaves threshold transport unresolved\n\n'
    body+=f'''At the source-locked threshold, Global/Hybrid specificity was {cg.specificity95:.3f}/{ch.specificity95:.3f} in Cristiano and {jg.specificity95:.3f}/{jh.specificity95:.3f} in Jiang. These achieved values differ substantially from the 0.950 source target despite patient-wise preprocessing independence.\n\n'''
    for study in [c,j]:
        for method in ['Global_Zscore','Hybrid_MultiScale']:
            r=thresholds[thresholds.scope.eq(PRIMARY)&thresholds.target.eq(study)&thresholds.method.eq(method)&thresholds.scheme.eq('PooledSource3Fold')].iloc[0]
            body+=f'''For {method.replace('_',' ')} targeting {study.split(' et al')[0]}, the source threshold was {r.threshold:.3f} (conditional bootstrap 95% CI {ci(r.threshold_low,r.threshold_high)}). Propagating source-control score uncertainty gave target specificity {ci(r.specificity_from_threshold_low,r.specificity_from_threshold_high)} with target patients fixed.\n\n'''
    ext=thresholds[thresholds.scope.eq(EXTERNAL)&thresholds.status.eq('complete')];t=[]
    for r in ext.itertuples():
        t.append({'Method':r.method,'Source calibration':'Mixed OOF' if r.scheme=='PooledSource3Fold' else 'Source-study held out','Threshold [95% CI]':f'{r.threshold:.3f} {ci(r.threshold_low,r.threshold_high)}','Cancer sensitivity':f'{r.target_sensitivity:.3f}','Healthy result (descriptive only)':f'{int(r.target_healthy_success)}/{int(r.target_healthy_n)}; exact 95% CI {ci(r.target_specificity_low_exact,r.target_specificity_high_exact)}'})
    body+='**Table 3B. External Snyder sensitivity and descriptive single-control outcomes.**\n\n'+pd.DataFrame(t).to_markdown(index=False,disable_numparse=True)+'\n\nSource-study-held-out EB-TargetMap calibration is unidentifiable because each calibration-training fold contains one study. A denominator of one supplies no stable population specificity estimate.\n\n'
    for method in ['Global_Zscore','Hybrid_MultiScale']:
        mixed=ext[ext.method.eq(method)&ext.scheme.eq('PooledSource3Fold')].iloc[0];lso=ext[ext.method.eq(method)&ext.scheme.eq('LeaveOneSourceStudyOut')].iloc[0]
        body+=f'''For external {method.replace('_',' ')}, source-study-held-out calibration changed the threshold from {mixed.threshold:.3f} to {lso.threshold:.3f} and Snyder cancer sensitivity from {mixed.target_sensitivity:.3f} to {lso.target_sensitivity:.3f}. This demonstrates calibration-design dependence; single-control outcomes remain descriptive.\n\n'''
    body+='### 3.3. Repeated anchors distinguish case fraction from companion identity\n\n'
    body+=f'''At fixed fraction 0.50, EB-TargetMap median absolute anchor shifts were {fc.median_abs_probability_shift:.3f} (Q1–Q3 {fc.q1_abs_probability_shift:.3f}–{fc.q3_abs_probability_shift:.3f}; IQR {fc.iqr_abs_probability_shift:.3f}) in Cristiano and {fj.median_abs_probability_shift:.3f} ({fj.q1_abs_probability_shift:.3f}–{fj.q3_abs_probability_shift:.3f}; IQR {fj.iqr_abs_probability_shift:.3f}) in Jiang. Decision-flip rates relative to the reference batch were {fc.threshold_crossing_flip_rate:.3f} and {fj.threshold_crossing_flip_rate:.3f}; the fractions of anchors ever changing decisions were {fc.fraction_anchors_ever_flip:.3f} and {fj.fraction_anchors_ever_flip:.3f}. Thus companion identity alone changes adapted scores at a fixed case fraction (Figure 3).\n\n'''
    for r in contrasts[contrasts.method.eq('EB-TargetMap')].itertuples():
        body+=f'''For {r.target.split(' et al')[0]}, the median absolute difference between mean anchor probabilities at fractions 0.75 and 0.25 was {r.median_abs_mean_075_minus_025:.3f}; the maximum per-anchor range across all tested batches was {r.max_probability_range:.3f}. The distribution table separates this mean case-fraction contrast from within-fraction identity ranges.\n\n'''
    body+='Raw, Global, Local and Hybrid preserved anchor probabilities under companion replacement to numerical tolerance, as required by the mathematical property and checked in unit tests. This invariant is a functional guarantee for frozen rules, not an empirical discrimination advantage. Full per-anchor predictions and all real batch memberships are supplied.\n\n'
    body+='### 3.4. Post-hoc exploratory sensitivity analyses\n\nThese post-hoc exploratory sensitivity analyses use the same public cohorts and alter target sample size and composition. They explore potential confounding directions and provide neither independent external validation nor confirmatory evidence.\n\n'
    for scenario,label in [('Depth_2to4x','In the post-hoc exploratory sensitivity analyses within the 2–4× common-depth interval'),('GI_vs_HCC','In the post-hoc exploratory sensitivity analyses under the gastrointestinal-versus-HCC restriction')]:
        g=sens[sens.scenario.eq(scenario)]
        for study in [c,j]:
            h=g[g.target.eq(study)].set_index('method');a,b=h.loc['Global_Zscore'],h.loc['Hybrid_MultiScale']
            body+=f'''{label}, {study.split(' et al')[0]} had {a.target_n} target participants. Global/Hybrid AUROCs were {a.auroc:.3f}/{b.auroc:.3f}, with achieved specificities {a.specificity:.3f}/{b.specificity:.3f}.\n\n'''
    body+='These post-hoc exploratory sensitivity analyses retained different source/target sizes and do not equalize stage or cancer composition. The gastrointestinal analysis contains bile duct and other digestive tumors rather than matched HCC cases. Both analyses retain the distinction between ranking and source-threshold transport (Supplementary Table S4).\n\n'
    body+='### 3.5. Spatial inference and dispersion qualify study-associated variation\n\n'
    for method in ['Raw_Ratio','Global_Zscore','LSN']:
        r=blocks[blocks.representation.eq(method)&blocks.annotation.eq('gc_fraction')].iloc[0]
        body+=f'''For {method.replace('_',' ')}, the GC–unique-study-variance Spearman correlation was {r.rho:.3f}, chromosome-bootstrap 95% CI {ci(r.block_ci_low,r.block_ci_high)}, block-wild p={r.block_permutation_p:.3f}, BH q={r.q_block:.3f}.\n\n'''
    rawdisp=disp[disp.family.eq('raw_ratio')].set_index('term')
    body+=f'''Raw-ratio PERMDISP detected study dispersion heterogeneity (F={rawdisp.loc['study','F']:.3f}, p={rawdisp.loc['study','p']:.3f}) and disease dispersion heterogeneity (F={rawdisp.loc['disease','F']:.3f}, p={rawdisp.loc['disease','p']:.3f}). The adjusted unique Type III PERMANOVA study fraction remains the recorded association measure, but its significance cannot be attributed exclusively to a centroid shift. Figure 1 and Supplementary Table S5 present dispersion alongside the genomic variance audit.\n\n'''
    alignment=pd.read_csv(root/'results/tables/Final_Alignment_Scale_Tradeoff.csv')
    effects=pd.read_csv(root/'results/tables/Disease_Effect_Sizes.csv').pivot(index=['bin_id','center'],columns='pipeline',values='cohen_d').dropna()
    arm=pd.read_csv(root/'results/tables/Final_Alignment_Scale_Tradeoff.csv').set_index('representation')
    enrichment=pd.read_csv(root/'results/tables/Filtered_Window_Enrichment.csv')
    body+=f'''The window-level audit measured median absolute disease contrasts of {effects.Raw_Ratio.abs().median():.3f} for Raw and {effects.LSN.abs().median():.3f} for local normalization. Broad arm-averaged effects were {arm.loc['Raw_Ratio','median_abs_arm_fragmentation_d']:.3f}, {arm.loc['Global_Zscore','median_abs_arm_fragmentation_d']:.3f} and {arm.loc['LSN','median_abs_arm_fragmentation_d']:.3f} for Raw, Global and Local, respectively. These fragmentation measurements describe the reference-scale trade-off; no independently measured CNA or tumor chromatin assay is inferred from them. The {len(enrichment)} reference-region enrichment tests yielded {int(enrichment.q.lt(.05).sum())} BH-significant associations.\n\n'''
    combined=metrics[metrics.cohort_scope.eq(COMBINED)&metrics['mode'].eq('loco')&metrics.model.eq('elasticnet')&metrics.center.eq('AllCenters')&metrics.status.eq('complete')]
    t3=[]
    for r in combined.itertuples():t3.append({'Method':display_method(r.pipeline),'Pooled AUROC':f'{r.auroc:.3f}','95% CI':ci(r.auroc_low,r.auroc_high),'AUPRC':f'{r.auprc:.3f}','Brier':f'{r.brier:.3f}','ECE':f'{r.ece:.3f}'})
    pairs=pd.read_csv(root/'results/tables/Step7_AUC_Pair_Decomposition.csv');pairs=pairs[pairs.pipeline.eq('InFold_ComBat')]
    within=pairs[pairs.pair_type.eq('within_study')];cross=pairs[pairs.pair_type.ne('within_study')]
    wauc=(within.correct_order_fraction*within.pair_n).sum()/within.pair_n.sum();cauc=(cross.correct_order_fraction*cross.pair_n).sum()/cross.pair_n.sum()
    body+=f'''### 3.6. Secondary pooled scores quantify cross-model comparability\n\nCancer–control pair decomposition gives within-study concordance {wauc:.3f} and cross-study concordance {cauc:.3f} for EB-TargetMap, with cross-study pairs accounting for {100*cross.pair_n.sum()/pairs.pair_n.sum():.3f}% of all pairs. Pooled three-study AUROC compares probabilities from separately fitted rules and is consequently a score-comparability diagnostic. The table retains those descriptive results separately from the primary direction-specific performance evidence (Supplementary Figure S5).\n\n'''
    body+='**Table 3C. Secondary Combined LOCO pooled-score comparability diagnostics.**\n\n'+pd.DataFrame(t3).to_markdown(index=False,disable_numparse=True)+'\n\n'
    body+=r'''## 4. Discussion

The study establishes a precise distinction between a mathematical property of preprocessing and the outcomes required for cross-study prediction. Within-sample normalization guarantees independence from companion patients for a frozen rule. Primary bidirectional results show that neither scale fusion nor that invariance establishes a consistent discrimination winner. Source thresholds calibrated at a nominal operating point can still produce substantially different target specificity. These endpoints answer different methodological questions and must be evaluated separately.

Global and local references redistribute the spatial contrasts available to a classifier. The local transform attenuates broad arm-level fragmentation effects along with slowly varying backgrounds. The retained block associations also show that a smaller median source-variance fraction need not remove the genomic ordering of residual variation. Hybrid exposes both scales to a regularized classifier; its paired target-specific comparisons quantify the dependence on evaluation context. Matched copy-number measurements, such as those enabled by ichorCNA [28], would directly test the genomic interpretation of the observed arm-level fragmentation contrasts.

Source calibration introduces a second transport problem. The threshold bootstrap measures finite-control uncertainty within a fixed calibration design, while source-study-held-out calibration changes the domain on which validation scores are produced. The resulting threshold and sensitivity changes demonstrate why calibration must be evaluated at the study level as well as the patient level. A higher transported threshold trades sensitivity for specificity; its utility cannot be selected using the target diagnosis.

Repeated anchor batches make the evaluated target-map dependence explicit. Holding case fraction fixed still permits substantial probability variability and changes in threshold decisions. Comparing mean anchor probabilities across fractions then separates a composition contrast from identity resampling. These are diagnostics of the EB-TargetMap rule under the tested real-patient mixtures; covariate-aware and reference-batch harmonization designs require their own evaluation. The retained within-study standardized disease effects and pair concordance distinguish target mapping and score comparability from a universal claim of biological signal removal.

### 4.1. Scope and limitations

The evidence comprises two primary sources with substantial control representation and one external source with 36 cancer participants and one healthy participant. The archived local FinaleDB/UNITE census yielded no additional human plasma WGS source with cancer and at least ten explicitly healthy participants. Cristiano lacked HCC labels. The post-hoc exploratory sensitivity analyses involving the shared-depth interval and the GI-versus-HCC restriction reuse the same public cohorts and change target sample size and composition. They explore potential confounding directions rather than providing independent external validation or confirmatory evidence. Stage was unrecorded for Jiang and other covariates remain incomplete. Study effects therefore combine technical and population differences.

Method development used the same public cohorts that motivated the representation comparison. Fold isolation addresses algorithmic leakage; independent confirmatory studies are needed to evaluate development-history and future-study uncertainty. Patient and threshold bootstraps condition on fixed fitted scores. Single-source calibration folds cannot identify multi-source harmonization or ICC filtering. The post-hoc exploratory shared-depth restriction is not exact distribution matching, and repeated anchors are dependent observations rather than additional independent patients.

Genomic block inference assumes independent chromosome contributions and sign symmetry for its null randomization. The unadjusted centroid-dispersion diagnostic and adjusted PERMANOVA address different hypotheses; detected heterogeneity precludes a pure-location explanation. Arm-level fragmentation and reference genomic annotations remain descriptive without complete matched CNA or tumor accessibility assays.

## 5. Conclusions

Within-sample normalization guarantees companion-independent target preprocessing for a frozen scoring rule. Its genomic reference scale produces target-dependent discrimination trade-offs, and the resulting classifier still requires explicit assessment of source-threshold transport. Bidirectional evaluation, source-study calibration and repeated real-patient anchors make these distinct properties testable in public cfDNA cohorts.

'''
    declarations=old.split('## Data and code availability\n')[1].split('## Figure legends\n')[0]
    body+='## Data and code availability\n'+declarations
    body=body.replace('All applicable roles were undertaken by the sole author, who approved the final manuscript.', 'All applicable roles were undertaken by the sole author.')
    body=body.replace('The author reviewed the code logic, statistical outputs and text, edited the material as needed, approved the manuscript, and takes full responsibility for its accuracy, integrity and interpretation.', 'The author takes responsibility for reviewing the code logic, statistical outputs and text, for the decision to submit, and for the accuracy, integrity and interpretation of the work.')
    source_description=declarations.strip().split('\n\n')[0]
    availability=source_description+'\n\nThe v0.2.0 software release package accompanies this manuscript and is prepared for deposit at https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit. It includes source-only calibration and anchor-resampling code, post-hoc exploratory sensitivity analyses, chromosome-block inference, PERMDISP, complete numeric CSV outputs and the current figures. No Zenodo DOI has yet been assigned to v0.2.0. The historical v0.1.0 archive remains available at https://doi.org/10.5281/zenodo.23653766; this DOI identifies the earlier archive rather than v0.2.0. Original software is MIT licensed; raw fragment streams, specimen feature matrices and fitted-model caches are excluded. Published data and third-party annotations retain their original terms.\n\n'
    release_file=root/'config/release.json'
    release=json.loads(release_file.read_text(encoding='utf-8')) if release_file.exists() else {}
    if release.get('version'):
        version=release['version']
        doi=release.get('zenodo_doi') or release.get('concept_zenodo_doi')
        scope='version-specific' if release.get('zenodo_doi') else 'collection (all versions)'
        availability=source_description+f'\n\nSoftware version v{version} and its publication files are identified by [{release["github_release_url"]}]({release["github_release_url"]}). The Zenodo {scope} DOI is [https://doi.org/{doi}](https://doi.org/{doi}). The Git tag identifies the frozen code, numeric CSV outputs and current figures for this version. The archive includes source-only calibration, real-patient anchor resampling, post-hoc exploratory sensitivity analyses, chromosome-block inference and PERMDISP. The repository is https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit. Original software is MIT licensed; raw fragment streams, specimen feature matrices and fitted-model caches are excluded. Published data and third-party annotations retain their original terms.\n\n'
    body=re.sub(r'## Data and code availability\n.*?(?=## CRediT authorship)', '## Data and code availability\n\n'+availability,body,flags=re.S)

    # One complete AI declaration in the entire manuscript.
    body=body.replace('## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process','## Declaration of generative AI and AI-assisted technologies in the writing process')
    legends='''## Figure legends

**Figure 1. Study-associated geometry and multivariate dispersion.** (a) Participant composition and observed healthy fragment-length profiles. The Snyder healthy profile is one descriptive observation. (b) Descriptive UMAP by study and disease; unique Type III R² is the drop-one, covariate-adjusted increment divided by total SS. (c) Unique study and disease window variance fractions. (d) Bias-adjusted Euclidean distances to each study centroid; study PERMDISP accompanies the adjusted association test and does not isolate a pure location effect.

**Figure 2. Primary bilateral discrimination and threshold transport.** (a) Target-specific ElasticNet AUROCs. (b) The same primary target-by-representation estimates with a sequential blue palette. (c) Achieved target specificity minus the nominal 95% source operating point. Exact operating-point intervals and conditional source-threshold bootstrap intervals are supplied in the accompanying tables. (d) Equal-weight mean AUROC across the two primary targets. Local denotes Local_LSN; source ICC filtering is unidentifiable in this one-source design. Secondary classifiers and random-CV comparisons appear in Supplementary Figure S6 and its table.

**Figure 3. Mathematical companion independence and empirical target-map dependence.** (a) The sample-wise functional property for frozen source parameters. (b) EB-TargetMap probability shifts of fixed real anchors across 50 companion draws at each case fraction, relative to replicate 0 at fraction 0.50; boxes show median and IQR, whiskers extend to 1.5 IQR. (c) The complete signed shift distribution at fraction 0.50 isolates companion-identity variability. (d) Anchor-evaluation decision-flip rates relative to the same reference. Observations repeat ten anchors per target and are not independent patients. Individual predictions, distributions and batch memberships are supplied.

**Figure 4. Spatial context and local feature ranks.** (a,b) Raw unique study-window variance against GC and telomere distance, with whole-chromosome bootstrap intervals on Spearman rho; block-wild null assumptions are specified in Methods. (c) Healthy median chr8 raw and local trajectories; gaps across excluded windows remain gaps. (d) Top 20 fixed Local_LSN_ICC coefficient ranks across the original three-study fits, shown with a blue palette and visible cell borders. Darker values indicate higher importance (lower rank), gray indicates unavailable/filtered; C/J/S denote held-out Cristiano/Jiang/Snyder. These coefficients are the fixed local model's descriptive feature context.

**Supplementary Figure S5. Secondary score comparability and source calibration.** (a) Within-study and cross-study pair contributions to pooled Combined LOCO AUROC. (b) Mixed-source OOF versus source-study-held-out thresholds fitted before testing Snyder. Pooled discrimination is a cross-model score diagnostic; the primary performance evidence is direction-specific.

**Supplementary Figure S6. Secondary classifier comparisons.** Bilateral mean target AUROC by representation and fixed classifier. Gray cells denote unidentifiable EB-TargetMap/DANN configurations. Local uses its Local_LSN branch in the primary one-source design. Full random and study-held-out metrics are in Supplementary Table S6.

'''
    body+=legends
    # Verify new dispersion reference with the archived Crossref works response.
    record=json.loads((base/'reports/PERMDISP_Crossref.json').read_text(encoding='utf-8'))['message']
    assert record['DOI'].lower()=='10.1111/j.1541-0420.2005.00440.x' and record['volume']=='62'
    refs=old.split('## References\n')[1].strip()+'\n\n29. Anderson MJ. Distance-Based Tests for Homogeneity of Multivariate Dispersions. Biometrics. 2006;62(1):245-253. https://doi.org/10.1111/j.1541-0420.2005.00440.x\n'
    # Renumber by first appearance, while retaining all 28 verified source works.
    lookup={int(m.group(1)):m.group(2) for m in re.finditer(r'(?m)^(\d+)\. (.*)$',refs)}
    pattern=re.compile(r'\[(\d+(?:,\d+)*)\]');order=[]
    for match in pattern.finditer(body):
        for n in map(int,match.group(1).split(',')):
            if n not in order:order.append(n)
    assert set(order)==set(range(1,30)),set(range(1,30))-set(order)
    mapping={n:i+1 for i,n in enumerate(order)}
    body=pattern.sub(lambda m:'['+','.join(str(mapping[int(n)]) for n in m.group(1).split(','))+']',body)
    full=body+'## References\n\n'+'\n\n'.join(f'{mapping[n]}. {lookup[n]}' for n in order)+'\n'
    assert all(ord(c)>=32 or c=='\n' for c in full), 'Unexpected control character in manuscript math'
    assert 'Proposed' not in full and 'InFold_ComBat' not in full
    assert full.count('Declaration of generative AI')==1 and '### 2.8.' not in full
    for path in [root/'results/MANUSCRIPT_FULL_TEXT.md',root/'results/reports/MANUSCRIPT_MASTER_DRAFT.md',out/'Manuscript_Main_Blind_and_Full.md']:
        path.write_text(full,encoding='utf-8')
    # Full numeric source tables accompany the publication, without relabelling history.
    metrics9=metrics.copy();metrics9['pipeline']=[display_method(m,s) for m,s in zip(metrics9.pipeline,metrics9.cohort_scope)]
    metrics9.to_csv(tables/'Table3_Benchmark_Metrics_Full_Revision.csv',index=False)
    pd.DataFrame([{'source_id':n,'reference':mapping[n],'citation':lookup[n]} for n in order]).to_csv(tables/'Revision_References.csv',index=False)
    revised_abstract=full.split('## Abstract\n')[1].split('## 1. Introduction')[0]
    (base/'reports/revision_abstract.md').write_text(revised_abstract,encoding='utf-8')
    # Supplementary tables are separate from the compact main narrative.
    supplement='# Supplementary material\n\n'+title+'\n\n'+legends.split('**Supplementary Figure S5.')[1].split('**Supplementary Figure S6.')[0].join(['**Supplementary Figure S5.',''])+'\n\n'
    supplement+=legends.split('**Supplementary Figure S6.')[1].join(['**Supplementary Figure S6.',''])+'\n\n'
    supplement+='## Statistical table guide\n\n'
    for title_,name,description in [
      ('S1','Primary_BH_Hypothesis_Family','All 38 primary hypotheses, models, target studies, DeLong p and BH q values. The complete family remains the historical prespecified family.'),
      ('S9','Threshold_Calibration_Comparison','Mixed-source OOF and source-study-held-out calibration, source threshold intervals and propagated operating points. Missing one-source EB-TargetMap cells are explicit.'),
      ('S3','Primary_Hybrid_Global_Delta','Paired Hybrid−Global ΔAUROC and 1000 diagnosis-stratified patient-bootstrap intervals.'),
      ('S4','Depth_CancerType_Sensitivity','Post-hoc exploratory sensitivity analyses of the shared-depth interval and GI-versus-HCC restriction, using the same public cohorts with altered target size/composition; not independent external validation or confirmatory evidence.'),
      ('S5','Chromosome_Block_Associations','Nine genomic block tests with full block bootstrap intervals and block-wild BH correction. PERMDISP_Dispersion.csv gives the four separate dispersion tests.'),
      ('S6','Table3_Benchmark_Metrics_Full_Revision','Complete fixed-classifier random and study-held-out metrics, with current method names and explicit unidentifiable cells.'),
      ('S7','Anchor_Distribution_Summary','Median, Q1, Q3, IQR, full range and threshold-crossing rates. Anchor_All_Predictions.csv gives all repeated observations and Anchor_Real_Batch_Members.csv supplies patient identities for reproducibility.'),
      ('S8','Table2_Study_Depth_Fragments_Composition','Exact depth/fragment ranges and cancer/stage categories for all sources. Source stage codes 0 and X remain as supplied.')]:
        supplement+=f'**Supplementary Table {title_}: {name}.csv.** {description}\n\n'
    versions=pd.read_csv(root/'results/tables/Methods_Software_Versions.csv')
    supplement+='**Supplementary Table S2. Software versions.**\n\n'+versions.to_markdown(index=False,disable_numparse=True)+'\n\n'
    compact=[]
    for r in thresholds[thresholds.scope.eq(PRIMARY)&thresholds.scheme.eq('PooledSource3Fold')].itertuples():
        compact.append({'Target':r.target.split(' et al')[0],'Method':r.method,'Source threshold [95% CI]':f'{r.threshold:.3f} {ci(r.threshold_low,r.threshold_high)}','Target specificity':f'{r.target_specificity:.3f}','Propagated specificity [95% CI]':ci(r.specificity_from_threshold_low,r.specificity_from_threshold_high)})
    supplement+='**Supplementary Table S9. Source threshold calibration (primary).**\n\n'+pd.DataFrame(compact).to_markdown(index=False,disable_numparse=True)+'\n\n'
    (out/'05_Supplementary_Material.md').write_text(supplement,encoding='utf-8')
    # A reviewer-facing map is generated from the actual results above.
    response='# Point-by-point major revision summary\n\n本次为已提交bioRxiv Version 1之后的本地方法学修订；尚未上传修订版，也未向期刊投稿。所有新增分析基于现有585名独立患者。原始v1矩阵、预测、表格和图件保留，新输出位于results/step09_major。\n\n'
    for heading,content in [
        ('1 主指标与命名',f'Abstract、Results 3.1和Table 3A改用primary双向结果。Cristiano ΔAUROC={dc.delta_auroc:.3f} {ci(dc.delta_low,dc.delta_high)}，q={dc.q:.3f}；Jiang ΔAUROC={dj.delta_auroc:+.3f} {ci(dj.delta_low,dj.delta_high)}，q={dj.q:.3f}。主稿移除Proposed；Local_LSN与Local_LSN_ICC按设计区分。Combined pooled AUROC仅作分数可比性诊断。'),
        ('2 阈值不确定性与源研究留一校准','Methods 2.4、Results 3.2、Table 3B及Threshold_Calibration_Comparison.csv给出1000次source-control bootstrap。LOSO校准为新拟合结果，EB-TargetMap因仅剩一个源研究无法估计；单健康对照改为0/1或1/1并附exact CI。'),
        ('3 固定锚点扩展',f'Methods 2.3/2.5、Results 3.3及Figure 3。零同伴依赖写为数学性质及单元测试。每目标10锚点、三种病例比例各50次同伴抽样。固定比例0.50时EB-TargetMap中位绝对变化Cristiano={fc.median_abs_probability_shift:.3f}、Jiang={fj.median_abs_probability_shift:.3f}，翻转率分别{fc.threshold_crossing_flip_rate:.3f}/{fj.threshold_crossing_flip_rate:.3f}。全分布、batch成员、IQR与均值比例对比均有CSV。'),
        ('4 队列异质性与敏感性','Table 2填补编号缺口，列深度、片段数、癌种和真实分期。Jiang分期未记录。Cristiano没有HCC标签，故执行明确的消化道肿瘤限制而不冒称同癌种。2–4×限制为共同深度范围，不是精确匹配；全部重拟合结果见Depth_CancerType_Sensitivity.csv。'),
        ('5 空间与多元统计','Methods 2.6及Results 3.5。22染色体整体bootstrap、block-wild rank-product sign randomization及其假设明确披露。9项BH家族；PERMDISP四项检验均p=0.001，因此不将PERMANOVA解释为纯中心位置效应。Primary_BH_Hypothesis_Family.csv逐项枚举38项主家族。'),
        ('6 方法细节与结构','Methods补齐C=1、l1_ratio=0.5、填充规则和所有scaler源样本范围；LSN、ICC、LOCO、OOF、ECE首次出现补全全称，次要分类器和Random-CV移入补充；引言删除具体结果；AI披露仅文末一处。'),
        ('7 图件与排版','所有主图及S5重新导出300 DPI PNG和矢量PDF；蓝色热图及清晰边框，最小绘图字体8.8 pt，按实际嵌入宽度核验；主表完整、表头重复、禁止拆行；最终Word/PDF逐页检查另记验证文件。')]:
        response+='## '+heading+'\n\n'+content+'\n\n'
    (out/'Point_by_Point_Revision_Summary.md').write_text(response,encoding='utf-8')
    return full
