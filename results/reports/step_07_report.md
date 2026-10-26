# Step 07：多尺度验证与科学叙事修订

本次只使用原有585名独立患者。主双向队列548人；外部高负荷队列37人。

## 主要结果

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

## 来源锁定的外部敏感性

| pipeline          |   sensitivity95 |   sensitivity_exact_low |   sensitivity_exact_high |   specificity95 |   source_threshold_min |
|:------------------|----------------:|------------------------:|-------------------------:|----------------:|-----------------------:|
| Global_Zscore     |           0.944 |                   0.813 |                    0.993 |           0.000 |                  0.771 |
| Hybrid_MultiScale |           0.944 |                   0.813 |                    0.993 |           0.000 |                  0.777 |
| InFold_ComBat     |           0.028 |                   0.001 |                    0.145 |           1.000 |                  0.658 |
| Proposed_LSN_ICC  |           0.917 |                   0.775 |                    0.982 |           0.000 |                  0.914 |
| Raw_Ratio         |           0.972 |                   0.855 |                    0.999 |           0.000 |                  0.567 |
| Tuned_LSN_ICC     |           0.972 |                   0.855 |                    0.999 |           0.000 |                  0.885 |

## 机制诊断

| pipeline          | target_study          |   anchors |   batch_n |   max_anchor_probability_range |   mean_anchor_probability_range | source_training_studies                  | scope                                                                                       |
|:------------------|:----------------------|----------:|----------:|-------------------------------:|--------------------------------:|:-----------------------------------------|:--------------------------------------------------------------------------------------------|
| Raw_Ratio         | Cristiano et al, 2019 |        10 |        32 |                       0.000000 |                        0.000000 | Jiang et al, 2015;Snyder et al, 2016     | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| Global_Zscore     | Cristiano et al, 2019 |        10 |        32 |                       0.000000 |                        0.000000 | Jiang et al, 2015;Snyder et al, 2016     | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| InFold_ComBat     | Cristiano et al, 2019 |        10 |        32 |                       0.075005 |                        0.053668 | Jiang et al, 2015;Snyder et al, 2016     | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| Proposed_LSN_ICC  | Cristiano et al, 2019 |        10 |        32 |                       0.000000 |                        0.000000 | Jiang et al, 2015;Snyder et al, 2016     | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| Hybrid_MultiScale | Cristiano et al, 2019 |        10 |        32 |                       0.000000 |                        0.000000 | Jiang et al, 2015;Snyder et al, 2016     | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| Raw_Ratio         | Jiang et al, 2015     |        10 |        32 |                       0.000000 |                        0.000000 | Cristiano et al, 2019;Snyder et al, 2016 | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| Global_Zscore     | Jiang et al, 2015     |        10 |        32 |                       0.000000 |                        0.000000 | Cristiano et al, 2019;Snyder et al, 2016 | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| InFold_ComBat     | Jiang et al, 2015     |        10 |        32 |                       0.214129 |                        0.105698 | Cristiano et al, 2019;Snyder et al, 2016 | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| Proposed_LSN_ICC  | Jiang et al, 2015     |        10 |        32 |                       0.000000 |                        0.000000 | Cristiano et al, 2019;Snyder et al, 2016 | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |
| Hybrid_MultiScale | Jiang et al, 2015     |        10 |        32 |                       0.000000 |                        0.000000 | Cristiano et al, 2019;Snyder et al, 2016 | diagnostic_target_labels_used_only_to_construct_real_patient_composition_not_to_tune_models |

| representation   |      rho |        p |   windows |   median_eta2_batch |   median_eta2_bio | scope                                                                            |        q |
|:-----------------|---------:|---------:|----------:|--------------------:|------------------:|:---------------------------------------------------------------------------------|---------:|
| Raw_Ratio        | 0.606342 | 0.000000 |       477 |            0.655450 |          0.008406 | descriptive_coverage_and_diagnosis_adjusted_unique_variance_not_CV_preprocessing | 0.000000 |
| Global_Zscore    | 0.158177 | 0.000525 |       477 |            0.025532 |          0.011189 | descriptive_coverage_and_diagnosis_adjusted_unique_variance_not_CV_preprocessing | 0.000525 |
| LSN              | 0.404597 | 0.000000 |       477 |            0.024555 |          0.009502 | descriptive_coverage_and_diagnosis_adjusted_unique_variance_not_CV_preprocessing | 0.000000 |

## 数据完整性

|   patients |   feature_columns |   overall_feature_missing_fraction |   max_column_missing_fraction |
|-----------:|------------------:|-----------------------------------:|------------------------------:|
| 585.000000 |       2739.000000 |                           0.000000 |                      0.000000 |

临床字段缺失率：

| field              |   missing_fraction |
|:-------------------|-------------------:|
| type               |           0.000000 |
| patient_id         |           0.000000 |
| specimen_id        |           0.000000 |
| timepoint          |           0.000000 |
| bam_id             |           0.000000 |
| study              |           0.000000 |
| sample_type        |           0.000000 |
| stage              |           0.565812 |
| library_kit        |           0.000000 |
| extraction_kit     |           0.000000 |
| seq_platform       |           0.000000 |
| age                |           0.292308 |
| gender             |           0.206838 |
| clinical_tf        |           1.000000 |
| data_source        |           0.000000 |
| seq_depth          |           0.000000 |
| mito_frac          |           1.000000 |
| ichorcna_tf        |           0.001709 |
| ichorcna_ploidy    |           0.001709 |
| ichorcna_gender    |           1.000000 |
| sample_id          |           0.000000 |
| center             |           0.000000 |
| label              |           0.000000 |
| cancer_type        |           0.000000 |
| fragment_url       |           0.000000 |
| center_source      |           0.000000 |
| label_source       |           0.000000 |
| study_id           |           0.000000 |
| publication        |           0.000000 |
| genome             |           0.000000 |
| genome_source      |           0.000000 |
| verified           |           0.000000 |
| sex                |           0.206838 |
| sequencing_depth   |           0.000000 |
| accepted_fragments |           0.000000 |
| center_id          |           0.000000 |
| disease_status     |           0.000000 |

ComBat pooled AUROC=0.543，加权中心内配对AUC=0.669，跨中心配对AUC=0.366。其中心内逐窗口Cohen d最大变化仅3.04e-08；不能据低pooled AUROC声称生物信号被摧毁。

## 设计边界与标题

单一来源无法估计多中心ICC/ComBat/DANN；ICC采用明确标记的LSN回退，ComBat和DANN相应单元格不估计。样本内变换不依赖目标批次构成，不等于分类器、训练患病率、PPV或临床泛化天然免疫混杂。

现稿贡献定位为Scale Trade-off系统审计，Hybrid是融合对照，不宣称击败现有基线。CNA保留及GC驱动漂移消除尚无独立因果验证；研究范围为公开cfDNA队列中的跨研究泛化。Study Source同时包含实验批次与人群构成混杂。

Combined LOCO的固定ElasticNet比较：Hybrid AUROC=0.731，Global=0.759，Proposed_LSN_ICC=0.724。这些点估计不构成融合算法的普遍显著优越性。

## 本地来源资格核查

| study                     |   runs |   healthy_patients |   cancer_patients | meets_cancer_and_10_healthy   |
|:--------------------------|-------:|-------------------:|------------------:|:------------------------------|
| Adalsteinsson et al, 2017 |   1508 |                  0 |               528 | False                         |
| Cristiano et al, 2019     |    538 |                262 |               231 | True                          |
| Jiang et al, 2015         |    225 |                 32 |                90 | True                          |
| Snyder et al, 2016        |     60 |                  4 |                48 | False                         |
| Sun et al, 2019           |     31 |                  0 |                12 | False                         |

本地未核验到额外符合癌症+至少10名健康人标准的来源；此结论仅限历史快照与存档来源汇总，不代表当前数据库或所有公开数据集的完整普查。两个双类别且健康人足量来源+一个病例主导来源是核心局限；不将双类别充足误称为严格1:1均衡。

完整手稿：../MANUSCRIPT_FULL_TEXT.md；所有统计量来自本地真实CSV，代码仓库、作者、资助和伦理确认项保持明确待作者填写。
