# 最终严谨性口径对齐

现稿题目：Scale Trade-offs in Within-Sample Self-Normalization for Cross-Study Generalization in Public cfDNA Cohorts

本轮只更新文稿、元数据资格核查和报告。未重复特征提取、分类器训练、参数搜索或图表重绘；核心输入、预测、Table3、配对检验和图件哈希全部保持一致。

## 贡献定位与对比口径

以下为同一固定ElasticNet、Combined LOCO、原有独立患者数据的结果，不混用主双向和三研究汇总口径：

| pipeline          |   n |   auroc |   auroc_low |   auroc_high |
|:------------------|----:|--------:|------------:|-------------:|
| Global_Zscore     | 585 |   0.759 |       0.721 |        0.791 |
| Hybrid_MultiScale | 585 |   0.731 |       0.694 |        0.767 |
| Proposed_LSN_ICC  | 585 |   0.724 |       0.685 |        0.761 |

Hybrid的点估计低于Global，高于固定LSN-ICC；后者的小幅数值差异不代表融合稳定显著优于基线。没有对不同已拟合模型拼接的pooled预测执行DeLong；逐来源的配对检验保持原限定范围：

| cohort_scope      | center                | comparator       |      delta |          p |          q |
|:------------------|:----------------------|:-----------------|-----------:|-----------:|-----------:|
| Primary_Bilateral | Cristiano et al, 2019 | Global_Zscore    |  -0.023726 |   0.001042 |   0.003301 |
| Primary_Bilateral | Cristiano et al, 2019 | Proposed_LSN_ICC |   0.001185 |   0.901590 |   0.951678 |
| Primary_Bilateral | Jiang et al, 2015     | Global_Zscore    |   0.018258 |   0.083436 |   0.144117 |
| Primary_Bilateral | Jiang et al, 2015     | Proposed_LSN_ICC |   0.061096 |   0.009054 |   0.022936 |
| Combined_LOCO     | Cristiano et al, 2019 | Global_Zscore    |  -0.049998 |   0.000277 |   0.000937 |
| Combined_LOCO     | Cristiano et al, 2019 | Proposed_LSN_ICC |  -0.023265 |   0.054122 |   0.108244 |
| Combined_LOCO     | Jiang et al, 2015     | Global_Zscore    |   0.014045 |   0.290383 |   0.440581 |
| Combined_LOCO     | Jiang et al, 2015     | Proposed_LSN_ICC |   0.048104 |   0.059434 |   0.113701 |
| Combined_LOCO     | Snyder et al, 2016    | Global_Zscore    | nan        | nan        | nan        |
| Combined_LOCO     | Snyder et al, 2016    | Proposed_LSN_ICC | nan        | nan        | nan        |

论文贡献定位为Global vs Local的Scale Trade-off系统审计，Hybrid作为融合探针。直接测得的是染色体臂平均片段化效应量及GC相关研究来源方差，不能将它们冒充经独立真值验证的CNA保留或GC因果批次效应消除：

| representation   |   median_abs_arm_fragmentation_d |   median_eta2_batch |      rho |           p |
|:-----------------|---------------------------------:|--------------------:|---------:|------------:|
| Raw_Ratio        |                         0.380998 |           0.65545   | 0.606342 | 3.20913e-49 |
| Global_Zscore    |                         0.270255 |           0.0255316 | 0.158177 | 0.000525459 |
| LSN              |                         0.244419 |           0.0245551 | 0.404597 | 3.25825e-20 |

Local并非在所有GC指标上都更好；没有为满足“权衡机制”叙事而反写实际关联方向。

## 临床口径

研究范围统一为cross-study generalization in public cfDNA cohorts，即公开cfDNA队列中的跨研究泛化。Study Source明确为包含实验批次、处理流程、招募和人群构成混杂的代理变量，不冒充实体医院或纯技术干预。来源论文的真实书名及临床元数据字段名称保留原意，不将研究结果外推为人群筛查有效性。

## 本地FinaleDB来源核查

本地历史快照含2362个运行记录、5个带逐例元数据的来源。计数使用同一来源内非空患者ID去重；Healthy须明确标记，不将肝硬化、移植等其他非癌状态改作健康人。

| study                     |   runs |   healthy_patients |   cancer_patients | meets_cancer_and_10_healthy   |
|:--------------------------|-------:|-------------------:|------------------:|:------------------------------|
| Adalsteinsson et al, 2017 |   1508 |                  0 |               528 | False                         |
| Cristiano et al, 2019     |    538 |                262 |               231 | True                          |
| Jiang et al, 2015         |    225 |                 32 |                90 | True                          |
| Snyder et al, 2016        |     60 |                  4 |                48 | False                         |
| Sun et al, 2019           |     31 |                  0 |                12 | False                         |

原库论文存档另列7个来源、2579个记录；其中尿液来源不满足血浆要求，产前单基因诊断来源没有本地可核验的癌症/健康资格。缺失逐例资料的健康/癌症人数保持未知，没有补成0。

本地未核验到其他同时具有癌症与至少10名健康人的候选。此结论仅限已存档材料；没有刷新实时数据库，也不能声称全球所有公开数据集都无更多合格队列。未纳入候选对象链接的实时可下载性本轮未检测。

核心局限已置于Discussion的Limitations小节：两个双类别且健康对照数足量的来源，加一个病例主导来源；这不等于两个来源都严格1:1均衡。严格双向验证仅有一个训练来源，故部分多来源去偏组件无法估计，不能据此提出普遍临床推广主张。

## 产物与验证

- 完整手稿：../MANUSCRIPT_FULL_TEXT.md；主稿：MANUSCRIPT_MASTER_DRAFT.md。
- 本地资格表：../tables/Final_Alignment_FinaleDB_Source_Census.csv。
- 原文来源清单：../tables/Final_Alignment_Original_Publication_Inventory.csv。
- 对比与尺度表：Final_Alignment_Baseline_Comparison.csv、Final_Alignment_Scale_Tradeoff.csv。
- 输入、模型结果和图件的冻结哈希见final_scope_alignment_integrity.json。

针对性验证：10项测试全部通过，0失败、0错误、0跳过；文稿生成在禁止进入诊断重拟合函数的保护下通过。测试证据见final_alignment_pytest.xml及final_alignment_validation.json。
