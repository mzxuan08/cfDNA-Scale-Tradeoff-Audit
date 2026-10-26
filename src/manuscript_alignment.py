"""Document-only final alignment; protect all fitted results and clinical inputs."""
import json
import re
import time
from pathlib import Path
import numpy as np
import pandas as pd
from .common import sha256,write_json
from .study import export_table


def finalize_scope_alignment(cfg):
    start=time.perf_counter();root=Path(cfg['paths']['results']);project=Path(cfg['_root']);tables=root/'tables'
    paths=[Path(cfg['paths']['metadata']),Path(cfg['paths']['processed'])/cfg['extraction']['final_file'],
           root/'reports/step7_fit_diagnostics.json',Path(cfg['_config_path'])]
    for name in ['Table3_Benchmark_Metrics_Full','Step7_Paired_Significance','Step7_Nested_Parameter_Search']:
        paths.extend([tables/(name+suffix) for suffix in ['.csv','.parquet']])
    paths.append(tables/'Step7_OOF_Predictions.parquet')
    paths.extend(p for p in (root/'figures').glob('Figure*') if p.suffix in ['.pdf','.png'])
    paths.extend(project/'src'/name for name in ['models.py','debiasing.py','evaluation.py','study.py','multiscale_benchmark.py'])
    frozen={str(p):sha256(p) for p in paths}
    from .multiscale_manuscript import compile_step7_manuscript,EVIDENCE_TITLE
    pred=pd.read_parquet(tables/'Step7_OOF_Predictions.parquet')
    compile_step7_manuscript(cfg,pred,None,None,None,None,None,reuse_existing=True)
    if any(sha256(p)!=digest for p,digest in frozen.items()):raise ValueError('Text alignment modified a frozen experiment artifact')
    metrics=pd.read_csv(tables/'Table3_Benchmark_Metrics_Full.csv')
    selected=metrics[metrics.cohort_scope.eq('Combined_LOCO')&metrics.model.eq('elasticnet')&metrics.center.eq('AllCenters')
                     &metrics.pipeline.isin(['Hybrid_MultiScale','Global_Zscore','Proposed_LSN_ICC'])].copy()
    export_table(selected,cfg,'Final_Alignment_Baseline_Comparison')
    gc=pd.read_csv(tables/'Step7_Residual_Batch_GC_Associations.csv')
    arms=pd.read_csv(tables/'Step7_Arm_Level_Disease_Contrasts.csv')
    effect=arms.groupby('representation').cohen_d.apply(lambda s:s.abs().median()).rename('median_abs_arm_fragmentation_d')
    scales=gc.merge(effect,on='representation',validate='one_to_one')
    scales['CNA_preservation_validated']=False;scales['GC_technical_causality_validated']=False
    export_table(scales,cfg,'Final_Alignment_Scale_Tradeoff')
    census=pd.read_csv(tables/'Final_Alignment_FinaleDB_Source_Census.csv')
    broad=census[census.screen.eq('plasma_snapshot_no_depth_restriction')]
    pub=pd.read_csv(tables/'Final_Alignment_Original_Publication_Inventory.csv')
    metadata=json.loads((root/'reports/final_alignment_metadata_audit.json').read_text())
    sig=pd.read_csv(tables/'Step7_Paired_Significance.csv')
    comparisons=sig[sig.model.eq('elasticnet')&sig['method'].eq('Hybrid_MultiScale')&sig.comparator.isin(['Global_Zscore','Proposed_LSN_ICC'])]
    export_table(comparisons,cfg,'Final_Alignment_Paired_Comparisons')
    headline=selected[['pipeline','n','auroc','auroc_low','auroc_high']]
    report=f'''# 最终严谨性口径对齐

现稿题目：{EVIDENCE_TITLE}

本轮只更新文稿、元数据资格核查和报告。未重复特征提取、分类器训练、参数搜索或图表重绘；核心输入、预测、Table3、配对检验和图件哈希全部保持一致。

## 贡献定位与对比口径

以下为同一固定ElasticNet、Combined LOCO、原有独立患者数据的结果，不混用主双向和三研究汇总口径：

{headline.to_markdown(index=False,floatfmt='.3f')}

Hybrid的点估计低于Global，高于固定LSN-ICC；后者的小幅数值差异不代表融合稳定显著优于基线。没有对不同已拟合模型拼接的pooled预测执行DeLong；逐来源的配对检验保持原限定范围：

{comparisons[['cohort_scope','center','comparator','delta','p','q']].to_markdown(index=False,floatfmt='.6f')}

论文贡献定位为Global vs Local的Scale Trade-off系统审计，Hybrid作为融合探针。直接测得的是染色体臂平均片段化效应量及GC相关研究来源方差，不能将它们冒充经独立真值验证的CNA保留或GC因果批次效应消除：

{scales[['representation','median_abs_arm_fragmentation_d','median_eta2_batch','rho','p']].to_markdown(index=False,floatfmt='.6g')}

Local并非在所有GC指标上都更好；没有为满足“权衡机制”叙事而反写实际关联方向。

## 临床口径

研究范围统一为cross-study generalization in public cfDNA cohorts，即公开cfDNA队列中的跨研究泛化。Study Source明确为包含实验批次、处理流程、招募和人群构成混杂的代理变量，不冒充实体医院或纯技术干预。来源论文的真实书名及临床元数据字段名称保留原意，不将研究结果外推为人群筛查有效性。

## 本地FinaleDB来源核查

本地历史快照含{metadata['snapshot_runs']}个运行记录、{metadata['run_level_studies']}个带逐例元数据的来源。计数使用同一来源内非空患者ID去重；Healthy须明确标记，不将肝硬化、移植等其他非癌状态改作健康人。

{broad[['study','runs','healthy_patients','cancer_patients','meets_cancer_and_10_healthy']].to_markdown(index=False)}

原库论文存档另列{metadata['archived_publication_studies']}个来源、{metadata['archived_publication_records']}个记录；其中尿液来源不满足血浆要求，产前单基因诊断来源没有本地可核验的癌症/健康资格。缺失逐例资料的健康/癌症人数保持未知，没有补成0。

本地未核验到其他同时具有癌症与至少10名健康人的候选。此结论仅限已存档材料；没有刷新实时数据库，也不能声称全球所有公开数据集都无更多合格队列。未纳入候选对象链接的实时可下载性本轮未检测。

核心局限已置于Discussion的Limitations小节：两个双类别且健康对照数足量的来源，加一个病例主导来源；这不等于两个来源都严格1:1均衡。严格双向验证仅有一个训练来源，故部分多来源去偏组件无法估计，不能据此提出普遍临床推广主张。

## 产物与验证

- 完整手稿：../MANUSCRIPT_FULL_TEXT.md；主稿：MANUSCRIPT_MASTER_DRAFT.md。
- 本地资格表：../tables/Final_Alignment_FinaleDB_Source_Census.csv。
- 原文来源清单：../tables/Final_Alignment_Original_Publication_Inventory.csv。
- 对比与尺度表：Final_Alignment_Baseline_Comparison.csv、Final_Alignment_Scale_Tradeoff.csv。
- 输入、模型结果和图件的冻结哈希见final_scope_alignment_integrity.json。
'''
    (root/'reports/final_scope_alignment.md').write_text(report,encoding='utf-8')
    acceptance=root/'reports/step7_acceptance_report.md'
    if acceptance.exists():
        text=acceptance.read_text(encoding='utf-8')
        text=re.sub(r'鉴于效能方向不一致.*?Step7_Scientific_Claim_Audit.csv。',
            '现稿已定位为Scale Trade-off系统审计，采用公开cfDNA队列中的跨研究泛化口径；CNA特异性保留与GC因果去偏均未被当作已验证机制。原候选题目及证据审计保留，最新来源资格核查见final_scope_alignment.md。',text,flags=re.S)
        note='> 最终口径已更新：贡献为Scale Trade-off审计；研究范围为公开cfDNA队列中的跨研究泛化。新增来源资格核查与证据边界见 [最终对齐报告](final_scope_alignment.md)。\n\n'
        if not text.startswith('> 最终口径已更新'):text=note+text
        acceptance.write_text(text,encoding='utf-8')
    manuscript=(root/'MANUSCRIPT_FULL_TEXT.md').read_text(encoding='utf-8')
    forbidden=re.compile(r'multi[- ]center clinical early detection|multi[- ]centre clinical early detection',re.I)
    assert not forbidden.search(manuscript)
    assert 'cross-study generalization in public cfdna cohorts' in manuscript.lower()
    assert 'Scale Trade-off' in manuscript and '### 4.1. Limitations' in manuscript
    assert manuscript in (root/'reports/MANUSCRIPT_MASTER_DRAFT.md').read_text(encoding='utf-8')
    trace=pd.read_csv(tables/'Step7_Manuscript_Numeric_Traceability.csv')
    for name,group in trace.groupby('source_csv'):
        data=pd.read_csv(tables/name)
        for column,rows in group.groupby('column'):
            np.testing.assert_allclose(data.loc[rows.row_index.astype(int),column].to_numpy(dtype=float),rows.value,rtol=1e-13,atol=1e-15)
    write_json(dict(status='passed',runtime_seconds=round(time.perf_counter()-start,3),
        feature_extraction_rerun=False,classifier_training_rerun=False,diagnostic_refitting_rerun=False,
        frozen_files=frozen,numeric_source_cells_verified=len(trace),metadata_scope=metadata['scope'],
        full_manuscript_word_count=len(manuscript.split()),
        abstract_word_count=len(manuscript.split('## Abstract')[1].split('**Keywords:**')[0].split()),
        additional_eligible_local_sources=metadata['additional_qualifying_count']),root/'reports/final_scope_alignment_integrity.json')
    from .manuscript import refresh_clinical_manifest
    refresh_clinical_manifest(cfg)
    return report
