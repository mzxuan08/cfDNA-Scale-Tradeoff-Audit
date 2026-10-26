"""Synchronize and seal the v0.2.0 publication package without model training."""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json,re,shutil,xml.etree.ElementTree as ET
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'submission_package';FINAL=ROOT/'results/step10_final'
RELEASE=OUT/'github_zenodo_release';VERSION='0.2.0'
def read(p):return Path(p).read_text(encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,value):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(value,indent=2,ensure_ascii=False),encoding='utf-8')
def copy(source,target):target=Path(target);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)

def validate_frozen():
    frozen=json.loads(read(FINAL/'reports/frozen_inputs_sha256.json'))
    for p,h in frozen.items():assert sha(ROOT/p)==h,p
    return frozen

def sync():
    validate_frozen()
    for folder in ['src','scripts','tests','config']:
        for p in (ROOT/folder).rglob('*'):
            if not p.is_file() or any(x in p.parts for x in ['__pycache__','.pytest_cache','.git']):continue
            if p.name=='submission_author.json' or p.suffix not in ['.py','.json','.yaml','.yml','.md','.toml']:continue
            copy(p,RELEASE/p.relative_to(ROOT))
    for base in [ROOT/'results/step09_major',FINAL]:
        for sub in ['tables','figures','reports']:
            for p in (base/sub).glob('*'):
                if p.is_file() and p.suffix in ['.csv','.pdf','.png','.json','.md','.xml']:
                    copy(p,RELEASE/p.relative_to(ROOT))
    for p in (FINAL/'figures').glob('*'):
        if p.suffix in ['.pdf','.png']:copy(p,RELEASE/'results/figures'/p.name)
    for p in [ROOT/'results/MANUSCRIPT_FULL_TEXT.md',ROOT/'results/reports/MANUSCRIPT_MASTER_DRAFT.md',ROOT/'results/reports/step_10_report.md']:
        copy(p,RELEASE/p.relative_to(ROOT))
    for sub in ['tables','figures','reports']:
        for p in (FINAL/sub).glob('*'):
            if p.is_file():copy(p,OUT/'figures_and_tables/step10_final'/sub/p.name)
    author=json.loads(read(ROOT/'config/submission_author.json'))
    citation=json.loads(read(RELEASE/'CITATION.cff'))
    citation.update(version=VERSION,message='Please cite this software version and its manuscript. The v0.2.0 package is prepared locally; its archival DOI is not yet assigned.')
    citation.pop('doi',None);dump(RELEASE/'CITATION.cff',citation);copy(RELEASE/'CITATION.cff',ROOT/'CITATION.cff')
    zen=json.loads(read(RELEASE/'.zenodo.json'));zen['version']=VERSION
    zen['description']='Software v0.2.0 for the normalization-scale audit of public cfDNA cohorts. Includes source-only calibration, repeated real-patient anchors, post-hoc exploratory shared-depth and GI-versus-HCC analyses, chromosome-block inference, PERMDISP, complete numeric CSV outputs and current figures. Third-party data retain their original terms.'
    zen['related_identifiers']=[{'identifier':'10.5281/zenodo.23653766','relation':'isNewVersionOf','resource_type':'software'}]
    dump(RELEASE/'.zenodo.json',zen)
    release_status={'version':VERSION,'tag':'v'+VERSION,'status':'prepared_locally_not_published','zenodo_doi':None,
        'historical_zenodo_doi':'10.5281/zenodo.23653766','classifier_fits_this_polish':0}
    dump(ROOT/'config/release.json',release_status);copy(ROOT/'config/release.json',RELEASE/'config/release.json')
    example=RELEASE/'config/submission_author.example.json'
    if example.exists():
        obj=json.loads(read(example));obj.update(software_version=VERSION,zenodo_doi=None,zenodo_status='v0.2.0 prepared; DOI not assigned')
        dump(example,obj)
    header='''# cfDNA Scale Tradeoff Audit v0.2.0

Scale Trade-offs in Within-Sample Self-Normalization for Cross-Study Generalization
in Public cfDNA Cohorts

The software package evaluates genomic normalization scale, companion dependence
and source-threshold transport as distinct properties. Frozen outputs represent
585 source-identified participants and a 548-participant bilateral benchmark.
The shared-depth and GI-versus-HCC comparisons are post-hoc exploratory
sensitivity analyses using the same cohorts with altered target size/composition;
they are not independent external validation or confirmatory evidence.

## Version and archival status

This package is version **v0.2.0**, prepared locally for GitHub/Zenodo publication.
No Zenodo DOI has been assigned to v0.2.0. The existing DOI
https://doi.org/10.5281/zenodo.23653766 identifies the historical v0.1.0 archive.
Public repository: https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit.
This package does not imply that a v0.2.0 tag, release or new DOI is already public.

## Package contents

- `src/major_revision.py`: source-only calibration, companion resampling,
  post-hoc exploratory sensitivity analyses, chromosome-block and PERMDISP code.
- `src/evaluation.py`: threshold uncertainty and paired fixed-score intervals.
- `src/revision_plotting.py` and `src/plotting.py`: current figure generation.
- `scripts/11_major_revision_analysis.py`: methodological analyses (includes
  fitting when deliberately run for research replication).
- `scripts/12_build_major_revision_package.py`: editable publication documents.
- `scripts/13_finalize_major_revision_package.py`: historical revision packager.
- `scripts/14_final_polish.py`: frozen-result terminology/stage checks, no training.
- `scripts/15_seal_v020_package.py`: release synchronization and integrity checks.
- `results/step09_major/tables/`: complete preceding revision numeric CSVs.
- `results/step10_final/tables/`: current complete numeric CSV set and stage audit.
- `results/step10_final/figures/`: current Figure 1–4 and Supplementary S5/S6,
  each supplied as vector PDF and 300 DPI PNG.
- `results/MANUSCRIPT_FULL_TEXT.md`: current manuscript and one AI-use disclosure.

Per-patient diagnostic CSVs contain repeated public source identifiers; repeated
anchor observations are not additional independent participants. Raw sequencing
streams, specimen feature matrices, trained-model caches, credentials and local
environment caches are excluded. Compact frozen score/coordinate provenance is
included for replay; original data-use terms continue to apply.

## Reproduction

Create the environment from `environment.yml` or `requirements.txt`, using
seed 2026 and the documented Python/CUDA versions. Restore published metadata
and features with the original streaming entry points before feature-dependent
analysis. The light package alone does not contain specimen feature matrices.
Legacy internal method IDs remain for checkpoint compatibility; publication
labels are Local_LSN, Local_LSN_ICC and EB-TargetMap.

For the final polish, run `python scripts/14_final_polish.py` against the complete
local frozen workspace. It reads existing scores, stage metadata and coordinates
and does not train or select a classifier. Use the document authoring runtime,
Pandoc, python-docx and LibreOffice to run script 12 and render all documents;
script 15 synchronizes the package and seals manifests. The final packaging
helpers require the complete local publication workspace and author configuration.
Scripts 01–11 are research-replication entry points, not part of the no-training
final-polish command sequence.

## Validation and provenance

The preceding major revision recorded 93 passing tests and one optional live
download test skipped. Final-polish tests are separate read-only publication
checks and do not rerun classifiers. See `results/step10_final/reports/` for
the current results and the frozen-input checksum manifest. Estimator code,
feature matrices, predictions and prior numeric analysis tables remain frozen.
Cristiano stage counts are traced to the selected source sample identifiers;
Stage III is 34, Stage 0 is 2 and Stage X is 1, with total cancer n=218.

## License and citation

Original software is MIT licensed. Cite the original public data papers and the
software version actually used. AI assistance is documented in the manuscript;
publication disclosures are distinct from research-code authorship metadata.
'''
    header+='\nAuthor: '+author['name']+'; '+author['affiliation']+'\n\nORCID: https://orcid.org/'+author['orcid']+'\n'
    (RELEASE/'README.md').write_text(header,encoding='utf-8');copy(RELEASE/'README.md',ROOT/'README.md')
    (OUT/'figures_and_tables/README_CURRENT_REVISION.md').write_text('# Current v0.2.0 materials\n\nUse step10_final/figures and step10_final/tables for this manuscript. Other subfolders retain historical scientific outputs and submitted-version provenance.\n',encoding='utf-8')
    for name in ['04_bioRxiv_Integrated_Preprint','05_Supplementary_Material']:
        p=ROOT/'work/step08/render'/name/(name+'.pdf')
        if p.exists():copy(p,OUT/(name+'.pdf'))
    return release_status

def finalize():
    state=sync();frozen=validate_frozen()
    spec=importlib.util.spec_from_file_location('publication',ROOT/'scripts/12_build_major_revision_package.py')
    builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
    validated=builder.validate_package()
    status=json.loads(read(ROOT/'work/step08/render_status.json'))
    for row in status:assert row['returncode']==0 and row['sha256']==sha(OUT/row['file'])
    text=read(ROOT/'results/MANUSCRIPT_FULL_TEXT.md');assert 'marginal' not in text.lower()
    titles=['Study: unique Type III R²=0.656; p=0.001','Disease: unique Type III R²=0.008; p=0.001']
    figuretext=' '.join(''.join(p.extract_text() or '' for p in PdfReader(FINAL/'figures/Figure1_Batch_Vulnerability_Audit.pdf').pages).split())
    # Text extraction can introduce kerning spaces inside a PDF word.
    for title in titles:assert re.sub(r'\s+','',title) in re.sub(r'\s+','',figuretext),(title,figuretext)
    effective=[]
    for record in json.loads(read(FINAL/'reports/figure_typography.json')):
        p=FINAL/'figures'/(record['figure']+'.png')
        with Image.open(p) as im:font=record['minimum_text_pt']*6.5/(im.width/300)
        assert font>=8.5,(p,font)
        effective.append(dict(record,minimum_embedded_text_pt=font))
    python_files=[p for folder in ['src','scripts','tests'] for p in (ROOT/folder).rglob('*.py')]
    # The legitimate disclosure values stay in author configuration and prose.
    author=json.loads(read(ROOT/'config/submission_author.json'))
    brands=[author['ai_tool'],author['ai_developer']]
    for p in python_files:
        assert not any(b.lower() in read(p).lower() for b in brands),p
    suites=list(ET.parse(FINAL/'reports/nontraining_tests.xml').getroot().iter('testsuite'))
    tests={k:sum(int(s.get(k,0)) for s in suites) for k in ['tests','failures','errors','skipped']}
    assert tests['failures']==tests['errors']==tests['skipped']==0
    summary=dict(state,check_status='passed',tests=tests,prior_full_test_run={'passed':93,'skipped':1,'scope':'preceding Step 9'},
        current_test_scope='Read-only publication checks; no classifier fitting',protected_files=len(frozen),
        rendered_documents=status,pdf_pages=len(PdfReader(OUT/'04_bioRxiv_Integrated_Preprint.pdf').pages),
        pdf_sha256=sha(OUT/'04_bioRxiv_Integrated_Preprint.pdf'),classifier_fits=0,figure_titles=titles,
        minimum_embedded_figure_font_pt=min(x['minimum_embedded_text_pt'] for x in effective),
        code_files_checked=len(python_files),public_release_created=False,biorxiv_revision_uploaded=False)
    dump(OUT/'verification/Build_Summary.json',summary)
    dump(OUT/'verification/Current_Test_Results.json',tests)
    dump(OUT/'verification/Current_Test_Summary.json',dict(tests,scope=summary['current_test_scope']))
    dump(OUT/'verification/Rendered_Documents.json',status)
    dump(OUT/'verification/Final_Seal_Validation.json',summary)
    dump(OUT/'verification/Render_QA_Review.json',{'status':'passed','scope':'Step 10','method':'Rendered all eight Word files; inspected final page contact sheets and full-size revised table, formula and Figure 1 pages.'})
    dump(OUT/'verification/Frozen_Input_Integrity.json',{'status':'passed','protected_files':frozen,'classifier_fits':0,'changed_scientific_inputs':[]})
    dump(FINAL/'reports/final_seal_validation.json',summary)
    copy(FINAL/'reports/nontraining_tests.xml',OUT/'verification/Final_Polish_Tests.xml')
    checklist='''# 最终五点核验 v0.2.0

1. **术语通过**：当前全文、Methods、Results、Figure 1 图注和绘图脚本统一为 unique Type III R²。Figure 1b 两个标题分别为 Study 0.656、Disease 0.008，p 均为 0.001；已重新导出 PDF/300 DPI PNG。
2. **解释边界通过**：Results 3.1 的 PERMANOVA 报告后紧接指定英文限定句；Discussion 明确研究来源兼含技术与人群构成，不作纯质心位置或纯技术效应解释。
3. **探索性标定通过**：共同深度范围和 GI-versus-HCC 结果在 Methods、Results、Discussion 中均标为 post-hoc exploratory sensitivity analyses，明确复用原队列、改变目标样本量与构成，不构成独立外部验证或确认性证据。
4. **分期加总通过**：依据冻结队列 ID 回查本地原始元数据，Cristiano 癌症分期为 Stage 0=2、I=44、II=102、III=34、IV=35、X=1、Unknown/Missing=0，合计 218。Stage III 原始值是 34，并非 3；没有补造 34 例未知患者。正文 Table 2、CSV、Word 已同步，逐患者核对见 Cristiano_Stage_Patient_Trace.csv。
5. **版本与完整性通过**：当前软件发布包版本为 v0.2.0；README、CITATION.cff、归档元数据和 Data and code availability 已对齐。全部校准、锚点、探索性分析、区组推断、PERMDISP 脚本、完整 CSV 及 Figure 1–4/S5/S6 均纳入发布目录。

代码已去掉硬编码的工具品牌与重复图注署名；论文必需的 AI 使用声明如实保留。历史已提交稿件和历史统计记录用于可追溯，不冒充当前版本。

**发布状态**：本轮完成的是本地 v0.2.0 封版与归档包；未创建远端 GitHub Release、未发布新 Zenodo 记录、未上传 bioRxiv 修订版。历史 DOI 10.5281/zenodo.23653766 仅对应 v0.1.0，未冒充新版 DOI。
'''
    checklist+=f'\n本轮非训练校验：{tests["tests"]} 项全部通过，{tests["skipped"]} 跳过；分类器训练次数为 0；{len(frozen)} 个受保护文件哈希一致；图文 PDF 共 {summary["pdf_pages"]} 页。\n'
    (OUT/'FINAL_5_POINTS_CHECK.md').write_text(checklist,encoding='utf-8')
    (OUT/'SUBMISSION_CHECKLIST_FOR_AUTHOR.md').write_text(checklist+'\n当前文件：01 Title Page、02 Highlights、03 双倍行距连续行号盲审稿、04 图文整合稿及 PDF、05 补充材料及 PDF、06/06b 两个目标期刊投稿信、07 利益冲突声明。使用 figures_and_tables/step10_final 中的图表。\n',encoding='utf-8')
    for p in (FINAL/'reports').glob('*'):
        if p.is_file():copy(p,RELEASE/p.relative_to(ROOT));copy(p,OUT/'figures_and_tables/step10_final/reports'/p.name)
    for p in [OUT/'FINAL_5_POINTS_CHECK.md',ROOT/'results/reports/step_10_report.md']:
        copy(p,RELEASE/'results/reports'/p.name)
    archives=builder.pack.finalize_archives()
    versioned=ROOT/f'cfDNA-Scale-Tradeoff-Audit-v{VERSION}.zip'
    copy(OUT/'github_zenodo_release.zip',versioned)
    archives.update(versioned_code_zip=str(versioned),versioned_code_zip_sha256=sha(versioned))
    dump(ROOT/'work/step10_final/final_package.json',dict(summary,archives=archives))
    print(json.dumps({'version':VERSION,'tests':tests,'classifier_fits':0,'protected':len(frozen),'pages':summary['pdf_pages'],'archives':archives},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--finalize',action='store_true');args=parser.parse_args()
    finalize() if args.finalize else sync()
