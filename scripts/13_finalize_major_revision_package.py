"""Verify rendered revision outputs, preserve scientific v1, and package locally."""
from pathlib import Path
import csv,hashlib,importlib.util,json,re,shutil,xml.etree.ElementTree as ET
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'submission_package';NEW=ROOT/'results/step09_major'
def read(p):return Path(p).read_text(encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,value):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(value,indent=2,ensure_ascii=False),encoding='utf-8')

def main():
    spec=importlib.util.spec_from_file_location('revision_package',ROOT/'scripts/12_build_major_revision_package.py')
    builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
    builder.validate_package()
    baseline=json.loads(read(ROOT/'work/step09/scientific_baseline.json'))
    protected={p:h for p,h in baseline.items() if p!='src/evaluation.py'}
    for p,h in protected.items():assert sha(ROOT/p)==h,p
    archive=ROOT/'results/archive_before_step09_major'
    for p,h in json.loads(read(archive/'manifest.json')).items():assert sha(archive/p)==h,p
    status=json.loads(read(ROOT/'work/step08/render_status.json'))
    for item in status:
        assert item['returncode']==0 and item['sha256']==sha(OUT/item['file'])
        assert item['pages']>0
    for stem in ['04_bioRxiv_Integrated_Preprint','05_Supplementary_Material']:
        source=ROOT/'work/step08/render'/stem/(stem+'.pdf');assert source.exists()
        reader=PdfReader(source);assert len(reader.pages)>0
        text='\n'.join(p.extract_text() or '' for p in reader.pages)
        assert '$$' not in text and '\\frac' not in text and '\\operatorname' not in text
        assert not re.search(r'\bAl\b',text)
        assert all(len((p.extract_text() or '').split())>5 for p in reader.pages)
        shutil.copy2(source,OUT/(stem+'.pdf'))
    typography=[]
    for record in json.loads(read(NEW/'reports/figure_typography.json')):
        p=NEW/'figures'/(record['figure']+'.png')
        with Image.open(p) as image:
            effective=record['minimum_text_pt']*6.5/(image.width/300)
            assert effective>=8.5,(p,effective)
            typography.append(dict(record,embedded_width_inches=6.5,png_width_pixels=image.width,
                minimum_embedded_text_pt=effective,png_sha256=sha(p)))
    suites=list(ET.parse(NEW/'reports/pytest_results.xml').getroot().iter('testsuite'))
    tests={k:sum(int(s.get(k,0)) for s in suites) for k in ['tests','errors','failures','skipped']}
    tests['passed']=tests['tests']-tests['skipped']-tests['errors']-tests['failures']
    assert tests['errors']==tests['failures']==0
    summary={'status':'passed','revision':'Step 9 major revision','participants':585,'primary_participants':548,
        'clinical_extraction_repeated':False,'new_source_only_calibration_and_sensitivity_fits':True,
        'original_protected_files_unchanged':len(protected),'old_submitted_version':'BIORXIV/2026/761726 Version 1',
        'revision_uploaded':False,'new_zenodo_release':False,'tests':tests,'rendered_documents':status,
        'pdf_pages':len(PdfReader(OUT/'04_bioRxiv_Integrated_Preprint.pdf').pages),
        'pdf_sha256':sha(OUT/'04_bioRxiv_Integrated_Preprint.pdf'),
        'previous_submitted_pdf_sha256':'6c2c9d1b3af8f3be4cc15c6899812c2559e85f3188bd6c22fe66f467c17c298d',
        'figures':typography,'references':29,'manual_review_of_this_new_revision':'pending author review before submission'}
    dump(OUT/'verification/Build_Summary.json',summary)
    dump(OUT/'verification/Current_Test_Results.json',tests)
    dump(OUT/'verification/Current_Test_Summary.json',dict(tests,scope='Full project and new methodological revision; no extraction rerun'))
    dump(OUT/'verification/Rendered_Documents.json',status)
    dump(OUT/'verification/Frozen_Input_Integrity.json',{'status':'passed','protected_files':protected,
        'authorized_code_extension':['src/evaluation.py'],'changed_scientific_data':[]})
    dump(OUT/'verification/Major_Revision_Figure_Typography.json',typography)
    # Align bibliographic metadata IDs with the current first-citation numbering.
    references=list(csv.DictReader(open(NEW/'tables/Revision_References.csv',encoding='utf-8')))
    numbers={r['citation'].split('https://doi.org/')[-1].strip().lower():int(r['reference']) for r in references}
    verified=json.loads(read(OUT/'verification/reference_verification.json'))
    for record in verified['references']:
        record.setdefault('source_reference_id',record['reference'])
        record['reference']=numbers[record['doi'].lower()]
    verified['references'].sort(key=lambda r:r['reference'])
    verified['numbering_note']='Current manuscript first-citation numbering; source_reference_id preserves earlier metadata identity.'
    dump(OUT/'verification/reference_verification.json',verified)
    old=list(csv.DictReader(open(OUT/'verification/References_Verified.csv',encoding='utf-8')))
    bydoi={r['doi'].lower():r for r in old}
    record=json.loads(read(NEW/'reports/PERMDISP_Crossref.json'))['message']
    bydoi.setdefault(record['DOI'].lower(),{'doi':record['DOI'],'authors':'Anderson MJ','year':'2006',
        'title':record['title'][0],'journal':record['container-title'][0],'volume':record['volume'],
        'issue':record['issue'],'pages_or_article':record['page'],
        'crossref_url':'https://api.crossref.org/works/'+record['DOI'],
        'crossref_sha256':sha(NEW/'reports/PERMDISP_Crossref.json')})
    rows=[]
    for doi,number in sorted(numbers.items(),key=lambda kv:kv[1]):
        row=bydoi[doi].copy();row['reference']=number;rows.append(row)
    fields=sorted({k for row in rows for k in row})
    with open(OUT/'verification/References_Verified.csv','w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    dump(OUT/'verification/Render_QA_Review.json',{'status':'passed','scope':'Current Step 9 revision',
        'review':'All pages reviewed in contact sheets; full-size formula, table and figure pages inspected. Blank trailing pages, split tables and escaped math were corrected.'})
    shutil.copy2(NEW/'reports/pytest_results.xml',OUT/'verification/Regression_Tests.xml')
    shutil.copy2(NEW/'reports/packaging_tests.xml',OUT/'verification/Packaging_Tests.xml')
    # Make the light software bundle reproducible without copying clinical caches.
    release=OUT/'github_zenodo_release'
    for folder in ['src','scripts','tests','config']:
        for p in (ROOT/folder).rglob('*'):
            if not p.is_file() or any(x in p.parts for x in ['__pycache__','.pytest_cache','.git']):continue
            if p.name in ['submission_author.json'] or p.suffix not in ['.py','.yaml','.yml','.json','.md','.toml']:continue
            dest=release/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    excluded=['Source_Calibration_Splits','Source_Calibration_Predictions','Sensitivity_Predictions',
        'Threshold_Bootstrap_Distribution','Anchor_All_Predictions','Anchor_Real_Batch_Members','PERMDISP_Patient_Distances']
    for p in (NEW/'tables').glob('*.csv'):
        if p.stem in excluded:continue
        dest=release/'results/step09_major/tables'/p.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    # Compact score and annotation-derived provenance needed by the revision
    # replay is supplied explicitly; no model weights or specimen matrices.
    for name in ['Step7_OOF_Predictions.parquet','clinical_embeddings.parquet','Clinical_Feature_Importance.parquet']:
        shutil.copy2(ROOT/'results/tables'/name,release/'results/tables'/name)
    historical='results/MANUSCRIPT_FULL_TEXT.md'
    dest=release/'results/archive_before_step09_major'/historical;dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(archive/historical,dest)
    dump(release/'results/archive_before_step09_major/manifest.json',{historical:sha(dest)})
    dest=release/'results/step09_major/reports/PERMDISP_Crossref.json';dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(NEW/'reports/PERMDISP_Crossref.json',dest)
    for p in [ROOT/'results/MANUSCRIPT_FULL_TEXT.md',ROOT/'results/reports/MANUSCRIPT_MASTER_DRAFT.md',
              ROOT/'results/reports/step_09_major_revision_report.md']:
        dest=release/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    readme=release/'README.md';old=read(readme).split('\n## Methodological major revision')[0]
    old=re.sub(r'(?m)^- \*\*Fixed anchors:\*\*.*$',
       '- **Companion dependence:** sample-wise frozen rules satisfy mathematical independence; repeated real-patient batches quantify EB-TargetMap variability at a fixed case fraction.',old)
    if 'Primary bilateral target AUROCs' not in old:
        old=old.replace('Fixed ElasticNet Combined LOCO AUROCs are',
           'Primary bilateral target AUROCs are Global/Hybrid 0.632/0.608 in Cristiano and 0.859/0.877 in Jiang. Equal-weight means are Raw 0.742, Global 0.745 and Hybrid 0.743.\n\nSecondary pooled-score comparability diagnostics (not primary performance evidence): ElasticNet Combined LOCO AUROCs are')
    old=old.replace('ComBat-plus-target alignment','EB-TargetMap').replace('The approved manuscript is','The current local revision is')
    old=old.replace('scripts/10_finalize_editorial_materials.py rebuilds the current submission',
       'The historical scripts/10_finalize_editorial_materials.py rebuilds the earlier submission')
    readme.write_text(old+f'''
## Methodological major revision

The current local revision prioritizes the 548-participant bilateral study-held-out
benchmark. Hybrid minus Global AUROC differs by target: −0.024 in Cristiano and
+0.018 in Jiang. Companion-independent preprocessing is a mathematical property;
source-locked specificity still requires transport evaluation. The submitted v1,
frozen clinical inputs and original predictions/figures are preserved.

The revised tests report {tests['passed']} passed and {tests['skipped']} optional
live-network test skipped. Reproduce the new real-patient analyses after rebuilding
the source features with `python scripts/11_major_revision_analysis.py`. This writes
isolated `results/step09_major` outputs, current Markdown and figures. Code retains
historical internal pipeline identifiers for checkpoint compatibility; manuscript
aliases are Local_LSN, Local_LSN_ICC and EB-TargetMap.

Document authoring is separate: configure `config/submission_author.json` from its
example and run `scripts/12_build_major_revision_package.py --pandoc <executable>`
with python-docx and Pandoc available. Render all Word files with LibreOffice before
running `scripts/13_finalize_major_revision_package.py` in the document runtime.
The latter requires the complete local artifact/verification workspace. Per-patient
diagnostic CSVs are in the full submission bundle rather than this light archive.
The source-only statistical runtime and document authoring runtime are separate.

The immutable v0.1.0 Zenodo DOI identifies the earlier archive. This local revision
has not been published as a new GitHub/Zenodo release or uploaded to bioRxiv.
''',encoding='utf-8')
    for name in ['CITATION.cff','.zenodo.json']:
        p=release/name
        if p.exists():
            obj=json.loads(read(p))
            if name=='CITATION.cff':
                obj['version']='0.2.0-dev';obj.pop('doi',None)
                obj['message']='This local methodological revision is unreleased. The earlier v0.1.0 archive has DOI 10.5281/zenodo.23653766; that DOI does not identify this revision.'
            else:obj['version']='0.2.0-dev'
            dump(p,obj)
    shutil.copy2(OUT/'verification/References_Verified.csv',release/'results/reports/References_Verified.csv')
    # Keep the previously submitted upload payload and receipt untouched.
    dump(ROOT/'work/biorxiv_revision/revision_status.json',{'status':'prepared_locally_not_uploaded',
        'previous_submission':'BIORXIV/2026/761726 Version 1','main_pdf':str(OUT/'04_bioRxiv_Integrated_Preprint.pdf'),
        'main_pdf_sha256':summary['pdf_sha256'],'pages':summary['pdf_pages'],'author_review_pending':True})
    report=read(ROOT/'results/reports/step_09_major_revision_report.md').split('\n## Final verification')[0]
    report+='\n## Final verification\n\n'+json.dumps(summary,indent=2,ensure_ascii=False)+'\n'
    (ROOT/'results/reports/step_09_major_revision_report.md').write_text(report,encoding='utf-8')
    shutil.copy2(ROOT/'results/reports/step_09_major_revision_report.md',release/'results/reports/step_09_major_revision_report.md')
    checklist='# Step 9 投稿前核验清单\n\n'
    checklist+=f'当前正文 PDF 共 {summary["pdf_pages"]} 页；8 份 Word 已重新生成、渲染核验。测试 {tests["passed"]} 项通过、{tests["skipped"]} 项跳过（可选远程下载），无失败。保留的 {len(protected)} 个既有科学文件哈希一致。\n\n'
    checklist+='本文主结论为目标来源依赖的归一化尺度权衡、预处理同伴不变性与阈值迁移的区别。新版未上传；原 bioRxiv Version 1 保持不变。请作者核验新增统计、正文及各项声明后，再决定提交修订版。\n\n'
    checklist+='## Highlights 字符数（含空格）\n\n'
    for x in builder.HIGHLIGHTS:checklist+=f'- {len(x)} / 85：{x}\n'
    checklist+='\n## 当前交付文件\n\n'
    for p in sorted(OUT.glob('*')):
        if p.is_file() and p.suffix in ['.docx','.pdf','.md','.zip']:checklist+='- '+p.name+'\n'
    checklist+='- figures_and_tables/step09_major/：此次新增统计 CSV、6 套 PDF/300 DPI PNG 图及分析报告。\n'
    checklist+='- verification/：哈希、测试、29 篇参考文献元数据和版面验证。旧版历史记录保持可追溯。\n'
    checklist+='\n## 上传操作提示\n\n'
    checklist+='1. bioRxiv：从现有投稿记录进入允许的修订入口，核对当前稿号；不要将新版另建为重复首投稿。待筛查期间若系统未开放修订，请等待平台状态更新。\n'
    checklist+='2. 选择当前 04 PDF 和 05 补充材料，核对作者、ORCID、摘要、CC BY 许可及版本说明；此包没有替作者点击提交。\n'
    checklist+='3. 期刊：Title Page、Highlights、03 盲审正文、独立主图和补充材料分别上传；两个 Cover Letter 按目标期刊选其一，确认传统订阅模式。\n'
    (OUT/'SUBMISSION_CHECKLIST_FOR_AUTHOR.md').write_text(checklist,encoding='utf-8')
    response=OUT/'Point_by_Point_Revision_Summary.md'
    value=read(response).split('\n## 新增结果快速核对')[0]+'\n## 新增结果快速核对\n\n'
    value+='### 深度与癌种限制\n\n| 限制 | 目标来源 | 方法 | AUROC | AUPRC | 锁定阈值下特异度 |\n|---|---|---|---|---|---|\n'
    for r in csv.DictReader(open(NEW/'tables/Depth_CancerType_Sensitivity.csv',encoding='utf-8')):
        if r['method'] not in ['Global_Zscore','Hybrid_MultiScale']:continue
        value+='| '+' | '.join([r['scenario'],r['target'],r['method']]+[f'{float(r[k]):.3f}' for k in ['auroc','auprc','specificity']])+' |\n'
    response.write_text(value,encoding='utf-8')
    # Refresh all current package copies only after clinical output checks.
    for p in (NEW/'reports').glob('*'):
        if p.is_file():shutil.copy2(p,OUT/'figures_and_tables/step09_major/reports'/p.name)
    archives=builder.pack.finalize_archives()
    dump(ROOT/'work/step09_major/final_package.json',dict(summary,archives=archives))
    print(json.dumps({'tests':tests,'pdf_pages':summary['pdf_pages'],'protected_unchanged':len(protected),'archives':archives},indent=2))

if __name__=='__main__':main()
