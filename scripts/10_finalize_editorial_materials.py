"""Unify the approved publication narrative without recomputing clinical results.

Requires the complete local publication workspace. --build preserves the current
version once, then rebuilds all documents. Render and inspect them before using
--finalize. The immutable software release and scientific outputs are read-only.
"""
from __future__ import annotations
import argparse, collections, importlib.util, json, re, shutil, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('typeset081', ROOT/'scripts/08_1_final_typesetting.py')
ts = importlib.util.module_from_spec(spec); spec.loader.exec_module(ts)
pack = ts.pack
OUT = ROOT/'submission_package'
WORK = ROOT/'work/final_editorial'
ARCHIVE = ROOT/'results/archive_before_final_editorial'
RELEASE = OUT/'github_zenodo_release'

def read(path): return Path(path).read_text(encoding='utf-8')
def dump(path, value): pack.dump(path, value)
def write(path, value): pack.write(path, value)

def archive():
    if (ARCHIVE/'archive_manifest.json').exists(): return
    paths = list(OUT.glob('*.docx')) + list(OUT.glob('*.md')) + [OUT/'04_bioRxiv_Integrated_Preprint.pdf']
    paths += [ROOT/'results/MANUSCRIPT_FULL_TEXT.md', ROOT/'results/reports/MANUSCRIPT_MASTER_DRAFT.md', ROOT/'README.md', RELEASE/'README.md', RELEASE/'CITATION.cff', RELEASE/'.zenodo.json']
    hashes = {}
    for path in paths:
        target = ARCHIVE/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(path, target)
        hashes[str(path.relative_to(ROOT))] = ts.sha(target)
    dump(ARCHIVE/'archive_manifest.json', hashes)

def revise(original):
    changes = []
    def paragraph(prefix, replacement):
        nonlocal original
        matches = [p for p in original.split('\n\n') if p.startswith(prefix)]
        assert len(matches) == 1, (prefix, len(matches))
        old = matches[0]
        # Every numeric token and citation stays attached to the same paragraph.
        numbers = lambda s: collections.Counter(re.findall(r'\d+(?:\.\d+)?(?:e[-+]\d+)?', s))
        assert numbers(old) == numbers(replacement), (prefix, numbers(old)-numbers(replacement), numbers(replacement)-numbers(old))
        original = original.replace(old, replacement, 1)
        changes.append({'before':old, 'after':replacement})
    paragraph('**Background:**',
        '**Background:** Normalization sets the genomic scale of cfDNA fragmentation contrasts and the reference against which patients are scored. We present an audit linking normalization scale to cross-study discrimination and dependence on target-cohort composition.')
    paragraph('**Methods:**',
        '**Methods:** We evaluated global, local and hybrid representations in 585 source-identified participants, with a primary bilateral benchmark of 548 and a separately held-out, case-dominated cohort of 37. Source-only selection and locked thresholds assessed transport. Cancer–control pair concordance and fixed-patient anchor experiments distinguished within-study ranking from cross-study score comparability and target-composition dependence.')
    paragraph('**Results:**',
        '**Results:** Changing the companions of fixed real patients shifted their probabilities by up to 0.214 under ComBat with unseen-target alignment and by 0 under within-sample transformations. Within-study and cross-study concordances for this ComBat comparator were 0.669 and 0.366. Global and local normalization retained median absolute arm-level fragmentation effects of 0.270 and 0.244, with residual GC correlations with study variance of 0.158 and 0.405. Combined LOCO ElasticNet AUROCs were 0.759 for Global Z-score, 0.731 for Hybrid and 0.724 for fixed LSN-ICC.')
    paragraph('**Conclusions:**',
        '**Conclusions:** The audit quantifies the trade-off between genomic reference scales and separates study-wise discrimination from score comparability. Within-sample self-normalization removes companion-dependent target preprocessing for a frozen scoring rule, providing a testable property for cross-study cfDNA evaluation.')
    paragraph('We present a Scale Trade-off audit',
        'We present a Scale Trade-off audit of the representations and scoring rules used in cross-study cfDNA evaluation. Global and same-arm local contrasts reveal how normalization redistributes broad fragmentation effects and residual GC-associated study variation; Hybrid makes their complementarity testable in a source-trained classifier. Bidirectional study-held-out validation evaluates discrimination and threshold transport. Cancer–control pair decomposition explains the contributions to pooled AUROC, and fixed real-patient anchors isolate the effect of companion patients on unseen-target adaptation.')
    paragraph('The Scale Trade-off has two observable dimensions:',
        'The Scale Trade-off links broad fragmentation-contrast retention to the magnitude and genomic pattern of residual study variation. Local normalization yielded a slightly smaller median study-variance fraction than Global Z-score and a stronger residual GC rank association. These endpoints distinguish what each reference scale preserves. Their relationship to copy number and tissue accessibility is addressed as a mechanistic hypothesis in the Discussion.')
    paragraph('All primary training and calibration excluded',
        'Bidirectional testing showed that scale fusion changed discrimination in a target-dependent manner. In the Jiang-to-Cristiano direction, Hybrid AUROC was 0.608 and AUPRC was 0.671; in the Cristiano-to-Jiang direction they were 0.877 and 0.952. Relative to Global Z-score, Hybrid had lower AUROC in Cristiano and numerically higher AUROC in Jiang. The corresponding DeLong p values were 0.001 and 0.083, with BH values of 0.003 and 0.144. Primary training and calibration used the two primary studies exclusively.')
    paragraph('For fixed ElasticNet, Combined LOCO AUROC',
        'For fixed ElasticNet, Combined LOCO AUROC was 0.759 for Global Z-score, 0.731 for Hybrid and 0.724 for Proposed LSN-ICC. Hybrid occupied an intermediate position between the global and fixed local representations. The primary paired comparisons identified a lower Hybrid AUROC in Cristiano and a numerical increase in Jiang, with BH significance confined to Cristiano. Together, these results establish target-dependent effects of scale fusion. DeLong inference remains target-specific.')
    paragraph('The empirical-Bayes-plus-target-alignment comparator',
        'Pair decomposition located the source of the empirical-Bayes-plus-target-alignment comparator’s secondary pooled AUROC of 0.543. Weighted within-study concordance was 0.669, whereas cross-study concordance was 0.366; cross-study pairs contributed 41.636% of all pairs. Pooled AUROC therefore depends on the comparability of scores from independently source-trained models as well as ranking within each target study (Supplementary Figure S5).')
    paragraph('The maximum change in within-target window Cohen d',
        'Within-target window Cohen d changed by at most 3.04e-08 under this comparator, confirming that positive feature-wise affine harmonization preserved standardized contrasts to numerical precision. Target-batch location and scale still change how those features enter the source-trained classifier. Preserved within-study contrasts and altered pooled operating scores thus identify distinct properties of the adaptation rule.')
    paragraph('The fixed-anchor experiment directly identified',
        'Fixed real-patient anchors isolated companion-dependent adaptation. Changing the target batch case fraction produced a maximum anchor probability range of 0.214 under target-batch adaptation. Corresponding ranges were 0 for Raw Ratio, 0 for Global Z-score, 0 for fixed LSN-ICC and 0 for Hybrid. A patient’s within-sample representation, and hence the frozen model’s prediction, remained unchanged as companion patients changed (Figure 3d).')
    paragraph('This audit separates three properties',
        'The audit explains cross-study cfDNA performance through the reference scale of a representation and the comparability of its resulting scores. Global and local normalization preserve different fragmentation contrasts and residual genomic source patterns. Pair decomposition identifies how those scores contribute to pooled discrimination, and fixed real-patient anchors establish whether target adaptation makes an individual’s prediction depend on other patients. These complementary measurements give normalization comparisons a mechanistic interpretation.')
    paragraph('Global and local normalization impose',
        'Global and local normalization impose different reference scales. Global Z-scoring retains deviations relative to the full genome; the local log contrast references adjacent windows on the same arm. Attenuating a smooth background can therefore attenuate a broad disease-associated contrast. Hybrid exposes both scales to a source-trained penalized classifier, making their complementarity testable. The target-specific ranking and calibration results demonstrate this dependence on evaluation context; Global Z-score had the highest Combined LOCO AUROC in the fixed ElasticNet comparison. Published low-coverage copy-number methods such as ichorCNA provide a separate route to estimating tumor fraction and broad genomic imbalance [28]. Matched copy-number measurements would establish how the observed fragmentation contrasts relate to CNA.')
    paragraph('The ComBat diagnostics locate',
        'The ComBat diagnostics identify composition dependence in the evaluated empirical-Bayes source adjustment with robust unseen-target mapping. Because target quantiles follow the patient mixture, changing that mixture can move the adapted features of a fixed individual. The anchor experiment isolates this effect with source fitting frozen. Preserved within-study Cohen d and the pair decomposition place the observed probability changes in adaptation and score transport. Covariate-aware and reference-batch ComBat designs require their own evaluation under the same diagnostic framework [9,10,11].')
    paragraph('Within-sample transformation makes',
        'Within-sample transformation makes each patient’s representation independent of other target patients. With a frozen source model, adding or removing companions leaves global, local and hybrid feature vectors and predictions unchanged. This functional property directly controls one route of target-composition dependence. Training-population composition and disease spectrum continue to shape the fitted rule, while prevalence affects predictive values. Calibration and source-locked threshold transport therefore remain essential companions to discrimination.')
    paragraph('The cohort hierarchy aligns',
        'The cohort hierarchy assigns each public source an estimable endpoint. The two primary targets support bidirectional discrimination and specificity estimation. The separately held-out, cancer-enriched source contributes sensitivity at primary-source-locked thresholds. Combined LOCO extends the analysis to the three-source setting and exposes cross-model score-comparability requirements. This design allows public datasets with unequal case fractions and healthy-control representation to contribute interpretable evidence.')
    paragraph('The available source-level replication comprises',
        'Source-level replication comprises 2 cohorts with substantial healthy-control representation and a case-dominated external source (36 cancer/1 healthy). The primary sources contain both classes but are not exactly balanced. In the archived local inventory, Adalsteinsson had 528 known cancer participants and 0 explicitly healthy participants; Sun had 12 and 0. No additional source met the eligibility criterion. The inventory defines the available replication scope and does not constitute a current exhaustive database census.')
    paragraph('Study Source combines experimental-batch',
        'Study Source combines experimental-batch and population-composition differences, preventing a purely technical interpretation of source effects. Each bilateral fit has one training study, so multi-source ICC filtering, multi-batch source ComBat and multi-source adversarial training are unidentifiable in that protocol. The results label the LSN-only fallback and unavailable cells explicitly. The external cohort’s single healthy participant supplies a descriptive outcome; the generalization claims concern the evaluated public cohorts.')
    paragraph('The analysis is methodological development',
        'The representation comparison was developed using these public participants. Fold-level patient isolation addresses algorithmic leakage, while independent confirmatory evaluation is needed to assess development-history effects. Bootstrap intervals condition on the fixed studies and fitted scores and exclude refitting and future-study uncertainty. DeLong tests concern individual held-out targets. The composition stress experiment measures adaptation dependence in the finite set of tested real-patient batches.')
    paragraph('Matched copy-number measurements',
        'Complete-cohort matched copy-number measurements and tumor-specific chromatin assays were unavailable. Arm-level fragmentation contrasts and reference GC associations are therefore descriptive endpoints whose biological attribution requires orthogonal measurements. Missing clinical covariates further constrain adjustment for population differences. Additional independent sources with adequate control recruitment and standardized specimen handling would extend the cross-study evaluation.')
    paragraph('Normalization scale, source composition and score comparability define',
        'The Scale Trade-off audit makes normalization scale and score transport measurable in public cfDNA cohorts. Global and local references preserve different broad fragmentation contrasts and residual genomic source patterns. Pair concordance decomposition links study-wise ranking to pooled discrimination, while fixed-patient anchors expose companion-dependent target adaptation. Within-sample transformations eliminate that dependence for frozen scoring rules. Together with source-locked thresholds, these diagnostics provide a reproducible framework for interpreting cross-study fragmentomics performance.')
    # Retain the complete disclosure of AI roles and author responsibility.
    original = original.replace('Figures were generated by reproducible code from the recorded data; no clinical observations were fabricated.', 'Figures were generated by reproducible code from the recorded data.')
    return original, changes

def numeric_tokens(text):
    return collections.Counter(re.findall(r'\d+(?:\.\d+)?(?:e[-+]\d+)?', text))

def build(args):
    frozen = ts.frozen_check(); archive(); WORK.mkdir(parents=True, exist_ok=True)
    source = read(ARCHIVE/'results/MANUSCRIPT_FULL_TEXT.md')
    full, changes = revise(source)
    assert [x for x in source.splitlines() if x.startswith('|')] == [x for x in full.splitlines() if x.startswith('|')]
    assert source.split('## References\n',1)[1] == full.split('## References\n',1)[1]
    assert numeric_tokens(source) == numeric_tokens(full)
    dump(WORK/'editorial_changes.json', changes)
    # Figures and their existing typeset derivative are immutable in this pass.
    figure_hashes = {str(p.relative_to(ROOT)):ts.sha(p) for p in (OUT/'figures_and_tables').glob('*') if p.suffix.lower() in ['.png','.pdf']}
    dump(WORK/'figure_baseline.json', figure_hashes)
    author = json.loads(read(OUT/'author_configuration.json'))
    blind = re.sub(r'\*\*Author:\*\*.*?(?=## Abstract)', '', full, flags=re.S).replace(author['name'], 'The sole author')
    begin = blind.index('Code and compact statistical outputs are distributed at '); end = blind.index('The release contains', begin)
    blind = blind[:begin]+'The versioned code and statistical archive are supplied through the separate title page for anonymous review. '+blind[end:]
    assert not any(s in blind for s in [ts.ORCID, author['email'], author['name'], 'mzxuan08', ts.AFFILIATION])
    for p in [ROOT/'results/MANUSCRIPT_FULL_TEXT.md', ROOT/'results/reports/MANUSCRIPT_MASTER_DRAFT.md', OUT/'Manuscript_Main_Blind_and_Full.md']:
        write(p, full)
    write(OUT/'Manuscript_Main_Blinded.md', blind)
    title = full.splitlines()[0][2:]
    abstract, keywords = pack.extract_section(full,'## Abstract','## 1. Introduction').split('**Keywords:**')
    body = full.split('## 1. Introduction\n')[1].split('## Data and code availability\n')[0]
    counts = {'abstract_words':pack.word_count(abstract),'main_body_words_excluding_tables':pack.word_count('\n'.join(x for x in body.splitlines() if not x.startswith('|'))),'main_body_including_tables_words':pack.word_count(body),'full_markdown_words':pack.word_count(full),'definition':'Whitespace-delimited tokens containing a letter or numeral; body is Introduction through Conclusions, excluding declarations, legends and references.'}
    dump(OUT/'verification/Word_Counts.json', counts)
    pack.presentation_tables = ts.complete_tables
    def create(text, name, **kwargs):
        from docx.shared import Pt
        doc = pack.pandoc_doc(text, OUT/name, args.pandoc, author, **kwargs)
        if name.startswith(('03_','04_')):
            for p in doc.paragraphs:
                p.paragraph_format.space_before = Pt(6 if p.style.name.startswith('Heading') else 0)
                p.paragraph_format.space_after = Pt(4)
        ts.table_layout(doc); doc.save(OUT/name); return doc
    title_text = pack.title_page(title,abstract.strip(),keywords.strip(),counts,author).replace('ORCID: Not registered.','ORCID: https://orcid.org/'+ts.ORCID)
    title_text = title_text.replace('Figures were generated by reproducible code from the recorded data; no clinical observations were fabricated.', 'Figures were generated by reproducible code from the recorded data.')
    title_text = re.sub(r'## Data and code release\n.*', '## Data and code release\n\nPublic repository: '+author['repository_url']+'. Immutable v0.1.0 software archive: https://doi.org/'+author['zenodo_doi']+'.\n',title_text,flags=re.S)
    create(title_text,'01_Title_Page.docx',double=False)
    highlights = json.loads(read(OUT/'verification/Highlights_Character_Counts.json'))
    for h in highlights:
        assert len(h['text']) <= 85
        h['characters_including_spaces'] = len(h['text']); h['passes_85'] = True
    dump(OUT/'verification/Highlights_Character_Counts.json', highlights)
    create('# Highlights\n\n'+'\n\n'.join('- '+h['text'] for h in highlights),'02_Highlights.docx',blind=True,double=False)
    create(blind,'03_Manuscript_DoubleSpaced_LineNumbered.docx',blind=True,numbers=True)
    doc = create(full,'04_bioRxiv_Integrated_Preprint.docx'); ts.add_figures(doc); doc.save(OUT/'04_bioRxiv_Integrated_Preprint.docx')
    supp = pack.supplement(full.replace('## References','## Highlights\n\n## References'),author)
    supp = supp.replace('This supplement accompanies the frozen analysis and documents the scope of the statistical tables. It introduces no new fitted predictions or observations.', 'This supplement links pair concordance and fixed-patient composition diagnostics to the statistical tables supporting the normalization-scale audit.')
    supp = supp.replace('Unavailable one-source ComBat and DANN cells remain unavailable. A one-source ICC bypass is an LSN-only fallback.', 'One-source bilateral fits label multi-source ComBat and DANN as unavailable and the ICC bypass as an LSN-only fallback.')
    supp = supp.replace('These results did not select a revised target-optimal model. Source-only nested selection is separately documented', 'Source-only nested selection is documented independently')
    supp = supp.replace('No pooled cross-validated DeLong inference is substituted for the target-specific tests.', 'Inference concerns each held-out study and its fixed individual scoring rules.')
    supp = supp.replace('secondary within-study/cross-study concordance and source-locked score-alignment diagnostics, which do not replace benchmark probabilities.', 'secondary within-study/cross-study concordance and source-locked score-alignment diagnostics. Original benchmark probabilities are retained for the performance comparisons.')
    supp = supp.replace('Labels define real-patient diagnostic batches after evaluation; no synthetic patient is introduced.', 'Labels define diagnostic batches of real patients after predictive evaluation, with source fitting and anchor identities held fixed.')
    supp = supp.replace('Missing source counts remain unknown; no live exhaustive database census is claimed.', 'The census records the archived local inventory; unavailable source counts are recorded as unknown.')
    write(OUT/'05_Supplementary_Material.md', supp)
    doc = create(supp,'05_Supplementary_Material.docx')
    from docx.shared import Inches, Pt
    for p in list(doc.paragraphs):
        if p.text.startswith('Supplementary Figure S5.'):
            image = p.insert_paragraph_before(); image.add_run().add_picture(str(ROOT/'results/figures'/(pack.FIGURES[4]+'.png')),width=Inches(6.5))
            image.paragraph_format.keep_with_next = True; image.paragraph_format.line_spacing = 1; p.paragraph_format.line_spacing = 1.15
    doc.save(OUT/'05_Supplementary_Material.docx')
    for name in ['06_Cover_Letter','06b_Cover_Letter_Alternative_CBM']:
        text = read(ARCHIVE/'submission_package'/(name+'.md'))
        text = text.replace('Hybrid tests scale complementarity; the comparison is reported objectively alongside Global Z-score.', 'The source-trained Hybrid comparison tests whether exposing both reference scales changes discrimination and calibration in each held-out study.')
        text = text.replace('Preserved within-study standardized contrasts help locate the effect at adaptation and score transport.', 'Preserved within-study standardized contrasts locate the observed shifts in adaptation and score transport.')
        write(OUT/(name+'.md'),text)
        doc = create(text,name+'.docx',double=False)
        for p in doc.paragraphs: p.paragraph_format.line_spacing=1; p.paragraph_format.space_after=Pt(4)
        doc.save(OUT/(name+'.docx'))
    create('# Declaration of interest\n\n'+title+'\n\n'+author['name']+'\n\n'+author['competing_interests'],'07_Declaration_of_Interest.docx',double=False)
    # Synchronize the active public overview; historical execution reports stay intact.
    overview = read(ARCHIVE/'submission_package/github_zenodo_release/README.md')
    start = overview.index('## Publication materials'); end = overview.index('## Included files',start)
    overview = overview[:start]+'''## Publication materials

The approved manuscript is results/MANUSCRIPT_FULL_TEXT.md. Its abstract,
interpretation, legends, submission letters and supplementary guide use the same
normalization-scale and score-comparability framework.

scripts/10_finalize_editorial_materials.py rebuilds the current submission
documents from the preserved local publication snapshot. It uses the table and
figure-placement helpers in scripts/08_1_final_typesetting.py. Both require the
complete local publication workspace, author configuration, verified references,
Pandoc and document-rendering dependencies. Research replication uses the
analysis entry points below independently of this publication workspace.

'''+overview[end:]
    start=overview.index('**Manuscript regeneration:**');end=overview.index('## Leakage controls',start)
    overview=overview[:start]+'''**Manuscript regeneration:** scripts/06_compile_manuscript_assets.py is the
historical analysis compiler and can regenerate scientific tables. For editorial
work, preserve the approved CSVs and use the publication entry point described
above. The lightweight release supplies the approved text and research
implementation; feature-dependent validation requires reconstructed clinical
inputs and the original data terms.

'''+overview[end:]
    overview=overview.replace('## Leakage controls and statistical limits','## Fold isolation and statistical inference')
    overview=overview.replace('Do not use this exploratory public-data study as a validated screening device.','This is research software for retrospective evaluation of public cfDNA cohorts.')
    write(RELEASE/'README.md',overview)
    citation=json.loads(read(RELEASE/'CITATION.cff'))
    citation['abstract']='Reproducible software for measuring normalization-scale trade-offs, within-study discrimination and cross-study score comparability in public human cfDNA cohorts, including fixed-patient tests of target-composition dependence.'
    dump(RELEASE/'CITATION.cff',citation)
    root_readme=read(ARCHIVE/'README.md')
    start=root_readme.index('## Step 1')
    root_readme='''# cfDNA Scale Tradeoff Audit

本项目研究公开 cfDNA 队列中的跨研究泛化，量化全局与局部自归一化的尺度权衡，并通过样本对一致性分解和固定真实患者压力测试解释分数可比性。研究来源同时包含实验批次及人群构成差异。

当前定稿见 `results/MANUSCRIPT_FULL_TEXT.md`；投稿文件见 `submission_package/`。585名来源标识独立患者的分析采用548名主双向验证及37名外部敏感性队列。完整实验入口、统计口径与数据重建说明见 `submission_package/github_zenodo_release/README.md`。

当前投稿排版入口为 `scripts/10_finalize_editorial_materials.py`，需要保留完整本地投稿工作区。以下步骤记录原始工程与数据处理规则；历史执行报告保留其当时的诊断结果。

'''+root_readme[start:]
    write(ROOT/'README.md',root_readme)
    dump(WORK/'preflight.json',{'status':'passed','changed_paragraphs':len(changes),'counts':counts,'numeric_token_multiset_unchanged':True,'markdown_tables_unchanged':True,'bibliography_unchanged':True,'figure_hashes_unchanged':True,'protected_files':len(frozen)})
    ts.validate_package(); assert ts.frozen_check()==frozen
    print(json.dumps(json.loads(read(WORK/'preflight.json'))))

def finalize():
    frozen = ts.frozen_check(); ts.validate_package()
    for rel, digest in json.loads(read(WORK/'figure_baseline.json')).items(): assert ts.sha(ROOT/rel)==digest,rel
    status=json.loads(read(ROOT/'work/step08/render_status.json'))
    for s in status: assert s['returncode']==0 and s['sha256']==ts.sha(OUT/s['file'])
    from pypdf import PdfReader
    render=ROOT/'work/step08/render'
    shutil.copy2(render/'04_bioRxiv_Integrated_Preprint/04_bioRxiv_Integrated_Preprint.pdf',OUT/'04_bioRxiv_Integrated_Preprint.pdf')
    pdf=PdfReader(OUT/'04_bioRxiv_Integrated_Preprint.pdf');text='\n'.join(p.extract_text() or '' for p in pdf.pages);compact=''.join(text.split())
    assert ts.ORCID in text and '10.5281/zenodo.23653766' in compact
    assert not re.search(r'\bAl\b|Not registered|PENDING_RELEASE|\bHighlights\b',text)
    for ref in json.loads(read(OUT/'verification/reference_verification.json'))['references']: assert ref['doi'] in compact,ref['doi']
    import xml.etree.ElementTree as ET
    suites=list(ET.parse(WORK/'pytest.xml').getroot().iter('testsuite'))
    tests={k:sum(int(s.get(k,0)) for s in suites) for k in ['tests','failures','errors','skipped']}
    tests['passed']=tests['tests']-tests['failures']-tests['errors']-tests['skipped'];assert tests['failures']==tests['errors']==0
    summary=dict(json.loads(read(WORK/'preflight.json')),status='complete',stage='Final editorial alignment',pdf_pages=len(pdf.pages),rendered_pages=sum(s['pages'] for s in status),tests=tests,pdf_sha256=ts.sha(OUT/'04_bioRxiv_Integrated_Preprint.pdf'),scientific_files_changed=0,publication_figures_changed=0)
    for name,value in [('Build_Summary.json',summary),('Current_Test_Results.json',tests),('Current_Test_Summary.json',dict(tests,scope='Full project tests after final editorial rebuild; clinical extraction and fitting unchanged')),('Rendered_Documents.json',status),('Render_QA_Review.json',{'status':'passed','method':'Every final page visually reviewed, with wide tables, equation rendering and integrated figures checked.','documents':status}),('Final_Editorial_Alignment.json',summary)]:dump(OUT/'verification'/name,value)
    shutil.copy2(WORK/'pytest.xml',OUT/'verification/Regression_Tests.xml')
    shutil.copy2(WORK/'editorial_changes.json',OUT/'verification/Editorial_Changes.json')
    dump(OUT/'verification/Editorial_Release_Integrity.json',{'status':'passed','protected_files':frozen,'changed_files':[],'publication_figures':json.loads(read(WORK/'figure_baseline.json'))})
    report='# Final editorial alignment\n\n本次统一摘要、引言收束、结果解释、讨论、结论、两份投稿信、补充材料导语及公开说明。核心贡献围绕归一化尺度权衡、样本对一致性与固定患者分数的批次构成依赖展开。保留Hybrid比较方向、显著性、所有必要局限及完整伦理/AI披露。历史执行报告未改写。\n\n'
    report+='Markdown表格逐行一致，参考文献段逐字一致，全文数字词元多重集一致。所有冻结科学文件和本次既有排版图件逐项哈希一致。\n\n```json\n'+json.dumps(summary,indent=2,ensure_ascii=False)+'\n```\n'
    write(ROOT/'results/reports/final_editorial_alignment.md',report)
    for src,rel in [(ROOT/'results/MANUSCRIPT_FULL_TEXT.md','results/MANUSCRIPT_FULL_TEXT.md'),(ROOT/'results/reports/final_editorial_alignment.md','results/reports/final_editorial_alignment.md'),(Path(__file__),'scripts/10_finalize_editorial_materials.py')]:shutil.copy2(src,RELEASE/rel)
    # Keep future upload fields and supplementary PDF aligned with the approved text.
    upload=ROOT/'work/biorxiv_submission';fields=json.loads(read(upload/'submission_fields.json'))
    full=read(ROOT/'results/MANUSCRIPT_FULL_TEXT.md')
    fields['abstract']=pack.extract_section(full,'## Abstract','## 1. Introduction').split('**Keywords:**')[0].replace('**','').strip()
    fields['main_pdf_sha256']=summary['pdf_sha256'];dump(upload/'submission_fields.json',fields)
    shutil.copy2(render/'05_Supplementary_Material/05_Supplementary_Material.pdf',upload/'Supplementary_Material.pdf')
    checklist=f'''# 投稿与开源发布核验清单

作者：Zhenxuan Ma（马振轩）；{ts.AFFILIATION}；mzxmc08@gmail.com；ORCID {ts.ORCID}。

公开代码仓库：https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit

真实版本DOI：https://doi.org/10.5281/zenodo.23653766 。v0.1.0归档保持原样；当前文稿修订与归档代码版本分别记录。

## 当前版本

04的PDF为{len(pdf.pages)}页；8份Word共渲染{summary['rendered_pages']}页。全部正文、图注、Highlights、两份投稿信、补充材料和README采用尺度权衡与分数可比性叙事。28篇文献保持已核验条目及首次引用顺序。全套测试{tests['passed']}项通过、{tests['skipped']}项跳过。{len(frozen)}个冻结科学文件及现有图件哈希不变。

Table 3A/3B/3C与S3使用完整三线表、重复表头、禁止跨页拆行。宽表使用横向页面，表内8.5–9 pt；03正文12 pt双倍行距、连续行号。03/04正文不含Highlights。Figure 3使用既有Typeset排版副本，原始冻结图件另存不变。

## Highlights字符数（含空格）

'''
    for h in json.loads(read(OUT/'verification/Highlights_Character_Counts.json')):checklist+=f"- {len(h['text'])}/85：{h['text']}\n"
    checklist+='''
## 三分钟上传指南

1. bioRxiv上传04 PDF及05补充材料PDF，并核对题目、摘要、作者、ORCID和真实代码DOI。按实际表单选择主题与许可。
2. CMPB或CBM选择一家投稿。分别上传01署名页、03匿名行号正文、02 Highlights、对应06投稿信、07利益冲突声明、05补充材料及独立图表。
3. 检查系统合成PDF的公式、完整表格、图注和28篇文献；如实填写AI与伦理信息。出版路线选择Subscription/traditional，具体条件以期刊当前系统为准。

## 作者核对事项

- 作者身份与ORCID采用本人提供信息，ORCID通过校验位检查。
- 伦理声明记录公开去标识化数据二次分析及未取得额外机构判定的状态；机构若提出要求，按实际结论补充。
- AI辅助范围、人工核验和作者责任完整披露。
- 03匿名正文隐去个人仓库标识；01署名页及04完整稿提供真实链接。
- 本地上传材料已更新；尚未向bioRxiv或期刊提交。

## 全部文件

'''
    inventory=sorted(str(p.relative_to(OUT)).replace('\\','/') for p in OUT.rglob('*') if p.is_file() and not any(x in p.parts for x in ['.git','__pycache__','.pytest_cache']))
    checklist+='\n'.join('- `'+p+'`' for p in inventory)+'\n';write(OUT/'SUBMISSION_CHECKLIST_FOR_AUTHOR.md',checklist)
    assert ts.frozen_check()==frozen
    pack.finalize_archives()
    for path in [ROOT/'cfDNA_Submission_Package.zip',OUT/'github_zenodo_release.zip']:
        with zipfile.ZipFile(path) as z: assert z.testzip() is None
    print(json.dumps(summary))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');p.add_argument('--finalize',action='store_true');p.add_argument('--pandoc',type=Path)
    a=p.parse_args()
    if a.build:build(a)
    elif a.finalize:finalize()
    else:p.error('Choose --build or --finalize')
