"""Build the locally revised submission files; never submit or refit clinical data.

Run 11_major_revision_analysis.py in the analysis environment first, then compile
src.revision_manuscript there. This authoring command uses the bundled document
runtime and the explicit --pandoc executable. Historical v1 remains archived.
"""
from pathlib import Path
import argparse,copy,csv,hashlib,importlib.util,json,re,shutil,zipfile
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt,Inches

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'submission_package'
NEW=ROOT/('results/step10_final' if (ROOT/'results/step10_final/reports/publication_manifest.json').exists() else 'results/step09_major');WORK=ROOT/'work/step09_major'
spec=importlib.util.spec_from_file_location('pack08',ROOT/'scripts/08_build_submission_package.py')
pack=importlib.util.module_from_spec(spec);spec.loader.exec_module(pack)
pack.presentation_tables=lambda text:text
pack.QA=WORK;WORK.mkdir(parents=True,exist_ok=True)
FIGURES=['Figure1_Batch_Vulnerability_Audit','Figure2_Primary_Ranking_Threshold_Transport',
 'Figure3_Companion_Dependence_and_Invariance','Figure4_Genomic_Context_and_Ranks',
 'FigureS5_Score_Comparability_and_Calibration','FigureS6_Secondary_Classifiers']
HIGHLIGHTS=[
 'Normalization-scale discrimination depends on the held-out cfDNA study.',
 'Source-locked thresholds shift despite companion-independent preprocessing.',
 'Real-patient anchor resampling separates composition and identity effects.',
 'Study-held-out calibration and genomic block inference expose distinct risks.']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return Path(p).read_text(encoding='utf-8')
def dump(p,data):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf-8')

def table_layout(doc):
    """Complete compact three-line tables, with deliberate landscape sections."""
    def section_boundary(element,landscape=False,before=False):
        p=OxmlElement('w:p');pr=OxmlElement('w:pPr');sect=copy.deepcopy(doc.sections[-1]._sectPr)
        size=sect.find(qn('w:pgSz'))
        size.set(qn('w:w'),'15840' if landscape else '12240');size.set(qn('w:h'),'12240' if landscape else '15840')
        size.set(qn('w:orient'),'landscape' if landscape else 'portrait')
        typ=sect.find(qn('w:type'))
        if typ is None:typ=OxmlElement('w:type');sect.insert(0,typ)
        typ.set(qn('w:val'),'nextPage');pr.append(sect);p.append(pr)
        element.addprevious(p) if before else element.addnext(p)
    for table in doc.tables:
        cap=table._tbl.getprevious();caption=''.join(cap.itertext()) if cap is not None else ''
        n=len(table.columns);wide=False
        width=648 if wide else 468
        if n==8:weights=[46,54,42,66,42,63,108,35]
        elif n==6 and 'Table 2.' in caption:weights=[58,80,120,145,145,100]
        elif n==5 and 'Table 3B.' in caption:weights=[85,110,120,75,180]
        elif n==4 and 'Table 2.' in caption:weights=[95,125,125,123]
        else:weights=[105]+[85]*(n-1)
        widths=[width*w/sum(weights) for w in weights]
        for c,w in zip(table.columns,widths):c.width=Pt(w)
        for old in list(table._tbl.tblPr.findall(qn('w:tblCellMar'))):table._tbl.tblPr.remove(old)
        margins=OxmlElement('w:tblCellMar')
        for side,value in [('top',35),('bottom',35),('left',42),('right',42)]:
            e=OxmlElement('w:'+side);e.set(qn('w:w'),str(value));e.set(qn('w:type'),'dxa');margins.append(e)
        table._tbl.tblPr.append(margins)
        for ri,row in enumerate(table.rows):
            pr=row._tr.get_or_add_trPr()
            if pr.find(qn('w:cantSplit')) is None:pr.append(OxmlElement('w:cantSplit'))
            if ri==0 and pr.find(qn('w:tblHeader')) is None:pr.append(OxmlElement('w:tblHeader'))
            for i,cell in enumerate(row.cells):
                cell.width=Pt(widths[i])
                for p in cell.paragraphs:
                    p.paragraph_format.line_spacing=1;p.paragraph_format.space_after=Pt(2)
                    p.paragraph_format.space_before=Pt(0);p.paragraph_format.keep_with_next=ri<len(table.rows)-1
                    for r in p.runs:pack.set_font(r,9);r.bold=ri==0
        if wide:
            section_boundary(cap,before=True);section_boundary(table._tbl,landscape=True)
    for i,sect in enumerate(doc._element.body.iter(qn('w:sectPr'))):
        line=sect.find(qn('w:lnNumType'))
        if line is not None and i:line.attrib.pop(qn('w:start'),None)
    return doc

def insert_figures(doc):
    from docx.text.paragraph import Paragraph
    mapping={'Figure 1.':FIGURES[0],'Figure 2.':FIGURES[1],'Figure 3.':FIGURES[2],
      'Figure 4.':FIGURES[3],'Supplementary Figure S5.':FIGURES[4],'Supplementary Figure S6.':FIGURES[5]}
    seen=[]
    for p in list(doc.paragraphs):
        key=next((key for key in mapping if p.text.startswith(key)),None)
        if key is None:continue
        e=OxmlElement('w:p');p._p.addprevious(e);imagep=Paragraph(e,p._parent)
        imagep.paragraph_format.page_break_before=not (doc.core_properties.subject=='Supplementary first figure' and not seen);imagep.paragraph_format.keep_with_next=True
        imagep.paragraph_format.line_spacing=1;imagep.paragraph_format.space_after=Pt(4)
        imagep.add_run().add_picture(str(NEW/'figures'/(mapping[key]+'.png')),width=Inches(6.5))
        p.paragraph_format.line_spacing=1.05;p.paragraph_format.space_after=Pt(4)
        p.paragraph_format.keep_together=True
        for r in p.runs:pack.set_font(r,10)
        seen.append(mapping[key])
    return seen

def refresh_figures_only():
    """Replace embedded rasters and their aspect ratios, leaving all text fixed."""
    from PIL import Image
    for name,names in [('04_bioRxiv_Integrated_Preprint',FIGURES),('05_Supplementary_Material',FIGURES[4:])]:
        p=OUT/(name+'.docx');doc=Document(p)
        assert len(doc.inline_shapes)==len(names)
        for shape,stem in zip(doc.inline_shapes,names):
            image=NEW/'figures'/(stem+'.png')
            rid=shape._inline.graphic.graphicData.pic.blipFill.blip.embed
            doc.part.related_parts[rid]._blob=image.read_bytes()
            with Image.open(image) as im:width,height=im.size
            shape.width=Inches(6.5);shape.height=Inches(6.5*height/width)
        doc.save(p)
    dump(OUT/'verification/Major_Revision_Document_Validation.json',validate_package())
    for p in (NEW/'figures').glob('*'):
        if p.suffix in ('.png','.pdf'):shutil.copy2(p,OUT/'figures_and_tables'/NEW.name/'figures'/p.name)

def author_text(full):
    return full.split('## Data and code availability\n')[1].split('## Figure legends\n')[0]

def validate_package():
    result={}
    for p in sorted(OUT.glob('*.docx')):
        doc=Document(p);text='\n'.join(x.text for x in doc.paragraphs)
        assert 'Not registered' not in text and 'PENDING_RELEASE' not in text,p
        assert 'Proposed_LSN' not in text and 'InFold_ComBat' not in text,p
        if p.name.startswith(('03_','04_')):
            assert not any(x.text.strip()=='Highlights' for x in doc.paragraphs)
            assert len(doc.tables)==5,(p,len(doc.tables))
            for label in ['Table 1.','Table 2.','Table 3A.','Table 3B.','Table 3C.']:assert label in text,label
            assert text.count('Declaration of generative AI')==1
            refs=text.split('References')[-1]
            assert len([line for line in refs.splitlines() if 'https://doi.org/' in line])==29
            assert '0/1; exact 95% CI [0.000, 0.975]' in '\n'.join(c.text for t in doc.tables for r in t.rows for c in r.cells)
        for t in doc.tables:
            assert t.rows[0]._tr.trPr.find(qn('w:tblHeader')) is not None
            assert all(r._tr.trPr.find(qn('w:cantSplit')) is not None for r in t.rows)
        if p.name.startswith('03_'):
            lines=list(doc._element.iter(qn('w:lnNumType')))
            assert lines and all(x.get(qn('w:restart'))=='continuous' for x in lines)
            assert 'Zhenxuan Ma' not in text and 'mzxmc08' not in text
        result[p.name]={'status':'passed','tables':len(doc.tables),'images':len(doc.inline_shapes),'sha256':sha(p)}
    assert len(result)==8,len(result)
    return result

def build(pandoc):
    full=read(ROOT/'results/MANUSCRIPT_FULL_TEXT.md');author=json.loads(read(ROOT/'config/submission_author.json'))
    title=full.splitlines()[0][2:];abstract=full.split('## Abstract\n')[1].split('## 1. Introduction')[0]
    body=full.split('## 1. Introduction')[1].split('## Data and code availability')[0]
    no_tables=re.sub(r'(?m)^\|.*\n?','',body)
    counts={'abstract_words':pack.word_count(abstract.split('**Keywords:**')[0]),
      'main_body_words_excluding_tables':pack.word_count(no_tables),'main_body_including_tables_words':pack.word_count(body),
      'convention':'Introduction through Conclusions including headings/equations; keywords excluded from abstract.'}
    dump(OUT/'verification/Word_Counts.json',counts)
    keywords=abstract.split('**Keywords:**')[1].strip()
    first=(f'# {title}\n\n{author["name"]}\n\n{author["affiliation"]}\n\nCorresponding author: {author["name"]}; {author["email"]}\n\n'
      f'ORCID: https://orcid.org/{author["orcid"]}\n\nAbstract words: {counts["abstract_words"]}\n\n'
      f'Main text words excluding tables: {counts["main_body_words_excluding_tables"]}\n\n'
      f'Main text words including tables: {counts["main_body_including_tables_words"]}\n\nKeywords: {keywords}\n\n'
      '## Data and code availability\n'+author_text(full))
    assert len(HIGHLIGHTS) in range(3,6) and all(len(x)<=85 for x in HIGHLIGHTS)
    texts={'01_Title_Page':first,'02_Highlights':'# Highlights\n\n'+'\n\n'.join('- '+x for x in HIGHLIGHTS),
      '03_Manuscript_DoubleSpaced_LineNumbered':'# '+title+'\n\n## Abstract\n'+full.split('## Abstract\n')[1],
      '04_bioRxiv_Integrated_Preprint':full,'05_Supplementary_Material':read(OUT/'05_Supplementary_Material.md')}
    # Blinded review copy carries factual declarations without author identifying URLs.
    blind=texts['03_Manuscript_DoubleSpaced_LineNumbered']
    blind=re.sub(r'## Data and code availability\n.*?(?=## Figure legends)',
      '## Declarations\n\nAuthor-identifying availability, authorship and compliance statements are supplied on the separate title page.\n\n'
      '## Declaration of generative AI and AI-assisted technologies in the writing process\n\n'
      +full.split('## Declaration of generative AI and AI-assisted technologies in the writing process\n')[1].split('## Figure legends')[0],blind,flags=re.S)
    texts['03_Manuscript_DoubleSpaced_LineNumbered']=blind
    delta=list(csv.DictReader(open(NEW/'tables/Primary_Hybrid_Global_Delta.csv',encoding='utf-8')))
    for name,journal in [('06_Cover_Letter','Computer Methods and Programs in Biomedicine'),('06b_Cover_Letter_Alternative_CBM','Computers in Biology and Medicine')]:
        texts[name]=(f'# Cover letter\n\nDear Editor-in-Chief,\n\nPlease consider my manuscript, “{title}”, for {journal}.\n\n'
         'This study distinguishes three properties that are often combined in public cfDNA evaluation: the genomic reference scale of normalization, dependence on companion patients, and transfer of a source-locked operating threshold. The primary evidence is a bidirectional study-held-out comparison of 548 participants; all analyses use an existing 585-participant public cohort.\n\n'
         f'Paired Hybrid−Global AUROC differences were {float(delta[0]["delta_auroc"]):.3f} in Cristiano and {float(delta[1]["delta_auroc"]):+.3f} in Jiang. '
         'The manuscript therefore presents a target-dependent scale trade-off. A mathematical invariance argument, repeated real-patient anchors, source-study-held-out calibration and chromosome-block inference establish distinct tests of preprocessing and operating-point transport. At fixed case fraction 0.50, the evaluated EB-TargetMap rule changed anchor decisions in 0.140 and 0.126 of repeated evaluations, whereas the sample-wise rules satisfy companion independence by construction.\n\n'
         'The accompanying reproducible code includes automated leakage and invariance tests. No new human recruitment or specimen collection was performed. Study source is treated as a proxy for technical and population differences.\n\n'
         'The work is original and is not under consideration at another journal. A previous version was submitted to bioRxiv as BIORXIV/2026/761726 Version 1 and is awaiting screening; this local revision has not been uploaded. I request the traditional subscription publishing model. I take responsibility for reviewing the manuscript and for its content at submission.\n\n'
         f'Sincerely,\n\n{author["name"]}\n\n{author["affiliation"]}\n\n{author["email"]}\n')
    texts['07_Declaration_of_Interest']='# Declaration of interest\n\n'+title+'\n\n'+author['name']+'\n\n'+author['competing_interests']+'\n'
    validation={}
    for name,text in texts.items():
        (OUT/(name+'.md')).write_text(text,encoding='utf-8')
        doc=pack.pandoc_doc(text,OUT/(name+'.docx'),pandoc,author,blind=name.startswith('03_'),double=name.startswith('03_'),numbers=name.startswith('03_'))
        table_layout(doc)
        aliases={'Global_Zscore':'Global','Hybrid_MultiScale':'Hybrid','Local_LSN_ICC':'Local LSN + ICC',
            'Local_LSN':'Local LSN','Raw_Ratio':'Raw','Tuned_Local_LSN_ICC':'Tuned Local'}
        for table in doc.tables:
            for row in table.rows:
                for cell in row.cells:
                    if cell.text in aliases:
                        cell.paragraphs[0].runs[0].text=aliases[cell.text]
        if name.startswith(('04_','05_')):
            inrefs=False
            for p in list(doc.paragraphs):
                if p.text=='Figure legends':p._p.getparent().remove(p._p);continue
                if p.text=='References':inrefs=True
                p.paragraph_format.line_spacing=1 if inrefs else 1.1
                p.paragraph_format.space_after=Pt(4)
            if name.startswith('05_'):doc.core_properties.subject='Supplementary first figure'
        if name.startswith(('04_','05_')):validation[name]=insert_figures(doc)
        for p in doc.paragraphs:
            if re.match(r'^(Table \d|Supplementary Table)',p.text):p.paragraph_format.keep_with_next=True
        doc.save(OUT/(name+'.docx'))
    refs=json.loads(read(OUT/'verification/reference_verification.json'))
    if len(refs['references'])==28:
        record=json.loads(read(NEW/'reports/PERMDISP_Crossref.json'))['message']
        refs['references'].append({'reference':29,'doi':record['DOI'],'year':2006,'volume':record['volume'],'issue':record['issue'],
            'pages_or_article':record['page'],'title':record['title'][0],'crossref_response':'results/step09_major/reports/PERMDISP_Crossref.json'})
    refs['numbering_note']='Metadata source IDs retain the original order; see Revision_References.csv for current manuscript numbers.'
    dump(OUT/'verification/reference_verification.json',refs)
    dump(OUT/'verification/Major_Revision_Document_Validation.json',validate_package())
    dump(OUT/'verification/Major_Revision_Highlights.json',[{'text':x,'characters':len(x)} for x in HIGHLIGHTS])
    dump(OUT/'verification/Major_Revision_Images.json',validation)
    for sub in ['figures','tables','reports']:
        dest=OUT/'figures_and_tables'/NEW.name/sub;dest.mkdir(parents=True,exist_ok=True)
        for p in (NEW/sub).glob('*'):
            if p.is_file() and p.suffix in ('.csv','.pdf','.png','.json','.md'):shutil.copy2(p,dest/p.name)
    (OUT/'figures_and_tables/README_CURRENT_REVISION.md').write_text(
      f'# Current revision\n\nThe current Word files embed figures from {NEW.name}/figures. Publication tables are in {NEW.name}/tables. Historical figure/table copies elsewhere in this folder belong to the submitted v1 and are retained for traceability.\n',encoding='utf-8')
    print(json.dumps({'documents':len(texts),'counts':counts,'highlights':[len(x) for x in HIGHLIGHTS]},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pandoc',type=Path);ap.add_argument('--refresh-figures-only',action='store_true');args=ap.parse_args()
    if args.refresh_figures_only:refresh_figures_only()
    else:
        if args.pandoc is None:ap.error('--pandoc is required for a full build')
        build(args.pandoc)
