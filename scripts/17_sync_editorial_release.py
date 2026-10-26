"""Prepare publication-only patch releases from existing frozen results.

Network publication is deliberately separate from this local assembly command.
The historical v0.2.0 packager remains available for reproducing that version.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import re
import shutil
from zipfile import ZipFile
from lxml import etree
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'submission_package'
RELEASE = OUT / 'github_zenodo_release'


def read(path):
    return path.read_text(encoding='utf-8')


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    assert sha(source) == sha(target)


def check():
    frozen = json.loads(read(ROOT / 'results/step10_final/reports/frozen_inputs_sha256.json'))
    for name, expected in frozen.items():
        assert sha(ROOT / name) == expected, name
    text = read(ROOT / 'results/MANUSCRIPT_FULL_TEXT.md')
    assert r'\overline{R}_{s}' in text
    assert r'\operatorname{sd}_{i}(R_{si})' in text
    assert 'Doebley et al. further link nucleosome protection profiles' in text
    ns = {'m': 'http://schemas.openxmlformats.org/officeDocument/2006/math'}
    for name in ['03_Manuscript_DoubleSpaced_LineNumbered', '04_bioRxiv_Integrated_Preprint']:
        with ZipFile(OUT / (name + '.docx')) as z:
            doc = etree.fromstring(z.read('word/document.xml'))
        equations = [m for m in doc.xpath('//m:oMath', namespaces=ns) if 'sd' in ''.join(m.itertext())]
        assert len(equations) == 1
        assert equations[0].xpath('.//m:bar[m:barPr/m:pos[@m:val="top"]]', namespaces=ns)
        assert '(' in ''.join(equations[0].itertext()) and ')' in ''.join(equations[0].itertext())
    return len(frozen)


def sync(pack=False):
    protected = check()
    state = json.loads(read(ROOT / 'config/release.json'))
    version = state['version']
    tag = 'v' + version
    url = state['github_release_url']
    doi = state.get('zenodo_doi')
    concept = state['concept_zenodo_doi']
    published = state.get('status') == 'published'
    citation = json.loads(read(RELEASE / 'CITATION.cff'))
    citation.update(version=version, message='Please cite the release tag and its Zenodo archive.', url=url)
    if doi:
        citation['doi'] = doi
    else:
        citation.pop('doi', None)
    dump(RELEASE / 'CITATION.cff', citation)
    copy(RELEASE / 'CITATION.cff', ROOT / 'CITATION.cff')
    zen = json.loads(read(RELEASE / '.zenodo.json'))
    zen.update(version=version, description=f'Software {tag}: editorial corrections to the sample-mean overbar and standard-deviation notation, explicit Doebley attribution, and organized current Word/PDF publication assets. Scientific results, predictions and figures remain unchanged from v0.2.0. Original software is MIT licensed; third-party data retain their original terms.')
    zen['related_identifiers'] = [{'identifier': state['previous_zenodo_doi'], 'relation': 'isNewVersionOf', 'resource_type': 'software'}]
    dump(RELEASE / '.zenodo.json', zen)
    example = RELEASE / 'config/submission_author.example.json'
    if example.exists():
        data = json.loads(read(example))
        data.update(software_version=version, zenodo_doi=doi, zenodo_status=state['status'])
        dump(example, data)
    for folder in ['src', 'scripts', 'tests', 'config']:
        for p in (ROOT / folder).rglob('*'):
            if not p.is_file() or any(part in p.parts for part in ['__pycache__', '.pytest_cache', '.git']):
                continue
            if p.name == 'submission_author.json' or p.suffix not in ['.py', '.json', '.yaml', '.yml', '.md', '.toml']:
                continue
            copy(p, RELEASE / p.relative_to(ROOT))
    for p in [ROOT / 'results/MANUSCRIPT_FULL_TEXT.md', ROOT / 'results/reports/MANUSCRIPT_MASTER_DRAFT.md',
              ROOT / 'results/reports/formula_citation_polish.md']:
        copy(p, RELEASE / p.relative_to(ROOT))
    for p in (ROOT / 'results/step10_final/reports').iterdir():
        if p.is_file():
            copy(p, RELEASE / p.relative_to(ROOT))
    header = read(RELEASE / 'README.md')
    header = re.sub(r'^# cfDNA Scale Tradeoff Audit v[^\n]+', '# cfDNA Scale Tradeoff Audit ' + tag, header)
    archive_line = f'Version-specific DOI: https://doi.org/{doi}.' if doi else f'Archive collection DOI: https://doi.org/{concept}. The exact version is identified by the {tag} Git tag; its version-specific DOI is assigned by Zenodo after publication.'
    section = f'## Version and archival status\n\nCurrent patch release: **{tag}**, {url}.\n{archive_line}\n\nThis patch corrects formula rendering and reference attribution, and includes current publication documents in `publication/`. Scientific inputs and numerical results are unchanged. Previous v0.2.0 remains immutable at https://doi.org/{state["previous_zenodo_doi"]}.\n\n'
    header = re.sub(r'## Version and archival status\n.*?(?=## Package contents)', section, header, flags=re.S)
    if '- `publication/`:' not in header:
        header = header.replace('- `results/MANUSCRIPT_FULL_TEXT.md`: current manuscript and one AI-use disclosure.', '- `results/MANUSCRIPT_FULL_TEXT.md`: current manuscript and one AI-use disclosure.\n- `publication/`: current editable Word manuscripts, integrated PDFs and file guide.\n- `scripts/17_sync_editorial_release.py`: publication-only patch synchronization; no classifier fitting.')
    header = header.replace('script 15 synchronizes the package and seals manifests.', 'script 15 seals the historical v0.2.0 package. For later editorial patches, use script 17 to synchronize and pack existing results.')
    (RELEASE / 'README.md').write_text(header, encoding='utf-8')
    copy(RELEASE / 'README.md', ROOT / 'README.md')
    guide = f'''# 当前定稿与发布入口

当前编辑修订版本：{tag}。
GitHub：{url}
{archive_line}

正文与投稿文件位于 `publication/`（公开仓库）及本地 `submission_package/`。
期刊正文：03_Manuscript_DoubleSpaced_LineNumbered.docx；预印本图文稿：04_bioRxiv_Integrated_Preprint.docx / .pdf。
01 为标题页，02 为 Highlights，05 为补充材料，06/06b 为两个备选期刊的投稿信，07 为利益冲突声明。实际投稿仅选择一个期刊。

当前图表：results/step10_final/figures 与 tables；历史阶段报告保留原有时间语境。
本版本已包含均值上横线、标准差括号和 Doebley 作者归属修订。科学结果及旧版本归档保持原样。
本轮仅同步 GitHub/Zenodo，未上传 bioRxiv 修订稿。
'''
    (OUT / 'START_HERE.md').write_text(guide, encoding='utf-8')
    for p in OUT.iterdir():
        if p.is_file() and p.suffix in ['.docx', '.pdf', '.md'] and (p.name[:2].isdigit() or p.name in ['START_HERE.md', 'Manuscript_Main_Blind_and_Full.md', 'Manuscript_Main_Blinded.md']):
            copy(p, RELEASE / 'publication' / p.name)
    for name in ['04_bioRxiv_Integrated_Preprint', '05_Supplementary_Material']:
        p = ROOT / 'work/step08/render' / name / (name + '.pdf')
        copy(p, OUT / p.name)
        copy(p, RELEASE / 'publication' / p.name)
    rendered = json.loads(read(ROOT / 'work/step08/render_status.json'))
    for item in rendered:
        assert item['returncode'] == 0 and item['sha256'] == sha(OUT / item['file'])
    old_summary = json.loads(read(OUT / 'verification/Build_Summary.json'))
    old_summary.update(state, rendered_documents=rendered,
                       pdf_pages=len(PdfReader(OUT / '04_bioRxiv_Integrated_Preprint.pdf').pages),
                       pdf_sha256=sha(OUT / '04_bioRxiv_Integrated_Preprint.pdf'),
                       public_release_created=published)
    dump(OUT / 'verification/Build_Summary.json', old_summary)
    dump(OUT / 'verification/Rendered_Documents.json', rendered)
    status = {'version': version, 'status': state['status'], 'zenodo_doi': doi,
              'concept_zenodo_doi': concept, 'protected_files': protected, 'classifier_fits': 0,
              'feature_extraction_runs': 0, 'publication_documents': 8, 'biorxiv_revision_uploaded': False}
    status.update(pdf_pages=old_summary['pdf_pages'], pdf_sha256=old_summary['pdf_sha256'])
    dump(OUT / 'verification/Editorial_Release_Status.json', status)
    copy(OUT / 'verification/Editorial_Release_Status.json', RELEASE / 'results/reports/Editorial_Release_Status.json')
    receipt = OUT / 'verification/Latest_Publication_Receipt.json'
    if receipt.exists():
        copy(receipt, RELEASE / 'results/reports/Latest_Publication_Receipt.json')
    (OUT / 'Manuscript_Main_Blinded.md').write_bytes((OUT / '03_Manuscript_DoubleSpaced_LineNumbered.md').read_bytes())
    copy(OUT / 'Manuscript_Main_Blinded.md', RELEASE / 'publication/Manuscript_Main_Blinded.md')
    if pack:
        spec = importlib.util.spec_from_file_location('publication', ROOT / 'scripts/12_build_major_revision_package.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.validate_package()
        archives = module.pack.finalize_archives()
        status['archives'] = archives
    print(json.dumps(status, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pack', action='store_true')
    sync(parser.parse_args().pack)
