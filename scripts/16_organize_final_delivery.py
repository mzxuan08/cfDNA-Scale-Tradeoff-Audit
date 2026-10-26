"""Assemble current publication assets without changing scientific outputs."""
from pathlib import Path
import hashlib
import json
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'submission_package'
DEST = ROOT / 'FINAL_DELIVERY'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy(source, relative):
    target = DEST / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    assert digest(source) == digest(target), relative


def main():
    frozen = json.loads((ROOT / 'results/step10_final/reports/frozen_inputs_sha256.json').read_text(encoding='utf-8'))
    for name, expected in frozen.items():
        assert digest(ROOT / name) == expected, name
    release_state = json.loads((ROOT / 'config/release.json').read_text(encoding='utf-8'))
    DEST.mkdir(exist_ok=True)
    for source in sorted(PACKAGE.glob('*')):
        if source.is_file() and source.suffix in ['.docx', '.pdf', '.md'] and source.name[:2].isdigit():
            copy(source, Path('01_Manuscripts') / source.name)
    copy(ROOT / 'results/MANUSCRIPT_FULL_TEXT.md', Path('01_Manuscripts/MANUSCRIPT_FULL_TEXT.md'))
    for folder, target, suffixes in [
        ('results/step10_final/figures', '02_Figures', {'.pdf', '.png'}),
        ('results/step10_final/tables', '03_Tables/current_revision', {'.csv'}),
        ('results/tables', '03_Tables/original_analysis_outputs', {'.csv'}),
        ('submission_package/verification', '04_Reports_and_Checks/verification', None),
        ('results/step10_final/reports', '04_Reports_and_Checks/final_seal', None),
        ('results/reports', '04_Reports_and_Checks/research_history', {'.md', '.json'}),
    ]:
        for source in sorted((ROOT / folder).iterdir()):
            if source.is_file() and (suffixes is None or source.suffix in suffixes):
                copy(source, Path(target) / source.name)
    for name in ['FINAL_5_POINTS_CHECK.md', 'SUBMISSION_CHECKLIST_FOR_AUTHOR.md', 'START_HERE.md']:
        copy(PACKAGE / name, Path('04_Reports_and_Checks') / name)
    copy(PACKAGE / 'github_zenodo_release.zip', Path('05_Code/current_working_copy.zip'))
    copy(ROOT / 'cfDNA-Scale-Tradeoff-Audit-v0.2.0.zip', Path('05_Code/published_v0.2.0_immutable.zip'))
    receipt = PACKAGE / 'verification/Latest_Publication_Receipt.json'
    copy(receipt if receipt.exists() else PACKAGE / 'verification/V020_Publication_Receipt.json', Path('05_Code/Publication_Receipt.json'))
    patch_archive = ROOT / ('cfDNA-Scale-Tradeoff-Audit-v' + release_state['version'] + '.zip')
    if release_state['version'] != '0.2.0' and patch_archive.exists():
        copy(patch_archive, Path('05_Code') / ('published_v' + release_state['version'] + '_immutable.zip'))
    note = '''# 最终交付文件导览

## 使用顺序

1. **01_Manuscripts**：当前完整英文正文、8 份 Word 投稿文件、图文整合 PDF 和补充材料 PDF。生物预印本使用 `04_bioRxiv_Integrated_Preprint.pdf`；期刊正文使用 `03_Manuscript_DoubleSpaced_LineNumbered.docx`，另交 Title Page、Highlights、补充材料和利益冲突声明。
2. **02_Figures**：当前 Figure 1–4、Supplementary Figure S5/S6，各含矢量 PDF 和 PNG。
3. **03_Tables/current_revision**：当前修订结果及分期核对 CSV。`original_analysis_outputs` 保留此前分析的完整数值输出，用于追溯；不能把不同验证范围的指标混作同一结果。
4. **04_Reports_and_Checks**：当前核验资料及研究历史。优先阅读 `START_HERE.md`、`FINAL_5_POINTS_CHECK.md` 和 `verification/Formula_Citation_Polish.md`。`research_history` 中的阶段报告描述各自当时的状态，当前状态以此导览及发布回执为准。
5. **05_Code**：当前本地代码目录压缩包、不可变的已发布 v0.2.0 代码包及发布回执。两者区别见下文。

## 版本与投稿状态

- GitHub v0.2.0：https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit/releases/tag/v0.2.0
- 已发布版本 DOI：https://doi.org/10.5281/zenodo.23659205
- `published_v0.2.0_immutable.zip` 对应已发布软件快照，不改写既有归档。
- `current_working_copy.zip` 包含归档后完成的公式排版、引文措辞修订与本次整理工具；它不是新发布版本，不代表远端主分支已同步这些修订。
- 当前 Word/PDF 已包含样本均值上横线、标准差参数括号，以及 Doebley 等人的明确作者归属。
- 两封 Cover Letter 对应两个备选期刊；实际投稿时仅选择一个，不同时投稿。
- 本地最新修订稿尚未另行上传 bioRxiv；本次整理不执行网站投稿或创建新归档版本。

## 完整性

本次仅整理与复制文件，不训练分类器、不重新提取特征、不修改统计 CSV 或图件。1,479 个受保护文件按既有清单校验。`DELIVERY_MANIFEST_SHA256.json` 记录本交付目录的文件校验值，`DELIVERY_STATUS.json` 记录文件数量与检查结果。

原始片段流、患者特征矩阵及模型缓存继续保留在研究工作区，不装入投稿交付包；代码与表格的公开范围沿用此前封版配置。
'''
    if release_state['version'] != '0.2.0':
        version = release_state['version']
        note = note.replace('GitHub v0.2.0：https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit/releases/tag/v0.2.0', 'GitHub v' + version + '：' + release_state['github_release_url'])
        note = note.replace('已发布版本 DOI：https://doi.org/10.5281/zenodo.23659205', '当前版本 DOI：https://doi.org/' + (release_state.get('zenodo_doi') or release_state['concept_zenodo_doi']))
        note = note.replace('`current_working_copy.zip` 包含归档后完成的公式排版、引文措辞修订与本次整理工具；它不是新发布版本，不代表远端主分支已同步这些修订。', '`current_working_copy.zip` 为当前工作目录副本；`published_v' + version + '_immutable.zip` 为正式发布快照。归档分配 DOI 后，主分支仅回填真实 DOI 和发布回执；既有版本标签不移动。')
        note = note.replace('本次整理不执行网站投稿或创建新归档版本。', '本次仅同步 GitHub/Zenodo，不执行网站投稿。')
        note += '\n当前发布元数据与回执位于 05_Code；历史 v0.2.0 仍单独保留。\n'
    (DEST / '00_READ_ME_FIRST.md').write_text(note, encoding='utf-8')
    files = sorted(p for p in DEST.rglob('*') if p.is_file() and p.name not in ['DELIVERY_STATUS.json', 'DELIVERY_MANIFEST_SHA256.json'])
    status = {'status': 'passed', 'protected_files_verified': len(frozen), 'classifier_fits': 0,
              'feature_extraction_runs': 0, 'github_or_biorxiv_upload_this_organization': False,
              'files_before_status_and_manifest': len(files),
              'documents': len(list((DEST / '01_Manuscripts').glob('*.docx'))),
              'figures': len(list((DEST / '02_Figures').glob('*'))),
              'current_revision_csv': len(list((DEST / '03_Tables/current_revision').glob('*.csv')))}
    (DEST / 'DELIVERY_STATUS.json').write_text(json.dumps(status, indent=2) + '\n', encoding='utf-8')
    manifest = {p.relative_to(DEST).as_posix(): digest(p) for p in sorted(DEST.rglob('*'))
                if p.is_file() and p.name != 'DELIVERY_MANIFEST_SHA256.json'}
    (DEST / 'DELIVERY_MANIFEST_SHA256.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    archive = ROOT / 'cfDNA_Final_Delivery.zip'
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(DEST.rglob('*')):
            if p.is_file():
                z.write(p, p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for name, expected in manifest.items():
            assert hashlib.sha256(z.read('FINAL_DELIVERY/' + name)).hexdigest() == expected, name
    (ROOT / 'cfDNA_Final_Delivery.zip.sha256').write_text(digest(archive) + '  cfDNA_Final_Delivery.zip\n', encoding='utf-8')
    print(json.dumps(dict(status, manifest_files=len(manifest), archive_bytes=archive.stat().st_size), indent=2))


if __name__ == '__main__':
    main()
