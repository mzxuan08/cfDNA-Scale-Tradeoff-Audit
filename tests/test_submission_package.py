"""Editorial-package safeguards; no production clinical fitting or synthetic patients."""
import importlib.util
import json
import re
from pathlib import Path

import pytest

pytest.importorskip('docx', reason='optional submission authoring dependency')

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('submission_builder', ROOT / 'scripts/08_build_submission_package.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)
typeset_spec = importlib.util.spec_from_file_location('typesetting', ROOT / 'scripts/08_1_final_typesetting.py')
typesetting = importlib.util.module_from_spec(typeset_spec)
typeset_spec.loader.exec_module(typesetting)


def test_complete_wide_table_keeps_joint_key_and_numeric_strings():
    # Use an actual frozen public-data table's first rows as the editorial fixture.
    records = builder.rows('Final_Alignment_Baseline_Comparison.csv') if (ROOT/'results/tables/Final_Alignment_Baseline_Comparison.csv').exists() else []
    if not records: pytest.skip('frozen editorial fixture not present')
    source = builder.csv_markdown(records, ['pipeline', 'center', 'auroc', 'auroc_low', 'auroc_high', 'auprc'])
    output = typesetting.complete_tables(source)
    assert output.count('| Method | Study |') == 1
    assert 'panel a' not in output and 'panel b' not in output
    for record in records:
        assert f'{float(record["auroc"]):.3f}' in output
        assert f'{float(record["auroc_high"]):.3f}' in output


def test_declarations_do_not_invent_orcid_irb_or_human_confirmation():
    config = ROOT/'config/submission_author.json'
    if not config.exists(): pytest.skip('private author confirmation excluded from light release')
    author = json.loads(config.read_text(encoding='utf-8'))
    assert author['manual_review_confirmed'] is True
    statements = builder.declarations(author)
    assert 'No additional institutional ethics approval or formal exemption determination was obtained' in statements['Ethics statement']
    assert 'takes full responsibility' in statements['Declaration of generative AI and AI-assisted technologies in the manuscript preparation process']
    assert author['orcid'] == '0009-0004-5432-4999'
    assert typesetting.orcid_valid(author['orcid'])


def test_citation_and_artifact_integrity_if_package_present():
    out = ROOT/'submission_package'
    if not (out/'verification/reference_verification.json').exists():
        pytest.skip('submission artifacts are not required for research unit tests')
    reference = json.loads((out/'verification/reference_verification.json').read_text(encoding='utf-8'))
    expected = 29 if (ROOT/'results/step09_major/reports/analysis_summary.json').exists() else 28
    assert len(reference['references']) == expected
    item = next(r for r in reference['references'] if r['doi']=='10.1126/sciadv.ady9432')
    assert (item['year'], item['volume'], item['issue'], item['pages_or_article']) == (2026,'12','28','eady9432')
    author = json.loads((ROOT/'config/submission_author.json').read_text(encoding='utf-8'))
    if expected == 29:
        revision_spec=importlib.util.spec_from_file_location('major_package',ROOT/'scripts/12_build_major_revision_package.py')
        revision=importlib.util.module_from_spec(revision_spec);revision_spec.loader.exec_module(revision)
        result=revision.validate_package()
    else:
        result = typesetting.validate_package()
    assert len(result) == 8
    assert result['04_bioRxiv_Integrated_Preprint.docx']['status']=='passed'


def test_typesetting_preserves_every_original_numeric_table_cell():
    from docx import Document
    manuscript=ROOT/'results/MANUSCRIPT_FULL_TEXT.md'
    if not manuscript.exists() or not (ROOT/'submission_package/04_bioRxiv_Integrated_Preprint.docx').exists():
        pytest.skip('submission artifacts not included in light software archive')
    source=manuscript.read_text(encoding='utf-8')
    blocks=re.findall(r'(?m)(?:^\|.*\n)+',source)
    tables=Document(ROOT/'submission_package/04_bioRxiv_Integrated_Preprint.docx').tables
    assert len(blocks)==len(tables)
    assert len(tables)==(5 if (ROOT/'results/step09_major/reports/analysis_summary.json').exists() else 7)
    for block,table in zip(blocks,tables):
        rows=[line.strip().strip('|').split('|') for line in block.strip().splitlines()]
        for expected,actual in zip(rows[2:],table.rows[1:]):
            for value,cell in zip(expected,actual.cells):
                value=value.strip()
                if re.fullmatch(r'-?\d+(?:\.\d+)?',value):assert cell.text==value
