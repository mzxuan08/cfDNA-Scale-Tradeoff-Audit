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


def test_wide_table_panels_keep_method_study_key_and_numeric_strings():
    # Use an actual frozen public-data table's first rows as the editorial fixture.
    records = builder.rows('Final_Alignment_Baseline_Comparison.csv') if (ROOT/'results/tables/Final_Alignment_Baseline_Comparison.csv').exists() else []
    if not records: pytest.skip('frozen editorial fixture not present')
    source = builder.csv_markdown(records, ['pipeline', 'center', 'auroc', 'auroc_low', 'auroc_high', 'auprc'])
    output = builder.presentation_tables(source)
    assert output.count('| Method | Study |') == 2
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
    assert author['orcid'] is None


def test_citation_and_artifact_integrity_if_package_present():
    out = ROOT/'submission_package'
    if not (out/'verification/reference_verification.json').exists():
        pytest.skip('submission artifacts are not required for research unit tests')
    reference = json.loads((out/'verification/reference_verification.json').read_text(encoding='utf-8'))
    assert len(reference['references']) == 8
    item = reference['references'][5]
    assert (item['year'], item['volume'], item['issue'], item['pages_or_article']) == (2026,'12','28','eady9432')
    author = json.loads((ROOT/'config/submission_author.json').read_text(encoding='utf-8'))
    result = builder.package_validation(author)
    assert len(result) == 8
    assert result['04_bioRxiv_Integrated_Preprint.docx']['all_five_embedded_figures_byte_identical']
