"""Execute Step 9 reviewer analyses on 585 existing real participants."""
from pathlib import Path
import sys,shutil,json,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.common import load_config,sha256,write_json
from src.major_revision import RevisionAnalysis

def main():
    cfg=load_config();root=Path(cfg['_root']);archive=root/'results/archive_before_step09_major'
    # Preserve the exact submitted publication files and source before revision.
    if not (archive/'manifest.json').exists():
        paths=[root/'results/MANUSCRIPT_FULL_TEXT.md',root/'results/reports/MANUSCRIPT_MASTER_DRAFT.md',root/'src/evaluation.py',root/'src/plotting.py',root/'tests/test_submission_package.py']
        paths+=list((root/'submission_package').glob('*.docx'))+list((root/'submission_package').glob('*.md'))+[root/'submission_package/04_bioRxiv_Integrated_Preprint.pdf']
        for p in paths:
            dest=archive/p.relative_to(root);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
        write_json({str(p.relative_to(root)):sha256(p) for p in paths},archive/'manifest.json')
    baseline_path=root/'work/step09/scientific_baseline.json'
    if baseline_path.exists():
        baseline=json.loads(baseline_path.read_text(encoding='utf-8'))
    else:
        # The portable software bundle excludes machine-specific work/ caches.
        # Freeze its reconstructed clinical inputs and supplied score provenance
        # before any new analysis; load_clinical independently validates stamps.
        paths=[Path(cfg['paths']['processed'])/cfg['extraction']['final_file'],
               Path(cfg['paths']['metadata'])]
        paths+=list((root/'data/ref').glob('*'))+list((root/'results/tables').glob('*'))
        baseline={str(p.relative_to(root)):sha256(p) for p in paths if p.is_file()}
        baseline_path.parent.mkdir(parents=True,exist_ok=True)
        write_json(baseline,baseline_path)
    # evaluation.py is expressly extended by this revision. Its historical copy
    # is supplied in the previous public tag and pre-revision archives.
    protected={p:h for p,h in baseline.items() if p!='src/evaluation.py'}
    for p,h in protected.items():assert sha256(root/p)==h,p
    analysis=RevisionAnalysis(cfg);analysis.run()
    for p,h in protected.items():assert sha256(root/p)==h,p
    write_json({'status':'passed','protected_unchanged':len(protected),'authorized_code_extension':['src/evaluation.py'],'changed_scientific_data':[]},analysis.out/'reports/input_integrity.json')
    report='# Step 9 major revision analysis\n\n'+(analysis.out/'reports/analysis_summary.json').read_text(encoding='utf-8')+'\n\nAll new outputs are isolated under results/step09_major/. Submitted v1 data, predictions, tables and original figures are retained.\n'
    (root/'results/reports/step_09_major_revision_report.md').write_text(report,encoding='utf-8')
    from src.revision_manuscript import compile_revision
    from src.plotting import major_revision_figures
    compile_revision(root)
    major_revision_figures(analysis)
if __name__=='__main__':main()
