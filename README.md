# cfDNA Scale Tradeoff Audit

**Scale Trade-offs in Within-Sample Self-Normalization for Cross-Study Generalization
in Public cfDNA Cohorts**

Zhenxuan Ma · School of Clinical Medicine, Hubei University of Medicine, Shiyan 442000, Hubei, China

ORCID: https://orcid.org/0009-0004-5432-4999

## What the audit resolves

A pooled cfDNA benchmark score combines representation choices, study-wise ranking and comparability across fitted scoring rules. This project separates those components and measures how normalization scale and target cohort composition affect cross-study evaluation.

- **Scale trade-offs:** Global, same-arm local and Hybrid representations quantify broad fragmentation-contrast retention and residual genomic source variation.
- **Pair concordance:** within-study and cross-study cancer–control pairs explain their separate contributions to pooled AUROC.
- **Fixed anchors:** changing companion patients produces a maximum probability shift of 0.214 for the evaluated ComBat-plus-target alignment and 0 for the within-sample transforms.
- **Locked thresholds:** source-only calibration makes target specificity drift visible alongside discrimination and calibration.

## Cohort and evaluation

The analysis contains 585 source-identified participants from 628 extracted specimens and uses 477 retained hg19 autosomal 5 Mb windows. The primary bilateral cohort contains Cristiano 2019 (427) and Jiang 2015 (121), total 548. Snyder 2016 (37; 36 cancer and one healthy) is a separately held-out sensitivity cohort. Secondary Combined LOCO evaluates all 585. Study Source combines experimental and population-composition differences.

Fixed ElasticNet Combined LOCO AUROCs are Global 0.759, Hybrid 0.731, fixed LSN-ICC 0.724 and ComBat-plus-target alignment 0.543. Hybrid is a scale-fusion comparator; the study’s contribution is the mechanistic evaluation framework. Use the exact CSV cells with cohort_scope for inference.

## Code and archive

Public repository: https://github.com/mzxuan08/cfDNA-Scale-Tradeoff-Audit

Immutable software release v0.1.0: https://doi.org/10.5281/zenodo.23653766

## Publication materials

The approved manuscript is results/MANUSCRIPT_FULL_TEXT.md. Its abstract,
interpretation, legends, submission letters and supplementary guide use the same
normalization-scale and score-comparability framework.

scripts/10_finalize_editorial_materials.py rebuilds the current submission
documents from the preserved local publication snapshot. It uses the table and
figure-placement helpers in scripts/08_1_final_typesetting.py. Both require the
complete local publication workspace, author configuration, verified references,
Pandoc and document-rendering dependencies. Research replication uses the
analysis entry points below independently of this publication workspace.

## Included files and exclusions

- src/, scripts/, tests/ and config/: original research implementation and tests.
- results/tables/: exact aggregate CSV outputs and full numeric traceability tables.
- results/reports/: selected frozen validation evidence and scope alignment.
- data/ref/: compact reference windows, annotations and interval masks.
- results/archive_step06/: minimal historical OOF predictions, split provenance,
  archived models.py and checksum manifest required by the revision engine.
- results/clinical_cv_checkpoints/selection_loco_*_Proposed_LSN_ICC.parquet:
  small fixed selection masks used to audit Figure 4, not fitted model weights.

Raw BAM/FASTQ, raw fragment streams, bulk metadata, derived specimen feature
matrices, fitted model caches and large third-party archives are excluded. The
minimal historical OOF bundle is intentional provenance, not a hidden fitting
cache. Some archived reports retain original absolute local paths as historical
evidence; those paths are not portable input configuration. Configured research
input paths are relative to the project root.

## Environment

Validated analysis: Python 3.10, Windows, seed 2026, RTX 3060 Ti. Use:

```console
conda env create -f environment.yml
conda activate cfdna-batchaudit
python -m pytest tests/ -q
```

Or create a Python 3.10 virtual environment and install requirements.txt. CUDA
12.1 PyTorch wheels are pinned for Windows/Linux. Neural batches are at most
128 and hidden widths at most 256. CUDA model availability affects whether all
locked model combinations can be reproduced; unavailable hardware must be
reported, not silently replaced with fabricated scores. Exact software versions and the validated Windows/CUDA environment are recorded
in the Methods tables.

## Rebuild from the public sources

Run from the repository root, in a separate analysis working copy. Rebuilding
will intentionally regenerate results in that working copy; preserve the frozen
release if comparing hashes. Network access, original data terms and available
remote objects are required. No hundreds-of-gigabytes BAM/FASTQ download is used.

```console
python scripts/01_setup_and_ref.py
python scripts/02_fetch_and_extract.py --help
python scripts/02_fetch_and_extract.py
python scripts/03_batch_vulnerability_audit.py
python scripts/04_run_cv_benchmarks.py --step7
python scripts/05_biological_interpretability.py
```

Read --help and config/config.yaml before full extraction. The discovery path
uses published UNITE metadata and source-explicit interval URLs, with six stream
workers and checkpoints. QC and source-identifier selection must reproduce the
documented cohort; do not invent missing source labels. The Step 7 engine reuses
the verified original four-pipeline Combined LOCO predictions supplied in the
compact archive and fits only the specified revision tasks. A de novo replication
of that original benchmark is distinct from replaying its preserved predictions.

The historical archive's SHA256 manifest is checked before prediction reuse,
and models.py must match its archived version. Research feature matrices must
first be reconstructed for scripts that load patient-level features. The light
release alone does not run feature-dependent validations without those inputs.

**Manuscript regeneration:** scripts/06_compile_manuscript_assets.py is the
historical analysis compiler and can regenerate scientific tables. For editorial
work, preserve the approved CSVs and use the publication entry point described
above. The lightweight release supplies the approved text and research
implementation; feature-dependent validation requires reconstructed clinical
inputs and the original data terms.

## Fold isolation and statistical inference

All fitted preprocessing, masks, nested parameter selection and source threshold
calibration are confined to source training folds. Primary bilateral training has
only one study; healthy-source ICC is bypassed explicitly, and multi-source
ComBat/DANN are unavailable. InFold ComBat includes an unlabeled unseen-target
median/IQR alignment and is explicitly transductive. Ordinary individual-score
DeLong inference is not reported for that target-batch-fitted adaptation.

Confidence intervals use 1000 fixed-score, study/diagnosis-stratified patient
bootstraps. They do not include model refitting or uncertainty over future
studies. Target-specific DeLong p values and BH q values are in the supplied
tables; no pooled CV DeLong test is asserted. Synthetic arrays occur only in
unit tests; production clinical outputs use real public participants.

## Verification and licensing

RELEASE_MANIFEST_SHA256.json records relative paths and hashes for this prepared
release. Compare those hashes before regeneration. The submission package's
separate frozen-input audit confirms the original matrices, CSVs and figure
files were unchanged during packaging. The Step 8.1 verification run had 80 passing tests and one optional
live-network test skipped. This verified software and publication formatting;
it did not repeat clinical extraction or benchmark fitting.

Original software is MIT licensed. Third-party reference/clinical resources
retain their own conditions; see DATA_PROVENANCE_AND_TERMS.md. Cite the original
data papers, the manuscript and the versioned code/Zenodo record once released.
This is research software for retrospective evaluation of public cfDNA cohorts.
