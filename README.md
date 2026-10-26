# cfDNA Scale Tradeoff Audit v0.2.0

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

Author: Zhenxuan Ma; School of Clinical Medicine, Hubei University of Medicine, Shiyan 442000, Hubei, China

ORCID: https://orcid.org/0009-0004-5432-4999
