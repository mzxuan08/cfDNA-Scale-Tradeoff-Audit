# Step 9 major revision analysis

{
  "participants": 585,
  "windows": 477,
  "ratio_missing": 0,
  "ratio_total": 279045,
  "length_missing": 0,
  "seconds": 151.23499999998603,
  "seed": 2026,
  "all_fits_converged": true,
  "stage_missing": 331,
  "matrix_sha256": "be7aed9e3ade19f3f492f6700251e54511c382e75f9781cf8f877f2f3b94c612",
  "clinical_refitting": "new source-only ElasticNet calibration and sensitivity analyses; old frozen predictions retained",
  "imputation": "per-feature source-training median; all-missing source columns use zero; post-transform medians; no target fit",
  "elasticnet": {
    "C": 1.0,
    "l1_ratio": 0.5,
    "solver": "saga",
    "max_iter": 20000,
    "tol": 0.001
  }
}

All new outputs are isolated under results/step09_major/. Submitted v1 data, predictions, tables and original figures are retained.

## Final verification

{
  "status": "passed",
  "revision": "Step 9 major revision",
  "participants": 585,
  "primary_participants": 548,
  "clinical_extraction_repeated": false,
  "new_source_only_calibration_and_sensitivity_fits": true,
  "original_protected_files_unchanged": 1411,
  "old_submitted_version": "BIORXIV/2026/761726 Version 1",
  "revision_uploaded": false,
  "new_zenodo_release": false,
  "tests": {
    "tests": 94,
    "errors": 0,
    "failures": 0,
    "skipped": 1,
    "passed": 93
  },
  "rendered_documents": [
    {
      "file": "01_Title_Page.docx",
      "sha256": "e021b9582e018b307195364d8684f20751134379c60866493cc93b6db20f6242",
      "returncode": 0,
      "pages": 2,
      "seconds": 10.76500000001397
    },
    {
      "file": "02_Highlights.docx",
      "sha256": "92cbf90cc6c034fcd38049696e8975649d84bb609df77540abf6595175cebda8",
      "returncode": 0,
      "pages": 1,
      "seconds": 10.5
    },
    {
      "file": "03_Manuscript_DoubleSpaced_LineNumbered.docx",
      "sha256": "6fe24a40b0eb377f5ee025a58b945ceb2313ba8fbf2964bb73a979d3e9f67236",
      "returncode": 0,
      "pages": 24,
      "seconds": 11.78100000001723
    },
    {
      "file": "04_bioRxiv_Integrated_Preprint.docx",
      "sha256": "d1da6bd3447b8c626fdc8c652610b3eb3c8487020774150c5851a3dca08378d5",
      "returncode": 0,
      "pages": 21,
      "seconds": 14.562999999965541
    },
    {
      "file": "05_Supplementary_Material.docx",
      "sha256": "7a199fe73a4fbf868ed83c9dff6605095360f252f8e9d0c06dd9137a70da18e8",
      "returncode": 0,
      "pages": 4,
      "seconds": 11.39100000000326
    },
    {
      "file": "06_Cover_Letter.docx",
      "sha256": "ca38b847e931bb16ebaac9a881f84cb4d2e39bf371b9fa4da71a19b9b23ed43b",
      "returncode": 0,
      "pages": 1,
      "seconds": 10.03100000001723
    },
    {
      "file": "06b_Cover_Letter_Alternative_CBM.docx",
      "sha256": "b092e6ba4b4e52a58d8e40365b43754be93266f761645736e9b14609bb6823e3",
      "returncode": 0,
      "pages": 1,
      "seconds": 9.59399999998277
    },
    {
      "file": "07_Declaration_of_Interest.docx",
      "sha256": "f7769c096ede290155814d582a77ce3eff9b23bca51740e828e607799df44f67",
      "returncode": 0,
      "pages": 1,
      "seconds": 9.202999999979511
    }
  ],
  "pdf_pages": 21,
  "pdf_sha256": "a8aed356c6f94d874710260d02ab6ecb1840b17b0fd4550adfe97c4b42a2ff09",
  "previous_submitted_pdf_sha256": "6c2c9d1b3af8f3be4cc15c6899812c2559e85f3188bd6c22fe66f467c17c298d",
  "figures": [
    {
      "figure": "Figure1_Batch_Vulnerability_Audit",
      "minimum_text_pt": 8.8,
      "width_inches": 6.5,
      "height_inches": 7.5,
      "png_dpi": 300,
      "embedded_width_inches": 6.5,
      "png_width_pixels": 1960,
      "minimum_embedded_text_pt": 8.755102040816327,
      "png_sha256": "5526cd7a308fa4b3a4a6532d7ee0bb26de649cf19cfd9dee4107d325df8a03dd"
    },
    {
      "figure": "Figure2_Primary_Ranking_Threshold_Transport",
      "minimum_text_pt": 8.8,
      "width_inches": 6.5,
      "height_inches": 6.7,
      "png_dpi": 300,
      "embedded_width_inches": 6.5,
      "png_width_pixels": 1960,
      "minimum_embedded_text_pt": 8.755102040816327,
      "png_sha256": "90f38eb72e58366f0fbfc4693ce62d9b07c6be6743449f36e5f7437ea50fe3b3"
    },
    {
      "figure": "Figure3_Companion_Dependence_and_Invariance",
      "minimum_text_pt": 8.8,
      "width_inches": 6.5,
      "height_inches": 7.0,
      "png_dpi": 300,
      "embedded_width_inches": 6.5,
      "png_width_pixels": 1960,
      "minimum_embedded_text_pt": 8.755102040816327,
      "png_sha256": "993f2a752fe2b52672ed60ff649cb87657f3a524bdcfb7c565d8bafc8a90fd9e"
    },
    {
      "figure": "Figure4_Genomic_Context_and_Ranks",
      "minimum_text_pt": 8.8,
      "width_inches": 6.5,
      "height_inches": 7.6,
      "png_dpi": 300,
      "embedded_width_inches": 6.5,
      "png_width_pixels": 1960,
      "minimum_embedded_text_pt": 8.755102040816327,
      "png_sha256": "73f4d75ce1033cb477c54e0ed67d9b4d6ed5c580552b8833b282fba163eaeb10"
    },
    {
      "figure": "FigureS5_Score_Comparability_and_Calibration",
      "minimum_text_pt": 8.8,
      "width_inches": 6.5,
      "height_inches": 6.5,
      "png_dpi": 300,
      "embedded_width_inches": 6.5,
      "png_width_pixels": 1960,
      "minimum_embedded_text_pt": 8.755102040816327,
      "png_sha256": "d88e043cb293a04dd9b44b95ae1b1cb188c352aec6899e1ae7b45fd5f476f7d9"
    },
    {
      "figure": "FigureS6_Secondary_Classifiers",
      "minimum_text_pt": 8.8,
      "width_inches": 6.5,
      "height_inches": 3.3,
      "png_dpi": 300,
      "embedded_width_inches": 6.5,
      "png_width_pixels": 1960,
      "minimum_embedded_text_pt": 8.755102040816327,
      "png_sha256": "dfcdbf5c098431b74ad38973ec8f1435cdcc5e139d9a34f28aeb4fb3ee6095e2"
    }
  ],
  "references": 29,
  "manual_review_of_this_new_revision": "pending author review before submission"
}
