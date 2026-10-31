# Reproducibility guide

## Fresh-start workflow

1. Install Python 3.10+ and R 4.3+.
2. Clone this repository and open PowerShell in the repository root.
3. Run `Set-ExecutionPolicy -Scope Process Bypass` if required by the local PowerShell policy.
4. Execute `./run_mosaic3.ps1`.
5. Inspect `results/QC/preflight_summary.json`, the result tables, warning logs and generated report.

The workflow downloads public TCGA/Xena matrices and GEO processed files, prepares matched cohorts, runs QC, models, cross-indication enrichment and regenerates the report. It preserves successfully downloaded files to support recovery after an interrupted run.

## Compute profile

The checked-in configuration is intended for a fuller reviewer rerun: three outer repetitions for both NSCLC and RA, five outer folds and three inner folds. NSCLC uses up to 600 features per layer and 12 lambda values; RA uses up to 120 features per layer and 10 lambda values. Random Survival Forest and linear SVM sensitivity settings are configurable. The run may be expensive on a laptop, especially during HM450 processing and nested model fitting.

## Archived run versus recommended rerun

The submitted metrics were generated under a constrained compute/time budget:

- NSCLC: one outer repetition, five outer folds, two inner folds, up to 250 features/layer and eight lambda values.
- RA: two outer repetitions, five outer folds, three inner folds and up to 120 features/layer.

The report explicitly presents these archived metrics. The new configuration increases the repetition count to three for reviewer execution. Do not describe the archived scores as if they were generated under the expanded configuration.

## Determinism and outputs

The project seed is recorded in `config/mosaic3_config.json`. Cohort manifests and acquisition manifests are regenerated during the workflow. Results are written beneath `results/`; reports are generated beneath `Reports/`. The main result files include fold-level and pooled out-of-fold predictions, model metrics, incremental layer contribution, selected features, stability summaries and warnings.

## Important statistical caveats

- The archived NSCLC run has only one outer repetition. Five-fold standard deviations describe variation across those five test folds, not between independent full CV repetitions.
- RA has 20 donors only. Inner folds can have small class counts and warnings are retained in `ra_model_warnings.csv`.
- In the archived RA run, TMM/logCPM was calculated once on the complete cohort without using labels. Although label-blind, this is cohort-wide preprocessing. Treat it as a potential transductive limitation and use a train-reference-only normalization scheme for a strict confirmatory analysis.
- The archived RA feature IDs did not preserve their modality prefixes. The submitted R source fixes this for future runs; the historical metrics are not recomputed.
- No independent external validation is claimed.

## Warning logs

The complete archived RA warning table is included at `results/RA/ra_model_warnings.csv`. The archived NSCLC console output reported numerical/convergence warnings at highly regularised lambda settings, but the full per-fold NSCLC warning log was not retained; a transparent provenance note is included at `results/NSCLC/nsclc_warning_summary.md`.
