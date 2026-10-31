# MOSAIC-3 | GeNeDis 2026

**Multi-Omics Stability, Integration and Clinical Insight**  
*A reproducible, disease-agnostic framework for predictive contribution, biomarker stability, and exploratory cross-indication comparison in non-small cell lung cancer (NSCLC) and rheumatoid arthritis (RA).*

## Project at a glance

MOSAIC-3 asks two linked questions:

1. How much does each molecular layer contribute beyond a simple clinical or transcriptomic baseline?
2. How stable are the genes and pathways prioritised by supervised models when training samples and available omics layers change?

The project analyses TCGA-LUAD and TCGA-LUSC survival and GEO synovial RA-versus-OA classification. A cross-indication step compares stable gene symbols and pathway signals as exploratory biological hypotheses; it does not imply causal relationships or equivalence between synovial fibroblasts and cancer-associated fibroblasts.

## Included deliverables

- Python and R source code for fresh-start acquisition, preprocessing, QC, modelling, biomarker stability, enrichment, and report generation.
- Public-source and accession manifest in `docs/DATA_PROVENANCE.md`.
- Archived cohort/model outputs, out-of-fold predictions, feature-selection tables, warning logs, and QC artefacts under `results/`.
- Full English scientific report in `Reports/MOSAIC3_final_report.html` and `Reports/MOSAIC3_final_report.pdf`, with publication-style figures in `Reports/figures/`.
- Hackathon checklist in `docs/HACKATHON_COMPLIANCE.md` and rerun guidance in `docs/REPRODUCIBILITY.md`.

**Data packaging policy:** raw inputs and large processed matrices are intentionally not committed. The acquisition scripts rebuild them from public UCSC Xena and GEO sources. The archived result tables and feature-selection outputs are included so reviewers can inspect the completed analysis without rerunning it.

## Archived cohort sizes

| Cohort | Complete samples | Events / RA cases | Censored / OA controls |
|---|---:|---:|---:|
| TCGA-LUAD | 440 | 157 events | 283 censored |
| TCGA-LUSC | 360 | 154 events | 206 censored |
| GEO RA/OA | 20 donors | 10 RA | 10 OA |

NSCLC includes RNA-seq, gene-level CNV and HumanMethylation450. MC3 gene-level SNV acquisition is retained as an exploratory data layer but is excluded from the primary model until callability can be established. RA uses RNA-seq, ATAC-seq, H3K27ac and H3K4me3. WGBS and raw FASTQ/SRA are excluded.

## Quick start (Windows)

Requirements: Python 3.10+, R 4.3+, internet access, and sufficient free disk space and memory. TCGA HM450 and RNA inputs are large; 16 GB RAM is a practical minimum and 32 GB is recommended.

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\run_mosaic3.ps1
```

Optional flags:

```powershell
.\run_mosaic3.ps1 -ForceDownloads
.\run_mosaic3.ps1 -SkipRDependencyInstall
```

The pipeline downloads TCGA/Xena and GEO files automatically, including processed RA ATAC and ChIP peak files. No manual ATAC folder copy is required for a fresh run. GEO WGBS, FASTQ/SRA and the full large TAR archive are not downloaded.

The runner executes dependency installation, acquisition, cohort preparation, preflight validation, NSCLC survival benchmarking, RA classification, cross-indication GO analysis, and report generation. It invokes the scripts from `src/scripts/`, stops on failed steps, and preserves completed downloads to support restart.

## Recommended reviewer configuration

`config/mosaic3_config.json` is set to a more complete rerun profile:

- NSCLC: 3 outer repetitions, 5 outer folds, 3 inner folds, up to 600 features per layer and 12 lambda values.
- RA: 3 outer repetitions, 5 outer folds, 3 inner folds, up to 120 features per layer and 10 lambda values.
- Random Survival Forest and linear-SVM sensitivity options remain configurable.

This setting is deliberately more thorough than the archived quick run and may take substantially longer or require more memory. The report's archived metrics will not change until a rerun is performed.

## Results: headline observations

- **LUAD:** RSF with clinical predictors reached a mean C-index of **0.680** (fold SD 0.070), compared with the Cox clinical baseline at 0.654 (SD 0.067).
- **LUSC:** RSF with clinical predictors reached a mean C-index of **0.565** (fold SD 0.093), compared with the Cox clinical baseline at 0.504 (SD 0.075).
- **RA:** RNA-only and several RNA-containing configurations reached mean AUC 1.000 over two outer repetitions, with n=20 donors. This is highly exploratory and must not be presented as validated diagnostic performance.
- **Cross-indication:** three mapped stable symbols were shared by RA and LUAD (NR0B1, PAX9, SPINK1); no mapped RA-LUSC overlap was observed in the archived comparison.

The NSCLC values are averages over five folds from one outer repetition, not repeated-CV estimates. Individual fold maxima are not used as headline claims.

## Statistical safeguards and caveats

- NSCLC feature filtering, imputation, scaling, clinical encoding and tuning are estimated from training partitions; survival time/event are not predictors.
- RA supervised feature selection, imputation, variance ranking and scaling are training-fold specific. However, the archived run calculated RNA TMM/logCPM normalization once across the cohort without labels before CV. This is disclosed as cohort-wide, label-blind preprocessing and a potential transductive limitation.
- The current RA R source now assigns modality-qualified column names to the combined design matrix and de-duplicates repeated selection rows within each model/repetition/fold. This improves future layer attribution without changing matrix values.
- The archived RA feature IDs did not preserve modality prefixes, so exact per-layer attribution for that historical run cannot be reconstructed.
- RA has only 20 donors; inner training class counts may be below eight. glmnet warnings are preserved in `results/RA/ra_model_warnings.csv`.
- No independent external validation is included. MOSAIC-3 is a research workflow, not a clinical diagnosis or prognosis product.

See the full report and documentation for more detail. Package hardening notes are in `CHANGELOG.md`.

## Repository map

```text
config/                  Analysis configuration
src/scripts/python/      Xena/GEO acquisition, TCGA preparation, QC, report
src/scripts/R/           RA preparation, survival/classification, cross-indication
results/                 Archived metrics, predictions, feature tables and QC
Reports/                 English HTML/PDF report and figures
docs/                    Sources, reproducibility, compliance and limitations
data/raw/                Created by the fresh-start downloader (not committed)
data/processed/          Created by preparation scripts (not committed)
```

## Citation and challenge

Project: MOSAIC-3, GeNeDis 2026 Satellite Hackathon - Predict the Future.  
Official challenge: <https://2026.genedis.org/hackathon/>

Please cite the original TCGA/Xena and GEO source studies/accessions when reusing any data or results. See `docs/DATA_PROVENANCE.md`.


To regenerate the PDF with the same HTML-based layout as the included file, optionally install `requirements_pdf_optional.txt` (WeasyPrint). The report generator falls back to ReportLab when WeasyPrint is not available.
