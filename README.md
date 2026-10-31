# MOSAIC-3 | GeNeDis 2026
> [!IMPORTANT]
> ## Processed Data Matrices — Google Drive Access
>
> The complete processed feature matrices used in MOSAIC-3 are available through the following Google Drive repository. They are hosted externally due to their file size and are not stored directly in this GitHub repository.
>
> **[https://drive.google.com/file/d/1NS8K4wEDslnMjTtaiYBk2xFqqHTK8BRJ/view)**
>
> The shared data include:
>
> - **TCGA-LUAD:** RNA-seq, copy number variation (CNV), and DNA methylation (HM450).
> - **TCGA-LUSC:** RNA-seq, copy number variation (CNV), and DNA methylation (HM450).
> - **Rheumatoid Arthritis (RA):** RNA-seq, ATAC-seq, H3K27ac, and H3K4me3.
>
> These processed matrices support the archived analyses and facilitate independent inspection and reproducibility.
>
> Please consult `docs/DATA_PROVENANCE.md` for dataset accessions, original sources, and processing details.

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

## Citation, Data Sources, and Methodological References

**Project:** MOSAIC-3 — Multi-Omics Stability, Integration and Clinical Insight.  
**Challenge:** GeNeDis 2026 Satellite Hackathon — Predict the Future.  
**Official challenge page:** https://2026.genedis.org/hackathon/

### 1. Data sources and cohort publications

1. The Cancer Genome Atlas Research Network. (2012). Comprehensive genomic characterization of squamous cell lung cancers. *Nature, 489*, 519–525. https://doi.org/10.1038/nature11404

2. The Cancer Genome Atlas Research Network. (2014). Comprehensive molecular profiling of lung adenocarcinoma. *Nature, 511*, 543–550. https://doi.org/10.1038/nature13385

3. Liu, J., Lichtenberg, T., Hoadley, K. A., et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. *Cell, 173*(2), 400–416.e11. https://doi.org/10.1016/j.cell.2018.02.052

4. Goldman, M. J., Craft, B., Hastie, M., et al. (2020). Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology, 38*, 675–678. https://doi.org/10.1038/s41587-020-0546-8

5. Ai, R., Laragione, T., Hammaker, D., et al. (2018). Comprehensive epigenetic landscape of rheumatoid arthritis fibroblast-like synoviocytes. *Nature Communications, 9*, 1921. https://doi.org/10.1038/s41467-018-04310-9

6. Clough, E., Barrett, T., Wilhite, S. E., et al. (2024). NCBI GEO: archive for gene expression and epigenomics data sets: 23-year update. *Nucleic Acids Research, 52*(D1), D138–D144. https://doi.org/10.1093/nar/gkad965

7. Bibikova, M., Barnes, B., Tsan, C., et al. (2011). High density DNA methylation array with single CpG site resolution. *Genomics, 98*(4), 288–295. https://doi.org/10.1016/j.ygeno.2011.07.007

8. Mermel, C. H., Schumacher, S. E., Hill, B., et al. (2011). GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. *Genome Biology, 12*, R41. https://doi.org/10.1186/gb-2011-12-4-r41

9. Ellrott, K., Bailey, M. H., Saksena, G., et al. (2018). Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. *Cell Systems, 6*(3), 271–281.e7. https://doi.org/10.1016/j.cels.2018.03.002

**GEO accessions used for rheumatoid arthritis synovial fibroblast data:**

- GSE112656 — RNA-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112656
- GSE112654 — ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112654
- GSE112655 — ChIP-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112655
- GSE112658 — Parent SuperSeries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112658

The RA datasets originate from the integrated epigenetic landscape study by Ai et al. (2018). The analysis uses RNA-seq, ATAC-seq, and the selected H3K27ac and H3K4me3 ChIP-seq marks.

**TCGA data access:** LUAD and LUSC molecular and clinical datasets were retrieved from the UCSC Xena public data hubs. The MC3 mutation dataset was retained as an exploratory resource and excluded from the primary benchmark pending a callability assessment. Exact dataset identifiers and acquisition details are documented in `docs/DATA_PROVENANCE.md`.

### 2. Statistical and computational methods

10. Cox, D. R. (1972). Regression Models and Life-Tables. *Journal of the Royal Statistical Society: Series B, 34*(2), 187–202. https://doi.org/10.1111/j.2517-6161.1972.tb00899.x

11. Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. *Journal of Statistical Software, 33*(1), 1–22. https://doi.org/10.18637/jss.v033.i01

12. Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2011). Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. *Journal of Statistical Software, 39*(5), 1–13. https://doi.org/10.18637/jss.v039.i05

13. Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random Survival Forests. *The Annals of Applied Statistics, 2*(3), 841–860. https://doi.org/10.1214/08-AOAS169

14. Wright, M. N., & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. *Journal of Statistical Software, 77*(1). https://doi.org/10.18637/jss.v077.i01

15. Cortes, C., & Vapnik, V. (1995). Support-vector networks. *Machine Learning, 20*, 273–297. https://doi.org/10.1007/BF00994018

16. Chang, C.-C., & Lin, C.-J. (2011). LIBSVM: A library for support vector machines. *ACM Transactions on Intelligent Systems and Technology, 2*(3), Article 27. https://doi.org/10.1145/1961189.1961199

17. Harrell, F. E., Lee, K. L., & Mark, D. B. (1996). Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Statistics in Medicine, 15*(4), 361–387. https://doi.org/10.1002/(SICI)1097-0258(19960229)15:4%3C361::AID-SIM168%3E3.0.CO;2-4

18. Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology, 11*, R25. https://doi.org/10.1186/gb-2010-11-3-r25

19. Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics, 26*(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

20. The Gene Ontology Consortium. (2021). The Gene Ontology resource: enriching a GOld mine. *Nucleic Acids Research, 49*(D1), D325–D334. https://doi.org/10.1093/nar/gkaa1113

21. Huber, W., Carey, V. J., Gentleman, R., et al. (2015). Orchestrating high-throughput genomic analysis with Bioconductor. *Nature Methods, 12*, 115–121. https://doi.org/10.1038/nmeth.3252

22. Lawrence, M., Huber, W., Pagès, H., et al. (2013). Software for computing and annotating genomic ranges. *PLoS Computational Biology, 9*(8), e1003118. https://doi.org/10.1371/journal.pcbi.1003118

23. Kent, W. J., Sugnet, C. W., Furey, T. S., et al. (2002). The human genome browser at UCSC. *Genome Research, 12*(6), 996–1006. https://doi.org/10.1101/gr.229102

### 3. Data provenance and citation policy

MOSAIC-3 uses publicly available secondary data from TCGA/Xena and NCBI GEO. Please cite the original cohort publications and dataset accessions when reusing data or derived results.

Exact dataset identifiers, source URLs, sample selection criteria, and acquisition details are provided in `docs/DATA_PROVENANCE.md`.

Raw source files and large processed matrices are not redistributed in this repository. They can be reconstructed using the included acquisition workflows.

### 4. Report regeneration

To regenerate the PDF using the same HTML-based layout as the included report, optionally install `requirements_pdf_optional.txt` (WeasyPrint). The report generator falls back to ReportLab when WeasyPrint is unavailable.



Please cite the original TCGA/Xena and GEO source studies/accessions when reusing any data or results. See `docs/DATA_PROVENANCE.md`.


To regenerate the PDF with the same HTML-based layout as the included file, optionally install `requirements_pdf_optional.txt` (WeasyPrint). The report generator falls back to ReportLab when WeasyPrint is not available.
