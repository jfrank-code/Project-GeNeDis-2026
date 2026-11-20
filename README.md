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

# MOSAIC-3 — Multi-Omics Stability, Integration and Clinical Insight

**GeNeDis 2026 Satellite Hackathon — Predict the Future**

MOSAIC-3 is a reproducible computational framework for studying how molecular data layers contribute to outcome modelling and how consistently supervised models select molecular features across patient partitions and omics configurations.

The project contains two distinct disease-specific analyses:

1. **NSCLC survival analysis:** TCGA lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC).
2. **Synovial classification:** Rheumatoid arthritis (RA) versus osteoarthritis (OA) using public GEO data.

A cross-disease analysis explores shared gene symbols and pathway signals. These overlaps are hypothesis-generating; they do not establish causality or biological equivalence between synovial tissue and lung tumours.

---

## Research Questions

MOSAIC-3 addresses two connected questions:

1. How much predictive information does each molecular layer add beyond a clinical or single-omics baseline?
2. How stable are the genes and pathways prioritised by supervised models when training samples or available omics layers change?

The two disease branches have different endpoints and are analysed separately: time-to-event survival in NSCLC and donor-level binary classification in RA versus OA.

---

## Cohorts and Data Sources

| Cohort | Complete samples / donors | Events / RA cases | Censored / OA controls |
|---|---:|---:|---:|
| TCGA-LUAD | 440 | 157 deaths | 283 censored |
| TCGA-LUSC | 360 | 154 deaths | 206 censored |
| GEO synovial RA/OA | 20 donors | 10 RA | 10 OA |

### 1. NSCLC: TCGA / UCSC Xena

The LUAD and LUSC cohorts were retrieved from public TCGA molecular and clinical resources through UCSC Xena. The model-complete cohorts contain primary tumour samples with aligned molecular profiles and valid overall-survival information.

The primary NSCLC analysis uses:

- **RNA:** RNA-seq-derived gene-expression profiles.
- **CNV_POS and CNV_NEG:** Positive and negative copy-number input-burden representations.
- **Methylation:** Illumina HumanMethylation450 (HM450) promoter methylation.

For pathway-level analyses, the project represents these assay families as four molecular layers: RNA, CNV_POS, CNV_NEG and methylation. Fifty Hallmark pathway features are evaluated per layer.

MC3 gene-level SNV data are retained as an exploratory resource but are excluded from the primary model because callability has not been established sufficiently for the primary benchmark.

### 2. Rheumatoid Arthritis: GEO Synovial Tissue

The RA branch uses public GEO datasets from synovial tissue:

- **GSE112656:** RNA-seq expression.
- **GSE112654:** ATAC-seq chromatin accessibility.
- **GSE112655:** ChIP-seq-derived histone-mark profiles.

The integrated cohort contains 10 RA and 10 OA donors, aligned across four molecular layers:

- RNA
- ATAC
- H3K27ac
- H3K4me3

WGBS and raw FASTQ/SRA inputs are not part of the benchmark described here.

---

## Study Design and Methods

### 1. NSCLC: Overall-Survival Modelling

The primary endpoint is overall survival, represented by follow-up time and an event indicator.

A clinical-only Cox model provides the reference baseline. Cox Elastic Net models evaluate RNA, CNV_POS, CNV_NEG and methylation individually, selected two-layer combinations and the full multi-omics configuration.

Harrell's C-index measures how well each model ranks patients by survival risk.

The updated NSCLC result package reports **20 outer repetitions × 5 folds**, yielding up to 100 outer-fold test assessments per cohort and model configuration when all folds complete successfully.

Hyperparameter tuning is performed within inner training/validation partitions. This is repeated internal validation, not independent external validation.

A separate pathway-characterisation analysis compares patients who died by 36 months with patients known to be alive at 36 months. Patients censored before 36 months are treated as unknown for this binary comparison, rather than being counted as survivors.

The full time-to-event endpoint remains the basis of Cox survival modelling.

### 2. RA vs. OA: Donor-Level Classification

The RA branch is a binary classification task, not a survival analysis.

Elastic Net logistic regression is the principal classifier, and linear SVM variants are used as sensitivity comparisons.

Model performance is described using:

- ROC-AUC
- Balanced accuracy
- Sensitivity
- Specificity
- Matthews correlation coefficient (MCC)

The updated RA figures report **10 outer repetitions**. The outer-fold count for that updated run is not specified in the report and should not be inferred.

Splitting is performed at donor level so that the molecular layers belonging to one donor remain together in the same partition.

### 3. Leakage Control and Preprocessing

The modelling design aims to prevent information from held-out test samples from influencing model training:

- For NSCLC, feature filtering, imputation, scaling, clinical encoding, supervised feature selection and hyperparameter tuning are fitted using training partitions.
- Overall-survival time and event status are outcomes, not predictors.
- For RA, supervised feature selection, variance ranking, imputation and scaling are performed within training partitions.
- All omics layers from a donor remain together during splitting.
- Hyperparameter selection is performed within the inner validation process.

**Historical RA preprocessing limitation:** An earlier RA run performed label-blind TMM/logCPM normalisation across the complete cohort before cross-validation. This is a transductive preprocessing limitation even though outcome labels were not used.

The revised RA benchmark code includes fold-specific normalisation. Archived results should only be described as fully fold-local when the run manifest confirms that they were generated with that revised workflow.

Repeated cross-validation assesses sensitivity to different partitions of the same cohort. It does not create new independent patients or replace external validation.

---

## Main Results

### 1. NSCLC Survival Prediction

The updated 20 × 5 result table reports the following highest-performing configurations:

| Cohort | Clinical baseline C-index | Best partial integration | Full multi-omics C-index |
|---|---:|---|---:|
| LUAD | 0.629 | RNA + methylation: 0.705 | **0.725** |
| LUSC | 0.550 | RNA + CNV_POS: 0.660 | **0.675** |

The full multi-omics model achieved the highest reported C-index in both cohorts.

The strongest partial integration differed by histology, suggesting that the relative contribution of molecular layers is disease-subtype dependent.

These are reported performance estimates; superiority over every alternative should not be described as statistically established without a paired inferential comparison.

### 2. NSCLC Pathway and Cross-Layer Findings

- In the v13 pathway audit, three LUAD RNA pathways — **Glycolysis, Mitotic Spindle and G2M Checkpoint** — passed both the three-year group comparison and adjusted Cox FDR criteria.
- In the same v13 audit, no LUSC pathway passed the Wilcoxon or adjusted-Cox FDR threshold.
- A consistent directional pattern in LUSC CNV_POS is descriptive, not confirmatory evidence.
- Pathway-score correlations across molecular layers were positive and mostly moderate (Spearman ρ = 0.26–0.46).
- RNA and methylation had the strongest reported correlation in both histologies.

Separate pathway-level hazard-ratio summaries report associations across RNA, CNV_POS and methylation. These are associations, not evidence of causal effects.

The report preserves distinct outputs where pathway HR estimates differ between analysis versions; they should not be averaged or silently substituted for one another.

### 3. RA Classification and Feature Stability

The updated RA figures report these approximate mean ROC-AUC values across ten outer repetitions:

| Model | Mean ROC-AUC |
|---|---:|
| RNA | **≈ 0.998** |
| RNA + ATAC | ≈ 0.984 |
| RNA + CHIP | ≈ 0.906 |
| All four layers | ≈ 0.904 |

RNA alone had the highest reported mean ROC-AUC. Adding more modalities did not improve the displayed mean AUC in this small cohort.

The RA feature-stability analysis reports the highest mean pairwise Jaccard similarity for RNA, although stability remains modest.

These findings are exploratory because there are only 20 independent donors.

### 4. Exploratory Gene Overlap

Cross-disease and cross-histology overlaps are retained as exploratory analyses, not validated biomarkers.

- The archived RA–LUAD comparison reports three shared mapped symbols: **NR0B1, PAX9 and SPINK1**.
- No mapped RA–LUSC overlap was observed in that archived comparison.
- A separate LUAD–LUSC appendix reports 16 shared symbols in an exploratory comparison of different model configurations.
- A stricter model-matched CLINICAL_RNA screen identified five shared candidates.

**Important provenance limitation:** The CSV exports used for that LUAD–LUSC appendix record five feature-selection opportunities. Those frequencies must not be presented as results from the NSCLC survival benchmark's separate 20 × 5 evaluation. A regenerated stability export is required before claiming that the overlap frequencies represent 50 or 100 selections.

---

## Limitations and Interpretation

- The RA/OA cohort contains only 20 independent donors (10 per class), creating a high feature-to-sample ratio and substantial overfitting risk.
- No independent external validation cohort is included for either disease branch.
- Internal cross-validation does not demonstrate transportability to new hospitals, populations or studies.
- Near-perfect internal RA AUC values are preliminary and must not be presented as validated diagnostic performance.
- Public observational data and pathway-level associations do not establish causal mechanisms or clinical utility.
- RA inner-validation warnings are retained for audit and require review alongside the run-level warning log.
- Historical cohort-wide RA normalisation should be distinguished from the revised fold-specific preprocessing code.
- Cross-disease and cross-histology gene overlaps are exploratory and depend on feature mapping, model configuration and selection-frequency provenance.

MOSAIC-3 is a research workflow for methodological and biological investigation. It is not a clinical diagnostic or prognostic product.

---

## Reproducibility and Repository Use

### Quick Start (Windows)

Requirements:

- Python 3.10+
- R 4.3+
- Internet access
- Sufficient free disk space and memory

TCGA HM450 and RNA inputs are large; 16 GB RAM is a practical minimum and 32 GB is recommended.

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\run_mosaic3.ps1
```

Optional flags:

```powershell
.\run_mosaic3.ps1 -ForceDownloads
.\run_mosaic3.ps1 -SkipRDependencyInstall
```

The runner is designed to acquire public TCGA/Xena and GEO inputs, prepare cohorts, run preflight quality control, execute the NSCLC survival and RA classification workflows, perform cross-indication analysis and generate reports.

It stops when a pipeline step fails and preserves completed downloads to support restart.

Raw inputs and large processed matrices are intentionally not committed. Acquisition scripts reconstruct them from public sources.

Archived metrics, out-of-fold predictions, feature-selection tables, warnings and quality-control artefacts are included for inspection.

### Reproducibility Notes

- The configuration in `config/mosaic3_config.json` controls analysis settings.
- The effective settings used by a particular result must be taken from that run's configuration and manifest, not assumed from a later edited configuration file.
- The NSCLC headline results above refer to the report's 20 × 5 result package.
- The updated RA performance figures refer to ten outer repetitions; the outer-fold count is not specified in the report.
- The RA–LUAD–LUSC feature-overlap appendix uses separate feature-stability CSV exports with five recorded selection opportunities.
- Do not conflate those exports with the NSCLC 20 × 5 survival-performance outputs.
- A rerun with a changed configuration can produce different results.
- Preserve the exact configuration, random seeds, software versions, fold-level outputs and warning logs for every reported run.

---

## Repository Structure

```text
config/                  Analysis configuration
src/scripts/python/      TCGA/Xena and GEO acquisition, preparation, QC and reporting
src/scripts/R/           NSCLC survival, RA classification and cross-indication analyses
results/                 Archived metrics, predictions, feature tables, warnings and QC
Reports/                 Scientific reports and figures
docs/                    Data provenance, reproducibility and compliance documentation
data/raw/                Created by acquisition scripts; not committed
data/processed/          Created by preparation scripts; not committed
```

---

## Data Sources and References

### Public Datasets and Cohort Studies

1. The Cancer Genome Atlas Research Network. (2012). Comprehensive genomic characterization of squamous cell lung cancers. *Nature, 489*, 519–525. https://doi.org/10.1038/nature11404

2. The Cancer Genome Atlas Research Network. (2014). Comprehensive molecular profiling of lung adenocarcinoma. *Nature, 511*, 543–550. https://doi.org/10.1038/nature13385

3. Liu, J., Lichtenberg, T., Hoadley, K. A., et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. *Cell, 173*(2), 400–416.e11. https://doi.org/10.1016/j.cell.2018.02.052

4. Goldman, M. J., Craft, B., Hastie, M., et al. (2020). Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology, 38*, 675–678. https://doi.org/10.1038/s41587-020-0546-8

5. Ai, R., Laragione, T., Hammaker, D., et al. (2018). Comprehensive epigenetic landscape of rheumatoid arthritis fibroblast-like synoviocytes. *Nature Communications, 9*, 1921. https://doi.org/10.1038/s41467-018-04310-9

6. Clough, E., Barrett, T., Wilhite, S. E., et al. (2024). NCBI GEO: archive for gene expression and epigenomics data sets: 23-year update. *Nucleic Acids Research, 52*(D1), D138–D144. https://doi.org/10.1093/nar/gkad965

### GEO Accessions Used for RA Synovial-Tissue Data

- GSE112656 — RNA-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112656
- GSE112654 — ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112654
- GSE112655 — ChIP-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112655
- GSE112658 — Parent SuperSeries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112658

### Methodological References

7. Cox, D. R. (1972). Regression Models and Life-Tables. *Journal of the Royal Statistical Society: Series B, 34*(2), 187–202. https://doi.org/10.1111/j.2517-6161.1972.tb00899.x

8. Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. *Journal of Statistical Software, 33*(1), 1–22. https://doi.org/10.18637/jss.v033.i01

9. Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2011). Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. *Journal of Statistical Software, 39*(5), 1–13. https://doi.org/10.18637/jss.v039.i05

10. Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random Survival Forests. *The Annals of Applied Statistics, 2*(3), 841–860. https://doi.org/10.1214/08-AOAS169

11. Wright, M. N., & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests in C++ and R. *Journal of Statistical Software, 77*(1). https://doi.org/10.18637/jss.v077.i01

12. Cortes, C., & Vapnik, V. (1995). Support-vector networks. *Machine Learning, 20*, 273–297. https://doi.org/10.1007/BF00994018

13. Chang, C.-C., & Lin, C.-J. (2011). LIBSVM: A library for support vector machines. *ACM Transactions on Intelligent Systems and Technology, 2*(3), Article 27. https://doi.org/10.1145/1961189.1961199

14. Harrell, F. E., Lee, K. L., & Mark, D. B. (1996). Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Statistics in Medicine, 15*(4), 361–387. https://doi.org/10.1002/(SICI)1097-0258(19960229)15:4%3C361::AID-SIM168%3E3.0.CO;2-4

15. Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology, 11*, R25. https://doi.org/10.1186/gb-2010-11-3-R25

16. Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics, 26*(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

17. The Gene Ontology Consortium. (2021). The Gene Ontology resource: enriching a GOld mine. *Nucleic Acids Research, 49*(D1), D325–D334. https://doi.org/10.1093/nar/gkaa1113

18. Huber, W., Carey, V. J., Gentleman, R., et al. (2015). Orchestrating high-throughput genomic analysis with Bioconductor. *Nature Methods, 12*, 115–121. https://doi.org/10.1038/nmeth.3252

19. Lawrence, M., Huber, W., Pagès, H., et al. (2013). Software for computing and annotating genomic ranges. *PLoS Computational Biology, 9*(8), e1003118. https://doi.org/10.1371/journal.pcbi.1003118

20. Kent, W. J., Sugnet, C. W., Furey, T. S., et al. (2002). The human genome browser at UCSC. *Genome Research, 12*(6), 996–1006. https://doi.org/10.1101/gr.229102

---

## Citation and Data Policy

MOSAIC-3 uses publicly available secondary data from TCGA/UCSC Xena and NCBI GEO. Cite the original cohort publications and dataset accessions when reusing data or derived results.

Exact dataset identifiers, source URLs, sample-selection criteria and acquisition details should be maintained in `docs/DATA_PROVENANCE.md`.

Raw source files and large processed matrices are not redistributed in this repository. They can be reconstructed using the included acquisition workflows.
