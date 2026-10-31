# Scientific scope and interpretation

## Central research question

How stable are the genes and pathways prioritised by supervised models when training samples and available omics layers change, and does the pattern differ between a cancer indication and an autoimmune disease?

## Design

- NSCLC: compare each molecular layer and combinations against clinical baselines for overall-survival modelling in LUAD and LUSC.
- RA: compare RNA, ATAC, H3K27ac and H3K4me3 layers and combinations for donor-level RA-versus-OA classification.
- Stability: estimate feature-selection frequency and pairwise Jaccard similarity across outer folds/repetitions.
- Cross-indication: compare stable mapped gene symbols and GO Biological Process terms. Treat results as hypothesis-generating.

## Interpretation boundaries

MOSAIC-3 is a computational benchmarking and biological prioritisation workflow. It does not establish causation, disease equivalence, clinical utility, diagnostic accuracy, treatment response, or independent prognostic validity. Cross-indication overlap is sensitive to feature identifiers, mapping choices, assay coverage and model-selection variability.
