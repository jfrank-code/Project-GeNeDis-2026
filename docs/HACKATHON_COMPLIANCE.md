# GeNeDis 2026 hackathon alignment checklist

The official challenge prioritises reproducible data workflows, public/legitimately accessible sources, explicit supervised labels, leakage checks, baseline comparisons, biological grounding, provenance, processed outputs, QC/model reports and clear limitations.

| Official evaluation/deliverable | MOSAIC-3 implementation | Status / caveat |
|---|---|---|
| Data organisation and open-source provenance | TCGA/Xena LUAD/LUSC and GEO GSE112656/GSE112654/GSE112655; source register and acquisition manifests | Implemented; raw inputs and large matrices are intentionally not committed |
| Explicit supervised targets | NSCLC overall-survival time/event; RA versus OA | Implemented |
| Baseline comparisons | Clinical-only Cox/RSF baselines for NSCLC; RNA-only baseline for RA | Implemented |
| Leakage safeguards | Training-fold preprocessing and nested tuning for NSCLC; training-fold supervised preprocessing for RA | Archived RA uses label-blind cohort-wide TMM/logCPM before CV; disclosed as a transductive limitation |
| Appropriate ML | Cox Elastic Net, Random Survival Forest sensitivity, Elastic Net logistic regression and linear SVM sensitivity | Implemented; warnings retained |
| Biological grounding | Stable gene-symbol overlap and GO Biological Process over-representation | Exploratory only; no causal or cellular-equivalence claim |
| Reproducible pipeline | Fresh-start PowerShell runner, Python/R scripts, config and install scripts | Implemented; full run may require substantial compute |
| Processed feature tables and key outputs | Cohort summaries, model metrics, OOF predictions, selected features, stability tables, enrichment and QC outputs | Included; raw and large full processed matrices are regenerated rather than bundled |
| QC/model report | English HTML and PDF report, figures, warnings and QC summaries | Included |
| Limitations | Small RA cohort, single-repeat NSCLC, normalization caveat, warnings, no external validation | Explicitly documented |
| Data ethics | Public, anonymised secondary datasets; no secrets or controlled-access data | No private keys or raw subject identifiers are intended to be included |

**Scientific integrity note:** all reported metrics are archived outputs, not synthetic values. The report distinguishes mean fold performance from single-fold maxima and does not call the findings clinically validated.

Official page: <https://2026.genedis.org/hackathon/>
