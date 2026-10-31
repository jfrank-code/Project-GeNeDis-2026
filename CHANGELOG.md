# Changelog

## Submission package - 2026-10-31

- Preserved the existing scientific workflow and archived model results.
- Corrected the PowerShell runner so it invokes scripts from the committed `src/scripts/` directory.
- Added modality-qualified feature names to the RA combined design matrix to prevent same-symbol features from different assays being conflated in future selections. Matrix values and the archived benchmark predictions were not changed.
- Audited archived RA selection counts by model/repetition/fold; duplicate records were removed from the derived stability table so selection frequency is bounded by 1.0. Historical modality attribution is explicitly marked unresolved.
- Increased the default jury rerun profile to three outer repetitions for each task; archived quick-run metrics remain unchanged and are labelled as such.
- Added an HTML-first PDF generation path with a ReportLab fallback.
- Normalised user-facing runner/log strings and removed local absolute paths from submitted provenance artefacts.
