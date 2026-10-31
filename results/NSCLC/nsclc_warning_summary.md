# NSCLC warning provenance

The archived LUAD/LUSC benchmark console output recorded `glmnet` numerical/convergence warnings at the most strongly regularised (final) lambda settings. The original full per-fold warning log was not retained in the uploaded result bundle, so this file intentionally does not invent warning counts or model/fold attribution. The RA warning records are available separately at `results/RA/ra_model_warnings.csv`.

Interpretation: the archived NSCLC fold-average metrics are reported as observed, but the affected penalised-model fits should be treated cautiously. A jury rerun will regenerate console warnings under the expanded configuration.
