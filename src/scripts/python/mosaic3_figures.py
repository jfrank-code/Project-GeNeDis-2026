#!/usr/bin/env python3
"""Generate reproducible MOSAIC-3 figures from the current CSV result artefacts.

The module never simulates or imputes model results. A figure is created only
when its corresponding result table contains the required columns and rows.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception as exc:
        warnings.warn(f"Could not read {path}: {exc}")
        return pd.DataFrame()


def numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    z = df.copy()
    for col in cols:
        if col in z.columns:
            z[col] = pd.to_numeric(z[col], errors="coerce")
    return z


def save(fig: plt.Figure, path: Path, dpi: int = 190) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def add_record(records: list[dict[str, str]], section: str, title: str,
               path: Path | None, source: str) -> None:
    if path is not None and path.exists():
        records.append({"section": section, "title": title,
                        "filename": path.name, "source": source})


def _clean_models(df: pd.DataFrame, metric: str = "model") -> pd.DataFrame:
    if df.empty or metric not in df.columns:
        return pd.DataFrame()
    z = df.copy()
    z[metric] = z[metric].astype(str)
    return z


def plot_nsclc_summary(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    rows = []
    for cohort in ("LUAD", "LUSC"):
        d = numeric(read_csv(root / "results" / "NSCLC" / cohort / "model_metrics.csv"),
                    ["mean_c_index", "sd_c_index"])
        if not d.empty and {"model", "mean_c_index"}.issubset(d.columns):
            d = d.dropna(subset=["mean_c_index"])
            d["cohort"] = cohort
            rows.append(d)
    if not rows:
        return
    z = pd.concat(rows, ignore_index=True)
    models = list(dict.fromkeys(z["model"].astype(str).tolist()))
    fig, ax = plt.subplots(figsize=(max(11, len(models) * 0.8), 5.8))
    x = np.arange(len(models)); width = 0.36
    for idx, cohort in enumerate(("LUAD", "LUSC")):
        q = z[z.cohort == cohort].set_index("model").reindex(models)
        vals = q["mean_c_index"].to_numpy(dtype=float)
        err = q["sd_c_index"].fillna(0).to_numpy(dtype=float) if "sd_c_index" in q else np.zeros(len(q))
        ax.bar(x + (idx - 0.5) * width, vals, width, yerr=err, capsize=3, label=cohort)
    ax.axhline(0.5, linestyle="--", linewidth=1, label="0.5 reference")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Harrell's C-index (mean ± SD)")
    ax.set_title("NSCLC survival performance by model", loc="left", weight="bold")
    ax.set_xticks(x); ax.set_xticklabels(models, rotation=35, ha="right")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=0.2)
    add_record(records, "NSCLC", "LUAD/LUSC C-index by model",
               save(fig, out / "nsclc_cindex_by_model.png"), "LUAD/LUSC model_metrics.csv")


def plot_nsclc_deltas(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    rows = []
    for cohort in ("LUAD", "LUSC"):
        d = numeric(read_csv(root / "results" / "NSCLC" / cohort / "layer_contribution_vs_clinical.csv"),
                    ["delta_c_index_vs_clinical", "sd_delta_c_index"])
        if not d.empty and {"model", "delta_c_index_vs_clinical"}.issubset(d.columns):
            d["cohort"] = cohort; rows.append(d)
    if not rows:
        return
    z = pd.concat(rows, ignore_index=True)
    models = list(dict.fromkeys(z["model"].astype(str).tolist()))
    fig, ax = plt.subplots(figsize=(max(10, len(models) * 0.8), 5.2))
    x = np.arange(len(models)); width = 0.36
    for idx, cohort in enumerate(("LUAD", "LUSC")):
        q = z[z.cohort == cohort].set_index("model").reindex(models)
        ax.bar(x + (idx - 0.5) * width, q.delta_c_index_vs_clinical, width,
               yerr=q.sd_delta_c_index.fillna(0) if "sd_delta_c_index" in q else None,
               capsize=2, label=cohort)
    ax.axhline(0, color="black", linewidth=1)
    ax.set_ylabel("Δ C-index vs clinical-only model")
    ax.set_title("Incremental contribution relative to clinical baseline", loc="left", weight="bold")
    ax.set_xticks(x); ax.set_xticklabels(models, rotation=35, ha="right")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=0.2)
    add_record(records, "NSCLC", "Incremental C-index over clinical baseline",
               save(fig, out / "nsclc_delta_cindex_vs_clinical.png"),
               "LUAD/LUSC layer_contribution_vs_clinical.csv")


def plot_nsclc_fold_distributions(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    for cohort in ("LUAD", "LUSC"):
        d = numeric(read_csv(root / "results" / "NSCLC" / cohort / "metrics_by_fold.csv"), ["c_index", "repetition", "fold"])
        if d.empty or not {"model", "c_index"}.issubset(d.columns):
            continue
        d = d.dropna(subset=["c_index"])
        if d.empty:
            continue
        means = d.groupby("model").c_index.mean().sort_values(ascending=False)
        models = means.index.tolist()
        values = [d.loc[d.model.astype(str) == str(m), "c_index"].to_numpy() for m in models]
        fig, ax = plt.subplots(figsize=(max(9, len(models) * 0.75), 5.2))
        ax.boxplot(values, labels=models, showmeans=True)
        ax.set_ylim(0, 1); ax.set_ylabel("Harrell's C-index per outer fold")
        ax.set_title(f"{cohort}: outer-fold C-index distribution", loc="left", weight="bold")
        ax.tick_params(axis="x", labelrotation=35); ax.grid(axis="y", alpha=0.2)
        add_record(records, "NSCLC", f"{cohort} outer-fold performance distribution",
                   save(fig, out / f"{cohort.lower()}_cindex_by_fold.png"),
                   f"results/NSCLC/{cohort}/metrics_by_fold.csv")


def plot_nsclc_feature_stability(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    for cohort in ("LUAD", "LUSC"):
        d = numeric(read_csv(root / "results" / "NSCLC" / cohort / "feature_stability.csv"), ["selection_frequency"])
        if d.empty or not {"model", "feature_id", "selection_frequency"}.issubset(d.columns):
            continue
        preference = ["CLINICAL_ALL", "OMICS_ALL", "CLINICAL_RNA", "RNA"]
        chosen = next((m for m in preference if (d.model.astype(str) == m).any()), str(d.model.iloc[0]))
        q = d[d.model.astype(str) == chosen].dropna(subset=["selection_frequency"])
        q = q.sort_values("selection_frequency", ascending=False).head(20).sort_values("selection_frequency")
        if q.empty:
            continue
        fig, ax = plt.subplots(figsize=(9.5, max(4.5, len(q) * 0.28)))
        ax.barh(q.feature_id.astype(str), q.selection_frequency)
        ax.set_xlim(0, 1); ax.set_xlabel("Selection frequency across completed outer folds")
        ax.set_title(f"{cohort}: top stable features ({chosen})", loc="left", weight="bold")
        ax.grid(axis="x", alpha=0.2)
        add_record(records, "NSCLC", f"{cohort} stable features ({chosen})",
                   save(fig, out / f"{cohort.lower()}_feature_stability_top20.png"),
                   f"results/NSCLC/{cohort}/feature_stability.csv")


def plot_nsclc_cohort_qc(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = numeric(read_csv(root / "results" / "NSCLC" / "nsclc_cohort_summary.csv"),
                ["n_model_complete", "events", "censored"])
    if d.empty or not {"cohort", "events", "censored"}.issubset(d.columns):
        return
    d = d.drop_duplicates("cohort")
    if d.empty:
        return
    x = np.arange(len(d)); width = 0.36
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(x - width / 2, d.events, width, label="Events")
    ax.bar(x + width / 2, d.censored, width, label="Censored")
    ax.set_xticks(x); ax.set_xticklabels(d.cohort.astype(str))
    ax.set_ylabel("Number of cases"); ax.set_title("NSCLC cohort event/censoring composition", loc="left", weight="bold")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=0.2)
    add_record(records, "NSCLC", "LUAD/LUSC event and censoring counts",
               save(fig, out / "nsclc_event_censoring_counts.png"),
               "results/NSCLC/nsclc_cohort_summary.csv")


def plot_ra_performance(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = numeric(read_csv(root / "results" / "RA" / "ra_model_metrics.csv"),
                ["mean_auc", "sd_auc", "mean_balanced_accuracy", "sd_balanced_accuracy"])
    if d.empty or not {"model", "mean_auc", "mean_balanced_accuracy"}.issubset(d.columns):
        return
    d = d.sort_values("mean_auc", ascending=False)
    x = np.arange(len(d)); width = 0.36
    fig, ax = plt.subplots(figsize=(max(10, len(d) * 0.85), 5.3))
    ax.bar(x - width / 2, d.mean_auc, width, yerr=d.sd_auc.fillna(0) if "sd_auc" in d else None,
           capsize=3, label="ROC-AUC")
    ax.bar(x + width / 2, d.mean_balanced_accuracy, width,
           yerr=d.sd_balanced_accuracy.fillna(0) if "sd_balanced_accuracy" in d else None,
           capsize=3, label="Balanced accuracy")
    ax.set_ylim(0, 1.05); ax.set_ylabel("Mean across outer repetitions")
    ax.set_title("RA vs OA: donor-level multi-omics performance", loc="left", weight="bold")
    ax.set_xticks(x); ax.set_xticklabels(d.model.astype(str), rotation=35, ha="right")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=0.2)
    add_record(records, "RA", "RA AUC and balanced accuracy by model",
               save(fig, out / "ra_auc_balanced_accuracy.png"), "results/RA/ra_model_metrics.csv")

    metric_cols = [c for c in ("mean_sensitivity", "mean_specificity", "mean_mcc") if c in d.columns]
    if metric_cols:
        q = d.set_index("model")[metric_cols]
        fig, ax = plt.subplots(figsize=(max(9, len(q) * 0.8), 5.0))
        q.plot(kind="bar", ax=ax)
        ax.set_ylim(-1, 1.05); ax.set_ylabel("Mean metric")
        ax.set_title("RA: sensitivity, specificity and MCC", loc="left", weight="bold")
        ax.tick_params(axis="x", labelrotation=35); ax.legend(frameon=False)
        ax.grid(axis="y", alpha=0.2)
        add_record(records, "RA", "RA sensitivity, specificity and MCC",
                   save(fig, out / "ra_sensitivity_specificity_mcc.png"), "results/RA/ra_model_metrics.csv")


def plot_ra_repeats(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = numeric(read_csv(root / "results" / "RA" / "ra_metrics_by_repeat.csv"),
                ["repetition", "mean_auc", "mean_balanced_accuracy"])
    if d.empty or not {"repetition", "model", "mean_auc"}.issubset(d.columns):
        return
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    for model, q in d.groupby("model", sort=False):
        q = q.sort_values("repetition")
        ax.plot(q.repetition, q.mean_auc, marker="o", linewidth=1.5, label=str(model))
    ax.set_ylim(0, 1.05); ax.set_xlabel("Outer repetition"); ax.set_ylabel("Pooled out-of-fold ROC-AUC")
    ax.set_title("RA: model performance across repeated outer CV", loc="left", weight="bold")
    ax.grid(alpha=0.2); ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    add_record(records, "RA", "RA AUC across outer repetitions",
               save(fig, out / "ra_auc_by_repetition.png"), "results/RA/ra_metrics_by_repeat.csv")


def plot_ra_contribution(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = numeric(read_csv(root / "results" / "RA" / "ra_layer_contribution.csv"),
                ["delta_auc_vs_RNA", "sd_delta_auc_vs_RNA"])
    if d.empty or not {"model", "delta_auc_vs_RNA"}.issubset(d.columns):
        return
    d = d.sort_values("delta_auc_vs_RNA")
    fig, ax = plt.subplots(figsize=(9.5, max(4.5, len(d) * 0.35)))
    ax.barh(d.model.astype(str), d.delta_auc_vs_RNA,
            xerr=d.sd_delta_auc_vs_RNA if "sd_delta_auc_vs_RNA" in d else None, capsize=2)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("Δ ROC-AUC vs RNA-only (paired by repetition)")
    ax.set_title("Incremental contribution beyond RNA-only baseline", loc="left", weight="bold")
    ax.grid(axis="x", alpha=0.2)
    add_record(records, "RA", "RA incremental contribution vs RNA-only",
               save(fig, out / "ra_delta_auc_vs_rna.png"), "results/RA/ra_layer_contribution.csv")


def plot_ra_stability(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = numeric(read_csv(root / "results" / "RA" / "ra_stability_jaccard.csv"),
                ["mean_jaccard", "sd_jaccard"])
    if not d.empty and {"model", "mean_jaccard"}.issubset(d.columns):
        d = d.sort_values("mean_jaccard")
        fig, ax = plt.subplots(figsize=(8.5, max(4.2, len(d) * 0.33)))
        ax.barh(d.model.astype(str), d.mean_jaccard,
                xerr=d.sd_jaccard if "sd_jaccard" in d else None, capsize=2)
        ax.set_xlim(0, 1); ax.set_xlabel("Mean pairwise Jaccard similarity")
        ax.set_title("RA selected-feature-set stability", loc="left", weight="bold")
        ax.grid(axis="x", alpha=0.2)
        add_record(records, "RA", "RA feature-set Jaccard stability",
                   save(fig, out / "ra_jaccard_stability.png"), "results/RA/ra_stability_jaccard.csv")

    f = numeric(read_csv(root / "results" / "RA" / "ra_stable_features.csv"), ["selection_frequency", "selected_count"])
    if f.empty or not {"model", "feature_id", "selection_frequency"}.issubset(f.columns):
        return
    model = "ALL" if (f.model.astype(str) == "ALL").any() else str(f.model.iloc[0])
    q = f[f.model.astype(str) == model].dropna(subset=["selection_frequency"])
    q = q.sort_values("selection_frequency", ascending=False).head(25).sort_values("selection_frequency")
    if not q.empty:
        fig, ax = plt.subplots(figsize=(9.5, max(5, len(q) * 0.29)))
        ax.barh(q.feature_id.astype(str), q.selection_frequency)
        ax.set_xlim(0, 1); ax.set_xlabel("Selection frequency across completed outer folds")
        ax.set_title(f"RA: top stable features ({model} model)", loc="left", weight="bold")
        ax.grid(axis="x", alpha=0.2)
        add_record(records, "RA", f"RA top stable features ({model})",
                   save(fig, out / "ra_top_stable_features.png"), "results/RA/ra_stable_features.csv")

    if "modality" in f.columns:
        q2 = f[(f.model.astype(str) == model) & (pd.to_numeric(f.selection_frequency, errors="coerce") > 0)]
        q2 = q2.groupby("modality").size().sort_values(ascending=False)
        if not q2.empty:
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.bar(q2.index.astype(str), q2.values)
            ax.set_ylabel("Number of selected features")
            ax.set_title(f"RA: selected-feature counts by modality ({model})", loc="left", weight="bold")
            ax.tick_params(axis="x", labelrotation=25); ax.grid(axis="y", alpha=0.2)
            add_record(records, "RA", f"RA selected-feature counts by modality ({model})",
                       save(fig, out / "ra_features_by_modality.png"), "results/RA/ra_stable_features.csv")


def roc_points(y: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    ok = np.isfinite(y) & np.isfinite(score)
    y = y[ok].astype(int); score = score[ok].astype(float)
    if len(y) == 0 or len(np.unique(y)) < 2:
        return None
    # Group tied thresholds so each ROC step is well-defined.
    thresholds = np.unique(score)[::-1]
    tpr = [0.0]; fpr = [0.0]
    positives = max(1, int((y == 1).sum())); negatives = max(1, int((y == 0).sum()))
    for threshold in thresholds:
        pred = score >= threshold
        tpr.append(float(((pred) & (y == 1)).sum()) / positives)
        fpr.append(float(((pred) & (y == 0)).sum()) / negatives)
    tpr.append(1.0); fpr.append(1.0)
    order = np.argsort(fpr, kind="stable")
    return np.asarray(fpr)[order], np.asarray(tpr)[order]


def plot_ra_roc(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = read_csv(root / "results" / "RA" / "ra_oof_predictions.csv")
    required = {"model", "repetition", "truth", "probability"}
    if d.empty or not required.issubset(d.columns):
        return
    d = d.copy()
    d["probability"] = pd.to_numeric(d.probability, errors="coerce")
    d["truth"] = pd.to_numeric(d.truth, errors="coerce")
    if "score_type" in d.columns:
        d = d[d.score_type.astype(str).str.lower().eq("probability")]
    d = d.dropna(subset=["probability", "truth", "repetition"])
    if d.empty:
        return
    preference = ["RNA", "RNA+ATAC", "RNA+CHIP", "ALL"]
    present = set(d.model.astype(str))
    models = [m for m in preference if m in present]
    if not models:
        models = list(dict.fromkeys(d.model.astype(str)))[:4]
    grid = np.linspace(0, 1, 101)
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    plotted = 0
    for model in models:
        curves = []
        for _, q in d[d.model.astype(str) == model].groupby("repetition"):
            pts = roc_points(q.truth.to_numpy(), q.probability.to_numpy())
            if pts is None:
                continue
            fpr, tpr = pts
            curves.append(np.interp(grid, fpr, tpr))
        if not curves:
            continue
        mean_tpr = np.mean(np.vstack(curves), axis=0)
        auc_value = float(np.trapz(mean_tpr, grid))
        ax.plot(grid, mean_tpr, linewidth=2, label=f"{model} (mean AUC≈{auc_value:.3f}; {len(curves)} reps)")
        plotted += 1
    if not plotted:
        plt.close(fig); return
    ax.plot([0, 1], [0, 1], linestyle="--", linewidth=1, label="Chance")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_xlabel("False-positive rate"); ax.set_ylabel("True-positive rate")
    ax.set_title("RA: repeated out-of-fold ROC curves", loc="left", weight="bold")
    ax.legend(frameon=False); ax.grid(alpha=0.2)
    add_record(records, "RA", "RA repeated out-of-fold ROC curves",
               save(fig, out / "ra_oof_roc_curves.png"), "results/RA/ra_oof_predictions.csv")


def plot_ra_feature_qc(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = numeric(read_csv(root / "results" / "RA" / "ra_feature_qc_by_fold.csv"),
                ["repetition", "outer_fold", "n_features_total", "n_nonzero_coefficients"])
    if d.empty or not {"model", "n_features_total", "n_nonzero_coefficients"}.issubset(d.columns):
        return
    models = list(dict.fromkeys(d.model.astype(str)))
    q = d.groupby("model")["n_nonzero_coefficients"].apply(list)
    order = [m for m in models if m in q.index]
    fig, ax = plt.subplots(figsize=(max(9, len(order) * 0.8), 5.0))
    ax.boxplot([q[m] for m in order], labels=order, showmeans=True)
    ax.set_ylabel("Non-zero Elastic Net coefficients per completed outer fold")
    ax.set_title("RA: selected-feature count by model and outer fold", loc="left", weight="bold")
    ax.tick_params(axis="x", labelrotation=35); ax.grid(axis="y", alpha=0.2)
    add_record(records, "RA", "RA non-zero feature counts by fold",
               save(fig, out / "ra_nonzero_features_by_fold.png"), "results/RA/ra_feature_qc_by_fold.csv")


def plot_ra_warnings(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = read_csv(root / "results" / "RA" / "ra_model_warnings.csv")
    if d.empty or "context" not in d.columns:
        return
    q = d.context.astype(str).value_counts().head(20).sort_values()
    if q.empty:
        return
    fig, ax = plt.subplots(figsize=(9, max(4.5, len(q) * 0.27)))
    ax.barh(q.index, q.values)
    ax.set_xlabel("Number of recorded warning rows")
    ax.set_title("RA: model warnings by analysis context", loc="left", weight="bold")
    ax.grid(axis="x", alpha=0.2)
    add_record(records, "RA", "RA warning counts by context",
               save(fig, out / "ra_model_warnings_by_context.png"), "results/RA/ra_model_warnings.csv")


def plot_cross_indication(root: Path, out: Path, records: list[dict[str, str]]) -> None:
    d = numeric(read_csv(root / "results" / "CROSS_INDICATION" / "shared_stable_genes.csv"),
                ["RA_selection_frequency", "NSCLC_selection_frequency"])
    if not d.empty and {"gene_symbol", "RA_selection_frequency", "NSCLC_selection_frequency"}.issubset(d.columns):
        fig, ax = plt.subplots(figsize=(7.8, 6.0))
        ax.scatter(d.RA_selection_frequency, d.NSCLC_selection_frequency)
        for _, r in d.iterrows():
            ax.annotate(str(r.gene_symbol), (r.RA_selection_frequency, r.NSCLC_selection_frequency),
                        xytext=(4, 3), textcoords="offset points", fontsize=8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("RA selection frequency"); ax.set_ylabel("NSCLC selection frequency")
        ax.set_title("Exploratory cross-indication stable-gene overlap", loc="left", weight="bold")
        ax.grid(alpha=0.2)
        add_record(records, "Cross-indication", "Shared stable genes: RA vs NSCLC frequency",
                   save(fig, out / "cross_indication_shared_genes.png"),
                   "results/CROSS_INDICATION/shared_stable_genes.csv")

    go = numeric(read_csv(root / "results" / "CROSS_INDICATION" / "GO_BP_enrichment.csv"), ["padj", "Count"])
    if not go.empty and {"context", "Description", "padj"}.issubset(go.columns):
        q = go[go.context.astype(str).eq("RA_stable")].dropna(subset=["padj"])
        q = q.sort_values("padj").head(12).copy()
        if not q.empty:
            q["minus_log10_padj"] = -np.log10(q.padj.clip(lower=1e-300))
            q = q.sort_values("minus_log10_padj")
            fig, ax = plt.subplots(figsize=(9.5, max(4.5, len(q) * 0.32)))
            ax.barh(q.Description.astype(str), q.minus_log10_padj)
            ax.set_xlabel("−log10(BH-adjusted p-value)")
            ax.set_title("RA stable-gene GO Biological Process enrichment", loc="left", weight="bold")
            ax.grid(axis="x", alpha=0.2)
            add_record(records, "Cross-indication", "GO Biological Process enrichment",
                       save(fig, out / "cross_indication_ra_go_bp.png"),
                       "results/CROSS_INDICATION/GO_BP_enrichment.csv")


def generate_all_figures(root: Path, figdir: Path) -> list[dict[str, str]]:
    figdir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, str]] = []
    plot_nsclc_summary(root, figdir, records)
    plot_nsclc_deltas(root, figdir, records)
    plot_nsclc_fold_distributions(root, figdir, records)
    plot_nsclc_feature_stability(root, figdir, records)
    plot_nsclc_cohort_qc(root, figdir, records)
    plot_ra_performance(root, figdir, records)
    plot_ra_repeats(root, figdir, records)
    plot_ra_contribution(root, figdir, records)
    plot_ra_stability(root, figdir, records)
    plot_ra_roc(root, figdir, records)
    plot_ra_feature_qc(root, figdir, records)
    plot_ra_warnings(root, figdir, records)
    plot_cross_indication(root, figdir, records)
    manifest = pd.DataFrame(records, columns=["section", "title", "filename", "source"])
    manifest.to_csv(figdir / "figure_manifest.csv", index=False)
    return records


def write_gallery(root: Path, records: list[dict[str, str]]) -> Path:
    reports = root / "Reports"
    reports.mkdir(parents=True, exist_ok=True)
    groups: dict[str, list[dict[str, str]]] = {}
    for rec in records:
        groups.setdefault(rec["section"], []).append(rec)
    sections = []
    for section, items in groups.items():
        cards = "\n".join(
            f'<figure><a href="figures/{r["filename"]}"><img src="figures/{r["filename"]}" alt="{r["title"]}"></a>'
            f'<figcaption><b>{r["title"]}</b><br><small>Source: {r["source"]}</small></figcaption></figure>'
            for r in items
        )
        sections.append(f"<h2>{section}</h2><div class='grid'>{cards}</div>")
    page = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>MOSAIC-3 figures</title><style>
body{font:15px/1.55 system-ui,Segoe UI,Arial;margin:0;background:#f5f7fa;color:#172b40}main{max-width:1200px;margin:auto;padding:28px}h1{font-size:32px}h2{margin-top:36px;border-bottom:1px solid #d8e0e8;padding-bottom:8px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:16px}figure{margin:0;background:white;border:1px solid #d8e0e8;border-radius:12px;padding:14px}img{width:100%;height:auto}figcaption{font-size:13px;color:#526272;margin-top:8px}small{font-size:11px}</style></head><body><main>
<h1>MOSAIC-3 | Figure gallery</h1><p>Figures are generated exclusively from the current result CSV files. Missing result tables are skipped; no values are simulated.</p>
""" + "\n".join(sections) + "</main></body></html>"
    path = reports / "MOSAIC3_figures_gallery.html"
    path.write_text(page, encoding="utf-8")
    return path

def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate MOSAIC-3 figures from existing result CSV files."
    )
    parser.add_argument(
        "--root", default=".",
        help="MOSAIC-3 project root (the directory containing results/)."
    )
    parser.add_argument(
        "--figdir", default=None,
        help="Output figure directory (default: <root>/Reports/figures)."
    )
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    if not root.exists():
        parser.error(f"Project root does not exist: {root}")

    results_dir = root / "results"
    if not results_dir.exists():
        parser.error(
            f"No 'results' directory found under {root}. "
            "Pass --root pointing to the MOSAIC-3 project directory."
        )

    figdir = (
        Path(args.figdir).expanduser().resolve()
        if args.figdir
        else root / "Reports" / "figures"
    )

    print(f"MOSAIC-3 project root: {root}")
    print(f"Figure output directory: {figdir}")
    records = generate_all_figures(root, figdir)

    if records:
        gallery = write_gallery(root, records)
        print(f"\nGenerated {len(records)} figure(s).")
        print(f"Figure manifest: {figdir / 'figure_manifest.csv'}")
        print(f"HTML gallery: {gallery}")
        for record in records:
            print(f"  [{record['section']}] {record['filename']}")
    else:
        print("\nNo figures were generated.")
        print("The script did not find result CSVs with the required columns/rows.")
        print(f"An empty figure manifest was written to: {figdir / 'figure_manifest.csv'}")
        print("Check that the benchmark CSV outputs exist under:")
        print(f"  {results_dir / 'RA'}")
        print(f"  {results_dir / 'NSCLC'}")
        print(f"  {results_dir / 'CROSS_INDICATION'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
