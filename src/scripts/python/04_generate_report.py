#!/usr/bin/env python3
"""Create the MOSAIC-3 figure gallery and dynamic HTML/PDF report.

Run after the NSCLC, RA and cross-indication analysis scripts:
    python src/scripts/python/04_generate_report.py --root .

All charts and reported metrics are read from the current result CSV artefacts.
No archived headline values are hard-coded into this report.
"""
from __future__ import annotations

import argparse
import html
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from mosaic3_figures import generate_all_figures, write_gallery


def load(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception as exc:
        warnings.warn(f"Could not read {path}: {exc}")
        return pd.DataFrame()


def fmt(value, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.{digits}f}"
    except Exception:
        pass
    return str(value)


def table_html(df: pd.DataFrame, columns: list[str] | None = None, max_rows: int = 25) -> str:
    if df.empty:
        return '<p class="muted">No result table is available yet. Run the corresponding analysis first.</p>'
    z = df.copy()
    if columns:
        use = [c for c in columns if c in z.columns]
        if use:
            z = z[use]
    z = z.head(max_rows)
    for col in z.columns:
        if pd.api.types.is_numeric_dtype(z[col]):
            z[col] = z[col].map(fmt)
    return z.to_html(index=False, escape=True, border=0, classes="data-table")


def figure_sections(records: list[dict[str, str]]) -> str:
    groups: dict[str, list[dict[str, str]]] = {}
    for rec in records:
        groups.setdefault(rec["section"], []).append(rec)
    parts = []
    for section, items in groups.items():
        cards = []
        for rec in items:
            title = html.escape(rec["title"])
            filename = html.escape(rec["filename"])
            source = html.escape(rec["source"])
            cards.append(
                f'<figure><a href="figures/{filename}"><img src="figures/{filename}" alt="{title}"></a>'
                f'<figcaption><b>{title}</b><br><small>Source: {source}</small></figcaption></figure>'
            )
        parts.append(f'<h3>{html.escape(section)}</h3><div class="grid">{"".join(cards)}</div>')
    if not parts:
        return '<p class="muted">No figures were generated because the expected result tables are missing or empty.</p>'
    return "\n".join(parts)


def get_ra_manifest(root: Path) -> pd.DataFrame:
    return load(root / "results/RA/ra_qc_model_manifest.csv")


def get_cancer_repeats(root: Path) -> dict[str, int]:
    result = {}
    for cohort in ("LUAD", "LUSC"):
        d = load(root / "results" / "NSCLC" / cohort / "metrics_by_fold.csv")
        if not d.empty and "repetition" in d.columns:
            values = pd.to_numeric(d.repetition, errors="coerce").dropna()
            result[cohort] = int(values.nunique())
    return result


def make_html(root: Path, records: list[dict[str, str]]) -> Path:
    luad = load(root / "results/NSCLC/LUAD/model_metrics.csv")
    lusc = load(root / "results/NSCLC/LUSC/model_metrics.csv")
    ra_metrics = load(root / "results/RA/ra_model_metrics.csv")
    ra_stable = load(root / "results/RA/ra_stable_features.csv")
    cross = load(root / "results/CROSS_INDICATION/shared_stable_genes.csv")
    go = load(root / "results/CROSS_INDICATION/GO_BP_enrichment.csv")
    cohort_qc = load(root / "results/NSCLC/nsclc_cohort_summary.csv")
    ra_manifest = get_ra_manifest(root)
    cancer_reps = get_cancer_repeats(root)

    def scalar(df: pd.DataFrame, col: str, default="Not available"):
        if df.empty or col not in df.columns:
            return default
        return fmt(df.iloc[0][col], 0)

    ra_donors = scalar(ra_manifest, "donors")
    ra_n = scalar(ra_manifest, "outer_repeats")
    ra_k = scalar(ra_manifest, "outer_folds")
    ra_inner = scalar(ra_manifest, "inner_folds")

    qc_by_cohort = {str(r.cohort): r for r in cohort_qc.itertuples()} if not cohort_qc.empty and "cohort" in cohort_qc else {}
    def qc_count(cohort: str, field: str) -> str:
        r = qc_by_cohort.get(cohort)
        return fmt(getattr(r, field)) if r is not None and hasattr(r, field) else "NA"

    feature_tables = []
    for cohort in ("LUAD", "LUSC"):
        d = load(root / f"results/NSCLC/{cohort}/feature_stability.csv")
        if not d.empty:
            model = next((m for m in ("CLINICAL_ALL", "OMICS_ALL") if "model" in d and (d.model.astype(str) == m).any()), None)
            q = d[d.model.astype(str) == model] if model else d
            q = q.sort_values("selection_frequency", ascending=False) if "selection_frequency" in q else q
            feature_tables.append(f'<div class="card"><h4>{cohort} stable features' + (f' ({html.escape(model)})' if model else '') + '</h4>' + table_html(q, ["model", "feature_id", "selected_count", "selection_frequency"], 15) + '</div>')
    if not feature_tables:
        feature_tables_html = '<p class="muted">NSCLC feature-stability tables are not available.</p>'
    else:
        feature_tables_html = '<div class="grid">' + ''.join(feature_tables) + '</div>'

    ra_feature_html = table_html(
        ra_stable.sort_values("selection_frequency", ascending=False) if not ra_stable.empty and "selection_frequency" in ra_stable else ra_stable,
        ["model", "feature_id", "modality", "selected_count", "completed_outer_folds", "selection_frequency"], 30
    )
    go_ra = go[go.context.astype(str).eq("RA_stable")] if not go.empty and "context" in go.columns else pd.DataFrame()

    html_doc = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>MOSAIC-3 | Current analysis report</title>
<style>
:root{{--ink:#172b40;--muted:#5e6b78;--teal:#147d7a;--blue:#2a6f97;--line:#dce5ec;--paper:#f5f7fa;--amber:#9b5b1f}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.65 Inter,Segoe UI,Arial,sans-serif}}
main{{max-width:1180px;margin:auto;padding:32px 24px 64px}}header{{background:linear-gradient(125deg,#10233f,#205a72 65%,#238b87);color:white;padding:44px;border-radius:16px}}
.eyebrow{{font-size:12px;letter-spacing:.16em;text-transform:uppercase;color:#c3e9e7;font-weight:700}}h1{{font-size:46px;line-height:1.1;margin:12px 0}}.subtitle{{font-size:18px;color:#e5eef3;max-width:900px}}
h2{{font-size:25px;margin:38px 0 12px;border-bottom:1px solid var(--line);padding-bottom:8px}}h3{{font-size:19px;margin:24px 0 12px}}h4{{margin:0 0 8px}}
.muted{{color:var(--muted)}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:15px;margin:14px 0}}
.card,figure{{background:white;border:1px solid var(--line);border-radius:12px;padding:16px;margin:0}}figure img{{width:100%;height:auto;display:block}}figcaption{{font-size:12px;color:var(--muted);margin-top:8px;overflow-wrap:anywhere}}
.data-table{{border-collapse:collapse;width:100%;font-size:12px;display:block;overflow-x:auto;background:white}}.data-table th{{background:#eaf1f5;text-align:left;color:#25435a;font-size:10px;text-transform:uppercase;letter-spacing:.05em}}.data-table td,.data-table th{{padding:8px 10px;border-bottom:1px solid var(--line);white-space:nowrap}}
.note{{background:#fff7ed;border-left:4px solid var(--amber);padding:14px 17px;margin:16px 0;border-radius:0 8px 8px 0}}.callout{{background:#eaf5f5;border-left:4px solid var(--teal);padding:14px 17px;margin:16px 0;border-radius:0 8px 8px 0}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:15px 0}}.kpi{{background:white;border:1px solid var(--line);padding:14px;border-radius:11px}}.kpi b{{display:block;font-size:26px;color:var(--teal)}}.kpi span{{font-size:12px;color:var(--muted)}}a{{color:var(--blue)}}.foot{{border-top:1px solid var(--line);margin-top:38px;padding-top:16px;color:var(--muted);font-size:12px}}
@media(max-width:650px){{main{{padding:14px}}header{{padding:28px 22px}}h1{{font-size:36px}}}}
@media print{{body{{background:white}}main{{max-width:none;padding:0}}header{{border-radius:0}}figure,.card{{break-inside:avoid}}}}
</style></head><body><main>
<header><div class="eyebrow">GeNeDis 2026 · MOSAIC-3</div><h1>Multi-Omics Stability, Integration and Clinical Insight</h1>
<div class="subtitle">Current, reproducible results for TCGA-LUAD, TCGA-LUSC and GEO synovial rheumatoid arthritis versus osteoarthritis. All figures and metric tables below are generated from the result files present at report-generation time.</div></header>
<h2>Run snapshot</h2><div class="kpis">
<div class="kpi"><b>{qc_count('LUAD','n_model_complete')}</b><span>LUAD complete cases</span></div>
<div class="kpi"><b>{qc_count('LUSC','n_model_complete')}</b><span>LUSC complete cases</span></div>
<div class="kpi"><b>{ra_donors}</b><span>RA/OA paired donors</span></div>
<div class="kpi"><b>{ra_n}</b><span>RA outer repetitions</span></div>
<div class="kpi"><b>{len(records)}</b><span>Figures generated</span></div></div>
<p class="muted">NSCLC outer repetitions in the current result files: LUAD {cancer_reps.get('LUAD','NA')}; LUSC {cancer_reps.get('LUSC','NA')}. RA validation design recorded by the benchmark: {ra_n} repetitions × {ra_k} folds, with {ra_inner} inner folds.</p>
<div class="note"><b>Interpretation limit.</b> Repeated cross-validation reuses the same participants across different partitions. Ten RA repetitions do not turn 20 donors into 200 independent donors and are not external validation. RA metrics remain exploratory and should not be described as clinically validated diagnostic performance.</div>
<h2>1. Cohorts and modelling</h2>
<div class="grid"><div class="card"><h4>NSCLC survival</h4><p>LUAD and LUSC are analysed separately. Models use a clinical baseline, Cox Elastic Net with RNA, gene-level CNV and HM450 methylation individually and in combinations, plus configured Random Survival Forest sensitivity models. Main performance metric: Harrell's C-index.</p>
{table_html(cohort_qc,["cohort","n_model_complete","events","censored","RNA_features","CNV_features","METHYLATION_features"],10)}</div>
<div class="card"><h4>RA versus OA classification</h4><p>Donor-level modelling uses RNA, ATAC, H3K27ac and H3K4me3. Elastic Net logistic regression is the primary model; linear SVM is a sensitivity analysis when enabled. Metrics include ROC-AUC, balanced accuracy, sensitivity, specificity, MCC and Brier score where probability predictions are available.</p>
{table_html(ra_manifest,["donors","RA","OA","outer_repeats","outer_folds","inner_folds","max_total_features_per_model","nlambda","svm_sensitivity_run"],5)}</div></div>
<h2>2. NSCLC results</h2><div class="grid"><div class="card"><h4>LUAD model metrics</h4>{table_html(luad,["model","mean_c_index","sd_c_index","folds"],30)}</div><div class="card"><h4>LUSC model metrics</h4>{table_html(lusc,["model","mean_c_index","sd_c_index","folds"],30)}</div></div>
{figure_sections([r for r in records if r['section']=='NSCLC'])}
<h3>NSCLC feature stability tables</h3>{feature_tables_html}
<h2>3. RA multi-omics results</h2>
{table_html(ra_metrics,["model","n_repetitions","mean_auc","sd_auc","mean_balanced_accuracy","sd_balanced_accuracy","mean_sensitivity","mean_specificity","mean_mcc","mean_brier"],30)}
{figure_sections([r for r in records if r['section']=='RA'])}
<h3>RA stable features</h3>{ra_feature_html}
<h2>4. Cross-indication results</h2>
{table_html(cross,["gene_symbol","context","RA_selection_frequency","NSCLC_selection_frequency"],30)}
{table_html(go_ra.sort_values('padj') if not go_ra.empty and 'padj' in go_ra else go_ra,["GO","Description","Count","pvalue","padj"],20)}
{figure_sections([r for r in records if r['section']=='Cross-indication'])}
<h2>5. Reproducibility and statistical safeguards</h2>
<div class="callout"><b>RA preprocessing.</b> The provided leakage-conscious R scripts retain raw RNA counts during preparation and estimate TMM/logCPM scaling within each training partition. Held-out donors are normalized relative to a fixed reference donor selected from that training partition. Expression/detection filtering, imputation, variance ranking, scaling and lambda selection are training-only. Outer-test labels are used for final evaluation only.</div>
<ul><li>For RA, the biological unit is the donor; technical files do not count as independent participants.</li><li>For NSCLC, clinical endpoints (survival time and event) are outcomes, not predictors. Preprocessing and feature selection are fitted within training partitions.</li><li>Feature selection frequency is calculated over completed model/fold runs; missing or unsuccessful folds must be inspected in the warning/QC outputs.</li><li>ROC curves are constructed from out-of-fold probability predictions and exclude SVM decision scores, which are not probabilities.</li><li>Cross-indication shared genes and GO enrichment are hypothesis-generating. They do not establish causality, identical effects or cell-state equivalence.</li><li>No independent external validation cohort is included in this workflow. No clinical utility or diagnostic validity is claimed.</li></ul>
<p><a href="MOSAIC3_figures_gallery.html">Open the standalone figure gallery</a>. The source CSV files remain under <code>results/</code>; figure provenance is listed in <code>Reports/figures/figure_manifest.csv</code>.</p>
<div class="foot">MOSAIC-3 · GeNeDis 2026 · Report generated from available local result artefacts. Missing tables are skipped rather than replaced by synthetic values.</div>
</main></body></html>'''
    out = root / "Reports" / "MOSAIC3_final_report.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html_doc, encoding="utf-8")
    return out


def make_pdf(root: Path, records: list[dict[str, str]]) -> Path | None:
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import mm
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak
    except Exception as exc:
        print(f"PDF skipped: ReportLab not available ({exc}).")
        return None
    out = root / "Reports" / "MOSAIC3_final_report.pdf"
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="MosaicHeading", parent=styles["Heading1"], textColor=colors.HexColor("#17324b"), spaceBefore=12))
    styles.add(ParagraphStyle(name="MosaicSmall", parent=styles["BodyText"], fontSize=8, leading=11))
    story = [Paragraph("MOSAIC-3", styles["Title"]),
             Paragraph("Multi-Omics Stability, Integration and Clinical Insight", styles["Heading2"]),
             Paragraph("Current report generated from local result CSV files. All findings are exploratory; no independent external validation is included.", styles["MosaicSmall"]),
             Spacer(1, 8)]
    last_section = None
    for rec in records:
        if rec["section"] != last_section:
            story.append(Paragraph(rec["section"], styles["MosaicHeading"]))
            last_section = rec["section"]
        story.append(Paragraph(html.escape(rec["title"]), styles["Heading3"]))
        p = root / "Reports" / "figures" / rec["filename"]
        try:
            im = Image(str(p))
            max_w, max_h = 170 * mm, 112 * mm
            scale = min(max_w / im.imageWidth, max_h / im.imageHeight, 1.0)
            im.drawWidth *= scale; im.drawHeight *= scale
            story.extend([im, Spacer(1, 7)])
        except Exception as exc:
            story.append(Paragraph(f"Figure unavailable: {html.escape(str(exc))}", styles["MosaicSmall"]))
    story.append(PageBreak())
    story.append(Paragraph("Interpretation and limitations", styles["MosaicHeading"]))
    story.append(Paragraph("RA repeated cross-validation reuses the same 20 donors and is not external validation. The NSCLC and RA analyses are exploratory research outputs, not clinical prediction tools. Cross-indication overlap is hypothesis-generating, not causal evidence.", styles["MosaicSmall"]))
    doc = SimpleDocTemplate(str(out), pagesize=A4, rightMargin=18 * mm, leftMargin=18 * mm,
                            topMargin=16 * mm, bottomMargin=18 * mm,
                            title="MOSAIC-3 Current Report", author="MOSAIC-3")
    doc.build(story)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="MOSAIC-3 repository root")
    parser.add_argument("--skip-pdf", action="store_true", help="Generate only HTML and PNG figures")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    reports = root / "Reports"
    figdir = reports / "figures"
    reports.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    records = generate_all_figures(root, figdir)
    gallery = write_gallery(root, records)
    report = make_html(root, records)
    print(f"Figures generated: {len(records)}")
    for rec in records:
        print(f"  [{rec['section']}] {rec['filename']} <- {rec['source']}")
    print("Figure manifest:", figdir / "figure_manifest.csv")
    print("Figure gallery:", gallery)
    print("HTML report:", report)

    if not args.skip_pdf:
        pdf = make_pdf(root, records)
        if pdf:
            print("PDF report:", pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
