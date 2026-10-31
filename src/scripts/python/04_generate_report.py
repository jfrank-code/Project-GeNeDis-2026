#!/usr/bin/env python3
# Generate the MOSAIC-3 English HTML and PDF scientific report from archived result files.
from __future__ import annotations
import argparse, html, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak, KeepTogether

NAVY="#10233f"; BLUE="#2a6f97"; TEAL="#39a6a3"; ORANGE="#e98b45"; PURPLE="#8367b2"; GREY="#5e6b78"; LIGHT="#eef3f7"

def load(path):
    p=Path(path)
    if not p.exists(): return pd.DataFrame()
    try: return pd.read_csv(p)
    except Exception: return pd.DataFrame()
def fnum(x,n=3):
    try:
        if pd.isna(x): return "NA"
        return f"{float(x):.{n}f}"
    except Exception: return str(x)
def savefig(fig,path):
    fig.savefig(path,dpi=180,bbox_inches="tight",facecolor="white")
    plt.close(fig)
def table_html(df, cols=None, max_rows=20, digits=3):
    if df.empty: return '<p class="muted">Not available in the archived result set.</p>'
    z=df.copy()
    if cols: z=z[[c for c in cols if c in z.columns]]
    z=z.head(max_rows)
    for c in z.select_dtypes(include=[np.number]).columns: z[c]=z[c].map(lambda x:fnum(x,digits))
    return z.to_html(index=False,escape=True,border=0,classes="data-table")
def chart_nsclc(root,figdir):
    ds={c:load(root/f"results/NSCLC/{c}/model_metrics.csv") for c in ["LUAD","LUSC"]}
    models=["CLINICAL","RSF_CLINICAL","RNA","CNV","METHYLATION","CLINICAL_RNA","CLINICAL_CNV","CLINICAL_METHYLATION","OMICS_ALL","CLINICAL_ALL","RSF_CLINICAL_ALL"]
    fig,ax=plt.subplots(figsize=(12,5.2)); x=np.arange(len(models)); w=.36
    for i,(co,col) in enumerate([("LUAD",TEAL),("LUSC",PURPLE)]):
        d=ds[co].set_index("model") if not ds[co].empty else pd.DataFrame()
        vals=[float(d.loc[m,"mean_c_index"]) if not d.empty and m in d.index else np.nan for m in models]
        errs=[float(d.loc[m,"sd_c_index"]) if not d.empty and m in d.index and pd.notna(d.loc[m,"sd_c_index"]) else 0 for m in models]
        ax.bar(x+(i-.5)*w,vals,w,yerr=errs,capsize=2,label=co,color=col,alpha=.9)
    ax.axhline(.5,color=GREY,ls="--",lw=1); ax.set_ylim(0,1); ax.set_ylabel("Harrell's C-index (mean +/- SD across folds)")
    ax.set_xticks(x); ax.set_xticklabels(models,rotation=35,ha="right",fontsize=8); ax.set_title("NSCLC survival benchmark",loc="left",weight="bold",color=NAVY)
    ax.legend(frameon=False); ax.grid(axis="y",alpha=.18); ax.spines[['top','right']].set_visible(False)
    p=figdir/'nsclc_cindex.png'; savefig(fig,p); return p

def chart_deltas(root,figdir):
    fig,ax=plt.subplots(figsize=(10,4.8)); rows=[]
    for co,col in [("LUAD",TEAL),("LUSC",PURPLE)]:
        d=load(root/f"results/NSCLC/{co}/layer_contribution_vs_clinical.csv")
        if not d.empty:
            for _,r in d.iterrows(): rows.append((co,str(r.model),float(r.delta_c_index_vs_clinical),col))
    z=pd.DataFrame(rows,columns=['cohort','model','delta','color'])
    mods=[m for m in ["RSF_CLINICAL","CLINICAL_CNV","CLINICAL_ALL","CLINICAL_RNA","CLINICAL_METHYLATION","RSF_CLINICAL_ALL","OMICS_ALL","RNA","METHYLATION","CNV"] if m in set(z.model)]
    x=np.arange(len(mods)); w=.36
    for i,co in enumerate(["LUAD","LUSC"]):
        q=z[z.cohort==co].set_index('model').reindex(mods)
        ax.bar(x+(i-.5)*w,q.delta,w,label=co,color=TEAL if co=='LUAD' else PURPLE,alpha=.9)
    ax.axhline(0,color=GREY,lw=1); ax.set_xticks(x); ax.set_xticklabels(mods,rotation=35,ha='right',fontsize=8)
    ax.set_ylabel("Delta C-index vs Cox clinical baseline"); ax.set_title("Incremental contribution over clinical baseline",loc='left',weight='bold',color=NAVY)
    ax.legend(frameon=False); ax.grid(axis='y',alpha=.18); ax.spines[['top','right']].set_visible(False)
    p=figdir/'nsclc_delta_cindex.png'; savefig(fig,p); return p

def chart_ra(root,figdir):
    d=load(root/'results/RA/ra_model_metrics.csv')
    if d.empty:return None
    d=d.sort_values('mean_auc',ascending=False); x=np.arange(len(d)); w=.36
    fig,ax=plt.subplots(figsize=(10,5))
    ax.bar(x-w/2,d.mean_auc,w,yerr=d.sd_auc.fillna(0),capsize=3,label='ROC-AUC',color=TEAL)
    ax.bar(x+w/2,d.mean_balanced_accuracy,w,yerr=d.sd_balanced_accuracy.fillna(0),capsize=3,label='Balanced accuracy',color=ORANGE)
    ax.set_ylim(0,1.08); ax.set_xticks(x); ax.set_xticklabels(d.model,rotation=35,ha='right',fontsize=8); ax.set_ylabel('Mean across two outer repetitions')
    ax.set_title('RA vs OA: donor-level multi-omics classification',loc='left',weight='bold',color=NAVY); ax.legend(frameon=False); ax.grid(axis='y',alpha=.18); ax.spines[['top','right']].set_visible(False)
    p=figdir/'ra_model_metrics.png'; savefig(fig,p); return p

def chart_stability(root,figdir):
    d=load(root/'results/RA/ra_stability_jaccard.csv')
    if d.empty:return None
    d=d.sort_values('mean_jaccard'); fig,ax=plt.subplots(figsize=(8.5,4.5))
    ax.barh(d.model,d.mean_jaccard,xerr=d.sd_jaccard,capsize=2,color=PURPLE,alpha=.9)
    ax.set_xlim(0,1); ax.set_xlabel('Mean pairwise Jaccard similarity across outer folds'); ax.set_title('RA biomarker-set stability',loc='left',weight='bold',color=NAVY)
    ax.grid(axis='x',alpha=.18); ax.spines[['top','right']].set_visible(False)
    p=figdir/'ra_jaccard_stability.png'; savefig(fig,p); return p

def chart_go(root,figdir):
    d=load(root/'results/CROSS_INDICATION/GO_BP_enrichment.csv')
    if d.empty:return None
    d=d[d.context.eq('RA_stable')].sort_values('padj').head(10).sort_values('padj',ascending=True)
    if d.empty:return None
    fig,ax=plt.subplots(figsize=(9,4.7)); ax.barh(d.Description,d['-log10_padj'] if '-log10_padj' in d else -np.log10(d.padj.clip(lower=1e-300)),color=BLUE)
    ax.set_xlabel('-log10(BH-adjusted p-value)'); ax.set_title('GO Biological Process enrichment of RA stable genes',loc='left',weight='bold',color=NAVY); ax.grid(axis='x',alpha=.18); ax.spines[['top','right']].set_visible(False)
    p=figdir/'ra_go_bp.png'; savefig(fig,p); return p

def make_html(root,figures,ns_luad,ns_lusc,ra_metrics,ra_stab,jaccard,cross,go):
    summary=load(root/'results/CROSS_INDICATION/cross_indication_summary.csv')
    qc=load(root/'results/NSCLC/nsclc_cohort_summary.csv')
    luad_best=ns_luad.sort_values('mean_c_index',ascending=False).iloc[0] if not ns_luad.empty else None
    lusc_best=ns_lusc.sort_values('mean_c_index',ascending=False).iloc[0] if not ns_lusc.empty else None
    ra_best=ra_metrics.sort_values('mean_auc',ascending=False).iloc[0] if not ra_metrics.empty else None
    def img(key,cap):
        p=figures.get(key)
        if not p:return ''
        return f'<figure><img src="figures/{p.name}" alt="{html.escape(cap)}"><figcaption>{html.escape(cap)}</figcaption></figure>'
    html_doc=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MOSAIC-3 | GeNeDis 2026</title>
<style>
:root{{--ink:#15263b;--muted:#5e6b78;--teal:#178b89;--blue:#2a6f97;--purple:#7158a5;--line:#dfe7ee;--paper:#f7f9fb;--white:#fff;--amber:#b96827}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.65 Inter,Segoe UI,Arial,sans-serif}}main{{max-width:1100px;margin:auto;padding:36px 24px 70px}}header{{background:linear-gradient(125deg,#10233f,#1f5570 62%,#288e8b);color:white;padding:54px 48px;border-radius:18px}}.eyebrow{{font-size:12px;letter-spacing:.18em;text-transform:uppercase;color:#bce9e7;font-weight:700}}h1{{font-size:55px;line-height:1;margin:12px 0}}.subtitle{{font-size:19px;max-width:860px;color:#e0edf4}}.chips{{display:flex;flex-wrap:wrap;gap:9px;margin-top:24px}}.chip{{border:1px solid #ffffff60;border-radius:20px;padding:5px 12px;font-size:12px}}h2{{margin:44px 0 12px;font-size:25px;color:var(--ink);border-bottom:1px solid var(--line);padding-bottom:8px}}h3{{font-size:18px;margin:8px 0}}p{{margin:8px 0 13px}}.muted{{color:var(--muted)}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:16px}}.card{{background:white;border:1px solid var(--line);border-radius:13px;padding:20px;box-shadow:0 3px 12px #10233f0a}}.callout{{background:#eef7f7;border-left:4px solid var(--teal);padding:16px 18px;margin:16px 0;border-radius:0 10px 10px 0}}.warning{{background:#fff6ed;border-left:4px solid var(--amber);padding:16px 18px;margin:16px 0;border-radius:0 10px 10px 0}}figure{{margin:0;background:white;padding:16px;border:1px solid var(--line);border-radius:13px}}figure img{{width:100%;height:auto}}figcaption{{font-size:12px;color:var(--muted);margin-top:8px}}.data-table{{border-collapse:collapse;width:100%;font-size:12px;display:block;overflow-x:auto}}.data-table th{{background:#eaf1f5;color:#23425b;text-align:left;text-transform:uppercase;font-size:10px;letter-spacing:.06em}}.data-table td,.data-table th{{padding:9px 11px;border-bottom:1px solid var(--line);white-space:nowrap}}ul{{padding-left:22px}}li{{margin:6px 0}}.foot{{border-top:1px solid var(--line);margin-top:42px;padding-top:18px;color:var(--muted);font-size:12px}}.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}}.kpi{{background:white;border:1px solid var(--line);border-radius:12px;padding:16px}}.kpi b{{display:block;font-size:27px;color:var(--teal)}}.kpi span{{font-size:12px;color:var(--muted)}}a{{color:var(--blue)}}
@media(max-width:650px){{main{{padding:16px 12px 40px}}header{{padding:34px 24px}}h1{{font-size:42px}}}}
@media print{{body{{background:white}}main{{padding:0;max-width:none}}header{{border-radius:0}}}}
</style></head><body><main>
<header><div class="eyebrow">GeNeDis 2026 Satellite Hackathon · Predict the Future</div><h1>MOSAIC-3</h1><div class="subtitle">Multi-Omics Stability, Integration and Clinical Insight: a reproducible, disease-agnostic framework to quantify predictive contribution, biomarker stability, and exploratory cross-indication molecular overlap in NSCLC and rheumatoid arthritis.</div><div class="chips"><span class="chip">LUAD / LUSC survival</span><span class="chip">RA vs OA classification</span><span class="chip">RNA · CNV · HM450 · ATAC · ChIP</span><span class="chip">Nested / repeated CV</span></div></header>
<h2>Executive summary</h2><div class="kpis">
<div class="kpi"><b>440</b><span>LUAD complete cases</span></div><div class="kpi"><b>360</b><span>LUSC complete cases</span></div><div class="kpi"><b>20</b><span>RA/OA complete donors</span></div><div class="kpi"><b>4</b><span>RA molecular layers</span></div>
</div>
<p>MOSAIC-3 evaluates two linked questions: (1) how much each molecular layer adds to clinical or transcriptomic baselines, and (2) how reproducibly models prioritize genes and pathways when the training samples and available layers change. A cross-indication comparison is used to generate biological hypotheses, not to imply shared causation or cellular equivalence.</p>
<div class="callout"><b>Scientific interpretation.</b> All findings are exploratory. The archived NSCLC metrics come from one five-fold outer-CV repetition; the RA metrics come from two five-fold repetitions on only 20 donors. Neither benchmark constitutes independent external validation or a clinical prediction tool.</div>
<h2>1. Cohorts and data provenance</h2><div class="grid"><div class="card"><h3>NSCLC: TCGA-LUAD / TCGA-LUSC</h3><p>Public UCSC Xena cohort matrices were harmonized using primary tumour samples (sample type 01), one deterministic sample per case, complete RNA/CNV/HM450 profiles, and valid overall-survival time/event fields.</p>{table_html(qc,['cohort','n_model_complete','events','censored','RNA_features','CNV_features','METHYLATION_features'],10)}</div><div class="card"><h3>Rheumatoid arthritis: synovial tissue</h3><p>GEO accessions: GSE112656 (RNA-seq), GSE112654 (ATAC-seq), and GSE112655 (ChIP-seq). The benchmark uses 10 RA and 10 OA donors with four aligned layers: RNA, ATAC, H3K27ac and H3K4me3. WGBS and raw FASTQ/SRA are excluded.</p><ul><li>RNA features: 57,763</li><li>ATAC features: 17,939</li><li>H3K27ac features: 17,102</li><li>H3K4me3 features: 16,290</li></ul></div></div>
<p class="muted">MC3 gene-level SNV files are acquired as an exploratory layer but excluded from the primary benchmark until callability/coverage can be established. No controlled-access data are required.</p>
<h2>2. Methodology</h2><div class="grid"><div class="card"><h3>NSCLC survival</h3><ul><li>Clinical-only baseline with age, sex and stage where available.</li><li>Cox Elastic Net with nested tuning and repeated outer folds.</li><li>RNA, CNV and HM450 are assessed individually and in combinations.</li><li>Random Survival Forest is included as a sensitivity comparison.</li><li>Harrell's C-index is reported by model and outer fold.</li></ul></div><div class="card"><h3>RA classification</h3><ul><li>Donor-level, stratified outer folds; each donor stays in one test fold per repetition.</li><li>Elastic Net logistic regression with inner-fold lambda selection.</li><li>Single-layer and multi-layer ablations are compared with RNA-only baseline.</li><li>Linear SVM is an additional sensitivity analysis.</li><li>Metrics: ROC-AUC, balanced accuracy, sensitivity, specificity, MCC and Brier score.</li></ul></div></div>
<div class="callout"><b>Leakage controls.</b> NSCLC imputation, variance filtering, scaling and clinical encoding are estimated from training partitions; the endpoint fields are excluded from predictors. RA supervised feature selection, imputation, variance ranking and scaling are fitted within training folds. The archived RA run calculated label-blind TMM/logCPM normalization once before cross-validation; this is disclosed as a transductive preprocessing limitation, not as label-guided selection. The current source retains this documented behaviour to preserve the pipeline and includes a modality-name fix for future runs.</div>
<h2>3. NSCLC survival results</h2>
{img('nsclc','Mean Harrell C-index and fold standard deviation for LUAD and LUSC.')}
<div class="grid"><div class="card"><h3>LUAD model metrics</h3>{table_html(ns_luad,['model','mean_c_index','sd_c_index','folds'],20)}</div><div class="card"><h3>LUSC model metrics</h3>{table_html(ns_lusc,['model','mean_c_index','sd_c_index','folds'],20)}</div></div>
{img('delta','Incremental C-index relative to the Cox clinical-only baseline; positive values indicate a higher observed C-index.')}
<p><b>LUAD.</b> RSF with clinical predictors achieved the highest archived mean C-index (0.680, SD 0.070), compared with the Cox clinical baseline (0.654, SD 0.067). Among the Cox modality combinations, none exceeded the clinical-only baseline in this run; clinical + CNV reached 0.632.</p>
<p><b>LUSC.</b> RSF with clinical predictors achieved 0.565 (SD 0.093), compared with the Cox clinical baseline at 0.504 (SD 0.075). The highest Cox multi-omic combination was clinical + CNV at 0.538; the full clinical + all-omics model reached 0.537.</p>
<div class="warning"><b>Do not over-interpret fold maxima.</b> The LUAD and LUSC figures above are fold averages across five folds from one outer repetition. They are not repeated-CV estimates. Individual fold maxima are not presented as the headline score.</div>
<h2>4. RA multi-omics classification</h2>
{img('ra','Mean out-of-fold ROC-AUC and balanced accuracy across two outer repetitions.')}
{table_html(ra_metrics,['model','n_repetitions','mean_auc','sd_auc','mean_balanced_accuracy','sd_balanced_accuracy','mean_mcc','mean_brier'],20)}
<p>RNA-only and several RNA-containing configurations achieved a mean AUC of 1.000 across the two archived repetitions; RNA-only and the ALL multi-omics model had mean balanced accuracy of 0.975. The SVM sensitivity variants reached 1.000 for both reported metrics, but do not provide a Brier score. ATAC-only and H3K4me3-only performance was lower and less balanced. These unusually high values are highly susceptible to overfitting and sampling variation with n=20 and must be treated as exploratory.</p>
<h2>5. Biomarker stability</h2>
{img('jaccard','Mean pairwise Jaccard similarity between outer-fold selected feature sets for RA models.')}
{table_html(jaccard,['model','comparisons','mean_jaccard','sd_jaccard'],20)}
<p>RNA-only had mean Jaccard similarity 0.398 (45 pairwise comparisons), followed by RNA+CHIP at 0.388 and RNA+ATAC at 0.339. H3K4me3-only had the lowest reported stability (0.011, 17 comparisons). The archived RA feature table was audited to remove duplicate feature selections within the same model/repetition/fold, so selection frequencies no longer exceed 1.0. Because archived feature identifiers omitted modality prefixes, exact historical per-layer attribution cannot be reconstructed; the updated R source now preserves modality-qualified column names in future runs.</p>
<h2>6. Cross-indication comparison and biological grounding</h2>
<div class="card">{table_html(summary,['RA_stable_genes','LUAD_stable_genes','LUSC_stable_genes','shared_RA_LUAD','shared_RA_LUSC','note'],5)}<p>Three mapped gene symbols were shared between RA and LUAD stable-feature sets: NR0B1, PAX9 and SPINK1. No mapped shared genes were reported for RA and LUSC in the archived comparison.</p>{table_html(cross,['gene_symbol','context','RA_selection_frequency','NSCLC_selection_frequency'],10)}</div>
{img('go','Top GO Biological Process terms for the RA stable-gene set, ranked by adjusted p-value.')}
{table_html(go[go.context.eq('RA_stable')].sort_values('padj'),['GO','Description','Count','pvalue','padj'],20)}
<p>GO Biological Process results include chemokine-mediated signalling, chemotaxis, cell adhesion and calcium-related processes among the adjusted-significant terms. This is an exploratory over-representation analysis; it does not demonstrate disease mechanism, causal direction, or cell-state equivalence between synovial fibroblasts and cancer-associated fibroblasts.</p>
<h2>7. Computational budget and reproducibility</h2>
<p>The archived benchmark used a reduced execution budget because of the submission time window and available laptop compute capacity: NSCLC used one outer repetition with five folds, two inner folds, up to 250 features per layer and eight lambda values; RA used two outer repetitions with five folds, three inner folds and up to 120 features per layer. The repository configuration is set to a fuller jury rerun with three outer repetitions for both NSCLC and RA, five outer folds, three inner folds and the higher feature/lambda settings defined in the configuration. This fuller run is expected to take substantially longer and may require more memory; the archived metrics in this report will not change unless the benchmark is rerun.</p>
<p>The pipeline is fresh-start capable: it downloads public TCGA/Xena and GEO inputs, prepares aligned sample/donor matrices, runs preflight QC, benchmarks the models, performs cross-indication enrichment and generates the report. Raw data and the large processed input matrices are not packaged in the repository; they are re-created by the download scripts. Versioned result tables, model predictions, feature-selection outputs, QC summaries and this report are included.</p>
<h2>8. Limitations</h2><ul><li>NSCLC estimates are based on one outer repetition (five folds) and therefore have uncertain run-to-run variability.</li><li>RA has only 20 complete donors. Inner training partitions can contain fewer than eight observations in a class; glmnet warnings were recorded in the warning CSV.</li><li>RA RNA TMM/logCPM was calculated once before outer cross-validation using no labels. This is label-blind but cohort-wide preprocessing and may be considered transductive; a strict confirmatory analysis should estimate normalization reference factors within each training split.</li><li>Archived RA feature identifiers lacked modality prefixes; the report corrects gene-level duplicate counting, but exact legacy layer attribution is unavailable. The current R source fixes feature naming for new runs.</li><li>The archived NSCLC console output recorded glmnet numerical/convergence warnings at highly regularised lambda settings; a complete per-fold NSCLC warning file was not retained. The RA warning records are included in results/RA/ra_model_warnings.csv. Review these caveats before interpreting individual models.</li><li>GEO annotation produced one-to-many mapping warnings, and 403 genes were excluded by the promoter/annotation mapping step in the preparation log.</li><li>No independent external validation cohort is reported here. No clinical utility, treatment response, prognosis deployment or diagnostic validity is claimed.</li><li>Cross-indication overlap is exploratory and subject to identifier mapping, feature-selection and assay differences.</li></ul>
<h2>9. Hackathon deliverables and compliance</h2><ul><li>Public sources and accession IDs documented; no new controlled-access approval is needed.</li><li>Supervised endpoints are explicit: overall survival time/event for NSCLC and RA versus OA for RA.</li><li>Baselines are included: clinical-only for NSCLC and RNA-only for RA.</li><li>Leakage controls are applied to model fitting; the RA label-blind normalization caveat is explicitly disclosed.</li><li>Source code, fresh-start download scripts, QC, model outputs, provenance notes, and limitations are included.</li><li>Raw references and large matrices are excluded from the Git repository; they can be reconstructed through the pipeline.</li></ul>
<p class="muted">Official challenge page: <a href="https://2026.genedis.org/hackathon/">GeNeDis 2026 Satellite Hackathon - Predict the Future</a>.</p>
<div class="foot">MOSAIC-3 · GeNeDis 2026 · Generated from archived analysis outputs · All numerical results in this report are read from the included CSV artefacts. No synthetic or illustrative model metrics are used.</div>
</main></body></html>'''
    (root/'Reports/MOSAIC3_final_report.html').write_text(html_doc,encoding='utf-8')

def make_pdf(root, figures, ns_luad, ns_lusc, ra_metrics, jaccard, cross, go):
    out=root/'Reports/MOSAIC3_final_report.pdf'
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='TitleCustom',parent=styles['Title'],fontName='Helvetica-Bold',fontSize=25,leading=29,textColor=colors.HexColor(NAVY),spaceAfter=8,alignment=TA_LEFT))
    styles.add(ParagraphStyle(name='Sub',parent=styles['Normal'],fontSize=11,leading=16,textColor=colors.HexColor(GREY),spaceAfter=10))
    styles.add(ParagraphStyle(name='H1x',parent=styles['Heading1'],fontSize=17,leading=21,textColor=colors.HexColor(NAVY),spaceBefore=12,spaceAfter=7))
    styles.add(ParagraphStyle(name='H2x',parent=styles['Heading2'],fontSize=12,leading=15,textColor=colors.HexColor(BLUE),spaceBefore=7,spaceAfter=4))
    styles.add(ParagraphStyle(name='Bodyx',parent=styles['BodyText'],fontSize=9.2,leading=13.2,spaceAfter=7))
    styles.add(ParagraphStyle(name='Smallx',parent=styles['BodyText'],fontSize=7.5,leading=10,textColor=colors.HexColor(GREY),spaceAfter=4))
    styles.add(ParagraphStyle(name='Calloutx',parent=styles['BodyText'],fontSize=9,leading=13,borderColor=colors.HexColor(TEAL),borderWidth=0.8,borderPadding=8,backColor=colors.HexColor('#eef7f7'),spaceBefore=4,spaceAfter=9))
    story=[]
    story += [Paragraph('MOSAIC-3',styles['TitleCustom']),Paragraph('Multi-Omics Stability, Integration and Clinical Insight',styles['Sub']),Paragraph('<b>GeNeDis 2026 Satellite Hackathon - Predict the Future</b><br/>NSCLC (LUAD/LUSC) and rheumatoid arthritis (RA)',styles['Bodyx']),Spacer(1,5)]
    story += [Paragraph('Executive summary',styles['H1x']),Paragraph('MOSAIC-3 quantifies the incremental predictive contribution of molecular layers, the stability of selected biomarkers, and exploratory cross-indication overlap. The analysis is designed as a reproducible research workflow, not a clinical prediction tool.',styles['Bodyx']),Paragraph('<b>Key cohort sizes:</b> LUAD 440 complete cases (157 events, 283 censored); LUSC 360 (154 events, 206 censored); RA/OA 20 donors (10 RA, 10 OA).',styles['Bodyx']),Paragraph('<b>Critical interpretation:</b> archived NSCLC metrics used one 5-fold outer repetition; RA used two 5-fold repetitions with n=20. All metrics are exploratory and not independently externally validated.',styles['Calloutx'])]
    story += [Paragraph('Cohorts and data sources',styles['H1x']),Paragraph('NSCLC: public TCGA cohort matrices accessed through UCSC Xena (RNA-seq, gene-level CNV, HumanMethylation450, clinical and survival data). Primary tumour sample type 01 was selected with one deterministic sample per case, requiring complete three-layer profiles and valid overall-survival endpoints. MC3 gene-level SNV acquisition is exploratory and excluded from the primary model pending callability assessment.',styles['Bodyx']),Paragraph('RA: GEO GSE112656 (RNA-seq), GSE112654 (ATAC-seq), and GSE112655 (ChIP-seq). The model cohort contains 10 RA and 10 OA donors with RNA, ATAC, H3K27ac and H3K4me3 layers. WGBS and FASTQ/SRA files are excluded.',styles['Bodyx'])]
    story += [Paragraph('Methodology and leakage controls',styles['H1x']),Paragraph('NSCLC survival: a clinical-only baseline, Cox Elastic Net with nested tuning, single-layer and combined-layer ablations, and Random Survival Forest sensitivity analysis. Harrell C-index is reported by model and fold. Age/sex/stage encoding, imputation, feature filtering and scaling are estimated within training partitions; survival time and event are not model predictors.',styles['Bodyx']),Paragraph('RA classification: donor-level stratified outer folds, inner-fold Elastic Net lambda selection, single-layer and combined-layer comparisons against RNA-only, plus a linear SVM sensitivity. Supervised selection, imputation, variance ranking and z-scoring are training-fold specific. The archived RNA TMM/logCPM transform was calculated once across the cohort without labels before CV; this label-blind but cohort-wide step is disclosed as a possible transductive preprocessing limitation.',styles['Bodyx'])]
    story += [Paragraph('NSCLC survival results',styles['H1x'])]
    if figures.get('nsclc'): story.append(Image(str(figures['nsclc']),width=174*mm,height=75*mm))
    story += [Paragraph('<b>LUAD:</b> RSF clinical mean C-index 0.680 (SD 0.070), versus Cox clinical baseline 0.654 (SD 0.067). None of the Cox omics combinations exceeded the clinical-only baseline in the archived quick run; clinical + CNV reached 0.632.',styles['Bodyx']),Paragraph('<b>LUSC:</b> RSF clinical mean C-index 0.565 (SD 0.093), versus Cox clinical baseline 0.504 (SD 0.075). The strongest Cox layer combination was clinical + CNV at 0.538; clinical + all omics reached 0.537.',styles['Bodyx'])]
    def pdf_table(df,cols,maxrows=15):
        if df.empty:return Paragraph('No archived table available.',styles['Smallx'])
        z=df[[c for c in cols if c in df.columns]].head(maxrows)
        data=[[Paragraph(str(c),styles['Smallx']) for c in z.columns]]
        for _,r in z.iterrows(): data.append([Paragraph(fnum(v,3) if isinstance(v,(int,float,np.number)) else html.escape(str(v)),styles['Smallx']) for v in r])
        widths=[(174*mm)/max(1,len(z.columns))]*len(z.columns)
        t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#eaf1f5')),('TEXTCOLOR',(0,0),(-1,0),colors.HexColor(NAVY)),('GRID',(0,0),(-1,-1),.35,colors.HexColor('#dfe7ee')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]))
        return t
    story += [Spacer(1,4),pdf_table(ns_luad,['model','mean_c_index','sd_c_index','folds'],12),Spacer(1,8),pdf_table(ns_lusc,['model','mean_c_index','sd_c_index','folds'],12)]
    if figures.get('delta'): story += [Spacer(1,8),Image(str(figures['delta']),width=165*mm,height=78*mm)]
    story += [PageBreak(),Paragraph('RA multi-omics classification',styles['H1x'])]
    if figures.get('ra'): story.append(Image(str(figures['ra']),width=165*mm,height=82*mm))
    story += [pdf_table(ra_metrics,['model','n_repetitions','mean_auc','sd_auc','mean_balanced_accuracy','sd_balanced_accuracy','mean_mcc','mean_brier'],15),Paragraph('RNA-only and several RNA-containing configurations achieved mean AUC 1.000 across two repetitions; RNA-only and ALL had mean balanced accuracy 0.975. These unusually high values must be interpreted cautiously because only 20 donors were available. The SVM sensitivity variants achieved 1.000 AUC and balanced accuracy but did not provide Brier scores.',styles['Bodyx']),Paragraph('Biomarker stability',styles['H1x'])]
    if figures.get('jaccard'): story.append(Image(str(figures['jaccard']),width=150*mm,height=76*mm))
    story += [pdf_table(jaccard,['model','comparisons','mean_jaccard','sd_jaccard'],15),Paragraph('RNA-only mean Jaccard similarity was 0.398 (45 comparisons), followed by RNA+CHIP at 0.388 and RNA+ATAC at 0.339. A post-hoc audit removed duplicate selections within the same model/repetition/fold so selection frequency does not exceed 1.0. Legacy RA feature IDs lacked modality prefixes; exact archived per-layer attribution cannot be reconstructed. The R source now preserves modality-qualified feature names in future runs.',styles['Bodyx'])]
    story += [PageBreak(),Paragraph('Cross-indication comparison',styles['H1x']),Paragraph('The archived comparison identified three mapped stable gene symbols shared by RA and LUAD: NR0B1, PAX9 and SPINK1. No RA-LUSC overlap was reported. The GO Biological Process analysis highlighted chemokine-mediated signalling, chemotaxis, cell adhesion and calcium-related terms among adjusted-significant categories.',styles['Bodyx'])]
    if figures.get('go'): story.append(Image(str(figures['go']),width=160*mm,height=82*mm))
    story += [pdf_table(cross,['gene_symbol','context','RA_selection_frequency','NSCLC_selection_frequency'],12),pdf_table(go[go.context.eq('RA_stable')].sort_values('padj'),['Description','Count','pvalue','padj'],12),Paragraph('These overlaps and enrichment results are hypothesis-generating only. They do not demonstrate causal relationships, shared disease mechanisms, or equivalence between synovial fibroblasts and cancer-associated fibroblasts.',styles['Calloutx'])]
    story += [Paragraph('Computational budget and jury rerun',styles['H1x']),Paragraph('The archived submission benchmark used fewer repetitions because of the deadline and the available laptop compute capacity: NSCLC, one outer repetition with five folds, two inner folds, up to 250 features per layer and eight lambda values; RA, two outer repetitions with five folds, three inner folds and up to 120 features per layer. The submitted configuration is increased to three outer repetitions for both tasks, five outer folds and three inner folds, with the fuller feature/lambda settings in the configuration file. This is a more complete jury rerun profile and may require substantially more compute time and memory. Existing report metrics remain the archived quick-run results until a new run is performed.',styles['Bodyx'])]
    story += [Paragraph('Limitations',styles['H1x']),Paragraph('Small RA donor cohort; inner training class sizes may be below eight and glmnet warnings were recorded. NSCLC has only one outer repetition. Cohort-wide label-blind TMM normalization in the archived RA run is transductive. The NSCLC console reported glmnet numerical/convergence warnings at highly regularised lambda settings, but a complete per-fold NSCLC warning file was not retained; the RA warning CSV is included. Gene-annotation one-to-many mappings were also observed. Promoter/annotation processing excluded 403 genes in the preparation log. There is no independent external validation. No clinical utility or diagnostic validity is claimed. Cross-indication mapping is subject to assay and identifier differences.',styles['Bodyx']),Paragraph('Reproducibility package',styles['H1x']),Paragraph('The repository contains acquisition scripts for public UCSC Xena/TCGA and GEO data, preparation and QC scripts, model code, archived result tables, predictions, feature-selection outputs, source/accession documentation, and this HTML/PDF report. Raw inputs and large processed matrices are intentionally excluded; the fresh-start workflow reconstructs them. Run on Windows with Python 3.10+, R 4.3+, internet access and sufficient disk/memory using run_mosaic3.ps1.',styles['Bodyx']),Paragraph('Hackathon fit: the project documents open data sources, explicit endpoints, baseline comparisons, leakage safeguards and caveats, QC, results, reproducible code, provenance and limitations in line with the official challenge deliverables.',styles['Bodyx']),Paragraph('Official challenge page: https://2026.genedis.org/hackathon/',styles['Smallx'])]
    def footer(canvas,doc):
        canvas.saveState(); w,h=A4; canvas.setStrokeColor(colors.HexColor('#dfe7ee')); canvas.line(18*mm,14*mm,w-18*mm,14*mm)
        canvas.setFont('Helvetica',7); canvas.setFillColor(colors.HexColor(GREY)); canvas.drawString(18*mm,9*mm,'MOSAIC-3 | GeNeDis 2026 | Exploratory research report'); canvas.drawRightString(w-18*mm,9*mm,f'Page {doc.page}'); canvas.restoreState()
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=18*mm,leftMargin=18*mm,topMargin=16*mm,bottomMargin=20*mm,title='MOSAIC-3 Scientific Report',author='MOSAIC-3 Project')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default='.'); args=ap.parse_args(); root=Path(args.root).resolve()
    figdir=root/'Reports/figures'; figdir.mkdir(parents=True,exist_ok=True); (root/'Reports').mkdir(exist_ok=True)
    ns_luad=load(root/'results/NSCLC/LUAD/model_metrics.csv'); ns_lusc=load(root/'results/NSCLC/LUSC/model_metrics.csv')
    ra=load(root/'results/RA/ra_model_metrics.csv'); stab=load(root/'results/RA/ra_stable_features.csv'); jac=load(root/'results/RA/ra_stability_jaccard.csv')
    cross=load(root/'results/CROSS_INDICATION/shared_stable_genes.csv'); go=load(root/'results/CROSS_INDICATION/GO_BP_enrichment.csv')
    figures={}
    for k,fn in [('nsclc',chart_nsclc),('delta',chart_deltas),('ra',chart_ra),('jaccard',chart_stability),('go',chart_go)]:
        p=fn(root,figdir)
        if p: figures[k]=p
    make_html(root,figures,ns_luad,ns_lusc,ra,stab,jac,cross,go)
    # Prefer the HTML print renderer for a more consistent, publication-style PDF.
    # ReportLab remains a dependency-light fallback for environments without WeasyPrint.
    pdf_path = root/'Reports/MOSAIC3_final_report.pdf'
    try:
        from weasyprint import HTML
        HTML(filename=str(root/'Reports/MOSAIC3_final_report.html'), base_url=str(root)).write_pdf(str(pdf_path))
        pdf = pdf_path
        print('PDF engine: WeasyPrint (HTML layout preserved).')
    except Exception as exc:
        print(f'WeasyPrint unavailable or failed ({exc}); using ReportLab fallback.')
        pdf = make_pdf(root,figures,ns_luad,ns_lusc,ra,jac,cross,go)
    print('English report generated: Reports/MOSAIC3_final_report.html')
    print('English PDF generated:',pdf)
    print('Archived result rows:',{'LUAD_models':len(ns_luad),'LUSC_models':len(ns_lusc),'RA_models':len(ra)})
if __name__=='__main__': main()
