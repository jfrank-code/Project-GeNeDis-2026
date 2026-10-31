#!/usr/bin/env python3
"""Prepare LUAD or LUSC Xena matrices and one primary tumour per case."""
from __future__ import annotations
import argparse, gzip, json, re
from pathlib import Path
import numpy as np
import pandas as pd

TCGA_SAMPLE_RE = re.compile(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}[A-Z]?)', re.I)
TCGA_CASE_RE = re.compile(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})', re.I)


def locate_layer(raw: Path, layer: str, cohort: str) -> Path:
    candidates = sorted((raw/layer).glob('*.tsv*'))
    if not candidates: raise FileNotFoundError(f'Downloaded file not found for {cohort}/{layer}')
    return candidates[0]


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    # pandas infers gzip from filename extension.
    return pd.read_csv(path, sep='\t', compression='infer', low_memory=False, **kwargs)


def norm_sample(x: str) -> str:
    m = TCGA_SAMPLE_RE.search(str(x).upper())
    return m.group(1) if m else ''


def case_id(x: str) -> str:
    m=TCGA_CASE_RE.search(str(x).upper())
    return m.group(1) if m else ''

def sample_type(x: str) -> str:
    # Anchor the sample-type field after the 12-character case barcode.
    # Searching for '-01' anywhere is incorrect: case barcodes can themselves
    # contain '-01' in their last four characters.
    m=re.match(r'^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-(\d{2})[A-Z]?',str(x).upper())
    return m.group(2) if m else ''


def select_primary_samples(ids: list[str]) -> dict[str,str]:
    by_case: dict[str,list[str]]={}
    for raw in ids:
        s=norm_sample(raw)
        if not s: continue
        if sample_type(s)!='01': continue
        by_case.setdefault(s[:12],[]).append(s)
    selected={}
    for c,vals in by_case.items():
        # Prefer vial A; then earliest deterministic aliquot label.
        vals=sorted(set(vals), key=lambda x:(0 if re.search(r'-01A(?:-|$)',x) else 1, x))
        selected[c]=vals[0]
    return selected


def load_xena_matrix(path: Path, layer: str) -> pd.DataFrame:
    print(f'Leyendo {layer}: {path.name}')
    # Read numeric feature columns as float32 to keep the large 450K methylation
    # matrix within reasonable RAM limits. Preserve the identifier column as text.
    header=pd.read_csv(path,sep='\t',compression='infer',nrows=0)
    if len(header.columns)<3: raise ValueError(f'Empty or non-tabular Xena matrix: {path}')
    dtype={str(c):np.float32 for c in header.columns[1:]}
    dtype[str(header.columns[0])]=str
    df=pd.read_csv(path,sep='\t',compression='infer',low_memory=False,dtype=dtype)
    if df.empty or df.shape[1]<3: raise ValueError(f'Empty or non-tabular Xena matrix: {path}')
    idcol=df.columns[0]
    firstvals=df.iloc[:min(30,len(df)),0].astype(str).tolist()
    n_sample_rows=sum(bool(norm_sample(x) or case_id(x)) for x in firstvals)
    # Xena genomic matrices normally have features in rows and sample barcodes
    # in columns. Some phenotypic-style exports put samples in rows instead.
    if n_sample_rows >= max(2, len(firstvals)//3):
        df=df.set_index(idcol).T
        df.index=df.index.astype(str)
    else:
        df=df.set_index(idcol)
    sample_cols=[c for c in df.columns if norm_sample(c)]
    # A transposed object may still expose IDs in the index after the first pass.
    if len(sample_cols)<3:
        sample_rows=[r for r in df.index.astype(str) if norm_sample(r)]
        if len(sample_rows)>=3:
            df=df.T
            sample_cols=[c for c in df.columns if norm_sample(c)]
    if len(sample_cols)<3:
        raise ValueError(f'TCGA columns were not detected in matrix {layer}. Example columns: {list(df.columns[:8])}')
    df=df.loc[:,sample_cols]
    df.index=df.index.astype(str)
    df=df.loc[~df.index.isin(['','nan','NA','ID_REF','_GENE_'])]
    df=df.apply(pd.to_numeric,errors='coerce').astype(np.float32).replace([np.inf,-np.inf],np.nan)
    df=df.groupby(level=0,sort=False).mean()
    # Remove entirely empty features before downstream variance ranking.
    df=df.loc[df.notna().any(axis=1)]
    return df


def map_matrix_to_rna_cases(mat: pd.DataFrame, rna_case_map: dict[str,str], layer: str) -> tuple[pd.DataFrame,pd.DataFrame]:
    colmap={}
    for c in mat.columns:
        s=norm_sample(c); case=s[:12] if s else ''
        # Only primary tumour sample-type 01.
        if case in rna_case_map and sample_type(s)=='01':
            colmap.setdefault(case,[]).append(s)
    chosen={}
    for case,vals in colmap.items():
        vals=sorted(set(vals),key=lambda x:(0 if re.search(r'-01A(?:-|$)',x) else 1,x))
        chosen[case]=vals[0]
    common=[case for case in rna_case_map if case in chosen]
    out=mat.loc[:,[next(c for c in mat.columns if norm_sample(c)==chosen[k]) for k in common]].copy()
    out.columns=[rna_case_map[k] for k in common]
    report=pd.DataFrame({'case_id':list(rna_case_map),'sample_id':[rna_case_map[k] for k in rna_case_map],
                         f'{layer}_present':[k in chosen for k in rna_case_map],
                         f'{layer}_source_sample':[chosen.get(k,'') for k in rna_case_map]})
    return out,report


def find_id_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        if c.lower() in {'sample','sample_id','submitter_id','bcr_patient_barcode','_patient'}: return c
    return df.columns[0]


def first_nonmissing(row: pd.Series, names: list[str]):
    lower={str(c).lower():c for c in row.index}
    for name in names:
        col=lower.get(name.lower())
        if col is not None and pd.notna(row[col]) and str(row[col]).strip().lower() not in {'','na','nan','[not available]','[not applicable]'}:
            return row[col]
    return np.nan


def prepare_clinical(raw: Path, rna_case_map: dict[str,str], proc: Path, cohort: str) -> pd.DataFrame:
    clin_path=locate_layer(raw,'CLINICAL',cohort)
    clin=read_table(clin_path)
    idcol=find_id_col(clin)
    clin['__id']=clin[idcol].astype(str)
    clin['case_id']=clin['__id'].map(case_id)
    # A separate PanCancer Atlas survival file, when available, takes precedence for endpoints.
    surv_files=list((raw/'CLINICAL').glob(f'*{cohort}_survival*.tsv*'))
    survival=None
    if surv_files:
        s=read_table(surv_files[0]); sid=find_id_col(s); s['case_id']=s[sid].astype(str).map(case_id)
        survival=s.drop_duplicates('case_id',keep='first').set_index('case_id')
    clin=clin[clin['case_id'].isin(rna_case_map)].copy()
    # Clinical Xena matrix is generally patient-level. Preserve a single row per case.
    clin=clin.drop_duplicates('case_id',keep='first').set_index('case_id')
    rows=[]
    for c,sample in rna_case_map.items():
        r=clin.loc[c] if c in clin.index else pd.Series(dtype=object)
        if isinstance(r,pd.DataFrame): r=r.iloc[0]
        sr=survival.loc[c] if survival is not None and c in survival.index else pd.Series(dtype=object)
        if isinstance(sr,pd.DataFrame): sr=sr.iloc[0]
        os_time=first_nonmissing(sr,['OS.time']) if len(sr) else np.nan
        os_event=first_nonmissing(sr,['OS']) if len(sr) else np.nan
        if pd.isna(os_time): os_time=first_nonmissing(r,['OS.time'])
        if pd.isna(os_event): os_event=first_nonmissing(r,['OS'])
        # Normalize event using only explicit OS/vital-status fields first.
        if isinstance(os_event,str):
            ev=1 if (os_event.strip().lower() in {'dead','deceased','1','yes','true'} or 'dead' in os_event.strip().lower()) else (0 if (os_event.strip().lower() in {'alive','0','no','false'} or 'alive' in os_event.strip().lower()) else np.nan)
        else:
            try:
                q=float(os_event) if pd.notna(os_event) else np.nan
                ev=int(q) if np.isfinite(q) and q in (0,1) else np.nan
            except Exception: ev=np.nan
        try: tm=float(os_time)
        except Exception: tm=np.nan
        # Fallback derived from TCGA follow-up fields; never cast days_to_death as event.
        dtd=first_nonmissing(r,['days_to_death']); dlf=first_nonmissing(r,['days_to_last_followup','days_to_last_follow_up'])
        vital=first_nonmissing(r,['vital_status'])
        if pd.isna(ev) and pd.notna(vital):
            v=str(vital).strip().lower(); ev=1 if v in {'dead','deceased'} else (0 if v=='alive' else np.nan)
        try:
            if (not np.isfinite(tm) or tm<=0) and pd.notna(dtd) and float(dtd)>0:
                tm=float(dtd); ev=1 if pd.isna(ev) else ev
            elif (not np.isfinite(tm) or tm<=0) and pd.notna(dlf) and float(dlf)>0:
                tm=float(dlf); ev=0 if pd.isna(ev) else ev
        except Exception: pass
        age=first_nonmissing(r,['age_at_initial_pathologic_diagnosis','age_at_diagnosis','age'])
        try: age=float(age)
        except Exception: age=np.nan
        if pd.notna(age) and age>200: age=age/365.25
        sex=first_nonmissing(r,['gender','sex'])
        stage=first_nonmissing(r,['pathologic_stage','ajcc_pathologic_tumor_stage','tumor_stage'])
        rows.append({'case_id':c,'sample_id':sample,'time_days':tm,'event':ev,
                     'age_years_candidate':age,'sex_candidate':sex,'stage_group':stage})
    out=pd.DataFrame(rows)
    out=out.replace([np.inf,-np.inf],np.nan)
    out.to_csv(proc/'clinical_all.csv',index=False)
    return out


def prepare_snv(raw: Path, canonical_cases: dict[str,str]) -> pd.DataFrame | None:
    files=list((raw/'SNV').glob('*.tsv*'))
    if not files: return None
    mat_raw=load_xena_matrix(files[0],'SNV_gene_level')
    mat, _ = map_matrix_to_rna_cases(mat_raw,canonical_cases,'SNV')
    # Xena gene-level mutation values encode called non-silent mutations. Keep as
    # an exploratory block only: absent data/callability cannot safely be set to WT.
    return mat


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default='.'); ap.add_argument('--cohort',choices=['LUAD','LUSC'],required=True); args=ap.parse_args(); root=Path(args.root).resolve()
    cohort=args.cohort
    raw=root/'data/raw'/cohort; proc=root/'data/processed'/cohort; proc.mkdir(parents=True,exist_ok=True)
    # RNA defines the canonical primary-tumour sample per case.
    rna_raw=load_xena_matrix(locate_layer(raw,'RNA',cohort),'RNA')
    rna_map=select_primary_samples(list(rna_raw.columns))
    if len(rna_map)<100: raise RuntimeError(f'{cohort}: has only {len(rna_map)} primary cases detected in RNA; check identifiers.')
    print(f'{cohort} RNA primary tumour: {len(rna_map)} casos.')
    canonical={case:sample for case,sample in rna_map.items()}
    rna=rna_raw.loc[:,[next(c for c in rna_raw.columns if norm_sample(c)==s) for s in canonical.values()]].copy()
    rna.columns=list(canonical.values())
    layer_mats={'RNA':rna}; align_rows=[]
    for layer in ['CNV','METHYLATION']:
        mat_raw=load_xena_matrix(locate_layer(raw,layer,cohort),layer)
        mat,rep=map_matrix_to_rna_cases(mat_raw,canonical,layer)
        layer_mats[layer]=mat
        align_rows.append(rep)
        print(f'{layer}: {mat.shape[1]} pacientes alineados, {mat.shape[0]} características.')
    # Clinical is case-level; it is explicitly aligned to the RNA-defined one-sample-per-case cohort.
    clin=prepare_clinical(raw,canonical,proc,cohort)
    # Optional somatic SNV feature block; it never shrinks the three-layer core cohort.
    snv=None; snv_error=''; snv_n_features=0
    try:
        snv=prepare_snv(raw,canonical)
    except Exception as exc:
        snv_error=str(exc)
        print(f'SNV exploratorio omitido (no bloquea el benchmark principal): {exc}')
    if snv is not None:
        # Secondary Xena gene-level mutation block. It is deliberately NOT in the
        # primary survival benchmark; callability/negative-call coverage is not inferred.
        snv.to_csv(proc/'SNV_exploratory_candidates.tsv.gz',sep='\t',compression='gzip',index_label='feature_id')
        snv_n_features=snv.shape[0]
        print(f'SNV exploratorio: {snv.shape[1]} pacientes; excluido del benchmark primario hasta verificar callability.')
    pd.DataFrame([{'status':'prepared' if snv is not None else ('unavailable_or_skipped' if not snv_error else 'parse_failed'),
                   'features':snv_n_features,'error':snv_error,
                   'interpretation':'Exploratory only; missing calls are not interpreted as wild type without callability.'}]).to_csv(proc/'SNV_exploratory_QC.csv',index=False)
    # Keep only primary tumour sample IDs and columns in common for main 3-layer models.
    core_common=set(clin.loc[clin.time_days.notna() & clin.event.notna() & (clin.time_days>0),'sample_id'])
    for layer in ['RNA','CNV','METHYLATION']: core_common &= set(layer_mats[layer].columns)
    common=[s for s in rna.columns if s in core_common]
    for layer in ['RNA','CNV','METHYLATION']:
        layer_mats[layer]=layer_mats[layer].loc[:,common]
    clin=clin[clin.sample_id.isin(common)].copy().sort_values('sample_id')
    # Write compact feature-by-sample TSVs; no binary formats required for the new subtype.
    for layer,mat in layer_mats.items():
        out=proc/f'{layer}.tsv.gz'
        mat.to_csv(out,sep='\t',compression={'method':'gzip','compresslevel':1},index_label='feature_id')
    clin.to_csv(proc/'clinical.csv',index=False)
    sample_qc=pd.DataFrame({'sample_id':common,'case_id':[s[:12] for s in common]})
    for layer,mat in layer_mats.items(): sample_qc[f'{layer}_present']=sample_qc.sample_id.isin(mat.columns)
    sample_qc.to_csv(proc/'sample_alignment.csv',index=False)
    clin[['case_id','sample_id','time_days','event']].sort_values('sample_id').to_csv(proc/'selected_sample_manifest.csv',index=False)
    pd.DataFrame([{'cohort':f'TCGA-{cohort}','n_primary_RNA_cases':len(rna_map),'n_model_complete_cases':len(common),
                   'n_events':int(clin.event.sum()),'n_censored':int((clin.event==0).sum()),
                   'RNA_features':layer_mats['RNA'].shape[0],'CNV_features':layer_mats['CNV'].shape[0],
                   'METHYLATION_features':layer_mats['METHYLATION'].shape[0],
                   'SNV_exploratory_features':snv_n_features}]).to_csv(proc/'cohort_qc.csv',index=False)
    print(f'\n{cohort} prepared. Complete cases across all three layers and survival: {len(common)}')
    print(f'Eventos: {int(clin.event.sum())}; censurados: {int((clin.event==0).sum())}')
    print(f'Se conservan RNA={layer_mats["RNA"].shape}, CNV={layer_mats["CNV"].shape}, Methylation={layer_mats["METHYLATION"].shape}.')
    return 0

if __name__=='__main__': raise SystemExit(main())
