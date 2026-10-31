#!/usr/bin/env python3
"""Cross-platform preflight: audit deliverables without changing sample cohorts."""
from __future__ import annotations
import argparse, hashlib, json, gzip
from pathlib import Path
import pandas as pd


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def tsv_columns(path: Path) -> list[str]:
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8',errors='replace') as f:
        return f.readline().rstrip('\r\n').split('\t')[1:]


def main() -> int:
    ap=argparse.ArgumentParser();ap.add_argument('--root',default='.');args=ap.parse_args();root=Path(args.root).resolve()
    records=[]; errors=[];summary={}
    (root/'results'/'QC').mkdir(parents=True,exist_ok=True)
    # LUAD and LUSC: verify same Xena-derived format and complete-case intersections.
    for cohort in ('LUAD','LUSC'):
        proc=root/'data/processed'/cohort
        files={x:proc/f'{x}.tsv.gz' for x in ['RNA','CNV','METHYLATION']}
        colsets={}
        for layer,p in files.items():
            ok=p.exists() and p.stat().st_size>0
            records.append({'cohort':cohort,'component':layer,'status':'READY' if ok else 'MISSING',
                            'path':str(p.relative_to(root)),'bytes':p.stat().st_size if ok else 0,
                            'sha256':sha(p) if ok else ''})
            if not ok: errors.append(f'{cohort} missing {layer}: {p}')
            else: colsets[layer]=set(tsv_columns(p))
        clin_path=proc/'clinical.csv'
        if not clin_path.exists():
            errors.append(f'{cohort} clinical.csv no existe.')
        else:
            d=pd.read_csv(clin_path)
            if not {'sample_id','time_days','event'}.issubset(d.columns):
                errors.append(f'{cohort} clinical requiere sample_id,time_days,event.')
            else:
                cl=set(d.sample_id.astype(str))
                inter=set.intersection(cl,*colsets.values()) if len(colsets)==3 else set()
                summary[cohort]={'clinical_rows':len(d),'unique_sample_ids':int(d.sample_id.nunique()),
                  'events':int(pd.to_numeric(d.event,errors='coerce').eq(1).sum()),
                  'censored':int(pd.to_numeric(d.event,errors='coerce').eq(0).sum()),
                  'common_samples_all_omics_and_clinical':len(inter),
                  'selection_rule':'primary tumour 01, one sample per case, complete RNA/CNV/HM450 and OS'}
                if d.sample_id.duplicated().any(): errors.append(f'{cohort} clinical contiene sample_id duplicados.')
                tm=pd.to_numeric(d.time_days,errors='coerce'); ev=pd.to_numeric(d.event,errors='coerce')
                if tm.isna().any() or (tm<=0).any(): errors.append(f'{cohort} has missing/non-positive OS times.')
                if ev.isna().any() or not ev.isin([0,1]).all(): errors.append(f'{cohort} has missing OS events or event values other than 0/1.')
                if len(inter)<100: errors.append(f'{cohort} RNA/CNV/methylation+survival intersection <100 ({len(inter)}).')
                if len(inter)!=len(cl): errors.append(f'{cohort} matriz/clínica no tienen exactamente los mismos sample_id.')
        manifest=root/'data/raw'/cohort/f'{cohort.lower()}_download_manifest.csv'
        if not manifest.exists(): errors.append(f'Download manifest missing {cohort}: {manifest}')
        else:
            dm=pd.read_csv(manifest)
            required_layers={'RNA','CNV','METHYLATION','CLINICAL','SURVIVAL'}
            ready=set(dm.loc[dm.status.isin(['downloaded','reused']),'layer'])
            missing=required_layers-ready
            if missing: errors.append(f'{cohort}: manifiesto sin capas obligatorias {sorted(missing)}')
    # RA fresh-start acquisition and processed four-layer intersection.
    ra_raw=root/'data/raw/RA'
    ra_manifest=ra_raw/'geo_download_manifest.csv'
    rna_path=ra_raw/'RNA'/'GSE112656_Counts_matix.txt.gz'
    ra_processed=root/'data/processed/RA'
    ra_rds=ra_processed/'ra_modalities.rds'
    ra_donors=ra_processed/'ra_manifest.csv'
    ra_layer_qc=ra_processed/'ra_layer_qc.csv'
    for label,path in [('RNA_counts',rna_path),('GEO_download_manifest',ra_manifest),
                       ('RA_paired_RDS',ra_rds),('RA_donor_manifest',ra_donors),
                       ('RA_layer_QC',ra_layer_qc)]:
        ok=path.exists() and path.stat().st_size>0
        records.append({'cohort':'RA','component':label,'status':'READY' if ok else 'MISSING',
                        'path':str(path.relative_to(root)),'bytes':path.stat().st_size if ok else 0,
                        'sha256':sha(path) if ok and path.is_file() else ''})
        if not ok: errors.append(f'RA missing {label}: {path}')
    if ra_manifest.exists():
        rd=pd.read_csv(ra_manifest)
        required_assays={'RNA','ATAC','H3K27ac','H3K4me3'}
        ready_assays=set(rd.loc[rd.status.isin(['downloaded','reused_project_download']),'assay']) if 'assay' in rd and 'status' in rd else set()
        missing_assays=required_assays-ready_assays
        if missing_assays: errors.append(f'RA GEO acquisition missing assays: {sorted(missing_assays)}')
        if 'path' in rd:
            absent=[str(x) for x in rd.loc[rd.assay.isin(required_assays)&rd.path.notna(),'path'] if str(x) and not Path(str(x)).exists()]
            if absent: errors.append(f'RA GEO manifest points to missing files: {absent[:5]}')
    if ra_donors.exists():
        rm=pd.read_csv(ra_donors)
        if not {'donor','group'}.issubset(rm.columns): errors.append('RA donor manifest requires donor,group.')
        else:
            if rm.donor.duplicated().any(): errors.append('RA donor manifest has duplicate donors.')
            if len(rm)<12 or set(rm.group.astype(str))!={'RA','OA'}: errors.append(f'RA complete matched donor set is insufficient or lacks both groups (n={len(rm)}).')
            summary['RA']={'complete_donors':len(rm),'RA':int((rm.group=='RA').sum()),'OA':int((rm.group=='OA').sum()),
                           'modalities':['RNA','ATAC','H3K27ac','H3K4me3']}
    # No WGBS, FASTQ/SRA or full GEO RAW TAR should be included in this compact workflow.
    forbidden=(list(ra_raw.glob('**/*WGBS*'))+list(ra_raw.glob('**/*.fastq*'))+
               list(ra_raw.glob('**/*.fq*'))+list(ra_raw.glob('**/*.sra'))+list(ra_raw.glob('**/*.tar'))+
               list(ra_raw.glob('**/*.tar.gz')))
    if forbidden: errors.append('RA raw folder contains excluded heavy files: '+', '.join(str(x) for x in forbidden[:5]))
    # Provenance digest list for raw download manifests.
    provenance=[]
    for p in [root/'data/raw/LUAD/luad_download_manifest.csv',root/'data/raw/LUSC/lusc_download_manifest.csv',root/'data/raw/RA/geo_download_manifest.csv']:
        if p.exists():
            df=pd.read_csv(p)
            for _,r in df.iterrows(): provenance.append({'manifest':str(p.relative_to(root)),**{k:(str(v) if pd.notna(v) else '') for k,v in r.items()}})
    pd.DataFrame(records).to_csv(root/'results/QC/component_status.csv',index=False)
    if provenance: pd.DataFrame(provenance).to_csv(root/'results/QC/provenance_aggregated.csv',index=False)
    summary['preflight']={'status':'PASS' if not errors else 'FAIL','errors':errors}
    (root/'results/QC/preflight_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(summary,indent=2,ensure_ascii=False))
    if errors:return 2
    print('Preflight passed. No samples were modified or selected.')
    return 0

if __name__=='__main__':raise SystemExit(main())
