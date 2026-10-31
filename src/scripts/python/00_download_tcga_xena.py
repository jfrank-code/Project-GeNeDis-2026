#!/usr/bin/env python3
"""Fresh-start acquisition of TCGA-LUAD and TCGA-LUSC processed multi-omics from UCSC Xena.
Downloads open/tabular processed datasets only; no GDC controlled-access files are needed.
"""
from __future__ import annotations
import argparse, hashlib, json, re, sys, time
from pathlib import Path
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

HEADERS = {"User-Agent": "MOSAIC3-GeNeDis/2.1 (open reproducible research workflow; contact: project README)"}
TIMEOUT = (30, 300)


def make_session() -> requests.Session:
    session = requests.Session()
    session.headers.update(HEADERS)
    retry = Retry(
        total=5, connect=5, read=3, backoff_factor=1.2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}), raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=4, pool_maxsize=4)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def download_one(session: requests.Session, urls: list[str], base: Path, force: bool=False, optional: bool=False) -> dict:
    base.parent.mkdir(parents=True,exist_ok=True)
    candidates=[base.with_suffix('.tsv.gz'),base.with_suffix('.tsv')]
    if not force:
        for p in candidates:
            if p.exists() and p.stat().st_size>1024:
                return {'path':str(p),'url':urls[0],'bytes':p.stat().st_size,'sha256':sha256(p),'status':'reused'}
    part=base.parent/(base.name+'.part')
    if force: part.unlink(missing_ok=True)
    errors=[]
    for url in urls:
        try:
            offset=part.stat().st_size if part.exists() else 0
            hdr={'Range':f'bytes={offset}-'} if offset else {}
            with session.get(url,headers=hdr,stream=True,timeout=TIMEOUT) as r:
                if r.status_code==404:
                    errors.append(f'404 {url}'); continue
                r.raise_for_status()
                ctype=r.headers.get('Content-Type','').lower()
                if 'text/html' in ctype:
                    errors.append(f'HTML returned by {url}'); continue
                append=bool(offset and r.status_code==206)
                if offset and not append: offset=0
                total=int(r.headers.get('Content-Length','0') or 0)+offset
                written=offset; last=time.time()
                with part.open('ab' if append else 'wb') as out:
                    for chunk in r.iter_content(1024*1024):
                        if chunk:
                            out.write(chunk); written+=len(chunk)
                            if time.time()-last>3:
                                if total: print(f'    {written/2**20:,.0f}/{total/2**20:,.0f} MiB ({100*written/total:.1f}%)',flush=True)
                                else: print(f'    {written/2**20:,.0f} MiB downloaded',flush=True)
                                last=time.time()
            if not part.exists() or part.stat().st_size<1024:
                part.unlink(missing_ok=True); errors.append(f'empty download: {url}'); continue
            with part.open('rb') as f: magic=f.read(2)
            if magic==b'\x1f\x8b':
                import gzip
                with gzip.open(part,'rb') as check: check.read(128)
                suffix='.tsv.gz'
            else:
                with part.open('rb') as check: head=check.read(4096)
                if b'\t' not in head:
                    raise ValueError('Downloaded object does not appear to be a valid TSV matrix.')
                suffix='.tsv'
            target=base.with_suffix(suffix)
            target.unlink(missing_ok=True); part.replace(target)
            return {'path':str(target),'url':url,'bytes':target.stat().st_size,'sha256':sha256(target),'status':'downloaded'}
        except Exception as e:
            errors.append(f'{url}: {e}')
            print(f'    Warning: {e}',file=sys.stderr)
            # Keep a partial transfer for the next execution, which retries the same first URL.
            if part.exists() and part.stat().st_size>0:
                break
    part.unlink(missing_ok=True)
    if optional:
        return {'path':'','url':urls[0],'bytes':0,'sha256':'','status':'not_available','note':' | '.join(errors)}
    raise RuntimeError(f'Download failed {urls[0]}: ' + ' | '.join(errors))


def dataset_urls(host: str, dataset: str, is_survival: bool=False, survival_host: str='') -> list[str]:
    hosts=[survival_host,host] if is_survival else [host]
    out=[]
    for h in hosts:
        h=h.rstrip('/')
        out.extend([h+'/download/'+dataset, h+'/download/'+dataset+'.gz'])
    return list(dict.fromkeys(out))


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='.')
    ap.add_argument('--force',action='store_true')
    args=ap.parse_args()
    root=Path(args.root).resolve(); config=json.loads((root/'config/mosaic3_config.json').read_text(encoding='utf-8'))
    sess=make_session()
    for cohort in ('LUAD','LUSC'):
        cfg=config[cohort.lower()]; raw=root/'data'/'raw'/cohort
        print(f'\n===== TCGA-{cohort}: adquisición desde UCSC Xena =====')
        rows=[]; host=cfg['xena_host']; survival_host=cfg.get('survival_host',host)
        for layer,dataset in cfg['datasets'].items():
            if layer=='SNV' and not cfg.get('download_snv',False):
                rows.append({'cohort':cohort,'layer':layer,'dataset':dataset,'path':'','bytes':0,'sha256':'','status':'skipped_optional'}); continue
            folder=raw/('CLINICAL' if layer in ('CLINICAL','SURVIVAL') else layer)
            base=folder/(f'{cohort}_survival' if layer=='SURVIVAL' else layer)
            is_surv=(layer=='SURVIVAL')
            urls=dataset_urls(host,dataset,is_surv,survival_host)
            print(f'[{cohort} / {layer}] {dataset}')
            try:
                rec=download_one(sess,urls,base,force=args.force,optional=(layer=='SNV'))
            except Exception as exc:
                rec={'path':'','url':urls[0],'bytes':0,'sha256':'','status':'failed','note':str(exc)}
                rec.update({'cohort':cohort,'layer':layer,'dataset':dataset})
                rows.append(rec)
                pd.DataFrame(rows).to_csv(raw/f'{cohort.lower()}_download_manifest.csv',index=False)
                raise
            rec.update({'cohort':cohort,'layer':layer,'dataset':dataset})
            rows.append(rec)
            # Persist the manifest after each asset so an interrupted fresh run is auditable.
            pd.DataFrame(rows).to_csv(raw/f'{cohort.lower()}_download_manifest.csv',index=False)
            print(f"    {rec['status']} | {rec['bytes']/2**20:,.1f} MiB")
        manifest=pd.DataFrame(rows)
        manifest.to_csv(raw/f'{cohort.lower()}_download_manifest.csv',index=False)
        required={'RNA','CNV','METHYLATION','CLINICAL','SURVIVAL'}
        missing=required-set(manifest.loc[manifest.status.isin(['downloaded','reused']),'layer'])
        if missing: raise SystemExit(f'TCGA-{cohort}: required datasets are missing: {sorted(missing)}')
        print(f'TCGA-{cohort}: required downloads complete.')
    print('\nTCGA acquisition complete. Files are stored in data/raw and reused on subsequent runs.')
    return 0

if __name__=='__main__': raise SystemExit(main())
