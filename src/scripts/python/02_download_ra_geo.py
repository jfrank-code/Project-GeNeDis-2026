#!/usr/bin/env python3
"""Download only processed, donor-level RA GEO matrices and peaks. WGBS/FASTQ/TAR excluded."""
from __future__ import annotations
import argparse, gzip, hashlib, json, re, sys, time
from pathlib import Path
from urllib.parse import urlparse
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

HEADERS = {'User-Agent': 'MOSAIC3-GeNeDis/2.1 (open reproducible research workflow)'}
RNA_DONORS = [f'RA_{i:02d}' for i in [1,2,3,4,5,6,7,8,10,11]] + [f'OA_{i:02d}' for i in range(2,12)]
DEFAULT_MAX_PEAK_MB = 35


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for c in iter(lambda:f.read(1024*1024),b''): h.update(c)
    return h.hexdigest()


def fetch(url: str, dest: Path, session: requests.Session, max_mb: float | None = None, force=False) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000 and not force:
        return {'url':url,'path':str(dest),'bytes':dest.stat().st_size,'sha256':sha256(dest),'status':'reused_project_download'}
    part=dest.with_suffix(dest.suffix+'.part')
    if force:
        for p in [dest,part]:
            if p.exists(): p.unlink()
    offset=part.stat().st_size if part.exists() else 0
    headers={'Range':f'bytes={offset}-'} if offset else {}
    with session.get(url,stream=True,headers=headers,timeout=(30,180)) as r:
        if r.status_code==404: raise FileNotFoundError(url)
        r.raise_for_status()
        if 'text/html' in r.headers.get('Content-Type','').lower(): raise ValueError('Unexpected HTML response: '+url)
        mode='ab' if offset and r.status_code==206 else 'wb'
        if offset and r.status_code!=206: offset=0
        with part.open(mode) as f:
            for chunk in r.iter_content(1024*1024):
                if chunk: f.write(chunk)
    if not part.exists() or part.stat().st_size<1000:
        raise IOError('Incomplete download: '+url)
    if max_mb and part.stat().st_size > max_mb*1024*1024:
        part.unlink(missing_ok=True)
        raise IOError(f'File exceeds the size limit of {max_mb} MB: {url}')
    # Validate gzip integrity where applicable.
    with part.open('rb') as f: magic=f.read(2)
    if magic==b'\x1f\x8b':
        with gzip.open(part,'rb') as gz: gz.read(64)
    part.replace(dest)
    return {'url':url,'path':str(dest),'bytes':dest.stat().st_size,'sha256':sha256(dest),'status':'downloaded'}


def read_soft(url: str, session: requests.Session, dest: Path, force=False) -> str:
    rec=fetch(url,dest,session,max_mb=15,force=force)
    if dest.suffix=='.gz':
        with gzip.open(dest,'rt',errors='replace') as f: return f.read()
    return dest.read_text(errors='replace')


def parse_family_soft(text: str) -> list[dict]:
    records=[]; current=None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):
            if current: records.append(current)
            current={'gsm':line.split('=',1)[1].strip(),'title':'','supplementary':[]}
        elif current and line.startswith('!Sample_title = '):
            current['title']=line.split('=',1)[1].strip().strip('"')
        elif current and line.startswith('!Sample_supplementary_file = '):
            val=line.split('=',1)[1].strip().strip('"')
            if val and val.lower() not in {'none','null'}: current['supplementary'].append(val)
    if current: records.append(current)
    return records


def normalize_url(url: str) -> str:
    return url.replace('ftp://','https://',1) if url.startswith('ftp://') else url


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='.')
    ap.add_argument('--force',action='store_true')
    args=ap.parse_args(); root=Path(args.root).resolve()
    cfg=json.loads((root/'config/mosaic3_config.json').read_text(encoding='utf-8'))['ra']
    session=requests.Session(); session.headers.update(HEADERS)
    retry=Retry(total=5,connect=5,read=3,backoff_factor=1.2,status_forcelist=(429,500,502,503,504),allowed_methods=frozenset({'GET'}),raise_on_status=False)
    adapter=HTTPAdapter(max_retries=retry,pool_connections=4,pool_maxsize=4)
    session.mount('https://',adapter); session.mount('http://',adapter)
    raw=root/'data/raw/RA'; (raw/'RNA').mkdir(parents=True,exist_ok=True)
    (raw/'ATAC').mkdir(parents=True,exist_ok=True); (raw/'CHIP').mkdir(parents=True,exist_ok=True)
    manifest=[]
    # GEO processed count matrix; do not retrieve any SRA/FASTQ or the 4.8-GB SuperSeries TAR.
    rna_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112656/suppl/GSE112656_Counts_matix.txt.gz'
    rna_path=raw/'RNA/GSE112656_Counts_matix.txt.gz'
    print('[RA] RNA-seq processed count matrix')
    manifest.append(fetch(rna_url,rna_path,session,max_mb=20,force=args.force)|{'assay':'RNA','donor':'ALL','mark':''})
    # Small family SOFT metadata gives sample titles and per-GSM processed file locations.
    soft_maps={}
    for accession,assay in [('GSE112654','ATAC'),('GSE112655','CHIP')]:
        url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:6]}nnn/{accession}/soft/{accession}_family.soft.gz'
        soft_path=raw/f'{accession}_family.soft.gz'
        print(f'[RA] Metadata {accession}')
        text=read_soft(url,session,soft_path,force=args.force)
        soft_maps[assay]=parse_family_soft(text)
    marks=set(cfg.get('selected_chip_marks',['H3K27ac','H3K4me3']))
    max_peak_mb=float(cfg.get('max_individual_peak_file_mb',DEFAULT_MAX_PEAK_MB))
    found={'ATAC':set(),'H3K27ac':set(),'H3K4me3':set()}
    for assay, records in soft_maps.items():
        for rec in records:
            title=rec['title']
            # GEO titles are e.g. RA_01_ATAC / RA_01_H3K27ac.
            m=re.search(r'((?:RA|OA)_\d{2})',title,re.I)
            if not m: continue
            donor=m.group(1).upper()
            if donor not in RNA_DONORS: continue
            if assay=='ATAC':
                if 'ATAC' not in title.upper(): continue
                label='ATAC'; outdir=raw/'ATAC'
            else:
                mark=next((x for x in marks if x.lower() in title.lower()),None)
                if not mark: continue
                label=mark; outdir=raw/'CHIP'
            candidates=[]
            for u in rec['supplementary']:
                ul=u.lower()
                if label=='ATAC' and ('.bed' in ul or '.narrowpeak' in ul or '.broadpeak' in ul): candidates.append(normalize_url(u))
                elif label in marks and label.lower() in ul and ('.narrowpeak' in ul or '.broadpeak' in ul): candidates.append(normalize_url(u))
            # GEO's sample-level URL convention is used only if SOFT has no link.
            if not candidates:
                gsm=rec['gsm']; folder=gsm[:7]+'nnn'
                suffix='_peak.bed.gz' if label=='ATAC' else f'_{label}.narrowPeak.gz'
                candidates=[f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{folder}/{gsm}/suppl/{gsm}_{donor}{suffix}']
            dest_name=f'{donor}_{label}_'+Path(urlparse(candidates[0]).path).name
            dest=outdir/dest_name
            try:
                recd=fetch(candidates[0],dest,session,max_mb=max_peak_mb,force=args.force)
                recd.update({'assay':label,'donor':donor,'gsm':rec['gsm'],'mark':label})
                manifest.append(recd); found[label].add(donor)
                print(f'  {label:8s} {donor:6s} {recd["bytes"]/1024**2:5.2f} MB')
            except Exception as exc:
                print(f'  WARNING: {label} {donor}: {exc}',file=sys.stderr)
    pd.DataFrame(manifest).to_csv(raw/'geo_download_manifest.csv',index=False)
    # Donor-level count matrix is the anchor; each layer is allowed to be incomplete,
    # but preparation will report the exact complete intersection rather than fake it.
    for layer in ['ATAC','H3K27ac','H3K4me3']:
        n=len(found[layer])
        if n==0: raise SystemExit(f'No processed files were retrieved for {layer}. Check GEO access and internet connectivity.')
        print(f'{layer}: {n} RNA-compatible donors found')
    print('\nRA data selectively downloaded. WGBS, SRA/FASTQ and full TAR archives excluded.')
    print('Next step: Rscript src/scripts/R/01_prepare_ra_features.R')
    return 0

if __name__=='__main__': raise SystemExit(main())
