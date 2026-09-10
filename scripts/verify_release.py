#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, hashlib, json

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',default='.'); args=ap.parse_args(); root=Path(args.repo).resolve()
    mf=root/'MANIFEST.tsv'; rows=list(csv.DictReader(mf.open(encoding='utf-8'),delimiter='\t'))
    missing=[]; bad=[]
    for r in rows:
        p=root/r['path']
        if not p.is_file(): missing.append(r['path']); continue
        got=sha256(p)
        if got!=r['sha256']: bad.append({'path':r['path'],'expected':r['sha256'],'got':got})
    print(json.dumps({'manifest_entries':len(rows),'missing':missing,'sha_mismatches':bad,'pass':not missing and not bad},indent=2))
    if missing or bad: raise SystemExit(2)
if __name__=='__main__': main()
