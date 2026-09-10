#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import pandas as pd

def parse_pdb(path:Path):
    chains=defaultdict(lambda:defaultdict(dict))
    for line in path.read_text(errors='ignore').splitlines():
        if not line.startswith(('ATOM  ','HETATM')): continue
        atom=line[12:16].strip(); alt=line[16:17]
        if alt not in (' ','A'): continue
        ch=line[21:22].strip() or '_'
        try: ri=int(line[22:26]); xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        except: continue
        chains[ch][ri][atom]=xyz
    return chains

def axis(resdict):
    keys=sorted(k for k,v in resdict.items() if 'CA' in v)
    if len(keys)<4:return None
    q=max(2,len(keys)//4); p0=np.mean([resdict[k]['CA'] for k in keys[:q]],0); p1=np.mean([resdict[k]['CA'] for k in keys[-q:]],0)
    v=p1-p0; n=np.linalg.norm(v); return None if n==0 else v/n

def orient_of(p):
    c=parse_pdb(p); ch=sorted(c)
    if len(ch)<2:return 'invalid',float('nan')
    a,b=axis(c[ch[0]]),axis(c[ch[1]])
    if a is None or b is None:return 'invalid',float('nan')
    dot=float(np.dot(a,b)); o='parallel' if dot>=0.5 else ('antiparallel' if dot<=-0.5 else 'oblique')
    return o,dot

def rank(p):
    m=re.search(r'rank[_-](\d+)',p.name); return int(m.group(1)) if m else 999

def query(p):
    m=re.match(r'(matchedcand_\d+__[^_]+__[^_]+)',p.name)
    if m:return m.group(1)
    # IDs can contain hyphens but no double underscores; stop before ColabFold suffix.
    x=p.name.split('_unrelaxed_rank_')[0].split('_relaxed_rank_')[0]
    return x if x.startswith('matchedcand_') else None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('project_root',type=Path); a=ap.parse_args(); root=a.project_root.resolve()
    out=root/'results/phase82d_matched_orientation'; pred=out/'colabfold_candidates'
    inv=pd.read_csv(out/'matched_orientation_candidate_inventory.csv')
    rows=[]
    pdbs=sorted(pred.glob('*unrelaxed*.pdb'))
    for p in pdbs:
        q=query(p)
        if not q:continue
        o,dot=orient_of(p); rows.append({'query':q,'rank':rank(p),'orientation':o,'axis_dot':dot,'pdb':str(p)})
    if not rows: raise SystemExit(f'No modeled candidate PDBs found under {pred}')
    geom=pd.DataFrame(rows).sort_values(['query','rank']); geom.to_csv(out/'matched_candidate_geometry.csv',index=False)
    sums=[]
    for q,g in geom.groupby('query'):
        # parse rank/partner from query header written by Phase82D0
        m=re.match(r'matchedcand_(\d+)__(.+)__(.+)$',q)
        fr=int(m.group(1)); shared=m.group(2); partner=m.group(3)
        cnt=g.orientation.value_counts(); n=len(g); par=int(cnt.get('parallel',0)); anti=int(cnt.get('antiparallel',0)); obl=int(cnt.get('oblique',0))
        top=g.sort_values('rank').head(5); mode=Counter(top.orientation).most_common(1)[0][0]
        invrow=inv[inv.frozen_score_rank.eq(fr)].iloc[0]
        sums.append({'frozen_score_rank':fr,'shared_member':shared,'partner':partner,'frozen_score':float(invrow.score),
                     'n_models':n,'parallel_n':par,'antiparallel_n':anti,'oblique_n':obl,'parallel_fraction':par/n,
                     'top5_orientation_mode':mode,'top5_orientation_consistency':float((top.orientation==mode).mean())})
    s=pd.DataFrame(sums).sort_values('frozen_score_rank'); s.to_csv(out/'matched_candidate_orientation_summary.csv',index=False)
    strict=s[(s.n_models>=25)&(s.parallel_fraction.eq(1.0))].copy()
    selected=strict.sort_values('frozen_score_rank').head(1)
    if len(selected):
        selected.to_csv(out/'selected_matched_parallel_control.csv',index=False); status='STRICT_MATCHED_PARALLEL_CONTROL_FOUND'
    else:
        pd.DataFrame(columns=s.columns).to_csv(out/'selected_matched_parallel_control.csv',index=False); status='NO_STRICT_MATCHED_PARALLEL_CONTROL_FOUND'
    gate={'phase':'Phase82D matched-orientation control','status':status,'candidate_count_modeled':int(len(s)),
          'selection_rule':'highest frozen-score heterodimer off-target sharing the fixed query chain with 25/25 parallel AF2-Multimer orientations',
          'primary_strongest_competitor_replaced':False,
          'interpretation':'Secondary sensitivity control only. The original strongest off-target remains the primary hard competitor.'}
    (out/'MATCHED_ORIENTATION_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2)); print('\n[ORIENTATION SUMMARY]'); print(s.to_string(index=False))
    if len(selected): print('\n[SELECTED MATCHED PARALLEL CONTROL]\n'+selected.to_string(index=False))

if __name__=='__main__':main()
