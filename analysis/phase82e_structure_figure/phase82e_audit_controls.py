#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, re
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

LABELS='abcdefg'
NEG3={'ASP','GLU'}; POS3={'LYS','ARG'}
CHARGE_ATOMS={'ASP':['OD1','OD2'],'GLU':['OE1','OE2'],'LYS':['NZ'],'ARG':['NH1','NH2','NE']}

def parse_pdb(path: Path):
    chains=defaultdict(lambda: defaultdict(dict)); rn={}
    for line in path.read_text(errors='ignore').splitlines():
        if not line.startswith(('ATOM  ','HETATM')): continue
        atom=line[12:16].strip(); alt=line[16:17]
        if alt not in (' ','A'): continue
        resn=line[17:20].strip(); ch=line[21:22].strip() or '_'
        try:
            ri=int(line[22:26]); xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])],float)
        except Exception:
            continue
        chains[ch][ri][atom]=xyz; rn[(ch,ri)]=resn
    return chains,rn

def axis(resdict):
    keys=sorted(k for k,v in resdict.items() if 'CA' in v)
    if len(keys)<4: return None
    q=max(2,len(keys)//4)
    p0=np.mean([resdict[k]['CA'] for k in keys[:q]],0); p1=np.mean([resdict[k]['CA'] for k in keys[-q:]],0)
    v=p1-p0; n=np.linalg.norm(v)
    return None if n==0 else v/n

def heavy_atoms(atoms): return [(n,c) for n,c in atoms.items() if not n.startswith('H')]
def min_dist(a,b):
    aa=heavy_atoms(a); bb=heavy_atoms(b)
    if not aa or not bb: return float('nan')
    return min(float(np.linalg.norm(ca-cb)) for _,ca in aa for _,cb in bb)

def charge_center(resn, atoms):
    pts=[atoms[n] for n in CHARGE_ATOMS.get(resn,[]) if n in atoms]
    return np.mean(pts,axis=0) if pts else None

def rank_from_name(path: Path):
    m=re.search(r'rank[_-](\d+)',path.name); return int(m.group(1)) if m else 999

def model_seed_key(path: Path):
    m=re.search(r'model_(\d+)_seed_(\d+)',path.name)
    return f'm{int(m.group(1))}_s{int(m.group(2))}' if m else path.name

def exact_score_json(pdb: Path):
    for tag in ('_unrelaxed_','_relaxed_'):
        if tag in pdb.name:
            j=pdb.parent/pdb.name.replace(tag,'_scores_').replace('.pdb','.json')
            if j.exists(): return j
    return None

def heptad_map(length:int,a_start:int):
    return {i:LABELS[(i-a_start)%7] for i in range(1,length+1)}

def analyze_pdb(pdb: Path, chainA_protein: str, chainB_protein: str, seqA: str, seqB: str, a_start: int, label: str):
    chains,rn=parse_pdb(pdb); chs=sorted(chains)
    if len(chs)<2: raise RuntimeError(f'<2 chains: {pdb}')
    A,B=chs[:2]; va,vb=axis(chains[A]),axis(chains[B])
    dot=float(np.dot(va,vb)); ori='parallel' if dot>=0.5 else ('antiparallel' if dot<=-0.5 else 'oblique')
    hmA=heptad_map(len(seqA),a_start); hmB=heptad_map(len(seqB),a_start)
    contacts=core=opp=same=0
    for ia,aa in chains[A].items():
        for ib,bb in chains[B].items():
            d=min_dist(aa,bb)
            if math.isnan(d) or d>6: continue
            if d<5:
                contacts += 1
                if hmA.get(ia,'') in 'ad' and hmB.get(ib,'') in 'ad': core += 1
            ra,rbb=rn.get((A,ia),'UNK'),rn.get((B,ib),'UNK')
            ca,cb=charge_center(ra,aa),charge_center(rbb,bb)
            if ca is None or cb is None: continue
            if hmA.get(ia,'') not in 'eg' or hmB.get(ib,'') not in 'eg': continue
            cd=float(np.linalg.norm(ca-cb))
            if cd>=6 or ra not in NEG3|POS3 or rbb not in NEG3|POS3: continue
            isopp=(ra in NEG3 and rbb in POS3) or (ra in POS3 and rbb in NEG3)
            if isopp: opp += 1
            else: same += 1
    score={}; sj=exact_score_json(pdb)
    if sj:
        try: score=json.loads(sj.read_text())
        except Exception: score={}
    return dict(label=label,pdb=str(pdb),rank=rank_from_name(pdb),model_seed_key=model_seed_key(pdb),
                chainA_protein=chainA_protein,chainB_protein=chainB_protein,orientation=ori,axis_dot=dot,
                contacts_lt5=contacts,ad_ad_lt5=core,eg_opp_lt6=opp,eg_same_lt6=same,
                iptm=float(score.get('iptm',np.nan)),ptm=float(score.get('ptm',np.nan)))

def resolve_path(p: str, root: Path, fallback_dir: Path):
    x=Path(str(p))
    if x.exists(): return x
    y=fallback_dir/x.name
    if y.exists(): return y
    hits=list(root.rglob(x.name))
    if len(hits)==1: return hits[0]
    raise FileNotFoundError(f'Cannot resolve {p}')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('project_root',type=Path)
    args=ap.parse_args(); root=args.project_root.resolve()
    cdir=root/'results/phase82c_modeled_structure'; ddir=root/'results/phase82d_matched_orientation'
    out=root/'results/phase82e_structure_figure'; out.mkdir(parents=True,exist_ok=True)
    contract=json.loads((cdir/'MODEL_INPUT_CONTRACT.json').read_text()); a_start=int(contract['heptad_a_start_1based'])
    selected=pd.read_csv(cdir/'selected_rank1_models.csv')
    geom_c=pd.read_csv(cdir/'model_geometry_by_prediction.csv')
    inv=pd.read_csv(ddir/'matched_orientation_candidate_inventory.csv')
    geom_d=pd.read_csv(ddir/'matched_candidate_geometry.csv')
    # Sequences.
    seqs={contract['shared_member']: str(inv.iloc[0].shared_sequence)}
    for r in inv.itertuples(index=False): seqs[str(r.partner)]=str(r.partner_sequence)
    # target partner / hard partner from contract
    target_partner=contract['target_model_order'][1]; hard_partner=contract['competitor_model_order'][1]
    # add sequences from Phase82C map if needed
    hmap=pd.read_csv(cdir/'sequence_heptad_map.csv')
    for prot,g in hmap.groupby('protein'):
        seqs[str(prot)]=str(g.sort_values('position').iloc[0].sequence)
    rows=[]
    # Phase82C all 50: resolve PDBs and recompute to keep one contact implementation.
    for r in geom_c.itertuples(index=False):
        p=resolve_path(r.pdb,root,cdir/'colabfold_main')
        label='target' if r.query=='main_target' else 'hard_competitor'
        partner=target_partner if label=='target' else hard_partner
        rows.append(analyze_pdb(p,contract['shared_member'],partner,seqs[contract['shared_member']],seqs[partner],a_start,label))
    # Phase82D candidates all models.
    for r in geom_d.itertuples(index=False):
        m=re.match(r'matchedcand_(\d+)__(.+)__(.+)$',str(r.query))
        if not m: continue
        fr=int(m.group(1)); shared=m.group(2); partner=m.group(3)
        p=resolve_path(r.pdb,root,ddir/'colabfold_candidates')
        label={'4H140-AA':'strict_parallel_control','4H1802-AA':'hard_parallel_favored_control'}.get(partner,f'candidate_{partner}')
        rows.append(analyze_pdb(p,shared,partner,seqs[shared],seqs[partner],a_start,label))
    m=pd.DataFrame(rows).sort_values(['label','rank','pdb']); m.to_csv(out/'phase82e_model_level_geometry.csv',index=False)
    # frozen pair scores
    score_map={'target':float(contract['target_score']),'hard_competitor':float(contract['competitor_score'])}
    for r in inv.itertuples(index=False):
        if r.partner=='4H140-AA': score_map['strict_parallel_control']=float(r.score)
        if r.partner=='4H1802-AA': score_map['hard_parallel_favored_control']=float(r.score)
    order=['target','hard_competitor','strict_parallel_control','hard_parallel_favored_control']
    summaries=[]
    for label in order:
        g=m[m.label.eq(label)].copy()
        if g.empty: continue
        cnt=g.orientation.value_counts(); top=g.sort_values('rank').head(5)
        summaries.append(dict(label=label,frozen_eol_xgb_score=score_map.get(label,np.nan),n_models=len(g),
            parallel_n=int(cnt.get('parallel',0)),antiparallel_n=int(cnt.get('antiparallel',0)),oblique_n=int(cnt.get('oblique',0)),
            top5_orientation_mode=Counter(top.orientation).most_common(1)[0][0],top5_orientation_consistency=float((top.orientation==Counter(top.orientation).most_common(1)[0][0]).mean()),
            mean_iptm=g.iptm.mean(),median_iptm=g.iptm.median(),mean_ptm=g.ptm.mean(),median_ptm=g.ptm.median(),
            mean_contacts=g.contacts_lt5.mean(),median_contacts=g.contacts_lt5.median(),mean_ad=g.ad_ad_lt5.mean(),median_ad=g.ad_ad_lt5.median(),
            mean_eg_opp=g.eg_opp_lt6.mean(),median_eg_opp=g.eg_opp_lt6.median(),mean_eg_same=g.eg_same_lt6.mean(),median_eg_same=g.eg_same_lt6.median(),
            rank1_pdb=str(g.sort_values('rank').iloc[0].pdb)))
    s=pd.DataFrame(summaries); s.to_csv(out/'phase82e_case_summary.csv',index=False)
    # 4H1802 parallel-only descriptive summary.
    hp=m[(m.label=='hard_parallel_favored_control')&(m.orientation=='parallel')]
    if len(hp):
        pd.DataFrame([dict(n_models=len(hp),mean_iptm=hp.iptm.mean(),median_iptm=hp.iptm.median(),mean_ptm=hp.ptm.mean(),median_ptm=hp.ptm.median(),
                           mean_contacts=hp.contacts_lt5.mean(),median_contacts=hp.contacts_lt5.median(),mean_ad=hp.ad_ad_lt5.mean(),median_ad=hp.ad_ad_lt5.median(),
                           mean_eg_opp=hp.eg_opp_lt6.mean(),median_eg_opp=hp.eg_opp_lt6.median(),mean_eg_same=hp.eg_same_lt6.mean(),median_eg_same=hp.eg_same_lt6.median())]).to_csv(out/'phase82e_4H1802_parallel_subset_summary.csv',index=False)
    # Paired model+seed contrasts versus target; descriptive only, not inferential statistics.
    base=m[m.label.eq('target')].set_index('model_seed_key')
    prs=[]
    for label in ['hard_competitor','strict_parallel_control','hard_parallel_favored_control']:
        other=m[m.label.eq(label)].set_index('model_seed_key'); j=base.join(other,lsuffix='_target',rsuffix='_control',how='inner')
        for metric in ['contacts_lt5','ad_ad_lt5','eg_opp_lt6','eg_same_lt6','iptm','ptm']:
            d=j[f'{metric}_target']-j[f'{metric}_control']
            prs.append(dict(control=label,metric=metric,n=len(d),target_mean=j[f'{metric}_target'].mean(),control_mean=j[f'{metric}_control'].mean(),delta_target_minus_control=d.mean(),target_gt_control=int((d>0).sum()),equal=int((d==0).sum()),target_lt_control=int((d<0).sum())))
    pd.DataFrame(prs).to_csv(out/'phase82e_paired_model_seed_contrasts.csv',index=False)
    gate={
      'phase':'Phase82E structural-case audit',
      'status':'READY_FOR_MATCHED_VIEW_RENDERING',
      'primary_case':'target versus strongest frozen E+O+L XGB shared-member off-target; do not replace after structure inspection',
      'strict_matched_orientation_control':'4H3651-AA + 4H140-AA (25/25 parallel); secondary sensitivity control only',
      'hard_parallel_favored_control':'4H3651-AA + 4H1802-AA (18/25 parallel; top5 5/5 parallel); supplementary sensitivity control',
      'experimental_structure':False,
      'inference_boundary':'AlphaFold/ColabFold coordinates are retrospective modeled illustrations. Model/seed replicates are descriptive ensemble samples, not independent experimental replicates.'
    }
    (out/'PHASE82E_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2)); print('\n[CASE SUMMARY]'); print(s.to_string(index=False))
    print(f'[OK] wrote {out}')

if __name__=='__main__': main()
