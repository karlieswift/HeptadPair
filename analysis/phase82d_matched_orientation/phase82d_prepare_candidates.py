#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
from typing import Iterable, Optional
import numpy as np
import pandas as pd

AA_RE=re.compile(r'[^A-Z]')

def norm_seq(x):
    if x is None or (isinstance(x,float) and np.isnan(x)): return ''
    return AA_RE.sub('',str(x).upper().replace('U','X'))

def read_table(path:Path):
    suf=''.join(path.suffixes).lower()
    if suf.endswith('.parquet'): return pd.read_parquet(path)
    if suf.endswith('.csv') or suf.endswith('.csv.gz'): return pd.read_csv(path,low_memory=False)
    if suf.endswith('.tsv') or suf.endswith('.tsv.gz'): return pd.read_csv(path,sep='\t',low_memory=False)
    if suf.endswith('.xlsx') or suf.endswith('.xls'): return pd.read_excel(path)
    raise ValueError(path)

def col_lookup(cols:Iterable[str], names:Iterable[str])->Optional[str]:
    lower={str(c).lower():c for c in cols}
    for n in names:
        if n.lower() in lower: return lower[n.lower()]
    return None

def add_seqmap(df:pd.DataFrame, seqmap:dict[str,str])->int:
    cols=list(df.columns); before=len(seqmap)
    layouts=[
      (["seq1_id","protein1_id","protein1","id1","name1"],["seq2_id","protein2_id","protein2","id2","name2"],
       ["seq1","sequence1","sequence_1","protein1_sequence","sequence_a"],["seq2","sequence2","sequence_2","protein2_sequence","sequence_b"]),
      (["seq_low_id","protein_low_id"],["seq_high_id","protein_high_id"],["seq_low","sequence_low","seq_low_sequence"],["seq_high","sequence_high","seq_high_sequence"]),
    ]
    for ac,bc,s1c,s2c in layouts:
        a=col_lookup(cols,ac); b=col_lookup(cols,bc); s1=col_lookup(cols,s1c); s2=col_lookup(cols,s2c)
        if all(x is not None for x in (a,b,s1,s2)):
            for aa,bb,x1,x2 in df[[a,b,s1,s2]].dropna().itertuples(index=False,name=None):
                q1,q2=norm_seq(x1),norm_seq(x2)
                if q1: seqmap.setdefault(str(aa),q1)
                if q2: seqmap.setdefault(str(bb),q2)
    ic=col_lookup(cols,["protein","protein_id","seq_id","name","id","cc_name"])
    sc=col_lookup(cols,["sequence","seq","aa_sequence","protein_sequence"])
    if ic is not None and sc is not None:
        for ii,ss in df[[ic,sc]].dropna().itertuples(index=False,name=None):
            q=norm_seq(ss)
            if q: seqmap.setdefault(str(ii),q)
    return len(seqmap)-before

def build_seqmap(root:Path, needed:set[str], out:Path):
    preferred=[
      root/'data/external/ccmax/ccmax_pairs.parquet',
      root/'data/external/ccmax/ccmax_pairs.csv',
      root/'results/phase7/ccmax/ccmax_frozen_features.parquet',
      root/'results/phase74/ccmax_supervised/ccmax_pair_augmented_features.parquet',
    ]
    candidates=[p for p in preferred if p.exists()]
    for base in [root/'data/external/ccmax',root/'data/external/ngb2h_official',root/'results/phase7/ccmax',root/'results/phase74/ccmax_supervised']:
        if not base.exists(): continue
        for p in base.rglob('*'):
            if p.is_file() and ''.join(p.suffixes).lower() in {'.csv','.csv.gz','.tsv','.tsv.gz','.parquet','.xlsx','.xls'} and p not in candidates:
                candidates.append(p)
    seqmap={}; audit=[]
    for p in candidates:
        if needed.issubset(seqmap): break
        try:
            df=read_table(p); added=add_seqmap(df,seqmap)
            audit.append({'path':str(p),'read_ok':True,'rows':len(df),'sequences_added':added,'needed_mapped_after':len(needed & set(seqmap))})
        except Exception as e:
            audit.append({'path':str(p),'read_ok':False,'error':repr(e),'sequences_added':0,'needed_mapped_after':len(needed & set(seqmap))})
    pd.DataFrame(audit).to_csv(out/'sequence_source_audit.csv',index=False)
    return seqmap, sorted(needed-set(seqmap))

def as_bool(x):
    if isinstance(x,(bool,np.bool_)): return bool(x)
    return str(x).strip().lower() in {'1','true','t','yes','y'}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('project_root',type=Path)
    args=ap.parse_args(); root=args.project_root.resolve()
    cdir=root/'results/phase82c_modeled_structure'
    adir=root/'results/phase82_structure_case_preflight'
    out=root/'results/phase82d_matched_orientation'; out.mkdir(parents=True,exist_ok=True)
    contract=json.loads((cdir/'MODEL_INPUT_CONTRACT.json').read_text())
    pairs=pd.read_csv(adir/'normalized_pair_scores_by_set.csv')
    sid=str(contract['official_set_id']); variant=str(contract['variant']); shared=str(contract['shared_member'])
    target_partner=str(contract['target_model_order'][1]); current_partner=str(contract['competitor_model_order'][1])
    g=pairs[(pairs.official_set_id.astype(str).eq(sid)) & (pairs.variant.astype(str).eq(variant))].copy()
    if g.empty: raise SystemExit('No Phase82A pair-score rows for fixed case set/variant')
    g['is_target_bool']=g['is_target'].map(as_bool)
    cand=g[(~g.is_target_bool) & ((g.pair_a.astype(str).eq(shared)) | (g.pair_b.astype(str).eq(shared)))].copy()
    cand['partner']=[b if str(a)==shared else a for a,b in zip(cand.pair_a.astype(str),cand.pair_b.astype(str))]
    # Matched-orientation control should remain a heterodimer; exclude self.
    cand=cand[cand.partner.astype(str).ne(shared)].copy()
    cand=cand.sort_values(['score','partner'],ascending=[False,True]).reset_index(drop=True)
    cand['frozen_score_rank']=np.arange(1,len(cand)+1)
    cand['is_current_strongest_competitor']=cand.partner.astype(str).eq(current_partner)
    target_row=g[g.is_target_bool & (((g.pair_a.astype(str).eq(shared)) & (g.pair_b.astype(str).eq(target_partner))) | ((g.pair_b.astype(str).eq(shared)) & (g.pair_a.astype(str).eq(target_partner))))]
    target_score=float(target_row.iloc[0].score) if len(target_row) else float(contract['target_score'])
    cand['target_score']=target_score; cand['target_minus_candidate']=target_score-pd.to_numeric(cand.score)
    needed={shared} | set(cand.partner.astype(str)) | {target_partner}
    seqmap,missing=build_seqmap(root,needed,out)
    cand['shared_sequence']=seqmap.get(shared,'')
    cand['partner_sequence']=cand.partner.astype(str).map(seqmap).fillna('')
    cand['shared_sequence_length']=cand.shared_sequence.map(len)
    cand['partner_sequence_length']=cand.partner_sequence.map(len)
    cand.to_csv(out/'matched_orientation_candidate_inventory.csv',index=False)
    # Model all remaining heterodimer off-targets except the already-modeled strongest competitor.
    model=cand[~cand.is_current_strongest_competitor].copy()
    fasta=out/'matched_parallel_candidate_batch.fasta'
    with fasta.open('w') as fh:
        for r in model.itertuples(index=False):
            if not r.shared_sequence or not r.partner_sequence: continue
            fh.write(f'>matchedcand_{int(r.frozen_score_rank):02d}__{shared}__{r.partner}\n{r.shared_sequence}:{r.partner_sequence}\n')
    status='MATCHED_ORIENTATION_CANDIDATES_READY' if not missing and len(model)>0 else 'REVIEW_SEQUENCE_MAPPING_OR_EMPTY_CANDIDATES'
    manifest={
      'phase':'Phase82D matched-orientation off-target control preparation','status':status,
      'official_set_id':sid,'variant':variant,'shared_member':shared,'fixed_target_partner':target_partner,
      'current_strongest_competitor_partner':current_partner,'target_score':target_score,
      'candidate_count_including_current':int(len(cand)),'new_candidates_to_model':int(len(model)),
      'missing_sequences':missing,
      'selection_rule':'Keep the Phase82C target fixed. Among heterodimer off-targets sharing the same query chain, order candidates by the existing frozen E+O+L XGB score before structure modeling. After identical AF2-Multimer modeling, choose the highest-scoring candidate with 25/25 parallel orientation as a secondary matched-orientation control. Do not replace the strongest overall competitor in the primary case.',
      'retrained_model':False,'retuned_model':False,'structure_used_for_primary_competitor_selection':False,
    }
    (out/'MATCHED_ORIENTATION_PREP.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
    print('\n[CANDIDATES IN FROZEN SCORE ORDER]')
    print(cand[['frozen_score_rank','partner','score','target_minus_candidate','is_current_strongest_competitor']].to_string(index=False))
    if status.startswith('REVIEW_'): raise SystemExit(3)

if __name__=='__main__': main()
