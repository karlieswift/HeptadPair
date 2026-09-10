#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, re
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

HYDRO3={"ALA","VAL","ILE","LEU","MET","PHE","TRP","TYR"}
NEG3={"ASP","GLU"}; POS3={"LYS","ARG"}
CHARGE_ATOMS={
 "ASP":["OD1","OD2"],"GLU":["OE1","OE2"],"LYS":["NZ"],"ARG":["NH1","NH2","NE"]
}


def json_safe(obj):
    """Recursively convert NumPy/pandas/path values to strict JSON-safe Python types."""
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if pd.isna(obj) if not isinstance(obj, (str, bytes, dict, list, tuple)) else False:
        return None
    return obj


def parse_pdb(path:Path):
    chains=defaultdict(lambda: defaultdict(dict))
    resname={}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM  ","HETATM")): continue
        atom=line[12:16].strip(); alt=line[16:17]
        if alt not in (" ","A"): continue
        rn=line[17:20].strip(); ch=line[21:22].strip() or "_"
        try: ri=int(line[22:26]); x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
        except: continue
        chains[ch][ri][atom]=np.array([x,y,z],float); resname[(ch,ri)]=rn
    return chains,resname


def ca_vector(resdict):
    keys=sorted(k for k,v in resdict.items() if "CA" in v)
    if len(keys)<4: return None
    # Robust axis from first/last quartile centroids.
    q=max(2,len(keys)//4)
    p0=np.mean([resdict[k]["CA"] for k in keys[:q]],axis=0)
    p1=np.mean([resdict[k]["CA"] for k in keys[-q:]],axis=0)
    v=p1-p0
    n=np.linalg.norm(v)
    return None if n==0 else v/n


def helix_like_fraction(resdict):
    keys=sorted(k for k,v in resdict.items() if "CA" in v)
    vals=[]
    for idx in range(len(keys)-4):
        i=keys[idx]
        if keys[idx+3]==i+3 and keys[idx+4]==i+4:
            d3=np.linalg.norm(resdict[i]["CA"]-resdict[i+3]["CA"])
            d4=np.linalg.norm(resdict[i]["CA"]-resdict[i+4]["CA"])
            vals.append(1.0 if 4.3<=d3<=6.3 and 5.2<=d4<=7.2 else 0.0)
    return float(np.mean(vals)) if vals else float("nan")


def heavy_atoms(atoms):
    return [(n,c) for n,c in atoms.items() if not n.startswith("H")]


def min_res_distance(a,b):
    aa=heavy_atoms(a); bb=heavy_atoms(b)
    if not aa or not bb: return float("nan"),None,None
    best=(1e9,None,None)
    for na,ca in aa:
        for nb,cb in bb:
            d=float(np.linalg.norm(ca-cb))
            if d<best[0]: best=(d,na,nb)
    return best


def charge_center(rn,atoms):
    names=CHARGE_ATOMS.get(rn,[])
    pts=[atoms[n] for n in names if n in atoms]
    return np.mean(pts,axis=0) if pts else None


def load_heptad_map(path):
    df=pd.read_csv(path)
    return {(str(r.protein),int(r.position)):str(r.heptad) for r in df.itertuples(index=False)}


def rank_from_name(p:Path):
    m=re.search(r"rank[_-](\d+)",p.name)
    return int(m.group(1)) if m else 999


def query_from_name(p:Path):
    n=p.name
    for q in ["main_target","main_competitor"]:
        if n.startswith(q) or q in n: return q
    return None


def find_score_json(pdb:Path):
    # ColabFold normally writes a score JSON with the exact same query/model/seed
    # suffix as the PDB. Prefer that exact mapping first. This is essential when
    # multiple queries share rank numbers in the same output directory.
    exact_names = []
    for tag in ("_unrelaxed_", "_relaxed_"):
        if tag in pdb.name:
            exact_names.append(pdb.name.replace(tag, "_scores_").replace(".pdb", ".json"))
    for name in exact_names:
        j = pdb.parent / name
        if j.exists():
            return j

    # Conservative fallback for ColabFold naming variants: require BOTH the
    # same query identity and the same rank. Never match rank alone because
    # target and competitor each have rank_001, rank_002, ... files.
    q = query_from_name(pdb)
    candidates=[]
    for j in pdb.parent.glob("*.json"):
        if "score" not in j.name and "scores" not in j.name:
            continue
        if q is not None and query_from_name(j) != q:
            continue
        if rank_from_name(j)==rank_from_name(pdb):
            candidates.append(j)
    return sorted(candidates,key=lambda x:(len(x.name),x.name))[0] if candidates else None


def analyze_one(pdb:Path, qmeta:dict, hmap):
    chains,rn=parse_pdb(pdb)
    chs=sorted(chains)
    if len(chs)<2: return {"pdb":str(pdb),"query":qmeta["query"],"status":"<2_CHAINS"},[]
    A,B=chs[:2]
    va,vb=ca_vector(chains[A]),ca_vector(chains[B])
    orient=float(np.dot(va,vb)) if va is not None and vb is not None else float("nan")
    orientation="parallel" if orient>=0.5 else ("antiparallel" if orient<=-0.5 else "oblique")
    contacts=[]
    residue_pairs=0; core_contacts=0; eg_opposite_close=0; eg_same_close=0
    for ia,aa_atoms in chains[A].items():
        for ib,bb_atoms in chains[B].items():
            d,atA,atB=min_res_distance(aa_atoms,bb_atoms)
            if math.isnan(d) or d>6.0: continue
            residue_pairs += int(d<5.0)
            ra,rbb=rn.get((A,ia),"UNK"),rn.get((B,ib),"UNK")
            ha=hmap.get((qmeta["chainA_protein"],ia),"")
            hb=hmap.get((qmeta["chainB_protein"],ib),"")
            ca=charge_center(ra,aa_atoms); cb=charge_center(rbb,bb_atoms)
            cd=float(np.linalg.norm(ca-cb)) if ca is not None and cb is not None else float("nan")
            charge_cat=""
            if ra in NEG3|POS3 and rbb in NEG3|POS3:
                opp=(ra in NEG3 and rbb in POS3) or (ra in POS3 and rbb in NEG3)
                charge_cat="opposite" if opp else "same"
                if not math.isnan(cd) and cd<6.0:
                    if opp: eg_opposite_close += int(ha in "eg" and hb in "eg")
                    else: eg_same_close += int(ha in "eg" and hb in "eg")
            if d<5.0 and ha in "ad" and hb in "ad": core_contacts += 1
            contacts.append({"query":qmeta["query"],"pdb":str(pdb),"chainA":A,"resA":ia,"resnA":ra,"heptadA":ha,
                             "chainB":B,"resB":ib,"resnB":rbb,"heptadB":hb,"min_heavy_A":atA,"min_heavy_B":atB,
                             "min_heavy_distance":d,"charge_center_distance":cd,"charge_category":charge_cat})
    score={}
    sj=find_score_json(pdb)
    if sj:
        try: score=json.loads(sj.read_text())
        except: score={}
    iptm=score.get("iptm",score.get("ipTM",float("nan")))
    ptm=score.get("ptm",score.get("pTM",float("nan")))
    # Some versions use 0-1, some presentation layers use 0-100. Keep raw.
    return {
      "query":qmeta["query"],"pdb":str(pdb),"rank":rank_from_name(pdb),"score_json":str(sj) if sj else "",
      "chainA":A,"chainB":B,"chainA_protein":qmeta["chainA_protein"],"chainB_protein":qmeta["chainB_protein"],
      "n_chainA_res":len(chains[A]),"n_chainB_res":len(chains[B]),"axis_dot":orient,"orientation":orientation,
      "helix_like_fraction_A":helix_like_fraction(chains[A]),"helix_like_fraction_B":helix_like_fraction(chains[B]),
      "interchain_residue_contacts_lt5A":residue_pairs,"ad_ad_contacts_lt5A":core_contacts,
      "eg_opposite_charge_centers_lt6A":eg_opposite_close,"eg_same_charge_centers_lt6A":eg_same_close,
      "iptm_raw":iptm,"ptm_raw":ptm,"status":"OK"
    },contacts


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("project_root",type=Path)
    args=ap.parse_args(); root=args.project_root.resolve()
    out=root/"results/phase82c_modeled_structure"; pred=out/"colabfold_main"
    contract=json.loads((out/"MODEL_INPUT_CONTRACT.json").read_text())
    hmap=load_heptad_map(out/"sequence_heptad_map.csv")
    qmeta={
      "main_target":{"query":"main_target","chainA_protein":contract["shared_member"],"chainB_protein":contract["target_model_order"][1]},
      "main_competitor":{"query":"main_competitor","chainA_protein":contract["shared_member"],"chainB_protein":contract["competitor_model_order"][1]},
    }
    pdbs=sorted(pred.rglob("*.pdb"))
    # Prefer unrelaxed if both relaxed and unrelaxed exist, to avoid duplicate counting.
    unrel=[p for p in pdbs if "unrelaxed" in p.name]
    if unrel: pdbs=unrel
    rows=[]; contacts=[]
    for p in pdbs:
        q=query_from_name(p)
        if not q: continue
        rr,cc=analyze_one(p,qmeta[q],hmap); rows.append(rr); contacts.extend(cc)
    if not rows: raise SystemExit(f"No main_target/main_competitor PDB files found under {pred}")
    m=pd.DataFrame(rows).sort_values(["query","rank","pdb"])
    m.to_csv(out/"model_geometry_by_prediction.csv",index=False)
    pd.DataFrame(contacts).to_csv(out/"interchain_contacts_by_prediction.csv",index=False)
    selected=[]; qsum=[]
    for q,g in m[m.status.eq("OK")].groupby("query"):
        gg=g.sort_values(["rank","pdb"])
        # Use top five ranked distinct predictions when available.
        top=gg.head(5)
        mode=Counter(top.orientation).most_common(1)[0][0]
        consistency=float((top.orientation==mode).mean())
        helix=float(np.nanmedian(np.minimum(top.helix_like_fraction_A,top.helix_like_fraction_B)))
        contacts_med=float(np.nanmedian(top.interchain_residue_contacts_lt5A))
        best=gg.iloc[0]
        selected.append(best)
        qsum.append({"query":q,"n_models":len(gg),"top5_orientation_mode":mode,"top5_orientation_consistency":consistency,
                     "top5_median_min_helix_like_fraction":helix,"top5_median_contacts_lt5A":contacts_med,
                     "rank1_pdb":best.pdb,"rank1_iptm_raw":best.iptm_raw,"rank1_ptm_raw":best.ptm_raw,
                     "rank1_ad_ad_contacts_lt5A":best.ad_ad_contacts_lt5A,"rank1_eg_opposite_charge_centers_lt6A":best.eg_opposite_charge_centers_lt6A,
                     "rank1_eg_same_charge_centers_lt6A":best.eg_same_charge_centers_lt6A})
    qs=pd.DataFrame(qsum); qs.to_csv(out/"model_ensemble_summary.csv",index=False)
    pd.DataFrame(selected).to_csv(out/"selected_rank1_models.csv",index=False)
    # Geometry gate intentionally does NOT call these structures experimentally validated.
    ready=True
    reasons=[]
    needed={"main_target","main_competitor"}
    if set(qs["query"])!=needed: ready=False; reasons.append("missing_target_or_competitor_models")
    for r in qs.itertuples(index=False):
        if r.top5_orientation_consistency < 0.6: ready=False; reasons.append(f"{r.query}:orientation_inconsistent")
        if r.top5_median_min_helix_like_fraction < 0.6: ready=False; reasons.append(f"{r.query}:weak_alpha_helix_geometry")
        if r.top5_median_contacts_lt5A < 8: ready=False; reasons.append(f"{r.query}:weak_interchain_contact")
    gate={"phase":"Phase-8.2C modeled structure analysis",
          "status":"QUALITATIVE_MODELED_STRUCTURE_READY_FOR_INTERFACE_ANNOTATION" if ready else "REVIEW_MODEL_GEOMETRY_BEFORE_FIGURE",
          "experimental_structure":False,"modeling_backend":"ColabFold/AlphaFold2-multimer as run by user",
          "wording":"Use modeled/predicted structural case study only; do not call these experimental structures.",
          "reasons":reasons,"ensemble_summary":qsum}
    gate=json_safe(gate)
    (out/"MODELED_STRUCTURE_GATE.json").write_text(json.dumps(gate,indent=2,allow_nan=False)+"\n")
    print(json.dumps(gate,indent=2,allow_nan=False))
    print("\n[ENSEMBLE SUMMARY]")
    print(qs.to_string(index=False))
    print(f"[OK] wrote {out}")

if __name__=="__main__": main()
