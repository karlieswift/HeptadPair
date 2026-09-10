#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import pandas as pd

LABELS = "abcdefg"
HYDRO = set("AVILMFWY")
CHARGED = set("DEKR")


def infer_a_start(sequences: dict[str, str]):
    rows=[]
    for a_start in range(1,8):
        score=0.0
        for name, seq in sequences.items():
            for pos, aa in enumerate(seq,1):
                lab = LABELS[(pos-a_start)%7]
                if lab in "ad":
                    score += 1.0 if aa in HYDRO else (0.25 if aa in "NQ" else (-1.0 if aa in CHARGED else 0.0))
                if lab in "eg":
                    score += 0.5 if aa in CHARGED else (-0.2 if aa in HYDRO else 0.0)
        rows.append({"a_start_1based":a_start,"phase_score":score})
    rows=sorted(rows,key=lambda r:r["phase_score"],reverse=True)
    best=rows[0]
    best["margin_over_second"] = best["phase_score"]-rows[1]["phase_score"]
    return best, rows


def label_positions(seq: str, a_start: int):
    out=[]
    for pos, aa in enumerate(seq,1):
        lab=LABELS[(pos-a_start)%7]
        role = "hydrophobic_core" if lab in "ad" else ("charge_specificity" if lab in "eg" else "sequence_path")
        out.append({"position":pos,"residue":aa,"heptad":lab,"role":role})
    return out


def charge_sign(aa: str):
    if aa in "KR": return 1
    if aa in "DE": return -1
    return 0


def simple_pair_audit(name: str, seq_a: str, seq_b: str, a_start: int):
    amap={x["position"]:x for x in label_positions(seq_a,a_start)}
    bmap={x["position"]:x for x in label_positions(seq_b,a_start)}
    # a-a and d-d aligned core categories; qualitative only.
    core=[]
    for lab in "ad":
        pa=[p for p,v in amap.items() if v["heptad"]==lab]
        pb=[p for p,v in bmap.items() if v["heptad"]==lab]
        for i,j in zip(pa,pb):
            aa,bb=seq_a[i-1],seq_b[j-1]
            if aa in HYDRO and bb in HYDRO: cat="hydrophobic_hydrophobic"
            elif aa in "NQ" and bb in "NQ": cat="polar_polar_core"
            elif (aa in HYDRO and bb in "NQ") or (bb in HYDRO and aa in "NQ"): cat="hydrophobic_polar_core"
            else: cat="other_core"
            core.append((lab,i,aa,j,bb,cat))
    # Same-heptad e<->g' bookkeeping. This is a sequence-level canonical-register audit, not a structural contact claim.
    eg=[]
    for la,lb in [("e","g"),("g","e")]:
        pa=[p for p,v in amap.items() if v["heptad"]==la]
        pb=[p for p,v in bmap.items() if v["heptad"]==lb]
        for i,j in zip(pa,pb):
            aa,bb=seq_a[i-1],seq_b[j-1]
            prod=charge_sign(aa)*charge_sign(bb)
            cat="opposite_charge" if prod==-1 else ("same_charge" if prod==1 else "neutral")
            eg.append((la+"-"+lb,i,aa,j,bb,cat))
    return {
        "pair": name,
        "core_hydrophobic_hydrophobic": sum(x[-1]=="hydrophobic_hydrophobic" for x in core),
        "core_polar_polar": sum(x[-1]=="polar_polar_core" for x in core),
        "core_hydrophobic_polar": sum(x[-1]=="hydrophobic_polar_core" for x in core),
        "eg_opposite_charge": sum(x[-1]=="opposite_charge" for x in eg),
        "eg_same_charge": sum(x[-1]=="same_charge" for x in eg),
        "eg_neutral": sum(x[-1]=="neutral" for x in eg),
        "note": "Canonical-register sequence audit only; spatial contacts must be checked on modeled coordinates."
    }, core, eg


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--case_csv", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=None)
    args=ap.parse_args()
    root=args.project_root.resolve()
    case_csv=(args.case_csv or root/"results/phase82_structure_case_preflight/selected_structure_cases.csv").resolve()
    out=(args.outdir or root/"results/phase82c_modeled_structure").resolve()
    inp=out/"inputs"; inp.mkdir(parents=True,exist_ok=True)
    if not case_csv.exists(): raise SystemExit(f"Missing {case_csv}")
    df=pd.read_csv(case_csv)
    if df.empty: raise SystemExit("selected_structure_cases.csv is empty")
    mainrow=df[df.case_role.eq("MAIN_success_tight_competitor")]
    if len(mainrow)!=1: raise SystemExit("Expected exactly one MAIN_success_tight_competitor")
    r=mainrow.iloc[0]
    seqs={}
    for prefix in ["target_a","target_b","competitor_a","competitor_b"]:
        seqs[str(r[prefix])]=str(r[prefix+"_sequence"])
    best, scores=infer_a_start(seqs)
    if best["margin_over_second"] < 2.0:
        raise SystemExit(f"Ambiguous heptad phase inference: {best}")
    a_start=int(best["a_start_1based"])
    pd.DataFrame(scores).to_csv(out/"heptad_phase_audit.csv",index=False)

    maps=[]
    for name,seq in seqs.items():
        for rr in label_positions(seq,a_start): maps.append({"protein":name,"sequence":seq,**rr})
    pd.DataFrame(maps).to_csv(out/"sequence_heptad_map.csv",index=False)

    ta,tb=str(r.target_a),str(r.target_b)
    ca,cb=str(r.competitor_a),str(r.competitor_b)
    # Put the shared member first in BOTH complexes to simplify direct structural comparison.
    shared=sorted(set([ta,tb]).intersection([ca,cb]))
    if len(shared)!=1: raise SystemExit(f"Expected exactly one shared member, got {shared}")
    shared=shared[0]
    target_other = tb if ta==shared else ta
    competitor_other = cb if ca==shared else ca
    seq_shared=seqs[shared]; seq_t=seqs[target_other]; seq_c=seqs[competitor_other]
    (inp/"main_target.fasta").write_text(f">main_target__{shared}__{target_other}\n{seq_shared}:{seq_t}\n")
    (inp/"main_competitor.fasta").write_text(f">main_competitor__{shared}__{competitor_other}\n{seq_shared}:{seq_c}\n")
    (inp/"main_pair_batch.fasta").write_text(
        f">main_target__{shared}__{target_other}\n{seq_shared}:{seq_t}\n"
        f">main_competitor__{shared}__{competitor_other}\n{seq_shared}:{seq_c}\n"
    )
    meta={
        "phase":"HeptadPair Phase-8.2C exact-sequence modeled structural case",
        "status":"MODEL_INPUTS_READY",
        "case_role":str(r.case_role),
        "official_set_id":str(r.official_set_id),
        "variant":str(r.variant),
        "target_original":[ta,tb],
        "competitor_original":[ca,cb],
        "shared_member":shared,
        "target_model_order":[shared,target_other],
        "competitor_model_order":[shared,competitor_other],
        "target_score":float(r.target_score),
        "competitor_score":float(r.competitor_score),
        "margin":float(r.margin),
        "reported_average_precision":float(r.reported_average_precision),
        "heptad_a_start_1based":a_start,
        "heptad_phase_score":float(best["phase_score"]),
        "heptad_phase_margin":float(best["margin_over_second"]),
        "structure_wording_policy":"No exact experimental PDB hit exists. Any coordinates produced here must be called modeled/predicted structures, not experimental structures.",
        "modeling_policy":"Primary: AlphaFold2-multimer-v3 via ColabFold, single-sequence/no-template, multiple seeds/models; structural illustration only, not a new performance benchmark."
    }
    (out/"MODEL_INPUT_CONTRACT.json").write_text(json.dumps(meta,indent=2)+"\n")

    audits=[]; core_rows=[]; eg_rows=[]
    for pname,x,y in [("main_target",seq_shared,seq_t),("main_competitor",seq_shared,seq_c)]:
        aud, core, eg=simple_pair_audit(pname,x,y,a_start); audits.append(aud)
        for t in core: core_rows.append({"pair":pname,"contact_class":t[0],"pos_a":t[1],"aa_a":t[2],"pos_b":t[3],"aa_b":t[4],"category":t[5]})
        for t in eg: eg_rows.append({"pair":pname,"contact_class":t[0],"pos_a":t[1],"aa_a":t[2],"pos_b":t[3],"aa_b":t[4],"category":t[5]})
    pd.DataFrame(audits).to_csv(out/"canonical_register_pair_audit.csv",index=False)
    pd.DataFrame(core_rows).to_csv(out/"canonical_core_pairs.csv",index=False)
    pd.DataFrame(eg_rows).to_csv(out/"canonical_eg_pairs.csv",index=False)
    print(json.dumps(meta,indent=2))
    print("\n[CANONICAL REGISTER SEQUENCE AUDIT]")
    print(pd.DataFrame(audits).to_string(index=False))
    print(f"[OK] wrote {out}")

if __name__=="__main__": main()
