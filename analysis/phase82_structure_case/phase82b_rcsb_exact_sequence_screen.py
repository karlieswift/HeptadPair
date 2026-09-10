#!/usr/bin/env python3
"""
HeptadPair Phase-8.2B v1.2.0 - conservative experimental RCSB PDB screen.

The RCSB Search API v2 sequence query uses sequence_type='protein', a 100%
identity cutoff, and experimental-only results. Every search hit is then
checked against the RCSB Data API canonical polymer-entity sequence. A common
PDB entry is only a candidate; biological-assembly and direct-interface
inspection remains mandatory before calling it an experimental dimer structure.
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

VERSION = "1.2.0"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/{entity}"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{entry}"
AA_RE = re.compile(r"[^A-Z]")


def norm_seq(x: Any) -> str:
    return AA_RE.sub("", str(x).upper())


def get_json(session: requests.Session, url: str, *, payload=None, timeout=40, retries=4):
    last = None
    for i in range(retries):
        try:
            r = session.get(url, timeout=timeout) if payload is None else session.post(url, json=payload, timeout=timeout)
            if r.status_code == 204:
                return None
            if r.status_code == 200:
                return r.json()
            last = RuntimeError(f"HTTP {r.status_code}: {r.text[:700]}")
            if r.status_code not in (429, 500, 502, 503, 504):
                break
        except Exception as e:
            last = e
        time.sleep(1.0 * (2 ** i))
    if last:
        raise last
    return None


def sequence_query_payload(seq: str, rows: int = 500) -> dict:
    # RCSB Search API v2 / rcsbsearchapi.SequenceQuery contract.
    return {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "value": seq,
                "sequence_type": "protein",
                "identity_cutoff": 1.0,
                "evalue_cutoff": 0.1,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "results_content_type": ["experimental"],
            "results_verbosity": "verbose",
            "paginate": {"start": 0, "rows": rows},
        },
    }


def search_100pct(session: requests.Session, seq: str, rows: int = 500) -> list[dict]:
    # Page through the complete 100%-identity experimental hit set so a valid
    # common entry cannot be missed simply because more than one page exists.
    all_hits = []
    start = 0
    while True:
        payload = sequence_query_payload(seq, rows)
        payload["request_options"]["paginate"]["start"] = start
        js = get_json(session, SEARCH_URL, payload=payload)
        if not js:
            break
        batch = list(js.get("result_set", []))
        all_hits.extend(batch)
        total = int(js.get("total_count", len(all_hits)) or len(all_hits))
        if not batch or len(all_hits) >= total:
            break
        start += len(batch)
    return all_hits


def entity_sequence(entity_js: dict) -> str:
    ep = entity_js.get("entity_poly", {}) or {}
    for key in ("pdbx_seq_one_letter_code_can", "pdbx_seq_one_letter_code"):
        if ep.get(key):
            return norm_seq(ep[key])
    return ""


def flatten_match_context(hit: dict) -> str:
    contexts = []
    def walk(x):
        if isinstance(x, dict):
            for k, v in x.items():
                if k in {"match_context", "matching_context"}:
                    contexts.append(v)
                else:
                    walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(hit)
    return json.dumps(contexts, ensure_ascii=False, sort_keys=True, default=str)


def entry_meta(session: requests.Session, entry: str) -> dict:
    try:
        js = get_json(session, ENTRY_URL.format(entry=entry)) or {}
    except Exception as e:
        return {"pdb_id": entry, "entry_meta_error": repr(e)}
    info = js.get("rcsb_entry_info", {}) or {}
    exptl = js.get("exptl", []) or []
    methods = sorted({str(x.get("method")) for x in exptl if x.get("method")})
    res = info.get("resolution_combined")
    if isinstance(res, list):
        res = min(res) if res else None
    return {
        "pdb_id": entry,
        "experimental_methods": " | ".join(methods),
        "resolution_A": res,
        "deposited_polymer_entity_instance_count": info.get("deposited_polymer_entity_instance_count"),
        "assembly_count": info.get("assembly_count"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--case-dir", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--sleep", type=float, default=0.25)
    ap.add_argument("--max-hits", type=int, default=500)
    args = ap.parse_args()

    root = args.project_root.expanduser().resolve()
    case_dir = (args.case_dir or root/"results/phase82_structure_case_preflight").resolve()
    out = (args.outdir or root/"results/phase82_rcsb_screen").resolve()
    out.mkdir(parents=True, exist_ok=True)
    case_path = case_dir/"selected_structure_cases.csv"
    if not case_path.exists():
        raise SystemExit(f"Missing {case_path}; Phase-8.2A must finish successfully first")
    cases = pd.read_csv(case_path)
    if cases.empty:
        raise SystemExit("selected_structure_cases.csv is empty")
    required = ["target_a","target_b","competitor_a","competitor_b"]
    for base in required:
        if base not in cases.columns or base+"_sequence" not in cases.columns:
            raise SystemExit(f"Missing {base} or {base}_sequence in {case_path}")
        if cases[base+"_sequence"].fillna("").astype(str).str.len().eq(0).any():
            raise SystemExit(f"Empty sequence in {base}_sequence; fix Phase-8.2A mapping before PDB screen")

    session = requests.Session()
    session.headers.update({"User-Agent": "HeptadPair-Phase82B/1.2.0 (+structural-case-audit)"})
    unique = {}
    for _, r in cases.iterrows():
        for base in required:
            name = str(r[base]); seq = norm_seq(r[base+"_sequence"])
            if seq:
                unique[(name,seq)] = None

    hit_rows = []
    raw_dir = out/"raw_search_json"; raw_dir.mkdir(exist_ok=True)
    query_dir = out/"query_payloads"; query_dir.mkdir(exist_ok=True)
    for idx, (name, seq) in enumerate(unique, start=1):
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
        payload = sequence_query_payload(seq, args.max_hits)
        (query_dir/f"{safe}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        if len(seq) < 25:
            hit_rows.append({"protein":name,"query_sequence":seq,"query_length":len(seq),"status":"QUERY_TOO_SHORT_FOR_RCSB_SEQUENCE_SEARCH"})
            continue
        print(f"[RCSB {idx}/{len(unique)}] {name} len={len(seq)}", flush=True)
        try:
            hits = search_100pct(session, seq, args.max_hits)
        except Exception as e:
            hit_rows.append({"protein":name,"query_sequence":seq,"query_length":len(seq),"status":"SEARCH_API_ERROR","error":repr(e)})
            continue
        (raw_dir/f"{safe}.json").write_text(json.dumps(hits, indent=2, ensure_ascii=False), encoding="utf-8")
        if not hits:
            hit_rows.append({"protein":name,"query_sequence":seq,"query_length":len(seq),"status":"NO_100PCT_IDENTITY_EXPERIMENTAL_HIT"})
            time.sleep(args.sleep); continue
        for h in hits:
            ident = str(h.get("identifier", ""))
            if "_" not in ident:
                hit_rows.append({"protein":name,"query_sequence":seq,"query_length":len(seq),"status":"UNEXPECTED_POLYMER_ENTITY_IDENTIFIER","identifier":ident})
                continue
            entry, entity = ident.split("_", 1)
            try:
                ej = get_json(session, ENTITY_URL.format(entry=entry, entity=entity)) or {}
                canon = entity_sequence(ej)
                exact_entity = bool(canon) and canon == seq
                exact_subsequence = bool(canon) and seq in canon
                status = "EXACT_ENTITY_SEQUENCE" if exact_entity else ("EXACT_QUERY_SUBSEQUENCE_IN_ENTITY" if exact_subsequence else "100PCT_IDENTITY_ALIGNMENT_ONLY")
                hit_rows.append({
                    "protein":name,"query_sequence":seq,"query_length":len(seq),
                    "pdb_id":entry.upper(),"polymer_entity_id":entity,
                    "entity_sequence":canon,"entity_length":len(canon),
                    "exact_entity_sequence":exact_entity,"exact_query_subsequence_in_entity":exact_subsequence,
                    "status":status,"search_match_context":flatten_match_context(h),
                })
            except Exception as e:
                hit_rows.append({"protein":name,"query_sequence":seq,"query_length":len(seq),"pdb_id":entry.upper(),"polymer_entity_id":entity,"status":"ENTITY_FETCH_ERROR","error":repr(e),"search_match_context":flatten_match_context(h)})
            time.sleep(args.sleep)

    hits = pd.DataFrame(hit_rows)
    hits.to_csv(out/"rcsb_sequence_hits.csv", index=False)
    if hits.empty or "status" not in hits.columns:
        evidence = pd.DataFrame()
    else:
        evidence = hits[hits.status.isin(["EXACT_ENTITY_SEQUENCE","EXACT_QUERY_SUBSEQUENCE_IN_ENTITY"])].copy()

    pair_rows = []; common_entries = set()
    for _, r in cases.iterrows():
        role = str(r.get("case_role", "case"))
        for pair_kind, aa, bb in [("target","target_a","target_b"),("competitor","competitor_a","competitor_b")]:
            a,b = str(r[aa]),str(r[bb])
            ha = evidence[evidence.protein.astype(str).eq(a)] if not evidence.empty else pd.DataFrame()
            hb = evidence[evidence.protein.astype(str).eq(b)] if not evidence.empty else pd.DataFrame()
            ea = set(ha.pdb_id.astype(str)) if not ha.empty else set(); eb = set(hb.pdb_id.astype(str)) if not hb.empty else set()
            common = sorted(ea & eb)
            if not common:
                pair_rows.append({"case_role":role,"pair_kind":pair_kind,"protein_a":a,"protein_b":b,"pdb_id":None,"both_exact_entity_sequences":False,"sequence_evidence":"NO_COMMON_EXPERIMENTAL_PDB_ENTRY_WITH_FULL_QUERY_SEQUENCE","assembly_gate":"NOT_REACHED"})
                continue
            for pdb in common:
                ah = ha[ha.pdb_id.astype(str).eq(pdb)]; bh = hb[hb.pdb_id.astype(str).eq(pdb)]
                both_entity = bool(ah.exact_entity_sequence.fillna(False).any() and bh.exact_entity_sequence.fillna(False).any())
                common_entries.add(pdb)
                pair_rows.append({"case_role":role,"pair_kind":pair_kind,"protein_a":a,"protein_b":b,"pdb_id":pdb,"both_exact_entity_sequences":both_entity,"sequence_evidence":"BOTH_EXACT_ENTITY" if both_entity else "BOTH_FULL_QUERY_EXACT_SUBSEQUENCE","assembly_gate":"NEEDS_BIOLOGICAL_ASSEMBLY_AND_DIRECT_INTERFACE_CHECK"})
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(out/"selected_pair_common_pdb_candidates.csv", index=False)
    meta = pd.DataFrame([entry_meta(session, e) for e in sorted(common_entries)])
    meta.to_csv(out/"common_entry_metadata.csv", index=False)

    query_errors = int(hits.status.eq("SEARCH_API_ERROR").sum()) if not hits.empty and "status" in hits.columns else 0
    target_candidates = pair_df[(pair_df.pair_kind == "target") & pair_df.pdb_id.notna()] if not pair_df.empty else pd.DataFrame()
    competitor_candidates = pair_df[(pair_df.pair_kind == "competitor") & pair_df.pdb_id.notna()] if not pair_df.empty else pd.DataFrame()
    if query_errors:
        status = "REVIEW_RCSB_API_ERRORS"
    elif not target_candidates.empty:
        status = "EXPERIMENTAL_PDB_SEQUENCE_CANDIDATES_FOUND_NEEDS_ASSEMBLY_INTERFACE_GATE"
    else:
        status = "NO_DIRECT_EXPERIMENTAL_PDB_SEQUENCE_CANDIDATE_USE_MODELED_STRUCTURE_IF_PROCEEDING"
    manifest = {
        "phase":"HeptadPair Phase-8.2B RCSB exact-sequence screen",
        "version":VERSION,
        "status":status,
        "selected_cases":int(len(cases)),
        "unique_protein_sequences":int(len(unique)),
        "search_api_errors":query_errors,
        "target_common_entry_candidates":int(len(target_candidates)),
        "competitor_common_entry_candidates":int(len(competitor_candidates)),
        "common_pdb_entries":sorted(common_entries),
        "hard_gate":"A common PDB entry is not sufficient. Verify both chains in the same relevant biological assembly and a direct coiled-coil interface before using experimental-structure wording.",
        "query_contract":{"service":"sequence","sequence_type":"protein","identity_cutoff":1.0,"evalue_cutoff":0.1,"results_content_type":["experimental"],"return_type":"polymer_entity"},
    }
    (out/"RCSB_STRUCTURE_SCREEN_GATE.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(f"[OK] wrote {out}")
    if query_errors:
        raise SystemExit(4)


if __name__ == "__main__":
    main()
