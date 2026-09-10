#!/usr/bin/env python3
"""
HeptadPair Phase-8.2A v1.2.0 - locked structural-case selector.

This hotfix deliberately avoids generic schema guessing. It reconstructs the
Section-3.4 partner-retrieval cases from the frozen Phase-7.4B outputs:
  * set_retrieval_metrics.csv
  * set_partner_retrieval.csv
  * family_repeat_averaged_predictions.parquet

No model is retrained and no score is retuned. Intended target pairs come from
set_partner_retrieval; candidate-pair scores come from the already-written
locked prediction table. The main case is the correctly recovered heterodimer
with the smallest positive margin to the strongest off-target sharing either
partner. A clear-margin successful case and an optional failure case are also
reported.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

VERSION = "1.2.0"
PREFERRED_VARIANTS = [
    "xgb_x2_1250_final_shared_2434",
    "hgb_c2_final_shared_2434",
    "xgb_x2_1250_final_dual_2434",
    "hgb_c2_final_dual_2434",
]

AA_RE = re.compile(r"[^A-Z]")


def norm_seq(x) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return AA_RE.sub("", str(x).upper().replace("U", "X"))


def read_table(path: Path, nrows: Optional[int] = None) -> pd.DataFrame:
    suf = "".join(path.suffixes).lower()
    if suf.endswith(".csv") or suf.endswith(".csv.gz"):
        return pd.read_csv(path, nrows=nrows, low_memory=False)
    if suf.endswith(".tsv") or suf.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", nrows=nrows, low_memory=False)
    if suf.endswith(".parquet"):
        df = pd.read_parquet(path)
        return df.head(nrows) if nrows else df
    if suf.endswith(".xlsx") or suf.endswith(".xls"):
        df = pd.read_excel(path)
        return df.head(nrows) if nrows else df
    raise ValueError(f"Unsupported table: {path}")


def require(path: Path, label: str) -> Path:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")
    return path


def pair_key(a, b) -> tuple[str, str]:
    a, b = str(a), str(b)
    return (a, b) if a <= b else (b, a)


def as_bool(x) -> bool:
    if pd.isna(x):
        return False
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, float, np.integer, np.floating)):
        return bool(float(x) > 0.5)
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y", "pass", "exact"}


def pick_variant(metrics: pd.DataFrame, partners: pd.DataFrame, preds: pd.DataFrame, requested: Optional[str]) -> str:
    common = set(metrics["variant"].dropna().astype(str))
    common &= set(partners["variant"].dropna().astype(str))
    common &= set(preds["variant"].dropna().astype(str))
    if not common:
        raise RuntimeError("No variant is shared by Phase-7.4B metrics, partner-retrieval and prediction tables")
    if requested:
        if requested not in common:
            raise RuntimeError(f"Requested variant {requested!r} is not present in all three Phase-7.4B tables. Shared={sorted(common)}")
        return requested
    for v in PREFERRED_VARIANTS:
        if v in common:
            return v
    # Fallback: maximize mean AP, then exact recovery, while staying within shared variants.
    mm = metrics[metrics.variant.astype(str).isin(common)].copy()
    score = mm.groupby("variant", as_index=False).agg(
        average_precision=("average_precision", "mean"),
        exact_topk_recovery=("exact_topk_recovery", "mean"),
    )
    score = score.sort_values(["average_precision", "exact_topk_recovery", "variant"], ascending=[False, False, True])
    return str(score.iloc[0].variant)


def build_targets(partners: pd.DataFrame, variant: str) -> tuple[dict[str, set[tuple[str, str]]], dict[str, set[str]], pd.DataFrame]:
    p = partners[partners.variant.astype(str).eq(variant)].copy()
    need = {"official_set_id", "protein", "intended_partner"}
    miss = need - set(p.columns)
    if miss:
        raise RuntimeError(f"Partner-retrieval table missing columns: {sorted(miss)}")
    target_map: dict[str, set[tuple[str, str]]] = {}
    member_map: dict[str, set[str]] = {}
    rows = []
    for sid, g in p.groupby("official_set_id", sort=False):
        sid = str(sid)
        targets = {pair_key(a, b) for a, b in zip(g.protein.astype(str), g.intended_partner.astype(str))}
        members = set()
        for a, b in targets:
            members.add(a); members.add(b)
        target_map[sid] = targets
        member_map[sid] = members
        rows.append({
            "official_set_id": sid,
            "variant": variant,
            "target_pair_count_from_partner_table": len(targets),
            "member_count_from_partner_table": len(members),
            "target_pairs": " | ".join(f"{a}::{b}" for a, b in sorted(targets)),
            "member_proteins": " | ".join(sorted(members)),
        })
    return target_map, member_map, pd.DataFrame(rows)


def normalize_predictions(preds: pd.DataFrame, variant: str) -> pd.DataFrame:
    need = {"variant", "seq_low_id", "seq_high_id", "y_pred"}
    miss = need - set(preds.columns)
    if miss:
        raise RuntimeError(f"Prediction table missing columns: {sorted(miss)}")
    p = preds[preds.variant.astype(str).eq(variant)].copy()
    p["score"] = pd.to_numeric(p["y_pred"], errors="coerce")
    p = p.dropna(subset=["score", "seq_low_id", "seq_high_id"])
    p["protein_a"] = p["seq_low_id"].astype(str)
    p["protein_b"] = p["seq_high_id"].astype(str)
    keys = [pair_key(a, b) for a, b in zip(p.protein_a, p.protein_b)]
    p["pair_a"] = [k[0] for k in keys]
    p["pair_b"] = [k[1] for k in keys]
    # Average only if the frozen output contains repeated views of exactly the same unordered pair.
    agg_spec = {
        "score": ("score", "mean"),
        "score_sd_across_written_rows": ("score", "std"),
        "written_prediction_rows": ("score", "size"),
    }
    if "interaction_score" in p.columns:
        p["interaction_score_num"] = pd.to_numeric(p["interaction_score"], errors="coerce")
        agg_spec["interaction_score"] = ("interaction_score_num", "mean")
    out = p.groupby(["pair_a", "pair_b"], as_index=False).agg(**agg_spec)
    out["score_sd_across_written_rows"] = out["score_sd_across_written_rows"].fillna(0.0)
    return out


def reconstruct_cases(metrics: pd.DataFrame, partners: pd.DataFrame, preds: pd.DataFrame, variant: str):
    target_map, member_map, target_inventory = build_targets(partners, variant)
    ps = normalize_predictions(preds, variant)
    m = metrics[metrics.variant.astype(str).eq(variant)].copy()
    if "official_set_id" not in m.columns or "exact_topk_recovery" not in m.columns:
        raise RuntimeError("Set metric table must contain official_set_id and exact_topk_recovery")
    m["official_set_id"] = m.official_set_id.astype(str)
    metric_by_set = {str(r.official_set_id): r for _, r in m.iterrows()}

    set_rows = []
    case_rows = []
    pair_rows = []
    for sid in sorted(target_map):
        targets = target_map[sid]
        members = member_map[sid]
        gm = ps[ps.pair_a.isin(members) & ps.pair_b.isin(members)].copy()
        gm = gm.drop_duplicates(["pair_a", "pair_b"])
        if gm.empty:
            set_rows.append({"official_set_id": sid, "variant": variant, "status": "NO_CANDIDATE_PREDICTIONS"})
            continue
        gm["is_target"] = [pair_key(a,b) in targets for a,b in zip(gm.pair_a, gm.pair_b)]
        target_rows = gm[gm.is_target].copy()
        off_rows = gm[~gm.is_target].copy()
        K = len(targets)
        target_covered = {pair_key(a,b) for a,b in zip(target_rows.pair_a, target_rows.pair_b)}
        all_targets_scored = target_covered == targets
        ranked = gm.sort_values(["score", "pair_a", "pair_b"], ascending=[False, True, True])
        derived_exact = bool(all_targets_scored and K > 0 and len(ranked) >= K and ranked.head(K).is_target.all())
        metric_row = metric_by_set.get(sid)
        reported_exact = as_bool(metric_row.exact_topk_recovery) if metric_row is not None else None
        exact_match = (reported_exact == derived_exact) if reported_exact is not None else None
        weakest_target = float(target_rows.score.min()) if len(target_rows) else np.nan
        strongest_off = float(off_rows.score.max()) if len(off_rows) else np.nan
        set_gap = weakest_target - strongest_off if np.isfinite(weakest_target) and np.isfinite(strongest_off) else np.nan
        set_rows.append({
            "official_set_id": sid,
            "variant": variant,
            "status": "OK" if all_targets_scored else "TARGET_PREDICTION_COVERAGE_INCOMPLETE",
            "member_count": len(members),
            "target_count": K,
            "target_count_scored": len(target_covered),
            "candidate_pair_count": len(gm),
            "reported_exact_topk_recovery": reported_exact,
            "derived_exact_topk_recovery": derived_exact,
            "reported_vs_derived_exact_match": exact_match,
            "weakest_target_score": weakest_target,
            "strongest_offtarget_score": strongest_off,
            "set_orthogonality_gap": set_gap,
            "reported_average_precision": float(metric_row.average_precision) if metric_row is not None and hasattr(metric_row, "average_precision") else np.nan,
            "reported_partner_top1_accuracy": float(metric_row.partner_top1_accuracy) if metric_row is not None and hasattr(metric_row, "partner_top1_accuracy") else np.nan,
        })
        for _, rr in gm.iterrows():
            pair_rows.append({"official_set_id": sid, "variant": variant, **rr.to_dict()})
        for a, b in sorted(targets):
            tr = gm[(gm.pair_a.eq(a)) & (gm.pair_b.eq(b))]
            if tr.empty:
                continue
            tscore = float(tr.iloc[0].score)
            off = off_rows[((off_rows.pair_a.eq(a)) | (off_rows.pair_b.eq(a)) | (off_rows.pair_a.eq(b)) | (off_rows.pair_b.eq(b)))].copy()
            shared = True
            if off.empty:
                off = off_rows.copy(); shared = False
            if off.empty:
                continue
            c = off.sort_values(["score", "pair_a", "pair_b"], ascending=[False, True, True]).iloc[0]
            case_rows.append({
                "official_set_id": sid,
                "variant": variant,
                "reported_exact_topk_recovery": reported_exact,
                "derived_exact_topk_recovery": derived_exact,
                "reported_vs_derived_exact_match": exact_match,
                "target_a": a,
                "target_b": b,
                "target_is_heterodimer": a != b,
                "target_score": tscore,
                "competitor_a": str(c.pair_a),
                "competitor_b": str(c.pair_b),
                "competitor_score": float(c.score),
                "margin": float(tscore - float(c.score)),
                "competitor_shares_target_member": shared,
                "reported_average_precision": float(metric_row.average_precision) if metric_row is not None and hasattr(metric_row, "average_precision") else np.nan,
            })
    return pd.DataFrame(pair_rows), pd.DataFrame(set_rows), pd.DataFrame(case_rows), target_inventory


def choose_cases(cases: pd.DataFrame) -> pd.DataFrame:
    if cases.empty:
        return pd.DataFrame()
    chosen = []
    # Main/backup require independently reconstructed exact recovery and a positive target-specific margin.
    success = cases[(cases.derived_exact_topk_recovery == True) & (cases.margin > 0)].copy()
    sh = success[success.target_is_heterodimer == True]
    pool = sh if not sh.empty else success
    if not pool.empty:
        chosen.append({"case_role": "MAIN_success_tight_competitor", **pool.sort_values(["margin", "official_set_id"]).iloc[0].to_dict()})
        backup = pool.sort_values(["margin", "official_set_id"], ascending=[False, True]).iloc[0]
        chosen.append({"case_role": "BACKUP_success_clear_margin", **backup.to_dict()})
    failure = cases[(cases.derived_exact_topk_recovery == False) | (cases.margin <= 0)].copy()
    fh = failure[failure.target_is_heterodimer == True]
    fpool = fh if not fh.empty else failure
    if not fpool.empty:
        chosen.append({"case_role": "OPTIONAL_failure_boundary", **fpool.sort_values(["margin", "official_set_id"]).iloc[0].to_dict()})
    out = pd.DataFrame(chosen)
    # Avoid duplicate main and backup rows when only one successful target is available.
    if len(out) >= 2:
        keycols = ["official_set_id", "target_a", "target_b", "competitor_a", "competitor_b"]
        if out.iloc[0][keycols].astype(str).tolist() == out.iloc[1][keycols].astype(str).tolist():
            out.loc[out.index[1], "case_role"] = "BACKUP_same_as_main_only_success_available"
    return out


def col_lookup(cols: Iterable[str], names: Iterable[str]) -> Optional[str]:
    lower = {str(c).lower(): c for c in cols}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def add_pair_sequence_map(df: pd.DataFrame, seqmap: dict[str, str]) -> int:
    cols = list(df.columns)
    layouts = [
        (["seq1_id", "protein1_id", "protein1", "id1", "name1"], ["seq2_id", "protein2_id", "protein2", "id2", "name2"],
         ["seq1", "sequence1", "sequence_1", "protein1_sequence", "sequence_a"], ["seq2", "sequence2", "sequence_2", "protein2_sequence", "sequence_b"]),
        (["seq_low_id", "protein_low_id"], ["seq_high_id", "protein_high_id"], ["seq_low", "sequence_low", "seq_low_sequence"], ["seq_high", "sequence_high", "seq_high_sequence"]),
    ]
    before = len(seqmap)
    for ac, bc, s1c, s2c in layouts:
        a = col_lookup(cols, ac); b = col_lookup(cols, bc); s1 = col_lookup(cols, s1c); s2 = col_lookup(cols, s2c)
        if all(x is not None for x in (a,b,s1,s2)):
            for aa, bb, x1, x2 in df[[a,b,s1,s2]].dropna().itertuples(index=False, name=None):
                q1, q2 = norm_seq(x1), norm_seq(x2)
                if q1: seqmap.setdefault(str(aa), q1)
                if q2: seqmap.setdefault(str(bb), q2)
    # Single-protein table fallback.
    ic = col_lookup(cols, ["protein", "protein_id", "seq_id", "name", "id", "cc_name"])
    sc = col_lookup(cols, ["sequence", "seq", "aa_sequence", "protein_sequence"])
    if ic is not None and sc is not None:
        for ii, ss in df[[ic,sc]].dropna().itertuples(index=False, name=None):
            q = norm_seq(ss)
            if q: seqmap.setdefault(str(ii), q)
    return len(seqmap) - before


def attach_sequences(selected: pd.DataFrame, root: Path, outdir: Path):
    if selected.empty:
        return selected, pd.DataFrame(), []
    needed = set()
    for c in ["target_a","target_b","competitor_a","competitor_b"]:
        needed.update(selected[c].dropna().astype(str))
    preferred = [
        root / "data/external/ccmax/ccmax_pairs.parquet",
        root / "data/external/ccmax/ccmax_pairs.csv",
        root / "results/phase7/ccmax/ccmax_frozen_features.parquet",
        root / "results/phase74/ccmax_supervised/ccmax_pair_augmented_features.parquet",
    ]
    candidates = [p for p in preferred if p.exists()]
    for base in [root/"data/external/ccmax", root/"data/external/ngb2h_official", root/"results/phase7/ccmax", root/"results/phase74/ccmax_supervised"]:
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if p.is_file() and "".join(p.suffixes).lower() in {".csv",".csv.gz",".tsv",".tsv.gz",".parquet",".xlsx",".xls"} and p not in candidates:
                candidates.append(p)
    seqmap: dict[str,str] = {}
    audit = []
    used = []
    for p in candidates:
        if needed.issubset(seqmap.keys()):
            break
        try:
            df = read_table(p)
            added = add_pair_sequence_map(df, seqmap)
            audit.append({"path":str(p), "read_ok":True, "rows":len(df), "sequences_added":added, "needed_mapped_after":len(needed & set(seqmap))})
            if added:
                used.append(str(p))
        except Exception as e:
            audit.append({"path":str(p), "read_ok":False, "error":repr(e), "sequences_added":0, "needed_mapped_after":len(needed & set(seqmap))})
    out = selected.copy()
    for c in ["target_a","target_b","competitor_a","competitor_b"]:
        out[c+"_sequence"] = out[c].astype(str).map(seqmap)
        out[c+"_sequence_length"] = out[c+"_sequence"].fillna("").map(len)
    pd.DataFrame(audit).to_csv(outdir/"sequence_source_audit.csv", index=False)
    missing = sorted(needed - set(seqmap))
    return out, pd.DataFrame(audit), missing


def write_fasta(selected: pd.DataFrame, path: Path):
    seen = set()
    with path.open("w") as fh:
        for _, r in selected.iterrows():
            for base in ["target_a","target_b","competitor_a","competitor_b"]:
                name = str(r[base]); seq = norm_seq(r.get(base+"_sequence", ""))
                key = (name, seq)
                if seq and key not in seen:
                    seen.add(key)
                    fh.write(f">{name}\n{seq}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--variant", default=None, help="Optional exact Phase-7.4B variant override")
    ap.add_argument("--metrics", type=Path, default=None)
    ap.add_argument("--partners", type=Path, default=None)
    ap.add_argument("--predictions", type=Path, default=None)
    args = ap.parse_args()

    root = args.project_root.expanduser().resolve()
    out = (args.outdir or root/"results/phase82_structure_case_preflight").resolve()
    out.mkdir(parents=True, exist_ok=True)

    metrics_path = require((args.metrics or root/"results/phase74b/locked_confirmation/family_stability/set_retrieval_metrics.csv").resolve(), "Phase-7.4B set retrieval metrics")
    partners_path = require((args.partners or root/"results/phase74b/locked_confirmation/family_stability/set_partner_retrieval.csv").resolve(), "Phase-7.4B set partner retrieval")
    preds_path = require((args.predictions or root/"results/phase74b/locked_confirmation/family_stability/family_repeat_averaged_predictions.parquet").resolve(), "Phase-7.4B family-averaged predictions")

    metrics = read_table(metrics_path)
    partners = read_table(partners_path)
    preds = read_table(preds_path)
    for label, df, cols in [
        ("metrics", metrics, ["variant","official_set_id","exact_topk_recovery","average_precision"]),
        ("partners", partners, ["variant","official_set_id","protein","intended_partner"]),
        ("predictions", preds, ["variant","seq_low_id","seq_high_id","y_pred"]),
    ]:
        miss = [c for c in cols if c not in df.columns]
        if miss:
            raise SystemExit(f"{label} table missing required columns {miss}")

    variant = pick_variant(metrics, partners, preds, args.variant)
    pair_scores, set_audit, cases, target_inventory = reconstruct_cases(metrics, partners, preds, variant)
    selected = choose_cases(cases)

    pair_scores.to_csv(out/"normalized_pair_scores_by_set.csv", index=False)
    set_audit.to_csv(out/"derived_set_recovery_audit.csv", index=False)
    cases.to_csv(out/"all_target_vs_competitor_cases.csv", index=False)
    target_inventory.to_csv(out/"target_pair_inventory.csv", index=False)

    selected, sequence_audit, missing_sequences = attach_sequences(selected, root, out)
    selected.to_csv(out/"selected_structure_cases.csv", index=False)
    if not selected.empty:
        write_fasta(selected, out/"selected_case_sequences.fasta")

    exact_matches = set_audit[set_audit.get("reported_vs_derived_exact_match", pd.Series(dtype=bool)) == True] if not set_audit.empty else pd.DataFrame()
    exact_mismatches = set_audit[set_audit.get("reported_vs_derived_exact_match", pd.Series(dtype=bool)) == False] if not set_audit.empty else pd.DataFrame()
    if selected.empty:
        status = "REVIEW_NO_CASES_DERIVED"
    elif missing_sequences:
        status = "STRUCTURE_CASES_SELECTED_SEQUENCE_MAPPING_INCOMPLETE"
    elif not exact_mismatches.empty:
        # Selected cases may still be valid if their own set reproduces, but flag any broader mismatch for audit.
        sel_sets = set(selected.official_set_id.astype(str))
        bad_sel = exact_mismatches[exact_mismatches.official_set_id.astype(str).isin(sel_sets)]
        status = "REVIEW_SELECTED_SET_RECOVERY_MISMATCH" if not bad_sel.empty else "STRUCTURE_CASES_SELECTED_WITH_NONSELECTED_SET_MISMATCHES"
    else:
        status = "STRUCTURE_CASES_SELECTED_READY_FOR_PDB_SCREEN"

    manifest = {
        "phase": "HeptadPair Phase-8.2A structural case selection",
        "version": VERSION,
        "status": status,
        "project_root": str(root),
        "variant": variant,
        "source_contract": {
            "metrics": str(metrics_path),
            "partners": str(partners_path),
            "predictions": str(preds_path),
            "retrained_model": False,
            "retuned_model": False,
            "selection_uses_existing_locked_predictions_only": True,
        },
        "set_count_audited": int(len(set_audit)),
        "reported_vs_derived_exact_matches": int(len(exact_matches)),
        "reported_vs_derived_exact_mismatches": int(len(exact_mismatches)),
        "candidate_target_cases": int(len(cases)),
        "selected_case_count": int(len(selected)),
        "missing_selected_protein_sequences": missing_sequences,
        "selection_policy": {
            "main": "derived exact-recovered heterodimer; smallest positive margin to strongest off-target sharing a target member",
            "backup": "derived exact-recovered heterodimer; largest positive margin",
            "failure": "non-exact or non-positive-margin heterodimer; smallest margin",
        },
        "interpretation": "Case selection is retrospective visualization of a frozen Section-3.4 retrieval result; it is not an additional performance benchmark.",
    }
    (out/"STRUCTURE_CASE_PREFLIGHT.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    if not selected.empty:
        show = [c for c in ["case_role","official_set_id","variant","target_a","target_b","target_score","competitor_a","competitor_b","competitor_score","margin"] if c in selected.columns]
        print("\n[SELECTED CASES]")
        print(selected[show].to_string(index=False))
    print(f"[OK] wrote {out}")
    if status.startswith("REVIEW_") or status.endswith("INCOMPLETE"):
        raise SystemExit(3)


if __name__ == "__main__":
    main()
