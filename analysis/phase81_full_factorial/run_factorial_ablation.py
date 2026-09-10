#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import itertools
import importlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ROOT_DEFAULT = Path.cwd()

BLOCK_LABELS = {
    "E": "engineered",
    "O": "ordered-register",
    "L": "local-register",
}

# All 2^3 combinations. The empty model is needed for exact factorial/Shapley decomposition.
SUBSETS: List[Tuple[str, ...]] = [
    tuple(),
    ("E",),
    ("O",),
    ("L",),
    ("E", "O"),
    ("E", "L"),
    ("O", "L"),
    ("E", "O", "L"),
]

VARIANT_NAME = {
    tuple(): "EMPTY_intercept_only_0",
    ("E",): "E_engineered_2214",
    ("O",): "O_ordered_118",
    ("L",): "L_local_102",
    ("E", "O"): "EO_engineered_plus_ordered_2332",
    ("E", "L"): "EL_engineered_plus_local_2316",
    ("O", "L"): "OL_ordered_plus_local_220",
    ("E", "O", "L"): "EOL_full_2434",
}

EXPECTED_BLOCK_COUNTS = {"E": 2214, "O": 118, "L": 102}
EXPECTED_SUBSET_COUNTS = {
    tuple(): 0,
    ("E",): 2214,
    ("O",): 118,
    ("L",): 102,
    ("E", "O"): 2332,
    ("E", "L"): 2316,
    ("O", "L"): 220,
    ("E", "O", "L"): 2434,
}

PRIMARY_METRIC = "spearman"
METRICS = ("pearson", "spearman", "r2", "rmse", "mae")
LOWER_IS_BETTER = {"rmse", "mae"}


def subset_key(subset: Sequence[str]) -> str:
    return "".join(subset) if subset else "EMPTY"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_table(path: Path) -> Path:
    """Resolve either an exact file path or a stem used by the HeptadPair project."""
    path = Path(path)
    if path.is_file():
        return path
    candidates = []
    if path.suffix:
        candidates.append(path)
    else:
        candidates.extend([
            path.with_suffix(".parquet"),
            path.with_suffix(".pkl"),
            path.with_suffix(".pickle"),
            path.with_suffix(".csv"),
        ])
    existing = [p for p in candidates if p.is_file()]
    if len(existing) == 1:
        return existing[0]
    if len(existing) > 1:
        # The project routinely uses stems with a parquet primary and occasional preview CSV.
        parquet = [p for p in existing if p.suffix.lower() == ".parquet"]
        if len(parquet) == 1:
            return parquet[0]
        raise RuntimeError(f"Ambiguous table stem {path}: {existing}")
    raise FileNotFoundError(f"Could not resolve table: {path}")


def read_table(path: Path) -> pd.DataFrame:
    actual = resolve_table(path)
    s = actual.suffix.lower()
    if s == ".parquet":
        return pd.read_parquet(actual)
    if s in {".pkl", ".pickle"}:
        return pd.read_pickle(actual)
    if s == ".csv":
        return pd.read_csv(actual)
    raise ValueError(f"Unsupported table format: {actual}")


def safe_corr(func, y: np.ndarray, pred: np.ndarray) -> float:
    if len(y) < 2 or np.nanstd(y) == 0 or np.nanstd(pred) == 0:
        return float("nan")
    try:
        return float(func(y, pred)[0])
    except Exception:
        return float("nan")


def regression_metrics(y: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    mask = np.isfinite(y) & np.isfinite(pred)
    y = y[mask]
    pred = pred[mask]
    if len(y) == 0:
        return {m: float("nan") for m in METRICS}
    return {
        "pearson": safe_corr(pearsonr, y, pred),
        "spearman": safe_corr(spearmanr, y, pred),
        "r2": float(r2_score(y, pred)) if len(y) >= 2 else float("nan"),
        "rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "mae": float(mean_absolute_error(y, pred)),
    }


def orient_gain(metric: str, delta: float) -> float:
    return -float(delta) if metric in LOWER_IS_BETTER else float(delta)


def find_partition_column(split: pd.DataFrame) -> str:
    preferred = ["partition", "split", "subset", "role", "set"]
    lower = {str(c).lower(): c for c in split.columns}
    for name in preferred:
        if name in lower:
            col = lower[name]
            values = set(split[col].dropna().astype(str).str.lower().unique())
            if values & {"train", "valid", "validation", "test", "discard_cross_partition"}:
                return col
    # Last resort: search any low-cardinality string column containing train/test.
    for col in split.columns:
        vals = split[col].dropna().astype(str).str.lower()
        if len(vals) and vals.nunique() <= 10:
            values = set(vals.unique())
            if "train" in values and "test" in values:
                return col
    raise RuntimeError(
        "Could not identify split partition column. Columns=" + ",".join(map(str, split.columns))
    )


def normalize_partition(value: object) -> str:
    s = str(value).strip().lower()
    if s == "validation":
        return "valid"
    return s


def choose_join_keys(keys_from_project: Sequence[str], split: pd.DataFrame, data: pd.DataFrame) -> List[str]:
    candidates = [k for k in keys_from_project if k in split.columns and k in data.columns]
    if candidates:
        return candidates
    fallback = [
        "pair_id", "pair_key", "id", "name",
        "seq1", "seq2", "sequence1", "sequence2",
        "protein1", "protein2", "protein_1", "protein_2",
    ]
    candidates = [k for k in fallback if k in split.columns and k in data.columns]
    # Prefer a compact identifier over two sequence columns when available.
    for single in ["pair_id", "pair_key", "id", "name"]:
        if single in candidates:
            return [single]
    seq_pairs = [
        ["seq1", "seq2"], ["sequence1", "sequence2"],
        ["protein1", "protein2"], ["protein_1", "protein_2"],
    ]
    for pair in seq_pairs:
        if all(x in candidates for x in pair):
            return pair
    if candidates:
        return candidates
    raise RuntimeError(
        "Could not identify join keys between split and merged CCNG1 table. "
        f"split columns={list(split.columns)[:30]}, data columns={list(data.columns)[:30]}"
    )


def merge_ccng1_tables(
    pairs: pd.DataFrame,
    engineered: pd.DataFrame,
    local: pd.DataFrame,
    ordered: pd.DataFrame,
    project_keys: Sequence[str],
    e_cols: Sequence[str],
    l_cols: Sequence[str],
    o_cols: Sequence[str],
    target_col: str,
) -> Tuple[pd.DataFrame, List[str]]:
    # Determine pair identity from the frozen project contract first.
    keys = [k for k in project_keys if all(k in frame.columns for frame in [pairs, engineered, local, ordered])]
    if not keys:
        keys = [k for k in project_keys if k in pairs.columns and k in engineered.columns]
    if not keys:
        common = set(pairs.columns) & set(engineered.columns) & set(local.columns) & set(ordered.columns)
        preferred = ["pair_id", "pair_key", "seq1", "seq2", "sequence1", "sequence2", "protein1", "protein2"]
        keys = [k for k in preferred if k in common]
        if "pair_id" in keys:
            keys = ["pair_id"]
        elif "pair_key" in keys:
            keys = ["pair_key"]
        elif all(k in keys for k in ["seq1", "seq2"]):
            keys = ["seq1", "seq2"]
    if not keys:
        raise RuntimeError(
            "Could not determine feature merge keys from heptadpair.phase67.KEYS or common columns. "
            f"project KEYS={list(project_keys)}"
        )

    def prep(frame: pd.DataFrame, cols: Sequence[str], label: str) -> pd.DataFrame:
        missing = [c for c in [*keys, *cols] if c not in frame.columns]
        if missing:
            raise RuntimeError(f"{label} table missing columns: {missing[:20]}")
        out = frame[[*keys, *cols]].copy()
        if out.duplicated(keys).any():
            dup = out.loc[out.duplicated(keys, keep=False), keys].head(10).to_dict("records")
            raise RuntimeError(f"{label} feature table has duplicate keys: {dup}")
        return out

    metadata_cols = list(keys)
    for c in [target_col, "seq1", "seq2", "sequence1", "sequence2", "protein1", "protein2"]:
        if c in pairs.columns and c not in metadata_cols:
            metadata_cols.append(c)
    base = pairs[metadata_cols].copy()
    if base.duplicated(keys).any():
        dup = base.loc[base.duplicated(keys, keep=False), keys].head(10).to_dict("records")
        raise RuntimeError(f"processed pairs table has duplicate keys: {dup}")

    merged = base.merge(prep(engineered, e_cols, "engineered"), on=keys, how="left", validate="one_to_one")
    merged = merged.merge(prep(local, l_cols, "local"), on=keys, how="left", validate="one_to_one")
    merged = merged.merge(prep(ordered, o_cols, "ordered"), on=keys, how="left", validate="one_to_one")

    if target_col not in merged.columns:
        # Some project versions keep target in the engineered table rather than data/processed/pairs.
        if target_col in engineered.columns:
            t = engineered[[*keys, target_col]].drop_duplicates(keys)
            merged = merged.merge(t, on=keys, how="left", validate="one_to_one")
        else:
            raise RuntimeError(f"Target column {target_col!r} not found in processed pairs or engineered features")

    all_features = [*e_cols, *l_cols, *o_cols]
    numeric = merged[all_features].apply(pd.to_numeric, errors="coerce")
    coverage = float(numeric.notna().all(axis=1).mean())
    if coverage < 0.95:
        bad = int((~numeric.notna().all(axis=1)).sum())
        raise RuntimeError(f"Final feature merge coverage below 95%: coverage={coverage:.6f}, bad_rows={bad}")
    merged = merged.loc[numeric.notna().all(axis=1)].copy().reset_index(drop=True)
    merged[all_features] = numeric.loc[numeric.notna().all(axis=1)].to_numpy(np.float32)
    merged[target_col] = pd.to_numeric(merged[target_col], errors="coerce")
    if not np.isfinite(merged[target_col].to_numpy(float)).all():
        raise RuntimeError("Non-finite CCNG1 target values after feature merge")
    return merged, keys


def exact_signflip_p(values: Sequence[float], alternative: str = "two-sided") -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    obs = float(np.mean(x))
    stats = []
    for signs in itertools.product([-1.0, 1.0], repeat=len(x)):
        stats.append(float(np.mean(x * np.asarray(signs))))
    stats = np.asarray(stats)
    if alternative == "greater":
        return float(np.mean(stats >= obs - 1e-15))
    if alternative == "less":
        return float(np.mean(stats <= obs + 1e-15))
    return float(np.mean(np.abs(stats) >= abs(obs) - 1e-15))


def bootstrap_mean_ci(values: Sequence[float], seed: int, n_boot: int = 20000) -> Tuple[float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    means = x[idx].mean(axis=1)
    lo, hi = np.quantile(means, [0.025, 0.975])
    return float(lo), float(hi)


def mean_std(values: Sequence[float]) -> Tuple[float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan"), float("nan")
    return float(np.mean(x)), float(np.std(x, ddof=1)) if len(x) > 1 else 0.0


def summarize_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (variant, subset, n_features), g in metrics.groupby(["variant", "subset", "n_features"], sort=False):
        row = {"variant": variant, "subset": subset, "n_features": int(n_features), "outer_splits": int(g["split_seed"].nunique())}
        for metric in METRICS:
            m, s = mean_std(g[metric])
            row[f"{metric}_mean"] = m
            row[f"{metric}_std"] = s
        rows.append(row)
    out = pd.DataFrame(rows)
    order = {subset_key(s): i for i, s in enumerate(SUBSETS)}
    out["_order"] = out["subset"].map(order)
    return out.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def leave_one_out_effects(metrics: pd.DataFrame, n_boot: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    full = metrics.loc[metrics["subset"].eq("EOL")].set_index("split_seed")
    without = {"E": "OL", "O": "EL", "L": "EO"}
    rows = []
    for block, subset in without.items():
        sub = metrics.loc[metrics["subset"].eq(subset)].set_index("split_seed")
        common = full.index.intersection(sub.index)
        for seed in common:
            for metric in METRICS:
                raw = float(full.loc[seed, metric] - sub.loc[seed, metric])
                rows.append({
                    "split_seed": int(seed),
                    "block": block,
                    "block_name": BLOCK_LABELS[block],
                    "comparison": f"EOL_minus_{subset}",
                    "metric": metric,
                    "raw_delta": raw,
                    "oriented_gain": orient_gain(metric, raw),
                })
    long = pd.DataFrame(rows)
    summary_rows = []
    for (block, block_name, comparison, metric), g in long.groupby(["block", "block_name", "comparison", "metric"], sort=False):
        vals = g["oriented_gain"].to_numpy(float)
        lo, hi = bootstrap_mean_ci(vals, seed=810000 + ord(block) + sum(map(ord, metric)), n_boot=n_boot)
        summary_rows.append({
            "block": block,
            "block_name": block_name,
            "comparison": comparison,
            "metric": metric,
            "n_splits": int(len(vals)),
            "oriented_gain_mean": float(np.mean(vals)),
            "oriented_gain_std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
            "bootstrap_ci95_low": lo,
            "bootstrap_ci95_high": hi,
            "fraction_positive": float(np.mean(vals > 0)),
            "exact_signflip_two_sided_p": exact_signflip_p(vals, "two-sided"),
            "exact_signflip_one_sided_greater_p": exact_signflip_p(vals, "greater"),
        })
    summary = pd.DataFrame(summary_rows)
    return long, summary


def standalone_effects(metrics: pd.DataFrame) -> pd.DataFrame:
    empty = metrics.loc[metrics["subset"].eq("EMPTY")].set_index("split_seed")
    rows = []
    for block in ["E", "O", "L"]:
        sub = metrics.loc[metrics["subset"].eq(block)].set_index("split_seed")
        common = empty.index.intersection(sub.index)
        for metric in METRICS:
            raw_vals = [float(sub.loc[s, metric] - empty.loc[s, metric]) for s in common]
            oriented = [orient_gain(metric, x) for x in raw_vals]
            rows.append({
                "block": block,
                "block_name": BLOCK_LABELS[block],
                "metric": metric,
                "standalone_oriented_gain_mean": float(np.nanmean(oriented)),
                "standalone_oriented_gain_std": float(np.nanstd(oriented, ddof=1)) if len(oriented) > 1 else 0.0,
                "fraction_positive": float(np.nanmean(np.asarray(oriented) > 0)),
            })
    return pd.DataFrame(rows)


def shapley_by_split(metrics: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Exact 3-player Shapley decomposition.

    For Spearman/Pearson, the intercept-only model is constant and therefore has undefined
    correlation. We use v(empty)=0 as an explicit no-signal convention for correlation-only
    Shapley reporting. The raw factorial metric table retains NaN for the empty-model correlation.
    R2/RMSE/MAE use the actual intercept-only predictions without a convention.
    """
    rows = []
    n = 3
    all_blocks = ("E", "O", "L")
    metric_lookup = {}
    for _, r in metrics.iterrows():
        metric_lookup[(int(r.split_seed), str(r.subset))] = r

    for seed in sorted(metrics["split_seed"].unique()):
        for metric in METRICS:
            value = {}
            for subset in SUBSETS:
                k = subset_key(subset)
                v = float(metric_lookup[(int(seed), k)][metric])
                if k == "EMPTY" and metric in {"pearson", "spearman"} and not np.isfinite(v):
                    v = 0.0
                # Convert error metrics to utility so positive Shapley = useful contribution.
                if metric in LOWER_IS_BETTER:
                    v = -v
                value[frozenset(subset)] = v

            for block in all_blocks:
                phi = 0.0
                others = [b for b in all_blocks if b != block]
                for r in range(len(others) + 1):
                    for S_tuple in itertools.combinations(others, r):
                        S = frozenset(S_tuple)
                        weight = math.factorial(len(S)) * math.factorial(n - len(S) - 1) / math.factorial(n)
                        phi += weight * (value[S | {block}] - value[S])
                rows.append({
                    "split_seed": int(seed),
                    "metric": metric,
                    "block": block,
                    "block_name": BLOCK_LABELS[block],
                    "shapley_utility": float(phi),
                    "correlation_empty_zero_convention": bool(metric in {"pearson", "spearman"}),
                })
    long = pd.DataFrame(rows)
    summary_rows = []
    for (metric, block, block_name), g in long.groupby(["metric", "block", "block_name"], sort=False):
        vals = g["shapley_utility"].to_numpy(float)
        summary_rows.append({
            "metric": metric,
            "block": block,
            "block_name": block_name,
            "n_splits": int(len(vals)),
            "shapley_mean": float(np.mean(vals)),
            "shapley_std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
            "fraction_positive": float(np.mean(vals > 0)),
            "correlation_empty_zero_convention": bool(metric in {"pearson", "spearman"}),
        })
    return long, pd.DataFrame(summary_rows)


def interaction_decomposition(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    lookup = {(int(r.split_seed), str(r.subset)): r for _, r in metrics.iterrows()}
    for seed in sorted(metrics["split_seed"].unique()):
        for metric in METRICS:
            def v(k: str) -> float:
                x = float(lookup[(int(seed), k)][metric])
                if k == "EMPTY" and metric in {"pearson", "spearman"} and not np.isfinite(x):
                    x = 0.0
                if metric in LOWER_IS_BETTER:
                    x = -x
                return x
            i_eo = v("EO") - v("E") - v("O") + v("EMPTY")
            i_el = v("EL") - v("E") - v("L") + v("EMPTY")
            i_ol = v("OL") - v("O") - v("L") + v("EMPTY")
            i_eol = (
                v("EOL") - v("EO") - v("EL") - v("OL")
                + v("E") + v("O") + v("L") - v("EMPTY")
            )
            for name, val in [("E:O", i_eo), ("E:L", i_el), ("O:L", i_ol), ("E:O:L", i_eol)]:
                rows.append({"split_seed": int(seed), "metric": metric, "interaction": name, "utility_interaction": float(val)})
    return pd.DataFrame(rows)


def replay_audit(current_summary: pd.DataFrame, reference: pd.DataFrame, tolerance: float) -> Tuple[pd.DataFrame, bool]:
    mapping = {
        "E": "engineered_hgb_2214",
        "EO": "engineered_plus_ordered_hgb_2332",
        "EOL": "final_local_plus_ordered_hgb_2434",
    }
    rows = []
    pass_all = True
    for subset, ref_variant in mapping.items():
        cur = current_summary.loc[current_summary["subset"].eq(subset)]
        ref = reference.loc[reference["variant"].astype(str).eq(ref_variant)]
        if len(cur) != 1 or len(ref) != 1:
            rows.append({"subset": subset, "reference_variant": ref_variant, "status": "MISSING_ROW"})
            pass_all = False
            continue
        c = cur.iloc[0]
        r = ref.iloc[0]
        row = {"subset": subset, "reference_variant": ref_variant, "status": "PASS"}
        for metric in METRICS:
            ccol = f"{metric}_mean"
            rcol = f"{metric}_mean"
            if rcol not in ref.columns:
                continue
            delta = float(c[ccol] - r[rcol])
            row[f"current_{ccol}"] = float(c[ccol])
            row[f"reference_{ccol}"] = float(r[rcol])
            row[f"delta_{ccol}"] = delta
            if abs(delta) > tolerance:
                row["status"] = "FAIL_NUMERIC_REPLAY"
                pass_all = False
        rows.append(row)
    return pd.DataFrame(rows), pass_all


def make_figures(summary: pd.DataFrame, loo_summary: pd.DataFrame, outdir: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        (outdir / "FIGURE_SKIPPED.txt").write_text(f"matplotlib unavailable: {e}\n", encoding="utf-8")
        return

    # Figure 1: full factorial Spearman means.
    order = [subset_key(s) for s in SUBSETS if s]
    sub = summary[summary["subset"].isin(order)].copy()
    sub["order"] = sub["subset"].map({k: i for i, k in enumerate(order)})
    sub = sub.sort_values("order")
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.bar(sub["subset"], sub["spearman_mean"], yerr=sub["spearman_std"], capsize=3)
    ax.set_xlabel("Feature-block subset")
    ax.set_ylabel("Test Spearman (mean +/- SD across 5 locked splits)")
    ax.set_title("Phase-8.1 full-factorial ablation")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "Figure_phase81_full_factorial_spearman.png", dpi=240)
    fig.savefig(outdir / "Figure_phase81_full_factorial_spearman.pdf")
    plt.close(fig)

    # Figure 2: leave-one-block-out contribution to the full model.
    s = loo_summary[loo_summary["metric"].eq("spearman")].copy()
    s = s.set_index("block").reindex(["E", "O", "L"]).reset_index()
    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.bar(s["block_name"], s["oriented_gain_mean"], yerr=s["oriented_gain_std"], capsize=3)
    ax.axhline(0.0, linewidth=1)
    ax.set_ylabel("Full-model Spearman loss when block is removed")
    ax.set_title("Conditional contribution within the 2,434-feature model")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "Figure_phase81_leave_one_block_out_spearman.png", dpi=240)
    fig.savefig(outdir / "Figure_phase81_leave_one_block_out_spearman.pdf")
    plt.close(fig)



def bootstrap_heptadpair_source(root: Path) -> Path:
    """Locate the project's local heptadpair source tree and put its parent on sys.path.

    Historical HeptadPair checkouts may store the package either as
    ``ROOT/heptadpair`` or ``ROOT/src/heptadpair``.  The Phase81 pack lives in a
    sibling directory, so relying on Python's script-directory sys.path is not
    sufficient.  This bootstrap deliberately resolves the source from the
    requested project root rather than requiring an editable pip install.
    """
    root = Path(root).resolve()
    direct = [
        root / "heptadpair" / "phase67.py",
        root / "src" / "heptadpair" / "phase67.py",
        root / "python" / "heptadpair" / "phase67.py",
    ]
    matches = [x for x in direct if x.is_file()]

    if not matches:
        excluded = {
            ".git", ".hg", ".svn", "__pycache__", "results", "data",
            "review_pack", "review_packs", "archive", "archives",
            "heptadpair_factorial_ablation_pack",
        }
        discovered = []
        try:
            for candidate in root.rglob("phase67.py"):
                if candidate.parent.name != "heptadpair" or not candidate.is_file():
                    continue
                rel = candidate.relative_to(root)
                if any(part in excluded for part in rel.parts[:-2]):
                    continue
                discovered.append(candidate)
        except OSError:
            discovered = []
        matches = discovered

    if not matches:
        raise RuntimeError(
            "Could not locate local HeptadPair source: expected heptadpair/phase67.py "
            "or src/heptadpair/phase67.py under project root " + str(root)
        )

    # Prefer the shallowest source tree; ties are deterministic.
    matches = sorted(matches, key=lambda x: (len(x.relative_to(root).parts), str(x)))
    phase67_file = matches[0].resolve()
    package_parent = phase67_file.parent.parent
    if str(package_parent) not in sys.path:
        sys.path.insert(0, str(package_parent))
    importlib.invalidate_caches()
    return phase67_file

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="HeptadPair Phase-8.1 exact 2^3 full-factorial feature-block ablation on frozen Phase-6.7 splits")
    p.add_argument("--root", type=Path, default=ROOT_DEFAULT)
    p.add_argument("--features", type=Path, default=None, help="Engineered feature table/stem")
    p.add_argument("--local", type=Path, default=None, help="Phase-6.4 local no-torsion feature table/stem")
    p.add_argument("--ordered", type=Path, default=None, help="Phase-6.5 ordered no-torsion feature table/stem")
    p.add_argument("--pairs", type=Path, default=None, help="Processed CCNG1 pairs table/stem")
    p.add_argument("--split-dir", type=Path, default=None)
    p.add_argument("--reference-summary", type=Path, default=None)
    p.add_argument("--outdir", type=Path, default=None)
    p.add_argument("--target", default="interaction_score")
    p.add_argument("--model-seed-base", type=int, default=810000)
    p.add_argument("--bootstrap", type=int, default=20000)
    p.add_argument("--replay-tol", type=float, default=0.005,
                   help="Maximum absolute difference in mean metric allowed for E/EO/EOL vs frozen Phase-6.7 summary. Default 0.005 is a strict protocol-replay guard while allowing HGB RNG/library micro-differences.")
    p.add_argument("--strict-replay", action="store_true", help="Exit nonzero if E/EO/EOL replay exceeds --replay-tol")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    phase67_source = bootstrap_heptadpair_source(root)
    print(f"[PYTHON SOURCE] phase67={phase67_source}")
    print(f"[PYTHON SOURCE] sys.path[0]={sys.path[0]}")

    try:
        from heptadpair.phase67 import KEYS, make_hgb, select_feature_blocks
    except Exception as e:
        raise RuntimeError(
            "Located local heptadpair/phase67.py but import still failed. "
            f"source={phase67_source}; sys.path[0]={sys.path[0]}; "
            f"python={sys.executable}. Original error: {type(e).__name__}: {e}"
        ) from e

    features_path = args.features or root / "results/features/heptadpair_mvp_features"
    local_path = args.local or root / "results/features/two_center_topology_features_v064"
    ordered_path = args.ordered or root / "results/features/ordered_two_center_topology_features_v065"
    pairs_path = args.pairs or root / "data/processed/pairs"
    split_dir = args.split_dir or root / "data/splits"
    reference_summary_path = args.reference_summary or root / "results/phase67/final_frozen/model_summary.csv"
    outdir = args.outdir or root / "results/phase81_factorial_ablation_v1_0_2"
    outdir.mkdir(parents=True, exist_ok=True)

    actual_paths = {
        "engineered": resolve_table(features_path),
        "local": resolve_table(local_path),
        "ordered": resolve_table(ordered_path),
        "pairs": resolve_table(pairs_path),
        "reference_summary": Path(reference_summary_path),
    }
    for name, pth in actual_paths.items():
        if not pth.is_file():
            raise FileNotFoundError(f"Missing required {name}: {pth}")

    print("[PHASE81 v1.0.2] exact paths")
    for name, pth in actual_paths.items():
        print(f"  {name}: {pth}")
    print(f"  project KEYS: {list(KEYS)}")

    engineered = read_table(actual_paths["engineered"])
    local = read_table(actual_paths["local"])
    ordered = read_table(actual_paths["ordered"])
    pairs = read_table(actual_paths["pairs"])

    blocks = select_feature_blocks(engineered, local, ordered)
    e_cols = list(blocks.engineered)
    l_cols = list(blocks.local_no_torsion)
    o_cols = list(blocks.ordered_no_torsion)
    observed_blocks = {"E": len(e_cols), "L": len(l_cols), "O": len(o_cols)}
    if observed_blocks != {"E": 2214, "L": 102, "O": 118}:
        raise RuntimeError(f"Frozen block-count contract failed: observed={observed_blocks}")

    data, merge_keys = merge_ccng1_tables(
        pairs=pairs,
        engineered=engineered,
        local=local,
        ordered=ordered,
        project_keys=list(KEYS),
        e_cols=e_cols,
        l_cols=l_cols,
        o_cols=o_cols,
        target_col=args.target,
    )
    print(f"[MERGE] rows={len(data):,} keys={merge_keys} target={args.target}")

    block_cols = {"E": e_cols, "O": o_cols, "L": l_cols}
    for subset in SUBSETS:
        n = sum(len(block_cols[b]) for b in subset)
        if n != EXPECTED_SUBSET_COUNTS[subset]:
            raise RuntimeError(f"Subset feature-count contract failed: {subset_key(subset)} expected={EXPECTED_SUBSET_COUNTS[subset]} observed={n}")

    metrics_rows = []
    pred_rows = []
    split_audit_rows = []
    start_time = time.time()

    for split_seed in [67, 68, 69, 70, 71]:
        split_path = split_dir / f"phase67_final_locked_both_chains_unseen_seed{split_seed}.csv"
        if not split_path.is_file():
            raise FileNotFoundError(f"Missing frozen Phase-6.7 split: {split_path}")
        split = pd.read_csv(split_path)
        partition_col = find_partition_column(split)
        split = split.copy()
        split["_partition_norm"] = split[partition_col].map(normalize_partition)
        join_keys = choose_join_keys(list(KEYS), split, data)

        membership_cols = [*join_keys, "_partition_norm"]
        membership = split[membership_cols].drop_duplicates(join_keys)
        if membership.duplicated(join_keys).any():
            raise RuntimeError(f"Split seed {split_seed}: duplicate membership keys after normalization")

        current = data.merge(membership, on=join_keys, how="inner", validate="one_to_one")
        train = current[current["_partition_norm"].eq("train")].copy()
        valid = current[current["_partition_norm"].eq("valid")].copy()
        test = current[current["_partition_norm"].eq("test")].copy()
        if len(train) == 0 or len(test) == 0:
            raise RuntimeError(
                f"Split seed {split_seed}: empty train/test after join; join_keys={join_keys}, partition_col={partition_col}"
            )

        split_audit_rows.append({
            "split_seed": split_seed,
            "split_file": str(split_path),
            "partition_column": str(partition_col),
            "join_keys": ";".join(join_keys),
            "train": len(train),
            "valid": len(valid),
            "test": len(test),
            "discard_cross_partition": int((split["_partition_norm"].eq("discard_cross_partition")).sum()),
            "merged_rows": len(current),
        })
        print(f"[SPLIT {split_seed}] train={len(train)} valid={len(valid)} test={len(test)} join={join_keys}")

        y_train = train[args.target].to_numpy(float)
        y_test = test[args.target].to_numpy(float)
        fit_seed = int(args.model_seed_base + split_seed)

        for subset in SUBSETS:
            skey = subset_key(subset)
            variant = VARIANT_NAME[subset]
            cols = [c for b in subset for c in block_cols[b]]
            if not cols:
                pred = np.full(len(test), float(np.mean(y_train)), dtype=float)
                model_kind = "intercept_only_train_mean"
            else:
                model = make_hgb(fit_seed)
                X_train = train[cols].to_numpy(np.float32, copy=False)
                X_test = test[cols].to_numpy(np.float32, copy=False)
                model.fit(X_train, y_train)
                pred = np.asarray(model.predict(X_test), dtype=float)
                model_kind = type(model).__name__

            met = regression_metrics(y_test, pred)
            metrics_rows.append({
                "split_seed": split_seed,
                "subset": skey,
                "variant": variant,
                "blocks": "+".join(subset) if subset else "none",
                "n_features": len(cols),
                "model_seed": fit_seed,
                "model_kind": model_kind,
                "n_train": len(train),
                "n_valid": len(valid),
                "n_test": len(test),
                **met,
            })
            for i, (_, row) in enumerate(test.iterrows()):
                rec = {
                    "split_seed": split_seed,
                    "subset": skey,
                    "variant": variant,
                    "row_in_test": i,
                    "y_true": float(y_test[i]),
                    "y_pred": float(pred[i]),
                }
                for k in join_keys:
                    rec[k] = row[k]
                pred_rows.append(rec)
            print(
                f"  [{skey:5s}] nfeat={len(cols):4d} "
                f"Pearson={met['pearson']:.6f} Spearman={met['spearman']:.6f} R2={met['r2']:.6f}"
            )

    metrics = pd.DataFrame(metrics_rows)
    predictions = pd.DataFrame(pred_rows)
    split_audit = pd.DataFrame(split_audit_rows)
    summary = summarize_metrics(metrics)

    # Historical replay audit for the three combinations that existed in Phase-6.7.
    reference_summary = pd.read_csv(actual_paths["reference_summary"])
    replay, replay_pass = replay_audit(summary, reference_summary, args.replay_tol)

    loo_long, loo_summary = leave_one_out_effects(metrics, n_boot=args.bootstrap)
    standalone = standalone_effects(metrics)
    shapley_long, shapley_summary = shapley_by_split(metrics)
    interactions = interaction_decomposition(metrics)

    metrics.to_csv(outdir / "factorial_metrics_by_split.csv", index=False)
    summary.to_csv(outdir / "factorial_summary.csv", index=False)
    predictions.to_csv(outdir / "factorial_test_predictions.csv", index=False)
    split_audit.to_csv(outdir / "phase67_split_reuse_audit.csv", index=False)
    replay.to_csv(outdir / "phase67_replay_audit.csv", index=False)
    loo_long.to_csv(outdir / "leave_one_block_out_effects_by_split.csv", index=False)
    loo_summary.to_csv(outdir / "leave_one_block_out_summary.csv", index=False)
    standalone.to_csv(outdir / "standalone_block_summary.csv", index=False)
    shapley_long.to_csv(outdir / "shapley_by_split.csv", index=False)
    shapley_summary.to_csv(outdir / "shapley_summary.csv", index=False)
    interactions.to_csv(outdir / "factorial_interactions_by_split.csv", index=False)

    make_figures(summary, loo_summary, outdir)

    # Primary ranking: conditional loss from removing each block from the full model.
    loo_primary = loo_summary[loo_summary["metric"].eq(PRIMARY_METRIC)].copy()
    loo_primary = loo_primary.sort_values("oriented_gain_mean", ascending=False).reset_index(drop=True)
    primary_ranking = loo_primary[["block", "block_name", "oriented_gain_mean", "bootstrap_ci95_low", "bootstrap_ci95_high", "fraction_positive", "exact_signflip_two_sided_p"]].to_dict("records")

    # Secondary ranking: exact Shapley for primary metric (empty correlation is explicitly defined as zero signal).
    shap_primary = shapley_summary[shapley_summary["metric"].eq(PRIMARY_METRIC)].copy()
    shap_primary = shap_primary.sort_values("shapley_mean", ascending=False).reset_index(drop=True)

    provenance = {
        name: {"path": str(path), "sha256": sha256_file(path)}
        for name, path in actual_paths.items()
    }
    provenance.update({
        f"split_seed_{s}": {
            "path": str(split_dir / f"phase67_final_locked_both_chains_unseen_seed{s}.csv"),
            "sha256": sha256_file(split_dir / f"phase67_final_locked_both_chains_unseen_seed{s}.csv"),
        }
        for s in [67, 68, 69, 70, 71]
    })

    gate = {
        "phase": "HeptadPair Phase-8.1 full-factorial ablation v1.0.2",
        "status": "PASS_FULL_FACTORIAL_ABLATION" if replay_pass else "REVIEW_REPLAY_MISMATCH",
        "replay_pass": bool(replay_pass),
        "replay_tolerance": float(args.replay_tol),
        "feature_blocks": {
            "E_engineered": 2214,
            "O_ordered_register": 118,
            "L_local_register": 102,
        },
        "evaluated_subsets": [subset_key(s) for s in SUBSETS],
        "outer_splits": [67, 68, 69, 70, 71],
        "primary_metric": PRIMARY_METRIC,
        "primary_contribution_definition": "leave-one-block-out loss from the full EOL model; positive means the full model performs better with the block present",
        "primary_block_ranking": primary_ranking,
        "secondary_shapley_ranking_spearman": shap_primary[["block", "block_name", "shapley_mean", "shapley_std", "fraction_positive", "correlation_empty_zero_convention"]].to_dict("records"),
        "shapley_caveat": "For Pearson/Spearman only, the intercept-only model has undefined correlation; v(empty)=0 is used as an explicit no-signal convention. R2/RMSE/MAE Shapley values use the actual intercept-only predictions.",
        "interpretation_rule": "Do not infer block importance from the previous cumulative E -> EO -> EOL increments alone. Use standalone, leave-one-out, full-factorial and Shapley results together.",
        "runtime_seconds": float(time.time() - start_time),
        "provenance": provenance,
    }
    (outdir / "FACTORIAL_ABLATION_GATE.json").write_text(json.dumps(gate, indent=2), encoding="utf-8")

    # Manuscript-ready compact text, generated only from the completed factorial table.
    lines = [
        "Phase-8.1 full-factorial ablation (draft reporting text)",
        "",
        "The three frozen HeptadPair feature blocks were evaluated under all 2^3 combinations on the same five Phase-6.7 both-chains-unseen splits, including an intercept-only reference. This design separates standalone predictive value from conditional contribution within the full model and avoids interpreting a cumulative E->EO->EOL ordering as a block-importance ranking.",
        "",
        f"Replay status: {'PASS' if replay_pass else 'REVIEW_REPLAY_MISMATCH'} (tolerance={args.replay_tol}).",
        "",
        "Primary leave-one-block-out Spearman ranking:",
    ]
    for rank, r in enumerate(primary_ranking, 1):
        lines.append(
            f"{rank}. {r['block_name']} ({r['block']}): mean full-model loss when removed = {r['oriented_gain_mean']:+.6f}; "
            f"95% bootstrap CI {r['bootstrap_ci95_low']:+.6f} to {r['bootstrap_ci95_high']:+.6f}; "
            f"positive in {r['fraction_positive']:.1%} of splits."
        )
    (outdir / "MANUSCRIPT_FACTORIAL_ABLATION_DRAFT.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n[PHASE81 SUMMARY]")
    print(summary[["subset", "n_features", "pearson_mean", "spearman_mean", "r2_mean", "rmse_mean", "mae_mean"]].to_string(index=False))
    print("\n[LEAVE-ONE-OUT PRIMARY SPEARMAN]")
    print(loo_primary[["block", "block_name", "oriented_gain_mean", "bootstrap_ci95_low", "bootstrap_ci95_high", "fraction_positive", "exact_signflip_two_sided_p"]].to_string(index=False))
    print(f"\n[REPLAY] pass={replay_pass} tolerance={args.replay_tol}")
    print(replay.to_string(index=False))
    print(f"\n[DONE] {outdir}")

    if args.strict_replay and not replay_pass:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
