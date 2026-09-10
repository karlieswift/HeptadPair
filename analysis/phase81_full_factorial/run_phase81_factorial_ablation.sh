#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
PACK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${PHASE81_OUTDIR:-$ROOT/results/phase81_factorial_ablation_v1_0_2}"
REPLAY_TOL="${PHASE81_REPLAY_TOL:-0.005}"
BOOTSTRAP="${PHASE81_BOOTSTRAP:-20000}"

cd "$ROOT"

for p in \
  results/features/heptadpair_mvp_features.parquet \
  data/splits/phase67_final_locked_both_chains_unseen_seed67.csv \
  data/splits/phase67_final_locked_both_chains_unseen_seed68.csv \
  data/splits/phase67_final_locked_both_chains_unseen_seed69.csv \
  data/splits/phase67_final_locked_both_chains_unseen_seed70.csv \
  data/splits/phase67_final_locked_both_chains_unseen_seed71.csv \
  results/phase67/final_frozen/model_summary.csv; do
  if [[ ! -f "$p" ]]; then
    echo "[FAIL] missing required project file: $ROOT/$p" >&2
    exit 2
  fi
done

# The local Python package is not necessarily installed into the active conda env.
# Historical checkouts may use ROOT/heptadpair or ROOT/src/heptadpair. Discover it
# from this exact project root and prepend the package parent to PYTHONPATH.
PHASE67_SOURCE=""
for candidate in \
  "$ROOT/heptadpair/phase67.py" \
  "$ROOT/src/heptadpair/phase67.py" \
  "$ROOT/python/heptadpair/phase67.py"; do
  if [[ -f "$candidate" ]]; then
    PHASE67_SOURCE="$candidate"
    break
  fi
done

if [[ -z "$PHASE67_SOURCE" ]]; then
  PHASE67_SOURCE="$(find "$ROOT" -maxdepth 6 -type f -path '*/heptadpair/phase67.py' \
    ! -path '*/results/*' ! -path '*/data/*' ! -path '*/.git/*' \
    ! -path '*/heptadpair_factorial_ablation_pack/*' -print 2>/dev/null | sort | head -n 1 || true)"
fi

if [[ -z "$PHASE67_SOURCE" ]]; then
  echo "[FAIL] could not locate heptadpair/phase67.py under $ROOT" >&2
  echo "[DEBUG] candidate package directories:" >&2
  find "$ROOT" -maxdepth 6 -type d -name heptadpair -print 2>/dev/null | head -n 30 >&2 || true
  exit 2
fi

PACKAGE_PARENT="$(cd "$(dirname "$(dirname "$PHASE67_SOURCE")")" && pwd)"
export PYTHONPATH="$PACKAGE_PARENT${PYTHONPATH:+:$PYTHONPATH}"

echo "[PHASE81 v1.0.2] root=$ROOT"
echo "[PHASE81 v1.0.2] python=$(command -v python)"
echo "[PHASE81 v1.0.2] phase67_source=$PHASE67_SOURCE"
echo "[PHASE81 v1.0.2] PYTHONPATH_head=$PACKAGE_PARENT"

python - <<'PY'
import sys
import heptadpair.phase67 as p67
print(f"[IMPORT PREFLIGHT PASS] python={sys.executable}")
print(f"[IMPORT PREFLIGHT PASS] heptadpair.phase67={p67.__file__}")
print(f"[IMPORT PREFLIGHT PASS] KEYS={list(p67.KEYS)}")
PY

python "$PACK_DIR/run_factorial_ablation.py" \
  --root "$ROOT" \
  --features "$ROOT/results/features/heptadpair_mvp_features" \
  --local "$ROOT/results/features/two_center_topology_features_v064" \
  --ordered "$ROOT/results/features/ordered_two_center_topology_features_v065" \
  --pairs "$ROOT/data/processed/pairs" \
  --split-dir "$ROOT/data/splits" \
  --reference-summary "$ROOT/results/phase67/final_frozen/model_summary.csv" \
  --outdir "$OUTDIR" \
  --bootstrap "$BOOTSTRAP" \
  --replay-tol "$REPLAY_TOL"

echo
printf '[DONE] Phase81 v1.0.2 outputs: %s\n' "$OUTDIR"
printf '[NEXT] bash %s/pack_results.sh %q\n' "$PACK_DIR" "$ROOT"
