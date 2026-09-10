#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
OUTDIR="${PHASE81_OUTDIR:-$ROOT/results/phase81_factorial_ablation_v1_0_2}"
STAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE="$ROOT/HeptadPair_Phase81_FullFactorial_Ablation_RESULTS_v1_0_2_${STAMP}.tar.gz"
if [[ ! -d "$OUTDIR" ]]; then
  echo "[FAIL] output directory not found: $OUTDIR" >&2
  exit 2
fi
tar -czf "$ARCHIVE" -C "$(dirname "$OUTDIR")" "$(basename "$OUTDIR")"
sha256sum "$ARCHIVE" | tee "$ARCHIVE.sha256"
echo "[DONE] $ARCHIVE"
