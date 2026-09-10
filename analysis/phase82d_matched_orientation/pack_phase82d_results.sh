#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="$ROOT/HeptadPair_Phase82D_MatchedOrientation_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" -C "$ROOT" results/phase82d_matched_orientation \
  2>/dev/null
sha256sum "$OUT" | tee "$OUT.sha256"
echo "[OK] $OUT"
echo "[OK] $OUT.sha256"
