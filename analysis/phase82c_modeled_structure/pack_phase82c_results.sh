#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
NAME="HeptadPair_Phase82C_ModeledStructure_RESULTS_${STAMP}.tar.gz"
cd "$ROOT"
tar -czf "$NAME" \
  results/phase82c_modeled_structure \
  run_phase82c0_prepare_and_backend_preflight.log \
  run_phase82c1_colabfold_main.log \
  run_phase82c2_analyze_models.log \
  run_phase82c3_make_pymol.log 2>/dev/null || \
tar -czf "$NAME" results/phase82c_modeled_structure
sha256sum "$NAME" > "$NAME.sha256"
echo "[OK] $ROOT/$NAME"
echo "[OK] $ROOT/$NAME.sha256"
