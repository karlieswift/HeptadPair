#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
PACK="$(cd "$(dirname "$0")" && pwd)"
OUT="$ROOT/results/phase82c_modeled_structure"
python "$PACK/phase82c_make_pymol_scripts.py" "$ROOT"
if command -v pymol >/dev/null 2>&1; then
  for p in "$OUT"/pymol/*.pml; do
    echo "[PYMOL] $p"
    pymol -cq "$p"
  done
else
  echo "[INFO] pymol not found; .pml scripts were generated but not rendered."
fi
