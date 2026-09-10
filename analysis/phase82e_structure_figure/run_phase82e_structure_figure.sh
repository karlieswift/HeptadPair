#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
PY="${PYTHON:-python}"
cd "$ROOT"
"$PY" phase82e_structure_figure_pack/phase82e_audit_controls.py "$ROOT"
# Prefer PyMOL's Python entry so the pymol module is definitely available.
if command -v pymol >/dev/null 2>&1; then
  PHASE82E_ROOT="$ROOT" PHASE82E_PYMOL_AUTORUN=1 pymol -cq phase82e_structure_figure_pack/phase82e_render_matched_views.py
else
  echo '[FAIL] pymol not found in PATH; activate bioh200 and rerun.' >&2
  exit 2
fi

# Hard gate: never enter composition unless all four matched-view renders exist.
OUT="$ROOT/results/phase82e_structure_figure"
REQ=(
  target_matched_view.png
  hard_matched_view.png
  strict_matched_view.png
  hardpar_matched_view.png
  matched_view_render_manifest.csv
)
for f in "${REQ[@]}"; do
  if [[ ! -s "$OUT/$f" ]]; then
    echo "[FAIL] matched-view renderer did not create: $OUT/$f" >&2
    echo "[HINT] PyMOL was found, but the Python renderer did not complete. Do not run the composer on missing renders." >&2
    exit 3
  fi
done
echo "[PASS] matched-view render gate: ${#REQ[@]} required outputs present"

"$PY" phase82e_structure_figure_pack/phase82e_compose_structure_panels.py "$ROOT"
echo "[DONE] $ROOT/results/phase82e_structure_figure"
