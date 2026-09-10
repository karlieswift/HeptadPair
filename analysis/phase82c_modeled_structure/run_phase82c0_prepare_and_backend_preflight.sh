#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
PACK="$(cd "$(dirname "$0")" && pwd)"
OUT="$ROOT/results/phase82c_modeled_structure"
mkdir -p "$OUT"
python "$PACK/phase82c_prepare_inputs.py" "$ROOT"
{
  echo "[PHASE82C BACKEND PREFLIGHT]"
  echo "date=$(date -Is)"
  echo "root=$ROOT"
  echo "python=$(command -v python || true)"
  python --version 2>&1 || true
  echo "colabfold_batch=$(command -v colabfold_batch || true)"
  if command -v colabfold_batch >/dev/null 2>&1; then
    colabfold_batch --help 2>&1 | sed -n '1,120p'
  fi
  echo "pymol=$(command -v pymol || true)"
  echo "chimerax=$(command -v chimerax || true)"
  echo "nvidia-smi=$(command -v nvidia-smi || true)"
  if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi -L || true
  fi
} | tee "$OUT/backend_preflight.txt"
python - "$OUT" <<'PY'
import json, shutil, sys
from pathlib import Path
out=Path(sys.argv[1])
status="COLABFOLD_READY" if shutil.which("colabfold_batch") else "COLABFOLD_NOT_FOUND"
d={"phase":"Phase-8.2C backend preflight","status":status,
   "colabfold_batch":shutil.which("colabfold_batch"),"pymol":shutil.which("pymol"),"chimerax":shutil.which("chimerax"),
   "next":"run_phase82c1_colabfold_main.sh" if status=="COLABFOLD_READY" else "Install/activate a ColabFold-capable environment, then rerun this preflight. Do not substitute experimental-structure wording."}
(out/"BACKEND_PREFLIGHT.json").write_text(json.dumps(d,indent=2)+"\n")
print(json.dumps(d,indent=2))
PY
